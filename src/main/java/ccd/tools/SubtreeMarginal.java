package ccd.tools;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Runnable;
import beastfx.app.tools.Application;
import beastfx.app.treeannotator.TreeAnnotator;
import beastfx.app.util.OutFile;
import beastfx.app.util.TreeFile;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.IdentityHashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import ccd.algorithms.LoadOrStoreTrees;
import ccd.model.AbstractCCD;
import ccd.model.CCD0;
import ccd.model.CCD1;
import ccd.model.CCDType;
import ccd.model.Clade;
import ccd.model.CladePartition;
import ccd.model.ITreeDistribution;
import ccd.model.KRegCCD;
import ccd.model.bitsets.BitSet;

/**
 * Compares, at the resolution of small (3- or 4-taxon) induced subtrees, the
 * distributions over topologies implied by the empirical posterior sample and by
 * CCD0, CCD1 and (optimised) regCCD built from it.
 * <p>
 * Purpose: to establish that <em>CCD1 is an adequate basis</em> to build on. CCD1
 * conditions each split on its parent clade and so carries the correlation structure
 * of the posterior; regularisation then smooths that basis to give it full support.
 * The three models are therefore not competitors here. CCD0, which factorises over
 * clades and cannot represent those correlations, is the negative control that shows
 * the comparison is not trivial; regCCD is CCD1 plus smoothing, and is expected to
 * sit slightly further from the posterior than CCD1 does -- a small, uniform
 * perturbation is what smoothing <em>is</em>, and it is the price of the coverage
 * that {@code miss} measures. The question this tool answers is whether CCD1's
 * deviation from the posterior is small relative to the posterior's own sampling
 * noise, not which model wins on total variation.
 * <p>
 * Rationale: the marginal posterior probability of a full topology is not
 * Monte-Carlo-estimable for non-trivial analyses (almost every sampled tree is
 * unique). The induced distribution over the topologies of a fixed small taxon
 * subset is, however, low-dimensional and reliably estimable, and interpolates
 * between single-clade marginals and the full topology. This tool estimates those
 * induced distributions by sampling from each model and restricting every sampled
 * tree to each subset. The empirical side additionally reports a split-half
 * "noise floor" (posterior first half vs second half) so that a model's deviation
 * can be judged against the irreducible sampling scatter.
 * <p>
 * Output is a tab-separated table with one row per (subset, induced shape):
 * {@code subset_id, k, taxa, shape, emp, emp_h1, emp_h2, ccd0, ccd1, regccd}
 * where each of the last six columns is a probability (frequency) under that
 * distribution for that subset. Shapes not seen under a distribution get 0.
 */
@Description("Compare induced 3-/4-taxon subtree topology distributions between the "
        + "posterior sample and CCD0/CCD1/regCCD, by sampling and restriction.")
public class SubtreeMarginal extends Runnable {

    final public Input<TreeFile> treeInput = new Input<>("trees",
            "trees file to construct the CCDs from and analyse", Input.Validate.REQUIRED);
    final public Input<OutFile> outputInput = new Input<>("out",
            "tab-separated output file", Input.Validate.REQUIRED);
    final public Input<Integer> burnInPercentageInput = new Input<>("burnin",
            "percentage of trees used as burn-in (ignored)", 10);
    final public Input<String> subsetSizesInput = new Input<>("k",
            "comma-separated induced-subtree sizes to analyse", "3,4");
    final public Input<Integer> subsetsPerSizeInput = new Input<>("subsets",
            "number of random taxon subsets drawn per size k", 200);
    final public Input<Integer> sampleSizeInput = new Input<>("length",
            "number of trees sampled from each CCD model", 100000);
    final public Input<Integer> foldsInput = new Input<>("folds",
            "number of cross-validation folds for regCCD (KRegCCD) parameter selection", 5);
    final public Input<Boolean> modelsInput = new Input<>("models",
            "also build CCD0/CCD1/regCCD and compare; if false, only the empirical posterior is "
                    + "summarised (calibration mode)", true);
    final public Input<Boolean> exactInput = new Input<>("exact",
            "compute EXACT induced-subtree probabilities for CCD0 and CCD1 (sum-product DP over the "
                    + "clade DAG) instead of sampling; gives true coverage-miss (regCCD is 0 by the "
                    + "full-coverage theorem). Overrides -models.", false);
    final public Input<Integer> maxTreesInput = new Input<>("maxTrees",
            "if > 0, thin the post-burn-in chain evenly down to this many trees. Applied while "
                    + "reading, so only the kept trees are ever parsed or held in memory -- the "
                    + "whole point on files too large to materialise (a 441-taxon, 36k-tree set "
                    + "exhausts a 40 GB heap otherwise). All three models are then built from this "
                    + "same subsample", 0);
    final public Input<Integer> reserveDepthInput = new Input<>("reserveDepth",
            "reserve depth k for regCCD, used for both the cross-validated parameter search and the "
                    + "fitted model. k = 2 (the default) solves eps from N_1 and N_2; k = 1 skips the "
                    + "boundary-4 match, which is where nearly all the cost is on large clade sets -- "
                    + "at 441 taxa the default did not complete one CV fold in 11 h of CPU. k = 1 shifts eps "
                    + "by a few percent, so results are not strictly comparable across depths",
            KRegCCD.DEFAULT_RESERVE_DEPTH);
    final public Input<Long> seedInput = new Input<>("seed",
            "random seed (subset choice and model sampling)");
    final public Input<Integer> minSubsetsInput = new Input<>("minSubsets",
            "minimum number of informative subsets (ones where the posterior itself spreads over "
                    + "more than one induced shape) required before the CCD models are built and "
                    + "sampled. Guards against spending hours on a posterior with too little "
                    + "topological uncertainty to measure; set to 0 to force the run", 1);

    // Distribution column indices in the per-shape count arrays.
    private static final int EMP = 0, EMP_H1 = 1, EMP_H2 = 2, CCD0 = 3, CCD1 = 4, REGCCD = 5;
    private static final int NDIST = 6;
    private static final String[] COLS = {"emp", "emp_h1", "emp_h2", "ccd0", "ccd1", "regccd"};

    @Override
    public void initAndValidate() {
    }

    @Override
    public void run() throws Exception {
        long seed = (seedInput.get() != null) ? seedInput.get() : System.currentTimeMillis();
        Random random = new Random(seed);

        Log.info.println("# Induced subtree-distribution comparison");
        Log.info.println("    trees file:  " + treeInput.get().getPath());
        Log.info.println("    burnin:      " + burnInPercentageInput.get());
        Log.info.println("    sizes k:     " + subsetSizesInput.get());
        Log.info.println("    subsets/k:   " + subsetsPerSizeInput.get());
        Log.info.println("    #samples:    " + sampleSizeInput.get());
        Log.info.println("    seed:        " + seed);

        // Load the posterior sample. With -maxTrees the thinning happens DURING the read, so only
        // the kept trees are parsed and retained: materialising a large posterior first and
        // subsampling afterwards defeats the purpose, and on the biggest sets it simply runs out
        // of heap before it can get there.
        List<Tree> posterior;
        if (maxTreesInput.get() > 0) {
            posterior = LoadOrStoreTrees.loadTrees(treeInput.get(),
                    burnInPercentageInput.get() / 100.0, maxTreesInput.get());
            Log.info.println("    posterior:   " + posterior.size()
                    + " trees (evenly thinned over the post-burn-in chain)");
        } else {
            posterior = CCDToolUtil.treesFromSet(
                    CCDToolUtil.getTreeSet(treeInput, burnInPercentageInput.get()));
            Log.info.println("    posterior:   " + posterior.size() + " trees");
        }
        if (posterior.isEmpty()) {
            throw new IllegalArgumentException("No trees left after burn-in.");
        }

        // Taxon names from the first tree.
        List<String> taxa = new ArrayList<>();
        for (Node leaf : posterior.get(0).getExternalNodes()) {
            taxa.add(leaf.getID());
        }
        Collections.sort(taxa);
        Log.info.println("    taxa:        " + taxa.size());

        // Model construction is deferred until after the empirical tally, so that a posterior with
        // no topological uncertainty can be detected and the (expensive) model work skipped -- see
        // the pre-flight check below. Their RNG seeds are drawn HERE, at the point in the sequence
        // where the models used to be built, so that deferring construction leaves subset selection
        // and model sampling bit-identical to before.
        AbstractCCD ccd0 = null, ccd1 = null, regccd = null;
        final boolean wantModels = modelsInput.get() && !exactInput.get();
        long ccd0Seed = 0, ccd1Seed = 0, regccdSeed = 0;
        if (wantModels) {
            ccd0Seed = random.nextLong();
            ccd1Seed = random.nextLong();
            regccdSeed = random.nextLong();
        } else {
            Log.info.println("> empirical-only calibration mode (no CCD models built)");
        }

        // Draw the taxon subsets for each size k.
        List<Subset> subsets = new ArrayList<>();
        for (String tok : subsetSizesInput.get().split(",")) {
            int k = Integer.parseInt(tok.trim());
            if (k < 2 || k > taxa.size()) {
                Log.warning("skipping k=" + k + " (out of range for " + taxa.size() + " taxa)");
                continue;
            }
            drawSubsets(taxa, k, subsetsPerSizeInput.get(), random, subsets);
        }
        Log.info.println("    subsets:     " + subsets.size());

        if (exactInput.get()) {
            runExact(posterior, subsets);
            return;
        }

        // Per-subset shape counts across the six distributions.
        int n = subsets.size();
        List<Map<String, double[]>> counts = new ArrayList<>(n);
        double[] totals = new double[NDIST];
        for (int i = 0; i < n; i++) {
            counts.add(new LinkedHashMap<>());
        }

        // Empirical side: iterate every posterior tree (exact frequencies) and, by index
        // parity, feed the two split-half distributions used for the noise floor.
        Log.info.println("> tallying empirical posterior (" + posterior.size() + " trees)...");
        int empReport = Math.max(1, posterior.size() / 10);
        for (int t = 0; t < posterior.size(); t++) {
            Tree tree = posterior.get(t);
            Node root = tree.getRoot();
            Map<String, Node> leaves = leafMap(tree);
            int half = (t % 2 == 0) ? EMP_H1 : EMP_H2;
            for (int i = 0; i < n; i++) {
                String shape = inducedShape(root, leaves, subsets.get(i).taxaSet);
                double[] c = counts.get(i).computeIfAbsent(shape, s -> new double[NDIST]);
                c[EMP]++;
                c[half]++;
            }
            if ((t + 1) % empReport == 0) {
                Log.info.println("    empirical: " + (t + 1) + "/" + posterior.size());
            }
        }
        totals[EMP] = posterior.size();
        totals[EMP_H1] = (posterior.size() + 1) / 2;
        totals[EMP_H2] = posterior.size() / 2;

        // Pre-flight: a subset carries signal only if the POSTERIOR itself spreads over more than
        // one induced shape. If the posterior is (near-)certain at this k, every subset is a point
        // mass, the downstream analysis discards all of them, and the model work -- regCCD
        // parameter selection plus three 50k-tree sampling runs, hours on a large dataset -- would
        // produce a table with nothing in it. Check before paying for it. (Encountered on
        // tornabene-2016: 1801 posterior trees, all the SAME topology, so 0 of 1000 subsets were
        // usable at k=4 and no larger k would have helped either.)
        int nonDegenerate = 0;
        for (Map<String, double[]> c : counts) {
            int shapesSeen = 0;
            for (double[] v : c.values()) {
                if (v[EMP] > 0) {
                    shapesSeen++;
                }
            }
            if (shapesSeen >= 2) {
                nonDegenerate++;
            }
        }
        Log.info.println("    distinct posterior topologies: " + countDistinctTopologies(posterior)
                + " of " + posterior.size() + " trees");
        Log.info.println("    informative subsets (posterior spreads over >1 shape): "
                + nonDegenerate + " / " + n);

        boolean buildModels = wantModels;
        if (wantModels && nonDegenerate < minSubsetsInput.get()) {
            Log.warning("SKIPPING the model comparison: only " + nonDegenerate + " of " + n
                    + " subsets are informative (need >= " + minSubsetsInput.get() + ").");
            Log.warning("The posterior has too little topological uncertainty at this k for the "
                    + "induced-subtree comparison to measure anything. Writing the empirical "
                    + "columns only; model columns will be 0. Raise -minSubsets to force the run.");
            buildModels = false;
        }

        if (buildModels) {
            Log.info.println("> building models...");
            // Built from the same tree list as regCCD, so that under -maxTrees all three models
            // see the identical subsample (the tree set would have re-read the whole file).
            ccd0 = new CCD0(posterior, 0.0);
            ccd1 = new CCD1(posterior, 0.0);
            Log.info.println("    selecting regCCD (KRegCCD) parameters by " + foldsInput.get()
                    + "-fold cross-validation...");
            regccd = KRegCCD.withOptimisedParameters(posterior, foldsInput.get(),
                    reserveDepthInput.get());
            Log.info.println("    " + regccd);
            ccd0.setRandom(new Random(ccd0Seed));
            ccd1.setRandom(new Random(ccd1Seed));
            regccd.setRandom(new Random(regccdSeed));
        }

        // Model side: sample from each CCD and restrict.
        if (buildModels) {
            tallyModel(ccd0, CCD0, sampleSizeInput.get(), subsets, counts, "CCD0");
            tallyModel(ccd1, CCD1, sampleSizeInput.get(), subsets, counts, "CCD1");
            tallyModel(regccd, REGCCD, sampleSizeInput.get(), subsets, counts, "regCCD");
            totals[CCD0] = totals[CCD1] = totals[REGCCD] = sampleSizeInput.get();
        }

        // Write the table.
        Log.info.println("> writing " + outputInput.get().getPath());
        try (PrintStream out = new PrintStream(outputInput.get())) {
            out.print("subset_id\tk\ttaxa\tshape");
            for (String col : COLS) {
                out.print("\t" + col);
            }
            out.println();
            for (int i = 0; i < n; i++) {
                Subset s = subsets.get(i);
                String taxaStr = String.join(",", s.taxaList);
                for (Map.Entry<String, double[]> e : counts.get(i).entrySet()) {
                    double[] c = e.getValue();
                    out.print(i + "\t" + s.taxaList.size() + "\t" + taxaStr + "\t" + e.getKey());
                    for (int d = 0; d < NDIST; d++) {
                        double p = (totals[d] > 0) ? c[d] / totals[d] : 0.0;
                        out.print("\t" + p);
                    }
                    out.println();
                }
            }
        }
        Log.info.println("... done.");
    }

    /**
     * Number of distinct topologies in the sample, ignoring branch lengths: a one-line diagnostic
     * for why a posterior yielded few informative subsets. A value of 1 means the chain returned a
     * single tree, in which case no choice of k can give this comparison anything to measure.
     */
    private static int countDistinctTopologies(List<Tree> posterior) {
        Set<String> seen = new HashSet<>();
        for (Tree t : posterior) {
            List<String> clades = new ArrayList<>();
            collectClades(t.getRoot(), clades);
            Collections.sort(clades);
            seen.add(String.join("|", clades));
        }
        return seen.size();
    }

    /** Appends a canonical string for every internal clade of the subtree at {@code v}. */
    private static List<Integer> collectClades(Node v, List<String> out) {
        List<Integer> taxa = new ArrayList<>();
        if (v.isLeaf()) {
            taxa.add(v.getNr());
        } else {
            for (Node child : v.getChildren()) {
                taxa.addAll(collectClades(child, out));
            }
            Collections.sort(taxa);
            out.add(taxa.toString());
        }
        return taxa;
    }

    /** Sample {@code length} trees from the distribution and tally induced shapes per subset. */
    private void tallyModel(ITreeDistribution dist, int col, int length,
                            List<Subset> subsets, List<Map<String, double[]>> counts, String label) {
        Log.info.println("> sampling " + length + " trees from " + label + "...");
        int report = Math.max(1, length / 10);
        for (int t = 0; t < length; t++) {
            Tree tree = dist.sampleTree();
            Node root = tree.getRoot();
            Map<String, Node> leaves = leafMap(tree);
            for (int i = 0; i < subsets.size(); i++) {
                String shape = inducedShape(root, leaves, subsets.get(i).taxaSet);
                double[] c = counts.get(i).computeIfAbsent(shape, s -> new double[NDIST]);
                c[col]++;
            }
            if ((t + 1) % report == 0) {
                Log.info.println("    " + label + ": " + (t + 1) + "/" + length);
            }
        }
    }

    /** Leaf-label to leaf-node map for one tree (built once, reused across all subsets). */
    private static Map<String, Node> leafMap(Tree tree) {
        Map<String, Node> map = new HashMap<>();
        for (Node leaf : tree.getExternalNodes()) {
            map.put(leaf.getID(), leaf);
        }
        return map;
    }

    /**
     * Canonical labelled rooted topology induced on {@code keep}, computed without copying or
     * mutating the tree. Marks the kept-leaf-to-root paths (O(k * height)), then reads the induced
     * topology off the marked nodes, suppressing any node that subtends only one kept leaf.
     */
    private static String inducedShape(Node root, Map<String, Node> leaves, Set<String> keep) {
        Map<Node, Integer> below = new IdentityHashMap<>();
        for (String taxon : keep) {
            Node x = leaves.get(taxon);
            while (x != null) {
                below.merge(x, 1, Integer::sum);
                x = x.getParent();
            }
        }
        return inducedKey(root, below);
    }

    /** Recurse through nodes subtending a kept leaf; a node with a single such child is suppressed. */
    private static String inducedKey(Node node, Map<Node, Integer> below) {
        if (node.isLeaf()) {
            return node.getID();
        }
        List<Node> relevant = new ArrayList<>();
        for (Node child : node.getChildren()) {
            if (below.getOrDefault(child, 0) > 0) {
                relevant.add(child);
            }
        }
        if (relevant.size() == 1) {
            return inducedKey(relevant.get(0), below);
        }
        List<String> parts = new ArrayList<>(relevant.size());
        for (Node child : relevant) {
            parts.add(inducedKey(child, below));
        }
        Collections.sort(parts);
        return "(" + String.join(",", parts) + ")";
    }

    /** Draw {@code count} distinct random size-{@code k} subsets of {@code taxa}. */
    private static void drawSubsets(List<String> taxa, int k, int count, Random random, List<Subset> out) {
        Set<String> seen = new HashSet<>();
        int n = taxa.size();
        int attempts = 0, maxAttempts = count * 50 + 1000;
        while (seen.size() < count && attempts++ < maxAttempts) {
            List<Integer> idx = new ArrayList<>();
            for (int i = 0; i < n; i++) {
                idx.add(i);
            }
            Collections.shuffle(idx, random);
            List<String> chosen = new ArrayList<>();
            for (int i = 0; i < k; i++) {
                chosen.add(taxa.get(idx.get(i)));
            }
            Collections.sort(chosen);
            String key = String.join(",", chosen);
            if (seen.add(key)) {
                out.add(new Subset(chosen));
            }
        }
    }

    /** A taxon subset, both as an ordered label list and as a lookup set. */
    private static final class Subset {
        final List<String> taxaList;
        final Set<String> taxaSet;

        Subset(List<String> taxaList) {
            this.taxaList = taxaList;
            this.taxaSet = new HashSet<>(taxaList);
        }
    }

    // ===== Exact induced-subtree probabilities for CCD0/CCD1 (sum-product DP over the clade DAG) =====

    /**
     * Exact mode. For each subset, tabulate the empirical induced-shape distribution and the EXACT
     * induced probability under CCD0 and CCD1. regCCD is not evaluated: its induced probability is
     * positive for every shape (full-coverage theorem), so its coverage-miss is 0 by construction.
     * Output: {@code subset_id, k, taxa, shape, emp, ccd0_exact, ccd1_exact, ccd0_massObs, ccd1_massObs}
     * where the mass columns are the total CCD probability on the subset's observed shapes, so the
     * analysis can form exact TV via the tail term (1 - massObs).
     */
    private void runExact(List<Tree> posterior, List<Subset> subsets) throws Exception {
        Log.info.println("> exact induced-probability mode (CCD0, CCD1)");
        AbstractCCD ccd0 = new CCD0(posterior, 0.0);
        AbstractCCD ccd1 = new CCD1(posterior, 0.0);
        String[] taxaNames = ccd1.getSomeBaseTree().getTaxaNames();
        Map<String, Integer> nameToIdx = new HashMap<>();
        for (int i = 0; i < taxaNames.length; i++) {
            nameToIdx.put(taxaNames[i], i);
        }

        int report = Math.max(1, subsets.size() / 10);
        try (java.io.PrintStream out = new java.io.PrintStream(outputInput.get())) {
            out.println("subset_id\tk\ttaxa\tshape\temp\tccd0_exact\tccd1_exact\tccd0_massObs\tccd1_massObs");
            double nTot = posterior.size();
            for (int si = 0; si < subsets.size(); si++) {
                Subset s = subsets.get(si);
                BitSet sBits = BitSet.newBitSet(taxaNames.length);
                for (String name : s.taxaList) {
                    sBits.set(nameToIdx.get(name));
                }
                // empirical induced distribution + one representative clade-map per distinct shape
                Map<String, Integer> shapeCount = new LinkedHashMap<>();
                Map<String, Map<BitSet, BitSet[]>> shapeTarget = new HashMap<>();
                for (Tree tree : posterior) {
                    Map<String, Node> leaves = leafMap(tree);
                    Map<Node, Integer> below = markBelow(tree.getRoot(), leaves, s.taxaSet);
                    Map<BitSet, BitSet[]> clades = new HashMap<>();
                    targetClades(tree.getRoot(), below, nameToIdx, taxaNames.length, clades);
                    String key = shapeKeyFromClades(clades);
                    shapeCount.merge(key, 1, Integer::sum);
                    shapeTarget.putIfAbsent(key, clades);
                }
                // exact CCD probabilities per observed shape (restricted-bits cache is target-independent)
                Map<Clade, BitSet> rb0 = new IdentityHashMap<>();
                Map<Clade, BitSet> rb1 = new IdentityHashMap<>();
                Map<String, double[]> exact = new HashMap<>();
                double mass0 = 0, mass1 = 0;
                for (Map.Entry<String, Map<BitSet, BitSet[]>> e : shapeTarget.entrySet()) {
                    double p0 = exactInduced(ccd0, sBits, e.getValue(), rb0);
                    double p1 = exactInduced(ccd1, sBits, e.getValue(), rb1);
                    exact.put(e.getKey(), new double[]{p0, p1});
                    mass0 += p0;
                    mass1 += p1;
                }
                String taxaStr = String.join(",", s.taxaList);
                for (Map.Entry<String, Integer> e : shapeCount.entrySet()) {
                    double[] px = exact.get(e.getKey());
                    out.println(si + "\t" + s.taxaList.size() + "\t" + taxaStr + "\t" + e.getKey()
                            + "\t" + (e.getValue() / nTot) + "\t" + px[0] + "\t" + px[1]
                            + "\t" + mass0 + "\t" + mass1);
                }
                if ((si + 1) % report == 0) {
                    Log.info.println("    exact: " + (si + 1) + "/" + subsets.size());
                }
            }
        }
        Log.info.println("... done.");
    }

    /** Mark, for every ancestor of a kept leaf, how many kept leaves lie below it. */
    private static Map<Node, Integer> markBelow(Node root, Map<String, Node> leaves, Set<String> keep) {
        Map<Node, Integer> below = new IdentityHashMap<>();
        for (String taxon : keep) {
            Node x = leaves.get(taxon);
            while (x != null) {
                below.merge(x, 1, Integer::sum);
                x = x.getParent();
            }
        }
        return below;
    }

    /**
     * Build the induced target topology as a map from each internal induced clade (as a BitSet over
     * CCD taxon indices) to its two child clades. Returns the clade below {@code node}.
     */
    private static BitSet targetClades(Node node, Map<Node, Integer> below, Map<String, Integer> nameToIdx,
                                       int nTaxa, Map<BitSet, BitSet[]> cladeChildren) {
        if (node.isLeaf()) {
            BitSet b = BitSet.newBitSet(nTaxa);
            b.set(nameToIdx.get(node.getID()));
            return b;
        }
        List<Node> rel = new ArrayList<>();
        for (Node child : node.getChildren()) {
            if (below.getOrDefault(child, 0) > 0) {
                rel.add(child);
            }
        }
        if (rel.size() == 1) {
            return targetClades(rel.get(0), below, nameToIdx, nTaxa, cladeChildren);
        }
        BitSet mine = BitSet.newBitSet(nTaxa);
        BitSet[] kids = new BitSet[rel.size()];
        for (int i = 0; i < rel.size(); i++) {
            kids[i] = targetClades(rel.get(i), below, nameToIdx, nTaxa, cladeChildren);
            mine.or(kids[i]);
        }
        cladeChildren.put(mine, new BitSet[]{kids[0], kids[1]});
        return mine;
    }

    /** Canonical key of an induced topology: the sorted set of its internal induced clades. */
    private static String shapeKeyFromClades(Map<BitSet, BitSet[]> cladeChildren) {
        List<String> keys = new ArrayList<>(cladeChildren.size());
        for (BitSet b : cladeChildren.keySet()) {
            keys.add(b.toString());
        }
        Collections.sort(keys);
        return String.join("|", keys);
    }

    /** Exact probability that the CCD induces {@code cladeChildren} on the subset {@code sBits}. */
    private static double exactInduced(AbstractCCD ccd, BitSet sBits,
                                       Map<BitSet, BitSet[]> cladeChildren, Map<Clade, BitSet> rbCache) {
        return gExact(ccd.getRootClade(), sBits, cladeChildren, rbCache, new IdentityHashMap<>());
    }

    /**
     * Sum-product recursion: probability that the CCD subtree below clade {@code c} induces, on
     * {@code c}'s S-restriction, the corresponding subtree of the target. Memoised per target.
     */
    private static double gExact(Clade c, BitSet sBits, Map<BitSet, BitSet[]> cladeChildren,
                                 Map<Clade, BitSet> rbCache, Map<Clade, Double> memo) {
        Double cached = memo.get(c);
        if (cached != null) {
            return cached;
        }
        BitSet a = restricted(c, sBits, rbCache);
        double result;
        if (a.cardinality() <= 1) {
            result = 1.0;                          // empty or singleton restriction: trivially induced
        } else {
            BitSet[] kids = cladeChildren.get(a);
            if (kids == null) {
                result = 0.0;                      // c's restriction is not a clade of the target
            } else {
                BitSet left = kids[0], right = kids[1];
                double total = 0.0;
                for (CladePartition part : c.getPartitions()) {
                    Clade[] cc = part.getChildClades();
                    BitSet a1 = restricted(cc[0], sBits, rbCache);
                    BitSet a2 = restricted(cc[1], sBits, rbCache);
                    double p = part.getCCP();
                    if (a1.isEmpty()) {
                        total += p * gExact(cc[1], sBits, cladeChildren, rbCache, memo);
                    } else if (a2.isEmpty()) {
                        total += p * gExact(cc[0], sBits, cladeChildren, rbCache, memo);
                    } else if ((a1.equals(left) && a2.equals(right))
                            || (a1.equals(right) && a2.equals(left))) {
                        total += p * gExact(cc[0], sBits, cladeChildren, rbCache, memo)
                                * gExact(cc[1], sBits, cladeChildren, rbCache, memo);
                    }
                }
                result = total;
            }
        }
        memo.put(c, result);
        return result;
    }

    /** S-restriction of a clade's taxon set, cached per clade (target-independent). */
    private static BitSet restricted(Clade c, BitSet sBits, Map<Clade, BitSet> cache) {
        BitSet r = cache.get(c);
        if (r == null) {
            r = (BitSet) c.getCladeInBitsTaxaOnly().clone();
            r.and(sBits);
            cache.put(c, r);
        }
        return r;
    }

    public static void main(String[] args) throws Exception {
        new Application(new SubtreeMarginal(), "Subtree Marginal Comparison", args);
    }
}
