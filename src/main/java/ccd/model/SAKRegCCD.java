package ccd.model;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator;
import ccd.model.bitsets.BitSet;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Full-support CCD for sampled-ancestor (SA) trees: {@link KRegCCD} on the topology, times a
 * smoothed model of the SA <em>split type</em> at every node.
 *
 * <h2>SA-ness as a split type</h2>
 * In the BEAST representation an SA tree is a binary tree in which a leaf {@code i} with a
 * zero-length branch is a sampled ancestor of its sibling subtree. So an SA tree is exactly a
 * (plain) rooted topology plus, at every internal node {@code v}, a <em>label</em>
 * <pre>
 *   t_v in {NONE} ∪ { SA(i) : i is a leaf child of v }
 * </pre>
 * (a cherry has three types: none, left SA, right SA; a node with one leaf child has two; a node
 * with no leaf child has only NONE). Labels never interact across nodes — each leaf is the child of
 * exactly one node — so any labelling of any topology is a valid SA tree, and
 * <pre>
 *   P(T) = P_topo(topology(T)) * prod_v P(t_v | context_v)
 * </pre>
 * is normalised as long as each per-node label distribution is normalised over that node's
 * feasible labels (which depend only on the node's split).
 *
 * <h2>Relation to CCD1EL / CCD1EI</h2>
 * <ul>
 *   <li><b>EL</b>: the context of {@code t_v} is the split {@code (C_v, A|B)}. With a plain CCD1
 *       topology and no smoothing this is <em>exactly</em> {@link CCD1EL}: an EL partition
 *       {@code ({i}_SA, B)} is the split {@code {i}|B} with label {@code SA(i)}, and the EL CCP
 *       factorises as {@code n(C,A|B)/n(C) * n(C,A|B,t)/n(C,A|B)}.</li>
 *   <li><b>EI</b>: in {@link CCD1EI} a clade's identity carries the SA flag of its own top node,
 *       so the parent picks its children's labels jointly. Here the parent {@code p} draws the pair
 *       {@code (t_A, t_B)} of its two children's labels conditional on
 *       {@code (C_p, t_p, A|B, split(A), split(B))}; the root's label is drawn from the root split.
 *       The topology is still scored by KRegCCD (it does not condition on {@code t_p}), so this is
 *       a "topology-first" EI. On the four-taxon worked example it reproduces CCD1EI's 27/27/27.</li>
 * </ul>
 *
 * <h2>Smoothing (full support)</h2>
 * Topology: KRegCCD (split expansion + additive {@code alpha} + novel-clade escape {@code mu}).
 * Labels: a back-off chain of additive smoothers, each level's pseudocount {@code beta} being
 * spread according to the next, coarser level:
 * <pre>
 *   P_EI(tA,tB | ctx)   = (n_EI(ctx,tA,tB) + beta * P_EL(tA|A) * P_EL(tB|B)) / (n_EI(ctx) + beta)
 *   P_EL(t | C, A|B)    = (n_EL(C,A|B,t)   + beta * q(t | A|B))              / (n_EL(C,A|B) + beta)
 *   q(t | A|B)          ∝ per-leaf independent prior, renormalised over the node's feasible labels:
 *                         q(NONE) ∝ prod_{leaf children j}(1 - pi_j),
 *                         q(SA(i)) ∝ pi_i * prod_{leaf children j != i}(1 - pi_j)
 *   pi_i                = (s_i + gamma * pibar) / (N + gamma)
 * </pre>
 * where {@code s_i} is the number of training trees in which leaf {@code i} is an SA, {@code N}
 * the number of trees, and {@code pibar = (sum_i s_i + 1) / (nN + 2)} the pooled SA rate (with a
 * Laplace pseudocount so it is never 0). Every {@code pi_i} lies strictly in (0, 1), hence every
 * label of every node — observed or novel — has positive probability: full support over SA trees.
 * Unobserved splits (including KRegCCD's novel clades) have {@code n = 0} and fall straight through
 * to {@code q}.
 *
 * <p>Notes: removing trees does not update the label counts (KRegCCD's held-out / parameter
 * optimisation applies to the topology part only). {@link #getNumberOfTrees()} is a coarse upper
 * bound on the SA-tree support: rooted topologies times {@code 2^n} labellings.
 */
public class SAKRegCCD extends KRegCCD {

    /** Which context the SA label of a node is conditioned on. */
    public enum SAContext {
        /** Label conditioned on the node's own split (CCD1EL-like). */
        EL,
        /** Children's labels drawn jointly at the parent, conditioned on the parent's labelled split (CCD1EI-like). */
        EI
    }

    public static final double DEFAULT_BETA = 1.0;
    public static final double DEFAULT_GAMMA = 1.0;

    /** Label value for "no sampled ancestor at this node". */
    static final int NONE = -1;

    /* Label statistics. NOT initialised in field declarations: the super constructor cladifies the
     * trees (through cladifyTree) before subclass field initialisers would run, and an initialiser
     * would then wipe the counts. They are created lazily in ensureLabelState(). */
    private Map<ElKey, Map<Integer, Double>> elCounts;
    private Map<EiKey, Map<Long, Double>> eiCounts;
    private int[] leafSACounts;
    private int numLabelledTrees;

    private final SAContext context;
    private final double beta;
    private final double gamma;

    /* Derived per-leaf prior pi_i, computed on first use (after all trees are in). */
    private volatile double[] leafPrior;

    /* -- KEYS -- */

    /** Split identity: parent clade and its canonical first child (taxa bits only). */
    private record ElKey(BitSet clade, BitSet firstChild) {
    }

    /**
     * EI context at parent {@code p}: p's split, p's own label, and both children's splits
     * ({@code null} first-child for a leaf child).
     */
    private record EiKey(BitSet clade, BitSet firstChild, int parentLabel,
                         BitSet firstChildSplit, BitSet secondChildSplit) {
    }

    /* -- CONSTRUCTORS -- */

    public SAKRegCCD(List<Tree> trees, double burnin, SAContext context) {
        this(trees, burnin, DEFAULT_MU, DEFAULT_ALPHA, DEFAULT_RESERVE_DEPTH, TailMode.BOUND,
                context, DEFAULT_BETA, DEFAULT_GAMMA);
    }

    /**
     * @param mu, alpha, k, tailMode KRegCCD topology parameters (see {@link KRegCCD})
     * @param context which SA-label context to use (see {@link SAContext})
     * @param beta    back-off pseudocount for the label smoothers (&gt; 0)
     * @param gamma   pseudocount shrinking each leaf's SA rate towards the pooled rate (&gt; 0)
     */
    public SAKRegCCD(List<Tree> trees, double burnin, double mu, double alpha, int k,
                     TailMode tailMode, SAContext context, double beta, double gamma) {
        super(trees, burnin, mu, alpha, k, tailMode);
        validate(beta, gamma);
        this.context = context;
        this.beta = beta;
        this.gamma = gamma;
    }

    public SAKRegCCD(TreeAnnotator.TreeSet treeSet, SAContext context) {
        this(treeSet, DEFAULT_MU, DEFAULT_ALPHA, DEFAULT_RESERVE_DEPTH, TailMode.BOUND,
                context, DEFAULT_BETA, DEFAULT_GAMMA);
    }

    public SAKRegCCD(TreeAnnotator.TreeSet treeSet, double mu, double alpha, int k,
                     TailMode tailMode, SAContext context, double beta, double gamma) {
        super(treeSet, mu, alpha, k, tailMode);
        validate(beta, gamma);
        this.context = context;
        this.beta = beta;
        this.gamma = gamma;
    }

    private static void validate(double beta, double gamma) {
        if (!(beta > 0) || !(gamma > 0)) {
            throw new IllegalArgumentException("beta and gamma must be > 0 for full support, got beta="
                    + beta + ", gamma=" + gamma);
        }
    }

    public SAContext getContext() {
        return context;
    }

    public double getBeta() {
        return beta;
    }

    public double getGamma() {
        return gamma;
    }

    /* -- TREE INSERTION: record label statistics alongside the topology -- */

    @Override
    protected void cladifyTree(Tree tree) {
        super.cladifyTree(tree);
        ensureLabelState();
        Map<Node, BitSet> bits = new HashMap<>();
        computeTaxaBits(tree.getRoot(), bits);
        recordLabels(tree.getRoot(), true, bits);
        for (Node leaf : tree.getExternalNodes()) {
            if (isSALeaf(leaf)) {
                leafSACounts[leaf.getNr()]++;
            }
        }
        numLabelledTrees++;
        leafPrior = null;
    }

    private void ensureLabelState() {
        if (elCounts == null) {
            elCounts = new HashMap<>();
            eiCounts = new HashMap<>();
            leafSACounts = new int[leafArraySize];
        }
    }

    /** Records the EL and EI label counts of {@code v} and its subtree. */
    private void recordLabels(Node v, boolean isRoot, Map<Node, BitSet> bits) {
        if (v.isLeaf()) {
            return;
        }
        int tv = observedLabel(v);
        add(elCounts, elKey(v, bits), tv);
        if (isRoot) {
            // root label in EI: no parent, condition on the root split only (same key as EL,
            // stored in the EI map under a sentinel parent label so the two stay separate)
            add(eiCounts, rootEiKey(v, bits), pairCode(tv, NONE));
        }
        Node[] kids = orderedChildren(v, bits);
        int tA = kids[0].isLeaf() ? NONE : observedLabel(kids[0]);
        int tB = kids[1].isLeaf() ? NONE : observedLabel(kids[1]);
        add(eiCounts, eiKey(v, tv, kids, bits), pairCode(tA, tB));
        recordLabels(kids[0], false, bits);
        recordLabels(kids[1], false, bits);
    }

    private static <K, E> void add(Map<K, Map<E, Double>> m, K key, E event) {
        m.computeIfAbsent(key, x -> new HashMap<>()).merge(event, 1.0, Double::sum);
    }

    /* -- TREE PROBABILITY -- */

    @Override
    public double getLogProbabilityOfTree(Tree tree) {
        return super.getLogProbabilityOfTree(tree) + labelLogProbability(tree);
    }

    @Override
    public double getLogProbabilityOfTree(Tree tree, double scoreMu) {
        return super.getLogProbabilityOfTree(tree, scoreMu) + labelLogProbability(tree);
    }

    @Override
    public double getLogProbabilityOfTree(Tree tree, Map<Node, BitSet> scratch) {
        return super.getLogProbabilityOfTree(tree, scratch) + labelLogProbability(tree);
    }

    @Override
    public BigInteger getNumberOfTrees() {
        // every topology admits at most 2^n labellings (each leaf SA or not), many infeasible;
        // a coarse upper bound on the SA-tree support
        return super.getNumberOfTrees().shiftLeft(getNumberOfLeaves());
    }

    /**
     * Log-probability of an SA tree with the labels scored under {@code labelContext} instead of
     * this model's own context. Both the EL and EI label tables are always recorded, so one model
     * can score under either (the topology part is identical).
     */
    public double getLogProbabilityOfTree(Tree tree, SAContext labelContext) {
        return super.getLogProbabilityOfTree(tree) + labelLogProbability(tree, labelContext);
    }

    /** {@code log P(labels | topology)} of an SA tree under this model. */
    public double labelLogProbability(Tree tree) {
        return labelLogProbability(tree, context);
    }

    /** {@code log P(labels | topology)} of an SA tree, with labels scored under {@code labelContext}. */
    public double labelLogProbability(Tree tree, SAContext labelContext) {
        Map<Node, BitSet> bits = new HashMap<>();
        computeTaxaBits(tree.getRoot(), bits);
        Node root = tree.getRoot();
        if (root.isLeaf()) {
            return 0.0;
        }
        if (labelContext == SAContext.EL) {
            double logp = 0.0;
            for (Node v : tree.getInternalNodes()) {
                logp += Math.log(elProb(v, observedLabel(v), bits));
            }
            return logp;
        }
        double logp = Math.log(eiRootProb(root, observedLabel(root), bits));
        for (Node p : tree.getInternalNodes()) {
            Node[] kids = orderedChildren(p, bits);
            int tA = kids[0].isLeaf() ? NONE : observedLabel(kids[0]);
            int tB = kids[1].isLeaf() ? NONE : observedLabel(kids[1]);
            logp += Math.log(eiProb(p, observedLabel(p), kids, tA, tB, bits));
        }
        return logp;
    }

    /* -- LABEL PROBABILITIES -- */

    /** Leaf-prior q(t | split of v), normalised over v's feasible labels. */
    private double priorProb(Node v, int t) {
        double[] pi = leafPrior();
        double wNone = 1.0;
        for (Node c : v.getChildren()) {
            if (c.isLeaf()) {
                wNone *= 1.0 - pi[c.getNr()];
            }
        }
        double z = wNone;
        double wt = (t == NONE) ? wNone : 0.0;
        for (Node c : v.getChildren()) {
            if (c.isLeaf()) {
                double w = wNone / (1.0 - pi[c.getNr()]) * pi[c.getNr()];
                z += w;
                if (c.getNr() == t) {
                    wt = w;
                }
            }
        }
        return wt / z;
    }

    private double elProb(Node v, int t, Map<Node, BitSet> bits) {
        Map<Integer, Double> counts = elCounts.get(elKey(v, bits));
        double q = priorProb(v, t);
        if (counts == null) {
            return q;
        }
        return (counts.getOrDefault(t, 0.0) + beta * q) / (total(counts) + beta);
    }

    /** Back-off for the joint child labels: product of the two children's EL probabilities. */
    private double elPairProb(Node[] kids, int tA, int tB, Map<Node, BitSet> bits) {
        double pA = kids[0].isLeaf() ? (tA == NONE ? 1.0 : 0.0) : elProb(kids[0], tA, bits);
        double pB = kids[1].isLeaf() ? (tB == NONE ? 1.0 : 0.0) : elProb(kids[1], tB, bits);
        return pA * pB;
    }

    private double eiProb(Node p, int tp, Node[] kids, int tA, int tB, Map<Node, BitSet> bits) {
        double backoff = elPairProb(kids, tA, tB, bits);
        Map<Long, Double> counts = eiCounts.get(eiKey(p, tp, kids, bits));
        if (counts == null) {
            return backoff;
        }
        return (counts.getOrDefault(pairCode(tA, tB), 0.0) + beta * backoff) / (total(counts) + beta);
    }

    private double eiRootProb(Node root, int t, Map<Node, BitSet> bits) {
        double backoff = elProb(root, t, bits);
        Map<Long, Double> counts = eiCounts.get(rootEiKey(root, bits));
        if (counts == null) {
            return backoff;
        }
        return (counts.getOrDefault(pairCode(t, NONE), 0.0) + beta * backoff) / (total(counts) + beta);
    }

    private static double total(Map<?, Double> counts) {
        double s = 0.0;
        for (double c : counts.values()) {
            s += c;
        }
        return s;
    }

    /** pi_i: per-leaf SA rate shrunk towards the pooled rate; strictly inside (0, 1). */
    private double[] leafPrior() {
        double[] pi = leafPrior;
        if (pi != null) {
            return pi;
        }
        ensureLabelState();
        int n = leafArraySize;
        long totalSA = 0;
        for (int s : leafSACounts) {
            totalSA += s;
        }
        double pibar = (totalSA + 1.0) / ((double) n * numLabelledTrees + 2.0);
        pi = new double[n];
        for (int i = 0; i < n; i++) {
            pi[i] = (leafSACounts[i] + gamma * pibar) / (numLabelledTrees + gamma);
            // s_i <= N, so pi_i < 1 because pibar < 1; guard anyway for the w/(1-pi) division
            pi[i] = Math.min(pi[i], 1.0 - 1e-12);
        }
        leafPrior = pi;
        return pi;
    }

    /** The fitted per-leaf SA prior pi_i (copy). */
    public double[] getLeafSAPrior() {
        return leafPrior().clone();
    }

    /* -- SAMPLING & MAP: topology from KRegCCD, then labels -- */

    @Override
    protected Tree getTreeBasedOnStrategy(SamplingStrategy samplingStrategy,
                                          HeightSettingStrategy heightStrategy) {
        // An SA is encoded by a zero-length leaf branch, which needs heights; without them every
        // leaf would read as an SA. Fall back to unit heights.
        HeightSettingStrategy hs = (heightStrategy == null || heightStrategy == HeightSettingStrategy.None)
                ? HeightSettingStrategy.One : heightStrategy;
        Tree tree = super.getTreeBasedOnStrategy(samplingStrategy, hs);
        Node root = tree.getRoot();
        if (root.isLeaf()) {
            return tree;
        }
        Map<Node, BitSet> bits = new HashMap<>();
        computeTaxaBits(root, bits);
        Map<Node, Integer> labels = new HashMap<>();
        if (samplingStrategy == SamplingStrategy.MAP) {
            assignMapLabels(tree, bits, labels);
        } else {
            assignSampledLabels(tree, bits, labels);
        }
        applyLabels(root, labels);

        double labelLog = labelLogProbabilityOf(tree, bits, labels);
        Object lp = root.getMetaData(LOG_PROB_SUBTREE_KEY);
        if (lp instanceof Double d) {
            root.setMetaData(LOG_PROB_SUBTREE_KEY, d + labelLog);
            root.setMetaData(PROB_SUBTREE_KEY, Math.exp(d + labelLog));
        }
        return tree;
    }

    /** Feasible labels of v: NONE, then each leaf child's number. */
    private static List<Integer> feasibleLabels(Node v) {
        List<Integer> f = new ArrayList<>(3);
        f.add(NONE);
        for (Node c : v.getChildren()) {
            if (c.isLeaf()) {
                f.add(c.getNr());
            }
        }
        return f;
    }

    private void assignSampledLabels(Tree tree, Map<Node, BitSet> bits, Map<Node, Integer> labels) {
        Node root = tree.getRoot();
        if (context == SAContext.EL) {
            for (Node v : tree.getInternalNodes()) {
                List<Integer> f = feasibleLabels(v);
                double[] w = new double[f.size()];
                for (int i = 0; i < w.length; i++) {
                    w[i] = elProb(v, f.get(i), bits);
                }
                labels.put(v, f.get(draw(w)));
            }
            return;
        }
        List<Integer> f = feasibleLabels(root);
        double[] w = new double[f.size()];
        for (int i = 0; i < w.length; i++) {
            w[i] = eiRootProb(root, f.get(i), bits);
        }
        labels.put(root, f.get(draw(w)));
        sampleChildrenEI(root, bits, labels);
    }

    private void sampleChildrenEI(Node p, Map<Node, BitSet> bits, Map<Node, Integer> labels) {
        Node[] kids = orderedChildren(p, bits);
        List<Integer> fA = kids[0].isLeaf() ? List.of(NONE) : feasibleLabels(kids[0]);
        List<Integer> fB = kids[1].isLeaf() ? List.of(NONE) : feasibleLabels(kids[1]);
        int tp = labels.get(p);
        double[] w = new double[fA.size() * fB.size()];
        for (int a = 0; a < fA.size(); a++) {
            for (int b = 0; b < fB.size(); b++) {
                w[a * fB.size() + b] = eiProb(p, tp, kids, fA.get(a), fB.get(b), bits);
            }
        }
        int idx = draw(w);
        if (!kids[0].isLeaf()) {
            labels.put(kids[0], fA.get(idx / fB.size()));
            sampleChildrenEI(kids[0], bits, labels);
        }
        if (!kids[1].isLeaf()) {
            labels.put(kids[1], fB.get(idx % fB.size()));
            sampleChildrenEI(kids[1], bits, labels);
        }
    }

    /** MAP labels given the topology (exact: per node for EL, tree DP for EI). */
    private void assignMapLabels(Tree tree, Map<Node, BitSet> bits, Map<Node, Integer> labels) {
        Node root = tree.getRoot();
        if (context == SAContext.EL) {
            for (Node v : tree.getInternalNodes()) {
                int best = NONE;
                double bestP = -1;
                for (int t : feasibleLabels(v)) {
                    double p = elProb(v, t, bits);
                    if (p > bestP) {
                        bestP = p;
                        best = t;
                    }
                }
                labels.put(v, best);
            }
            return;
        }
        Map<Node, Map<Integer, Double>> best = new HashMap<>();
        int bestT = NONE;
        double bestScore = Double.NEGATIVE_INFINITY;
        for (int t : feasibleLabels(root)) {
            double s = Math.log(eiRootProb(root, t, bits)) + mapScore(root, t, bits, best);
            if (s > bestScore) {
                bestScore = s;
                bestT = t;
            }
        }
        labels.put(root, bestT);
        mapBacktrack(root, bits, best, labels);
    }

    /** Best log-prob of the labels strictly below p, given p's label tp (memoised). */
    private double mapScore(Node p, int tp, Map<Node, BitSet> bits, Map<Node, Map<Integer, Double>> memo) {
        Map<Integer, Double> m = memo.computeIfAbsent(p, x -> new HashMap<>());
        Double cached = m.get(tp);
        if (cached != null) {
            return cached;
        }
        double bestS = Double.NEGATIVE_INFINITY;
        for (int[] pair : childLabelPairs(p, bits)) {
            double s = pairScore(p, tp, pair, bits, memo);
            bestS = Math.max(bestS, s);
        }
        m.put(tp, bestS);
        return bestS;
    }

    private double pairScore(Node p, int tp, int[] pair, Map<Node, BitSet> bits,
                             Map<Node, Map<Integer, Double>> memo) {
        Node[] kids = orderedChildren(p, bits);
        double s = Math.log(eiProb(p, tp, kids, pair[0], pair[1], bits));
        if (!kids[0].isLeaf()) {
            s += mapScore(kids[0], pair[0], bits, memo);
        }
        if (!kids[1].isLeaf()) {
            s += mapScore(kids[1], pair[1], bits, memo);
        }
        return s;
    }

    private void mapBacktrack(Node p, Map<Node, BitSet> bits, Map<Node, Map<Integer, Double>> memo,
                              Map<Node, Integer> labels) {
        int tp = labels.get(p);
        Node[] kids = orderedChildren(p, bits);
        int[] bestPair = null;
        double bestS = Double.NEGATIVE_INFINITY;
        for (int[] pair : childLabelPairs(p, bits)) {
            double s = pairScore(p, tp, pair, bits, memo);
            if (s > bestS) {
                bestS = s;
                bestPair = pair;
            }
        }
        for (int i = 0; i < 2; i++) {
            if (!kids[i].isLeaf()) {
                labels.put(kids[i], bestPair[i]);
                mapBacktrack(kids[i], bits, memo, labels);
            }
        }
    }

    private List<int[]> childLabelPairs(Node p, Map<Node, BitSet> bits) {
        Node[] kids = orderedChildren(p, bits);
        List<Integer> fA = kids[0].isLeaf() ? List.of(NONE) : feasibleLabels(kids[0]);
        List<Integer> fB = kids[1].isLeaf() ? List.of(NONE) : feasibleLabels(kids[1]);
        List<int[]> out = new ArrayList<>();
        for (int a : fA) {
            for (int b : fB) {
                out.add(new int[]{a, b});
            }
        }
        return out;
    }

    /** Label log-probability from an explicit label map (used to stamp sampled trees). */
    private double labelLogProbabilityOf(Tree tree, Map<Node, BitSet> bits, Map<Node, Integer> labels) {
        if (context == SAContext.EL) {
            double logp = 0.0;
            for (Node v : tree.getInternalNodes()) {
                logp += Math.log(elProb(v, labels.get(v), bits));
            }
            return logp;
        }
        Node root = tree.getRoot();
        double logp = Math.log(eiRootProb(root, labels.get(root), bits));
        for (Node p : tree.getInternalNodes()) {
            Node[] kids = orderedChildren(p, bits);
            int tA = kids[0].isLeaf() ? NONE : labels.get(kids[0]);
            int tB = kids[1].isLeaf() ? NONE : labels.get(kids[1]);
            logp += Math.log(eiProb(p, labels.get(p), kids, tA, tB, bits));
        }
        return logp;
    }

    /**
     * Encodes labels as zero-length leaf branches: an SA leaf gets its parent's height. Heights of
     * non-SA leaves are left as the topology sampler set them.
     */
    private static void applyLabels(Node v, Map<Node, Integer> labels) {
        if (v.isLeaf()) {
            return;
        }
        int t = labels.get(v);
        for (Node c : v.getChildren()) {
            if (c.isLeaf() && c.getNr() == t) {
                c.setHeight(v.getHeight());
            }
            applyLabels(c, labels);
        }
    }

    private int draw(double[] w) {
        double total = 0.0;
        for (double x : w) {
            total += x;
        }
        double r = random().nextDouble() * total;
        double cum = 0.0;
        for (int i = 0; i < w.length; i++) {
            cum += w[i];
            if (r < cum) {
                return i;
            }
        }
        return w.length - 1;
    }

    /* -- HELPERS -- */

    /** Whether a leaf is a sampled ancestor (zero-length branch), as in CCD1EL/CCD1EI. */
    private static boolean isSALeaf(Node leaf) {
        return leaf.isLeaf() && !leaf.isRoot() && leaf.getLength() == 0;
    }

    /** Observed label of an internal node: the number of its SA leaf child, or NONE. */
    static int observedLabel(Node v) {
        int t = NONE;
        for (Node c : v.getChildren()) {
            if (isSALeaf(c)) {
                if (t != NONE) {
                    throw new IllegalArgumentException("Node has two zero-length leaf children ("
                            + t + ", " + c.getNr() + "); not a valid sampled-ancestor tree");
                }
                t = c.getNr();
            }
        }
        return t;
    }

    private BitSet computeTaxaBits(Node v, Map<Node, BitSet> bits) {
        BitSet b = BitSet.newBitSet(leafArraySize);
        if (v.isLeaf()) {
            b.set(v.getNr());
        } else {
            b.or(computeTaxaBits(v.getChildren().get(0), bits));
            b.or(computeTaxaBits(v.getChildren().get(1), bits));
        }
        bits.put(v, b);
        return b;
    }

    /** Children in canonical order (by smallest taxon), so keys are invariant to child order. */
    private static Node[] orderedChildren(Node v, Map<Node, BitSet> bits) {
        Node a = v.getChildren().get(0);
        Node b = v.getChildren().get(1);
        return (bits.get(a).nextSetBit(0) < bits.get(b).nextSetBit(0))
                ? new Node[]{a, b} : new Node[]{b, a};
    }

    private static BitSet firstChildBits(Node v, Map<Node, BitSet> bits) {
        return v.isLeaf() ? null : bits.get(orderedChildren(v, bits)[0]);
    }

    private static ElKey elKey(Node v, Map<Node, BitSet> bits) {
        return new ElKey(bits.get(v), firstChildBits(v, bits));
    }

    private static EiKey eiKey(Node p, int tp, Node[] kids, Map<Node, BitSet> bits) {
        return new EiKey(bits.get(p), bits.get(kids[0]), tp,
                firstChildBits(kids[0], bits), firstChildBits(kids[1], bits));
    }

    /** Sentinel parent label distinct from NONE and every taxon, for the root's own label. */
    private static final int ROOT_SENTINEL = Integer.MIN_VALUE;

    private static EiKey rootEiKey(Node root, Map<Node, BitSet> bits) {
        return new EiKey(bits.get(root), firstChildBits(root, bits), ROOT_SENTINEL, null, null);
    }

    private static long pairCode(int tA, int tB) {
        return ((long) tA << 32) ^ (tB & 0xffffffffL);
    }

    @Override
    public String toString() {
        return "SAKRegCCD-" + context + " [beta = " + beta + ", gamma = " + gamma + "] "
                + super.toString();
    }
}
