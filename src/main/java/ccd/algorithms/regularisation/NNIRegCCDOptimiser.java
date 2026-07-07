package ccd.algorithms.regularisation;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator;
import ccd.algorithms.NNICladeExpansion.PairingMode;
import ccd.model.AbstractCCD;
import ccd.model.Clade;
import ccd.model.CladePartition;
import ccd.model.NNIRegCCD;
import ccd.model.bitsets.BitSet;
import ccd.tools.CCDToolUtil;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

/**
 * Selects the two-level smoothing pseudocounts ({@code alpha} for observed and
 * split-expanded splits, {@code beta} for NNI-derived splits) of an
 * {@link NNIRegCCD} by maximising the mean held-out log probability under
 * batched (strided) leave-one-tree-out.
 *
 * <p>
 * The pseudocounts do not change the clade/split graph, only the smoothing, so
 * the expensive per-fold graph build is done once. For each held-out tree we
 * extract a lightweight <em>scoring trace</em> (the raw counts and NNI-flag at
 * each multi-split node it traverses) and discard the graph; the objective at any
 * {@code (alpha, beta)} is then a closed-form sum over those traces:
 * </p>
 * <pre>
 *   log CCP(S) = log( f(S) + w )  -  log( f(C) + alpha*|R(C)| + beta*|N(C)| ),
 *               w = beta if S is NNI-derived else alpha,
 * </pre>
 * with single-split clades contributing a deterministic factor of 1. Coverage is
 * independent of {@code (alpha, beta)} (any positive values keep every split
 * positive), so the objective is a fixed-denominator mean over the covered
 * held-out trees and the optimisation cannot trade coverage for sharpness.
 *
 * @author Claude (CCD-Sophie)
 */
public class NNIRegCCDOptimiser {

    /** The per-node contribution of one multi-split clade to a tree's log probability. */
    public record Factor(int fS, int fC, int regCount, int nniCount, boolean nni) {
        double logCCP(double alpha, double beta) {
            double w = nni ? beta : alpha;
            return Math.log(fS + w) - Math.log(fC + alpha * regCount + beta * nniCount);
        }
    }

    /** Result of an (alpha, beta) optimisation for one pairing mode. */
    public record OptResult(PairingMode mode, double alpha, double beta,
                            double meanLogProb, int covered, int total) {
        @Override
        public String toString() {
            return String.format(
                    "%-12s best alpha=%.4f beta=%.4f   mean held-out logP=%9.4f   covered %d/%d (%.1f%%)",
                    mode, alpha, beta, meanLogProb, covered, total,
                    100.0 * covered / total);
        }
    }

    /* -- SCORING TRACES -- */

    /**
     * Extracts a scoring trace per test tree from a built CCD graph. A trace is
     * the list of multi-split-clade {@link Factor}s the tree traverses; it is
     * {@code null} if the tree is not covered (a clade or split is missing).
     *
     * @param ccd       a built (NNI-)regCCD graph (its CCPs are not used)
     * @param testTrees trees to trace
     * @return one trace (or null) per test tree
     */
    public static List<List<Factor>> extractTraces(AbstractCCD ccd, List<Tree> testTrees) {
        Map<Clade, int[]> weightCache = new HashMap<>();
        List<List<Factor>> traces = new ArrayList<>(testTrees.size());
        for (Tree tree : testTrees) {
            List<Factor> trace = new ArrayList<>();
            boolean[] uncovered = {false};
            extract(tree.getRoot(), ccd, trace, uncovered, weightCache);
            traces.add(uncovered[0] ? null : trace);
        }
        return traces;
    }

    /* Recursive trace extraction; returns the clade of this vertex or null if missing. */
    private static Clade extract(Node vertex, AbstractCCD ccd, List<Factor> trace,
                                 boolean[] uncovered, Map<Clade, int[]> weightCache) {
        if (uncovered[0]) {
            return null;
        }
        if (vertex.isLeaf()) {
            BitSet bits = BitSet.newBitSet(ccd.getSizeOfLeavesArray());
            bits.set(vertex.getNr());
            Clade leaf = ccd.getClade(bits);
            if (leaf == null) {
                uncovered[0] = true;
            }
            return leaf;
        }

        Clade c1 = extract(vertex.getChildren().get(0), ccd, trace, uncovered, weightCache);
        Clade c2 = extract(vertex.getChildren().get(1), ccd, trace, uncovered, weightCache);
        if (uncovered[0] || c1 == null || c2 == null) {
            uncovered[0] = true;
            return null;
        }

        BitSet bits = (BitSet) c1.getCladeInBits().clone();
        bits.or(c2.getCladeInBits());
        Clade clade = ccd.getClade(bits);
        if (clade == null) {
            uncovered[0] = true;
            return null;
        }
        CladePartition partition = clade.getCladePartition(c1, c2);
        if (partition == null) {
            uncovered[0] = true;
            return null;
        }

        // single-split clades are deterministic (CCP = 1), contributing nothing
        if (clade.getNumberOfPartitions() > 1) {
            int[] w = weightCache.computeIfAbsent(clade, NNIRegCCDOptimiser::splitWeights);
            trace.add(new Factor(partition.getNumberOfOccurrences(),
                    clade.getNumberOfOccurrences(), w[0], w[1],
                    CCD1Regularisor.isNNISplit(clade, partition)));
        }
        return clade;
    }

    /* {regular split count, NNI split count} for a clade. */
    private static int[] splitWeights(Clade clade) {
        int reg = 0, nni = 0;
        for (CladePartition partition : clade.getPartitions()) {
            if (CCD1Regularisor.isNNISplit(clade, partition)) {
                nni++;
            } else {
                reg++;
            }
        }
        return new int[]{reg, nni};
    }

    /* -- OBJECTIVE -- */

    /** Log probability of one covered tree's trace at the given pseudocounts. */
    public static double logProb(List<Factor> trace, double alpha, double beta) {
        double lp = 0.0;
        for (Factor factor : trace) {
            lp += factor.logCCP(alpha, beta);
        }
        return lp;
    }

    /** Mean log probability over the covered traces (null traces are skipped). */
    public static double meanLogProb(List<List<Factor>> traces, double alpha, double beta) {
        double sum = 0.0;
        int n = 0;
        for (List<Factor> trace : traces) {
            if (trace != null) {
                sum += logProb(trace, alpha, beta);
                n++;
            }
        }
        return (n == 0) ? Double.NaN : sum / n;
    }

    /* -- OPTIMISATION -- */

    /**
     * Builds the strided CV fold graphs once and extracts a held-out scoring
     * trace per tree. This is the expensive step; the cheap {@code (alpha, beta)}
     * grid search then reuses these traces.
     *
     * @param trees all trees (post burn-in)
     * @param folds number of strided CV folds (B)
     * @param mode  NNI pairing mode
     * @return one trace (or null if uncovered) per tree
     */
    public static List<List<Factor>> buildFoldTraces(List<Tree> trees, int folds, PairingMode mode) {
        int n = trees.size();
        List<List<Factor>> allTraces = new ArrayList<>(n);
        for (int f = 0; f < folds; f++) {
            List<Tree> train = new ArrayList<>(n);
            List<Tree> test = new ArrayList<>();
            for (int i = 0; i < n; i++) {
                if (i % folds == f) {
                    test.add(trees.get(i));
                } else {
                    train.add(trees.get(i));
                }
            }
            NNIRegCCD foldModel = new NNIRegCCD(train, 0.0, mode, NNIRegCCD.DEFAULT_ALPHA);
            allTraces.addAll(extractTraces(foldModel, test));
            System.out.println("  [" + mode + "] fold " + (f + 1) + "/" + folds + " traced");
        }
        return allTraces;
    }

    /** Grid search for the maximum mean held-out log probability over the traces. */
    public static OptResult argmax(List<List<Factor>> traces, PairingMode mode,
                                   double[] alphas, double[] betas) {
        int covered = 0;
        for (List<Factor> trace : traces) {
            if (trace != null) {
                covered++;
            }
        }
        double bestAlpha = alphas[0], bestBeta = betas[0], best = Double.NEGATIVE_INFINITY;
        for (double alpha : alphas) {
            for (double beta : betas) {
                double mean = meanLogProb(traces, alpha, beta);
                if (mean > best) {
                    best = mean;
                    bestAlpha = alpha;
                    bestBeta = beta;
                }
            }
        }
        return new OptResult(mode, bestAlpha, bestBeta, best, covered, traces.size());
    }

    /**
     * Builds the strided CV fold graphs once and grid searches {@code (alpha,
     * beta)} for the maximum mean held-out log probability.
     *
     * @param trees  all trees (post burn-in)
     * @param folds  number of strided CV folds (B)
     * @param mode   NNI pairing mode
     * @param alphas alpha grid
     * @param betas  beta grid
     * @return the best {@code (alpha, beta)} and its mean held-out log probability
     */
    public static OptResult optimise(List<Tree> trees, int folds, PairingMode mode,
                                     double[] alphas, double[] betas) {
        return argmax(buildFoldTraces(trees, folds, mode), mode, alphas, betas);
    }

    private static final double[] DEFAULT_ALPHAS = {0.05, 0.1, 0.2, 0.4, 0.8, 1.6};
    private static final double[] DEFAULT_BETAS = {0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.4, 0.8};

    /* Warn if an optimum sits on the edge of the search grid (true optimum may lie beyond). */
    private static void warnIfOnBoundary(OptResult r, double[] alphas, double[] betas) {
        if (r.alpha() == alphas[0] || r.alpha() == alphas[alphas.length - 1]) {
            System.out.println("  warning: best alpha at grid boundary (" + r.alpha() + "); widen the grid");
        }
        if (r.beta() == betas[0] || r.beta() == betas[betas.length - 1]) {
            System.out.println("  warning: best beta at grid boundary (" + r.beta() + "); widen the grid");
        }
    }

    private static List<Tree> loadTrees(String path, int burninPercentage) throws Exception {
        TreeAnnotator.MemoryFriendlyTreeSet treeSet = CCDToolUtil.getTreeSet(path, burninPercentage);
        List<Tree> trees = new ArrayList<>();
        treeSet.reset();
        while (treeSet.hasNext()) {
            trees.add(treeSet.next());
        }
        return trees;
    }

    private static List<PairingMode> parseModes(String sel) {
        return switch (sel.toLowerCase()) {
            case "co", "co_occurring" -> List.of(PairingMode.CO_OCCURRING);
            case "gw", "graph_wide" -> List.of(PairingMode.GRAPH_WIDE);
            default -> List.of(PairingMode.values());
        };
    }

    /* batch <out.csv> <folds> <burnin%> <co|gw|both> <trees...> : optimise each dataset, append a CSV row. */
    private static void batchMain(String[] args) throws Exception {
        String out = args[1];
        int folds = Integer.parseInt(args[2]);
        int burnin = Integer.parseInt(args[3]);
        List<PairingMode> modes = parseModes(args[4]);

        File csv = new File(out);
        boolean writeHeader = !csv.exists() || csv.length() == 0;
        try (PrintWriter w = new PrintWriter(new FileWriter(csv, true))) {
            if (writeHeader) {
                w.println("dataset,taxa,trees,mode,alpha,beta,meanLogProb,covered,total");
                w.flush();
            }
            for (int a = 5; a < args.length; a++) {
                String path = args[a];
                String name = new File(path).getName().replaceAll("\\.trees$", "");
                List<Tree> trees;
                try {
                    trees = loadTrees(path, burnin);
                } catch (Exception e) {
                    System.out.println("SKIP " + name + " (" + e.getMessage() + ")");
                    continue;
                }
                int taxa = trees.get(0).getLeafNodeCount();
                System.out.println("=== " + name + " : " + taxa + " taxa, " + trees.size() + " trees ===");
                for (PairingMode mode : modes) {
                    OptResult r = optimise(trees, folds, mode, DEFAULT_ALPHAS, DEFAULT_BETAS);
                    w.printf(Locale.US, "%s,%d,%d,%s,%.4f,%.4f,%.6f,%d,%d%n",
                            name, taxa, trees.size(), mode, r.alpha(), r.beta(),
                            r.meanLogProb(), r.covered(), r.total());
                    w.flush();
                    System.out.println(r);
                    warnIfOnBoundary(r, DEFAULT_ALPHAS, DEFAULT_BETAS);
                }
            }
        }
        System.out.println("Wrote " + out);
    }

    /* Parse "start:end:step" into an inclusive grid (e.g. "0.15:0.25:0.01" → 11 values). */
    private static double[] parseGrid(String spec) {
        String[] parts = spec.split(":");
        if (parts.length != 3) {
            throw new IllegalArgumentException("grid spec must be start:end:step, got: " + spec);
        }
        double start = Double.parseDouble(parts[0]);
        double end = Double.parseDouble(parts[1]);
        double step = Double.parseDouble(parts[2]);
        if (step <= 0 || end < start) {
            throw new IllegalArgumentException("invalid grid: " + spec);
        }
        int n = (int) Math.round((end - start) / step) + 1;
        double[] grid = new double[n];
        for (int i = 0; i < n; i++) {
            grid[i] = start + i * step;
        }
        return grid;
    }

    /* heldout <train> <test> <burnin%> <co|gw|both> <alpha-spec> <beta-spec> : fit on train, score on test. */
    private static void heldOutMain(String[] args) throws Exception {
        if (args.length < 7) {
            System.err.println("Usage: NNIRegCCDOptimiser heldout <train.trees> <test.trees> "
                    + "<burnin%> <co|gw|both> <alpha-start:end:step> <beta-start:end:step>");
            System.exit(1);
        }
        String trainPath = args[1];
        String testPath = args[2];
        int burnin = Integer.parseInt(args[3]);
        List<PairingMode> modes = parseModes(args[4]);
        double[] alphas = parseGrid(args[5]);
        double[] betas = parseGrid(args[6]);

        List<Tree> train = loadTrees(trainPath, burnin);
        List<Tree> test = loadTrees(testPath, burnin);
        System.out.printf("train: %d trees (%s), test: %d trees (%s)%n",
                train.size(), new File(trainPath).getName(),
                test.size(), new File(testPath).getName());
        System.out.printf("grid: %d alphas [%g .. %g] × %d betas [%g .. %g] = %d points%n",
                alphas.length, alphas[0], alphas[alphas.length - 1],
                betas.length, betas[0], betas[betas.length - 1],
                alphas.length * betas.length);

        for (PairingMode mode : modes) {
            System.out.println("building NNIRegCCD[" + mode + "] on training set...");
            NNIRegCCD model = new NNIRegCCD(train, 0.0, mode, NNIRegCCD.DEFAULT_ALPHA);
            List<List<Factor>> traces = extractTraces(model, test);
            OptResult r = argmax(traces, mode, alphas, betas);
            System.out.println(r);
            warnIfOnBoundary(r, alphas, betas);
        }
    }

    /* surface <out.csv> <trees> <folds> <burnin%> <co|gw> : dump the full (alpha,beta) objective grid. */
    private static void surfaceMain(String[] args) throws Exception {
        String out = args[1];
        int folds = Integer.parseInt(args[3]);
        int burnin = Integer.parseInt(args[4]);
        PairingMode mode = parseModes(args[5]).get(0);
        List<Tree> trees = loadTrees(args[2], burnin);
        List<List<Factor>> traces = buildFoldTraces(trees, folds, mode);
        try (PrintWriter w = new PrintWriter(new FileWriter(out))) {
            w.println("alpha,beta,meanLogProb");
            for (double alpha : DEFAULT_ALPHAS) {
                for (double beta : DEFAULT_BETAS) {
                    w.printf(Locale.US, "%.4f,%.4f,%.6f%n", alpha, beta, meanLogProb(traces, alpha, beta));
                }
            }
        }
        System.out.println("Wrote surface to " + out);
    }

    /**
     * CLI. Forms:
     * <ul>
     * <li>{@code <trees> [folds] [burnin%]} --- optimise (alpha, beta) for each
     * pairing mode on one dataset;</li>
     * <li>{@code batch <out.csv> <folds> <burnin%> <co|gw|both> <trees...>} ---
     * optimise each dataset and append a row per (dataset, mode) to a CSV;</li>
     * <li>{@code surface <out.csv> <trees> <folds> <burnin%> <co|gw>} --- dump the
     * full (alpha, beta) objective grid for one dataset.</li>
     * </ul>
     *
     * @param args see above
     * @throws Exception if a tree file cannot be read
     */
    public static void main(String[] args) throws Exception {
        if (args.length >= 1 && args[0].equals("batch")) {
            batchMain(args);
            return;
        }
        if (args.length >= 1 && args[0].equals("surface")) {
            surfaceMain(args);
            return;
        }
        if (args.length >= 1 && args[0].equals("heldout")) {
            heldOutMain(args);
            return;
        }
        if (args.length < 1) {
            System.err.println("Usage: NNIRegCCDOptimiser <trees> [folds] [burnin%]");
            System.err.println("   or: NNIRegCCDOptimiser batch <out.csv> <folds> <burnin%> <co|gw|both> <trees...>");
            System.err.println("   or: NNIRegCCDOptimiser surface <out.csv> <trees> <folds> <burnin%> <co|gw>");
            System.err.println("   or: NNIRegCCDOptimiser heldout <train> <test> <burnin%> <co|gw|both> <alpha-spec> <beta-spec>");
            System.exit(1);
        }
        int folds = (args.length >= 2) ? Integer.parseInt(args[1]) : 20;
        int burnin = (args.length >= 3) ? Integer.parseInt(args[2]) : 10;
        List<Tree> trees = loadTrees(args[0], burnin);
        System.out.println("(alpha,beta) optimisation: " + folds + "-fold strided CV, "
                + trees.size() + " trees");
        for (PairingMode mode : PairingMode.values()) {
            OptResult result = optimise(trees, folds, mode, DEFAULT_ALPHAS, DEFAULT_BETAS);
            System.out.println(result);
            warnIfOnBoundary(result, DEFAULT_ALPHAS, DEFAULT_BETAS);
        }
    }
}
