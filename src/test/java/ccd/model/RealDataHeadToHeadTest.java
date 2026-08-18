package ccd.model;

import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator;
import ccd.tools.CCDToolUtil;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * Head-to-head held-out predictive comparison of the regularised CCD variants on a real posterior
 * tree sample.
 *
 * <p>Protocol follows the manuscript's RSV2 comparison: the model is trained on trees drawn from the
 * first half of the chain and scored on trees from the second half, so the test trees are genuinely
 * out of training. Each model's hyperparameters are selected on an inner fit/validation split of the
 * training half alone, then the model is rebuilt on the whole training half and scored on the test
 * half. Reported per model: support coverage (fraction of test trees with positive probability), mean
 * log probability over all test trees, and -- for the full-support models -- the paired per-tree
 * comparison against KRegCCD.
 *
 * <p>Point the test at a tree file with {@code -Dccd.trees=/path/to/x.trees}; it is skipped when no
 * file is given. {@code -Dccd.n=1000} sets the number of training and test trees.
 */
public class RealDataHeadToHeadTest {

    private static final String PATH = System.getProperty("ccd.trees", "");
    private static final int N = Integer.parseInt(System.getProperty("ccd.n", "1000"));
    private static final double BURNIN_PERCENT = 10;

    private interface Scorer {
        double logP(Tree t);
    }

    private static List<String> newickCache;
    private static List<String> taxaCache;

    /** Parses the tree file once, keeping the topologies as newick strings. */
    private static void load() throws Exception {
        if (newickCache != null) {
            return;
        }
        TreeAnnotator.TreeSet ts = CCDToolUtil.getTreeSet(PATH, (int) BURNIN_PERCENT);
        ts.reset();
        List<String> nwk = new ArrayList<>();
        List<String> taxa = null;
        while (ts.hasNext()) {
            Tree t = ts.next();
            if (taxa == null) {
                taxa = new ArrayList<>();
                String[] byNr = new String[t.getLeafNodeCount()];
                for (beast.base.evolution.tree.Node leaf : t.getExternalNodes()) {
                    byNr[leaf.getNr()] = leaf.getID();
                }
                taxa.addAll(List.of(byNr));
            }
            nwk.add(t.getRoot().toNewick() + ";");
        }
        newickCache = nwk;
        taxaCache = taxa;
        System.out.printf("loaded %d trees, %d taxa from %s%n",
                nwk.size(), taxa.size(), new File(PATH).getName());
    }

    /**
     * {@code count} evenly spaced trees from the chain segment {@code [from, to)} (as fractions of the
     * post-burn-in chain), freshly parsed on every call because the CCD constructors take ownership of
     * the trees they are given.
     *
     * <p>Segments must be disjoint for hyperparameter selection to be honest: scoring trees that are
     * also in the fitted set drives every regularisation parameter to zero, because the backbone
     * already fits its own training trees and any reserved mass is then pure loss.
     */
    private static List<Tree> read(int count, double from, double to) throws Exception {
        load();
        int lo = (int) (from * newickCache.size());
        int hi = (int) (to * newickCache.size());
        List<String> pool = newickCache.subList(lo, hi);
        List<Tree> out = new ArrayList<>();
        double step = Math.max(1.0, pool.size() / (double) count);
        for (int i = 0; i < count && (int) (i * step) < pool.size(); i++) {
            out.add(new beast.base.evolution.tree.TreeParser(taxaCache, pool.get((int) (i * step)), 0, false));
        }
        return out;
    }

    /** Trees for the final models and for scoring: first half of the chain trains, second half tests. */
    private static List<Tree> train(int count) throws Exception {
        return read(count, 0.0, 0.5);
    }

    /** Disjoint inner split of the training half: [0, 0.25) fits, [0.25, 0.5) validates. */
    private static List<Tree> fitSet(int count) throws Exception {
        return read(count, 0.0, 0.25);
    }

    private static List<Tree> valSet(int count) throws Exception {
        return read(count, 0.25, 0.5);
    }

    private static double[] score(Scorer s, List<Tree> test) {
        int covered = 0;
        double sum = 0.0;
        for (Tree t : test) {
            double lp = s.logP(t);
            if (Double.isFinite(lp)) {
                covered++;
                sum += lp;
            }
        }
        return new double[]{covered, covered == 0 ? Double.NEGATIVE_INFINITY : sum / covered};
    }

    /** Mean log probability over every test tree, treating unsupported trees as -infinity. */
    private static double meanAll(Scorer s, List<Tree> test) {
        double sum = 0.0;
        for (Tree t : test) {
            double lp = s.logP(t);
            if (!Double.isFinite(lp)) {
                return Double.NEGATIVE_INFINITY;
            }
            sum += lp;
        }
        return sum / test.size();
    }

    @Test
    public void headToHeadOnRealData() throws Exception {
        Assumptions.assumeTrue(!PATH.isEmpty() && new File(PATH).exists(),
                "set -Dccd.trees=/path/to/x.trees to run this comparison");

        List<Tree> probe = train(N);
        int nTaxa = probe.get(0).getLeafNodeCount();
        System.out.printf("%n=== %s: %d taxa, %d train / %d test trees ===%n",
                new File(PATH).getName(), nTaxa, probe.size(), N);

        List<Tree> test = read(N, 0.5, 1.0);
        List<Tree> val = valSet(N / 2);
        int nFit = N / 2;

        // ---- CCD1 ----
        CCD1 ccd1 = new CCD1(train(N), 0.0);
        report("CCD1", "-", ccd1::getLogProbabilityOfTree, test);

        // ---- RegCCD: alpha on validation ----
        double bestAlpha = 0.4;
        double bestAlphaScore = Double.NEGATIVE_INFINITY;
        for (double alpha : new double[]{0.01, 0.05, 0.1, 0.2, 0.4, 0.8, 1.0}) {
            RegCCD m = new RegCCD(fitSet(nFit), 0.0, alpha);
            double sc = score(m::getLogProbabilityOfTree, val)[1];
            if (sc > bestAlphaScore) {
                bestAlphaScore = sc;
                bestAlpha = alpha;
            }
        }
        RegCCD reg = new RegCCD(train(N), 0.0, bestAlpha);
        report("RegCCD", String.format("alpha=%.2f", bestAlpha), reg::getLogProbabilityOfTree, test);

        // ---- KRegCCD: mu on validation at alpha = 0.4 (the manuscript's setting) ----
        KRegCCD kFit = new KRegCCD(fitSet(nFit), 0.0, 0.005, 0.4);
        double bestMu = 0.005;
        double bestMuScore = Double.NEGATIVE_INFINITY;
        for (double mu : new double[]{0.00002, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05}) {
            final double m = mu;
            double sc = score(t -> kFit.getLogProbabilityOfTree(t, m), val)[1];
            if (sc > bestMuScore) {
                bestMuScore = sc;
                bestMu = mu;
            }
        }
        KRegCCD kreg = new KRegCCD(train(N), 0.0, bestMu, 0.4);
        report("KRegCCD", String.format("alpha=0.4, mu=%.5f", bestMu),
                kreg::getLogProbabilityOfTree, test);

        // ---- MRegCCD: mu on validation ----
        MRegCCD mFit = new MRegCCD(fitSet(nFit), 0.0, MRegCCD.DEFAULT_MU);
        double bestMMu = MRegCCD.DEFAULT_MU;
        double bestMMuScore = Double.NEGATIVE_INFINITY;
        for (double mu : new double[]{0.00005, 0.0002, 0.001, 0.002, 0.008, 0.0159, 0.05, 0.1}) {
            final double m = mu;
            double sc = score(t -> mFit.getLogProbabilityOfTree(t, m), val)[1];
            if (sc > bestMMuScore) {
                bestMMuScore = sc;
                bestMMu = mu;
            }
        }
        MRegCCD mreg = new MRegCCD(train(N), 0.0, bestMMu);
        report("MRegCCD", String.format("mu=%.5f", bestMMu), mreg::getLogProbabilityOfTree, test);

        // ---- CRegCCD: (a2, a3, a4) on validation, a1 = 0 ----
        CRegCCD cFit = new CRegCCD(fitSet(nFit), 0.0);
        double[] grid = {0.002, 0.01, 0.05, 0.2, 0.4, 1.0, 2.0, 5.0, 12.0};
        double[] bestC = {0.0, 0.4, 0.4, 0.4};
        double bestCScore = Double.NEGATIVE_INFINITY;
        for (double b2 : grid) {
            for (double b3 : grid) {
                for (double b4 : grid) {
                    double sc = score(t -> cFit.getLogProbabilityOfTree(t, b2, b3, b4), val)[1];
                    if (sc > bestCScore) {
                        bestCScore = sc;
                        bestC = new double[]{0.0, b2, b3, b4};
                    }
                }
            }
        }
        CRegCCD creg = new CRegCCD(train(N), 0.0, bestC[1], bestC[2], bestC[3]);
        report("CRegCCD", String.format("alpha=%.3f, alpha1=%.3f, alpha2=%.3f", bestC[1], bestC[2], bestC[3]),
                creg::getLogProbabilityOfTree, test);

        // ---- paired comparison of the full-support models against KRegCCD ----
        System.out.printf("%npaired per-tree comparison against KRegCCD (n = %d test trees):%n", test.size());
        paired("MRegCCD", mreg::getLogProbabilityOfTree, kreg::getLogProbabilityOfTree, test);
        paired("CRegCCD", creg::getLogProbabilityOfTree, kreg::getLogProbabilityOfTree, test);
    }

    private static void report(String name, String params, Scorer s, List<Tree> test) {
        double[] sc = score(s, test);
        double all = meanAll(s, test);
        System.out.printf("%-9s %-27s coverage %6.1f%%   mean logP(covered) %10.2f   mean logP(all) %s%n",
                name, params, 100.0 * sc[0] / test.size(), sc[1],
                Double.isFinite(all) ? String.format("%10.2f", all) : "      -inf");
    }

    private static void paired(String name, Scorer a, Scorer baseline, List<Tree> test) {
        int wins = 0;
        double sumDiff = 0.0;
        double sumSq = 0.0;
        for (Tree t : test) {
            double d = a.logP(t) - baseline.logP(t);
            if (d > 0) {
                wins++;
            }
            sumDiff += d;
            sumSq += d * d;
        }
        int n = test.size();
        double mean = sumDiff / n;
        double se = Math.sqrt(Math.max(0, sumSq / n - mean * mean) / n);
        System.out.printf("  %-9s mean log-ratio %+8.2f +/- %.2f nats/tree, better on %d/%d trees%n",
                name, mean, se, wins, n);
    }

    /** Scale check: sampling from a real 129-taxon posterior must be fast, valid and
     *  self-consistent with the scorer. */
    @Test
    public void samplingOnRealData() throws Exception {
        Assumptions.assumeTrue(!PATH.isEmpty() && new File(PATH).exists(),
                "set -Dccd.trees=/path/to/x.trees to run this comparison");
        CRegCCD ccd = new CRegCCD(train(N), 0.0);
        int nTaxa = ccd.getSizeOfLeavesArray();
        long t0 = System.nanoTime();
        int draws = 200;
        for (int i = 0; i < draws; i++) {
            Tree t = ccd.sampleTree();
            if (t.getLeafNodeCount() != nTaxa || t.getNodeCount() != 2 * nTaxa - 1) {
                throw new AssertionError("invalid sampled tree");
            }
            double stamped = (Double) t.getRoot().getMetaData(CCD1.LOG_PROB_SUBTREE_KEY);
            double scored = ccd.getLogProbabilityOfTree(t);
            if (Math.abs(stamped - scored) > 1e-9) {
                throw new AssertionError("stamped " + stamped + " != scored " + scored);
            }
        }
        double secs = (System.nanoTime() - t0) / 1e9;
        System.out.printf("CRegCCD sampling on %s: %d taxa, %d trees in %.2f s (%.1f ms/tree), "
                + "all valid and self-consistent%n",
                new File(PATH).getName(), nTaxa, draws, secs, 1000 * secs / draws);
    }

    /** MAP and entropy on a real posterior: correctness certificate and wall-clock cost. */
    @Test
    public void mapAndEntropyOnRealData() throws Exception {
        Assumptions.assumeTrue(!PATH.isEmpty() && new File(PATH).exists(),
                "set -Dccd.trees=/path/to/x.trees to run this comparison");
        CRegCCD creg = new CRegCCD(train(N), 0.0, 2.0, 0.4, 0.05);

        long t0 = System.nanoTime();
        double maxLog = creg.getMaxLogTreeProbability();
        boolean certified = creg.isMAPCertifiedGlobal();
        double bound = creg.getOffBackboneBound();
        double mapSecs = (System.nanoTime() - t0) / 1e9;

        t0 = System.nanoTime();
        Tree map = creg.getMAPTree();
        double treeSecs = (System.nanoTime() - t0) / 1e9;
        double scored = creg.getLogProbabilityOfTree(map);

        t0 = System.nanoTime();
        double[] h = creg.getEntropyMonteCarlo(20_000);
        double entSecs = (System.nanoTime() - t0) / 1e9;

        System.out.printf("%n=== %s: CRegCCD MAP and entropy ===%n", new File(PATH).getName());
        System.out.printf("MAP DP + certificate : %.2f s, max logP = %.4f, "
                        + "off-backbone bound = %.4f, certified global = %s%n",
                mapSecs, maxLog, bound, certified);
        System.out.printf("MAP tree traceback   : %.2f s, scored logP = %.4f (matches: %s)%n",
                treeSecs, scored, Math.abs(scored - maxLog) < 1e-9);
        System.out.printf("entropy (20k draws)  : %.2f s, H = %.3f +/- %.3f nats%n",
                entSecs, h[0], h[1]);
    }

    /** Deterministic entropy recursion vs the unbiased Monte-Carlo estimator on a real posterior. */
    @Test
    public void entropyRecursionVersusMonteCarloOnRealData() throws Exception {
        Assumptions.assumeTrue(!PATH.isEmpty() && new File(PATH).exists(),
                "set -Dccd.trees=/path/to/x.trees to run this comparison");
        System.out.printf("%n=== %s: CRegCCD entropy, recursion vs Monte Carlo ===%n",
                new File(PATH).getName());
        System.out.printf("%-22s %-11s %-9s %-22s %-9s %-10s%n",
                "pseudocounts", "recursion", "rec (s)", "Monte Carlo", "MC (s)", "difference");
        for (double[] p : new double[][]{{2.0, 0.4, 0.05}, {5.0, 2.0, 0.4}, {0.4, 0.4, 0.4}}) {
            CRegCCD ccd = new CRegCCD(train(N), 0.0, p[0], p[1], p[2]);
            long t0 = System.nanoTime();
            double rec = ccd.getEntropyRecursive();
            double recSecs = (System.nanoTime() - t0) / 1e9;
            t0 = System.nanoTime();
            double[] mc = ccd.getEntropyMonteCarlo(200_000);
            double mcSecs = (System.nanoTime() - t0) / 1e9;
            double diff = rec - mc[0];
            System.out.printf("a=(%.2f,%.2f,%.2f)%-6s %-11.4f %-9.2f %8.4f +/-%.4f %-9.2f %+.4f (%+.3f%%)%n",
                    p[0], p[1], p[2], "", rec, recSecs, mc[0], mc[1], mcSecs, diff, 100 * diff / mc[0]);
        }
    }

    /** Does the exact A_0/A_1 MAP search stay tractable on a real posterior? */
    @Test
    public void exactMapOnRealData() throws Exception {
        Assumptions.assumeTrue(!PATH.isEmpty() && new File(PATH).exists(),
                "set -Dccd.trees=/path/to/x.trees to run this comparison");
        CRegCCD creg = new CRegCCD(train(N), 0.0, 0.4, 0.4, 0.05);
        double backbone = creg.getMaxLogTreeProbability();
        System.out.printf("%n=== %s: MAP search by allowed A_1 depth ===%n", new File(PATH).getName());
        System.out.printf("backbone (A_0 only) max logP = %.4f%n", backbone);
        System.out.printf("%-6s %-12s %-12s %-12s %-10s %-8s%n",
                "maxA1", "max logP", "improvement", "bound", "certified", "states");
        for (int k = 0; k <= 3; k++) {
            long t0 = System.nanoTime();
            CRegCCD.MapResult r = creg.solveMAP(k);
            double secs = (System.nanoTime() - t0) / 1e9;
            if (!r.complete()) {
                System.out.printf("%-6d exceeded the state budget after %d states (%.1f s)%n",
                        k, r.statesExplored(), secs);
                break;
            }
            System.out.printf("%-6d %-12.4f %-12.4f %-12.4f %-10s %-8d (%.1f s)%n",
                    k, r.maxLogProbability(), r.maxLogProbability() - backbone,
                    r.offBackboneBound(), r.a2Excluded(), r.statesExplored(), secs);
        }
    }
}
