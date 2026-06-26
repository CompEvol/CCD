package test.ccd.model;

import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import ccd.model.KRegCCD;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Compares the two KRegCCD entropy estimators -- the Monte-Carlo plug-in
 * ({@link KRegCCD#getEntropyMonteCarlo(int)}) and the generalised Lewis recursion
 * ({@link KRegCCD#getEntropyRecursive()}) -- against the <em>exact</em> entropy obtained by brute
 * force over every rooted topology.
 *
 * <p>On four taxa with {@link KRegCCD.TailMode#NONE} the full-support model normalises to 1
 * exactly (every blue-region boundary part is a leaf or cherry, so nothing reserves), which makes
 * brute force a true gold standard: the recursion must match it to numerical precision and the
 * Monte-Carlo estimate to within its standard error. This simultaneously validates both estimators
 * and the score/sampler consistency they each rely on.
 */
public class KRegCCDEntropyTest {

    /* --------------------------------------------------------------------- *
     * four-taxon fixtures (model normalises exactly -> exact entropy available)
     * --------------------------------------------------------------------- */

    private static final List<String> TAXA4 = Arrays.asList("A", "B", "C", "D");

    /** Observed clades AB, CD, ABC -- leaves room for novel clades, so escape entropy is nonzero. */
    private static List<Tree> training4() {
        List<Tree> trees = new ArrayList<>();
        trees.add(new TreeParser(TAXA4, "(((A:1,B:1):1,C:1):1,D:1):0;", 1, false));
        trees.add(new TreeParser(TAXA4, "((A:1,B:1):1,(C:1,D:1):1):0;", 1, false));
        return trees;
    }

    /* --------------------------------------------------------------------- *
     * five-taxon fixtures (boundary parts may reserve -> exact via brute force still,
     * but the sampler/score split P != Q becomes observable)
     * --------------------------------------------------------------------- */

    private static final List<String> TAXA5 = Arrays.asList("A", "B", "C", "D", "E");

    private static List<Tree> training5() {
        List<Tree> trees = new ArrayList<>();
        for (String nwk : new String[]{
                "((((A:1,B:1):1,C:1):1,D:1):1,E:1):0;",
                "((((A:1,B:1):1,C:1):1,D:1):1,E:1):0;",
                "(((A:1,B:1):1,(C:1,D:1):1):1,E:1):0;",
                "(((A:1,B:1):1,(C:1,D:1):1):1,E:1):0;",
                "((((A:1,C:1):1,B:1):1,D:1):1,E:1):0;",
                "(((A:1,B:1):1,C:1):1,(D:1,E:1):1):0;",
        }) {
            trees.add(new TreeParser(TAXA5, nwk, 1, false));
        }
        return trees;
    }

    /* --------------------------------------------------------------------- */

    /** Exact entropy (and total mass) by brute force over all topologies: {@code -sum P log P}. */
    private static double[] exactEntropyAndMass(KRegCCD ccd, List<Tree> topologies) {
        double mass = 0.0;
        double entropy = 0.0;
        for (Tree t : topologies) {
            double logp = ccd.getLogProbabilityOfTree(t);
            double p = Math.exp(logp);
            mass += p;
            entropy -= p * logp;
        }
        return new double[]{entropy, mass};
    }

    /* --------------------------------------------------------------------- */

    @Test
    public void testRecursionMatchesExactOnFourTaxa() {
        for (KRegCCD.NovelMode mode : KRegCCD.NovelMode.values()) {
            KRegCCD ccd = new KRegCCD(training4(), 0.0, 0.05, 0.4, 2, KRegCCD.TailMode.NONE, mode);
            List<Tree> topologies = allRootedTopologies(TAXA4);

            double[] exact = exactEntropyAndMass(ccd, topologies);
            assertEquals(1.0, exact[1], 1e-9,
                    "model must normalise on four taxa [" + mode + "]");

            double recursive = ccd.getEntropyRecursive();
            assertEquals(exact[0], recursive, 1e-9,
                    "recursion must equal the exact entropy [" + mode + "]");
        }
    }

    @Test
    public void testMonteCarloMatchesExactOnFourTaxa() {
        for (KRegCCD.NovelMode mode : KRegCCD.NovelMode.values()) {
            KRegCCD ccd = new KRegCCD(training4(), 0.0, 0.05, 0.4, 2, KRegCCD.TailMode.NONE, mode);
            ccd.setRandom(new Random(20260614L));
            List<Tree> topologies = allRootedTopologies(TAXA4);

            double exact = exactEntropyAndMass(ccd, topologies)[0];
            double[] mc = ccd.getEntropyMonteCarlo(300_000);

            assertEquals(exact, mc[0], 5 * mc[1] + 1e-3,
                    "Monte-Carlo estimate " + mc[0] + " +/- " + mc[1]
                            + " must bracket the exact entropy " + exact + " [" + mode + "]");
        }
    }

    @Test
    public void testRecursionAndMonteCarloAgreeOnFiveTaxa() {
        // On five taxa some boundary parts reserve, so the sampler's Q and the scored P diverge by
        // O(mu) per reserving boundary part: the plug-in -E_Q[log P] sits slightly above the
        // recursion's H(Q). They should still agree to within a few standard errors plus that small
        // bias at the operating mu.
        KRegCCD ccd = new KRegCCD(training5(), 0.0, 0.02, 0.4, 2, KRegCCD.TailMode.NONE);
        ccd.setRandom(new Random(7L));

        double recursive = ccd.getEntropyRecursive();
        double[] mc = ccd.getEntropyMonteCarlo(300_000);

        assertTrue(Double.isFinite(recursive) && recursive > 0, "recursion must give a positive entropy");
        assertEquals(recursive, mc[0], 5 * mc[1] + 0.02,
                "recursion " + recursive + " and Monte-Carlo " + mc[0] + " +/- " + mc[1]
                        + " should agree up to the O(mu) sampler reservation");
    }

    @Test
    public void testGetEntropyFollowsSamplingFidelity() {
        // SELF_CONSISTENT: getEntropy() is exactly the deterministic recursion.
        KRegCCD sc = new KRegCCD(training4(), 0.0, 0.05, 0.4, 2, KRegCCD.TailMode.NONE);
        sc.setSamplingFidelity(KRegCCD.SamplingFidelity.SELF_CONSISTENT);
        assertEquals(sc.getEntropyRecursive(), sc.getEntropy(), 1e-12,
                "SELF_CONSISTENT getEntropy() must equal the exact recursion");

        // FULL_SUPPORT: getEntropy() uses the Monte-Carlo plug-in. On four taxa there is no reserve
        // tail, so the full-support distribution equals the self-consistent one and getEntropy()
        // must match the exact entropy within Monte-Carlo error.
        double exact = exactEntropyAndMass(
                new KRegCCD(training4(), 0.0, 0.05, 0.4, 2, KRegCCD.TailMode.NONE),
                allRootedTopologies(TAXA4))[0];
        KRegCCD fs = new KRegCCD(training4(), 0.0, 0.05, 0.4, 2, KRegCCD.TailMode.SAMPLED);
        fs.setSamplingFidelity(KRegCCD.SamplingFidelity.FULL_SUPPORT);
        fs.setRandom(new Random(3L));
        assertEquals(exact, fs.getEntropy(), 0.03,
                "FULL_SUPPORT getEntropy() must match the exact entropy on four taxa (no tail)");
    }

    @Test
    public void testLewisEntropyRejected() {
        KRegCCD ccd = new KRegCCD(training4(), 0.0, 0.05, 0.4, 2, KRegCCD.TailMode.NONE);
        assertThrows(UnsupportedOperationException.class, ccd::getEntropyLewis,
                "the observed-only Lewis entropy must be rejected for a full-support model");
    }

    /* --------------------------------------------------------------------- *
     * topology enumeration (all rooted binary trees over the taxa)
     * --------------------------------------------------------------------- */

    private static List<Tree> allRootedTopologies(List<String> taxa) {
        List<Tree> trees = new ArrayList<>();
        for (String shape : shapes(taxa)) {
            trees.add(new TreeParser(taxa, shape + ";", 1, false));
        }
        return trees;
    }

    /** Newick fragments of every rooted binary tree over {@code taxa}; the first taxon is forced
     * into the left side so each unordered split is counted once. */
    private static List<String> shapes(List<String> taxa) {
        List<String> out = new ArrayList<>();
        if (taxa.size() == 1) {
            out.add(taxa.get(0) + ":1");
            return out;
        }
        String first = taxa.get(0);
        List<String> rest = taxa.subList(1, taxa.size());
        int n = rest.size();
        for (int mask = 0; mask < (1 << n); mask++) {
            List<String> left = new ArrayList<>();
            left.add(first);
            List<String> right = new ArrayList<>();
            for (int i = 0; i < n; i++) {
                if ((mask & (1 << i)) != 0) {
                    left.add(rest.get(i));
                } else {
                    right.add(rest.get(i));
                }
            }
            if (right.isEmpty()) {
                continue;
            }
            for (String l : shapes(left)) {
                for (String r : shapes(right)) {
                    out.add("(" + l + "," + r + "):1");
                }
            }
        }
        return out;
    }

    /* --------------------------------------------------------------------- *
     * printed comparison (run directly: not part of the JUnit assertions)
     * --------------------------------------------------------------------- */

    public static void main(String[] args) {
        System.out.printf("%-8s %-7s %-6s %-14s %-22s %-12s %-10s%n",
                "taxa", "novel", "mu", "exact H", "Monte-Carlo H (+/-SE)", "recursion H", "totalMass");
        compareRow(TAXA4, training4(), 0.05, KRegCCD.NovelMode.FLAT);
        compareRow(TAXA4, training4(), 0.05, KRegCCD.NovelMode.SHARED);
        compareRow(TAXA4, training4(), 0.01, KRegCCD.NovelMode.FLAT);
        compareRow(TAXA5, training5(), 0.02, KRegCCD.NovelMode.FLAT);
        compareRow(TAXA5, training5(), 0.05, KRegCCD.NovelMode.FLAT);
    }

    private static void compareRow(List<String> taxa, List<Tree> training, double mu,
                                   KRegCCD.NovelMode mode) {
        KRegCCD ccd = new KRegCCD(training, 0.0, mu, 0.4, 2, KRegCCD.TailMode.NONE, mode);
        ccd.setRandom(new Random(1234L));
        double[] exact = exactEntropyAndMass(ccd, allRootedTopologies(taxa));
        double[] mc = ccd.getEntropyMonteCarlo(300_000);
        double rec = ccd.getEntropyRecursive();
        System.out.printf("%-8d %-7s %-6.2f %-14.6f %10.6f +/-%8.6f   %-12.6f %-10.6f%n",
                taxa.size(), mode, mu, exact[0], mc[0], mc[1], rec, exact[1]);
    }
}
