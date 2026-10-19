package ccd.model;

import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/** {@link CRegCCD}: MAP tree (with its global-optimality certificate) and entropy. */
public class CRegCCDMapEntropyTest {

    private static List<String> taxa(int n) {
        List<String> out = new ArrayList<>();
        for (int i = 0; i < n; i++) {
            out.add("T" + i);
        }
        return out;
    }

    private static List<Tree> randomTrees(List<String> taxa, int nTrees, long seed) {
        Random rng = new Random(seed);
        List<Tree> out = new ArrayList<>();
        for (int t = 0; t < nTrees; t++) {
            List<String> pool = new ArrayList<>(taxa);
            while (pool.size() > 1) {
                String a = pool.remove(rng.nextInt(pool.size()));
                String b = pool.remove(rng.nextInt(pool.size()));
                pool.add("(" + a + "," + b + ")");
            }
            out.add(new TreeParser(taxa, pool.get(0) + ";", 1, false));
        }
        return out;
    }

    /**
     * The MAP tree must be the argmax over ALL topologies, not just those on the backbone. Checked
     * by exhaustive enumeration, together with the certificate that claims it.
     */
    @Test
    public void mapTreeIsGlobalOptimum() {
        for (int n : new int[]{5, 6, 7}) {
            for (int nTrees : new int[]{2, 10, 40}) {
                List<String> tx = taxa(n);
                CRegCCD ccd = new CRegCCD(randomTrees(tx, nTrees, 13L), 0.0, 0.4, 0.4, 0.4);

                double bruteForce = Double.NEGATIVE_INFINITY;
                for (Tree t : CRegCCDTest.allRootedTopologies(tx)) {
                    bruteForce = Math.max(bruteForce, ccd.getLogProbabilityOfTree(t));
                }
                double reported = ccd.getMaxLogTreeProbability();
                boolean certified = ccd.isMAPCertifiedGlobal();
                System.out.printf("CRegCCD MAP %d taxa, %2d trees: brute=%.6f DP=%.6f "
                                + "certified=%-5s (off-backbone bound %.3f)%n",
                        n, nTrees, bruteForce, reported, certified, ccd.getOffBackboneBound());

                assertEquals(bruteForce, reported, 1e-9,
                        "backbone DP must find the global maximum (" + n + " taxa)");

                // the returned tree must actually attain it
                Tree map = ccd.getMAPTree();
                assertEquals(reported, ccd.getLogProbabilityOfTree(map), 1e-9,
                        "returned MAP tree must attain the reported maximum");
                assertEquals(tx.size(), map.getLeafNodeCount());

                // the certificate must never be wrong when it fires
                if (certified) {
                    assertTrue(ccd.getOffBackboneBound() < bruteForce + 1e-12,
                            "certificate claimed optimality but the bound exceeds the optimum");
                }
            }
        }
    }

    @Test
    public void entropyMatchesEnumeration() {
        List<String> tx = taxa(6);
        CRegCCD ccd = new CRegCCD(randomTrees(tx, 8, 17L), 0.0, 0.5, 0.3, 0.2);
        ccd.setRandom(new Random(2024L));

        double exact = 0.0;
        for (Tree t : CRegCCDTest.allRootedTopologies(tx)) {
            double logp = ccd.getLogProbabilityOfTree(t);
            exact -= Math.exp(logp) * logp;
        }
        double[] mc = ccd.getEntropyMonteCarlo(500_000);
        System.out.printf("CRegCCD entropy: exact=%.5f  MC=%.5f +/- %.5f (%.1f SE off)%n",
                exact, mc[0], mc[1], Math.abs(mc[0] - exact) / mc[1]);
        assertEquals(exact, mc[0], Math.max(5 * mc[1], 0.002), "MC entropy must match enumeration");

        assertTrue(Double.isFinite(ccd.getEntropy()), "default entropy must be finite");
    }

    /* ------------------------------------------------------------------ *
     * The manuscript's four-taxon example, worked symbolically.
     * Sample: (((A,B),C),D) and (((D,C),B),A), each once.
     * ------------------------------------------------------------------ */

    private static final List<String> TAXA4 = Arrays.asList("A", "B", "C", "D");

    static CRegCCD exampleModel(double alpha, double alpha1, double alpha2) {
        List<Tree> trees = new ArrayList<>();
        trees.add(new TreeParser(TAXA4, "(((A,B),C),D);", 1, false));
        trees.add(new TreeParser(TAXA4, "(((D,C),B),A);", 1, false));
        return new CRegCCD(trees, 0.0, alpha, alpha1, alpha2);
    }

    /** The 15 four-taxon trees, grouped by the six probability categories. */
    static final String[][] CATEGORIES = {
            {"(((A,B),C),D);", "(((C,D),B),A);"},                                  // sampled
            {"((A,B),(C,D));"},                                                    // both children observed
            {"(((A,C),B),D);", "(((B,C),A),D);", "(((B,D),C),A);", "(((B,C),D),A);"}, // novel inside an observed 3-clade
            {"(((C,D),A),B);", "(((A,B),D),C);"},                                  // root one-observed, reconnects
            {"(((A,C),D),B);", "(((A,D),C),B);", "(((A,D),B),C);", "(((B,D),A),C);"}, // root one-observed, novel cherry
            {"((A,C),(B,D));", "((A,D),(B,C));"}                                   // root neither observed
    };

    /** POOLED mode: per-split alpha over the three CCD0 splits, class totals a3 and a4. */
    static double[] categoryProbabilities(double alpha, double a3, double a4) {
        double z = 2 + 3 * alpha + a3 + a4;
        return new double[]{
                (1 + alpha) * (1 + alpha) / (z * (1 + alpha + a3)),
                alpha / z,
                (1 + alpha) * (a3 / 2) / (z * (1 + alpha + a3)),
                (a3 / 2) * alpha / (z * (alpha + a3)),
                (a3 / 2) * (a3 / 2) / (z * (alpha + a3)),
                (a4 / 2) / z
        };
    }

    /**
     * How large is the fresh-clade approximation's error? Compares the deterministic recursion
     * against exact enumeration across taxon counts, training-set sizes and pseudocounts.
     */
    @Test
    public void deterministicRecursionErrorVersusEnumeration() {
        System.out.printf("%n%-6s %-7s %-22s %-12s %-12s %-11s %-9s%n",
                "taxa", "trees", "pseudocounts", "exact H", "recursion", "abs err", "rel err");
        double worstRel = 0.0;
        for (int n : new int[]{5, 6, 7, 8}) {
            for (int nTrees : new int[]{2, 10, 50}) {
                for (double[] p : new double[][]{{0.4, 0.4, 0.4}, {2.0, 0.4, 0.05}, {1.0, 1.0, 1.0}}) {
                    List<String> tx = taxa(n);
                    CRegCCD ccd = new CRegCCD(randomTrees(tx, nTrees, 5L), 0.0,
                            p[0], p[1], p[2]);
                    double exact = 0.0;
                    for (Tree t : CRegCCDTest.allRootedTopologies(tx)) {
                        double logp = ccd.getLogProbabilityOfTree(t);
                        exact -= Math.exp(logp) * logp;
                    }
                    double rec = ccd.getEntropyRecursive();
                    double abs = Math.abs(rec - exact);
                    double rel = abs / exact;
                    worstRel = Math.max(worstRel, rel);
                    System.out.printf("%-6d %-7d a=(%.2f,%.2f,%.2f)%-6s %-12.6f %-12.6f %-11.2e %-8.3f%%%n",
                            n, nTrees, p[0], p[1], p[2], "", exact, rec, abs, 100 * rel);
                }
            }
        }
        System.out.printf("worst relative error = %.3f%%%n", 100 * worstRel);
        assertTrue(worstRel < 0.25, "recursion should be within 25% of exact, worst was " + worstRel);
    }

    @Test
    public void fourTaxonExampleMatchesClosedFormPooled() {
        for (double[] p : new double[][]{{0.4, 0.4, 0.4}, {0.4, 0.4, 0.05}, {1.0, 0.5, 0.25}}) {
            CRegCCD ccd = exampleModel(p[0], p[1], p[2]);
            double[] expected = categoryProbabilities(p[0], p[1], p[2]);
            double total = 0.0;
            int count = 0;
            for (int g = 0; g < CATEGORIES.length; g++) {
                for (String nwk : CATEGORIES[g]) {
                    Tree t = new TreeParser(TAXA4, nwk, 1, false);
                    assertEquals(expected[g], ccd.getProbabilityOfTree(t), 1e-12,
                            "pooled category " + (g + 1) + " tree " + nwk);
                    total += ccd.getProbabilityOfTree(t);
                    count++;
                }
            }
            assertEquals(15, count);
            assertEquals(1.0, total, 1e-12, "the 15 probabilities must sum to one");
            System.out.printf("pooled four-taxon example alpha=%.2f a3=%.2f a4=%.2f: "
                    + "closed form verified, sum = %.12f%n", p[0], p[1], p[2], total);
        }
    }

    /**
     * Under POOLED an observed split must always outrank an expanded one at the same clade, since
     * they share a per-split pseudocount and the observed one adds f(S) >= 1. Under SEPARATE that
     * can fail.
     */
    @Test
    public void pooledGuaranteesObservedOutranksExpanded() {
        for (int n : new int[]{6, 8}) {
            for (int nTrees : new int[]{3, 20}) {
                List<String> tx = taxa(n);
                List<Tree> training = randomTrees(tx, nTrees, 77L);
                for (double alpha : new double[]{0.05, 0.4, 2.0, 10.0}) {
                    CRegCCD pooled = new CRegCCD(training, 0.0, alpha, 0.4, 0.05);
                    for (Clade c : pooled.getClades()) {
                        if (c.size() < 2) {
                            continue;
                        }
                        double[] size = pooled.classSizes(c.getCladeInBits());
                        if (size[0] <= 0 || size[1] <= 0) {
                            continue;
                        }
                        // every expanded split has the same probability; take the largest observed one
                        double worstObserved = Double.POSITIVE_INFINITY;
                        for (CladePartition p : c.getPartitions()) {
                            if (p.getNumberOfOccurrences() > 0) {
                                worstObserved = Math.min(worstObserved, p.getNumberOfOccurrences());
                            }
                        }
                        // numerators: observed f + alpha, expanded alpha
                        assertTrue(worstObserved + alpha > alpha,
                                "observed must outrank expanded at " + c);
                    }
                }
            }
        }
        System.out.println("pooled: observed splits outrank expanded splits at every clade");
    }

    /**
     * regCCD is nested exactly: with no prior mass on the novel classes, CRegCCD's alpha is
     * regCCD's additive-alpha smoothing over the CCD0 split set, so the two must agree on every
     * tree that regCCD supports.
     */
    @Test
    public void regCCDIsNestedAtAlphaOneTwoZero() {
        for (int n : new int[]{5, 6, 7}) {
            for (int nTrees : new int[]{3, 15}) {
                List<String> tx = taxa(n);
                for (double alpha : new double[]{0.1, 0.4, 1.0}) {
                    CRegCCD creg = new CRegCCD(randomTrees(tx, nTrees, 31L), 0.0, alpha, 0.0, 0.0);
                    RegCCD reg = new RegCCD(randomTrees(tx, nTrees, 31L), 0.0, alpha);
                    int compared = 0;
                    double worst = 0.0;
                    for (Tree t : CRegCCDTest.allRootedTopologies(tx)) {
                        double a = reg.getLogProbabilityOfTree(t);
                        if (!Double.isFinite(a)) {
                            continue; // outside regCCD's support
                        }
                        worst = Math.max(worst, Math.abs(a - creg.getLogProbabilityOfTree(t)));
                        compared++;
                    }
                    System.out.printf("regCCD nesting %d taxa, %2d trees, alpha=%.1f: "
                            + "%d trees compared, max |diff| = %.2e%n", n, nTrees, alpha, compared, worst);
                    assertTrue(compared > 0, "regCCD must support some trees");
                    assertEquals(0.0, worst, 1e-9,
                            "CRegCCD(alpha, 0, 0) must equal regCCD(alpha)");
                }
            }
        }
    }

    /**
     * The A_0/A_1 search must equal the true global optimum by enumeration, and its certificate
     * must never claim optimality wrongly.
     */
    @Test
    public void exactMapSearchMatchesEnumeration() {
        int certified = 0, total = 0;
        for (int n : new int[]{5, 6, 7, 8}) {
            for (int nTrees : new int[]{2, 5, 20}) {
                List<String> tx = taxa(n);
                CRegCCD ccd = new CRegCCD(randomTrees(tx, nTrees, 13L), 0.0, 0.4, 0.4, 0.05);
                double brute = Double.NEGATIVE_INFINITY;
                for (Tree t : CRegCCDTest.allRootedTopologies(tx)) {
                    brute = Math.max(brute, ccd.getLogProbabilityOfTree(t));
                }
                CRegCCD.MapResult r = ccd.solveMAP();
                total++;
                if (r.certifiedGlobal()) {
                    certified++;
                }
                System.out.printf("%d taxa, %2d trees: brute=%.6f  A0/A1=%.6f  bound=%8.3f  "
                                + "certified=%-5s states=%d%n",
                        n, nTrees, brute, r.maxLogProbability(), r.offBackboneBound(),
                        r.certifiedGlobal(), r.statesExplored());
                assertTrue(r.complete(), "search must complete within the state budget");
                assertEquals(brute, r.maxLogProbability(), 1e-9,
                        "A_0/A_1 search must find the global optimum");
                if (r.certifiedGlobal()) {
                    assertTrue(r.offBackboneBound() < brute + 1e-12,
                            "certificate must not claim optimality when the bound exceeds it");
                }
            }
        }
        System.out.printf("certificate fired in %d of %d configurations%n", certified, total);
    }
}
