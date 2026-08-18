package ccd.model;

import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import ccd.model.bitsets.BitSet;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/** {@link CRegCCD}: exact normalisation, class-size arithmetic, and full support. */
public class CRegCCDTest {

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

    private static double totalMass(CRegCCD ccd, List<String> taxa) {
        double sum = 0.0;
        for (Tree t : allRootedTopologies(taxa)) {
            sum += Math.exp(ccd.getLogProbabilityOfTree(t));
        }
        return sum;
    }

    @Test
    public void exactlyNormalisedOverTreeSpace() {
        for (int n : new int[]{4, 5, 6, 7}) {
            List<String> tx = taxa(n);
            for (int nTrees : new int[]{1, 3, 20}) {
                List<Tree> training = randomTrees(tx, nTrees, 7L);
                CRegCCD ccd = new CRegCCD(training, 0.0, 0.4, 0.4, 0.4);
                double mass = totalMass(ccd, tx);
                System.out.printf("CRegCCD %d taxa, %2d training trees: totalMass = %.12f%n",
                        n, nTrees, mass);
                assertEquals(1.0, mass, 1e-9,
                        "CRegCCD must be exactly normalised (" + n + " taxa, " + nTrees + " trees)");
            }
        }
    }

    @Test
    public void normalisedAcrossPseudocountChoices() {
        List<String> tx = taxa(6);
        List<Tree> training = randomTrees(tx, 10, 11L);
        CRegCCD ccd = new CRegCCD(training, 0.0);
        for (double[] p : new double[][]{
                {0.0, 0.4, 0.4, 0.4},
                {1.0, 1.0, 1.0, 1.0},
                {0.0, 2.0, 0.5, 0.05},
                {0.3, 0.01, 5.0, 0.2}}) {
            double sum = 0.0;
            for (Tree t : allRootedTopologies(tx)) {
                sum += Math.exp(ccd.getLogProbabilityOfTree(t, p[1], p[2], p[3]));
            }
            System.out.printf("CRegCCD 6 taxa a=(%.2f,%.2f,%.2f,%.2f): totalMass = %.12f%n",
                    p[0], p[1], p[2], p[3], sum);
            assertEquals(1.0, sum, 1e-9, "normalisation must hold for any pseudocounts");
        }
    }

    @Test
    public void classSizesPartitionAllBipartitions() {
        List<String> tx = taxa(8);
        List<Tree> training = randomTrees(tx, 50, 3L);
        CRegCCD ccd = new CRegCCD(training, 0.0);
        CCD0 ccd0 = new CCD0(training, 0);
        int checked = 0;
        for (Clade c : ccd.getClades()) {
            if (c.size() < 2) {
                continue;
            }
            double[] s = ccd.classSizes(c.getCladeInBits());
            double total = Math.pow(2.0, c.size() - 1) - 1.0;
            assertEquals(total, s[0] + s[1] + s[2] + s[3], 1e-6,
                    "class sizes must partition all bipartitions of " + c);

            // |A_1| + |A_2| is exactly the CCD0 split count (observed + expanded)
            Clade c0 = ccd0.getClade(c.getCladeInBits());
            if (c0 != null) {
                assertEquals(c0.getPartitions().size(), (int) Math.round(s[0] + s[1]),
                        "|A_1|+|A_2| must equal the CCD0 partition count for " + c);
            }
            checked++;
        }
        System.out.printf("class-size arithmetic verified on %d clades%n", checked);
        assertTrue(checked > 0);
    }

    @Test
    public void fullSupportOnHeldOutTrees() {
        List<String> tx = taxa(30);
        List<Tree> training = randomTrees(tx, 100, 5L);
        List<Tree> heldOut = randomTrees(tx, 100, 99L);
        CRegCCD ccd = new CRegCCD(training, 0.0);
        int covered = 0;
        double sum = 0.0;
        for (Tree t : heldOut) {
            double lp = ccd.getLogProbabilityOfTree(t);
            if (Double.isFinite(lp)) {
                covered++;
                sum += lp;
            }
        }
        System.out.printf("CRegCCD 30 taxa: %d/%d held-out trees with finite logP, mean = %.2f%n",
                covered, heldOut.size(), sum / covered);
        assertEquals(heldOut.size(), covered, "every held-out tree must have positive probability");
    }

    /**
     * If the simulator draws from q AND reports the correct log q, then the mean sampled -log q
     * equals the entropy computed by enumeration with the scorer. This checks the sampler and the
     * scorer are the same distribution without needing per-topology counts.
     */
    @Test
    public void samplerEntropyMatchesScorer() {
        List<String> tx = taxa(6);
        List<Tree> training = randomTrees(tx, 8, 17L);
        CRegCCD ccd = new CRegCCD(training, 0.0, 0.5, 0.3, 0.2);
        ccd.setRandom(new Random(31337L));

        double mass = 0.0;
        double h = 0.0;
        for (Tree t : allRootedTopologies(tx)) {
            double logp = ccd.getLogProbabilityOfTree(t);
            double p = Math.exp(logp);
            mass += p;
            h -= p * logp;
        }
        assertEquals(1.0, mass, 1e-9, "scored distribution must be normalised");

        int n = 1_000_000;
        double s1 = 0.0;
        double s2 = 0.0;
        for (int i = 0; i < n; i++) {
            double logp = ccd.sampleTreeLogProbability();
            s1 += -logp;
            s2 += logp * logp;
        }
        double hHat = s1 / n;
        double se = Math.sqrt(Math.max(0, s2 / n - hHat * hHat) / n);
        System.out.printf("CRegCCD sampler: H_enum=%.5f  H_MC=%.5f +/- %.5f (%.1f SE off)%n",
                h, hHat, se, Math.abs(hHat - h) / se);
        assertEquals(h, hHat, Math.max(5 * se, 0.005),
                "sampler entropy must match the scorer's enumerated entropy");
    }

    /**
     * The direct check: sampled topology frequencies must match the scored probabilities. Compares
     * every topology whose expected count is large enough for a normal approximation.
     */
    @Test
    public void sampledFrequenciesMatchScoredProbabilities() {
        List<String> tx = taxa(5);
        List<Tree> training = randomTrees(tx, 5, 23L);
        CRegCCD ccd = new CRegCCD(training, 0.0, 0.5, 0.3, 0.2);
        ccd.setRandom(new Random(4242L));

        java.util.Map<String, Double> expected = new java.util.HashMap<>();
        for (Tree t : allRootedTopologies(tx)) {
            expected.put(canonical(t), Math.exp(ccd.getLogProbabilityOfTree(t)));
        }

        int n = 1_000_000;
        java.util.Map<String, Integer> counts = new java.util.HashMap<>();
        for (int i = 0; i < n; i++) {
            counts.merge(canonical(ccd.sampleTree()), 1, Integer::sum);
        }

        int checked = 0;
        double worst = 0.0;
        String worstKey = null;
        for (var e : expected.entrySet()) {
            double p = e.getValue();
            if (n * p < 30) {
                continue; // too rare for a normal approximation
            }
            int obs = counts.getOrDefault(e.getKey(), 0);
            double z = Math.abs(obs - n * p) / Math.sqrt(n * p * (1 - p));
            if (z > worst) {
                worst = z;
                worstKey = e.getKey();
            }
            checked++;
        }
        System.out.printf("CRegCCD frequencies: %d topologies checked, worst |z| = %.2f (%s)%n",
                checked, worst, worstKey);
        assertTrue(checked >= 10, "expected a reasonable number of comparable topologies");
        assertTrue(worst < 4.5, "sampled frequencies must match scored probabilities, worst z = " + worst);

        // no sampled topology may fall outside the enumerated support
        for (String key : counts.keySet()) {
            assertTrue(expected.containsKey(key), "sampler produced an unrecognised topology " + key);
        }
    }

    @Test
    public void sampledTreesAreValidAndSelfConsistent() {
        List<String> tx = taxa(20);
        List<Tree> training = randomTrees(tx, 50, 8L);
        CRegCCD ccd = new CRegCCD(training, 0.0);
        ccd.setRandom(new Random(5L));
        for (int i = 0; i < 200; i++) {
            Tree t = ccd.sampleTree();
            assertEquals(tx.size(), t.getLeafNodeCount(), "sampled tree must have every taxon");
            assertEquals(2 * tx.size() - 1, t.getNodeCount(), "sampled tree must be binary");
            double stamped = (Double) t.getRoot().getMetaData(CCD1.LOG_PROB_SUBTREE_KEY);
            assertEquals(ccd.getLogProbabilityOfTree(t), stamped, 1e-9,
                    "stamped log probability must equal the scorer's");
        }
        System.out.println("CRegCCD: 200 sampled 20-taxon trees valid and self-consistent");
    }

    /** Canonical topology key: nested sorted taxon-index sets. */
    private static String canonical(Tree t) {
        return canonical(t.getRoot());
    }

    private static String canonical(beast.base.evolution.tree.Node v) {
        if (v.isLeaf()) {
            return String.valueOf(v.getNr());
        }
        String a = canonical(v.getChild(0));
        String b = canonical(v.getChild(1));
        return (a.compareTo(b) <= 0) ? "(" + a + "," + b + ")" : "(" + b + "," + a + ")";
    }

    /* --------------------------------------------------------------------- */

    static List<Tree> allRootedTopologies(List<String> taxa) {
        List<Tree> trees = new ArrayList<>();
        for (String shape : shapes(taxa)) {
            trees.add(new TreeParser(taxa, shape + ";", 1, false));
        }
        return trees;
    }

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
}
