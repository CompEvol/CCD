package ccd.model;

import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import ccd.model.bitsets.BitSet;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;

/** MRegCCD must reproduce MRegCCDSlow exactly: same boundary counts, same tree probabilities. */
public class MRegCCDAgreementTest {

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

    @Test
    public void boundaryCountsAndProbabilitiesAgree() {
        int cladesChecked = 0;
        for (int n : new int[]{6, 8, 10, 14}) {
            for (int nTrees : new int[]{5, 25}) {
                for (int depth : new int[]{2, 3, 4}) {
                    List<String> tx = taxa(n);
                    MRegCCDSlow slow = new MRegCCDSlow(randomTrees(tx, nTrees, 21L), 0.0, 0.02, depth, true);
                    MRegCCD fast = new MRegCCD(randomTrees(tx, nTrees, 21L), 0.0, 0.02, depth, true);

                    for (Clade c : slow.getClades()) {
                        BitSet cb = c.getCladeInBits();
                        assertArrayEquals(slow.countsFor(cb), fast.countsFor(cb),
                                "boundary counts differ at " + n + " taxa, depth " + depth
                                        + ", clade " + cb);
                        cladesChecked++;
                    }

                    List<Tree> probe = randomTrees(tx, 40, 99L);
                    for (Tree t : probe) {
                        assertEquals(slow.getLogProbabilityOfTree(t), fast.getLogProbabilityOfTree(t),
                                1e-9, "tree probability differs at " + n + " taxa, depth " + depth);
                    }
                }
            }
        }
        System.out.printf("MRegCCD agrees with MRegCCDSlow on %d clades and every probed tree%n",
                cladesChecked);
    }

    /**
     * The exact-normalisation guarantee is stated for the model, and MRegCCDTest checks it on the
     * reference implementation. Since MRegCCD is the class callers get, and its fast path computes
     * the counts that the reserve is solved from, check the property directly on it too.
     */
    @Test
    public void fastImplementationIsExactlyNormalisedAtFullDepth() {
        for (int n : new int[]{5, 6}) {
            for (double mu : new double[]{0.02, 0.1, 0.25}) {
                List<String> tx = taxa(n);
                // full reserve depth: no omitted tail, so the model must normalise exactly
                MRegCCD m = new MRegCCD(randomTrees(tx, 6, 5L), 0.0, mu, tx.size(), false);
                double sum = 0.0;
                for (Tree t : allRootedTopologies(tx)) {
                    sum += Math.exp(m.getLogProbabilityOfTree(t));
                }
                System.out.printf("MRegCCD %d taxa mu=%.2f full-depth SUM = %.12f%n", n, mu, sum);
                assertEquals(1.0, sum, 1e-9,
                        "MRegCCD at full reserve depth must be exactly normalised");
            }
        }
    }

    private static List<Tree> allRootedTopologies(List<String> taxa) {
        List<Tree> out = new ArrayList<>();
        for (String shape : shapes(taxa)) {
            out.add(new TreeParser(taxa, shape + ";", 1, false));
        }
        return out;
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
