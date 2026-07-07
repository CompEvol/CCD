package test.ccd.model;

import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import ccd.model.HeightSettingStrategy;
import ccd.model.KRegCCD;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Tests that {@link KRegCCD}'s discounted MAP DP returns the true full-support maximum.
 *
 * <p>The inherited DP maximises only the red-only log-CCP sum and omits the per-clade
 * {@code log(1 - mu - tail(C))} escape discount, so it both overstates the maximum and can pick the
 * wrong topology (the discount is incurred once per reservable internal clade, and reservability
 * differs across topologies). On a taxon set small enough to enumerate every rooted topology we can
 * pin the DP against a brute-force maximum of {@link KRegCCD#getLogProbabilityOfTree}.
 */
public class KRegCCDMapTreeTest {

    private static final List<String> TAXA = Arrays.asList("A", "B", "C", "D", "E");

    private static Tree parse(String newick) {
        return new TreeParser(TAXA, newick, 1, false);
    }

    /** Training trees giving several size-3/4 observed clades (so reservable clades exist and the
     * discount actually bites), with enough split variety that the MAP is non-trivial. */
    private static List<Tree> trainingTrees() {
        List<Tree> trees = new ArrayList<>();
        trees.add(parse("((((A:1,B:1):1,C:1):1,D:1):1,E:1):0;"));
        trees.add(parse("((((A:1,B:1):1,C:1):1,D:1):1,E:1):0;"));
        trees.add(parse("(((A:1,B:1):1,(C:1,D:1):1):1,E:1):0;"));
        trees.add(parse("(((A:1,B:1):1,C:1):1,(D:1,E:1):1):0;"));
        return trees;
    }

    private static KRegCCD build(double mu, KRegCCD.TailMode tailMode) {
        return new KRegCCD(trainingTrees(), 0.0, mu, 0.4, 2, tailMode);
    }

    /**
     * The discounted DP value equals the brute-force maximum of {@code getLogProbabilityOfTree} over
     * every rooted topology, the returned MAP tree attains it, and {@code getMaxTreeProbability} is
     * its exp -- across several mu and both tail modes (so the discount {@code 1 - mu - tail} varies).
     */
    @Test
    public void testDiscountedMapMatchesBruteForce() {
        for (KRegCCD.TailMode tail : new KRegCCD.TailMode[]{KRegCCD.TailMode.NONE, KRegCCD.TailMode.BOUND}) {
            for (double mu : new double[]{0.005, 0.02, 0.05}) {
                KRegCCD ccd = build(mu, tail);

                double bruteMax = Double.NEGATIVE_INFINITY;
                for (Tree t : allRootedTopologies()) {
                    bruteMax = Math.max(bruteMax, ccd.getLogProbabilityOfTree(t));
                }

                double dpMax = ccd.getMaxLogTreeProbability();
                assertEquals(bruteMax, dpMax, 1e-9,
                        "discounted DP must equal the brute-force max log-prob (mu=" + mu + ", tail=" + tail + ")");

                // The MAP tree the DP reports must actually attain that maximum.
                Tree map = ccd.getMAPTree(HeightSettingStrategy.None);
                assertEquals(dpMax, ccd.getLogProbabilityOfTree(map), 1e-9,
                        "getLogProbabilityOfTree(getMAPTree()) must equal getMaxLogTreeProbability() (mu="
                                + mu + ", tail=" + tail + ")");

                assertEquals(Math.exp(dpMax), ccd.getMaxTreeProbability(), 1e-12,
                        "getMaxTreeProbability must be exp(getMaxLogTreeProbability) (mu=" + mu + ", tail=" + tail + ")");

                // No tree can beat the reported maximum.
                for (Tree t : allRootedTopologies()) {
                    assertTrue(ccd.getLogProbabilityOfTree(t) <= dpMax + 1e-9,
                            "no topology may exceed the reported max (mu=" + mu + ", tail=" + tail + ")");
                }
            }
        }
    }

    /* --------------------------------------------------------------------- */

    private static List<Tree> allRootedTopologies() {
        List<Tree> trees = new ArrayList<>();
        for (String shape : shapes(TAXA)) {
            trees.add(parse(shape + ";"));
        }
        return trees;
    }

    /** Newick fragments (with dummy branch lengths) of all rooted binary trees over {@code taxa};
     * each unordered bipartition is counted once by forcing the first taxon into the left side. */
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
                continue; // left always holds `first`; require right non-empty too
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
