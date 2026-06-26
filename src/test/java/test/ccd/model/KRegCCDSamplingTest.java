package test.ccd.model;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import ccd.model.AbstractCCD;
import ccd.model.HeightSettingStrategy;
import ccd.model.KRegCCD;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Tests that {@link KRegCCD} can <em>sample</em> from the full-support distribution it scores
 * (not just the truncated red-only backbone the inherited sampler produces). The sampler escapes
 * into blue regions that contain novel clades, and its empirical distribution matches
 * {@code exp(getLogProbabilityOfTree)}.
 *
 * <p>On four taxa every blue region's boundary parts are leaves or cherries, which reserve no
 * escape mass, so the sampler is <em>exact</em>: the model normalises to 1 and sampled-tree
 * frequencies match the score.
 */
public class KRegCCDSamplingTest {

    private static final List<String> TAXA = Arrays.asList("A", "B", "C", "D");

    private static Tree parse(String newick) {
        return new TreeParser(TAXA, newick, 1, false);
    }

    /** Observed clades: AB, CD (size 2) and ABC (size 3); leaves room for novel clades like
     * {A,C}, {A,B,D}, etc. */
    private static List<Tree> trainingTrees() {
        List<Tree> trees = new ArrayList<>();
        trees.add(parse("(((A:1,B:1):1,C:1):1,D:1):0;"));
        trees.add(parse("((A:1,B:1):1,(C:1,D:1):1):0;"));
        return trees;
    }

    private static KRegCCD build(double mu, KRegCCD.TailMode tailMode) {
        // k = 2 (default reserve depth); on four taxa this already covers every region order.
        return new KRegCCD(trainingTrees(), 0.0, mu, 0.4, 2, tailMode);
    }

    /* --------------------------------------------------------------------- */

    /**
     * With the tail set to zero the full-support model is exactly normalised on four taxa:
     * summing {@code exp(getLogProbabilityOfTree)} over every rooted binary topology gives 1.
     * This is the strongest correctness check — it confirms the red discount, the blue-region
     * weights, and the pathcount normalisation all line up.
     */
    @Test
    public void testModelNormalisesOnFourTaxa() {
        KRegCCD ccd = build(0.05, KRegCCD.TailMode.NONE);
        double total = 0.0;
        for (Tree t : allRootedTopologies()) {
            double p = Math.exp(ccd.getLogProbabilityOfTree(t));
            assertTrue(p > 0.0, "full support: every topology must have positive probability");
            total += p;
        }
        assertEquals(1.0, total, 1e-6,
                "exp(score) should sum to 1 over all topologies (exact on four taxa)");
    }

    /**
     * Every sampled tree's stamped log-probability (accumulated in the node metadata while
     * sampling) must equal {@link KRegCCD#getLogProbabilityOfTree} re-scoring that same tree.
     * This validates the metadata stamping independently of Monte-Carlo noise, in both modes.
     */
    @Test
    public void testSampleLogProbabilityMatchesScore() {
        for (KRegCCD.SamplingFidelity fidelity : KRegCCD.SamplingFidelity.values()) {
            KRegCCD ccd = build(0.05, KRegCCD.TailMode.NONE);
            ccd.setSamplingFidelity(fidelity);
            ccd.setRandom(new Random(42));
            for (int i = 0; i < 2000; i++) {
                Tree sampled = ccd.sampleTree(HeightSettingStrategy.None);
                double stamped = (double) sampled.getRoot().getMetaData(AbstractCCD.LOG_PROB_SUBTREE_KEY);
                double rescored = ccd.getLogProbabilityOfTree(sampled);
                assertEquals(rescored, stamped, 1e-9,
                        "stamped log-prob must equal getLogProbabilityOfTree [" + fidelity + "]");
            }
        }
    }

    /**
     * Per-topology empirical frequencies match {@code exp(score)} on four taxa (exact target).
     * At least one novel-clade topology must be both expected and observed, proving the sampler
     * actually escapes the observed support.
     */
    @Test
    public void testEmpiricalFrequenciesMatchScore() {
        double mu = 0.05;
        int n = 300_000;
        KRegCCD ccd = build(mu, KRegCCD.TailMode.NONE);
        ccd.setRandom(new Random(20260605L));

        Map<String, Integer> counts = new HashMap<>();
        for (int i = 0; i < n; i++) {
            counts.merge(topologyKey(ccd.sampleTree(HeightSettingStrategy.None)), 1, Integer::sum);
        }

        boolean sawNovelTopology = false;
        for (Tree t : allRootedTopologies()) {
            double expected = Math.exp(ccd.getLogProbabilityOfTree(t));
            String key = topologyKey(t);
            double observed = counts.getOrDefault(key, 0) / (double) n;
            // binomial standard error; only assert where we have enough samples to be meaningful
            if (expected * n >= 100) {
                double se = Math.sqrt(expected * (1 - expected) / n);
                assertEquals(expected, observed, 5 * se + 0.002,
                        "frequency of topology " + key + " should match exp(score)");
            }
            if (containsNovelClade(ccd, t) && expected * n >= 100) {
                sawNovelTopology = true;
                assertTrue(counts.getOrDefault(key, 0) > 0,
                        "a novel-clade topology with non-trivial mass must actually be sampled");
            }
        }
        assertTrue(sawNovelTopology,
                "the sampler must produce novel-clade topologies (the whole point of full support)");
    }

    /**
     * Aggregate calibration: the fraction of sampled trees that contain a novel clade matches the
     * model's closed-form novel-clade probability {@link KRegCCD#getNovelCladeProbability()}.
     */
    @Test
    public void testEmpiricalNovelFractionMatchesModel() {
        double mu = 0.05;
        int n = 300_000;
        KRegCCD ccd = build(mu, KRegCCD.TailMode.NONE);
        ccd.setRandom(new Random(7L));

        double modelNovel = ccd.getNovelCladeProbability();
        assertTrue(modelNovel > 0 && modelNovel < 1, "novel-clade probability should be interior");

        int novel = 0;
        for (int i = 0; i < n; i++) {
            if (containsNovelClade(ccd, ccd.sampleTree(HeightSettingStrategy.None))) {
                novel++;
            }
        }
        double empirical = novel / (double) n;
        double se = Math.sqrt(modelNovel * (1 - modelNovel) / n);
        assertEquals(modelNovel, empirical, 5 * se + 0.003,
                "empirical novel-clade fraction should match getNovelCladeProbability()");
    }

    /**
     * Sampled trees are structurally valid (correct leaf count, positive branch lengths under
     * height strategies that set heights), some contain novel clades, and all have positive
     * probability under the model. Exercises the CommonAncestorHeights path that {@code CCDSampler}
     * uses, plus the One strategy.
     */
    @Test
    public void testSampledTreesAreValid() {
        for (HeightSettingStrategy strategy :
                List.of(HeightSettingStrategy.One, HeightSettingStrategy.CommonAncestorHeights)) {
            KRegCCD ccd = build(0.1, KRegCCD.TailMode.BOUND);
            ccd.setRandom(new Random(99));
            boolean anyNovel = false;
            for (int i = 0; i < 500; i++) {
                Tree sampled = ccd.sampleTree(strategy);
                assertEquals(TAXA.size(), sampled.getLeafNodeCount(),
                        "sampled tree must have all taxa [" + strategy + "]");
                assertPositiveBranches(sampled.getRoot(), strategy);
                assertTrue(Math.exp(ccd.getLogProbabilityOfTree(sampled)) > 0.0,
                        "sampled tree must have positive probability [" + strategy + "]");
                anyNovel |= containsNovelClade(ccd, sampled);
            }
            assertTrue(anyNovel, "some sampled trees should contain a novel clade [" + strategy + "]");
        }
    }

    /* --------------------------------------------------------------------- *
     * helpers
     * --------------------------------------------------------------------- */

    private static boolean containsNovelClade(KRegCCD ccd, Tree t) {
        return ccd.containsNovelClade(t);
    }

    /** Topology-invariant key: the sorted set of leaf-index sets at the internal nodes. */
    private static String topologyKey(Tree t) {
        List<String> clades = new ArrayList<>();
        cladeKeys(t.getRoot(), clades);
        Collections.sort(clades);
        return clades.toString();
    }

    private static java.util.BitSet cladeKeys(Node v, List<String> clades) {
        java.util.BitSet b = new java.util.BitSet();
        if (v.isLeaf()) {
            b.set(v.getNr());
        } else {
            for (Node c : v.getChildren()) {
                b.or(cladeKeys(c, clades));
            }
            clades.add(b.toString());
        }
        return b;
    }

    private static void assertPositiveBranches(Node v, HeightSettingStrategy strategy) {
        for (Node c : v.getChildren()) {
            double bl = v.getHeight() - c.getHeight();
            assertTrue(!Double.isNaN(bl) && bl > -1e-9,
                    "branch length must be non-negative and finite [" + strategy + "], got " + bl);
            assertPositiveBranches(c, strategy);
        }
    }

    /** Every rooted binary topology over {@link #TAXA} (15 of them for four taxa). */
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
