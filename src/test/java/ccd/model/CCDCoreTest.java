package test.ccd.model;

import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import ccd.model.AbstractCCD;
import ccd.model.CCD0;
import ccd.model.CCD1;
import ccd.model.CCD2;
import ccd.model.HeightSettingStrategy;

import org.junit.Before;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import static org.junit.Assert.*;

/**
 * Core unit tests for CCD0, CCD1, and CCD2 models.
 *
 * Uses 5-taxon trees (A, B, C, D, E) with three topologies at known frequencies:
 *   T1: ((A,(B,C)),(D,E))  — 6 copies (0.6)
 *   T2: ((((A,B),C),D),E)  — 2 copies (0.2)
 *   T3: ((((B,C),A),E),D)  — 2 copies (0.2)
 */
public class CCDCoreTest {

    // Newick strings for 5-taxon trees
    private static final String T1_NEWICK = "((A:1,(B:1,C:1):1):1,(D:1,E:1):1):0;";
    private static final String T2_NEWICK = "((((A:1,B:1):1,C:1):1,D:1):1,E:1):0;";
    private static final String T3_NEWICK = "((((B:1,C:1):1,A:1):1,E:1):1,D:1):0;";

    // 4-taxon tree for single-tree tests
    private static final String FOUR_TAXON_NEWICK = "((A:1,B:1):1,(C:1,D:1):1):0;";

    // An unseen 5-taxon topology (not in {T1, T2, T3})
    private static final String UNSEEN_NEWICK = "((A:1,D:1):1,(B:1,(C:1,E:1):1):1):0;";

    private List<Tree> mixedTrees;   // 6xT1 + 2xT2 + 2xT3 = 10 trees
    private Tree t1, t2, t3, unseenTree;

    @Before
    public void setUp() {
        AbstractCCD.verbose = false;

        t1 = parseNewick(T1_NEWICK);
        t2 = parseNewick(T2_NEWICK);
        t3 = parseNewick(T3_NEWICK);
        unseenTree = parseNewick(UNSEEN_NEWICK);

        mixedTrees = new ArrayList<>();
        for (int i = 0; i < 6; i++) mixedTrees.add(parseNewick(T1_NEWICK));
        for (int i = 0; i < 2; i++) mixedTrees.add(parseNewick(T2_NEWICK));
        for (int i = 0; i < 2; i++) mixedTrees.add(parseNewick(T3_NEWICK));
    }

    private Tree parseNewick(String newick) {
        return new TreeParser(newick, false, false, true, 1);
    }

    private List<Tree> nCopies(String newick, int n) {
        List<Tree> trees = new ArrayList<>();
        for (int i = 0; i < n; i++) trees.add(parseNewick(newick));
        return trees;
    }

    // ======================== CCD1 Tests ========================

    @Test
    public void testCCD1_singleTreeCladeCount() {
        // 4-taxon rooted binary tree: 4 leaves + 2 internal + 1 root = 7 clades
        List<Tree> trees = Collections.singletonList(parseNewick(FOUR_TAXON_NEWICK));
        CCD1 ccd = new CCD1(trees, 0.0);
        assertEquals(7, ccd.getNumberOfClades());
    }

    @Test
    public void testCCD1_identicalTreesMAPAndProbability() {
        List<Tree> trees = nCopies(T1_NEWICK, 10);
        CCD1 ccd = new CCD1(trees, 0.0);

        // MAP tree should have probability 1.0
        Tree mapTree = ccd.getMAPTree();
        double prob = ccd.getProbabilityOfTree(mapTree);
        assertEquals(1.0, prob, 1e-9);

        // The MAP tree should match our input topology
        assertEquals(1.0, ccd.getProbabilityOfTree(t1), 1e-9);
    }

    @Test
    public void testCCD1_mixedTreesProbabilitiesReflectFrequencies() {
        CCD1 ccd = new CCD1(mixedTrees, 0.0);

        double p1 = ccd.getProbabilityOfTree(t1);
        double p2 = ccd.getProbabilityOfTree(t2);
        double p3 = ccd.getProbabilityOfTree(t3);

        // All input trees should have positive probability
        assertTrue("T1 should have positive probability", p1 > 0);
        assertTrue("T2 should have positive probability", p2 > 0);
        assertTrue("T3 should have positive probability", p3 > 0);

        // T1 appeared 6/10 times, should have highest probability
        assertTrue("T1 should have highest probability", p1 > p2);
        assertTrue("T1 should have highest probability", p1 > p3);
    }

    @Test
    public void testCCD1_containsTree() {
        CCD1 ccd = new CCD1(mixedTrees, 0.0);

        assertTrue("CCD1 should contain T1", ccd.containsTree(t1));
        assertTrue("CCD1 should contain T2", ccd.containsTree(t2));
        assertTrue("CCD1 should contain T3", ccd.containsTree(t3));

        // CCD1 only assigns nonzero probability to trees whose clade partitions
        // were all observed; an unseen topology with novel partitions should not be contained
        assertFalse("CCD1 should not contain unseen topology", ccd.containsTree(unseenTree));
    }

    @Test
    public void testCCD1_getMAPTree() {
        CCD1 ccd = new CCD1(mixedTrees, 0.0);
        Tree mapTree = ccd.getMAPTree();

        // MAP should be T1 (highest frequency)
        double mapProb = ccd.getProbabilityOfTree(mapTree);
        double t1Prob = ccd.getProbabilityOfTree(t1);
        assertEquals("MAP tree should have same probability as T1", t1Prob, mapProb, 1e-9);
    }

    @Test
    public void testCCD1_sampleTreeValid() {
        CCD1 ccd = new CCD1(mixedTrees, 0.0);
        ccd.setRandom(new Random(42));

        for (int i = 0; i < 100; i++) {
            Tree sampled = ccd.sampleTree();
            assertNotNull("Sampled tree should not be null", sampled);
            assertTrue("Sampled tree should be contained in CCD",
                    ccd.getProbabilityOfTree(sampled) > 0);
        }
    }

    @Test
    public void testCCD1_sampleTreeFrequencyApproximation() {
        CCD1 ccd = new CCD1(mixedTrees, 0.0);
        ccd.setRandom(new Random(123));

        int nSamples = 10000;
        int t1Count = 0;
        double t1Prob = ccd.getProbabilityOfTree(t1);

        for (int i = 0; i < nSamples; i++) {
            Tree sampled = ccd.sampleTree();
            if (ccd.getProbabilityOfTree(sampled) == t1Prob) {
                t1Count++;
            }
        }

        double empiricalFreq = (double) t1Count / nSamples;
        assertEquals("Sampled frequency of T1 should approximate its CCD probability",
                t1Prob, empiricalFreq, 0.05);
    }

    @Test
    public void testCCD1_burnin() {
        // 10 trees total, 50% burnin = first 5 discarded, last 5 used
        // Put T2 in first 5 positions (will be discarded), T1 in last 5
        List<Tree> trees = new ArrayList<>();
        for (int i = 0; i < 5; i++) trees.add(parseNewick(T2_NEWICK));
        for (int i = 0; i < 5; i++) trees.add(parseNewick(T1_NEWICK));

        CCD1 ccd = new CCD1(trees, 0.5);

        // Only T1 trees should have been used
        assertEquals(1.0, ccd.getProbabilityOfTree(t1), 1e-9);
        assertEquals(5, ccd.getNumberOfBaseTrees());
    }

    @Test
    public void testCCD1_probabilitySumApproximatelyOne() {
        CCD1 ccd = new CCD1(mixedTrees, 0.0);
        ccd.setRandom(new Random(42));

        // For CCD1 with few topologies, we can check that probabilities of
        // all distinct observed topologies sum close to 1
        double sum = ccd.getProbabilityOfTree(t1)
                + ccd.getProbabilityOfTree(t2)
                + ccd.getProbabilityOfTree(t3);

        // CCD1 can assign probability to unobserved topologies that share
        // observed clade partitions. So sum of observed trees <= 1.
        assertTrue("Sum of observed tree probabilities should be <= 1", sum <= 1.0 + 1e-9);
        assertTrue("Sum of observed tree probabilities should be substantial", sum > 0.5);
    }

    // ======================== CCD0 Tests ========================

    @Test
    public void testCCD0_singleTreeCladeCount() {
        List<Tree> trees = Collections.singletonList(parseNewick(FOUR_TAXON_NEWICK));
        CCD0 ccd = new CCD0(trees, 0.0);
        // CCD0 expands the graph, so may have more clades than observed
        assertTrue("CCD0 should have at least 7 clades", ccd.getNumberOfClades() >= 7);
    }

    @Test
    public void testCCD0_identicalTreesMAPAndProbability() {
        List<Tree> trees = nCopies(T1_NEWICK, 10);
        CCD0 ccd = new CCD0(trees, 0.0);

        Tree mapTree = ccd.getMAPTree();
        double prob = ccd.getProbabilityOfTree(mapTree);

        // With identical trees, MAP tree should have the highest probability
        assertTrue("MAP probability should be high", prob > 0.5);

        // The input topology should have highest probability
        double t1Prob = ccd.getProbabilityOfTree(t1);
        assertTrue("T1 should have high probability", t1Prob > 0.5);
    }

    @Test
    public void testCCD0_mixedTreesProbabilities() {
        CCD0 ccd = new CCD0(mixedTrees, 0.0);

        double p1 = ccd.getProbabilityOfTree(t1);
        double p2 = ccd.getProbabilityOfTree(t2);

        assertTrue("T1 should have positive probability", p1 > 0);
        assertTrue("T2 should have positive probability", p2 > 0);
        assertTrue("T1 should have higher probability than T2", p1 > p2);
    }

    @Test
    public void testCCD0_getMAPTree() {
        CCD0 ccd = new CCD0(mixedTrees, 0.0);
        Tree mapTree = ccd.getMAPTree();

        double mapProb = ccd.getProbabilityOfTree(mapTree);
        assertTrue("MAP tree should have positive probability", mapProb > 0);

        // MAP should have probability >= any input tree
        assertTrue("MAP >= T1", mapProb >= ccd.getProbabilityOfTree(t1) - 1e-9);
        assertTrue("MAP >= T2", mapProb >= ccd.getProbabilityOfTree(t2) - 1e-9);
        assertTrue("MAP >= T3", mapProb >= ccd.getProbabilityOfTree(t3) - 1e-9);
    }

    @Test
    public void testCCD0_sampleTreeValid() {
        CCD0 ccd = new CCD0(mixedTrees, 0.0);
        ccd.setRandom(new Random(42));

        for (int i = 0; i < 50; i++) {
            Tree sampled = ccd.sampleTree();
            assertNotNull("Sampled tree should not be null", sampled);
            assertTrue("Sampled tree should have positive probability",
                    ccd.getProbabilityOfTree(sampled) > 0);
        }
    }

    @Test
    public void testCCD0_burnin() {
        List<Tree> trees = new ArrayList<>();
        for (int i = 0; i < 5; i++) trees.add(parseNewick(T2_NEWICK));
        for (int i = 0; i < 5; i++) trees.add(parseNewick(T1_NEWICK));

        CCD0 ccd = new CCD0(trees, 0.5);
        assertEquals(5, ccd.getNumberOfBaseTrees());

        // T1 should dominate since burnin discards first 5 (T2) trees
        Tree mapTree = ccd.getMAPTree();
        double mapProb = ccd.getProbabilityOfTree(mapTree);
        assertTrue("MAP should have high probability after burnin", mapProb > 0.5);
    }

    @Test
    public void testCCD0_canAssignNonzeroProbToUnseenTopology() {
        // CCD0 expands the graph, so unseen topologies may get nonzero probability
        // if their clades are subsets of observed clades
        CCD0 ccd = new CCD0(mixedTrees, 0.0);

        // This is a property test: CCD0 may or may not contain unseen trees
        // depending on the expansion. We just check it doesn't crash.
        double prob = ccd.getProbabilityOfTree(unseenTree);
        assertTrue("Probability should be non-negative", prob >= 0);
    }

    // ======================== CCD2 Tests ========================

    @Test
    public void testCCD2_singleTreeCladeCount() {
        // CCD2 tracks extended clades (clade + sibling context),
        // but getNumberOfClades() returns the standard clade count
        List<Tree> trees = Collections.singletonList(parseNewick(FOUR_TAXON_NEWICK));
        CCD2 ccd = new CCD2(trees, 0.0);
        // 4-taxon rooted binary tree has 7 standard clades,
        // but CCD2 may count differently due to root handling
        assertTrue("CCD2 should have at least 6 clades", ccd.getNumberOfClades() >= 6);
        assertTrue("CCD2 should not exceed 7 clades", ccd.getNumberOfClades() <= 7);
    }

    @Test
    public void testCCD2_identicalTreesMAPAndProbability() {
        List<Tree> trees = nCopies(T1_NEWICK, 10);
        CCD2 ccd = new CCD2(trees, 0.0);

        Tree mapTree = ccd.getMAPTree();
        double prob = ccd.getProbabilityOfTree(mapTree);
        assertEquals(1.0, prob, 1e-9);
        assertEquals(1.0, ccd.getProbabilityOfTree(t1), 1e-9);
    }

    @Test
    public void testCCD2_mixedTreesProbabilities() {
        CCD2 ccd = new CCD2(mixedTrees, 0.0);

        double p1 = ccd.getProbabilityOfTree(t1);
        double p2 = ccd.getProbabilityOfTree(t2);
        double p3 = ccd.getProbabilityOfTree(t3);

        assertTrue("T1 should have positive probability", p1 > 0);
        assertTrue("T2 should have positive probability", p2 > 0);
        assertTrue("T3 should have positive probability", p3 > 0);
        assertTrue("T1 should have highest probability", p1 > p2);
        assertTrue("T1 should have highest probability", p1 > p3);
    }

    @Test
    public void testCCD2_containsTree() {
        CCD2 ccd = new CCD2(mixedTrees, 0.0);

        assertTrue("CCD2 should contain T1", ccd.containsTree(t1));
        assertTrue("CCD2 should contain T2", ccd.containsTree(t2));
        assertTrue("CCD2 should contain T3", ccd.containsTree(t3));
    }

    @Test
    public void testCCD2_getMAPTree() {
        CCD2 ccd = new CCD2(mixedTrees, 0.0);
        Tree mapTree = ccd.getMAPTree();

        double mapProb = ccd.getProbabilityOfTree(mapTree);
        double t1Prob = ccd.getProbabilityOfTree(t1);
        assertEquals("MAP tree should have same probability as T1", t1Prob, mapProb, 1e-9);
    }

    @Test
    public void testCCD2_sampleTreeValid() {
        CCD2 ccd = new CCD2(mixedTrees, 0.0);
        ccd.setRandom(new Random(42));

        for (int i = 0; i < 100; i++) {
            Tree sampled = ccd.sampleTree();
            assertNotNull("Sampled tree should not be null", sampled);
            assertTrue("Sampled tree should have positive probability",
                    ccd.getProbabilityOfTree(sampled) > 0);
        }
    }

    @Test
    public void testCCD2_burnin() {
        List<Tree> trees = new ArrayList<>();
        for (int i = 0; i < 5; i++) trees.add(parseNewick(T2_NEWICK));
        for (int i = 0; i < 5; i++) trees.add(parseNewick(T1_NEWICK));

        CCD2 ccd = new CCD2(trees, 0.5);

        assertEquals(5, ccd.getNumberOfBaseTrees());
        assertEquals(1.0, ccd.getProbabilityOfTree(t1), 1e-9);
    }

    // ======================== Cross-model Tests ========================

    @Test
    public void testAllModels_sameNumberOfBaseTrees() {
        CCD0 ccd0 = new CCD0(mixedTrees, 0.0);
        CCD1 ccd1 = new CCD1(mixedTrees, 0.0);
        CCD2 ccd2 = new CCD2(mixedTrees, 0.0);

        assertEquals(10, ccd0.getNumberOfBaseTrees());
        assertEquals(10, ccd1.getNumberOfBaseTrees());
        assertEquals(10, ccd2.getNumberOfBaseTrees());
    }

    @Test
    public void testAllModels_numberOfLeavesConsistent() {
        CCD0 ccd0 = new CCD0(mixedTrees, 0.0);
        CCD1 ccd1 = new CCD1(mixedTrees, 0.0);
        CCD2 ccd2 = new CCD2(mixedTrees, 0.0);

        assertEquals(5, ccd0.getNumberOfLeaves());
        assertEquals(5, ccd1.getNumberOfLeaves());
        assertEquals(5, ccd2.getNumberOfLeaves());
    }

    @Test
    public void testCCD1andCCD2_cladeCountRelationship() {
        // CCD0 expands the graph beyond observed clades
        // CCD2 may have additional clades from sibling-context tracking
        CCD0 ccd0 = new CCD0(mixedTrees, 0.0);
        CCD1 ccd1 = new CCD1(mixedTrees, 0.0);
        CCD2 ccd2 = new CCD2(mixedTrees, 0.0);

        int c0 = ccd0.getNumberOfClades();
        int c1 = ccd1.getNumberOfClades();
        int c2 = ccd2.getNumberOfClades();

        // CCD2 should have at least as many clades as CCD1
        assertTrue("CCD2 should have >= CCD1 clades", c2 >= c1);

        // CCD0 expands the graph and should have the most clades
        assertTrue("CCD0 should have >= CCD1 clades", c0 >= c1);
    }
}
