package test.ccd.model;

import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import ccd.algorithms.NNICladeExpansion.PairingMode;
import ccd.model.Clade;
import ccd.model.CladePartition;
import ccd.model.NNIRegCCD;
import ccd.model.RegCCD;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Tests that {@link NNIRegCCD} assigns nonzero probability to held-out trees
 * containing an unsampled but NNI-reachable clade --- the case plain regCCD
 * cannot cover --- while remaining a valid distribution.
 */
public class NNIRegCCDTest {

    private static final List<String> TAXA = Arrays.asList("A", "B", "C", "D", "E");

    private static Tree parse(String newick) {
        return new TreeParser(TAXA, newick, 1, false);
    }

    private static List<Tree> trainingTrees() {
        List<Tree> trees = new ArrayList<>();
        trees.add(parse("((((A:1,B:1):1,C:1):1,D:1):1,E:1):0;"));
        trees.add(parse("((((A:1,C:1):1,B:1):1,E:1):1,D:1):0;"));
        return trees;
    }

    // one NNI from training tree 1: clade {A,B,D} is not in the sample
    private static Tree heldOutWithUnsampledClade() {
        return parse("((((A:1,B:1):1,D:1):1,C:1):1,E:1):0;");
    }

    /**
     * The held-out tree contains the clade {A,B,D}, which no training tree has.
     * Plain regCCD gives it zero probability; NNIRegCCD (either pairing mode)
     * gives it nonzero probability because the NNI clade expansion adds {A,B,D}
     * and the split {A,B,C,D} -> {{A,B,D}, C}.
     */
    @Test
    public void testHeldOutUnsampledCladeBecomesNonzero() {
        List<Tree> training = trainingTrees();
        Tree heldOut = heldOutWithUnsampledClade();

        RegCCD reg = new RegCCD(training, 0.0);
        double pReg = reg.getProbabilityOfTree(heldOut);
        assertEquals(0.0, pReg, 0.0, "plain regCCD cannot cover the unsampled clade");

        for (PairingMode mode : PairingMode.values()) {
            NNIRegCCD nni = new NNIRegCCD(training, 0.0, mode);
            double p = nni.getProbabilityOfTree(heldOut);
            assertTrue(p > 0.0, "NNIRegCCD[" + mode + "] should cover the unsampled clade, got " + p);
            assertTrue(p <= 1.0, "probability must be <= 1, got " + p);
            assertTrue(nni.getNumberOfNNIClades() > 0, "expansion should add clades");
        }
    }

    /** The training trees themselves must keep nonzero probability. */
    @Test
    public void testTrainingTreesRemainSupported() {
        List<Tree> training = trainingTrees();
        for (PairingMode mode : PairingMode.values()) {
            NNIRegCCD nni = new NNIRegCCD(training, 0.0, mode);
            for (Tree t : training) {
                assertTrue(nni.getProbabilityOfTree(t) > 0.0,
                        "training tree should stay supported under NNIRegCCD[" + mode + "]");
            }
        }
    }

    /**
     * NNIRegCCD must remain a valid CCD: the conditional clade probabilities of
     * every non-leaf clade's splits sum to 1.
     */
    @Test
    public void testConditionalDistributionsNormalised() {
        List<Tree> training = trainingTrees();
        for (PairingMode mode : PairingMode.values()) {
            List<NNIRegCCD> models = List.of(
                    new NNIRegCCD(training, 0.0, mode),                 // single alpha
                    new NNIRegCCD(training, 0.0, mode, 0.4, 0.02));     // two-level alpha/beta
            for (NNIRegCCD nni : models) {
                for (Clade clade : nni.getClades()) {
                    if (clade.isLeaf()) {
                        continue;
                    }
                    double sum = 0.0;
                    for (CladePartition partition : clade.getPartitions()) {
                        sum += partition.getCCP();
                    }
                    assertEquals(1.0, sum, 1e-9,
                            "CCPs of clade " + clade + " should sum to 1 under NNIRegCCD[" + mode + "]");
                }
            }
        }
    }

    /**
     * Two-level smoothing with beta == alpha must reproduce single-alpha
     * smoothing exactly (the AdditiveXY path collapses to AdditiveX).
     */
    @Test
    public void testTwoLevelEqualsSingleAlphaWhenBetaEqualsAlpha() {
        List<Tree> training = trainingTrees();
        Tree heldOut = heldOutWithUnsampledClade();
        double alpha = 0.4;

        for (PairingMode mode : PairingMode.values()) {
            NNIRegCCD single = new NNIRegCCD(training, 0.0, mode, alpha);
            NNIRegCCD two = new NNIRegCCD(training, 0.0, mode, alpha, alpha);

            assertEquals(single.getProbabilityOfTree(heldOut), two.getProbabilityOfTree(heldOut), 1e-12,
                    "held-out prob should match for beta == alpha [" + mode + "]");
            for (Tree t : training) {
                assertEquals(single.getProbabilityOfTree(t), two.getProbabilityOfTree(t), 1e-12,
                        "training prob should match for beta == alpha [" + mode + "]");
            }
        }
    }

    /**
     * A smaller beta gives NNI-derived splits less mass: the held-out tree that
     * relies on an NNI parent split loses probability, while a training tree
     * built from observed splits at the same parent gains probability. Both stay
     * positive.
     */
    @Test
    public void testSmallerBetaShiftsMassToObservedSplits() {
        List<Tree> training = trainingTrees();
        Tree t1 = training.get(0);
        Tree heldOut = heldOutWithUnsampledClade();
        double alpha = 0.4;

        NNIRegCCD big = new NNIRegCCD(training, 0.0, PairingMode.CO_OCCURRING, alpha, alpha);
        NNIRegCCD small = new NNIRegCCD(training, 0.0, PairingMode.CO_OCCURRING, alpha, 0.01);

        double pHeldBig = big.getProbabilityOfTree(heldOut);
        double pHeldSmall = small.getProbabilityOfTree(heldOut);
        assertTrue(pHeldSmall > 0, "held-out tree still covered with small beta");
        assertTrue(pHeldSmall < pHeldBig, "smaller beta gives the NNI-reliant tree less mass");

        double pTrainBig = big.getProbabilityOfTree(t1);
        double pTrainSmall = small.getProbabilityOfTree(t1);
        assertTrue(pTrainSmall > pTrainBig, "smaller beta gives the observed-split tree more mass");
    }

    /**
     * Graph-wide pairing should cover at least as many clades as co-occurring,
     * since the strict NNI neighbourhood is contained in the graph-wide one.
     */
    @Test
    public void testGraphWideCoversAtLeastCoOccurring() {
        List<Tree> training = trainingTrees();
        NNIRegCCD broad = new NNIRegCCD(training, 0.0, PairingMode.GRAPH_WIDE);
        NNIRegCCD strict = new NNIRegCCD(training, 0.0, PairingMode.CO_OCCURRING);
        assertTrue(broad.getNumberOfNNIClades() >= strict.getNumberOfNNIClades(),
                "graph-wide should add at least as many clades as co-occurring");
    }
}
