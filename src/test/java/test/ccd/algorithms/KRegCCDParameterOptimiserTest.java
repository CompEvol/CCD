package test.ccd.algorithms;

import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import ccd.algorithms.regularisation.KRegCCDParameterOptimiser;
import ccd.algorithms.regularisation.NNIHeldOutComparison.FoldAssignment;
import ccd.model.KRegCCD;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Tests selecting KRegCCD's hyperparameters by held-out tree probability:
 * {@link KRegCCD#getLogProbabilityOfTree(Tree, double)} re-scores at an arbitrary mu, and
 * {@link KRegCCDParameterOptimiser} maximises the cross-validated held-out log-probability over
 * {@code (alpha, mu)}.
 */
public class KRegCCDParameterOptimiserTest {

    private static final List<String> TAXA = Arrays.asList("A", "B", "C", "D", "E");

    private static Tree parse(String newick) {
        return new TreeParser(TAXA, newick, 1, false);
    }

    /** A small, varied posterior-like sample with topological turnover, so cross-validation
     * folds genuinely contain clades the training folds miss. */
    private static List<Tree> trees() {
        List<Tree> trees = new ArrayList<>();
        trees.add(parse("((((A:1,B:1):1,C:1):1,D:1):1,E:1):0;"));
        trees.add(parse("((((A:1,B:1):1,C:1):1,D:1):1,E:1):0;"));
        trees.add(parse("(((A:1,B:1):1,(C:1,D:1):1):1,E:1):0;"));
        trees.add(parse("(((A:1,B:1):1,(C:1,D:1):1):1,E:1):0;"));
        trees.add(parse("((((A:1,C:1):1,B:1):1,D:1):1,E:1):0;"));
        trees.add(parse("((((A:1,B:1):1,D:1):1,C:1):1,E:1):0;"));
        trees.add(parse("(((A:1,B:1):1,C:1):1,(D:1,E:1):1):0;"));
        trees.add(parse("((((A:1,B:1):1,E:1):1,C:1):1,D:1):0;"));
        return trees;
    }

    /**
     * Scoring at an overridden mu on one backbone equals a freshly built NONE-tail model at that
     * mu — this is what lets the optimiser sweep mu without rebuilding.
     */
    @Test
    public void testMuOverrideMatchesFreshNoneModel() {
        List<Tree> trees = trees();
        // backbone built at one mu / tail mode; scoring overrides mu and uses tail = 0
        KRegCCD built = new KRegCCD(trees, 0.0, 0.02, 0.4, 2, KRegCCD.TailMode.BOUND);
        for (double mu : new double[]{0.005, 0.05, 0.2}) {
            KRegCCD fresh = new KRegCCD(trees, 0.0, mu, 0.4, 2, KRegCCD.TailMode.NONE);
            for (Tree t : trees) {
                assertEquals(fresh.getLogProbabilityOfTree(t), built.getLogProbabilityOfTree(t, mu), 1e-9,
                        "override scoring at mu=" + mu + " must match a fresh NONE model");
            }
        }
    }

    /**
     * The optimiser returns valid parameters, its reported objective is the genuine cross-validated
     * held-out log-probability at (alpha*, mu*), and that mu* is at least as good as the extreme mu
     * values at the same alpha* — i.e. it really maximises held-out probability over mu.
     */
    @Test
    public void testOptimiserMaximisesHeldOutOverMu() {
        List<Tree> trees = trees();
        KRegCCDParameterOptimiser.Params p = KRegCCDParameterOptimiser.optimise(trees, 4);

        assertTrue(p.mu() > 0 && p.mu() < 1, "mu* must be a valid escape probability, got " + p.mu());
        assertTrue(p.alpha() > 0, "alpha* must be positive, got " + p.alpha());
        assertTrue(Double.isFinite(p.heldOutLogProb()), "held-out logP must be finite");

        double atOptimum = KRegCCDParameterOptimiser.crossValidatedLogProb(
                trees, 4, FoldAssignment.CONTIGUOUS, p.alpha(), p.mu());
        assertEquals(p.heldOutLogProb(), atOptimum, 1e-6, "reported objective must be reproducible");

        // compare against the two ends of the searched (valid) mu range
        double atTinyMu = KRegCCDParameterOptimiser.crossValidatedLogProb(
                trees, 4, FoldAssignment.CONTIGUOUS, p.alpha(), 1e-3);
        double atCeilingMu = KRegCCDParameterOptimiser.crossValidatedLogProb(
                trees, 4, FoldAssignment.CONTIGUOUS, p.alpha(), KRegCCD.MU_RELIABLE_MAX);
        assertTrue(atOptimum >= atTinyMu - 1e-9,
                "mu* should be at least as good as a near-zero mu (" + atOptimum + " vs " + atTinyMu + ")");
        assertTrue(atOptimum >= atCeilingMu - 1e-9,
                "mu* should be at least as good as the ceiling mu (" + atOptimum + " vs " + atCeilingMu + ")");
        assertTrue(p.mu() <= KRegCCD.MU_RELIABLE_MAX + 1e-9,
                "mu* must not exceed the reliability ceiling, got " + p.mu());
    }

    /**
     * Held-out probability genuinely prefers some escape mass here: an interior mu beats mu -> 0,
     * because the cross-validation folds contain clades their training folds miss (which a
     * zero-escape model would assign probability zero / -inf).
     */
    @Test
    public void testHeldOutPrefersNonzeroEscape() {
        List<Tree> trees = trees();
        double alpha = 0.4;
        double atSmall = KRegCCDParameterOptimiser.crossValidatedLogProb(
                trees, 4, FoldAssignment.STRIDED, alpha, 1e-3);
        double atModerate = KRegCCDParameterOptimiser.crossValidatedLogProb(
                trees, 4, FoldAssignment.STRIDED, alpha, 0.05);
        assertTrue(atModerate > atSmall,
                "moderate escape should beat near-zero escape when held-out trees have novel clades ("
                        + atModerate + " vs " + atSmall + ")");
    }

    /** The factory builds a usable model at the optimised parameters: training trees score finite
     * and it can sample. */
    @Test
    public void testFactoryBuildsUsableModel() {
        List<Tree> trees = trees();
        KRegCCD ccd = KRegCCD.withOptimisedParameters(trees, 4);
        for (Tree t : trees) {
            assertTrue(Double.isFinite(ccd.getLogProbabilityOfTree(t)),
                    "training trees must score finite under the optimised model");
        }
        ccd.setRandom(new java.util.Random(1));
        Tree sampled = ccd.sampleTree();
        assertEquals(TAXA.size(), sampled.getLeafNodeCount(), "sampled tree must have all taxa");
    }
}
