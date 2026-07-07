package test.ccd.algorithms;

import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import ccd.algorithms.NNICladeExpansion.PairingMode;
import ccd.algorithms.regularisation.NNIRegCCDOptimiser;
import ccd.algorithms.regularisation.NNIRegCCDOptimiser.Factor;
import ccd.algorithms.regularisation.NNIRegCCDOptimiser.OptResult;
import ccd.model.NNIRegCCD;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Tests for {@link NNIRegCCDOptimiser}: the closed-form trace scorer must agree
 * with the baked model, and the grid search must run.
 */
public class NNIRegCCDOptimiserTest {

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

    /**
     * The trace-based log probability at (alpha, beta) must equal the baked
     * NNIRegCCD's getLogProbabilityOfTree for the same (alpha, beta), for both
     * covered and uncovered trees, in both pairing modes.
     */
    @Test
    public void testTraceScorerMatchesBakedModel() {
        List<Tree> training = trainingTrees();
        List<Tree> testTrees = new ArrayList<>(training);
        testTrees.add(parse("((((A:1,B:1):1,D:1):1,C:1):1,E:1):0;")); // unsampled clade {A,B,D}

        double[][] params = {{0.4, 0.4}, {0.4, 0.02}, {0.25, 0.1}, {0.1, 0.5}};

        for (PairingMode mode : PairingMode.values()) {
            for (double[] ab : params) {
                double alpha = ab[0], beta = ab[1];
                NNIRegCCD model = new NNIRegCCD(training, 0.0, mode, alpha, beta);
                List<List<Factor>> traces = NNIRegCCDOptimiser.extractTraces(model, testTrees);

                for (int i = 0; i < testTrees.size(); i++) {
                    double modelLP = model.getLogProbabilityOfTree(testTrees.get(i));
                    if (modelLP == Double.NEGATIVE_INFINITY) {
                        assertNull(traces.get(i),
                                "uncovered tree should have a null trace [" + mode + "]");
                    } else {
                        assertEquals(modelLP, NNIRegCCDOptimiser.logProb(traces.get(i), alpha, beta), 1e-9,
                                "trace logP must match baked model [" + mode
                                        + ", a=" + alpha + ", b=" + beta + "]");
                    }
                }
            }
        }
    }

    /** The grid search runs and returns a finite optimum within the grid. */
    @Test
    public void testOptimiseRuns() {
        List<Tree> training = trainingTrees();
        double[] alphas = {0.1, 0.4, 1.0};
        double[] betas = {0.01, 0.1, 0.4};

        OptResult result = NNIRegCCDOptimiser.optimise(training, 2, PairingMode.CO_OCCURRING, alphas, betas);

        assertTrue(Double.isFinite(result.meanLogProb()), "optimum should be finite");
        assertTrue(result.alpha() >= 0.1 && result.alpha() <= 1.0, "alpha within grid");
        assertTrue(result.beta() >= 0.01 && result.beta() <= 0.4, "beta within grid");
        assertTrue(result.covered() > 0, "some held-out trees covered");
    }
}
