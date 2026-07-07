package ccd.algorithms.regularisation;

import beast.base.evolution.tree.Tree;
import ccd.algorithms.regularisation.NNIHeldOutComparison.FoldAssignment;
import ccd.model.MRegCCD;

import java.util.ArrayList;
import java.util.List;

/**
 * Selects the one hyperparameter {@code mu} of {@link MRegCCD} by maximising <em>cross-validated
 * held-out</em> tree log-probability — the single-parameter counterpart of
 * {@link KRegCCDParameterOptimiser}.
 *
 * <p>The objective is the total log-probability of held-out trees under {@code b}-fold
 * cross-validation: each fold is held out in turn while the backbone is rebuilt on the rest and
 * scores the fold. Scoring training trees instead would drive {@code mu -> 0} (the backbone fits them
 * exactly, so any reserved escape mass is pure loss), so the scored trees must be out of training.
 *
 * <p>Folds are <strong>CONTIGUOUS</strong> blocks by default, not strided: posterior trees are
 * autocorrelated, so strided folds leave each held-out tree's near-neighbours in its own training
 * fold and make the CV optimistic. A contiguous held-out block is far less correlated with its
 * training set.
 *
 * <p>{@code mu} does not touch the clade counts — only {@code eps(C)} re-solves — so the whole grid is
 * evaluated on one set of per-fold backbones via {@link MRegCCD#getLogProbabilityOfTree(Tree, double)}.
 * The models are built with the tail correction ON: unlike the un-tailed score the tail keeps each
 * candidate {@code mu} properly normalised, so the held-out objective is honest (a tail-free score
 * super-normalises and would bias the search toward large {@code mu}).
 *
 * @author Claude (CCD-Sophie)
 */
public class MRegCCDParameterOptimiser {

    /** The selected escape probability and the cross-validated held-out log-probability at it. */
    public record MuResult(double mu, double heldOutLogProb) {
    }

    /** Default number of cross-validation folds. */
    public static final int DEFAULT_FOLDS = 5;

    /* mu grid: log-spaced. The lower bound is negligible escape; the upper bound is generous (the
     * per-new-split model spends more escape mass per region than KRegCCD, so its optimum sits
     * higher). MRegCCD self-normalises via the geometric tail at every mu, so there is no hard
     * reliability ceiling as in KRegCCD's tail-free search; floor/ceiling warnings still flag a
     * non-interior optimum. */
    private static final double MU_LO = 1e-4;
    private static final double MU_HI = 0.2;
    private static final int MU_GRID = 61;

    public static MuResult optimiseMu(List<Tree> trees) {
        return optimiseMu(trees, DEFAULT_FOLDS, FoldAssignment.CONTIGUOUS);
    }

    public static MuResult optimiseMu(List<Tree> trees, int folds, FoldAssignment assignment) {
        int n = trees.size();
        if (n < 2) {
            throw new IllegalArgumentException("need at least 2 trees to cross-validate, got " + n);
        }
        int b = Math.max(2, Math.min(folds, n));

        List<List<Tree>> trainByFold = new ArrayList<>(b);
        List<List<Tree>> testByFold = new ArrayList<>(b);
        for (int f = 0; f < b; f++) {
            List<Tree> train = new ArrayList<>();
            List<Tree> test = new ArrayList<>();
            for (int i = 0; i < n; i++) {
                boolean isTest = switch (assignment) {
                    case STRIDED -> (i % b == f);
                    case CONTIGUOUS -> (i >= (int) ((long) f * n / b)) && (i < (int) ((long) (f + 1) * n / b));
                };
                (isTest ? test : train).add(trees.get(i));
            }
            trainByFold.add(train);
            testByFold.add(test);
        }

        // One backbone per fold (tail ON so each mu is properly normalised); reused across the grid.
        MRegCCD[] models = new MRegCCD[b];
        for (int f = 0; f < b; f++) {
            models[f] = new MRegCCD(trainByFold.get(f), 0.0, MRegCCD.DEFAULT_MU,
                    MRegCCD.DEFAULT_RESERVE_DEPTH, true);
        }

        double[] muGrid = logspace(MU_LO, MU_HI, MU_GRID);
        double bestMu = muGrid[0];
        double bestLogProb = Double.NEGATIVE_INFINITY;
        for (double mu : muGrid) {
            double total = 0.0;
            for (int f = 0; f < b; f++) {
                for (Tree t : testByFold.get(f)) {
                    total += models[f].getLogProbabilityOfTree(t, mu);
                }
            }
            if (total > bestLogProb) {
                bestLogProb = total;
                bestMu = mu;
            }
        }
        if (bestMu >= MU_HI * (1 - 1e-9)) {
            System.err.println("WARNING: MRegCCD held-out logP still rising at the mu ceiling "
                    + MU_HI + "; treat as a capped estimate.");
        }
        if (bestMu <= MU_LO * (1 + 1e-9)) {
            System.err.println("WARNING: MRegCCD held-out logP still rising as mu falls to the floor "
                    + MU_LO + "; treat as a capped estimate.");
        }
        return new MuResult(bestMu, bestLogProb);
    }

    private static double[] logspace(double lo, double hi, int n) {
        double[] grid = new double[n];
        double logLo = Math.log(lo);
        double logHi = Math.log(hi);
        for (int i = 0; i < n; i++) {
            grid[i] = Math.exp(logLo + (logHi - logLo) * i / (n - 1));
        }
        return grid;
    }
}
