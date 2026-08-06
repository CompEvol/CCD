package ccd.algorithms.regularisation;

import beast.base.evolution.tree.Tree;
import ccd.algorithms.regularisation.NNIHeldOutComparison.FoldAssignment;
import ccd.model.KRegCCD;
import org.apache.commons.math4.legacy.analysis.UnivariateFunction;
import org.apache.commons.math4.legacy.optim.MaxEval;
import org.apache.commons.math4.legacy.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math4.legacy.optim.univariate.BrentOptimizer;
import org.apache.commons.math4.legacy.optim.univariate.SearchInterval;
import org.apache.commons.math4.legacy.optim.univariate.UnivariateObjectiveFunction;
import org.apache.commons.math4.legacy.optim.univariate.UnivariatePointValuePair;

import java.util.ArrayList;
import java.util.List;

/**
 * Selects {@link KRegCCD}'s hyperparameters by maximising <em>cross-validated held-out</em>
 * tree log-probability — the counterpart of {@link RegCCDParameterOptimiser} (which tunes the
 * backbone smoothing {@code alpha} alone) for the full-support model's two parameters.
 *
 * <p>
 * The objective is the total log-probability of held-out trees under batched leave-one-tree-out
 * (B-fold cross-validation): each fold is held out in turn while the backbone is rebuilt on the
 * rest and scores the fold. Scoring training trees instead would drive {@code mu -> 0} (the
 * backbone fits them exactly, so any reserved escape mass is pure loss), so the trees scored must
 * be out of training.
 *
 * <p>
 * The two parameters are optimised jointly but with different search strategies that exploit their
 * structure:
 * <ul>
 *   <li>{@code alpha} (backbone smoothing) changes the clade/split graph, so each value needs a
 *       full rebuild per fold; it has a sharp interior optimum and is found by Brent's method
 *       (as in {@link RegCCDParameterOptimiser}).</li>
 *   <li>{@code mu} (escape probability) does not touch the counts — only {@code eps(C)} re-solves
 *       from the cached {@code N_j} — so for a fixed {@code alpha} the whole {@code mu} sweep is
 *       evaluated on one set of trained backbones via
 *       {@link KRegCCD#getLogProbabilityOfTree(Tree, double)}. Its held-out curve is shallow, so
 *       it is searched on a (log-spaced) grid and taken at the argmax rather than by a
 *       derivative-based method that would wander on the flat ridge.</li>
 * </ul>
 * For each candidate {@code alpha}, {@code mu} is profiled out (best over the grid); Brent then
 * maximises that profiled objective over {@code alpha}. When {@link #FIXED_ALPHA} is set, the Brent
 * search is skipped and only {@code mu} is profiled at that pinned {@code alpha} — a large speedup
 * (one set of per-fold backbones instead of one per alpha evaluation).
 *
 * @author Claude (CCD-Sophie)
 */
public class KRegCCDParameterOptimiser {

    /** The selected hyperparameters and the cross-validated held-out log-probability at them. */
    public record Params(double alpha, double mu, double heldOutLogProb) {
    }

    /** Default number of cross-validation folds. */
    public static final int DEFAULT_FOLDS = 5;

    /** Default reserve depth for the candidate models (matches {@link KRegCCD#DEFAULT_RESERVE_DEPTH}). */
    private static final int RESERVE_DEPTH = KRegCCD.DEFAULT_RESERVE_DEPTH;

    /* alpha search interval (strictly positive; mirrors RegCCDParameterOptimiser's range). */
    private static final double ALPHA_LO = 1e-3;
    private static final double ALPHA_HI = 2.0;
    private static final double ALPHA_START = KRegCCD.DEFAULT_ALPHA;

    /* Fixed backbone-smoothing alpha. The Yule50 sweep found the fitted alpha clustered tightly
     * around 0.4 (per-size medians 0.46/0.40/0.36) and the held-out objective is shallow in alpha,
     * so the Brent search buys little but costs a full backbone rebuild per alpha evaluation (up to
     * MAX_ALPHA_EVALS x folds rebuilds). Pinning alpha = 0.4 and profiling only mu collapses that to
     * one set of per-fold backbones -- the dominant speedup for the 100-rep re-run. Set to null to
     * restore the Brent search over [ALPHA_LO, ALPHA_HI]. */
    private static final Double FIXED_ALPHA = 0.4;

    /* mu grid: log-spaced from negligible escape up to the reliability ceiling. The upper bound is
     * KRegCCD.MU_RELIABLE_MAX, not larger: above it the depth-k reserve approximation breaks down
     * (the omitted tail grows ~mu^3), and because the held-out objective is scored with tail = 0
     * (NONE) that breakdown inflates the score ~mu^3 -- which would bias the search toward large
     * mu. Searching beyond the ceiling is therefore both invalid and misleading.
     *
     * The lower bound is 1e-4 (was 1e-3): on large posterior samples the held-out optimum routinely
     * sat on the old 1e-3 floor (once most clades are seen the data wants very little escape mass),
     * so the floor was binding rather than interior. MU_GRID is widened from 41 to 61 to keep the
     * per-decade grid density roughly constant over the now-wider range. A floor warning (below)
     * flags fits that still pin at the lowered floor. */
    private static final double MU_LO = 1e-4;
    private static final double MU_HI = KRegCCD.MU_RELIABLE_MAX;
    private static final int MU_GRID = 61;

    private static final int MAX_ALPHA_EVALS = 60;

    public static Params optimise(List<Tree> trees) {
        return optimise(trees, DEFAULT_FOLDS);
    }

    public static Params optimise(List<Tree> trees, int folds) {
        // CONTIGUOUS (block) folds, not STRIDED: posterior trees are autocorrelated, so strided folds
        // leave each held-out tree's near-neighbours in its own training fold, making the CV optimistic
        // and selecting under-regularised (too-small alpha/mu) parameters. A contiguous held-out block
        // is far less correlated with its training set, giving an honest held-out objective.
        return optimise(trees, folds, FoldAssignment.CONTIGUOUS);
    }

    public static Params optimise(List<Tree> trees, int folds, FoldAssignment assignment) {
        return optimise(trees, folds, assignment, RESERVE_DEPTH);
    }

    /**
     * As {@link #optimise(List, int, FoldAssignment)} but with an explicit reserve depth {@code k}.
     * The per-clade reserve solve dominates cross-validation on large clade sets, and almost all of
     * it is the boundary-4 (N_2) match that {@code k = 1} skips outright -- on a 441-taxon set the
     * default k = 2 did not finish a single fold in 11 h of CPU. k = 1 solves eps from N_1 alone, which
     * shifts eps by a few percent, so prefer the default unless the search is otherwise intractable.
     */
    public static Params optimise(List<Tree> trees, int folds, FoldAssignment assignment,
                                  int reserveDepth) {
        int n = trees.size();
        if (n < 2) {
            throw new IllegalArgumentException("need at least 2 trees to cross-validate, got " + n);
        }
        int b = Math.max(2, Math.min(folds, n));

        // Fixed fold partition (only the backbones change with alpha, not the split).
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

        double[] muGrid = logspace(MU_LO, MU_HI, MU_GRID);

        double alphaStar;
        if (FIXED_ALPHA != null) {
            // alpha pinned (see FIXED_ALPHA): skip the Brent search and profile only mu.
            alphaStar = FIXED_ALPHA;
        } else {
            UnivariateFunction profiled = alpha -> {
                double[] best = bestMuLogProb(trainByFold, testByFold, alpha, muGrid, reserveDepth);
                System.out.println(String.format(
                        "  testing alpha = %.5f -> best mu = %.5f, held-out logP = %.5f",
                        alpha, best[0], best[1]));
                return best[1];
            };

            BrentOptimizer optimiser = new BrentOptimizer(1e-4, 1e-4);
            UnivariatePointValuePair sol = optimiser.optimize(
                    new MaxEval(MAX_ALPHA_EVALS),
                    new UnivariateObjectiveFunction(profiled),
                    GoalType.MAXIMIZE,
                    new SearchInterval(ALPHA_LO, ALPHA_HI, ALPHA_START));
            alphaStar = sol.getPoint();
        }

        double[] best = bestMuLogProb(trainByFold, testByFold, alphaStar, muGrid, reserveDepth);
        System.out.println(String.format(
                "KRegCCD optimised: alpha = %.5f, mu = %.5f, held-out logP = %.5f",
                alphaStar, best[0], best[1]));
        if (best[0] >= MU_HI * (1 - 1e-9)) {
            System.err.println("WARNING: held-out probability is still increasing in mu at the "
                    + "reliability ceiling mu = " + MU_HI + " (KRegCCD.MU_RELIABLE_MAX); the data "
                    + "wants more escape mass than the depth-k approximation can validly provide. "
                    + "Treat mu = " + MU_HI + " as a capped estimate, not an interior optimum.");
        }
        if (best[0] <= MU_LO * (1 + 1e-9)) {
            System.err.println("WARNING: held-out probability is still increasing as mu decreases at "
                    + "the search floor mu = " + MU_LO + "; the data wants even less escape mass than "
                    + "this. Treat mu = " + MU_LO + " as a capped (non-interior) estimate.");
        }
        return new Params(alphaStar, best[0], best[1]);
    }

    /**
     * Optimise mu alone, given an alpha value
     */
    public static MuResult optimiseMu(List<Tree> trees, int folds, FoldAssignment assignment, double alpha) {
        int n = trees.size();
        if (n < 2) throw new IllegalArgumentException("need at least 2 trees to cross-validate, got " + n);
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
        double[] muGrid = logspace(MU_LO, MU_HI, MU_GRID);
        double[] best = bestMuLogProb(trainByFold, testByFold, alpha, muGrid, RESERVE_DEPTH);
        return new MuResult(best[0], best[1]);
    }

    public record MuResult(double mu, double heldOutLogProb) {
    }

    /**
     * The cross-validated held-out total log-probability at a fixed {@code (alpha, mu)} — the
     * objective the optimiser maximises. Rebuilds the per-fold backbones, so prefer the internal
     * grid path for a full search; this is for evaluating or comparing individual points.
     */
    public static double crossValidatedLogProb(List<Tree> trees, int folds, FoldAssignment assignment,
                                               double alpha, double mu) {
        int n = trees.size();
        int b = Math.max(2, Math.min(folds, n));
        double total = 0.0;
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
            KRegCCD model = new KRegCCD(train, 0.0, KRegCCD.DEFAULT_MU, alpha,
                    RESERVE_DEPTH, KRegCCD.TailMode.NONE);
            for (Tree t : test) {
                total += model.getLogProbabilityOfTree(t, mu);
            }
        }
        return total;
    }

    /**
     * For a fixed {@code alpha}: builds one NONE-tail backbone per fold and returns
     * {@code {argmax mu over the grid, that held-out total log-probability}}. All grid values reuse
     * the same backbones — {@code mu} only re-solves {@code eps} from the cached counts.
     */
    private static double[] bestMuLogProb(List<List<Tree>> trainByFold, List<List<Tree>> testByFold,
                                          double alpha, double[] muGrid, int reserveDepth) {
        int b = trainByFold.size();
        KRegCCD[] models = new KRegCCD[b];
        for (int f = 0; f < b; f++) {
            // mu here is a placeholder: scoring overrides it via getLogProbabilityOfTree(tree, mu).
            models[f] = new KRegCCD(trainByFold.get(f), 0.0, KRegCCD.DEFAULT_MU, alpha,
                    reserveDepth, KRegCCD.TailMode.NONE);
            // Solve the reserves up front, in parallel. Scoring the held-out trees would otherwise
            // fill regCache lazily on one thread during the first mu of the grid, and that solve --
            // not the scoring arithmetic -- is what dominates on a large clade set. Same values
            // either way (a clade's reserve depends only on its own observed subclades, and the
            // N_j counts are mu-independent), so this changes cost only.
            models[f].precomputeReserves();
        }
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
        return new double[]{bestMu, bestLogProb};
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
