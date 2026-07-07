package ccd.algorithms.regularisation;

import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator;
import ccd.algorithms.NNICladeExpansion.PairingMode;
import ccd.model.KRegCCD;
import ccd.model.NNIRegCCD;
import ccd.model.RegCCD;
import ccd.tools.CCDToolUtil;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;

/**
 * Compares the out-of-sample coverage of regCCD and NNI-expanded regCCD: how many
 * held-out trees receive nonzero probability (i.e. are <em>covered</em>) and the
 * mean log probability over the covered trees.
 *
 * <p>
 * This targets the residual regCCD limitation that trees containing an unsampled
 * clade receive zero probability. The NNI clade expansion adds plausible unseen
 * clades, so it should reduce the zero-probability count. Scoring uses the
 * regularised model directly ({@code getLogProbabilityOfTree}); held-out trees
 * are genuinely out-of-training, so no leave-one-out count adjustment is needed.
 * </p>
 *
 * <p>
 * Two evaluation protocols are provided:
 * <ul>
 * <li>{@link #compare} --- train on one set, score an independent held-out set
 * (e.g.\ a separate chain, or one file split in two);</li>
 * <li>{@link #crossValidate} --- <em>batched leave-one-tree-out</em>: B-fold
 * cross-validation that rebuilds the (NNI-)regCCD on each training fold and
 * scores the held-out fold. This approximates true leave-one-tree-out with only
 * B rebuilds instead of N. A full rebuild per fold is required because the NNI
 * clade set depends on which trees are present, so the memoised count-subtraction
 * shortcut used for plain regCCD does not apply to the NNI model.</li>
 * </ul>
 * </p>
 *
 * @author Claude (CCD-Sophie)
 */
public class NNIHeldOutComparison {

    /**
     * Coverage / mean-log-prob summary for one model on a held-out set.
     *
     * @param model              model description
     * @param total              number of held-out trees
     * @param covered            number assigned nonzero probability
     * @param meanLogProbCovered mean log prob over this model's covered trees
     * @param meanLogProbCommon  mean log prob over trees covered by all models
     *                           (the fair, same-denominator sharpness metric)
     */
    public record Result(String model, int total, int covered,
                         double meanLogProbCovered, double meanLogProbCommon) {
        public double coverageFraction() {
            return (total == 0) ? 0 : covered / (double) total;
        }

        @Override
        public String toString() {
            return String.format(
                    "%-24s covered %5d / %-5d (%5.1f%%)   logP/own = %9.4f   logP/common = %9.4f",
                    model, covered, total, 100 * coverageFraction(),
                    meanLogProbCovered, meanLogProbCommon);
        }
    }

    /** A named recipe for building a (Reg)CCD from a set of training trees. */
    private record ModelSpec(String label, Function<List<Tree>, RegCCD> builder) {
    }

    /**
     * How trees are assigned to cross-validation folds.
     * <ul>
     * <li>{@link #STRIDED}: tree {@code i} goes to fold {@code i mod B}. With an
     * autocorrelated chain this keeps each held-out tree's neighbours in
     * training, so it tracks true leave-one-tree-out (the default).</li>
     * <li>{@link #CONTIGUOUS}: fold {@code f} is a contiguous block. This removes
     * whole correlated stretches from training, giving a more pessimistic
     * (blocked-CV) estimate that does <em>not</em> approximate LOTO.</li>
     * </ul>
     */
    public enum FoldAssignment {STRIDED, CONTIGUOUS}

    /** The models compared: baseline regCCD, NNI (both modes), and a co-occurring beta sweep. */
    private static List<ModelSpec> modelSpecs(double alpha) {
        List<ModelSpec> specs = new ArrayList<>();
        specs.add(new ModelSpec("regCCD", tr -> new RegCCD(tr, 0.0, alpha)));
        for (PairingMode mode : PairingMode.values()) {
            specs.add(new ModelSpec("NNIRegCCD[" + mode + "]",
                    tr -> new NNIRegCCD(tr, 0.0, mode, alpha)));
        }
        for (double beta : new double[]{alpha / 4, alpha / 20}) {
            specs.add(new ModelSpec(String.format("NNIRegCCD[CO_OCCURRING,b=%.3f]", beta),
                    tr -> new NNIRegCCD(tr, 0.0, PairingMode.CO_OCCURRING, alpha, beta)));
        }
        // full-support KRegCCD (reserve depth 2); compare the three tail modes to
        // show the normalisation effect on mean logP: NONE (super-normalised),
        // BOUND (upper bound, sub-normalised), SAMPLED (Knuth, near-exact)
        for (double mu : new double[]{0.01, 0.05, 0.1, 0.2}) {
            specs.add(new ModelSpec(String.format("KRegCCD[rd=2,none,mu=%.2f]", mu),
                    tr -> new KRegCCD(tr, 0.0, mu, alpha, 2, KRegCCD.TailMode.NONE)));
            specs.add(new ModelSpec(String.format("KRegCCD[rd=2,bound,mu=%.2f]", mu),
                    tr -> new KRegCCD(tr, 0.0, mu, alpha, 2, KRegCCD.TailMode.BOUND)));
            specs.add(new ModelSpec(String.format("KRegCCD[rd=2,sampled,mu=%.2f]", mu),
                    tr -> new KRegCCD(tr, 0.0, mu, alpha, 2, KRegCCD.TailMode.SAMPLED)));
        }
        return specs;
    }

    /* Per-tree log probabilities aligned with the given list. */
    private static double[] logProbs(RegCCD ccd, List<Tree> trees) {
        double[] lp = new double[trees.size()];
        for (int i = 0; i < trees.size(); i++) {
            lp[i] = ccd.getLogProbabilityOfTree(trees.get(i));
        }
        return lp;
    }

    /* Turn per-model held-out log probabilities into coverage / sharpness results. */
    private static List<Result> aggregate(List<String> labels, List<double[]> logs) {
        int n = logs.isEmpty() ? 0 : logs.get(0).length;

        boolean[] common = new boolean[n];
        for (int i = 0; i < n; i++) {
            common[i] = true;
            for (double[] lp : logs) {
                if (lp[i] == Double.NEGATIVE_INFINITY) {
                    common[i] = false;
                    break;
                }
            }
        }

        List<Result> results = new ArrayList<>();
        for (int m = 0; m < labels.size(); m++) {
            double[] lp = logs.get(m);
            int covered = 0, commonCount = 0;
            double sumOwn = 0.0, sumCommon = 0.0;
            for (int i = 0; i < lp.length; i++) {
                if (lp[i] != Double.NEGATIVE_INFINITY) {
                    covered++;
                    sumOwn += lp[i];
                }
                if (common[i]) {
                    sumCommon += lp[i];
                    commonCount++;
                }
            }
            double meanOwn = (covered == 0) ? Double.NaN : sumOwn / covered;
            double meanCommon = (commonCount == 0) ? Double.NaN : sumCommon / commonCount;
            results.add(new Result(labels.get(m), n, covered, meanOwn, meanCommon));
        }
        return results;
    }

    /**
     * Trains each model on {@code training} and scores {@code heldOut}.
     *
     * @param training training trees
     * @param heldOut  independent held-out trees
     * @param alpha    additive-smoothing pseudocount used for all models
     * @return one report line per model
     */
    public static List<Result> compare(List<Tree> training, List<Tree> heldOut, double alpha) {
        List<ModelSpec> specs = modelSpecs(alpha);
        List<String> labels = new ArrayList<>();
        List<double[]> logs = new ArrayList<>();
        for (ModelSpec spec : specs) {
            labels.add(spec.label());
            logs.add(logProbs(spec.builder().apply(training), heldOut));
        }
        return aggregate(labels, logs);
    }

    /**
     * Batched leave-one-tree-out using {@link FoldAssignment#STRIDED} folds.
     *
     * @param trees all trees (post burn-in)
     * @param folds number of cross-validation folds (B)
     * @param alpha additive-smoothing pseudocount used for all models
     * @return one aggregated report line per model
     */
    public static List<Result> crossValidate(List<Tree> trees, int folds, double alpha) {
        return crossValidate(trees, folds, alpha, FoldAssignment.STRIDED);
    }

    /**
     * Batched leave-one-tree-out: B-fold cross-validation. Each fold is held out
     * in turn while every model is rebuilt on the remaining trees and scores the
     * held-out fold, so each tree is scored exactly once, out-of-training. A full
     * rebuild per fold is required because the NNI clade set depends on which
     * trees are present.
     *
     * <p>
     * To approximate true leave-one-tree-out use {@link FoldAssignment#STRIDED}
     * (the default): with an autocorrelated sample it keeps each held-out tree's
     * neighbours in training, as LOTO does, and converges to the LOTO estimate at
     * much smaller {@code folds} than contiguous blocks. Either way the estimate
     * converges to LOTO as {@code folds -> N}.
     * </p>
     *
     * @param trees      all trees (post burn-in)
     * @param folds      number of cross-validation folds (B)
     * @param alpha      additive-smoothing pseudocount used for all models
     * @param assignment how trees are assigned to folds
     * @return one aggregated report line per model
     */
    public static List<Result> crossValidate(List<Tree> trees, int folds, double alpha,
                                             FoldAssignment assignment) {
        int n = trees.size();
        List<ModelSpec> specs = modelSpecs(alpha);
        double[][] logs = new double[specs.size()][n];

        for (int f = 0; f < folds; f++) {
            List<Tree> train = new ArrayList<>(n);
            List<Integer> testIdx = new ArrayList<>();
            for (int i = 0; i < n; i++) {
                boolean isTest = switch (assignment) {
                    case STRIDED -> (i % folds == f);
                    case CONTIGUOUS -> (i >= (int) ((long) f * n / folds))
                            && (i < (int) ((long) (f + 1) * n / folds));
                };
                if (isTest) {
                    testIdx.add(i);
                } else {
                    train.add(trees.get(i));
                }
            }

            for (int m = 0; m < specs.size(); m++) {
                RegCCD model = specs.get(m).builder().apply(train);
                for (int i : testIdx) {
                    logs[m][i] = model.getLogProbabilityOfTree(trees.get(i));
                }
            }
            System.out.println("  fold " + (f + 1) + "/" + folds + " done ("
                    + testIdx.size() + " held out, " + train.size() + " train)");
        }

        List<String> labels = new ArrayList<>();
        List<double[]> logList = new ArrayList<>();
        for (int m = 0; m < specs.size(); m++) {
            labels.add(specs.get(m).label());
            logList.add(logs[m]);
        }
        return aggregate(labels, logList);
    }

    private static List<Tree> loadTrees(String path, int burninPercentage) throws Exception {
        TreeAnnotator.MemoryFriendlyTreeSet treeSet = CCDToolUtil.getTreeSet(path, burninPercentage);
        List<Tree> trees = new ArrayList<>();
        treeSet.reset();
        while (treeSet.hasNext()) {
            trees.add(treeSet.next());
        }
        return trees;
    }

    /**
     * CLI. Forms:
     * <ul>
     * <li>{@code <train.trees> <heldout.trees> [burnin%] [alpha]} --- train on
     * one file, evaluate on an independent file;</li>
     * <li>{@code <all.trees> [burnin%] [alpha]} --- split one file into first
     * (train) and second (held-out) half;</li>
     * <li>{@code cv <all.trees> <folds> [burnin%] [alpha]} --- batched
     * leave-one-tree-out (B-fold cross-validation).</li>
     * </ul>
     *
     * @param args see above
     * @throws Exception if a tree file cannot be read
     */
    public static void main(String[] args) throws Exception {
        if (args.length >= 1 && args[0].equals("calib")) {
            // calib <train.trees> [<heldout.trees>] [burnin%] [alpha]
            // Compares KRegCCD's model mass on novel-clade trees against the
            // empirical fraction of held-out trees that contain a novel clade.
            boolean twoFiles = args.length >= 3 && args[2].endsWith(".trees");
            int ai = twoFiles ? 3 : 2;
            int burnin = (args.length > ai) ? Integer.parseInt(args[ai]) : 10;
            double alpha = (args.length > ai + 1) ? Double.parseDouble(args[ai + 1]) : NNIRegCCD.DEFAULT_ALPHA;
            List<Tree> train, heldOut;
            if (twoFiles) {
                train = loadTrees(args[1], burnin);
                heldOut = loadTrees(args[2], burnin);
            } else {
                List<Tree> all = loadTrees(args[1], burnin);
                int half = all.size() / 2;
                train = new ArrayList<>(all.subList(0, half));
                heldOut = new ArrayList<>(all.subList(half, all.size()));
            }
            // empirical novel-clade fraction (depends only on the expanded clade set)
            KRegCCD ref = new KRegCCD(train, 0.0, 0.01, alpha, 2, KRegCCD.TailMode.NONE);
            int novel = 0;
            for (Tree t : heldOut) {
                if (ref.containsNovelClade(t)) {
                    novel++;
                }
            }
            double empirical = novel / (double) heldOut.size();
            System.out.printf("calibration: train = %d, held-out = %d, alpha = %.3f%n",
                    train.size(), heldOut.size(), alpha);
            System.out.printf("empirical novel-clade fraction = %.4f (%d/%d)%n",
                    empirical, novel, heldOut.size());
            System.out.println("  mu       model P(novel)  empirical  ratio   mean logP/tree");
            for (double mu : new double[]{0.002, 0.003, 0.0035, 0.004, 0.0045, 0.005, 0.006, 0.01}) {
                KRegCCD m = new KRegCCD(train, 0.0, mu, alpha, 2, KRegCCD.TailMode.NONE);
                double pNovel = m.getNovelCladeProbability();
                double sumLogP = 0.0;
                for (Tree t : heldOut) {
                    sumLogP += m.getLogProbabilityOfTree(t);
                }
                double meanLogP = sumLogP / heldOut.size();
                System.out.printf("  %-7.3f  %12.4f   %8.4f  %.2fx   %10.4f%n",
                        mu, pNovel, empirical, pNovel / empirical, meanLogP);
            }
            return;
        }
        if (args.length >= 1 && args[0].equals("cv")) {
            if (args.length < 3) {
                System.err.println("Usage: NNIHeldOutComparison cv <all.trees> <folds> [burnin%] [alpha] [strided|contiguous]");
                System.exit(1);
            }
            int folds = Integer.parseInt(args[2]);
            int burnin = (args.length >= 4) ? Integer.parseInt(args[3]) : 10;
            double alpha = (args.length >= 5) ? Double.parseDouble(args[4]) : NNIRegCCD.DEFAULT_ALPHA;
            FoldAssignment assignment = (args.length >= 6)
                    ? FoldAssignment.valueOf(args[5].toUpperCase()) : FoldAssignment.STRIDED;
            List<Tree> all = loadTrees(args[1], burnin);
            System.out.println(folds + "-fold cross-validation (" + assignment
                    + "), alpha = " + alpha + ", " + all.size() + " trees");
            for (Result result : crossValidate(all, folds, alpha, assignment)) {
                System.out.println(result);
            }
            return;
        }

        if (args.length < 1) {
            System.err.println("Usage: NNIHeldOutComparison <train.trees> <heldout.trees> [burnin%] [alpha]");
            System.err.println("   or: NNIHeldOutComparison <all.trees> [burnin%] [alpha]   (split in half)");
            System.err.println("   or: NNIHeldOutComparison cv <all.trees> <folds> [burnin%] [alpha]");
            System.exit(1);
        }

        List<Tree> training;
        List<Tree> heldOut;
        double alpha;
        int burnin;

        boolean twoFiles = args.length >= 2 && args[1].endsWith(".trees");
        if (twoFiles) {
            burnin = (args.length >= 3) ? Integer.parseInt(args[2]) : 10;
            alpha = (args.length >= 4) ? Double.parseDouble(args[3]) : NNIRegCCD.DEFAULT_ALPHA;
            training = loadTrees(args[0], burnin);
            heldOut = loadTrees(args[1], burnin);
        } else {
            burnin = (args.length >= 2) ? Integer.parseInt(args[1]) : 10;
            alpha = (args.length >= 3) ? Double.parseDouble(args[2]) : NNIRegCCD.DEFAULT_ALPHA;
            List<Tree> all = loadTrees(args[0], burnin);
            int half = all.size() / 2;
            training = new ArrayList<>(all.subList(0, half));
            heldOut = new ArrayList<>(all.subList(half, all.size()));
            System.out.println("Split one file: " + training.size() + " train / " + heldOut.size() + " held-out");
        }

        System.out.println("alpha = " + alpha + ", training = " + training.size()
                + " trees, held-out = " + heldOut.size() + " trees");
        for (Result result : compare(training, heldOut, alpha)) {
            System.out.println(result);
        }
    }
}
