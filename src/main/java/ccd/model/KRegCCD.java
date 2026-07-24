package ccd.model;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator;
import ccd.algorithms.regularisation.KRegCCDParameterOptimiser;
import ccd.algorithms.regularisation.NNIHeldOutComparison;
import ccd.model.bitsets.BitSet;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.ConcurrentHashMap;
import java.util.List;
import java.util.Map;

/**
 * Remco's "blue-region" CCD regularisation with <em>full support</em> (every tree
 * gets nonzero probability), parameterised by the reserve depth {@code k}.
 *
 * <p>
 * The split set is first widened by {@link RegCCD}'s CCD split expansion and
 * additive-{@code alpha} smoothing, so a split among observed clades is always
 * representable ("red", priced by its smoothed CCP). The only genuinely unseen
 * structures are then <em>novel clades</em>. To score a tree we colour each
 * internal node red (split in the expanded CCD) or blue (not — i.e. it involves a
 * novel clade), cluster blue nodes into maximal regions, and price:
 * <ul>
 *   <li>each red node at clade {@code C} by {@code (1 - mu) * p(split|C)};</li>
 *   <li>each blue region with top clade {@code C} and {@code j} novel clades by
 *       {@code eps(C)^j / pathcount} — one {@code log eps} per novel clade.
 *       {@code pathcount} is the number of all-novel resolutions of {@code C} into
 *       the region's boundary subclades (a cheap subset-DP over the boundary).</li>
 * </ul>
 * A region with a boundary of {@code m} observed subclades has {@code m - 2} novel
 * clades. Regions of <em>any</em> depth are scored, so any tree gets nonzero
 * probability.
 *
 * <p>
 * The single hyperparameter is {@code mu}, the per-clade escape probability.
 * {@code eps(C)} is the root of {@code R(C) = sum_{j>=1} N_j(C) eps^j = mu}, where
 * {@code N_j(C)} counts the region boundaries of {@code C} with {@code j} novel
 * clades. Computing the exact reserve is #P-hard (a set-partition count over a
 * non-laminar clade family), so we <em>decouple coverage from normalisation</em>:
 * {@code eps} is solved from only the cheap low orders {@code j = 1..k} (boundaries
 * {@code 3..k+2}; {@code N_1} is {@code O(m^2)}, {@code N_2} is {@code O(m^3)} in
 * the number {@code m} of observed subclades), while regions of any depth are still
 * scored with that {@code eps}. The omitted deep orders contribute {@code O(mu^2/m)}
 * to the reserve — a rounding error — so this barely affects probabilities while
 * keeping the model tractable and full-support. Reserve depth {@code k = 2} (the
 * default) is recommended; an op-budget (scaled to the largest clade's N_1, overridable via
 * {@code -Dkreg.enumOps})
 * bounds even the low-order enumeration on pathological clades, and a clade with
 * only deep regions climbs to its first nonzero order so {@code eps} stays positive.
 *
 * @author Claude (CCD-Sophie)
 */
public class KRegCCD extends RegCCD {

    /**
     * Per-clade enumeration-op budget: a clade's N_j enumeration is capped at this many steps so a
     * genuine combinatorial blow-up (deep orders of a large clade) is bounded, while the low orders
     * that set eps still complete.
     *
     * <p>The budget {@link #opsBudget} scales with the dataset: {@code OPS_BUDGET_FACTOR * m^2}, where
     * {@code m} is the largest clade's observed-subclade count (the root's) and {@code m^2} is the cost
     * of its order-1 (N_1) enumeration. So N_1 — which dominates eps — always completes regardless of
     * dataset size, while the expensive deep orders stay capped. A flat {@code 1e9} default was ~1000x
     * larger than needed on 129-taxon data (the reserve is unchanged from {@code ~m^2} upward). An
     * explicit {@code -Dkreg.enumOps} value (accepts {@code 1e7} or {@code 10000000}) overrides the scaling.
     */
    private static final long OPS_BUDGET_FACTOR = 8L;
    private static final long MIN_OPS_BUDGET = 2_000_000L;
    private final long opsBudget;

    /** Thrown internally when a single clade's N_j enumeration exceeds {@link #opsBudget}. */
    private static final class BudgetExceeded extends RuntimeException {
        BudgetExceeded() {
            super(null, null, false, false);
        }
    }

    private static final BudgetExceeded BUDGET_EXCEEDED = new BudgetExceeded();

    /** Per-thread enumeration-step counter for the op budget. Thread-local so {@link #precomputeReserves}
     *  can run clades in parallel: each clade's enumeration resets and accumulates its own thread's
     *  counter (every enumeration routine is reset-then-walk within a single thread). */
    private final ThreadLocal<long[]> enumOps = ThreadLocal.withInitial(() -> new long[1]);

    /**
     * Penalty exponent reduction: a boundary-{@code m} region (with {@code m-2}
     * novel clades) is priced {@code eps^(m - EXP_REDUCTION)} = {@code eps^(#novel
     * clades)}, i.e. one factor of {@code eps} per novel clade.
     */
    private static final int EXP_REDUCTION = 2;

    /** Per-clade escape probability (reserved mass for unseen resolutions). */
    private final double mu;

    /** Additive-smoothing pseudocount for the split-expanded backbone (the fitted alpha). */
    private final double alpha;

    private final double log1mMu;

    /**
     * Reserve depth: include novel-clade orders 1..k (boundaries 3..k+2) when
     * solving eps. Scoring is NOT truncated at this depth — every region is scored
     * (full support), so coverage is independent of k.
     */
    private final int reserveBoundary; // = k + 2

    /**
     * How to correct for the omitted reserve tail (orders &gt; reserve depth) when
     * discounting red splits by {@code 1 - mu - tail(C)}:
     * <ul>
     *   <li>{@link #NONE}: no correction ({@code tail = 0}); slightly
     *       super-normalised (inflates held-out scores by the tail).</li>
     *   <li>{@link #BOUND}: geometric upper bound on the tail; provably
     *       sub-normalised (never inflates).</li>
     *   <li>{@link #SAMPLED}: Knuth estimate of the actual tail; near-exactly
     *       normalised (up to Monte-Carlo noise).</li>
     * </ul>
     */
    public enum TailMode {NONE, BOUND, SAMPLED}

    private final TailMode tailMode;

    /**
     * How novel-clade mass is weighted across the distinct novel resolutions of a blue region.
     * <ul>
     *   <li>{@link #SHARED}: the original scheme. Each observed-subclade boundary reserves a total
     *       {@code eps^(m-2)}, shared evenly among its {@code pathcount} all-novel resolutions, so
     *       a tree is priced {@code eps^(m-2) / pathcount} and the reserve sums boundary counts
     *       {@code N_j}. This penalises novel clades clustered into one region (pathcount grows
     *       ~{@code (2k+1)!!}).</li>
     *   <li>{@link #FLAT}: every distinct novel resolution is priced {@code eps^(m-2)} directly
     *       (no {@code /pathcount}), and the reserve sums tree counts {@code M_j = sum of pathcount
     *       over boundaries}. This penalises each novel clade equally, independent of clustering.</li>
     * </ul>
     * Both are self-consistent, properly normalised priors (the reserve count and the per-tree
     * normaliser are a matched pair); they agree at the single-novel-clade order and diverge only
     * for clustered multi-clade regions.
     */
    public enum NovelMode {SHARED, FLAT}

    private final NovelMode novelMode;

    /** Default novel-clade weighting (clustering-neutral). */
    public static final NovelMode DEFAULT_NOVEL_MODE = NovelMode.FLAT;

    /**
     * Largest {@code mu} at which the depth-{@code k} reserve approximation stays reliable. Above
     * it the omitted reserve tail grows ~{@code mu^3} and the un-corrected (NONE) inflation becomes
     * material (&gt; ~0.01 nat/tree), so the model should not be operated — nor its parameters
     * searched (see {@link ccd.algorithms.regularisation.KRegCCDParameterOptimiser}) — beyond it.
     */
    public static final double MU_RELIABLE_MAX = 0.05;

    /** mu above which the un-corrected inflation becomes material (&gt; ~0.01 nat). */
    private static final double MU_WARN_THRESHOLD = MU_RELIABLE_MAX;

    /** Knuth samples per order when estimating the tail (SAMPLED mode). */
    private static final int TAIL_SAMPLES = 1000;

    /** Below this, the bounded tail is negligible, so SAMPLED skips the sampling. */
    private static final double TAIL_SAMPLE_THRESHOLD = 1e-5;

    private final java.util.Random rng = new java.util.Random(12345);

    /**
     * Fidelity of {@link #sampleTree()} relative to the full-support score
     * ({@link #getLogProbabilityOfTree}).
     * <ul>
     *   <li>{@link #SELF_CONSISTENT}: escape with probability {@code mu} and draw blue regions
     *       only up to the computed novelty orders (1..k, whose weights {@code N_j eps^j} sum to
     *       {@code mu} by the eps-solve). This is exactly {@code exp(score)} with the tail set to
     *       zero; regions deeper than {@code k} (mass ~{@code mu^3}) are never produced.</li>
     *   <li>{@link #FULL_SUPPORT}: escape with probability {@code mu + tail} and additionally draw
     *       the leading tail orders ({@code k+1, k+2}) via the Knuth estimator, so frequencies
     *       approach {@code exp(getLogProbabilityOfTree)} up to an {@code O(mu^4)} truncation of
     *       orders beyond {@code k+2}.</li>
     * </ul>
     * In both modes a region's boundary parts are resolved with an observed (red) split, which
     * makes the sampled blue regions <em>maximal</em> — matching the score's region decomposition
     * (and avoiding the over-counting a free recursion would cause). The sampling distribution is
     * exact whenever no boundary part is itself a reserving clade (always true for trees up to ~5
     * taxa); for larger trees it deviates from {@code exp(score)} by the boundary parts' own
     * {@code (1 - mu - tail)} reservation, an {@code O(mu)}-per-reserving-boundary-part effect.
     */
    public enum SamplingFidelity {SELF_CONSISTENT, FULL_SUPPORT}

    private SamplingFidelity samplingFidelity = SamplingFidelity.SELF_CONSISTENT;

    /** Reservoir state for {@link #sampleBoundary}: count of valid boundaries seen so far. */
    private int boundarySeen;
    private double boundaryWeightSeen; // weighted-reservoir accumulator for FLAT boundary sampling
    /** Reservoir state for {@link #sampleBoundary}: the boundary currently selected. */
    private BitSet[] boundaryPick;

    /**
     * Strictly-positive fallback increment for novel-node heights when the region top is not
     * strictly above a boundary part (e.g. common-ancestor-height non-monotonicity).
     */
    private static final double NOVEL_HEIGHT_EPS = 1e-8;

    /**
     * Per-clade reserve info: whether C reserves mu, log eps(C) (NEGATIVE_INFINITY
     * if nothing computable), and the tail correction to subtract from (1 - mu).
     *
     * <p>{@code nCounts} caches the per-boundary partition counts used to solve eps:
     * {@code nCounts[m]} is N_{m-2}(C) (the number of boundary-m partitions admitting an
     * all-novel resolution), for boundaries {@code 3..lastBoundary}. The sampler reuses these
     * to draw a region's novelty order j with probability proportional to N_j eps^j, exactly
     * matching the score (the same counts that pinned eps). {@code null}/{@code 0} when C is
     * not reservable.
     */
    private record CladeReg(boolean reservable, double logEps, double tail,
                            int[] nCounts, int lastBoundary) {
    }

    private final Map<Clade, CladeReg> regCache = new ConcurrentHashMap<>();

    /**
     * Default per-clade escape probability used when a tool builds a KRegCCD without
     * specifying one (e.g. via {@link CCDType#KRegCCD}); the RSV2 operating point.
     */
    public static final double DEFAULT_MU = 0.01;
    /**
     * Default additive-smoothing pseudocount for the split-expanded backbone (shared with
     * {@link RegCCD}'s default).
     */
    public static final double DEFAULT_ALPHA = 0.4;
    /** Default reserve depth k. */
    public static final int DEFAULT_RESERVE_DEPTH = 2;

    /**
     * Default variant: reserve depth k = 2 (eps from the one- and two-novel-clade
     * orders), full-support scoring.
     *
     * @param trees  training trees (post burn-in)
     * @param burnin fraction discarded as burn-in (0.0 if already removed)
     * @param mu     per-clade escape probability, in (0, 1)
     * @param alpha  additive-smoothing pseudocount for the split-expanded backbone
     */
    public KRegCCD(List<Tree> trees, double burnin, double mu, double alpha) {
        this(trees, burnin, mu, alpha, 2);
    }

    /**
     * @param k reserve depth (see below); tail bound correction on.
     */
    public KRegCCD(List<Tree> trees, double burnin, double mu, double alpha, int k) {
        this(trees, burnin, mu, alpha, k, TailMode.BOUND);
    }

    /**
     * @param trees    training trees (post burn-in)
     * @param burnin   fraction discarded as burn-in (0.0 if already removed)
     * @param mu       per-clade escape probability, in (0, 1)
     * @param alpha    additive-smoothing pseudocount for the split-expanded backbone
     * @param k        reserve depth: include novel-clade orders {@code 1..k} when
     *                 solving {@code eps} (>= 1). Scoring is NOT truncated at
     *                 {@code k}, so coverage is full regardless; {@code k} only
     *                 affects accuracy.
     * @param tailMode how to correct for the omitted reserve tail (see {@link TailMode})
     */
    public KRegCCD(List<Tree> trees, double burnin, double mu, double alpha, int k,
                   TailMode tailMode) {
        this(trees, burnin, mu, alpha, k, tailMode, DEFAULT_NOVEL_MODE);
    }

    /**
     * @param novelMode how novel-clade mass is weighted across resolutions (see {@link NovelMode})
     */
    public KRegCCD(List<Tree> trees, double burnin, double mu, double alpha, int k,
                   TailMode tailMode, NovelMode novelMode) {
        super(trees, burnin, false); // plain CCD1; split expansion + smoothing applied below
        validateParams(mu, k, alpha);
        this.mu = mu;
        this.alpha = alpha;
        this.log1mMu = Math.log(1.0 - mu);
        this.reserveBoundary = k + 2; // boundary size = novel clades + 2
        this.tailMode = tailMode;
        this.novelMode = novelMode;
        expandRegCCD();
        regulariseRegCCD(alpha);
        this.opsBudget = computeOpsBudget();
    }

    /**
     * Constructor from a burn-in-free tree set (used by the CCD tools, e.g. {@link CCDType#KRegCCD}
     * via {@code CCDToolUtil}). Mirrors the {@link List}-based constructor: builds the plain CCD1
     * backbone, then split-expands and additively smooths it, leaving eps to be solved per clade
     * on demand.
     *
     * @param treeSet  trees without burn-in whose distribution is approximated
     * @param mu       per-clade escape probability, in (0, 1)
     * @param alpha    additive-smoothing pseudocount for the split-expanded backbone
     * @param k        reserve depth (see the {@link List}-based constructor)
     * @param tailMode how to correct for the omitted reserve tail (see {@link TailMode})
     */
    public KRegCCD(TreeAnnotator.TreeSet treeSet, double mu, double alpha, int k, TailMode tailMode) {
        this(treeSet, mu, alpha, k, tailMode, DEFAULT_NOVEL_MODE);
    }

    /**
     * @param novelMode how novel-clade mass is weighted across resolutions (see {@link NovelMode})
     */
    public KRegCCD(TreeAnnotator.TreeSet treeSet, double mu, double alpha, int k, TailMode tailMode,
                   NovelMode novelMode) {
        super(treeSet, false); // plain CCD1 from the tree set; expansion + smoothing applied below
        validateParams(mu, k, alpha);
        this.mu = mu;
        this.alpha = alpha;
        this.log1mMu = Math.log(1.0 - mu);
        this.reserveBoundary = k + 2;
        this.tailMode = tailMode;
        this.novelMode = novelMode;
        expandRegCCD();
        regulariseRegCCD(alpha);
        this.opsBudget = computeOpsBudget();
    }

    /**
     * Builds a KRegCCD on {@code trees} with {@code (alpha, mu)} selected by maximising
     * cross-validated held-out tree log-probability (see {@link KRegCCDParameterOptimiser}),
     * rather than using the fixed {@link #DEFAULT_MU}/{@link #DEFAULT_ALPHA}. Reserve depth is
     * {@link #DEFAULT_RESERVE_DEPTH} and the tail correction is {@link TailMode#BOUND}.
     *
     * @param trees the trees the model is built on and cross-validated over
     * @param folds number of cross-validation folds
     */
    public static KRegCCD withOptimisedParameters(List<Tree> trees, int folds) {
        KRegCCDParameterOptimiser.Params p = KRegCCDParameterOptimiser.optimise(trees, folds, NNIHeldOutComparison.FoldAssignment.CONTIGUOUS);
        return new KRegCCD(trees, 0.0, p.mu(), p.alpha(), DEFAULT_RESERVE_DEPTH, TailMode.BOUND);
    }

    /** As {@link #withOptimisedParameters(List, int)} with {@link KRegCCDParameterOptimiser#DEFAULT_FOLDS} folds. */
    public static KRegCCD withOptimisedParameters(List<Tree> trees) {
        return withOptimisedParameters(trees, KRegCCDParameterOptimiser.DEFAULT_FOLDS);
    }

    public static KRegCCD withOptimisedMu(List<Tree> trees) {
        KRegCCDParameterOptimiser.MuResult r =
                KRegCCDParameterOptimiser.optimiseMu(trees, KRegCCDParameterOptimiser.DEFAULT_FOLDS,
                        NNIHeldOutComparison.FoldAssignment.CONTIGUOUS, DEFAULT_ALPHA);

        System.out.println("Optimised mu = " + r.mu() + ", held-out logP = " + r.heldOutLogProb());

        return new KRegCCD(trees, 0.0, r.mu(), DEFAULT_ALPHA, DEFAULT_RESERVE_DEPTH, TailMode.BOUND);
    }

    /** The additive-smoothing pseudocount alpha this model was built/fitted with. */
    public double getAlpha() {
        return alpha;
    }

    /** The per-clade escape probability mu this model was built/fitted with. */
    public double getMu() {
        return mu;
    }

    private static void validateParams(double mu, int k, double alpha) {
        if (mu <= 0 || mu >= 1) {
            throw new IllegalArgumentException("mu must be in (0, 1), got " + mu);
        }
        if (k < 1) {
            throw new IllegalArgumentException("reserve depth k must be >= 1, got " + k);
        }
        if (Double.isNaN(alpha)) {
            throw new IllegalArgumentException("KRegCCD requires split expansion (finite alpha)");
        }
        if (mu > MU_WARN_THRESHOLD) {
            System.err.println("WARNING: KRegCCD mu = " + mu + " > " + MU_WARN_THRESHOLD
                    + "; the reserve-tail inflation grows ~mu^3 and becomes material here "
                    + "(~" + String.format("%.2g", Math.pow(mu / 0.05, 3) * 0.01)
                    + " nat/tree at this mu). Use a smaller mu, or TailMode.BOUND/SAMPLED.");
        }
    }

    @Override
    public String toString() {
        return "KRegCCD [mu = " + mu + ", reserve depth k = " + (reserveBoundary - 2)
                + ", tail=" + tailMode + ", novel=" + novelMode + ", full support, split-expanded]";
    }

    /* ----------------------------------------------------------------------
     * Tree scoring
     * ------------------------------------------------------------------- */

    @Override
    public double getLogProbabilityOfTree(Tree tree) {
        Map<Node, BitSet> bits = new HashMap<>();
        computeBits(tree.getRoot(), bits);
        double[] logp = new double[]{0.0};
        scoreFresh(tree.getRoot(), bits, logp, mu, true);
        return logp[0];
    }

    /**
     * Full-support log-probability of a tree at an arbitrary escape probability {@code scoreMu},
     * reusing this model's ({@code mu}-independent) clade/split backbone and cached per-clade
     * {@code N_j} counts. Scored with no reserve-tail correction (tail = 0, i.e.
     * {@link TailMode#NONE} semantics), which is the objective used to select {@code mu}. This
     * lets a parameter optimiser evaluate many {@code mu} values on one trained backbone without
     * rebuilding — {@code eps(C)} is just re-solved from the cached counts. For
     * {@code scoreMu == mu} and a {@code NONE}-built model it equals {@link #getLogProbabilityOfTree(Tree)}.
     */
    public double getLogProbabilityOfTree(Tree tree, double scoreMu) {
        if (scoreMu <= 0 || scoreMu >= 1) {
            throw new IllegalArgumentException("scoreMu must be in (0, 1), got " + scoreMu);
        }
        Map<Node, BitSet> bits = new HashMap<>();
        computeBits(tree.getRoot(), bits);
        double[] logp = new double[]{0.0};
        scoreFresh(tree.getRoot(), bits, logp, scoreMu, false);
        return logp[0];
    }

    /**
     * Buffer-reusing variant of {@link #getLogProbabilityOfTree(Tree)} for repeated scoring (e.g. an SPR
     * operator scoring every re-attachment of a pruned subtree): the caller supplies a scratch map whose
     * per-node {@link BitSet}s are CLEARED and reused rather than reallocated, removing the per-call
     * {@code HashMap} and per-node {@code BitSet} allocation that dominates such loops. Entries for nodes not in
     * {@code tree} are harmless (only nodes reachable from the root are written and read). Returns exactly the
     * same value as {@link #getLogProbabilityOfTree(Tree)}.
     */
    public double getLogProbabilityOfTree(Tree tree, Map<Node, BitSet> scratch) {
        computeBitsReusing(tree.getRoot(), scratch);
        double[] logp = new double[]{0.0};
        scoreFresh(tree.getRoot(), scratch, logp, mu, true);
        return logp[0];
    }

    /** As {@link #computeBits} but reuses each node's {@link BitSet} from {@code bits} (clear + recompute)
     *  instead of allocating a fresh one, so repeated calls on stable node objects do not allocate. */
    private BitSet computeBitsReusing(Node v, Map<Node, BitSet> bits) {
        BitSet b = bits.get(v);
        if (b == null) {
            b = BitSet.newBitSet(leafArraySize);
            bits.put(v, b);
        } else {
            b.clear();
        }
        if (v.isLeaf()) {
            b.set(v.getNr());
        } else {
            b.or(computeBitsReusing(v.getChildren().get(0), bits));
            b.or(computeBitsReusing(v.getChildren().get(1), bits));
        }
        return b;
    }

    /**
     * Incrementally score EVERY single-leaf re-graft of leaf {@code leafIndex} onto {@code reduced} (a tree that
     * does NOT contain that leaf): for each target node {@code x}, the full-support log-probability of the tree
     * obtained by attaching the leaf on the branch above {@code x}. This is the engine for a CCD-guided SPR/Gibbs
     * operator that must score all {@code ~2n} re-attachments per move; scoring each independently with
     * {@link #getLogProbabilityOfTree(Tree)} is {@code O(n)} each ({@code O(n^2)} total), dominated by re-walking
     * the unchanged backbone. Here the backbone is priced ONCE: {@code subScore[v]} = the self-contained
     * full-support score of {@code v}'s subtree, then each re-graft is a MEMOISED {@code scoreFresh} that returns
     * the cached {@code subScore} for any subtree not containing the leaf, so the recursion only descends the
     * leaf's insertion path. Result is identical to {@link #getLogProbabilityOfTree(Tree)} of the spliced tree.
     *
     * @param reduced   a tree on the taxon set MINUS the leaf (leaves carry this CCD's canonical numbering)
     * @param leafIndex the canonical index of the leaf being re-grafted
     * @param targets   the nodes above which to score an attachment (must exclude the root); {@code null} = all
     *                  non-root nodes
     * @return map from each target node to the log-probability of attaching the leaf above it
     */
    public Map<Node, Double> logProbabilityOfAllRegrafts(Tree reduced, int leafIndex, List<Node> targets) {
        Map<Node, BitSet> bits = new HashMap<>();
        computeBits(reduced.getRoot(), bits);

        // price the backbone once: self-contained subtree score for every node (only OBSERVED-clade nodes are ever
        // read back, but the post-order recurses through novel-clade interiors to reach their observed boundaries)
        Map<Node, Double> subScore = new HashMap<>();
        computeSubtreeScores(reduced.getRoot(), bits, subScore);

        // spare nodes spliced in to form each candidate tree (reused across candidates)
        BitSet iBits = BitSet.newBitSet(leafArraySize);
        iBits.set(leafIndex);
        Node iLeaf = new Node();
        iLeaf.setNr(leafIndex);
        bits.put(iLeaf, iBits);
        Node p = new Node();

        List<Node> tgts = targets;
        if (tgts == null) {
            tgts = new ArrayList<>();
            collectNonRoot(reduced.getRoot(), reduced.getRoot(), tgts);   // getNodesAsArray() may be uninitialised
        }
        Map<Node, Double> out = new HashMap<>(tgts.size() * 2);
        for (Node x : tgts) {
            out.put(x, scoreOneRegraft(reduced.getRoot(), x, leafIndex, iBits, iLeaf, p, bits, subScore));
        }
        return out;
    }

    /** Score the single tree formed by attaching the leaf above {@code x}: splice {@code p=(leaf, x)} into x's
     *  edge, override the bits of {@code p} and the path {@code x.parent..root} to carry the leaf, run the
     *  memoised scoreFresh, then restore everything (so {@code bits}/the tree are unchanged for the next target). */
    private double scoreOneRegraft(Node root, Node x, int leafIndex, BitSet iBits, Node iLeaf, Node p,
                                   Map<Node, BitSet> bits, Map<Node, Double> subScore) {
        Node xp = x.getParent();
        xp.removeChild(x);
        xp.addChild(p);
        p.addChild(iLeaf);
        p.addChild(x);
        BitSet pBits = (BitSet) bits.get(x).clone();
        pBits.or(iBits);
        bits.put(p, pBits);
        List<Node> pathNodes = new ArrayList<>();
        List<BitSet> pathSaved = new ArrayList<>();
        for (Node a = xp; a != null; a = a.getParent()) {
            pathNodes.add(a);
            pathSaved.add(bits.get(a));
            BitSet ab = (BitSet) bits.get(a).clone();
            ab.or(iBits);
            bits.put(a, ab);
        }
        double[] logp = {0.0};
        scoreFreshMemo(root, bits, logp, leafIndex, subScore);
        for (int q = 0; q < pathNodes.size(); q++) bits.put(pathNodes.get(q), pathSaved.get(q));
        bits.remove(p);
        xp.removeChild(p);
        p.removeChild(iLeaf);
        p.removeChild(x);
        xp.addChild(x);
        return logp[0];
    }

    /** {@code subScore[v]} = self-contained full-support score of v's subtree, computed bottom-up (mirrors a
     *  standalone {@link #scoreFresh} of v). Only OBSERVED-clade nodes' values are ever read; novel-clade nodes get
     *  0 (they only ever appear as region interior, never as a memoised recursion root). */
    private double computeSubtreeScores(Node v, Map<Node, BitSet> bits, Map<Node, Double> subScore) {
        if (v.isLeaf()) {
            subScore.put(v, 0.0);
            return 0.0;
        }
        Node c1 = v.getChildren().get(0), c2 = v.getChildren().get(1);
        computeSubtreeScores(c1, bits, subScore);
        computeSubtreeScores(c2, bits, subScore);
        BitSet vb = bits.get(v), c1b = bits.get(c1), c2b = bits.get(c2);
        double s;
        if (isRed(vb, c1b, c2b)) {
            Clade c = getClade(vb);
            CladeReg reg = computeReg(c);
            s = (reg.reservable() ? Math.log(1.0 - mu - reg.tail()) : 0.0)
                    + rawLogCCP(c, c1b, c2b) + subScore.get(c1) + subScore.get(c2);
        } else if (getClade(vb) != null) {
            // observed clade, novel split: standalone region top
            List<Node> boundary = new ArrayList<>();
            collectRegion(v, bits, boundary);
            s = regionScore(vb, getClade(vb), boundary.size(), boundary, bits);
            for (Node b : boundary) s += subScore.get(b);
        } else {
            s = 0.0;   // novel clade: never read back
        }
        subScore.put(v, s);
        return s;
    }

    /** Memoised {@link #scoreFresh}: identical pricing, but recursion into a subtree NOT containing the re-grafted
     *  leaf returns the cached {@code subScore} instead of descending. Only ever entered on observed-clade nodes. */
    private void scoreFreshMemo(Node v, Map<Node, BitSet> bits, double[] logp, int leafIndex,
                                Map<Node, Double> subScore) {
        if (v.isLeaf()) {
            return;
        }
        Node c1 = v.getChildren().get(0), c2 = v.getChildren().get(1);
        BitSet vb = bits.get(v), c1b = bits.get(c1), c2b = bits.get(c2);
        if (isRed(vb, c1b, c2b)) {
            Clade c = getClade(vb);
            CladeReg reg = computeReg(c);
            if (reg.reservable()) {
                logp[0] += Math.log(1.0 - mu - reg.tail());
            }
            logp[0] += rawLogCCP(c, c1b, c2b);
            memoRecurse(c1, bits, logp, leafIndex, subScore);
            memoRecurse(c2, bits, logp, leafIndex, subScore);
        } else {
            List<Node> boundary = new ArrayList<>();
            collectRegion(v, bits, boundary);
            logp[0] += regionScore(vb, getClade(vb), boundary.size(), boundary, bits);
            for (Node b : boundary) {
                memoRecurse(b, bits, logp, leafIndex, subScore);
            }
        }
    }

    private void collectNonRoot(Node v, Node root, List<Node> out) {
        if (v != root) out.add(v);
        for (Node c : v.getChildren()) collectNonRoot(c, root, out);
    }

    private void memoRecurse(Node child, Map<Node, BitSet> bits, double[] logp, int leafIndex,
                             Map<Node, Double> subScore) {
        if (bits.get(child).get(leafIndex)) {
            scoreFreshMemo(child, bits, logp, leafIndex, subScore);
        } else {
            logp[0] += subScore.get(child);
        }
    }

    /** The blue-region contribution {@code (m - 2) * logEps(top) - logPath} (fallback {@code log(mu) - logPath}),
     *  exactly as {@link #scoreFresh}'s else-branch (useStoredReg path). */
    private double regionScore(BitSet vb, Clade top, int m, List<Node> boundary, Map<Node, BitSet> bits) {
        double logPath = 0.0;
        if (novelMode != NovelMode.FLAT) {
            BitSet[] parts = new BitSet[m];
            for (int i = 0; i < m; i++) parts[i] = bits.get(boundary.get(i));
            logPath = Math.log(countAllNovelResolutions(vb, parts));
        }
        double logEps = computeReg(top).logEps();
        return (logEps == Double.NEGATIVE_INFINITY) ? Math.log(mu) - logPath : (m - EXP_REDUCTION) * logEps - logPath;
    }

    /**
     * Probability of the given tree under the full-support model, {@code exp(getLogProbabilityOfTree)}.
     *
     * <p>Overrides {@link AbstractCCD#getProbabilityOfTree}, whose red-only CCP product is wrong here:
     * it returns 0 for any tree containing a novel clade (this model gives every tree nonzero mass) and
     * omits the per-clade {@code (1 - mu)} escape discount even for fully observed trees. Note this can
     * underflow to 0 for large trees; use {@link #getLogProbabilityOfTree(Tree)} when that matters.
     */
    @Override
    public double getProbabilityOfTree(Tree tree) {
        return Math.exp(getLogProbabilityOfTree(tree));
    }

    /**
     * Always {@code true}: KRegCCD is full support, so every tree on this taxon set has nonzero
     * probability. Overrides the inherited {@code getProbabilityOfTree(tree) > 0} test, which would
     * both inherit the red-only bug above and report a false negative when a large tree's probability
     * underflows to 0.
     */
    @Override
    public boolean containsTree(Tree tree) {
        return true;
    }

    /**
     * Number of internal clades of {@code tree} that are not present in this CCD's observed clade
     * set (the novel-clade count = the total {@code m-2} over the tree's blue regions). This is the
     * exponent that the per-region {@code eps} factors are raised to, and is mode-independent (it
     * depends only on the observed backbone, not on SHARED/FLAT weighting). Useful for bucketing
     * held-out trees by how much novel (clustered) structure they carry.
     */
    public int novelCladeCount(Tree tree) {
        Map<Node, BitSet> bits = new HashMap<>();
        computeBits(tree.getRoot(), bits);
        int count = 0;
        for (Node v : tree.getNodesAsArray()) {
            if (v.isLeaf()) {
                continue;
            }
            if (getClade(bits.get(v)) == null) {
                count++;
            }
        }
        return count;
    }

    /* ----------------------------------------------------------------------
     * MAP tree (discounted DP)
     *
     * The inherited max-probability DP (AbstractCCD/Clade.getMaxSubtreeLogCCP) maximises only the
     * sum of observed log-CCPs -- the red-only backbone. It is wrong for a full-support model in two
     * ways: (1) it OMITS the per-clade red discount log(1 - mu - tail(C)) that scoreFresh charges at
     * every reservable clade, so its returned max probability is inflated; and (2) because that
     * discount is incurred once per reservable internal clade and reservability differs between
     * topologies, the discount can change which resolution is optimal -- so the max-CCP topology need
     * not be the full-support MAP. We replace it with a discounted DP over the same observed-clade DAG:
     *
     *   D(c)  = disc(c) + max_p [ logCCP(c -> L,R) + D(L) + D(R) ],   D(leaf) = 0
     *   disc(c) = reservable(c) ? log(1 - mu - tail(c)) : 0
     *
     * which is exactly scoreFresh's score of an all-red tree, maximised. Escapes are never taken: a
     * blue region contributes eps^(m-2)/pathcount (with eps << mu), strictly less than the red
     * alternative's (1 - mu) * ccp at the operating mu (<= MU_RELIABLE_MAX), so the global MAP is the
     * best all-red tree this DP finds. getMAPTree follows the same argmax partitions (via the MAP
     * branch of getPartitionBasedOnStrategy), so getMaxLogTreeProbability() ==
     * getLogProbabilityOfTree(getMAPTree()) by construction.
     * ------------------------------------------------------------------- */

    private Map<Clade, Double> mapLogProbMemo;
    private Map<Clade, CladePartition> mapPartitionMemo;

    /**
     * Computes the discounted MAP DP over all clades (children before parents) and caches, per
     * clade, the best all-red subtree log-probability and the argmax partition. Idempotent.
     */
    private void computeDiscountedMap() {
        if (mapLogProbMemo != null) {
            return;
        }
        Map<Clade, Double> logProb = new HashMap<>();
        Map<Clade, CladePartition> argmax = new HashMap<>();
        List<Clade> clades = new ArrayList<>(getClades());
        clades.sort((a, b) -> Integer.compare(a.size(), b.size()));
        for (Clade c : clades) {
            if (c.isLeaf()) {
                logProb.put(c, 0.0);
                continue;
            }
            double best = Double.NEGATIVE_INFINITY;
            CladePartition bestPartition = null;
            for (CladePartition p : c.getPartitions()) {
                Clade[] ch = p.getChildClades();
                double v = p.getLogCCP() + logProb.get(ch[0]) + logProb.get(ch[1]);
                if (v > best) { // first-wins on ties: getPartitions() order is deterministic
                    best = v;
                    bestPartition = p;
                }
            }
            CladeReg reg = computeReg(c);
            double disc = reg.reservable() ? Math.log(1.0 - mu - reg.tail()) : 0.0;
            logProb.put(c, disc + best);
            argmax.put(c, bestPartition);
        }
        mapLogProbMemo = logProb;
        mapPartitionMemo = argmax;
    }

    /**
     * Log-probability of the most likely tree under the full-support model: the discounted DP value
     * at the root (see the section comment above). Overrides the inherited red-only maximum, which
     * omits the per-clade {@code (1 - mu - tail)} escape discount and so overstates the maximum.
     */
    @Override
    public double getMaxLogTreeProbability() {
        computeDiscountedMap();
        return mapLogProbMemo.get(getRootClade());
    }

    /**
     * Routes the MAP partition selection through the discounted DP so {@link #getMAPTree()} returns
     * the full-support MAP tree (whose log-probability is {@link #getMaxLogTreeProbability()}), not
     * the red-only max-CCP tree. Other strategies keep the inherited behaviour.
     */
    @Override
    protected CladePartition getPartitionBasedOnStrategy(Clade clade, SamplingStrategy samplingStrategy) {
        if (samplingStrategy == SamplingStrategy.MAP) {
            computeDiscountedMap();
            return mapPartitionMemo.get(clade);
        }
        return super.getPartitionBasedOnStrategy(clade, samplingStrategy);
    }

    /* ----------------------------------------------------------------------
     * Entropy
     *
     * The inherited CCD entropy (AbstractCCD.getEntropy / getEntropyLewis, the
     * Lewis et al. 2016 recursion) sums only over observed clade partitions and
     * treats their CCPs as a complete distribution. That is the entropy of the
     * renormalised red-only backbone (the smoothed, split-expanded CCD1) -- it
     * silently drops the escape mass mu, the blue regions, and the whole
     * novel-clade support, so it is the WRONG quantity for the full-support
     * KRegCCD. We provide two estimators of the real thing instead:
     *
     *  - getEntropyMonteCarlo: a plug-in average -mean(log P) over trees drawn
     *    from the model's own sampler. Simple, public-API only, with a reported
     *    standard error; bias is the O(mu)-per-reserving-boundary-part gap the
     *    sampler already carries (negligible at the operating mu <= 0.05).
     *
     *  - getEntropyRecursive: a generalisation of the Lewis recursion to the
     *    escape mixture. Exact for the self-consistent escape distribution
     *    (tail = 0, escape mass exactly mu, region orders 1..k), modulo the same
     *    op-budget truncation the rest of the model uses; far lower variance.
     *
     * {@link #getEntropy()} mirrors {@link #sampleTree()}: it reports the entropy of whichever
     * distribution the sampler would draw at the current {@link SamplingFidelity} -- the exact
     * recursion for SELF_CONSISTENT, the Monte-Carlo plug-in for FULL_SUPPORT (whose reserve tail
     * has no cheap closed form). Under SELF_CONSISTENT both estimators target the same normalised,
     * tail-excluded distribution, so they are directly comparable; on small trees, where the model
     * normalises exactly, they agree with the brute-force entropy.
     * ------------------------------------------------------------------- */

    /**
     * Number of sampled trees used by {@link #getEntropy()} (the default Monte-Carlo path is
     * not used; the recursion is). Kept for the convenience overload.
     */
    public static final int DEFAULT_ENTROPY_SAMPLES = 100_000;

    /**
     * The inherited Lewis-recursion entropy is invalid for a full-support model (it ignores escape
     * mass and novel clades). Use {@link #getEntropy()} (the generalised recursion) or
     * {@link #getEntropyMonteCarlo(int)} instead.
     */
    @Override
    public double getEntropyLewis() {
        throw new UnsupportedOperationException(
                "The Lewis-recursion entropy only covers observed clade partitions and so omits "
                        + "KRegCCD's escape mass and novel-clade support. Use getEntropy() (the "
                        + "generalised recursion) or getEntropyMonteCarlo(samples).");
    }

    /**
     * Entropy (in nats) of the KRegCCD topology distribution that {@link #sampleTree()} would draw
     * from at the current {@link SamplingFidelity} -- so the reported entropy always matches the
     * distribution the sampler produces. Like the sampler, the mode is read from
     * {@link #setSamplingFidelity}:
     * <ul>
     *   <li>{@link SamplingFidelity#SELF_CONSISTENT}: the normalised distribution with escape mass
     *       exactly {@code mu} and blue regions only up to reserve depth {@code k}. Computed by the
     *       exact, deterministic recursion {@link #getEntropyRecursive()}.</li>
     *   <li>{@link SamplingFidelity#FULL_SUPPORT}: the full model, including the reserve tail
     *       (regions deeper than {@code k}). The tail is the same #P-hard quantity the sampler
     *       Knuth-estimates, so there is no cheap exact recursion for it; this falls back to the
     *       Monte-Carlo plug-in {@link #getEntropyMonteCarlo(int)} ({@link #DEFAULT_ENTROPY_SAMPLES}
     *       samples), exactly as full-support <em>sampling</em> relies on the Knuth tail.</li>
     * </ul>
     * For the two estimators to be mutually consistent (and the Monte-Carlo path unbiased), pair the
     * fidelity with the matching {@link TailMode} the sampler expects -- {@code SELF_CONSISTENT} with
     * {@link TailMode#NONE}, {@code FULL_SUPPORT} with {@link TailMode#SAMPLED} -- so the sampler's
     * red discount {@code (1 - mu - tail)} matches its escape mass. Overrides the inherited
     * observed-only formula, which is wrong for a full-support model.
     */
    @Override
    public double getEntropy() {
        return switch (samplingFidelity) {
            case SELF_CONSISTENT -> getEntropyRecursive();
            case FULL_SUPPORT -> {
                if (tailMode == TailMode.NONE) {
                    System.err.println("WARNING: getEntropy() in FULL_SUPPORT fidelity on a "
                            + "TailMode.NONE model is biased: the sampler escapes with mu + tail but "
                            + "the red discount omits the tail (super-normalised). Build with "
                            + "TailMode.SAMPLED for an unbiased full-support entropy.");
                }
                yield getEntropyMonteCarlo(DEFAULT_ENTROPY_SAMPLES)[0];
            }
        };
    }

    /**
     * Monte-Carlo plug-in estimate of the entropy: draw {@code samples} trees from the model's
     * sampler and average their negative log-probability,
     * {@code H = -mean(log P(T_s))}. Returns {@code {entropy, standardError}} (both in nats), the
     * standard error being {@code sd(log P) / sqrt(samples)}.
     *
     * <p>This estimates {@code -E_Q[log P]} where {@code Q} is the sampler's distribution and
     * {@code P} the scored distribution; the two coincide except for the sampler's
     * {@code O(mu)}-per-reserving-boundary-part reservation (the same gap documented on
     * {@link SamplingFidelity}), so on trees small enough that no boundary part reserves (P == Q,
     * exactly normalised) this is an unbiased estimate of the true entropy. Set the RNG via
     * {@link #setRandom} for reproducibility; uses the current {@link SamplingFidelity}.
     */
    public double[] getEntropyMonteCarlo(int samples) {
        if (samples < 2) {
            throw new IllegalArgumentException("need at least 2 samples for a standard error");
        }
        double sumLogP = 0.0;
        double sumLogP2 = 0.0;
        for (int s = 0; s < samples; s++) {
            Tree t = sampleTree(HeightSettingStrategy.None);
            double logp = (Double) t.getRoot().getMetaData(LOG_PROB_SUBTREE_KEY);
            sumLogP += logp;
            sumLogP2 += logp * logp;
        }
        double meanLogP = sumLogP / samples;
        double varLogP = Math.max(0.0, sumLogP2 / samples - meanLogP * meanLogP);
        double stdErr = Math.sqrt(varLogP / samples);
        return new double[]{-meanLogP, stdErr};
    }

    /**
     * Generalised Lewis recursion for the full-support entropy. The KRegCCD distribution is
     * generated by the same recursive branching process the sampler walks, so entropy decomposes by
     * the chain rule over the observed-clade DAG. At a reservable clade {@code C} the immediate
     * decision is a mixture: with probability {@code 1 - mu} an observed (red) split (weighted by
     * its CCP), and with probability {@code mu} an escape into a blue region. Writing
     * {@code H_free(C)} for the entropy of the subtree generated when {@code C} may escape and
     * {@code H_red(C)} for when {@code C} is forced to take an observed split (a blue region's
     * boundary part), with escape mass {@code mu} and {@code eps = eps(C)}:
     * <pre>
     *   H_red(C)  = sum_(L,R) ccp(L,R) * ( -log ccp(L,R) + H_free(L) + H_free(R) )
     *   H_free(C) = (1 - mu) * ( -log(1 - mu) + H_red(C) )                       // red branch
     *             + (-log eps) * sum_m (m-2) * N_m * eps^(m-2)                   // blue self-information
     *             + sum_(regions) P(region) * sum_i H_red(B_i)                   // blue children (forced red)
     * </pre>
     * The red part and the blue self-information are exact (finite sums over the computed orders).
     * The blue-children term is computed exactly by enumerating the same canonical boundaries the
     * reserve uses (FLAT weights each boundary by its pathcount, SHARED by one), bounded by the same
     * op-budget; for SHARED the per-boundary {@code log pathcount} self-information is accumulated in
     * the same pass. Clades are processed in increasing size so every child and boundary part is
     * already memoised when a clade is reached (the enumeration then does no nested re-entry).
     *
     * @return the entropy (in nats) of the self-consistent (tail = 0, escape mass {@code mu})
     * KRegCCD distribution
     */
    public double getEntropyRecursive() {
        Map<Clade, Double> freeMemo = new ConcurrentHashMap<>();
        Map<Clade, Double> redMemo = new ConcurrentHashMap<>();
        // Both entropyRedForced(c) and entropyFree(c) read only strictly-smaller clades (a clade's
        // partition children and its blue-region boundary parts are all smaller than it), so once
        // every smaller size level is memoised, all clades of a given size are independent and can
        // be computed in parallel -- exactly as precomputeReserves parallelises computeReg. Each
        // task writes only its own clade's freeMemo/redMemo entries (a concurrent map), and the
        // per-clade arithmetic is unchanged, so the entropy is identical to the serial computation.
        List<List<Clade>> bySize = new ArrayList<>(leafArraySize + 1);
        for (int s = 0; s <= leafArraySize; s++) {
            bySize.add(new ArrayList<>());
        }
        for (Clade c : getClades()) {
            bySize.get(c.size()).add(c);
        }
        for (int s = 1; s <= leafArraySize; s++) {
            List<Clade> level = bySize.get(s);
            if (!level.isEmpty()) {
                level.parallelStream().forEach(c -> {
                    entropyRedForced(c, freeMemo, redMemo);
                    entropyFree(c, freeMemo, redMemo);
                });
            }
        }
        return freeMemo.getOrDefault(getRootClade(), 0.0);
    }

    /**
     * Entropy of the subtree generated at {@code c} when {@code c} may escape (the root, or a red
     * child of another node). Assumes all strictly-smaller clades are already memoised.
     */
    private double entropyFree(Clade c, Map<Clade, Double> freeMemo, Map<Clade, Double> redMemo) {
        if (c.isLeaf()) {
            return 0.0;
        }
        Double cached = freeMemo.get(c);
        if (cached != null) {
            return cached;
        }
        CladeReg reg = computeReg(c);
        double h;
        if (!reg.reservable()) {
            h = entropyRedForced(c, freeMemo, redMemo); // no escape mass: pure observed-split node
        } else {
            double eps = Math.exp(reg.logEps());
            int[] n = reg.nCounts();
            int last = reg.lastBoundary();
            // self-consistent order weights w[m] = N_m * eps^(m-2), m = 3..last; sum == mu by the
            // eps-solve. "blueSelfInfoSum" accumulates sum_m (m-2) w[m] for the blue self-information.
            double escapeMass = 0.0;
            double blueSelfInfoSum = 0.0;
            for (int m = 3; m <= last; m++) {
                if (n[m] > 0) {
                    double w = n[m] * Math.pow(eps, m - EXP_REDUCTION);
                    escapeMass += w;
                    blueSelfInfoSum += (m - EXP_REDUCTION) * w;
                }
            }
            double redMass = 1.0 - escapeMass;
            double redContribution = redMass * (-Math.log(redMass) + entropyRedForced(c, freeMemo, redMemo));
            BlueTerms blue = blueContribution(c, reg, eps, freeMemo, redMemo);
            // FLAT: blue self-information is exactly -log(eps) * sum_m (m-2) w[m]. SHARED adds the
            // per-region -log(1/pathcount) = +log(pathcount) self-information collected by the enumeration.
            double blueSelfInfo = -reg.logEps() * blueSelfInfoSum + blue.extraLogPc();
            h = redContribution + blueSelfInfo + blue.childEntropy();
        }
        freeMemo.put(c, h);
        return h;
    }

    /**
     * Entropy of the subtree generated at {@code c} when {@code c} is forced to take an observed
     * (red) split -- i.e. as a blue region's boundary part: split ~ CCP, children free. Assumes all
     * strictly-smaller clades are already memoised.
     */
    private double entropyRedForced(Clade c, Map<Clade, Double> freeMemo, Map<Clade, Double> redMemo) {
        if (c.isLeaf()) {
            return 0.0;
        }
        Double cached = redMemo.get(c);
        if (cached != null) {
            return cached;
        }
        double h = 0.0;
        for (CladePartition p : c.getPartitions()) {
            double ccp = p.getCCP();
            Clade[] ch = p.getChildClades();
            h += ccp * (-p.getLogCCP()
                    + entropyFree(ch[0], freeMemo, redMemo) + entropyFree(ch[1], freeMemo, redMemo));
        }
        redMemo.put(c, h);
        return h;
    }

    /**
     * The two enumeration-derived parts of a clade's blue (escape) entropy contribution:
     * {@code childEntropy} = sum over regions of P(region) * sum_i H_red(boundary part i), and
     * {@code extraLogPc} = the SHARED-only per-region +log(pathcount) self-information (0 for FLAT).
     */
    private record BlueTerms(double childEntropy, double extraLogPc) {
    }

    /**
     * Enumerates the canonical boundaries of {@code c} (orders 3..lastBoundary) into observed
     * subclades admitting an all-novel resolution -- the same enumeration the reserve uses -- and
     * accumulates the blue-region entropy terms exactly. Bounded by the per-clade op-budget.
     */
    private BlueTerms blueContribution(Clade c, CladeReg reg, double eps,
                                       Map<Clade, Double> freeMemo, Map<Clade, Double> redMemo) {
        List<Clade> subs = observedSubclades(c);
        BitSet cBits = c.getCladeInBits();
        int lastBoundary = reg.lastBoundary();
        double[] acc = new double[2]; // [0] childEntropy, [1] extraLogPc (SHARED)
        enumOps.get()[0] = 0;
        try {
            // Orders 3 and 4 (the default reserve depth) via the fast sum-arithmetic pair pass.
            blueContributionFast(subs, cBits, eps, lastBoundary, acc, freeMemo, redMemo);
            // Deeper orders (only reached when boundaries 3 and 4 are both empty) keep the old walk.
            for (int m = 5; m <= lastBoundary; m++) {
                blueWalk(cBits, subs, m, 0, BitSet.newBitSet(leafArraySize), new ArrayList<>(m),
                        eps, acc, freeMemo, redMemo);
            }
        } catch (BudgetExceeded e) {
            // deeper boundaries omitted (negligible, like the reserve tail); same guard as countNj
        }
        return new BlueTerms(acc[0], acc[1]);
    }

    /**
     * Accumulates the blue-region entropy terms for boundary orders 3 and 4 in one disjoint-pair
     * pass, the entropy counterpart of {@link #countN1N2}: it enumerates the same boundaries by
     * weighted-sum complement lookups, but at each boundary adds {@code boundaryMass * sum of the
     * boundary parts' forced-red entropy} to {@code acc[0]} (and, in SHARED mode, the per-region
     * {@code eps^(m-2) * log pathcount} self-information to {@code acc[1]}). Emits exactly the same
     * boundaries as the old {@link #blueWalk} for these orders.
     */
    private void blueContributionFast(List<Clade> subs, BitSet cBits, double eps, int lastBoundary,
                                      double[] acc, Map<Clade, Double> freeMemo,
                                      Map<Clade, Double> redMemo) {
        int m = subs.size();
        BitSet[] sb = new BitSet[m];
        long[] partSum = new long[m];
        int[] partCard = new int[m];
        double[] redH = new double[m]; // forced-red entropy of each boundary part (memoised lookup)
        for (int i = 0; i < m; i++) {
            sb[i] = subs.get(i).getCladeInBits();
            partSum[i] = weightedSum(sb[i]);
            partCard[i] = sb[i].cardinality();
            redH[i] = entropyRedForced(subs.get(i), freeMemo, redMemo);
        }
        long sumC = weightedSum(cBits);
        int cardC = cBits.cardinality();
        boolean flat = novelMode == NovelMode.FLAT;
        double eps1 = eps;         // eps^(3-2)
        double eps2 = eps * eps;   // eps^(4-2)
        boolean do4 = lastBoundary >= 4;
        final long[] ops = enumOps.get();

        int stb = Integer.highestOneBit(Math.max(1, m)) << 1;
        int smask = stb - 1;
        int[] subHead = new int[stb];
        java.util.Arrays.fill(subHead, -1);
        int[] subNxt = new int[m];
        for (int i = 0; i < m; i++) {
            int b = mixHash(partSum[i]) & smask;
            subNxt[i] = subHead[b];
            subHead[b] = i;
        }

        int cap = 256, P = 0;
        int[] ei = new int[cap], ej = new int[cap], eh = new int[cap];
        long[] sU = new long[cap];
        BitSet ue = BitSet.newBitSet(leafArraySize);
        BitSet uf = BitSet.newBitSet(leafArraySize);
        for (int i = 0; i < m; i++) {
            BitSet a = sb[i];
            for (int j = i + 1; j < m; j++) {
                if (++ops[0] > opsBudget) {
                    throw BUDGET_EXCEEDED;
                }
                if (a.intersects(sb[j])) {
                    continue;
                }
                long s = partSum[i] + partSum[j];
                long sR = sumC - s;
                boolean ueBuilt = false;
                for (int d = subHead[mixHash(sR) & smask]; d >= 0; d = subNxt[d]) {
                    if (d <= j || partSum[d] != sR) {
                        continue;
                    }
                    if (partCard[i] + partCard[j] + partCard[d] != cardC) {
                        continue;
                    }
                    if (!ueBuilt) {
                        ue.clear();
                        ue.or(a);
                        ue.or(sb[j]);
                        ueBuilt = true;
                    }
                    if (ue.intersects(sb[d])) {
                        continue;
                    }
                    BitSet[] parts = {sb[i], sb[j], sb[d]};
                    int pc = countAllNovelResolutions(cBits, parts);
                    if (pc > 0) {
                        double childH = redH[i] + redH[j] + redH[d];
                        acc[0] += (flat ? pc * eps1 : eps1) * childH;
                        if (!flat) {
                            acc[1] += eps1 * Math.log(pc);
                        }
                    }
                }
                if (do4) {
                    if (P == cap) {
                        cap <<= 1;
                        ei = java.util.Arrays.copyOf(ei, cap);
                        ej = java.util.Arrays.copyOf(ej, cap);
                        eh = java.util.Arrays.copyOf(eh, cap);
                        sU = java.util.Arrays.copyOf(sU, cap);
                    }
                    ei[P] = i;
                    ej[P] = j;
                    sU[P] = s;
                    eh[P] = mixHash(s);
                    P++;
                }
            }
        }
        if (!do4 || P == 0) {
            return;
        }
        int tb = Integer.highestOneBit(P) << 1;
        int mask = tb - 1;
        int[] head = new int[tb];
        java.util.Arrays.fill(head, -1);
        int[] nxt = new int[P];
        for (int e = 0; e < P; e++) {
            int b = eh[e] & mask;
            nxt[e] = head[b];
            head[b] = e;
        }
        for (int e = 0; e < P; e++) {
            int i = ei[e], j = ej[e];
            long sR = sumC - sU[e];
            boolean ueBuilt = false;
            for (int f = head[mixHash(sR) & mask]; f >= 0; f = nxt[f]) {
                if (++ops[0] > opsBudget) {
                    throw BUDGET_EXCEEDED;
                }
                int k = ei[f], l = ej[f];
                if (k <= j) {
                    continue;
                }
                if (partCard[i] + partCard[j] + partCard[k] + partCard[l] != cardC) {
                    continue;
                }
                if (!ueBuilt) {
                    ue.clear();
                    ue.or(sb[i]);
                    ue.or(sb[j]);
                    ueBuilt = true;
                }
                uf.clear();
                uf.or(sb[k]);
                uf.or(sb[l]);
                if (ue.intersects(uf)) {
                    continue;
                }
                BitSet[] parts = {sb[i], sb[j], sb[k], sb[l]};
                int pc = countAllNovelResolutions(cBits, parts);
                if (pc > 0) {
                    double childH = redH[i] + redH[j] + redH[k] + redH[l];
                    acc[0] += (flat ? pc * eps2 : eps2) * childH;
                    if (!flat) {
                        acc[1] += eps2 * Math.log(pc);
                    }
                }
            }
        }
    }

    private void blueWalk(BitSet cBits, List<Clade> subs, int m, int startIdx, BitSet used,
                          List<BitSet> chosen, double eps, double[] acc,
                          Map<Clade, Double> freeMemo, Map<Clade, Double> redMemo) {
        if (++enumOps.get()[0] > opsBudget) {
            throw BUDGET_EXCEEDED;
        }
        if (chosen.size() == m - 1) {
            BitSet last = (BitSet) cBits.clone();
            last.andNot(used);
            if (last.cardinality() == 0 || getClade(last) == null) {
                return;
            }
            if (compareBitSets(chosen.get(chosen.size() - 1), last) >= 0) {
                return;
            }
            BitSet[] parts = new BitSet[m];
            for (int i = 0; i < m - 1; i++) {
                parts[i] = chosen.get(i);
            }
            parts[m - 1] = last;
            int pc = countAllNovelResolutions(cBits, parts);
            if (pc <= 0) {
                return;
            }
            double epsPow = Math.pow(eps, m - EXP_REDUCTION);
            double sumChildH = 0.0;
            for (BitSet pb : parts) {
                sumChildH += entropyRedForced(getClade(pb), freeMemo, redMemo); // memoised lookup
            }
            // mass on this boundary: pathcount * eps^(m-2) [FLAT] or eps^(m-2) [SHARED]
            double boundaryMass = (novelMode == NovelMode.FLAT) ? pc * epsPow : epsPow;
            acc[0] += boundaryMass * sumChildH;
            if (novelMode == NovelMode.SHARED) {
                acc[1] += epsPow * Math.log(pc);
            }
            return;
        }
        for (int i = startIdx; i < subs.size(); i++) {
            BitSet pb = subs.get(i).getCladeInBits();
            if (pb.intersects(used)) {
                continue;
            }
            chosen.add(pb);
            BitSet newUsed = (BitSet) used.clone();
            newUsed.or(pb);
            blueWalk(cBits, subs, m, i + 1, newUsed, chosen, eps, acc, freeMemo, redMemo);
            chosen.remove(chosen.size() - 1);
        }
    }

    /* ----------------------------------------------------------------------
     * Calibration: total model mass on trees with at least one novel clade
     * ------------------------------------------------------------------- */

    /**
     * The model's total probability mass on trees containing at least one novel
     * clade, {@code P(novel) = 1 - M(root)}, where
     * {@code M(c) = disc(c) * sum_{observed (L,R)} ccp(L,R) * M(L) * M(R)} is the
     * mass on subtrees of {@code c} that stay entirely within the observed clade
     * set, and {@code disc(c) = 1 - mu} for clades that reserve escape mass (the
     * tail correction is negligible for this aggregate). Compare against the
     * empirical fraction of held-out trees that {@link #containsNovelClade}.
     *
     * @return model probability that a tree contains a novel clade
     */
    public double getNovelCladeProbability() {
        return getNovelCladeProbability(mu);
    }

    /**
     * Closed-form P(tree contains a novel clade) at an arbitrary escape probability {@code atMu},
     * reusing this model's (mu-independent) alpha-smoothed backbone CCPs and reservability. This is the
     * P(novel) analogue of {@link #getLogProbabilityOfTree(Tree, double)}: it lets a sweep evaluate many
     * mu on one trained backbone without rebuilding. For {@code atMu == mu} it equals the no-arg form.
     */
    public double getNovelCladeProbability(double atMu) {
        return 1.0 - massNoNovelClade(getRootClade(), new HashMap<>(), atMu);
    }

    private double massNoNovelClade(Clade c, Map<Clade, Double> memo, double atMu) {
        if (c.isLeaf()) {
            return 1.0;
        }
        Double cached = memo.get(c);
        if (cached != null) {
            return cached;
        }
        double sum = 0.0;
        for (CladePartition p : c.getPartitions()) {
            Clade[] ch = p.getChildClades();
            sum += p.getCCP() * massNoNovelClade(ch[0], memo, atMu) * massNoNovelClade(ch[1], memo, atMu);
        }
        double m = (cheapReservable(c) ? (1.0 - atMu) : 1.0) * sum;
        memo.put(c, m);
        return m;
    }

    /**
     * Cheap reservability test (early-exit over the reserve-depth orders), used
     * for the aggregate {@link #getNovelCladeProbability} where the exact tail is
     * not needed.
     */
    private boolean cheapReservable(Clade c) {
        if (c.size() < 3) {
            return false;
        }
        List<Clade> subs = observedSubclades(c);
        enumOps.get()[0] = 0;
        try {
            for (int j = 3; j <= reserveBoundary && j <= c.size(); j++) {
                if (countNj(c, subs, j, true) > 0) {
                    return true;
                }
            }
        } catch (BudgetExceeded e) {
            // boundary 3 (O(m^2)) always completes
        }
        return false;
    }

    /** Whether the given tree contains a clade that is not in the (expanded) CCD. */
    public boolean containsNovelClade(Tree tree) {
        Map<Node, BitSet> bits = new HashMap<>();
        computeBits(tree.getRoot(), bits);
        for (Map.Entry<Node, BitSet> e : bits.entrySet()) {
            if (!e.getKey().isLeaf() && getClade(e.getValue()) == null) {
                return true;
            }
        }
        return false;
    }

    private BitSet computeBits(Node v, Map<Node, BitSet> bits) {
        BitSet b = BitSet.newBitSet(leafArraySize);
        if (v.isLeaf()) {
            b.set(v.getNr());
        } else {
            b.or(computeBits(v.getChildren().get(0), bits));
            b.or(computeBits(v.getChildren().get(1), bits));
        }
        bits.put(v, b);
        return b;
    }

    /**
     * Scores the subtree at {@code v} at escape probability {@code scoreMu}.
     * {@code useStoredReg} controls the source of the per-clade reserve: when {@code true}
     * (the {@code scoreMu == mu} path) it uses the cached {@code eps}/{@code tail} solved at
     * construction; when {@code false} it re-solves {@code eps} from the cached {@code N_j} at
     * {@code scoreMu} with tail = 0 (the held-out / mu-search path).
     */
    private void scoreFresh(Node v, Map<Node, BitSet> bits, double[] logp,
                            double scoreMu, boolean useStoredReg) {
        if (v.isLeaf() || logp[0] == Double.NEGATIVE_INFINITY) {
            return;
        }
        Node c1 = v.getChildren().get(0);
        Node c2 = v.getChildren().get(1);
        BitSet vb = bits.get(v);

        if (isRed(vb, bits.get(c1), bits.get(c2))) {
            Clade c = getClade(vb);
            CladeReg reg = computeReg(c);
            if (reg.reservable()) {
                // tail = 0 for NONE / the mu-search path -> plain (1 - mu); otherwise (1 - mu - tail)
                double tail = useStoredReg ? reg.tail() : 0.0;
                logp[0] += Math.log(1.0 - scoreMu - tail);
            }
            logp[0] += rawLogCCP(c, bits.get(c1), bits.get(c2));
            scoreFresh(c1, bits, logp, scoreMu, useStoredReg);
            scoreFresh(c2, bits, logp, scoreMu, useStoredReg);
        } else {
            // v is the top of a maximal blue region. No truncation: every region
            // is scored, so any tree gets nonzero probability (full support).
            List<Node> boundary = new ArrayList<>();
            collectRegion(v, bits, boundary);
            int m = boundary.size();
            Clade c = getClade(vb); // region top is always an observed clade
            BitSet[] parts = new BitSet[m];
            for (int i = 0; i < m; i++) {
                parts[i] = bits.get(boundary.get(i));
            }
            // SHARED shares a boundary's eps^(m-2) over its pathcount resolutions (-> -log pathcount);
            // FLAT prices each distinct resolution eps^(m-2) directly (no pathcount normaliser).
            double logPath = (novelMode == NovelMode.FLAT)
                    ? 0.0 : Math.log(countAllNovelResolutions(vb, parts));
            double logEps = useStoredReg ? computeReg(c).logEps() : logEpsAtMu(c, scoreMu);
            if (logEps == Double.NEGATIVE_INFINITY) {
                // C reserves nothing computable: give this region fallback mass
                // mu (SHARED: shared over its pathcount resolutions)
                logp[0] += Math.log(scoreMu) - logPath;
            } else {
                // m - 2 = number of novel clades in this region: one log eps each
                logp[0] += (m - EXP_REDUCTION) * logEps - logPath;
            }
            for (Node b : boundary) {
                scoreFresh(b, bits, logp, scoreMu, useStoredReg);
            }
        }
    }

    /**
     * Re-solves {@code log eps(C)} from the cached ({@code mu}-independent) {@code N_j} counts at
     * an arbitrary {@code scoreMu}; {@code NEGATIVE_INFINITY} if {@code C} reserves nothing.
     */
    private double logEpsAtMu(Clade c, double scoreMu) {
        CladeReg reg = computeReg(c);
        if (!reg.reservable()) {
            return Double.NEGATIVE_INFINITY;
        }
        return solveLogEps(reg.nCounts(), reg.lastBoundary(), scoreMu);
    }

    private void collectRegion(Node v, Map<Node, BitSet> bits, List<Node> boundary) {
        for (Node child : v.getChildren()) {
            if (child.isLeaf()) {
                boundary.add(child);
            } else if (isRed(bits.get(child), bits.get(child.getChildren().get(0)),
                    bits.get(child.getChildren().get(1)))) {
                boundary.add(child);
            } else {
                collectRegion(child, bits, boundary);
            }
        }
    }

    /* ----------------------------------------------------------------------
     * Tree sampling (full support)
     *
     * The inherited AbstractCCD sampler only ever picks an observed clade partition, so it
     * draws from the truncated red-only distribution and can never produce a novel clade. We
     * override the recursion so that at each reservable observed clade we either take an
     * observed (red) split with probability 1 - mu, or escape (probability mu, or mu + tail in
     * FULL_SUPPORT) into a blue region: a binary resolution of the clade through novel
     * intermediate clades down to a boundary of observed subclades. The escape draw reproduces
     * the score's blue-region weight exactly:
     *   P(region) = escapeMass * (N_j eps^j / escapeMass) * (1/N_j) * (1/pathcount)
     *             = eps^(m-2) / pathcount.
     * Boundary parts are resolved red (an observed split), which makes regions maximal and the
     * decomposition unique, matching scoreFresh's collectRegion. The PROB/LOG_PROB metadata is
     * stamped with the same factors scoreFresh charges, so the materialised tree's stamped
     * (log-)probability equals getLogProbabilityOfTree of that tree.
     * ------------------------------------------------------------------- */

    /**
     * Selects how faithfully {@link #sampleTree()} reproduces the full-support score; see
     * {@link SamplingFidelity}. Default {@link SamplingFidelity#SELF_CONSISTENT}.
     */
    public void setSamplingFidelity(SamplingFidelity fidelity) {
        this.samplingFidelity = fidelity;
    }

    public SamplingFidelity getSamplingFidelity() {
        return samplingFidelity;
    }

    @Override
    protected Node getVertexBasedOnStrategy(Clade clade, SamplingStrategy samplingStrategy,
                                            HeightSettingStrategy heightStrategy) {
        // Escape sampling is meaningful only for random Sampling; point estimates (MAP,
        // MaxSumCladeCredibility) and leaves keep the inherited observed-only behaviour.
        if (samplingStrategy != SamplingStrategy.Sampling || clade.isLeaf()) {
            return super.getVertexBasedOnStrategy(clade, samplingStrategy, heightStrategy);
        }

        CladeReg reg = computeReg(clade);
        if (reg.reservable()) {
            double eps = Math.exp(reg.logEps());
            double[] orderWeight = escapeOrderWeights(clade, reg, eps);
            double escapeMass = 0.0;
            for (double w : orderWeight) {
                escapeMass += w;
            }
            if (escapeMass > 0 && random.nextDouble() < escapeMass) {
                Node region = sampleBlueRegion(clade, reg, orderWeight, escapeMass, heightStrategy);
                if (region != null) {
                    return region;
                }
                // region sampling hit the op-budget; fall through to an observed split
            }
        }
        return resolveRedRoot(clade, samplingStrategy, heightStrategy);
    }

    /**
     * Resolves {@code clade} with an observed (red) split and recurses into its children
     * through the overriding sampler (so descendant clades may still escape). Stamps the same
     * {@code (1 - mu - tail) * theta} factor scoreFresh charges at a red node. Leaves delegate
     * to the inherited leaf construction.
     */
    private Node resolveRedRoot(Clade clade, SamplingStrategy samplingStrategy,
                                HeightSettingStrategy heightStrategy) {
        if (clade.isLeaf()) {
            return super.getVertexBasedOnStrategy(clade, samplingStrategy, heightStrategy);
        }
        CladePartition partition = getPartitionBasedOnStrategy(clade, samplingStrategy);
        if (partition == null) {
            throw new AssertionError("Unsuccessful to find clade partition of clade: "
                    + clade.getCladeInBits());
        }
        Node firstChild = getVertexBasedOnStrategy(partition.getChildClades()[0],
                samplingStrategy, heightStrategy);
        Node secondChild = getVertexBasedOnStrategy(partition.getChildClades()[1],
                samplingStrategy, heightStrategy);
        Node vertex = buildInternalVertex(clade, partition, firstChild, secondChild, heightStrategy);
        CladeReg reg = computeReg(clade);
        if (reg.reservable()) {
            applyLogFactor(vertex, Math.log(1.0 - mu - reg.tail()));
        }
        return vertex;
    }

    /**
     * Escape weight per boundary size: {@code w[m] = N_m * eps^(m-2)}. SELF_CONSISTENT uses the
     * computed orders (whose sum is mu); FULL_SUPPORT additionally Knuth-estimates the leading
     * tail orders ({@code k+1, k+2}).
     */
    private double[] escapeOrderWeights(Clade c, CladeReg reg, double eps) {
        int maxBoundary = reg.lastBoundary();
        int[] n = reg.nCounts();
        int hi = maxBoundary;
        if (samplingFidelity == SamplingFidelity.FULL_SUPPORT) {
            hi = Math.min(c.size(), reserveBoundary + 2);
        }
        double[] w = new double[Math.max(maxBoundary, hi) + 1];
        for (int m = 3; m <= maxBoundary; m++) {
            if (n[m] > 0) {
                w[m] = n[m] * Math.pow(eps, m - EXP_REDUCTION);
            }
        }
        if (samplingFidelity == SamplingFidelity.FULL_SUPPORT && hi > maxBoundary) {
            List<Clade> subs = observedSubclades(c);
            for (int m = maxBoundary + 1; m <= hi; m++) {
                double nm = estimateNj(c, subs, m, TAIL_SAMPLES);
                if (nm > 0) {
                    w[m] = nm * Math.pow(eps, m - EXP_REDUCTION);
                }
            }
        }
        return w;
    }

    /**
     * Samples a maximal blue region rooted at observed clade {@code c}: draws a boundary size,
     * a uniform valid boundary, a uniform all-novel resolution shape, resolves the boundary parts
     * red, and stamps the region factor {@code eps^(m-2)/pathcount}. Returns {@code null} if the
     * boundary enumeration exceeds the op-budget (the caller then takes an observed split).
     */
    private Node sampleBlueRegion(Clade c, CladeReg reg, double[] orderWeight, double escapeMass,
                                  HeightSettingStrategy heightStrategy) {
        double target = random.nextDouble() * escapeMass;
        double acc = 0.0;
        int m = -1;
        for (int mm = 3; mm < orderWeight.length; mm++) {
            if (orderWeight[mm] <= 0) {
                continue;
            }
            acc += orderWeight[mm];
            if (target < acc) {
                m = mm;
                break;
            }
        }
        if (m < 0) { // numerical guard: fall back to the largest available order
            for (int mm = orderWeight.length - 1; mm >= 3; mm--) {
                if (orderWeight[mm] > 0) {
                    m = mm;
                    break;
                }
            }
        }
        if (m < 0) {
            return null;
        }

        List<Clade> subs = observedSubclades(c);
        BitSet[] parts;
        try {
            parts = sampleBoundary(c, subs, m);
        } catch (BudgetExceeded e) {
            return null;
        }
        if (parts == null) {
            return null;
        }

        Resolution res = sampleResolution(c.getCladeInBits(), parts);

        // Resolve every boundary part with an observed split (red) so the region is maximal.
        Node[] boundaryNodes = new Node[m];
        for (int i = 0; i < m; i++) {
            Clade part = getClade(parts[i]);
            boundaryNodes[i] = resolveRedRoot(part, SamplingStrategy.Sampling, heightStrategy);
        }

        double topHeight = regionTopHeight(c, heightStrategy);
        Node region = assembleRegion(res.shape(), boundaryNodes, c, topHeight, heightStrategy);
        // One factor of eps per novel clade (m - 2 of them). SHARED spreads it over the pathcount
        // resolutions (-log pathcount); FLAT prices the resolution directly (boundary already drawn
        // proportional to pathcount above), matching the FLAT scorer.
        double regionFactor = (m - EXP_REDUCTION) * reg.logEps();
        if (novelMode == NovelMode.SHARED) {
            regionFactor -= Math.log(res.pathcount());
        }
        applyLogFactor(region, regionFactor);
        return region;
    }

    /**
     * Uniformly samples one valid boundary of {@code c} into {@code m} observed subclades
     * (admitting an all-novel resolution) by reservoir sampling over the same canonical
     * enumeration {@link #countNj} counts. Returns the chosen parts, or {@code null} if none.
     */
    private BitSet[] sampleBoundary(Clade c, List<Clade> subs, int m) {
        boundarySeen = 0;
        boundaryWeightSeen = 0.0;
        boundaryPick = null;
        enumOps.get()[0] = 0;
        sampleBoundaryWalk(c.getCladeInBits(), subs, m, 0,
                BitSet.newBitSet(leafArraySize), new ArrayList<>(m));
        return boundaryPick;
    }

    private void sampleBoundaryWalk(BitSet cBits, List<Clade> subs, int m, int startIdx,
                                    BitSet used, List<BitSet> chosen) {
        if (++enumOps.get()[0] > opsBudget) {
            throw BUDGET_EXCEEDED;
        }
        if (chosen.size() == m - 1) {
            BitSet last = (BitSet) cBits.clone();
            last.andNot(used);
            if (last.cardinality() == 0 || getClade(last) == null) {
                return;
            }
            if (compareBitSets(chosen.get(chosen.size() - 1), last) >= 0) {
                return;
            }
            BitSet[] parts = new BitSet[m];
            for (int i = 0; i < m - 1; i++) {
                parts[i] = chosen.get(i);
            }
            parts[m - 1] = last;
            int pc = countAllNovelResolutions(cBits, parts);
            if (pc <= 0) {
                return;
            }
            if (novelMode == NovelMode.FLAT) {
                // FLAT samples a boundary in proportion to its pathcount (so that, with the
                // uniform resolution pick below, every distinct novel tree is equiprobable).
                boundaryWeightSeen += pc;
                if (random.nextDouble() * boundaryWeightSeen < pc) { // weighted reservoir
                    boundaryPick = parts;
                }
            } else {
                boundarySeen++;
                if (random.nextInt(boundarySeen) == 0) { // uniform reservoir: keep with prob 1/seen
                    boundaryPick = parts;
                }
            }
            return;
        }
        for (int i = startIdx; i < subs.size(); i++) {
            BitSet pb = subs.get(i).getCladeInBits();
            if (pb.intersects(used)) {
                continue;
            }
            chosen.add(pb);
            BitSet newUsed = (BitSet) used.clone();
            newUsed.or(pb);
            sampleBoundaryWalk(cBits, subs, m, i + 1, newUsed, chosen);
            chosen.remove(chosen.size() - 1);
        }
    }

    /**
     * A sampled all-novel resolution of a clade into its boundary parts, with the total number
     * of such resolutions (the pathcount) that normalises the region.
     */
    private record Resolution(Shape shape, int pathcount) {
    }

    /**
     * Node of a sampled resolution shape: a boundary part (leaf, {@code partIndex >= 0}) or a
     * novel internal split with two children.
     */
    private static final class Shape {
        final BitSet union;
        final int partIndex;
        final Shape left;
        final Shape right;

        Shape(BitSet union, int partIndex) {
            this.union = union;
            this.partIndex = partIndex;
            this.left = null;
            this.right = null;
        }

        Shape(BitSet union, Shape left, Shape right) {
            this.union = union;
            this.partIndex = -1;
            this.left = left;
            this.right = right;
        }

        boolean isLeaf() {
            return partIndex >= 0;
        }
    }

    /**
     * Uniformly samples one all-novel binary resolution of {@code cBits} into {@code parts} by
     * top-down sampling from the same subset DP {@link #countAllNovelResolutions} builds (at each
     * mask a split is chosen with probability proportional to {@code f[s1]*f[s2]}).
     */
    private Resolution sampleResolution(BitSet cBits, BitSet[] parts) {
        int k = parts.length;
        if (k == 1) {
            return new Resolution(new Shape(parts[0], 0), 1);
        }
        int full = (1 << k) - 1;
        BitSet[] unionOf = new BitSet[1 << k];
        unionOf[0] = BitSet.newBitSet(leafArraySize);
        for (int mask = 1; mask <= full; mask++) {
            int low = Integer.numberOfTrailingZeros(mask);
            BitSet u = (BitSet) unionOf[mask & (mask - 1)].clone();
            u.or(parts[low]);
            unionOf[mask] = u;
        }
        int[] f = new int[1 << k];
        for (int mask = 1; mask <= full; mask++) {
            if (Integer.bitCount(mask) == 1) {
                f[mask] = 1;
                continue;
            }
            int low = mask & (-mask);
            int rest = mask ^ low;
            int count = 0;
            for (int sub = rest; ; sub = (sub - 1) & rest) {
                int s1 = sub | low;
                int s2 = mask ^ s1;
                if (s2 != 0 && !isSplitObserved(unionOf[mask], unionOf[s1], unionOf[s2])) {
                    count += f[s1] * f[s2];
                }
                if (sub == 0) {
                    break;
                }
            }
            f[mask] = count;
        }
        Shape shape = sampleShape(full, parts, unionOf, f);
        return new Resolution(shape, f[full]);
    }

    private Shape sampleShape(int mask, BitSet[] parts, BitSet[] unionOf, int[] f) {
        if (Integer.bitCount(mask) == 1) {
            int idx = Integer.numberOfTrailingZeros(mask);
            return new Shape(parts[idx], idx);
        }
        int low = mask & (-mask);
        int rest = mask ^ low;
        int targetCount = random.nextInt(f[mask]); // f[mask] > 0 at every visited mask
        int acc = 0;
        int chosenS1 = -1;
        int chosenS2 = -1;
        for (int sub = rest; ; sub = (sub - 1) & rest) {
            int s1 = sub | low;
            int s2 = mask ^ s1;
            if (s2 != 0 && !isSplitObserved(unionOf[mask], unionOf[s1], unionOf[s2])) {
                acc += f[s1] * f[s2];
                if (targetCount < acc) {
                    chosenS1 = s1;
                    chosenS2 = s2;
                    break;
                }
            }
            if (sub == 0) {
                break;
            }
        }
        Shape left = sampleShape(chosenS1, parts, unionOf, f);
        Shape right = sampleShape(chosenS2, parts, unionOf, f);
        return new Shape((BitSet) unionOf[mask].clone(), left, right);
    }

    /**
     * Materialises a sampled resolution shape into nodes: boundary parts reuse their already
     * resolved (red) nodes, novel internal nodes get fresh nodes with identity probability factors
     * (the region factor is applied once at the top by the caller) and interpolated heights.
     */
    private Node assembleRegion(Shape shape, Node[] boundaryNodes, Clade c, double topHeight,
                                HeightSettingStrategy heightStrategy) {
        if (shape.isLeaf()) {
            return boundaryNodes[shape.partIndex];
        }
        Node left = assembleRegion(shape.left, boundaryNodes, c, topHeight, heightStrategy);
        Node right = assembleRegion(shape.right, boundaryNodes, c, topHeight, heightStrategy);

        Node vertex = new Node();
        vertex.setNr(nextRunningInnerIndex());
        vertex.addChild(left);
        vertex.addChild(right);

        boolean isTop = shape.union.equals(c.getCladeInBits());
        Clade obs = getClade(shape.union); // observed at the top (= c) and at any observed interior
        double support = (obs != null) ? obs.getProbability() : 0.0;
        vertex.setMetaData(CLADE_SUPPORT_KEY, support);
        String posteriorSupport = CLADE_SUPPORT_KEY + "=" + support;
        vertex.metaDataString = (vertex.metaDataString != null)
                ? vertex.metaDataString + "," + posteriorSupport : posteriorSupport;

        // identity factors: the per-region eps^(m-2)/pathcount is stamped once at the region top.
        vertex.setMetaData(PROB_SUBTREE_KEY,
                (Double) left.getMetaData(PROB_SUBTREE_KEY) * (Double) right.getMetaData(PROB_SUBTREE_KEY));
        vertex.setMetaData(LOG_PROB_SUBTREE_KEY,
                (Double) left.getMetaData(LOG_PROB_SUBTREE_KEY) + (Double) right.getMetaData(LOG_PROB_SUBTREE_KEY));

        setNovelHeight(vertex, left, right, isTop, topHeight, heightStrategy);
        return vertex;
    }

    /**
     * Region-top height under the given strategy ({@code NaN} for One/None, which the novel-height
     * scheme handles without a fixed top).
     */
    private double regionTopHeight(Clade c, HeightSettingStrategy heightStrategy) {
        if (heightStrategy == HeightSettingStrategy.MeanOccurredHeights) {
            return c.getMeanOccurredHeight();
        }
        if (heightStrategy == HeightSettingStrategy.CommonAncestorHeights) {
            return c.getCommonAncestorHeight();
        }
        return Double.NaN;
    }

    /**
     * Sets a height for a novel (or region-top) internal node that keeps branch lengths positive:
     * One uses {@code max(child) + 1}; Mean/CommonAncestor put the top at its clade height and each
     * interior node at the midpoint between its children and the top, falling back to
     * {@code max(child) + eps} if the top is not strictly above the children; None leaves it unset.
     */
    private void setNovelHeight(Node vertex, Node left, Node right, boolean isTop,
                                double topHeight, HeightSettingStrategy heightStrategy) {
        if (heightStrategy == null || heightStrategy == HeightSettingStrategy.None) {
            return;
        }
        double maxChild = Math.max(left.getHeight(), right.getHeight());
        double height;
        if (heightStrategy == HeightSettingStrategy.One) {
            height = maxChild + 1.0;
        } else if (isTop) {
            height = (topHeight > maxChild) ? topHeight : maxChild + NOVEL_HEIGHT_EPS;
        } else {
            height = (topHeight > maxChild) ? 0.5 * (maxChild + topHeight) : maxChild + NOVEL_HEIGHT_EPS;
        }
        vertex.setHeight(height);
    }

    private static void applyLogFactor(Node vertex, double logFactor) {
        vertex.setMetaData(LOG_PROB_SUBTREE_KEY,
                (Double) vertex.getMetaData(LOG_PROB_SUBTREE_KEY) + logFactor);
        vertex.setMetaData(PROB_SUBTREE_KEY,
                (Double) vertex.getMetaData(PROB_SUBTREE_KEY) * Math.exp(logFactor));
    }

    /* ----------------------------------------------------------------------
     * Observed-split / observed-clade queries
     * ------------------------------------------------------------------- */

    private boolean isSplitObserved(BitSet parentBits, BitSet c1Bits, BitSet c2Bits) {
        Clade parent = getClade(parentBits);
        if (parent == null) {
            return false;
        }
        Clade a = getClade(c1Bits);
        if (a == null) {
            return false;
        }
        Clade b = getClade(c2Bits);
        if (b == null) {
            return false;
        }
        return parent.getCladePartition(a, b) != null;
    }

    private boolean isRed(BitSet parentBits, BitSet c1Bits, BitSet c2Bits) {
        return isSplitObserved(parentBits, c1Bits, c2Bits);
    }

    private double rawLogCCP(Clade parent, BitSet c1Bits, BitSet c2Bits) {
        Clade a = getClade(c1Bits);
        Clade b = getClade(c2Bits);
        CladePartition p = parent.getCladePartition(a, b);
        return p.getLogCCP();
    }

    /* ----------------------------------------------------------------------
     * pathcount: number of all-novel resolutions of C into the given parts,
     * via a subset DP over the parts (handles any boundary size).
     * ------------------------------------------------------------------- */

    private int countAllNovelResolutions(BitSet cBits, BitSet[] parts) {
        int k = parts.length;
        if (k == 1) {
            return 1;
        }
        int full = (1 << k) - 1;
        BitSet[] unionOf = new BitSet[1 << k];
        unionOf[0] = BitSet.newBitSet(leafArraySize);
        for (int mask = 1; mask <= full; mask++) {
            int low = Integer.numberOfTrailingZeros(mask);
            BitSet u = (BitSet) unionOf[mask & (mask - 1)].clone();
            u.or(parts[low]);
            unionOf[mask] = u;
        }
        int[] f = new int[1 << k];
        for (int mask = 1; mask <= full; mask++) {
            if (Integer.bitCount(mask) == 1) {
                f[mask] = 1;
                continue;
            }
            int low = mask & (-mask); // lowest set bit, fixed in S1 to avoid double counting
            int rest = mask ^ low;
            int count = 0;
            // iterate non-empty proper submasks S1 of mask that contain `low`
            for (int sub = rest; ; sub = (sub - 1) & rest) {
                int s1 = sub | low;
                int s2 = mask ^ s1;
                if (s2 != 0 && !isSplitObserved(unionOf[mask], unionOf[s1], unionOf[s2])) {
                    count += f[s1] * f[s2];
                }
                if (sub == 0) {
                    break;
                }
            }
            f[mask] = count;
        }
        return f[full];
    }

    /**
     * Eagerly compute and cache the per-clade reserve eps for every clade in this CCD. The reserve solve
     * runs the op-capped N_j enumeration once per clade ({@link #computeReg}); doing it up front means
     * subsequent {@code sampleTree}/{@code getLogProbabilityOfTree} calls hit {@code regCache} instead of
     * paying that enumeration lazily — and unevenly — during the first samples. Idempotent (cache hits
     * after the first call). Useful before a large sampling loop (e.g. variational ELBO estimation).
     */
    public void precomputeReserves() {
        // Clades are independent: each computeReg only reads immutable clade structure + writes its own
        // entry in the (concurrent) regCache, and the op-budget counter is thread-local. With tailMode
        // BOUND/NONE there is no random tail estimation, so this is deterministic and thread-safe.
        cladeMapping.values().parallelStream().forEach(this::computeReg);
    }

    /** The per-clade enumeration op budget: an explicit {@code -Dkreg.enumOps} override (accepts
     *  {@code 1e7} or {@code 10000000}), else {@link #OPS_BUDGET_FACTOR} times the largest clade's N_1
     *  cost (root subclade count squared), floored at {@link #MIN_OPS_BUDGET}. Scales with dataset size
     *  so N_1 always completes while deep orders stay capped. */
    private long computeOpsBudget() {
        String override = System.getProperty("kreg.enumOps");
        if (override != null && !override.isBlank()) {
            return (long) Double.parseDouble(override.trim());
        }
        long m = observedSubclades(getRootClade()).size();   // largest clade's observed-subclade count
        return Math.max(MIN_OPS_BUDGET, OPS_BUDGET_FACTOR * m * m);
    }

    /** The per-clade enumeration op budget in effect (see {@link #computeOpsBudget}). */
    public long getOpsBudget() {
        return opsBudget;
    }

    /* ----------------------------------------------------------------------
     * Per-clade reservation: N_j and eps
     * ------------------------------------------------------------------- */

    /* Per-clade reserve: solve eps from R(C) = sum_{j>=1} N_{j+2} eps^j = mu using
     * orders up to the reserve depth (boundaries 3..reserveBoundary; climbing to the
     * first nonzero order if those are all zero, so eps stays positive), and bound
     * the omitted tail (orders > k) geometrically from the two computed orders:
     *   rho = (N_2/N_1) eps,  tail <= N_2 eps^2 * rho/(1-rho)
     * (rigorous if the N_j are log-concave; clamped to [0, mu]). Computed once per
     * clade and cached. */
    private CladeReg computeReg(Clade c) {
        CladeReg cached = regCache.get(c);
        if (cached != null) {
            return cached;
        }
        CladeReg reg;
        if (c.size() < 3) {
            reg = new CladeReg(false, Double.NEGATIVE_INFINITY, 0.0, null, 0);
            regCache.put(c, reg);
            return reg;
        }
        List<Clade> subs = observedSubclades(c);
        int[] n = new int[c.size() + 1];
        int last = 2;
        boolean anyNonzero = false;
        enumOps.get()[0] = 0;
        if (reserveBoundary == 3 || reserveBoundary == 4) {
            // Fast path for the practical reserve depths: k = 1 (N_1 only) or k = 2 (N_1 and N_2),
            // sharing one weighted-sum disjoint-pair pass. N_2 is skipped entirely when k = 1.
            boolean computeN2 = reserveBoundary == 4 && c.size() >= 4;
            try {
                long[] n12 = countN1N2(c, subs, computeN2);
                n[3] = (int) n12[0];
                last = 3;
                if (n[3] > 0) {
                    anyNonzero = true;
                }
                if (computeN2) {
                    n[4] = (int) n12[1];
                    last = 4;
                    if (n[4] > 0) {
                        anyNonzero = true;
                    }
                }
            } catch (BudgetExceeded e) {
                // pathological pair blow-up: leave N_1/N_2 = 0 (same capped-out reserve as before)
            }
            // Only climb past the reserve depth if all computed orders were empty (keeps eps positive).
            for (int j = last + 1; j <= c.size() && !anyNonzero; j++) {
                try {
                    n[j] = countNj(c, subs, j, false);
                } catch (BudgetExceeded e) {
                    break;
                }
                last = j;
                if (n[j] > 0) {
                    anyNonzero = true;
                }
            }
        } else {
            for (int j = 3; j <= c.size(); j++) {
                if (j > reserveBoundary && anyNonzero) {
                    break; // reached reserve depth and have a usable (positive) reserve
                }
                try {
                    n[j] = countNj(c, subs, j, false); // accumulates enumOps across j
                } catch (BudgetExceeded e) {
                    break;
                }
                last = j;
                if (n[j] > 0) {
                    anyNonzero = true;
                }
            }
        }
        double logEps = anyNonzero ? solveLogEps(n, last, mu) : Double.NEGATIVE_INFINITY;

        double tail = 0.0;
        int n1 = n[3];
        int n2 = (n.length > 4) ? n[4] : 0;
        if (tailMode != TailMode.NONE && logEps != Double.NEGATIVE_INFINITY && n1 > 0 && n2 > 0) {
            double eps = Math.exp(logEps);
            // geometric upper bound from the two computed orders (rigorous if log-concave)
            double rho = ((double) n2 / n1) * eps;
            double tailUB = (rho > 0 && rho < 1) ? n2 * eps * eps * rho / (1 - rho) : mu;
            tailUB = Math.min(tailUB, mu);

            if (tailMode == TailMode.BOUND || tailUB < TAIL_SAMPLE_THRESHOLD) {
                // BOUND mode, or a negligible bound: use the (cheap) bound directly
                tail = tailUB;
            } else {
                // SAMPLED mode with a non-negligible bound: Knuth-estimate the actual
                // tail (orders > reserve depth), summing the two leading orders
                tail = 0.0;
                for (int m = reserveBoundary + 1; m <= reserveBoundary + 2 && m <= c.size(); m++) {
                    double nm = estimateNj(c, subs, m, TAIL_SAMPLES);
                    tail += nm * Math.pow(eps, m - EXP_REDUCTION);
                }
                tail = Math.min(tail, mu);
            }
        }

        reg = new CladeReg(anyNonzero, logEps, tail, anyNonzero ? n : null, last);
        regCache.put(c, reg);
        return reg;
    }

    /* Knuth (1975) backtrack-tree estimator of N_m (number of valid boundary-m
     * partitions): average over random root-to-leaf descents of the product of
     * branching factors, times the leaf-valid indicator. Each descent is O(m * |subs|);
     * unbiased. */
    private double estimateNj(Clade c, List<Clade> subs, int m, int samples) {
        BitSet cBits = c.getCladeInBits();
        double sum = 0.0;
        for (int s = 0; s < samples; s++) {
            sum += knuthDescent(cBits, subs, m);
        }
        return sum / samples;
    }

    private double knuthDescent(BitSet cBits, List<Clade> subs, int m) {
        BitSet used = BitSet.newBitSet(leafArraySize);
        double product = 1.0;
        int startIdx = 0;
        BitSet prev = null;
        BitSet[] chosen = new BitSet[m];
        for (int level = 0; level < m - 1; level++) {
            // children = candidates at index >= startIdx disjoint from `used`
            int d = 0;
            for (int i = startIdx; i < subs.size(); i++) {
                if (!subs.get(i).getCladeInBits().intersects(used)) {
                    d++;
                }
            }
            if (d == 0) {
                return 0.0; // dead end
            }
            product *= d;
            int target = rng.nextInt(d);
            int pickIdx = -1, seen = 0;
            for (int i = startIdx; i < subs.size(); i++) {
                if (!subs.get(i).getCladeInBits().intersects(used)) {
                    if (seen == target) {
                        pickIdx = i;
                        break;
                    }
                    seen++;
                }
            }
            BitSet pb = subs.get(pickIdx).getCladeInBits();
            chosen[level] = pb;
            used = (BitSet) used.clone();
            used.or(pb);
            prev = pb;
            startIdx = pickIdx + 1;
        }
        BitSet last = (BitSet) cBits.clone();
        last.andNot(used);
        if (last.cardinality() == 0 || getClade(last) == null) {
            return 0.0;
        }
        if (compareBitSets(prev, last) >= 0) {
            return 0.0; // canonical order violated
        }
        chosen[m - 1] = last;
        return countAllNovelResolutions(cBits, chosen) > 0 ? product : 0.0;
    }

    /* Observed clades strictly contained in C, sorted in canonical (lexicographic) order. */
    private List<Clade> observedSubclades(Clade c) {
        List<Clade> subs = new ArrayList<>();
        int cSize = c.size();
        for (Clade x : getClades()) {
            int s = x.size();
            if (s >= 1 && s < cSize && c.contains(x.getCladeInBitsTaxaOnly())) {
                subs.add(x);
            }
        }
        subs.sort((a, b) -> compareBitSets(a.getCladeInBits(), b.getCladeInBits()));
        return subs;
    }

    /* Count m-part boundaries of C into observed subclades admitting an all-novel
     * resolution. earlyExit returns 1 as soon as one is found. */
    private int countNj(Clade c, List<Clade> subs, int m, boolean earlyExit) {
        BitSet cBits = c.getCladeInBits();
        if (m == 4) {
            // Boundary-4 (N_2) dominates construction: the generic enumerate is O(m^3) (choose 3
            // parts, derive the 4th as the complement), so on the big clades it burns the whole
            // op-budget. A 4-part boundary is exactly two disjoint observed PAIRS whose unions are
            // complementary within C, so we can meet in the middle: enumerate the O(m^2) disjoint
            // pairs once, bucket them by their union bitset, then match each union U with C\U. This
            // visits only combinations that actually tile C (no rejected triples) -- O(m^2) build
            // plus one visit per valid boundary -- and returns exactly the same count as enumerate
            // for both FLAT (sum of pathcounts) and SHARED (admissible-boundary indicator).
            return countN2ViaPairs(cBits, subs, earlyExit);
        }
        BitSet used = BitSet.newBitSet(leafArraySize);
        List<BitSet> chosen = new ArrayList<>(m);
        return enumerate(cBits, subs, m, 0, used, chosen, earlyExit);
    }

    /**
     * Meet-in-the-middle count of boundary-4 (N_2) partitions of {@code cBits} into observed
     * subclades (see {@link #countNj}). A 4-part boundary {@code P1<P2<P3<P4} is the pairing
     * {P1,P2} | {P3,P4} of two disjoint observed pairs whose unions are complementary in C. We
     * enumerate the O(m^2) disjoint observed pairs, then for each pair {@code (i,j)} look up the
     * pairs {@code (k,l)} whose union is the complement {@code C∖union(i,j)} with {@code k > j}, so
     * each boundary is emitted exactly once (indices {@code i<j<k<l}).
     *
     * <p>Complements are found by <em>weighted-sum arithmetic</em>, not bitset materialisation.
     * {@link Clade#getCladeInBits()}'s {@code hashCode} weights word {@code w} by {@code w+1}, so
     * {@code weightedSum} is additive over disjoint sets: a pair's union sum is the sum of its two
     * parts' precomputed sums, and the complement's sum is {@code weightedSum(C) - union sum}. Pairs
     * are hashed by that {@code long} sum (no per-pair union bitset or {@code hashCode}); a chain hit
     * is confirmed <em>exactly</em> — cheaply and immune to sum wraparound — by the cardinality
     * identity {@code |Pi| + |Pj| + |Pk| + |Pl| == |C|} together with disjointness, which for parts
     * that are subclades of C is equivalent to tiling C. This is the dominant cost of building the
     * reserves and is allocation-free per pair. Shares the per-clade {@link #enumOps}/{@link
     * #opsBudget} guard with {@link #enumerate}, so a pathological clade still degrades gracefully to
     * the same capped-out {@code N_2 = 0}.
     */
    private int countN2ViaPairs(BitSet cBits, List<Clade> subs, boolean earlyExit) {
        int m = subs.size();
        BitSet[] sb = new BitSet[m];
        long[] partSum = new long[m];
        int[] partCard = new int[m];
        for (int i = 0; i < m; i++) {
            sb[i] = subs.get(i).getCladeInBits();
            partSum[i] = weightedSum(sb[i]);
            partCard[i] = sb[i].cardinality();
        }
        long sumC = weightedSum(cBits);
        int cardC = cBits.cardinality();
        final long[] ops = enumOps.get();

        // Enumerate disjoint observed pairs into flat arrays; entry e is pair (ei[e], ej[e]) with
        // union weighted-sum sU[e]. No per-pair bitset is built (sums are added from the parts).
        int cap = 256, P = 0;
        int[] ei = new int[cap], ej = new int[cap], eh = new int[cap];
        long[] sU = new long[cap];
        for (int i = 0; i < m; i++) {
            BitSet a = sb[i];
            for (int j = i + 1; j < m; j++) {
                if (++ops[0] > opsBudget) {
                    throw BUDGET_EXCEEDED;
                }
                if (a.intersects(sb[j])) {
                    continue;
                }
                if (P == cap) {
                    cap <<= 1;
                    ei = java.util.Arrays.copyOf(ei, cap);
                    ej = java.util.Arrays.copyOf(ej, cap);
                    eh = java.util.Arrays.copyOf(eh, cap);
                    sU = java.util.Arrays.copyOf(sU, cap);
                }
                long s = partSum[i] + partSum[j];
                ei[P] = i;
                ej[P] = j;
                sU[P] = s;
                eh[P] = mixHash(s);
                P++;
            }
        }
        if (P == 0) {
            return 0;
        }

        // Chained open hash over the union sums, so a pair can be found by its complement's sum.
        int tb = Integer.highestOneBit(P) << 1; // power of two in [P, 2P)
        int mask = tb - 1;
        int[] head = new int[tb];
        java.util.Arrays.fill(head, -1);
        int[] nxt = new int[P];
        for (int e = 0; e < P; e++) {
            int b = eh[e] & mask;
            nxt[e] = head[b];
            head[b] = e;
        }

        // For each pair (i,j) match pairs (k,l) with union sum == weightedSum(C∖union(i,j)) and k>j;
        // confirm exactly with the cardinality/disjointness tiling test.
        int count = 0;
        BitSet ue = BitSet.newBitSet(leafArraySize);
        BitSet uf = BitSet.newBitSet(leafArraySize);
        for (int e = 0; e < P; e++) {
            int i = ei[e], j = ej[e];
            long sR = sumC - sU[e];
            boolean ueBuilt = false;
            for (int f = head[mixHash(sR) & mask]; f >= 0; f = nxt[f]) {
                if (++ops[0] > opsBudget) {
                    throw BUDGET_EXCEEDED;
                }
                int k = ei[f], l = ej[f];
                if (k <= j) {
                    continue; // canonical i<j<k<l, so each boundary is counted once
                }
                if (partCard[i] + partCard[j] + partCard[k] + partCard[l] != cardC) {
                    continue; // parts cannot tile C (also filters sum-hash collisions)
                }
                if (!ueBuilt) {
                    ue.clear();
                    ue.or(sb[i]);
                    ue.or(sb[j]);
                    ueBuilt = true;
                }
                uf.clear();
                uf.or(sb[k]);
                uf.or(sb[l]);
                if (ue.intersects(uf)) {
                    continue; // disjoint + matching cardinality  <=>  the four parts tile C exactly
                }
                BitSet[] parts = {sb[i], sb[j], sb[k], sb[l]};
                int res = countAllNovelResolutions(cBits, parts);
                if (res > 0) {
                    count += (novelMode == NovelMode.FLAT) ? res : 1;
                    if (earlyExit) {
                        return count;
                    }
                }
            }
        }
        return count;
    }

    /**
     * Computes the boundary-3 ({@code N_1}) and boundary-4 ({@code N_2}) counts of {@code c}
     * together in a <em>single</em> disjoint-observed-pair pass (returns {@code {N_1, N_2}}). This is
     * the default reserve depth (k = 2) and by far the hottest part of building the reserves, so the
     * two orders share the one O(m^2) pair enumeration instead of scanning it twice.
     *
     * <p>Both boundaries are matched by weighted-sum arithmetic (see {@link #countN2ViaPairs}): for a
     * disjoint pair {@code (i,j)} the boundary-3 complement is the single observed subclade whose
     * {@code weightedSum} equals {@code weightedSum(C) - union sum} (looked up in a sum-keyed hash of
     * the subclades), and the boundary-4 complement is another disjoint pair with that sum. Every hit
     * is confirmed exactly with the cardinality/disjointness tiling test, so sum collisions and
     * wraparound are harmless. Honours the same {@link #enumOps}/{@link #opsBudget} guard.
     *
     * @param computeN2 when {@code false} (reserve depth k = 1, {@code eps} from {@code N_1} alone)
     *                  the boundary-4 stash and match are skipped and {@code N_2} is returned as 0.
     */
    private long[] countN1N2(Clade c, List<Clade> subs, boolean computeN2) {
        int m = subs.size();
        BitSet[] sb = new BitSet[m];
        long[] partSum = new long[m];
        int[] partCard = new int[m];
        for (int i = 0; i < m; i++) {
            sb[i] = subs.get(i).getCladeInBits();
            partSum[i] = weightedSum(sb[i]);
            partCard[i] = sb[i].cardinality();
        }
        BitSet cBits = c.getCladeInBits();
        long sumC = weightedSum(cBits);
        int cardC = c.size();
        final long[] ops = enumOps.get();

        // Sum-keyed chained hash of the subclades, for the boundary-3 complement (a single subclade).
        int stb = Integer.highestOneBit(Math.max(1, m)) << 1;
        int smask = stb - 1;
        int[] subHead = new int[stb];
        java.util.Arrays.fill(subHead, -1);
        int[] subNxt = new int[m];
        for (int i = 0; i < m; i++) {
            int b = mixHash(partSum[i]) & smask;
            subNxt[i] = subHead[b];
            subHead[b] = i;
        }

        // One disjoint-pair pass: score boundary-3 inline, stash pairs for the boundary-4 match.
        long n1 = 0;
        int cap = 256, P = 0;
        int[] ei = new int[cap], ej = new int[cap], eh = new int[cap];
        long[] sU = new long[cap];
        BitSet ue = BitSet.newBitSet(leafArraySize);
        BitSet uf = BitSet.newBitSet(leafArraySize);
        for (int i = 0; i < m; i++) {
            BitSet a = sb[i];
            for (int j = i + 1; j < m; j++) {
                if (++ops[0] > opsBudget) {
                    throw BUDGET_EXCEEDED;
                }
                if (a.intersects(sb[j])) {
                    continue;
                }
                long s = partSum[i] + partSum[j];
                long sR = sumC - s;
                // boundary-3: complement is one observed subclade d, canonical d > j
                boolean ueBuilt = false;
                for (int d = subHead[mixHash(sR) & smask]; d >= 0; d = subNxt[d]) {
                    if (d <= j || partSum[d] != sR) {
                        continue;
                    }
                    if (partCard[i] + partCard[j] + partCard[d] != cardC) {
                        continue;
                    }
                    if (!ueBuilt) {
                        ue.clear();
                        ue.or(a);
                        ue.or(sb[j]);
                        ueBuilt = true;
                    }
                    if (ue.intersects(sb[d])) {
                        continue;
                    }
                    BitSet[] parts = {sb[i], sb[j], sb[d]};
                    int res = countAllNovelResolutions(cBits, parts);
                    if (res > 0) {
                        n1 += (novelMode == NovelMode.FLAT) ? res : 1;
                    }
                }
                if (computeN2) {
                    if (P == cap) {
                        cap <<= 1;
                        ei = java.util.Arrays.copyOf(ei, cap);
                        ej = java.util.Arrays.copyOf(ej, cap);
                        eh = java.util.Arrays.copyOf(eh, cap);
                        sU = java.util.Arrays.copyOf(sU, cap);
                    }
                    ei[P] = i;
                    ej[P] = j;
                    sU[P] = s;
                    eh[P] = mixHash(s);
                    P++;
                }
            }
        }

        // boundary-4: match each pair with the pair whose union sum is the complement's, k > j.
        long n2 = 0;
        if (computeN2 && P > 0) {
            int tb = Integer.highestOneBit(P) << 1;
            int mask = tb - 1;
            int[] head = new int[tb];
            java.util.Arrays.fill(head, -1);
            int[] nxt = new int[P];
            for (int e = 0; e < P; e++) {
                int b = eh[e] & mask;
                nxt[e] = head[b];
                head[b] = e;
            }
            for (int e = 0; e < P; e++) {
                int i = ei[e], j = ej[e];
                long sR = sumC - sU[e];
                boolean ueBuilt = false;
                for (int f = head[mixHash(sR) & mask]; f >= 0; f = nxt[f]) {
                    if (++ops[0] > opsBudget) {
                        throw BUDGET_EXCEEDED;
                    }
                    int k = ei[f], l = ej[f];
                    if (k <= j) {
                        continue;
                    }
                    if (partCard[i] + partCard[j] + partCard[k] + partCard[l] != cardC) {
                        continue;
                    }
                    if (!ueBuilt) {
                        ue.clear();
                        ue.or(sb[i]);
                        ue.or(sb[j]);
                        ueBuilt = true;
                    }
                    uf.clear();
                    uf.or(sb[k]);
                    uf.or(sb[l]);
                    if (ue.intersects(uf)) {
                        continue;
                    }
                    BitSet[] parts = {sb[i], sb[j], sb[k], sb[l]};
                    int res = countAllNovelResolutions(cBits, parts);
                    if (res > 0) {
                        n2 += (novelMode == NovelMode.FLAT) ? res : 1;
                    }
                }
            }
        }
        return new long[]{n1, n2};
    }

    /** Weighted bit-sum matching {@link Clade#getCladeInBits()}'s {@code hashCode} word weighting
     *  (word {@code w} weighted by {@code w+1}); additive over disjoint sets, so it lets the
     *  boundary matcher find a pair's complement by {@code long} arithmetic. Wraparound mod 2^64 is
     *  harmless — the matcher confirms hits with an exact cardinality/disjointness test. */
    private static long weightedSum(BitSet b) {
        long s = 0;
        for (int bit = b.nextSetBit(0); bit >= 0; bit = b.nextSetBit(bit + 1)) {
            s += ((long) ((bit >>> 6) + 1)) << (bit & 63);
        }
        return s;
    }

    /** Scalar mixer spreading a {@code long} union sum to a hash-table bucket. */
    private static int mixHash(long x) {
        long z = x * 0x9E3779B97F4A7C15L;
        return (int) (z ^ (z >>> 32));
    }

    /* Recursive canonical enumeration: pick (m-1) parts at strictly increasing
     * indices (= increasing canonical order), then derive the last part as the
     * complement and require it to be the canonical-largest. */
    private int enumerate(BitSet cBits, List<Clade> subs, int m, int startIdx,
                          BitSet used, List<BitSet> chosen, boolean earlyExit) {
        if (++enumOps.get()[0] > opsBudget) {
            throw BUDGET_EXCEEDED;
        }
        if (chosen.size() == m - 1) {
            BitSet last = (BitSet) cBits.clone();
            last.andNot(used);
            if (last.cardinality() == 0 || getClade(last) == null) {
                return 0;
            }
            BitSet prev = chosen.get(chosen.size() - 1);
            // canonical: last must be strictly greater than the previous (largest) chosen
            if (compareBitSets(prev, last) >= 0) {
                return 0;
            }
            BitSet[] parts = new BitSet[m];
            for (int i = 0; i < m - 1; i++) {
                parts[i] = chosen.get(i);
            }
            parts[m - 1] = last;
            // SHARED counts admissible boundaries (N_j, an indicator); FLAT counts distinct novel
            // trees (M_j = sum of pathcount over boundaries) by returning the pathcount itself.
            int res = countAllNovelResolutions(cBits, parts);
            return novelMode == NovelMode.FLAT ? res : (res > 0 ? 1 : 0);
        }
        int count = 0;
        for (int i = startIdx; i < subs.size(); i++) {
            BitSet pb = subs.get(i).getCladeInBits();
            if (pb.intersects(used)) {
                continue;
            }
            chosen.add(pb);
            BitSet newUsed = (BitSet) used.clone();
            newUsed.or(pb);
            count += enumerate(cBits, subs, m, i + 1, newUsed, chosen, earlyExit);
            chosen.remove(chosen.size() - 1);
            if (earlyExit && count > 0) {
                return count;
            }
        }
        return count;
    }

    /* Monotone bisection solve of sum_{m>=3} n[m] eps^(m-2) = mu on eps >= 0
     * (eps^(m-2) = eps^(#novel clades)). */
    private static double solveLogEps(int[] n, int maxDeg, double mu) {
        boolean any = false;
        for (int j = 3; j <= maxDeg; j++) {
            if (n[j] > 0) {
                any = true;
                break;
            }
        }
        if (!any) {
            return Double.NEGATIVE_INFINITY; // not reservable; eps unused
        }
        double lo = 0.0;
        double hi = 1.0;
        while (evalPoly(n, maxDeg, hi) < mu) {
            hi *= 2.0;
        }
        for (int it = 0; it < 100; it++) {
            double mid = 0.5 * (lo + hi);
            if (evalPoly(n, maxDeg, mid) < mu) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        return Math.log(0.5 * (lo + hi));
    }

    /* Consistent total order on clade bitsets: compare sorted set-bit indices
     * lexicographically (a proper prefix sorts before its extension). */
    private static int compareBitSets(BitSet a, BitSet b) {
        int ia = a.nextSetBit(0);
        int ib = b.nextSetBit(0);
        while (ia >= 0 && ib >= 0) {
            if (ia != ib) {
                return Integer.compare(ia, ib);
            }
            ia = a.nextSetBit(ia + 1);
            ib = b.nextSetBit(ib + 1);
        }
        // whichever still has set bits is "greater"; -1 (exhausted) sorts first
        return Integer.compare(ia, ib);
    }

    /* sum_{m>=3} n[m] * x^(m - 2). */
    private static double evalPoly(int[] n, int maxDeg, double x) {
        double s = 0.0;
        for (int m = 3; m <= maxDeg; m++) {
            if (n[m] != 0) {
                s += n[m] * Math.pow(x, m - EXP_REDUCTION);
            }
        }
        return s;
    }
}
