package ccd.model;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator.TreeSet;
import ccd.model.bitsets.BitSet;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * CRegCCD -- the <em>class-based</em> regularised CCD (Jonathan's proposal): a full-support tree
 * distribution obtained by additive smoothing over <em>all</em> bipartitions of every clade, with the
 * pseudocount depending only on which of four classes a bipartition falls into.
 *
 * <p>At a clade {@code C} of {@code m} taxa, each of the {@code 2^(m-1) - 1} bipartitions
 * {@code {A, B}} belongs to exactly one class:
 * <ol>
 *   <li>{@code A_1}: the split was observed in the sample;</li>
 *   <li>{@code A_2}: unobserved split, both {@code A} and {@code B} are observed clades
 *       (this is exactly the CCD0 split expansion);</li>
 *   <li>{@code A_3}: unobserved split, exactly one of {@code A}, {@code B} is an observed clade;</li>
 *   <li>{@code A_4}: neither {@code A} nor {@code B} is an observed clade.</li>
 * </ol>
 * Only classes 1--3 are ever materialised; {@code |A_4|} follows in closed form as
 * {@code (2^(m-1) - 1) - |A_1| - |A_2| - |A_3|}, so the whole of tree space is represented without
 * enumerating it. There is no escape probability, no reserve equation and no region decomposition.
 *
 * <p><b>Per-class totals, not per-split constants.</b> Each {@code a_j} is the <em>total</em>
 * pseudocount mass of its class, so the per-split pseudocount is {@code a_j / |A_j|} and
 * <pre>
 *   theta(S) = (f(S) + a_{j(S)} / |A_{j(S)}|) / (f(C) + sum over non-empty classes of a_j).
 * </pre>
 * This matters: {@code |A_4|} is essentially {@code 2^(m-1)}, so a constant per-split pseudocount
 * would give class 4 all the mass and the data none (on 40 taxa the observed splits retain about
 * {@code 6e-9} of the probability at a root clade; see {@code SplitClassSizeAnalysis}). With per-class
 * totals the retained mass is independent of taxon count. Equivalently: draw a class, then draw
 * uniformly within it.
 *
 * <p>Consequences, all by construction rather than by correction:
 * <ul>
 *   <li><b>Exactly normalised.</b> Each clade -- observed or not -- carries a normalised categorical
 *       over all of its bipartitions, and a tree's probability is the chain-rule product over its
 *       internal nodes, so the distribution sums to one over tree space by induction. No partition
 *       function is needed, and the {@code Theta(mu^2)} maximality deficit of {@link KRegCCD}
 *       cannot arise because there are no regions.</li>
 *   <li><b>Full support.</b> Every bipartition of every clade has positive probability whenever
 *       {@code a_2, a_3, a_4 > 0}.</li>
 *   <li><b>Consistent.</b> The denominator contains {@code f(C)}, so the mass held back for unobserved
 *       splits shrinks as the evidence for {@code C} accumulates.</li>
 *   <li><b>No mass reserved where nothing can escape.</b> A clade all of whose bipartitions were
 *       observed has classes 2--4 empty, so those {@code a_j} never enter the denominator.</li>
 * </ul>
 *
 * @author Claude
 */
public class CRegCCD extends CCD1 {

    /**
     * Per-split pseudocount on the CCD0 split set (splits introducing no novel clade). This is
     * regCCD's {@code alpha}; the fitted value was 0.4 on every real data set tested.
     */
    public static final double DEFAULT_ALPHA = 0.4;
    /** Total prior mass for splits introducing one novel clade. */
    public static final double DEFAULT_ALPHA1 = 0.4;
    /** Total prior mass for splits introducing two novel clades. */
    public static final double DEFAULT_ALPHA2 = 0.05;

    private final double alpha;
    private final double alpha1;
    private final double alpha2;

    /** Class sizes are parameter-independent, so they are computed once and reused across a search.
     *  Concurrent because {@link #sampleTrees} fans draws out over threads. */
    private final Map<BitSet, double[]> sizeCache = new java.util.concurrent.ConcurrentHashMap<>();
    private volatile List<BitSet> sortedCladeBits;

    /**
     * Strictly-positive height increment for a novel internal node whose clade has no recorded
     * height, so that branch lengths stay positive.
     */
    private static final double NOVEL_HEIGHT_EPS = 1e-8;

    /**
     * Draw a class-4 split by rejection while at least this fraction of bipartitions are class 4,
     * which bounds the expected number of attempts by its reciprocal; below it, enumerate instead.
     * Since classes 1-3 are only polynomially large, a small acceptance rate implies a small
     * {@code 2^(m-1)}, so the enumeration branch is always cheap.
     */
    private static final double MIN_REJECTION_ACCEPTANCE = 0.02;

    public CRegCCD(List<Tree> trees, double burnin) {
        this(trees, burnin, DEFAULT_ALPHA, DEFAULT_ALPHA1, DEFAULT_ALPHA2);
    }

    public CRegCCD(List<Tree> trees, double burnin, double alpha, double alpha1, double alpha2) {
        super(trees, burnin);
        validate(alpha, alpha1, alpha2);
        this.alpha = alpha;
        this.alpha1 = alpha1;
        this.alpha2 = alpha2;
    }

    public CRegCCD(TreeSet treeSet) {
        this(treeSet, DEFAULT_ALPHA, DEFAULT_ALPHA1, DEFAULT_ALPHA2);
    }

    public CRegCCD(TreeSet treeSet, double alpha, double alpha1, double alpha2) {
        super(treeSet, false);
        validate(alpha, alpha1, alpha2);
        this.alpha = alpha;
        this.alpha1 = alpha1;
        this.alpha2 = alpha2;
    }

    private static void validate(double alpha, double alpha1, double alpha2) {
        if (alpha <= 0) {
            throw new IllegalArgumentException("alpha must be > 0, got " + alpha);
        }
        if (alpha1 < 0 || alpha2 < 0) {
            throw new IllegalArgumentException(
                    "alpha1 and alpha2 must be >= 0, got " + alpha1 + ", " + alpha2);
        }
    }

    /** regCCD's per-split pseudocount on the CCD0 split set. */
    public double getAlpha() {
        return alpha;
    }

    /** Total prior mass for splits introducing one novel clade. */
    public double getAlpha1() {
        return alpha1;
    }

    /** Total prior mass for splits introducing two novel clades. */
    public double getAlpha2() {
        return alpha2;
    }

    @Override
    public String toString() {
        return "CRegCCD(alpha=" + alpha + ", alpha1=" + alpha1 + ", alpha2=" + alpha2 + ")";
    }

    /**
     * Per-split pseudocount of each split class at {@code cBits} (as logs, so that
     * {@code alpha2/|A_2|} with an exponentially large class cannot underflow), plus the normaliser
     * {@code Z} in the last slot. Indices 0 and 1 are the two halves of the CCD0 split set and share
     * the per-split {@code alpha}; indices 2 and 3 are the one- and two-novel-clade classes, whose
     * class totals are spread over their members.
     */
    private double[] logWeightsAndZ(BitSet cBits, double alpha, double alpha1, double alpha2) {
        double[] size = classSizes(cBits);
        Clade c = getClade(cBits);
        double z = (c != null) ? c.getNumberOfOccurrences() : 0.0;
        double[] w = {Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY,
                Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY};
        double n0 = size[0] + size[1];
        if (n0 > 0) {
            w[0] = Math.log(alpha);
            w[1] = w[0];
            z += alpha * n0;
        }
        if (size[2] > 0) {
            w[2] = Math.log(alpha1) - Math.log(size[2]);
            z += alpha1;
        }
        if (size[3] > 0) {
            w[3] = Math.log(alpha2) - Math.log(size[3]);
            z += alpha2;
        }
        return new double[]{w[0], w[1], w[2], w[3], z};
    }

    /* ----------------------------------------------------------------------
     * Scoring
     * ------------------------------------------------------------------- */

    @Override
    public double getLogProbabilityOfTree(Tree tree) {
        return scoreTree(tree, alpha, alpha1, alpha2);
    }

    /**
     * Log probability at pseudocounts other than this model's own, reusing the cached (parameter-free)
     * class sizes. Lets a cross-validation sweep evaluate many parameter vectors on one trained model.
     */
    public double getLogProbabilityOfTree(Tree tree, double b0, double b1, double b2) {
        validate(b0, b1, b2);
        return scoreTree(tree, b0, b1, b2);
    }

    @Override
    public double getProbabilityOfTree(Tree tree) {
        return Math.exp(getLogProbabilityOfTree(tree));
    }

    /** Always true: CRegCCD has full support over the trees on its taxon set. */
    @Override
    public boolean containsTree(Tree tree) {
        return true;
    }

    private double scoreTree(Tree tree, double b0, double b1, double b2) {
        Map<Node, BitSet> bits = new HashMap<>();
        computeBits(tree.getRoot(), bits);
        double logp = 0.0;
        for (Node v : tree.getNodesAsArray()) {
            if (v.isLeaf()) {
                continue;
            }
            logp += logSplitProbability(bits.get(v),
                    bits.get(v.getChildren().get(0)),
                    bits.get(v.getChildren().get(1)),
                    b0, b1, b2);
        }
        return logp;
    }

    /**
     * Log conditional probability of the bipartition {@code {aBits, bBits}} of clade {@code cBits},
     * for a clade that need not be observed. This is the whole model: every internal node of a tree
     * contributes exactly one such factor.
     */
    double logSplitProbability(BitSet cBits, BitSet aBits, BitSet bBits,
                               double b0, double b1, double b2) {
        double[] wz = logWeightsAndZ(cBits, b0, b1, b2);
        int cls = splitClass(cBits, aBits, bBits);
        double fS = (cls == 0)
                ? observedPartition(getClade(cBits), aBits, bBits).getNumberOfOccurrences() : 0.0;
        double logNumerator = (fS > 0) ? Math.log(fS + Math.exp(wz[cls])) : wz[cls];
        return logNumerator - Math.log(wz[4]);
    }

    /**
     * Which of the four classes the bipartition {@code {aBits, bBits}} of {@code cBits} belongs to,
     * as a 0-based index (0 = observed split, 3 = neither child observed).
     */
    int splitClass(BitSet cBits, BitSet aBits, BitSet bBits) {
        if (observedPartition(getClade(cBits), aBits, bBits) != null) {
            return 0;
        }
        boolean aObs = getClade(aBits) != null;
        boolean bObs = getClade(bBits) != null;
        return (aObs && bObs) ? 1 : ((aObs || bObs) ? 2 : 3);
    }

    private CladePartition observedPartition(Clade c, BitSet aBits, BitSet bBits) {
        if (c == null) {
            return null;
        }
        Clade ca = getClade(aBits);
        Clade cb = getClade(bBits);
        if (ca == null || cb == null) {
            return null;
        }
        return c.getCladePartition(ca, cb);
    }

    /* ----------------------------------------------------------------------
     * Class sizes
     * ------------------------------------------------------------------- */

    /**
     * Sizes {@code {|A_1|, |A_2|, |A_3|, |A_4|}} of the four split classes of clade {@code cBits}.
     *
     * <p>{@code |A_1|} is the number of observed splits; a single pass over the observed subclades of
     * {@code C} yields {@code |A_3|} (observed subclade whose complement is not observed) and the
     * number of observed-clade pairs, from which {@code |A_2|} follows; {@code |A_4|} is the
     * remainder of {@code 2^(m-1) - 1}. Cached, and independent of the pseudocounts.
     */
    double[] classSizes(BitSet cBits) {
        double[] cached = sizeCache.get(cBits);
        if (cached != null) {
            return cached;
        }
        int m = cBits.cardinality();
        Clade c = getClade(cBits);

        double n1 = 0.0;
        if (c != null) {
            for (CladePartition p : c.getPartitions()) {
                if (p.getNumberOfOccurrences() > 0) {
                    n1++;
                }
            }
        }

        // one pass over observed clades strictly inside C
        int bothObservedEnds = 0; // counts each both-observed bipartition twice (once per side)
        double n3 = 0.0;
        for (BitSet d : sortedCladeBits()) {
            if (d.cardinality() >= m || !subset(d, cBits)) {
                continue;
            }
            BitSet complement = BitSet.newBitSet(cBits);
            complement.andNot(d);
            if (getClade(complement) != null) {
                bothObservedEnds++;
            } else {
                n3++;
            }
        }
        double n2 = Math.max(0.0, bothObservedEnds / 2.0 - n1);

        double total = Math.pow(2.0, m - 1) - 1.0;
        double n4 = Math.max(0.0, total - n1 - n2 - n3);

        double[] size = {n1, n2, n3, n4};
        sizeCache.put(BitSet.newBitSet(cBits), size);
        return size;
    }

    private synchronized List<BitSet> sortedCladeBits() {
        if (sortedCladeBits == null) {
            List<BitSet> all = new ArrayList<>();
            for (Clade c : getClades()) {
                all.add(c.getCladeInBits());
            }
            all.sort(CRegCCD::compareBitSets);
            sortedCladeBits = all;
        }
        return sortedCladeBits;
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

    private static boolean subset(BitSet a, BitSet c) {
        BitSet tmp = BitSet.newBitSet(a);
        tmp.andNot(c);
        return tmp.isEmpty();
    }

    /* ----------------------------------------------------------------------
     * MAP tree
     *
     * The maximum over all of tree space is a DP over the subset lattice, so instead we run the DP
     * over the observed-clade DAG using only the both-children-observed splits (classes 1 and 2 --
     * exactly CCD0's split set) and then *certify* that no off-backbone tree can beat it.
     *
     * The certificate is a pair of upper-bound DPs over the same DAG. U(C) bounds the best subtree
     * log-probability over ALL trees on C, and V(C) bounds it over trees that use at least one
     * off-backbone (class 3 or 4) split. Any subtree contributes at most 0, so a novel child is
     * bounded by 0; every class-3 split at C shares one theta, as does every class-4 split, because
     * the model is uniform within a class. If best(root) > V(root), no tree using an off-backbone
     * split anywhere can beat the backbone optimum, so the backbone MAP is the global MAP.
     * ------------------------------------------------------------------- */

    private volatile Map<BitSet, Double> mapBest;
    private volatile Map<BitSet, BitSet[]> mapArg;
    private volatile double offBackboneBound = Double.NaN;

    /** All both-children-observed bipartitions of {@code cBits} (classes 1 and 2), each once. */
    private List<BitSet[]> backboneSplits(BitSet cBits) {
        int m = cBits.cardinality();
        List<BitSet[]> out = new ArrayList<>();
        for (BitSet d : sortedCladeBits()) {
            if (d.cardinality() >= m || !subset(d, cBits)) {
                continue;
            }
            BitSet complement = BitSet.newBitSet(cBits);
            complement.andNot(d);
            if (getClade(complement) == null || compareBitSets(d, complement) >= 0) {
                continue;
            }
            out.add(new BitSet[]{d, complement});
        }
        return out;
    }

    /** Log theta shared by every split of the given off-backbone class at {@code cBits}. */
    private double logThetaOfClass(BitSet cBits, int cls) {
        double[] size = classSizes(cBits);
        if (size[cls] <= 0) {
            return Double.NEGATIVE_INFINITY;
        }
        double[] wz = logWeightsAndZ(cBits, alpha, alpha1, alpha2);
        return wz[cls] - Math.log(wz[4]);
    }

    private synchronized void computeMAP() {
        if (mapBest != null) {
            return;
        }
        List<Clade> clades = new ArrayList<>(getClades());
        clades.sort(java.util.Comparator.comparingInt(Clade::size));

        Map<BitSet, Double> best = new HashMap<>();
        Map<BitSet, BitSet[]> arg = new HashMap<>();
        Map<BitSet, Double> upper = new HashMap<>();      // U: best over all trees
        Map<BitSet, Double> upperOff = new HashMap<>();   // V: best over trees using an off-backbone split

        for (Clade c : clades) {
            BitSet cb = c.getCladeInBits();
            if (c.size() == 1) {
                best.put(cb, 0.0);
                upper.put(cb, 0.0);
                upperOff.put(cb, Double.NEGATIVE_INFINITY);
                continue;
            }
            double bBest = Double.NEGATIVE_INFINITY;
            BitSet[] bArg = null;
            double bUpper = Double.NEGATIVE_INFINITY;
            double bOff = Double.NEGATIVE_INFINITY;

            for (BitSet[] s : backboneSplits(cb)) {
                Double l = best.get(s[0]);
                Double r = best.get(s[1]);
                if (l == null || r == null) {
                    continue;
                }
                double theta = logSplitProbability(cb, s[0], s[1], alpha, alpha1, alpha2);
                double v = theta + l + r;
                if (v > bBest) {
                    bBest = v;
                    bArg = s;
                }
                double ul = upper.get(s[0]);
                double ur = upper.get(s[1]);
                bUpper = Math.max(bUpper, theta + ul + ur);
                double vl = upperOff.get(s[0]);
                double vr = upperOff.get(s[1]);
                bOff = Math.max(bOff, theta + Math.max(vl + ur, ul + vr));
            }

            // class 3: one child observed (bounded above by U of that child, novel side by 0)
            double t3 = logThetaOfClass(cb, 2);
            if (t3 > Double.NEGATIVE_INFINITY) {
                double bestObservedSide = Double.NEGATIVE_INFINITY;
                int m = cb.cardinality();
                for (BitSet d : sortedCladeBits()) {
                    if (d.cardinality() >= m || !subset(d, cb)) {
                        continue;
                    }
                    BitSet complement = BitSet.newBitSet(cb);
                    complement.andNot(d);
                    if (getClade(complement) == null) { // exactly one side observed
                        Double u = upper.get(d);
                        if (u != null) {
                            bestObservedSide = Math.max(bestObservedSide, u);
                        }
                    }
                }
                if (bestObservedSide > Double.NEGATIVE_INFINITY) {
                    bUpper = Math.max(bUpper, t3 + bestObservedSide);
                    bOff = Math.max(bOff, t3 + bestObservedSide);
                }
            }

            // class 4: both children novel, each bounded by 0
            double t4 = logThetaOfClass(cb, 3);
            if (t4 > Double.NEGATIVE_INFINITY) {
                bUpper = Math.max(bUpper, t4);
                bOff = Math.max(bOff, t4);
            }

            best.put(cb, bBest);
            arg.put(cb, bArg);
            upper.put(cb, bUpper);
            upperOff.put(cb, bOff);
        }

        this.offBackboneBound = upperOff.get(getRootClade().getCladeInBits());
        this.mapArg = arg;
        this.mapBest = best;
    }

    /**
     * Exact MAP over all trees that use no two-novel-clade split, by memoised recursion over
     * classes {@code A_0} and {@code A_1}.
     *
     * <p>Restricting to {@code A_0} keeps the recursion on the observed-clade DAG. Admitting
     * {@code A_1} as well -- peel off an observed clade, leave a novel remainder -- widens the state
     * space to clades of the form {@code root} minus a union of disjoint observed clades. That set
     * can in principle be large, so the search is capped by {@link #MAP_STATE_BUDGET} distinct
     * clades; in practice it stays small because an {@code A_1} split is expensive and the recursion
     * only ever descends.
     *
     * <p>Returns {@code {best, viaA2}} for the clade: the best log probability using only
     * {@code A_0}/{@code A_1} splits, and an upper bound on any subtree that uses an {@code A_2}
     * split somewhere. The second is the certificate: if {@code best > viaA2} at the root, no tree
     * containing a two-novel-clade split can reach the optimum, so the answer is the global MAP.
     */
    private static final long MAP_STATE_BUDGET =
            Long.getLong("creg.mapStates", 4_000_000L);

    private static final class BudgetExhausted extends RuntimeException {
        BudgetExhausted() {
            super(null, null, false, false);
        }
    }

    private double[] solveFull(BitSet cBits, int a1Budget, List<Map<BitSet, double[]>> memos,
                               long[] states) {
        Map<BitSet, double[]> memo = memos.get(a1Budget);
        double[] cached = memo.get(cBits);
        if (cached != null) {
            return cached;
        }
        if (++states[0] > MAP_STATE_BUDGET) {
            throw new BudgetExhausted();
        }
        int m = cBits.cardinality();
        if (m == 1) {
            double[] leaf = {0.0, Double.NEGATIVE_INFINITY};
            memo.put(BitSet.newBitSet(cBits), leaf);
            return leaf;
        }
        double best = Double.NEGATIVE_INFINITY;
        double viaA2 = Double.NEGATIVE_INFINITY;

        for (BitSet d : sortedCladeBits()) {
            if (d.cardinality() >= m || !subset(d, cBits)) {
                continue;
            }
            BitSet complement = BitSet.newBitSet(cBits);
            complement.andNot(d);
            boolean complementObserved = getClade(complement) != null;
            if (complementObserved && compareBitSets(d, complement) >= 0) {
                continue; // both-observed bipartitions are reached from their smaller side only
            }
            int childBudget = complementObserved ? a1Budget : a1Budget - 1;
            if (childBudget < 0) {
                continue;   // no A_1 split allowance left on this path
            }
            double theta = logSplitProbability(cBits, d, complement, alpha, alpha1, alpha2);
            double[] left = solveFull(d, childBudget, memos, states);
            double[] right = solveFull(complement, childBudget, memos, states);
            best = Math.max(best, theta + left[0] + right[0]);
            double ul = Math.max(left[0], left[1]);
            double ur = Math.max(right[0], right[1]);
            viaA2 = Math.max(viaA2, theta + Math.max(left[1] + ur, ul + right[1]));
        }

        // an A_2 split taken here; both children are novel and bounded above by zero
        double[] size = classSizes(cBits);
        if (size[3] > 0) {
            viaA2 = Math.max(viaA2, logThetaOfClass(cBits, 3));
        }

        double[] value = {best, viaA2};
        memo.put(BitSet.newBitSet(cBits), value);
        return value;
    }

    /**
     * Result of the MAP search.
     *
     * <p>{@code a2Excluded} says only that no {@code A_2} split can improve the optimum <em>within
     * the searched class of trees</em>, i.e. among trees using at most {@code a1Depth} one-novel-clade
     * splits per path. It upgrades to a genuine global certificate ({@code certifiedGlobal}) only
     * when that depth was not binding, so that every {@code A_0}/{@code A_1} tree was considered.
     */
    public record MapResult(double maxLogProbability, double offBackboneBound,
                            boolean a2Excluded, int a1Depth, boolean exhaustiveInA1,
                            long statesExplored, boolean complete) {

        /** True only when the search covered every A_0/A_1 tree and excluded A_2 as well. */
        public boolean certifiedGlobal() {
            return complete && exhaustiveInA1 && a2Excluded;
        }
    }

    /**
     * Runs the {@code A_0}/{@code A_1} search and reports whether the optimum it found is provably
     * the global MAP. {@code complete} is false when the state budget was exhausted, in which case
     * the backbone DP result should be used instead.
     */
    public MapResult solveMAP() {
        return solveMAP(Integer.MAX_VALUE / 2);
    }

    /**
     * As {@link #solveMAP()} but allowing at most {@code maxA1} one-novel-clade splits on any
     * root-to-leaf path. {@code maxA1 = 0} is the backbone DP; raising it enlarges the search until
     * either the optimum stops improving or the state budget is exhausted.
     */
    public MapResult solveMAP(int maxA1) {
        int cap = Math.min(maxA1, getSizeOfLeavesArray());
        List<Map<BitSet, double[]>> memos = new ArrayList<>();
        for (int i = 0; i <= cap; i++) {
            memos.add(new HashMap<>());
        }
        long[] states = {0};
        boolean exhaustive = cap >= getSizeOfLeavesArray() - 2;
        try {
            double[] root = solveFull(getRootClade().getCladeInBits(), cap, memos, states);
            return new MapResult(root[0], root[1], root[0] > root[1], cap, exhaustive,
                    states[0], true);
        } catch (BudgetExhausted e) {
            return new MapResult(getMaxLogTreeProbability(), getOffBackboneBound(),
                    false, cap, exhaustive, states[0], false);
        }
    }

    /** Log probability of the backbone MAP tree. */
    @Override
    public double getMaxLogTreeProbability() {
        computeMAP();
        return mapBest.get(getRootClade().getCladeInBits());
    }

    /**
     * Whether the backbone MAP tree is provably the global MAP over all of tree space: true when no
     * tree using a class-3 or class-4 split anywhere can reach the backbone optimum.
     */
    public boolean isMAPCertifiedGlobal() {
        computeMAP();
        return getMaxLogTreeProbability() > offBackboneBound;
    }

    /** The certificate's upper bound on any tree that uses an off-backbone split. */
    public double getOffBackboneBound() {
        computeMAP();
        return offBackboneBound;
    }

    /* ----------------------------------------------------------------------
     * Entropy
     *
     * The sampler draws from exactly the scored distribution, so E[-log q] is an unbiased estimate
     * of H(q) with no truncation to correct for. A deterministic recursion would instead have to
     * approximate the subtree entropy of novel clades, so the Monte-Carlo estimator is both simpler
     * and more accurate here; only its standard error stands between it and the exact value.
     * ------------------------------------------------------------------- */

    /* ----------------------------------------------------------------------
     * Deterministic entropy recursion
     *
     * H(C) = H_split(C) + sum_S theta(S) [H(A_S) + H(B_S)], with H(leaf) = 0.
     *
     * The local term is closed form: within class j >= 2 every member has the same
     * theta_j = a_j / (|A_j| Z), so that class contributes -(a_j/Z) log theta_j as a single term --
     * the exponentially large class 4 is never enumerated. The expectation term is exact for
     * classes 1 and 2 (both children observed, so the recursion stays on the clade DAG) and needs a
     * value for the novel child of a class-3 split and for both children of a class-4 split.
     *
     * APPROXIMATION: a novel clade is treated as *fresh*, i.e. as containing no observed clades
     * other than its singletons, so its subtree entropy depends only on its size and is given by a
     * universal g(k) computed once by an O(n^2) recursion. Real novel clades usually do contain
     * observed clades, so g overestimates their structure-free entropy; the error enters only
     * through the class-3 and class-4 branches, whose total weight at a clade is (a_3 + a_4)/Z.
     * Class-4 splits are grouped by the sizes of the two sides, whose counts follow from the
     * binomials minus the enumerable classes, keeping the whole pass O(K + m) per clade.
     * ------------------------------------------------------------------- */

    /** g(k): subtree entropy of a fresh (no observed subclades but singletons) clade of size k. */
    private volatile double[] freshEntropy;

    private synchronized double[] freshEntropy() {
        if (freshEntropy != null) {
            return freshEntropy;
        }
        int n = getSizeOfLeavesArray();
        double[] g = new double[Math.max(3, n + 1)];
        g[1] = 0.0;
        if (g.length > 2) {
            g[2] = 0.0; // the single split of a novel cherry has probability one
        }
        for (int k = 3; k <= n; k++) {
            double total = Math.pow(2.0, k - 1) - 1.0;
            double n3 = k;                 // {leaf, rest}, rest unobserved since k-1 >= 2
            double n4 = total - n3;        // both sides of size >= 2, so both unobserved
            double z = alpha1 + (n4 > 0 ? alpha2 : 0.0);

            double logT3 = Math.log(alpha1) - Math.log(n3) - Math.log(z);
            double h = -(alpha1 / z) * logT3;
            double e = (alpha1 / z) * g[k - 1]; // class-3 children are {1, k-1}
            if (n4 > 0) {
                double logT4 = Math.log(alpha2) - Math.log(n4) - Math.log(z);
                h -= (alpha2 / z) * logT4;
                double weighted = 0.0;
                for (int j = 2; j <= k / 2; j++) {
                    double cnt = binomial(k, j);
                    if (j == k - j) {
                        cnt /= 2.0;
                    }
                    weighted += cnt * (g[j] + g[k - j]);
                }
                e += Math.exp(logT4) * weighted;
            }
            g[k] = h + e;
        }
        freshEntropy = g;
        return g;
    }

    /** Local split entropy at {@code cBits}: {@code -sum_S theta(S) log theta(S)}, closed form. */
    private double localSplitEntropy(BitSet cBits) {
        double[] size = classSizes(cBits);
        double[] wz = logWeightsAndZ(cBits, alpha, alpha1, alpha2);
        double z = wz[4];
        Clade c = getClade(cBits);
        double h = 0.0;
        if (size[0] > 0) { // class 1 is explicit: theta varies with the split count
            double perSplit = Math.exp(wz[0]);
            for (CladePartition p : c.getPartitions()) {
                if (p.getNumberOfOccurrences() == 0) {
                    continue;
                }
                double theta = (p.getNumberOfOccurrences() + perSplit) / z;
                h -= theta * Math.log(theta);
            }
        }
        for (int j = 1; j < 4; j++) { // classes 2-4 are uniform within the class
            if (size[j] > 0) {
                double logTheta = wz[j] - Math.log(z);
                h -= size[j] * Math.exp(logTheta) * logTheta;
            }
        }
        return h;
    }

    /**
     * Deterministic entropy in nats, using the fresh-clade approximation for novel subclades. Exact
     * whenever no class-3 or class-4 split leads to a novel clade that contains an observed clade.
     */
    public double getEntropyRecursive() {
        double[] g = freshEntropy();
        List<Clade> clades = new ArrayList<>(getClades());
        clades.sort(java.util.Comparator.comparingInt(Clade::size));
        Map<BitSet, Double> entropy = new HashMap<>();

        for (Clade c : clades) {
            BitSet cb = c.getCladeInBits();
            int m = c.size();
            if (m == 1) {
                entropy.put(cb, 0.0);
                continue;
            }
            double[] size = classSizes(cb);
            double[] wz = logWeightsAndZ(cb, alpha, alpha1, alpha2);
            double z = wz[4];

            double e = 0.0;

            // class 1: observed splits, both children observed
            if (size[0] > 0) {
                double perSplit = Math.exp(wz[0]);
                for (CladePartition p : c.getPartitions()) {
                    if (p.getNumberOfOccurrences() == 0) {
                        continue;
                    }
                    double theta = (p.getNumberOfOccurrences() + perSplit) / z;
                    e += theta * (entropy.get(p.getChildClades()[0].getCladeInBits())
                            + entropy.get(p.getChildClades()[1].getCladeInBits()));
                }
            }

            // classes 2 and 3, plus the size profile of everything that is not class 4
            double t2 = (size[1] > 0) ? Math.exp(wz[1]) / z : 0.0;
            double t3 = (size[2] > 0) ? Math.exp(wz[2]) / z : 0.0;
            double[] nonClass4 = new double[m / 2 + 1];
            for (BitSet d : sortedCladeBits()) {
                if (d.cardinality() >= m || !subset(d, cb)) {
                    continue;
                }
                BitSet complement = BitSet.newBitSet(cb);
                complement.andNot(d);
                int j = Math.min(d.cardinality(), complement.cardinality());
                if (getClade(complement) != null) {
                    if (compareBitSets(d, complement) < 0) {
                        nonClass4[j]++;
                        if (observedPartition(c, d, complement) == null) { // class 2
                            e += t2 * (entropy.get(d) + entropy.get(complement));
                        }
                    }
                } else { // class 3: d observed, complement novel
                    nonClass4[j]++;
                    e += t3 * (entropy.get(d) + g[complement.cardinality()]);
                }
            }

            // class 4, grouped by the sizes of the two sides
            if (size[3] > 0) {
                double logT4 = wz[3] - Math.log(z);
                double weighted = 0.0;
                for (int j = 1; j <= m / 2; j++) {
                    double totalPairs = binomial(m, j);
                    if (j == m - j) {
                        totalPairs /= 2.0;
                    }
                    double count4 = totalPairs - nonClass4[j];
                    if (count4 > 0) {
                        weighted += count4 * (g[j] + g[m - j]);
                    }
                }
                e += Math.exp(logT4) * weighted;
            }

            entropy.put(cb, localSplitEntropy(cb) + e);
        }
        return entropy.get(getRootClade().getCladeInBits());
    }

    private static double binomial(int n, int k) {
        double r = 1.0;
        for (int i = 1; i <= k; i++) {
            r = r * (n - k + i) / i;
        }
        return r;
    }

    /** Draws used by {@link #getEntropy()}. */
    public static final int DEFAULT_ENTROPY_SAMPLES = 100_000;

    /**
     * Unbiased Monte-Carlo entropy in nats.
     *
     * @param samples number of draws
     * @return {@code {estimate, standard error}}
     */
    public double[] getEntropyMonteCarlo(int samples) {
        double s1 = 0.0;
        double s2 = 0.0;
        for (int i = 0; i < samples; i++) {
            double logp = sampleTreeLogProbability();
            s1 += -logp;
            s2 += logp * logp;
        }
        double mean = s1 / samples;
        double se = Math.sqrt(Math.max(0.0, s2 / samples - mean * mean) / samples);
        return new double[]{mean, se};
    }

    /** Monte-Carlo entropy at {@link #DEFAULT_ENTROPY_SAMPLES} draws. */
    @Override
    public double getEntropy() {
        return getEntropyMonteCarlo(DEFAULT_ENTROPY_SAMPLES)[0];
    }

    /** Not applicable: the Lewis recursion assumes the distribution is supported on the CCD graph. */
    @Override
    public double getEntropyLewis() {
        throw new UnsupportedOperationException(
                "CRegCCD has support outside the CCD graph; use getEntropyMonteCarlo(samples).");
    }

    /* ----------------------------------------------------------------------
     * Sampling
     *
     * The generative process is the model read forwards: at each clade draw a class with
     * probability proportional to (its observed count + a_j) over the non-empty classes, then a
     * member uniformly within that class, then recurse into both children. Classes 1-3 are
     * explicitly enumerable in one pass over the observed clades; class 4 is drawn by rejection
     * from uniform bipartitions, which accepts with probability |A_4| / (2^(m-1) - 1) -- close to
     * 1 for any clade large enough for that to matter.
     *
     * Every drawn split is scored with the same {@link #logSplitProbability} the scorer uses, so
     * the sampling distribution equals exp(getLogProbabilityOfTree) by construction: there is no
     * truncation and no separate sampling fidelity to choose.
     * ------------------------------------------------------------------- */

    /** Simulates one draw and returns its log probability, without materialising a tree. */
    @Override
    public double sampleTreeLogProbability() {
        return simulate(getRootClade().getCladeInBits());
    }

    private double simulate(BitSet cBits) {
        if (cBits.cardinality() == 1) {
            return 0.0;
        }
        BitSet[] split = sampleSplit(cBits);
        return logSplitProbability(cBits, split[0], split[1], alpha, alpha1, alpha2)
                + simulate(split[0]) + simulate(split[1]);
    }

    /**
     * The inherited sampler only ever picks an observed clade partition, so it would draw from the
     * observed-splits-only distribution and could never produce a novel clade. Random sampling is
     * therefore overridden, as is MAP, which uses this model's own backbone DP rather than the
     * inherited CCD1 conditional clade probabilities.
     */
    @Override
    protected Node getVertexBasedOnStrategy(Clade clade, SamplingStrategy samplingStrategy,
                                            HeightSettingStrategy heightStrategy) {
        if (clade.isLeaf()) {
            return super.getVertexBasedOnStrategy(clade, samplingStrategy, heightStrategy);
        }
        if (samplingStrategy == SamplingStrategy.Sampling) {
            return sampleVertex(clade.getCladeInBits(), heightStrategy);
        }
        if (samplingStrategy == SamplingStrategy.MAP) {
            computeMAP();
            return mapVertex(clade.getCladeInBits(), heightStrategy);
        }
        return super.getVertexBasedOnStrategy(clade, samplingStrategy, heightStrategy);
    }

    /** Traceback of the backbone MAP DP. Every clade it visits is observed, by construction. */
    private Node mapVertex(BitSet cBits, HeightSettingStrategy heightStrategy) {
        if (cBits.cardinality() == 1) {
            return super.getVertexBasedOnStrategy(getClade(cBits),
                    SamplingStrategy.MAP, heightStrategy);
        }
        BitSet[] split = mapArg.get(cBits);
        if (split == null) {
            throw new AssertionError("no backbone split for clade " + cBits);
        }
        Node left = mapVertex(split[0], heightStrategy);
        Node right = mapVertex(split[1], heightStrategy);
        double logFactor = logSplitProbability(cBits, split[0], split[1], alpha, alpha1, alpha2);
        return buildVertex(cBits, left, right, logFactor, heightStrategy);
    }

    private Node sampleVertex(BitSet cBits, HeightSettingStrategy heightStrategy) {
        if (cBits.cardinality() == 1) {
            // leaves are always observed clades, so the inherited leaf construction applies
            return super.getVertexBasedOnStrategy(getClade(cBits),
                    SamplingStrategy.Sampling, heightStrategy);
        }
        BitSet[] split = sampleSplit(cBits);
        Node left = sampleVertex(split[0], heightStrategy);
        Node right = sampleVertex(split[1], heightStrategy);
        double logFactor = logSplitProbability(cBits, split[0], split[1], alpha, alpha1, alpha2);
        return buildVertex(cBits, left, right, logFactor, heightStrategy);
    }

    /** Assembles an internal node from two resolved children, stamping the subtree probability. */
    private Node buildVertex(BitSet cBits, Node left, Node right, double logFactor,
                             HeightSettingStrategy heightStrategy) {
        Node vertex = new Node();
        vertex.setNr(nextRunningInnerIndex());
        vertex.addChild(left);
        vertex.addChild(right);

        Clade observed = getClade(cBits);
        double support = (observed != null) ? observed.getProbability() : 0.0;
        vertex.setMetaData(CLADE_SUPPORT_KEY, support);
        String posteriorSupport = CLADE_SUPPORT_KEY + "=" + support;
        vertex.metaDataString = (vertex.metaDataString != null)
                ? vertex.metaDataString + "," + posteriorSupport : posteriorSupport;

        double logP = (Double) left.getMetaData(LOG_PROB_SUBTREE_KEY)
                + (Double) right.getMetaData(LOG_PROB_SUBTREE_KEY) + logFactor;
        vertex.setMetaData(LOG_PROB_SUBTREE_KEY, logP);
        vertex.setMetaData(PROB_SUBTREE_KEY, Math.exp(logP));

        setSampledHeight(vertex, left, right, observed, heightStrategy);
        return vertex;
    }

    /** Heights: {@code One} stacks by one; the height strategies use the clade's recorded height
     *  when it is available and strictly above both children, else a minimal positive increment. */
    private void setSampledHeight(Node vertex, Node left, Node right, Clade observed,
                                  HeightSettingStrategy heightStrategy) {
        if (heightStrategy == null || heightStrategy == HeightSettingStrategy.None) {
            return;
        }
        double maxChild = Math.max(left.getHeight(), right.getHeight());
        if (heightStrategy == HeightSettingStrategy.One) {
            vertex.setHeight(maxChild + 1.0);
            return;
        }
        double recorded = Double.NaN;
        if (observed != null) {
            recorded = (heightStrategy == HeightSettingStrategy.CommonAncestorHeights)
                    ? observed.getCommonAncestorHeight() : observed.getMeanOccurredHeight();
        }
        vertex.setHeight(recorded > maxChild ? recorded : maxChild + NOVEL_HEIGHT_EPS);
    }

    /** Draws one bipartition of {@code cBits} from this model's conditional distribution. */
    private BitSet[] sampleSplit(BitSet cBits) {
        double[] size = classSizes(cBits);
        Clade c = getClade(cBits);
        double[] wz = logWeightsAndZ(cBits, alpha, alpha1, alpha2);

        // total mass of each class: counts plus its share of the prior
        double[] weight = new double[4];
        for (int j = 0; j < 4; j++) {
            weight[j] = (size[j] > 0) ? Math.exp(wz[j]) * size[j] : 0.0;
        }
        if (size[0] > 0) {
            weight[0] += c.getNumberOfOccurrences();
        }

        double total = weight[0] + weight[1] + weight[2] + weight[3];
        double target = random().nextDouble() * total;
        int cls = -1;
        double acc = 0.0;
        for (int j = 0; j < 4; j++) {
            if (weight[j] <= 0) {
                continue;
            }
            acc += weight[j];
            if (target < acc) {
                cls = j;
                break;
            }
        }
        if (cls < 0) { // numerical guard: fall back to the last non-empty class
            for (int j = 3; j >= 0; j--) {
                if (size[j] > 0) {
                    cls = j;
                    break;
                }
            }
        }

        return switch (cls) {
            case 0 -> sampleObservedSplit(c, size[0]);
            case 1, 2 -> sampleEnumerableSplit(cBits, cls);
            default -> sampleNovelSplit(cBits, size[3]);
        };
    }

    /** Class 1: an observed split, with weight {@code f(S) + a_1/|A_1|}. */
    private BitSet[] sampleObservedSplit(Clade c, double n1) {
        double perSplit = Math.exp(logWeightsAndZ(c.getCladeInBits(), alpha, alpha1, alpha2)[0]);
        double totalWeight = c.getNumberOfOccurrences() + perSplit * n1;
        double target = random().nextDouble() * totalWeight;
        double acc = 0.0;
        CladePartition last = null;
        for (CladePartition p : c.getPartitions()) {
            if (p.getNumberOfOccurrences() == 0) {
                continue;
            }
            last = p;
            acc += p.getNumberOfOccurrences() + perSplit;
            if (target < acc) {
                return childBits(p);
            }
        }
        return childBits(last); // numerical guard
    }

    private static BitSet[] childBits(CladePartition p) {
        return new BitSet[]{p.getChildClades()[0].getCladeInBits(),
                p.getChildClades()[1].getCladeInBits()};
    }

    /**
     * Classes 2 and 3, drawn uniformly by reservoir sampling over the same single pass across the
     * observed clades inside {@code cBits} that produced the class sizes. A both-observed
     * bipartition is reached from either side, so it is only considered from its canonically
     * smaller side; a one-observed bipartition is reached exactly once, from its observed side.
     */
    private BitSet[] sampleEnumerableSplit(BitSet cBits, int cls) {
        int m = cBits.cardinality();
        int seen = 0;
        BitSet[] pick = null;
        for (BitSet d : sortedCladeBits()) {
            if (d.cardinality() >= m || !subset(d, cBits)) {
                continue;
            }
            BitSet complement = BitSet.newBitSet(cBits);
            complement.andNot(d);
            boolean complementObserved = getClade(complement) != null;
            boolean candidate;
            if (complementObserved) {
                candidate = cls == 1
                        && compareBitSets(d, complement) < 0
                        && observedPartition(getClade(cBits), d, complement) == null;
            } else {
                candidate = cls == 2;
            }
            if (candidate) {
                seen++;
                if (random().nextInt(seen) == 0) {
                    pick = new BitSet[]{d, complement};
                }
            }
        }
        return pick;
    }

    /**
     * Class 4, drawn uniformly among bipartitions with neither side an observed clade. Uniform
     * bipartitions are generated by pinning the lowest taxon to one side and flipping a fair coin
     * for the rest, and rejected unless both sides are novel. When acceptance would be poor the
     * bipartition set is necessarily small, so it is enumerated instead.
     */
    private BitSet[] sampleNovelSplit(BitSet cBits, double n4) {
        int m = cBits.cardinality();
        int[] idx = new int[m];
        int k = 0;
        for (int b = cBits.nextSetBit(0); b >= 0; b = cBits.nextSetBit(b + 1)) {
            idx[k++] = b;
        }
        double total = Math.pow(2.0, m - 1) - 1.0;

        if (n4 / total >= MIN_REJECTION_ACCEPTANCE) {
            while (true) {
                BitSet left = BitSet.newBitSet(leafArraySize);
                BitSet right = BitSet.newBitSet(leafArraySize);
                left.set(idx[0]);
                for (int i = 1; i < m; i++) {
                    if (random().nextBoolean()) {
                        left.set(idx[i]);
                    } else {
                        right.set(idx[i]);
                    }
                }
                if (right.isEmpty()) {
                    continue;
                }
                if (getClade(left) == null && getClade(right) == null) {
                    return new BitSet[]{left, right};
                }
            }
        }

        // low acceptance => 2^(m-1) is small; enumerate and reservoir-sample
        int seen = 0;
        BitSet[] pick = null;
        for (int mask = 0; mask < (1 << (m - 1)); mask++) {
            BitSet left = BitSet.newBitSet(leafArraySize);
            BitSet right = BitSet.newBitSet(leafArraySize);
            left.set(idx[0]);
            for (int i = 1; i < m; i++) {
                if ((mask & (1 << (i - 1))) != 0) {
                    left.set(idx[i]);
                } else {
                    right.set(idx[i]);
                }
            }
            if (right.isEmpty()) {
                continue;
            }
            if (getClade(left) == null && getClade(right) == null) {
                seen++;
                if (random().nextInt(seen) == 0) {
                    pick = new BitSet[]{left, right};
                }
            }
        }
        return pick;
    }

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
        return Integer.compare(ia, ib);
    }
}
