package ccd.model;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * GRegCCD -- the Geometrically-Regularised CCD: a one-parameter, full-support tree distribution in
 * which every <em>unobserved</em> clade-split is priced by a single escape rate {@code eps} (observed
 * splits keep their sample count). A tree's probability is
 * <pre>
 *   P(T) = eps^(s(T)) / Z,    s(T) = number of (clade -> split) pairs of T not seen in the sample,
 * </pre>
 * i.e. a geometric law in split-novelty distance. This unifies RegCCD's {@code alpha} split-expansion
 * (a recombination of observed clades is one new split) and KRegCCD's {@code mu} escape (a novel clade
 * costs its creating split plus its own split, i.e. >= 2 new splits) into a single rate, so there is no
 * {@code alpha} and no per-clade reserve depth.
 *
 * <p>Unlike KRegCCD it is <em>exactly</em> normalised, because {@code Z} is the true partition function
 * of the standard CCD sum-product
 * <pre>
 *   W(S) = sum over bipartitions {L,R} of S of  w(S -> L|R) * W(L) * W(R),   Z = W(rootClade),
 * </pre>
 * with {@code w} the observed count or {@code eps} for any unseen bipartition. A subset {@code S} that
 * contains no observed clade is <em>fresh</em>: every bipartition is priced {@code eps}, so {@code W(S)}
 * depends only on {@code |S|} and equals the universal {@code g(|S|)} computed once by
 * {@code g(m) = eps * (1/2) * sum_{j=1}^{m-1} C(m,j) g(j) g(m-j)}.
 *
 * <p>{@code Z} is computed exactly by a subset DP, which is {@code O(3^n)} and so feasible only up to
 * {@link #MAX_EXACT_TAXA} taxa. A tractable large-{@code n} approximation of {@code Z} (fresh-set term
 * plus low-order observed-clade corrections, mirroring KRegCCD's reserve decoupling) is future work; the
 * scoring of a given tree, {@code s(T)}, is always cheap.
 */
public class GRegCCD {

    /** Largest taxon count for which the exact O(3^n) partition function is computed. */
    public static final int MAX_EXACT_TAXA = 16;

    private final int nTaxa;
    private final double eps;
    private final double logEps;

    /** (cladeMask, minChildMask) -> observed split count. */
    private final Map<Long, Double> obsWeight = new HashMap<>();
    /** Masks of observed clades (split parents); used to detect fresh subsets. */
    private final long[] observedClades;

    private final double[] g;          // g[m] = W of a fresh m-set
    private final Map<Long, Double> wMemo = new HashMap<>();
    private final double logZ;

    public GRegCCD(List<Tree> trees, double eps) {
        if (eps <= 0 || eps >= 1) {
            throw new IllegalArgumentException("eps must be in (0, 1), got " + eps);
        }
        this.eps = eps;
        this.logEps = Math.log(eps);
        this.nTaxa = trees.get(0).getLeafNodeCount();

        Set<Long> cladeSet = new HashSet<>();
        for (Tree t : trees) {
            recordSplits(t.getRoot(), cladeSet);
        }
        this.observedClades = cladeSet.stream().mapToLong(Long::longValue).sorted().toArray();

        this.g = freshSetWeights(nTaxa, eps);

        long full = (nTaxa == 64) ? -1L : (1L << nTaxa) - 1;
        if (nTaxa > MAX_EXACT_TAXA) {
            throw new UnsupportedOperationException(
                    "GRegCCD exact normalisation is O(3^n); n=" + nTaxa + " exceeds MAX_EXACT_TAXA="
                            + MAX_EXACT_TAXA + ". A tractable approximate Z is not yet implemented.");
        }
        this.logZ = Math.log(W(full));
    }

    /* ---- construction: observed split weights ---- */

    private long recordSplits(Node v, Set<Long> cladeSet) {
        if (v.isLeaf()) {
            return 1L << v.getNr();
        }
        long l = recordSplits(v.getChildren().get(0), cladeSet);
        long r = recordSplits(v.getChildren().get(1), cladeSet);
        long mask = l | r;
        cladeSet.add(mask);
        obsWeight.merge(splitKey(mask, l, r), 1.0, Double::sum);
        return mask;
    }

    private static long splitKey(long mask, long childA, long childB) {
        return (mask << 32) | Math.min(childA, childB);
    }

    /* ---- universal fresh-set partition function g(m) ---- */

    private static double[] freshSetWeights(int n, double eps) {
        long[][] c = new long[n + 1][n + 1];
        for (int m = 0; m <= n; m++) {
            c[m][0] = 1;
            for (int k = 1; k <= m; k++) {
                c[m][k] = c[m - 1][k - 1] + (k <= m - 1 ? c[m - 1][k] : 0);
            }
        }
        double[] g = new double[n + 1];
        if (n >= 1) g[1] = 1.0;
        for (int m = 2; m <= n; m++) {
            double s = 0;
            for (int j = 1; j <= m - 1; j++) {
                s += c[m][j] * g[j] * g[m - j];
            }
            g[m] = eps * 0.5 * s;
        }
        return g;
    }

    /* ---- exact partition function over all rooted trees on `mask` ---- */

    private double W(long mask) {
        int size = Long.bitCount(mask);
        if (size == 1) {
            return 1.0;
        }
        if (isFresh(mask)) {
            return g[size];
        }
        Double memo = wMemo.get(mask);
        if (memo != null) {
            return memo;
        }
        long low = mask & (-mask);
        double sum = 0.0;
        for (long sub = mask; sub != 0; sub = (sub - 1) & mask) {
            if ((sub & low) == 0 || sub == mask) {
                continue; // L contains the lowest bit; R = mask ^ sub nonempty
            }
            long L = sub, R = mask ^ sub;
            double w = obsWeight.getOrDefault(splitKey(mask, L, R), eps);
            sum += w * W(L) * W(R);
        }
        wMemo.put(mask, sum);
        return sum;
    }

    /** A subset is fresh (W == g[size]) iff it contains no observed clade as a subset. */
    private boolean isFresh(long mask) {
        for (long c : observedClades) {
            if ((c & ~mask) == 0) { // c subset of mask
                return false;
            }
        }
        return true;
    }

    /* ---- scoring ---- */

    public double getLogProbabilityOfTree(Tree tree) {
        double[] logProd = {0.0};
        scoreSplits(tree.getRoot(), logProd);
        return logProd[0] - logZ;
    }

    public double getProbabilityOfTree(Tree tree) {
        return Math.exp(getLogProbabilityOfTree(tree));
    }

    /** Always true: every tree on this taxon set has positive probability. */
    public boolean containsTree(Tree tree) {
        return true;
    }

    /** Number of (clade -> split) pairs of the tree that were never seen in the sample. */
    public int novelSplitCount(Tree tree) {
        int[] count = {0};
        scoreSplits(tree.getRoot(), new double[1], count);
        return count[0];
    }

    private long scoreSplits(Node v, double[] logProd) {
        return scoreSplits(v, logProd, null);
    }

    private long scoreSplits(Node v, double[] logProd, int[] novel) {
        if (v.isLeaf()) {
            return 1L << v.getNr();
        }
        long l = scoreSplits(v.getChildren().get(0), logProd, novel);
        long r = scoreSplits(v.getChildren().get(1), logProd, novel);
        long mask = l | r;
        Double w = obsWeight.get(splitKey(mask, l, r));
        if (w == null) {
            logProd[0] += logEps;
            if (novel != null) novel[0]++;
        } else {
            logProd[0] += Math.log(w);
        }
        return mask;
    }

    public double getEpsilon() {
        return eps;
    }

    public double getLogPartitionFunction() {
        return logZ;
    }

    @Override
    public String toString() {
        return "GRegCCD [eps = " + eps + ", taxa = " + nTaxa + ", logZ = " + logZ + "]";
    }
}
