package ccd.model;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import ccd.model.bitsets.BitSet;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Tractable approximation of the GRegCCD partition function Z = W(rootClade) for large taxon sets,
 * where the exact O(3^n) sum is infeasible. Recurses over the OBSERVED clade DAG, pricing every
 * escape into a fresh remainder by the universal fresh-set weight g(m) and keeping recombinations
 * into two observed clades exact:
 * <pre>
 *   W(C) = g(|C|) + sum over observed clades A strict-subset of C of
 *          [ w(C->{A, C\A}) * Wtilde(A) * Wtilde(C\A)  -  eps * g(|A|) * g(|C\A|) ],
 * </pre>
 * with w = observed split count (else eps), and Wtilde(X) = W(X) if X is an observed clade else
 * g(|X|). Exact for observed trees and one-recombination ("one new split") trees; deeper novelty
 * (subsets that contain an observed clade but are not themselves observed, e.g. an observed cherry
 * plus a stray taxon) is approximated as fresh, an O(eps^2)-per-clade error. All terms are
 * non-negative, so the recursion is carried in log space (log-sum-exp) without cancellation.
 *
 * <p>Validated against {@link GRegCCD}'s exact Z on small enumerable sets ({@code GRegZApproxTest}).
 */
public class GRegZApprox {

    private final int nTaxa;
    private final List<BitSet> clades;                 // observed clades (incl. root + singletons)
    private final java.util.Set<BitSet> cladeSet;
    private final Map<BitSet, Map<BitSet, Double>> obsSplit; // parent -> canonChild -> count
    private final BitSet root;

    private double[] logG;
    private double eps, logEps;
    private Map<BitSet, Double> memo;

    private GRegZApprox(int n, List<BitSet> clades, java.util.Set<BitSet> cladeSet,
                        Map<BitSet, Map<BitSet, Double>> obsSplit, BitSet root) {
        this.nTaxa = n;
        this.clades = clades;
        this.cladeSet = cladeSet;
        this.obsSplit = obsSplit;
        this.root = root;
    }

    public static GRegZApprox fromTrees(List<Tree> trees) {
        int n = trees.get(0).getLeafNodeCount();
        Map<BitSet, Map<BitSet, Double>> obsSplit = new HashMap<>();
        java.util.Set<BitSet> cladeSet = new java.util.HashSet<>();
        for (Tree t : trees) {
            Map<Node, BitSet> bits = new HashMap<>();
            computeBits(t.getRoot(), bits, n);
            for (Node v : t.getNodesAsArray()) {
                BitSet pb = bits.get(v);
                cladeSet.add(pb);
                if (v.isLeaf()) {
                    continue;
                }
                BitSet canon = canonChild(pb, bits.get(v.getChildren().get(0)), bits.get(v.getChildren().get(1)));
                obsSplit.computeIfAbsent(pb, k -> new HashMap<>()).merge(canon, 1.0, Double::sum);
            }
        }
        BitSet root = BitSet.newBitSet(n);
        root.set(0, n);
        List<BitSet> clades = new ArrayList<>(cladeSet);
        clades.sort((a, b) -> Integer.compare(a.cardinality(), b.cardinality()));
        return new GRegZApprox(n, clades, cladeSet, obsSplit, root);
    }

    /** Approximate log Z at the given escape rate. */
    public double logZ(double eps) {
        this.eps = eps;
        this.logEps = Math.log(eps);
        this.logG = freshLogWeights(nTaxa, logEps);
        this.memo = new HashMap<>();
        return logW(root);
    }

    private double logW(BitSet C) {
        int size = C.cardinality();
        if (size == 1) {
            return 0.0;
        }
        Double m = memo.get(C);
        if (m != null) {
            return m;
        }
        int low = C.nextSetBit(0);
        double acc = logG[size];
        Map<BitSet, Double> splits = obsSplit.get(C);
        for (BitSet A : clades) {
            int aCard = A.cardinality();
            if (aCard >= size) {
                continue; // strict subset only (clades sorted ascending, but cheap to test)
            }
            if (!subset(A, C)) {
                continue;
            }
            BitSet B = BitSet.newBitSet(C);
            B.andNot(A);
            if (B.isEmpty()) {
                continue;
            }
            boolean aHasLow = A.get(low);
            boolean bInD = cladeSet.contains(B);
            if (!aHasLow && bInD) {
                continue; // bipartition will be (or was) counted from B's side
            }
            BitSet canon = aHasLow ? A : B;
            Double cnt = (splits == null) ? null : splits.get(canon);
            double w = (cnt != null) ? cnt : eps;

            int bCard = B.cardinality();
            double logWA = logW(A);
            double logWB = bInD ? logW(B) : logG[bCard];
            double logA = Math.log(w) + logWA + logWB;
            double logB = logEps + logG[aCard] + logG[bCard];
            if (logA > logB) {
                // log(exp(logA) - exp(logB)) stably for logA >= logB (expm1(logA-logB) overflows)
                double logTerm = logA + Math.log1p(-Math.exp(logB - logA));
                acc = logSumExp(acc, logTerm);
            }
        }
        memo.put(C, acc);
        return acc;
    }

    /* ---- helpers ---- */

    /** log g(m): g(m) = eps^(m-1) * (2m-3)!! (the fresh / uniform partition function on m taxa). */
    private static double[] freshLogWeights(int n, double logEps) {
        double[] logG = new double[n + 1];
        double logDoubleFact = 0.0; // log of (2m-3)!! accumulated
        if (n >= 1) logG[1] = 0.0;
        for (int m = 2; m <= n; m++) {
            logDoubleFact += Math.log(2 * m - 3); // multiply in the next odd factor
            logG[m] = (m - 1) * logEps + logDoubleFact;
        }
        return logG;
    }

    private static boolean subset(BitSet a, BitSet c) {
        BitSet tmp = BitSet.newBitSet(a);
        tmp.andNot(c);
        return tmp.isEmpty();
    }

    private static double logSumExp(double a, double b) {
        if (a == Double.NEGATIVE_INFINITY) return b;
        if (b == Double.NEGATIVE_INFINITY) return a;
        double max = Math.max(a, b);
        return max + Math.log(Math.exp(a - max) + Math.exp(b - max));
    }

    private static BitSet canonChild(BitSet parent, BitSet c0, BitSet c1) {
        int lb = parent.nextSetBit(0);
        return c0.get(lb) ? c0 : c1;
    }

    private static BitSet computeBits(Node v, Map<Node, BitSet> bits, int n) {
        BitSet b = BitSet.newBitSet(n);
        if (v.isLeaf()) {
            b.set(v.getNr());
        } else {
            b.or(computeBits(v.getChildren().get(0), bits, n));
            b.or(computeBits(v.getChildren().get(1), bits, n));
        }
        bits.put(v, b);
        return b;
    }
}
