package ccd.model;

import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator.TreeSet;
import ccd.model.bitsets.BitSet;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

/**
 * The one-parameter "per-new-split" regularised CCD, with the boundary counts computed the way
 * {@link KRegCCD} computes its reserve rather than by direct recursive enumeration.
 *
 * <p>The model is unchanged: this overrides only {@code countsFor}, so every probability, sample and
 * point estimate is defined exactly as in {@link MRegCCDSlow}, which this overrides in a single
 * method. The two must therefore agree wherever the reference implementation's op budget does not
 * truncate its enumeration, which is what {@code MRegCCDAgreementTest} checks.
 *
 * <p>The reference implementation walks every ordered choice of boundary parts, which costs {@code O(m^(k-1))} in
 * the number {@code m} of observed subclades of a clade and blows through a flat op budget on large
 * analyses -- silently, since {@code countsFor} catches the overflow and leaves the remaining orders
 * at zero, which also disables the tail correction that reads them. Here the boundaries are found by
 * indexing instead:
 * <ul>
 *   <li>the disjoint pairs of observed subclades are enumerated once, {@code O(m^2)}, and grouped by
 *       the bitset they cover;</li>
 *   <li>a boundary of 2 is an observed subclade whose complement in {@code C} is observed;</li>
 *   <li>a boundary of 3 is an observed subclade whose complement is covered by a pair;</li>
 *   <li>a boundary of 4 is a pair whose complement is covered by another pair.</li>
 * </ul>
 * Each lookup is a hash probe rather than a search, so orders 2 to 4 cost {@code O(m^2)} in total.
 * Orders beyond 4 fall back to the inherited enumeration, so this is a strict speed-up of the
 * practical depths and never changes what is computed.
 */
public class MRegCCD extends MRegCCDSlow {

    private final Map<BitSet, int[]> fastCounts = new ConcurrentHashMap<>();

    public MRegCCD(List<Tree> trees, double burnin, double mu) {
        super(trees, burnin, mu);
    }

    public MRegCCD(List<Tree> trees, double burnin, double mu, int reserveDepth, boolean useTail) {
        super(trees, burnin, mu, reserveDepth, useTail);
    }

    public MRegCCD(TreeSet treeSet, double mu) {
        super(treeSet, mu);
    }

    public MRegCCD(TreeSet treeSet, double mu, int reserveDepth, boolean useTail) {
        super(treeSet, mu, reserveDepth, useTail);
    }

    @Override
    int[] countsFor(BitSet C) {
        int[] cached = fastCounts.get(C);
        if (cached != null) {
            return cached;
        }
        int card = C.cardinality();
        int depth = Math.min(card, getReserveDepth());
        // The pair index only pays for itself once order 4 needs it; below that the inherited
        // enumeration is cheaper, and identical by construction.
        if (depth < 4) {
            return super.countsFor(C);
        }
        int[] n = new int[depth + 1];
        if (card >= 2 && depth >= 2) {
            List<BitSet> subs = subclades(C);

            // Every disjoint pair of observed subclades, grouped by the bitset it covers. Only
            // orders 3 and 4 consult this, so at depth 2 the index is not worth building: the O(m^2)
            // pass would cost more than the order-2 scan it would serve.
            Map<BitSet, List<BitSet[]>> pairsByUnion = new HashMap<>();
            for (int i = 0; i < subs.size(); i++) {
                BitSet a = subs.get(i);
                for (int j = i + 1; j < subs.size(); j++) {
                    BitSet b = subs.get(j);
                    if (a.intersects(b)) {
                        continue;
                    }
                    BitSet union = BitSet.newBitSet(a);
                    union.or(b);
                    pairsByUnion.computeIfAbsent(union, k -> new ArrayList<>())
                            .add(new BitSet[]{a, b});
                }
            }

            for (BitSet d : subs) {
                BitSet rest = BitSet.newBitSet(C);
                rest.andNot(d);
                if (rest.isEmpty()) {
                    continue;
                }
                // boundary 2: {d, rest}, counted once from its canonically smaller side
                if (depth >= 2 && isObs(rest) && compareBitSets(d, rest) < 0) {
                    n[2] += countAllNovelResolutions(C, new BitSet[]{d, rest});
                }
                // boundary 3: {d} plus a pair covering the remainder
                if (depth >= 3) {
                    for (BitSet[] p : pairsByUnion.getOrDefault(rest, List.of())) {
                        if (compareBitSets(d, p[0]) < 0) {   // d must be the canonically first part
                            n[3] += countAllNovelResolutions(C, new BitSet[]{d, p[0], p[1]});
                        }
                    }
                }
            }

            // boundary 4: a pair whose complement is covered by another pair
            if (depth >= 4) {
                for (Map.Entry<BitSet, List<BitSet[]>> e : pairsByUnion.entrySet()) {
                    BitSet rest = BitSet.newBitSet(C);
                    rest.andNot(e.getKey());
                    if (rest.isEmpty() || !subset(e.getKey(), C)) {
                        continue;
                    }
                    List<BitSet[]> others = pairsByUnion.get(rest);
                    if (others == null) {
                        continue;
                    }
                    for (BitSet[] p : e.getValue()) {
                        for (BitSet[] q : others) {
                            // A 4-part boundary splits into two pairs in three ways, so count only
                            // the one whose first pair holds the two canonically smallest parts:
                            // for parts w < x < y < z that is {w,x}|{y,z} and no other.
                            if (compareBitSets(p[1], q[0]) < 0) {
                                n[4] += countAllNovelResolutions(C,
                                        new BitSet[]{p[0], p[1], q[0], q[1]});
                            }
                        }
                    }
                }
            }

            // orders beyond 4 are rare in practice; defer to the inherited enumeration
            if (depth >= 5) {
                int[] slow = super.countsFor(C);
                for (int m = 5; m < n.length && m < slow.length; m++) {
                    n[m] = slow[m];
                }
            }
        }
        fastCounts.put(BitSet.newBitSet(C), n);
        return n;
    }

    private static boolean subset(BitSet a, BitSet c) {
        BitSet tmp = BitSet.newBitSet(a);
        tmp.andNot(c);
        return tmp.isEmpty();
    }
}
