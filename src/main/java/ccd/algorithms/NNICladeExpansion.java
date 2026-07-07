package ccd.algorithms;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator;
import ccd.model.AbstractCCD;
import ccd.model.CCD0;
import ccd.model.Clade;
import ccd.model.CladePartition;
import ccd.model.bitsets.BitSet;
import ccd.tools.CCDToolUtil;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Clade-level NNI expansion of a CCD: the analogue of the CCD0 split expansion
 * ({@link ccd.model.CCD0#expand()}), but operating on the <em>clades</em>
 * themselves rather than on clade partitions.
 *
 * <p>
 * The CCD0 expansion adds clade partitions {@code P -> {L, R}} that were never
 * observed but whose parent {@code P} and children {@code L}, {@code R} all are.
 * This class instead enumerates novel <em>clades</em> that a single
 * nearest-neighbour-interchange (NNI) move can produce from the splits already
 * stored in the CCD, without ever building a tree.
 * </p>
 *
 * <p>
 * <b>The move at the clade level.</b> An internal edge of a tree is a stored
 * split {@code P -> {L, R}}: the parent clade {@code P} branches into children
 * {@code L} and {@code R}. An NNI at that edge needs the grandchildren: pick a
 * child to descend into, say {@code L -> {A, B}}, with sibling {@code R}. The
 * two NNI moves recombine {@code {A, B, R}} so that the child of {@code P} is no
 * longer {@code L = A ∪ B} but either {@code A ∪ R} or {@code B ∪ R}. Those
 * unions are exactly the novel clades; each comes with the NNI-implied seed
 * split into {@code {kept, sibling}} (e.g. {@code (A ∪ R) -> {A, R}}), whose two
 * constituents are themselves existing clades. We read {@code A}, {@code B},
 * {@code R} straight out of the split table.
 * </p>
 *
 * <p>
 * <b>Pairing modes.</b> The CCD graph mixes splits drawn from many posterior
 * trees, so pairing a split {@code P -> {L, R}} with a grandchild split
 * {@code L -> {A, B}} can combine splits that never co-occurred in one tree.
 * The {@link PairingMode} controls which pairs are allowed:
 * <ul>
 * <li>{@link PairingMode#GRAPH_WIDE} pairs <em>any</em> stored split on
 * {@code P} with <em>any</em> stored split on its child {@code L}. This is the
 * cheapest, broadest reading; it generates clades that no single sampled tree is
 * one NNI move away from.</li>
 * <li>{@link PairingMode#CO_OCCURRING} only recombines grandchild splits that
 * actually appeared together under {@code P} in some single base tree, i.e. the
 * true NNI neighbourhood of the posterior sample. It walks the stored base trees
 * and therefore requires the CCD to have been built with stored base trees
 * (currently non-sampled-ancestor CCDs only).</li>
 * </ul>
 * </p>
 *
 * <p>
 * The pass is read-only: it does not modify the CCD. Cost is proportional to the
 * size of the CCD itself ({@code O(|splits| × d̄)}, with {@code d̄} the average
 * number of splits per clade), not to the number of trees the CCD represents.
 * </p>
 *
 * @author Claude (CCD-Sophie)
 */
public class NNICladeExpansion {

    /** How splits may be paired to form NNI recombinations. */
    public enum PairingMode {
        /** Pair any stored split on P with any stored split on its child L. */
        GRAPH_WIDE,
        /** Only pair splits that co-occurred under P in some single base tree. */
        CO_OCCURRING
    }

    /**
     * One way a novel clade can be produced by an NNI move. The novel clade's
     * taxa are {@code kept ∪ sibling}; its NNI-implied seed split is
     * {@code {kept, sibling}}; it appears as a new child of {@code parent} in
     * place of the original child {@code kept ∪ dropped}.
     *
     * @param parent  the parent clade P at the rearranged edge
     * @param sibling the sibling clade R that joins the kept grandchild
     * @param kept    the grandchild (A or B) that joins the sibling
     * @param dropped the other grandchild, left behind under P
     */
    public record Provenance(Clade parent, Clade sibling, Clade kept, Clade dropped) {
    }

    /** A candidate clade not present in the CCD, with all the ways it arises. */
    public static final class NovelClade {
        /** Taxa of the novel clade (taxa-only bits). */
        public final BitSet taxa;
        /** Number of taxa in the novel clade. */
        public final int size;
        /** Every (parent, sibling, kept, dropped) that produces this clade. */
        public final List<Provenance> provenances = new ArrayList<>(2);

        private NovelClade(BitSet taxa) {
            this.taxa = taxa;
            this.size = taxa.cardinality();
        }
    }

    private final AbstractCCD ccd;
    private final PairingMode mode;

    private final Map<BitSet, NovelClade> novelClades = new HashMap<>();
    private long candidatesEmitted = 0;
    private long candidatesAlreadyPresent = 0;
    private boolean computed = false;

    /** Taxa-only bits of every clade in the CCD, for fast presence tests. */
    private Set<BitSet> existingTaxa;

    /**
     * @param ccd  the CCD whose clade set is expanded (not modified)
     * @param mode how splits may be paired into NNI recombinations
     */
    public NNICladeExpansion(AbstractCCD ccd, PairingMode mode) {
        this.ccd = ccd;
        this.mode = mode;
    }

    /** Runs the expansion if not already done; idempotent. */
    public NNICladeExpansion compute() {
        if (computed) {
            return this;
        }

        existingTaxa = new HashSet<>(ccd.getNumberOfClades() * 2);
        for (Clade clade : ccd.getClades()) {
            existingTaxa.add(clade.getCladeInBitsTaxaOnly());
        }

        switch (mode) {
            case GRAPH_WIDE -> expandGraphWide();
            case CO_OCCURRING -> expandCoOccurring();
        }

        computed = true;
        return this;
    }

    /* -- GRAPH-WIDE EXPANSION (pair any stored adjacent splits) -- */

    private void expandGraphWide() {
        for (Clade parent : ccd.getClades()) {
            // need at least 3 taxa to have a (cherry-or-larger child) + sibling
            if (parent.size() < 3) {
                continue;
            }
            for (CladePartition split : parent.getPartitions()) {
                Clade childA = split.getChildClades()[0];
                Clade childB = split.getChildClades()[1];
                // descend into each child, with the other child as the sibling
                recombineAllGrandchildSplits(parent, childA, childB);
                recombineAllGrandchildSplits(parent, childB, childA);
            }
        }
    }

    /* For parent P, descend-into child L and sibling R; recombine every L-split. */
    private void recombineAllGrandchildSplits(Clade parent, Clade descend, Clade sibling) {
        for (CladePartition grandSplit : descend.getPartitions()) {
            Clade a = grandSplit.getChildClades()[0];
            Clade b = grandSplit.getChildClades()[1];
            emit(parent, sibling, a, b); // novel = A ∪ R, drop B
            emit(parent, sibling, b, a); // novel = B ∪ R, drop A
        }
    }

    /* -- CO-OCCURRING EXPANSION (only splits seen together in one base tree) -- */

    private void expandCoOccurring() {
        List<Tree> baseTrees = ccd.getBaseTrees();
        if (baseTrees == null) {
            throw new IllegalStateException(
                    "CO_OCCURRING pairing needs all base trees stored; "
                            + "build the CCD with storeBaseTrees = true.");
        }
        for (Tree tree : baseTrees) {
            Map<Node, Clade> nodeToClade = new HashMap<>();
            mapNodesToClades(tree.getRoot(), nodeToClade);
            recombineTree(tree.getRoot(), nodeToClade);
        }
    }

    /* Map every node of a base tree to its CCD clade via taxa bits. */
    private BitSet mapNodesToClades(Node vertex, Map<Node, Clade> map) {
        BitSet bits = BitSet.newBitSet(ccd.getSizeOfLeavesArray());
        if (vertex.isLeaf()) {
            bits.set(vertex.getNr());
        } else {
            for (Node child : vertex.getChildren()) {
                bits.or(mapNodesToClades(child, map));
            }
        }
        Clade clade = ccd.getClade(bits);
        if (clade == null) {
            throw new UnsupportedOperationException(
                    "No clade found for a base-tree node; CO_OCCURRING pairing "
                            + "currently supports non-sampled-ancestor CCDs only.");
        }
        map.put(vertex, clade);
        return bits;
    }

    /* Walk internal edges of one base tree, emitting their NNI recombinations. */
    private void recombineTree(Node vertex, Map<Node, Clade> map) {
        if (vertex.isLeaf()) {
            return;
        }
        List<Node> children = vertex.getChildren();
        if (children.size() == 2) {
            Node left = children.get(0);
            Node right = children.get(1);
            Clade parent = map.get(vertex);
            // edge above an internal child = NNI site; the other child is the sibling
            if (!left.isLeaf()) {
                recombineAt(parent, left, right, map);
            }
            if (!right.isLeaf()) {
                recombineAt(parent, right, left, map);
            }
        }
        for (Node child : children) {
            recombineTree(child, map);
        }
    }

    /* The two NNI moves at the edge (parent, descend), with given sibling. */
    private void recombineAt(Clade parent, Node descend, Node siblingNode, Map<Node, Clade> map) {
        if (descend.getChildren().size() != 2) {
            return;
        }
        Clade sibling = map.get(siblingNode);
        Clade a = map.get(descend.getChildren().get(0));
        Clade b = map.get(descend.getChildren().get(1));
        emit(parent, sibling, a, b);
        emit(parent, sibling, b, a);
    }

    /* -- SHARED EMIT -- */

    /* Record the candidate clade kept ∪ sibling; novel iff not already in C. */
    private void emit(Clade parent, Clade sibling, Clade kept, Clade dropped) {
        BitSet taxa = (BitSet) kept.getCladeInBitsTaxaOnly().clone();
        taxa.or(sibling.getCladeInBitsTaxaOnly());

        candidatesEmitted++;
        if (existingTaxa.contains(taxa)) {
            candidatesAlreadyPresent++;
            return;
        }
        NovelClade novel = novelClades.computeIfAbsent(taxa, NovelClade::new);
        novel.provenances.add(new Provenance(parent, sibling, kept, dropped));
    }

    /* -- RESULTS -- */

    /** @return the distinct novel candidate clades (not present in the CCD) */
    public Collection<NovelClade> getNovelClades() {
        compute();
        return novelClades.values();
    }

    /** @return number of distinct novel clades found */
    public int getNumberOfNovelClades() {
        compute();
        return novelClades.size();
    }

    /** @return total candidate emissions, counting multiplicity */
    public long getNumberOfCandidatesEmitted() {
        compute();
        return candidatesEmitted;
    }

    /** @return emissions whose taxa already form a clade in the CCD */
    public long getNumberOfCandidatesAlreadyPresent() {
        compute();
        return candidatesAlreadyPresent;
    }

    /**
     * Number of tree topologies the CCD would gain if this NNI expansion were
     * applied. Computed without mutating the CCD.
     *
     * <p>
     * Equivalent to running {@link ccd.algorithms.regularisation.NNICladeExpander}
     * and then asking the resulting CCD for {@link AbstractCCD#getNumberOfTrees()}
     * minus the pre-expansion count, but skipping the graph mutation. The
     * read-only enumeration must still run (via {@link #compute()}); the saving
     * is in the mutation step itself.
     * </p>
     *
     * @return T_after - T_before, the count of new trees the expansion adds
     */
    public BigInteger getNumberOfNovelTrees() {
        return getNumberOfTreesAfterExpansion().subtract(ccd.getNumberOfTrees());
    }

    /**
     * Number of tree topologies the CCD would contain after this NNI expansion
     * were applied, computed without mutating the CCD.
     *
     * <p>
     * The expansion only ever (a) adds partitions {@code P -> {K, B}} on
     * existing parents {@code P}, with {@code K} a novel clade and {@code B}
     * existing; and (b) gives each novel {@code K} seed splits {@code K ->
     * {kept, sibling}} whose children are both existing. So the dependency
     * graph is still a DAG and a single size-ordered pass settles all
     * topology counts.
     * </p>
     *
     * @return the post-expansion topology count at the root clade
     */
    public BigInteger getNumberOfTreesAfterExpansion() {
        compute();

        // (existing parent P) -> (novel K -> set of distinct dropped clades B);
        // mirrors NNICladeExpander's per-partition de-duplication
        Map<Clade, Map<NovelClade, Set<Clade>>> novelChildrenByParent = new HashMap<>();
        // (novel K) -> set of distinct unordered {kept, sibling} seed splits
        Map<NovelClade, Set<Set<Clade>>> seedSplitsByNovel = new HashMap<>();

        for (NovelClade k : novelClades.values()) {
            for (Provenance prov : k.provenances) {
                Set<Clade> seed = new HashSet<>(2);
                seed.add(prov.kept());
                seed.add(prov.sibling());
                seedSplitsByNovel.computeIfAbsent(k, x -> new HashSet<>()).add(seed);

                novelChildrenByParent
                        .computeIfAbsent(prov.parent(), p -> new HashMap<>())
                        .computeIfAbsent(k, x -> new HashSet<>())
                        .add(prov.dropped());
            }
        }

        // T_new(c) for every existing clade c, and T_new(K) for every novel K.
        // We recompute the topology product over the *existing* partitions of c
        // using updated child counts - we cannot reuse c.getNumberOfTopologies()
        // here, because a novel partition deeper in the tree can change the
        // count of an existing child of c.
        Map<Clade, BigInteger> tExisting = new HashMap<>(ccd.getNumberOfClades() * 2);
        Map<NovelClade, BigInteger> tNovel = new HashMap<>(novelClades.size() * 2);

        List<Clade> existingBySize = new ArrayList<>(ccd.getClades());
        existingBySize.sort(Comparator.comparingInt(Clade::size));
        List<NovelClade> novelBySize = new ArrayList<>(novelClades.values());
        novelBySize.sort(Comparator.comparingInt(n -> n.size));

        // merged size-ascending walk; within a size no novel/existing pair
        // depends on the other, so equal-size ordering is arbitrary
        int ei = 0;
        int ni = 0;
        while (ei < existingBySize.size() || ni < novelBySize.size()) {
            int eSize = (ei < existingBySize.size())
                    ? existingBySize.get(ei).size() : Integer.MAX_VALUE;
            int nSize = (ni < novelBySize.size())
                    ? novelBySize.get(ni).size : Integer.MAX_VALUE;
            if (nSize <= eSize) {
                NovelClade k = novelBySize.get(ni++);
                BigInteger t = BigInteger.ZERO;
                for (Set<Clade> seed : seedSplitsByNovel.get(k)) {
                    Iterator<Clade> it = seed.iterator();
                    Clade a = it.next();
                    Clade b = it.next();
                    t = t.add(tExisting.get(a).multiply(tExisting.get(b)));
                }
                tNovel.put(k, t);
            } else {
                Clade c = existingBySize.get(ei++);
                BigInteger t;
                if (c.isLeaf()) {
                    t = BigInteger.ONE;
                } else {
                    t = BigInteger.ZERO;
                    for (CladePartition p : c.getPartitions()) {
                        BigInteger left = tExisting.get(p.getChildClades()[0]);
                        BigInteger right = tExisting.get(p.getChildClades()[1]);
                        t = t.add(left.multiply(right));
                    }
                }
                Map<NovelClade, Set<Clade>> kids = novelChildrenByParent.get(c);
                if (kids != null) {
                    for (Map.Entry<NovelClade, Set<Clade>> e : kids.entrySet()) {
                        BigInteger tK = tNovel.get(e.getKey());
                        for (Clade dropped : e.getValue()) {
                            t = t.add(tK.multiply(tExisting.get(dropped)));
                        }
                    }
                }
                tExisting.put(c, t);
            }
        }

        return tExisting.get(ccd.getRootClade());
    }

    /**
     * Number of tree topologies the CCD would contain after this NNI clade
     * expansion <em>and</em> a subsequent {@link
     * ccd.algorithms.regularisation.CCDExpansion} pass (i.e. the
     * {@link ccd.model.NNIRegCCD} pipeline graph), computed without mutating
     * the CCD.
     *
     * <p>
     * After both passes every clade {@code P} (observed or NNI-novel) carries
     * every compatible split {@code (L, R)} - {@code L ∪ R = P},
     * {@code L ∩ R = ∅}, both {@code L} and {@code R} clades in the universe.
     * So we enumerate compatible splits directly over the union of existing
     * and novel taxa-sets.
     * </p>
     *
     * @return the post-(NNI + split-expansion) topology count at the root
     */
    public BigInteger getNumberOfTreesAfterFullExpansion() {
        compute();

        // universe (existing + novel), bucketed by clade size for the split scan
        Map<Integer, List<BitSet>> bySize = new HashMap<>();
        int maxSize = 0;
        for (Clade c : ccd.getClades()) {
            bySize.computeIfAbsent(c.size(), k -> new ArrayList<>())
                    .add(c.getCladeInBitsTaxaOnly());
            if (c.size() > maxSize) maxSize = c.size();
        }
        for (NovelClade nk : novelClades.values()) {
            bySize.computeIfAbsent(nk.size, k -> new ArrayList<>()).add(nk.taxa);
            if (nk.size > maxSize) maxSize = nk.size;
        }

        Map<BitSet, BigInteger> tMap = new HashMap<>(
                (ccd.getNumberOfClades() + novelClades.size()) * 2);

        for (int size = 1; size <= maxSize; size++) {
            List<BitSet> here = bySize.get(size);
            if (here == null) continue;
            for (BitSet pBits : here) {
                if (size == 1) {
                    tMap.put(pBits, BigInteger.ONE);
                    continue;
                }
                BigInteger t = BigInteger.ZERO;
                int halfSize = size / 2;
                for (int j = 1; j <= halfSize; j++) {
                    List<BitSet> bucket = bySize.get(j);
                    if (bucket == null) continue;
                    int complementSize = size - j;
                    for (BitSet lBits : bucket) {
                        if (!pBits.contains(lBits)) continue;
                        BitSet rBits = (BitSet) pBits.clone();
                        rBits.andNot(lBits);
                        BigInteger tR = tMap.get(rBits);
                        if (tR == null) continue;
                        if (j < complementSize) {
                            t = t.add(tMap.get(lBits).multiply(tR));
                        } else if (lBits.nextSetBit(0) < rBits.nextSetBit(0)) {
                            // even split: L and R both appear in this bucket -
                            // use first-bit ordering to count each unordered
                            // partition exactly once
                            t = t.add(tMap.get(lBits).multiply(tR));
                        }
                    }
                }
                tMap.put(pBits, t);
            }
        }

        return tMap.get(ccd.getRootClade().getCladeInBitsTaxaOnly());
    }

    /**
     * @return T_full - T_before, the number of trees that the full
     *         NNI + split expansion pipeline adds on top of the input CCD
     */
    public BigInteger getNumberOfNovelTreesAfterFullExpansion() {
        return getNumberOfTreesAfterFullExpansion().subtract(ccd.getNumberOfTrees());
    }

    /**
     * Total probability mass that the {@link ccd.model.NNIRegCCD} pipeline
     * (NNI clade expansion + {@link ccd.algorithms.regularisation.CCDExpansion}
     * + additive-{@code alpha} regularisation) would assign to trees containing
     * at least one novel clade. Computed without mutating the CCD.
     *
     * <p>
     * Equivalent to {@code getProbabilityOfNNIExpandedTrees(alpha, alpha)}
     * (single-level {@code AdditiveX} smoothing, matching {@code
     * regulariseRegCCD(alpha)}).
     * </p>
     *
     * @param alpha additive-smoothing pseudocount used on every partition
     * @return P(tree contains a novel clade) under the post-pipeline distribution
     */
    public double getProbabilityOfNNIExpandedTrees(double alpha) {
        return getProbabilityOfNNIExpandedTrees(alpha, alpha);
    }

    /**
     * Two-level additive smoothing variant: existing-and-split-expanded
     * partitions get pseudocount {@code alpha}; partitions touching a novel
     * clade get pseudocount {@code beta} (matching {@code
     * CCDRegularisationStrategy.AdditiveXY}).
     *
     * <p>
     * Probability is computed via {@code 1 - M(root)}, where {@code M(c)} is
     * the post-regularisation total CCP-product mass on trees rooted at
     * {@code c} that use only existing clades. {@code M} obeys
     * {@code M(leaf) = 1}, {@code M(novel) = 0}, and for existing {@code c}
     * with at least 2 partitions in the full graph:
     * {@code M(c) = Σ_p CCP(p) · M(L) · M(R)} over the existing-only
     * partitions {@code p = (L, R)}, with CCPs read off the post-regularisation
     * formula. Single-partition clades retain {@code CCP = 1} (the regulariser
     * skips them).
     * </p>
     *
     * @param alpha pseudocount for existing/split-expanded partitions
     * @param beta  pseudocount for NNI-derived (novel-touching) partitions
     * @return P(tree contains a novel clade) under the post-pipeline distribution
     */
    public double getProbabilityOfNNIExpandedTrees(double alpha, double beta) {
        compute();

        // universe (existing + novel) for the compatible-split count
        Map<Integer, List<BitSet>> bySize = new HashMap<>();
        Set<BitSet> universe = new HashSet<>();
        for (Clade c : ccd.getClades()) {
            BitSet bits = c.getCladeInBitsTaxaOnly();
            bySize.computeIfAbsent(c.size(), k -> new ArrayList<>()).add(bits);
            universe.add(bits);
        }
        for (NovelClade nk : novelClades.values()) {
            bySize.computeIfAbsent(nk.size, k -> new ArrayList<>()).add(nk.taxa);
            universe.add(nk.taxa);
        }

        // M(c) bottom-up over existing clades (novel clades have M = 0)
        Map<Clade, Double> mExisting = new HashMap<>(ccd.getNumberOfClades() * 2);
        List<Clade> existingBySize = new ArrayList<>(ccd.getClades());
        existingBySize.sort(Comparator.comparingInt(Clade::size));

        for (Clade c : existingBySize) {
            if (c.isLeaf()) {
                mExisting.put(c, 1.0);
                continue;
            }
            int nExisting = c.getNumberOfPartitions();
            int nTotal = countCompatibleSplits(c.getCladeInBitsTaxaOnly(),
                    c.size(), bySize, universe);
            int nNNI = nTotal - nExisting;

            double m = 0.0;
            if (nTotal <= 1 || c.isCherry()) {
                // regulariser skips: CCP of the single partition is 1
                CladePartition only = c.getPartitions().get(0);
                m = mExisting.get(only.getChildClades()[0])
                        * mExisting.get(only.getChildClades()[1]);
            } else {
                double denom = c.getNumberOfOccurrences()
                        + nExisting * alpha + nNNI * beta;
                for (CladePartition p : c.getPartitions()) {
                    double ccp = (p.getNumberOfOccurrences() + alpha) / denom;
                    Clade l = p.getChildClades()[0];
                    Clade r = p.getChildClades()[1];
                    m += ccp * mExisting.get(l) * mExisting.get(r);
                }
            }
            mExisting.put(c, m);
        }

        return 1.0 - mExisting.get(ccd.getRootClade());
    }

    private static int countCompatibleSplits(BitSet pBits, int pSize,
            Map<Integer, List<BitSet>> bySize, Set<BitSet> universe) {
        int count = 0;
        int halfSize = pSize / 2;
        for (int j = 1; j <= halfSize; j++) {
            List<BitSet> bucket = bySize.get(j);
            if (bucket == null) continue;
            int complementSize = pSize - j;
            for (BitSet lBits : bucket) {
                if (!pBits.contains(lBits)) continue;
                BitSet rBits = (BitSet) pBits.clone();
                rBits.andNot(lBits);
                if (!universe.contains(rBits)) continue;
                if (j < complementSize) {
                    count++;
                } else if (lBits.nextSetBit(0) < rBits.nextSetBit(0)) {
                    count++;
                }
            }
        }
        return count;
    }

    /** Where a held-out tree sits relative to this CCD's NNI universe. */
    public enum TreeZone {
        /** Every internal clade of the tree is observed in the input CCD. */
        OBSERVED,
        /** All internal clades are in the universe (existing ∪ NNI-novel),
         *  and at least one is NNI-novel. */
        NNI_EXPANDED,
        /** At least one internal clade is outside the universe; even the
         *  NNI-expanded graph cannot give this tree positive probability. */
        OUT_OF_UNIVERSE
    }

    /** Counts of held-out trees in each zone. */
    public record ZoneCounts(int observed, int nniExpanded, int outOfUniverse) {
        public int total() {
            return observed + nniExpanded + outOfUniverse;
        }
    }

    /**
     * Classify each held-out tree by whether it sits inside the input CCD's
     * observed-clade set, only inside the NNI-expanded universe, or contains a
     * clade outside the universe altogether.
     *
     * <p>
     * Assumes the held-out trees use the same leaf indexing as the CCD
     * (consistent BEAST translate blocks across independent runs of the same
     * analysis - the usual case).
     * </p>
     *
     * @param trees independent trees to classify against this NNI universe
     * @return counts of trees in each {@link TreeZone}
     */
    public ZoneCounts classifyTrees(Iterable<Tree> trees) {
        compute();
        Set<BitSet> novelTaxa = new HashSet<>(novelClades.size() * 2);
        for (NovelClade nk : novelClades.values()) {
            novelTaxa.add(nk.taxa);
        }

        int observed = 0;
        int nniExpanded = 0;
        int outOfUniverse = 0;
        int leafArraySize = ccd.getSizeOfLeavesArray();

        for (Tree tree : trees) {
            List<BitSet> internals = new ArrayList<>(tree.getNodeCount() / 2);
            collectInternalBits(tree.getRoot(), leafArraySize, internals);
            boolean anyNovel = false;
            boolean anyOut = false;
            for (BitSet bits : internals) {
                if (existingTaxa.contains(bits)) continue;
                if (novelTaxa.contains(bits)) {
                    anyNovel = true;
                    continue;
                }
                anyOut = true;
                break;
            }
            if (anyOut) outOfUniverse++;
            else if (anyNovel) nniExpanded++;
            else observed++;
        }
        return new ZoneCounts(observed, nniExpanded, outOfUniverse);
    }

    private static BitSet collectInternalBits(Node v, int size, List<BitSet> internals) {
        BitSet bits = BitSet.newBitSet(size);
        if (v.isLeaf()) {
            bits.set(v.getNr());
        } else {
            for (Node c : v.getChildren()) {
                bits.or(collectInternalBits(c, size, internals));
            }
            internals.add(bits);
        }
        return bits;
    }

    /** @return novel clades found per size; index i holds the count of size i */
    public int[] novelSizeHistogram() {
        compute();
        int[] hist = new int[ccd.getSizeOfLeavesArray() + 1];
        for (NovelClade novel : novelClades.values()) {
            hist[novel.size]++;
        }
        return hist;
    }

    /* -- REPORTING -- */

    /** @return a human-readable summary of this single expansion */
    public String report() {
        compute();
        int numClades = ccd.getNumberOfClades();
        int numNovel = novelClades.size();
        double blowup = (numClades == 0) ? 0 : (numNovel / (double) numClades);
        double presentFrac = (candidatesEmitted == 0) ? 0
                : (candidatesAlreadyPresent / (double) candidatesEmitted);

        StringBuilder sb = new StringBuilder();
        sb.append("NNI clade expansion [").append(mode).append("]\n");
        sb.append(String.format("  taxa                 : %d%n", ccd.getNumberOfLeaves()));
        sb.append(String.format("  base trees           : %d%n", ccd.getNumberOfBaseTrees()));
        sb.append(String.format("  clades in C          : %d%n", numClades));
        sb.append(String.format("  candidates emitted   : %d%n", candidatesEmitted));
        sb.append(String.format("  already in C         : %d (%.1f%%)%n",
                candidatesAlreadyPresent, 100 * presentFrac));
        sb.append(String.format("  novel clades         : %d%n", numNovel));
        sb.append(String.format("  blow-up (novel/|C|)  : %.3f%n", blowup));
        BigInteger treesBefore = ccd.getNumberOfTrees();
        BigInteger treesAfterNNI = getNumberOfTreesAfterExpansion();
        BigInteger treesAfterFull = getNumberOfTreesAfterFullExpansion();
        sb.append(String.format("  trees in C           : %s%n", treesBefore));
        sb.append(String.format("  trees after NNI      : %s%n", treesAfterNNI));
        sb.append(String.format("  trees after NNI+split: %s%n", treesAfterFull));
        sb.append(String.format("  novel trees (NNI)    : %s%n", treesAfterNNI.subtract(treesBefore)));
        sb.append(String.format("  novel trees (full)   : %s%n", treesAfterFull.subtract(treesBefore)));
        // P(NNI-expanded tree) varies with both α and β under AdditiveXY;
        // show a small grid covering the single-level default and a typical
        // fitted two-level setting where β << α dampens speculative mass
        sb.append("  P(novel tree) at (α, β):\n");
        sb.append(String.format("    (0.4 , 0.4 )  single-level default : %.6f%n",
                getProbabilityOfNNIExpandedTrees(0.4, 0.4)));
        sb.append(String.format("    (0.4 , 0.04)  β 10× damped         : %.6f%n",
                getProbabilityOfNNIExpandedTrees(0.4, 0.04)));
        sb.append(String.format("    (0.2 , 0.02)  recalled fit         : %.6f%n",
                getProbabilityOfNNIExpandedTrees(0.2, 0.02)));
        sb.append(String.format("    (0.265, 0.0014) GW held-out fit    : %.6f%n",
                getProbabilityOfNNIExpandedTrees(0.265, 0.0014)));
        sb.append(String.format("    (0.260, 0.0084) CO held-out fit    : %.6f%n",
                getProbabilityOfNNIExpandedTrees(0.260, 0.0084)));
        sb.append("  novel by size        :");
        int[] hist = novelSizeHistogram();
        for (int s = 2; s < hist.length; s++) {
            if (hist[s] > 0) {
                sb.append(' ').append(s).append(':').append(hist[s]);
            }
        }
        sb.append('\n');
        return sb.toString();
    }

    /**
     * Runs both pairing modes on the given CCD and reports each, plus how the
     * graph-wide novel set relates to the strict co-occurring one.
     *
     * @param ccd CCD to analyse (must have base trees stored for the strict mode)
     * @return a combined report
     */
    public static String compareModes(AbstractCCD ccd) {
        return compareModes(ccd, null);
    }

    /**
     * As {@link #compareModes(AbstractCCD)}, plus a held-out classification
     * section reporting the fraction of {@code heldOut} trees that fall inside
     * the observed-clade zone, the NNI-expanded zone, or outside the universe
     * altogether - in each pairing mode.
     *
     * @param ccd     CCD to analyse (training)
     * @param heldOut independent trees to score against the NNI universe;
     *                if {@code null} or empty, no classification section is added
     * @return combined report
     */
    public static String compareModes(AbstractCCD ccd, List<Tree> heldOut) {
        NNICladeExpansion broad = new NNICladeExpansion(ccd, PairingMode.GRAPH_WIDE).compute();
        StringBuilder sb = new StringBuilder();
        sb.append(broad.report());

        NNICladeExpansion strict = null;
        if (ccd.getBaseTrees() == null) {
            sb.append("NNI clade expansion [CO_OCCURRING]: skipped (base trees not stored)\n");
        } else {
            strict = new NNICladeExpansion(ccd, PairingMode.CO_OCCURRING).compute();
            sb.append(strict.report());
        }

        if (heldOut != null && !heldOut.isEmpty()) {
            sb.append(heldOutClassification(broad, strict, heldOut));
        }

        if (strict == null) {
            return sb.toString();
        }

        Set<BitSet> broadSet = broad.novelClades.keySet();
        Set<BitSet> strictSet = strict.novelClades.keySet();
        int onlyBroad = 0;
        for (BitSet b : broadSet) {
            if (!strictSet.contains(b)) {
                onlyBroad++;
            }
        }
        int onlyStrict = 0;
        for (BitSet b : strictSet) {
            if (!broadSet.contains(b)) {
                onlyStrict++;
            }
        }
        sb.append("comparison\n");
        sb.append(String.format("  shared novel clades  : %d%n", strictSet.size() - onlyStrict));
        sb.append(String.format("  graph-wide only      : %d%n", onlyBroad));
        sb.append(String.format("  co-occurring only    : %d%n", onlyStrict));
        return sb.toString();
    }

    /* Format the held-out classification block for one or both pairing modes. */
    private static String heldOutClassification(NNICladeExpansion broad,
                                                NNICladeExpansion strict,
                                                List<Tree> heldOut) {
        StringBuilder sb = new StringBuilder();
        sb.append(String.format("held-out classification (%d trees)%n", heldOut.size()));
        appendZoneLine(sb, "GRAPH_WIDE", broad.classifyTrees(heldOut));
        if (strict != null) {
            appendZoneLine(sb, "CO_OCCURRING", strict.classifyTrees(heldOut));
        }
        return sb.toString();
    }

    private static void appendZoneLine(StringBuilder sb, String label, ZoneCounts c) {
        int total = c.total();
        double pct = total == 0 ? 0 : 100.0 / total;
        sb.append(String.format("  [%s] observed: %d (%.2f%%)  NNI-expanded: %d (%.2f%%)  out-of-universe: %d (%.2f%%)%n",
                label,
                c.observed(), c.observed() * pct,
                c.nniExpanded(), c.nniExpanded() * pct,
                c.outOfUniverse(), c.outOfUniverse() * pct));
    }

    /**
     * CLI: report NNI clade-expansion properties of a tree file, optionally with
     * a held-out tree file classified by zone (observed / NNI-expanded /
     * out-of-universe).
     * <p>
     * Usage: {@code NNICladeExpansion <trees-file> [burnin-percentage] [held-out-trees-file]}
     *
     * @param args trees file path, optional burn-in percentage (default 10),
     *             and optional held-out trees file
     * @throws Exception if a tree file cannot be read
     */
    public static void main(String[] args) throws Exception {
        if (args.length < 1) {
            System.err.println("Usage: NNICladeExpansion <trees-file> [burnin-percentage] [held-out-trees-file]");
            System.exit(1);
        }
        int burnin = (args.length >= 2) ? Integer.parseInt(args[1]) : 10;
        TreeAnnotator.MemoryFriendlyTreeSet treeSet = CCDToolUtil.getTreeSet(args[0], burnin);
        CCD0 ccd = new CCD0(treeSet, true);

        List<Tree> heldOut = null;
        if (args.length >= 3) {
            TreeAnnotator.MemoryFriendlyTreeSet heldOutSet = CCDToolUtil.getTreeSet(args[2], burnin);
            heldOut = new ArrayList<>();
            heldOutSet.reset();
            while (heldOutSet.hasNext()) {
                heldOut.add(heldOutSet.next());
            }
        }

        System.out.print(compareModes(ccd, heldOut));
    }
}
