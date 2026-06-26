package ccd.model;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import ccd.model.bitsets.BitSet;

import java.util.ArrayList;

/**
 * Shared logic for the extended-split (SJ) sampled-ancestor clade identity,
 * factored out so {@link CCD1SJ} and {@link CCD0SJ} share one implementation
 * (Java's single-inheritance forces them to extend different parents).
 *
 * <p>The SJ identity (per {@code alexei_c_four_sa_models.pdf}, Table 1):
 * each internal clade carries SA flags for its direct leaf children only,
 * with no propagation up the tree. Identity bitset width is {@code 2n+1}:
 * <ul>
 *   <li>bits {@code [0, n)}: taxa membership
 *   <li>bits {@code [n, 2n)}: bit {@code n+i} set iff leaf {@code i} is a
 *       direct SA leaf child of this clade
 *   <li>bit {@code 2n}: sentinel for the placeholder root clade
 * </ul>
 *
 * <p>The placeholder {@code rootClade} carries one
 * {@code (variant, emptyClade)} partition per observed root variant; the
 * shared {@code emptyClade} is a neutral sister with no taxa.
 *
 * <p>All methods are static and package-private. State is passed in
 * explicitly via the {@link AbstractCCD} parameter (whose protected fields
 * are accessible from this same package) and the {@code emptyClade}
 * argument.
 */
final class SJSupport {

    private SJSupport() {
    }

    /**
     * Bitset width used for all clades under the SJ identity: {@code 2n + 1}
     * (taxa bits, SA-flag bits, sentinel bit).
     */
    static int extendedBitWidth(int leafArraySize) {
        return 2 * leafArraySize + 1;
    }

    /* -- KEY CONSTRUCTION -- */

    /**
     * Bitset key for a leaf clade: only the taxon bit is set (SA flags live
     * on the parent under the SJ identity, not on the leaf itself).
     */
    static BitSet leafKey(int taxonIndex, int leafArraySize) {
        BitSet key = BitSet.newBitSet(extendedBitWidth(leafArraySize));
        key.set(taxonIndex);
        return key;
    }

    /**
     * Bitset key for an internal clade under the SJ identity: union of
     * descendant taxa bits, plus an SA-flag bit for each direct leaf child of
     * this vertex with branch length zero. SA flags do not propagate up.
     */
    static BitSet extendedSelfKey(Node vertex, int leafArraySize) {
        BitSet key = BitSet.newBitSet(extendedBitWidth(leafArraySize));
        collectTaxaBits(vertex, key);
        for (Node child : vertex.getChildren()) {
            if (child.isLeaf() && child.getLength() == 0) {
                key.set(leafArraySize + child.getNr());
            }
        }
        return key;
    }

    private static void collectTaxaBits(Node vertex, BitSet key) {
        if (vertex.isLeaf()) {
            key.set(vertex.getNr());
        } else {
            for (Node child : vertex.getChildren()) {
                collectTaxaBits(child, key);
            }
        }
    }

    /* -- ROOT INITIALISATION -- */

    /**
     * Sets up the placeholder root clade (sentinel bit only) and the shared
     * empty sister, registers both in {@code ccd.cladeMapping}, and returns
     * the empty sister so the caller can store the field reference. Must be
     * invoked from the subclass's {@code initializeRootClade} override.
     */
    static Clade initializeSJRoot(AbstractCCD ccd, int numLeaves) {
        ccd.leafArraySize = numLeaves;
        int width = extendedBitWidth(numLeaves);

        BitSet rootBitSet = BitSet.newBitSet(width);
        rootBitSet.set(2 * numLeaves);
        ccd.rootClade = new Clade(rootBitSet, ccd);
        ccd.cladeMapping.put(ccd.rootClade.getCladeInBits(), ccd.rootClade);

        BitSet emptyKey = BitSet.newBitSet(width);
        return ccd.addNewClade(emptyKey);
    }

    /* -- TREE INSERTION -- */

    /**
     * Top-level cladify entry: walks the tree under the SJ identity and
     * attaches the resulting root variant to {@code ccd.rootClade} via a
     * {@code (variant, emptyClade)} partition. Increments the placeholder
     * root's occurrence count once per tree.
     */
    static Clade cladifyTreeRoot(AbstractCCD ccd, Clade emptyClade, Node root) {
        return cladifySJVertex(ccd, emptyClade, root, true);
    }

    private static Clade cladifySJVertex(AbstractCCD ccd, Clade emptyClade, Node vertex, boolean isRoot) {
        if (vertex.isLeaf()) {
            BitSet leafKey = leafKey(vertex.getNr(), ccd.leafArraySize);
            Clade leafClade = ccd.cladeMapping.get(leafKey);
            if (leafClade == null) {
                leafClade = ccd.addNewClade(leafKey);
            }
            leafClade.increaseOccurrenceCount(vertex.getHeight());
            return leafClade;
        }

        Clade firstChildClade = cladifySJVertex(ccd, emptyClade, vertex.getChildren().get(0), false);
        Clade secondChildClade = cladifySJVertex(ccd, emptyClade, vertex.getChildren().get(1), false);

        BitSet selfKey = extendedSelfKey(vertex, ccd.leafArraySize);
        Clade currentClade = ccd.cladeMapping.get(selfKey);
        if (currentClade == null) {
            currentClade = ccd.addNewClade(selfKey);
        }
        currentClade.increaseOccurrenceCount(vertex.getHeight());

        CladePartition partition = currentClade.getCladePartition(firstChildClade, secondChildClade);
        if (partition == null) {
            partition = currentClade.createCladePartition(firstChildClade, secondChildClade);
        }
        partition.increaseOccurrenceCount(vertex.getHeight());

        if (isRoot) {
            // Attach this root variant under the placeholder rootClade via a
            // (variant, emptyClade) partition. rootClade's own occurrence
            // count ticks up once per tree.
            ccd.rootClade.increaseOccurrenceCount(vertex.getHeight());
            CladePartition rootPart = ccd.rootClade.getCladePartition(currentClade, emptyClade);
            if (rootPart == null) {
                rootPart = ccd.rootClade.createCladePartition(currentClade, emptyClade);
            }
            rootPart.increaseOccurrenceCount(vertex.getHeight());
        }

        return currentClade;
    }

    /* -- SA DETECTION & DISPLAY -- */

    /**
     * True iff the clade has at least one direct sampled-ancestor leaf
     * child (any bit set in the SA region {@code [n, 2n)}).
     */
    static boolean isSampledAncestor(Clade clade, int leafArraySize) {
        BitSet bits = clade.getCladeInBits();
        for (int i = leafArraySize; i < 2 * leafArraySize; i++) {
            if (bits.get(i)) return true;
        }
        return false;
    }

    /**
     * Renders the clade's SA-flag region as a human-readable string for
     * inclusion in {@link Clade#toString()}.
     */
    static String saInfoString(Clade clade, int leafArraySize) {
        BitSet bits = clade.getCladeInBits();
        java.util.List<Integer> sa = new ArrayList<>();
        for (int i = leafArraySize; i < 2 * leafArraySize; i++) {
            if (bits.get(i)) sa.add(i - leafArraySize);
        }
        if (sa.isEmpty()) return "sampled ancestor = none";
        return "sampled ancestor taxa = " + sa;
    }

    /* -- SJ-FLAG COMPATIBILITY (used by CCD0SJ.expand) -- */

    /**
     * True iff the partition {@code (A, B)} is consistent with
     * {@code parent}'s SA-flag identity: every SA-flagged leaf {@code l} of
     * {@code parent} must appear as the singleton clade {@code {l}} in
     * {@code A} or {@code B}. Used to filter expanded partitions so the
     * model never assigns probability to trees that violate the identity.
     */
    static boolean isSJCompatible(Clade parent, Clade A, Clade B, int leafArraySize) {
        BitSet bits = parent.getCladeInBits();
        for (int i = leafArraySize; i < 2 * leafArraySize; i++) {
            if (!bits.get(i)) continue;
            int taxon = i - leafArraySize;
            boolean aIsSingletonForTaxon = A.size() == 1 && A.getCladeInBitsTaxaOnly().get(taxon);
            boolean bIsSingletonForTaxon = B.size() == 1 && B.getCladeInBitsTaxaOnly().get(taxon);
            if (!aIsSingletonForTaxon && !bIsSingletonForTaxon) return false;
        }
        return true;
    }

    /* -- TREE PROBABILITY -- */

    /**
     * Computes {@code P(T)} (or its log) by walking {@code root} under the
     * SJ identity and folding in CCPs along the path, including the
     * placeholder root's {@code (variant, emptyClade)} factor. Returns the
     * variant clade for the tree's root, or {@code null} if any vertex's
     * extended clade or partition isn't in the model.
     */
    static Clade computeTreeProbability(AbstractCCD ccd, Clade emptyClade, Node root,
                                        double[] runningProbability, boolean computeLog) {
        Clade variant = computeProbSJVertex(ccd, root, runningProbability, computeLog);
        if (variant == null) return null;
        CladePartition rootPart = ccd.rootClade.getCladePartition(variant, emptyClade);
        if (rootPart == null) {
            AbstractCCD.setComputedNoProbability(runningProbability, computeLog);
            return null;
        }
        if (computeLog) {
            runningProbability[0] += rootPart.getLogCCP();
        } else {
            runningProbability[0] *= rootPart.getCCP();
        }
        return variant;
    }

    private static Clade computeProbSJVertex(AbstractCCD ccd, Node vertex,
                                             double[] runningProbability, boolean computeLog) {
        if (vertex.isLeaf()) {
            BitSet leafKey = leafKey(vertex.getNr(), ccd.leafArraySize);
            Clade leafClade = ccd.cladeMapping.get(leafKey);
            if (leafClade == null) {
                AbstractCCD.setComputedNoProbability(runningProbability, computeLog);
                return null;
            }
            return leafClade;
        }

        Clade firstChildClade = computeProbSJVertex(ccd, vertex.getChildren().get(0), runningProbability, computeLog);
        Clade secondChildClade = computeProbSJVertex(ccd, vertex.getChildren().get(1), runningProbability, computeLog);

        if (firstChildClade == null || secondChildClade == null) {
            AbstractCCD.setComputedNoProbability(runningProbability, computeLog);
            return null;
        }

        BitSet selfKey = extendedSelfKey(vertex, ccd.leafArraySize);
        Clade currentClade = ccd.cladeMapping.get(selfKey);
        if (currentClade == null) {
            AbstractCCD.setComputedNoProbability(runningProbability, computeLog);
            return null;
        }

        CladePartition partition = currentClade.getCladePartition(firstChildClade, secondChildClade);
        if (partition == null) {
            AbstractCCD.setComputedNoProbability(runningProbability, computeLog);
            return null;
        }

        if (computeLog) {
            runningProbability[0] += partition.getLogCCP();
        } else {
            runningProbability[0] *= partition.getCCP();
        }
        return currentClade;
    }

    /* -- HELD-OUT (LEAVE-ONE-TREE-OUT) TREE PROBABILITY -- */

    /**
     * Leave-one-tree-out counterpart of {@link #computeTreeProbability}: folds
     * in the held-out CCPs {@code (count-1+alpha)/(parentCount+alpha*k-1)}
     * instead of the in-sample CCPs, so the returned value is
     * {@code P(T | CCD trained on all trees except T)}. Mirrors the EL/CP
     * held-out traversal in {@link CCD1#getLogProbOfHeldOutTree}: a clade or
     * partition seen only in the held-out tree contributes no probability.
     */
    static Clade computeTreeProbabilityHeldOut(AbstractCCD ccd, Clade emptyClade, Node root, double[] runningProbability, double alpha, boolean computeLog) {
        Clade variant = computeProbSJVertexHeldOut(ccd, root, runningProbability, alpha, computeLog);
        if (variant == null) return null;
        CladePartition rootPart = ccd.rootClade.getCladePartition(variant, emptyClade);
        if (rootPart == null) {
            AbstractCCD.setComputedNoProbability(runningProbability, computeLog);
            return null;
        }
        if (computeLog) {
            runningProbability[0] += rootPart.getLogCCPInHeldOutTree(alpha);
        } else {
            runningProbability[0] *= rootPart.getCCPInHeldOutTree(alpha);
        }
        return variant;
    }

    private static Clade computeProbSJVertexHeldOut(AbstractCCD ccd, Node vertex, double[] runningProbability, double alpha, boolean computeLog) {
        if (vertex.isLeaf()) {
            BitSet leafKey = leafKey(vertex.getNr(), ccd.leafArraySize);
            Clade leafClade = ccd.cladeMapping.get(leafKey);
            if (leafClade == null) {
                AbstractCCD.setComputedNoProbability(runningProbability, computeLog);
                return null;
            }
            return leafClade;
        }

        Clade firstChildClade = computeProbSJVertexHeldOut(ccd, vertex.getChildren().get(0), runningProbability, alpha, computeLog);
        Clade secondChildClade = computeProbSJVertexHeldOut(ccd, vertex.getChildren().get(1), runningProbability, alpha, computeLog);

        if (firstChildClade == null || secondChildClade == null) {
            AbstractCCD.setComputedNoProbability(runningProbability, computeLog);
            return null;
        }

        BitSet selfKey = extendedSelfKey(vertex, ccd.leafArraySize);
        Clade currentClade = ccd.cladeMapping.get(selfKey);
        // A clade seen in exactly one tree exists only because of the held-out
        // tree; the held-out model has never seen it, so contribute no probability
        // (also avoids a zero denominator in the held-out CCP at alpha=0).
        if (currentClade == null || currentClade.getNumberOfOccurrences() == 1) {
            AbstractCCD.setComputedNoProbability(runningProbability, computeLog);
            return null;
        }

        CladePartition partition = currentClade.getCladePartition(firstChildClade, secondChildClade);
        if (partition == null) {
            AbstractCCD.setComputedNoProbability(runningProbability, computeLog);
            return null;
        }

        if (computeLog) {
            runningProbability[0] += partition.getLogCCPInHeldOutTree(alpha);
        } else {
            runningProbability[0] *= partition.getCCPInHeldOutTree(alpha);
        }
        return currentClade;
    }

    /* -- SAMPLING & MAP -- */

    /**
     * Builds a tree under the SJ identity: pick a root variant via the
     * placeholder root's {@code (variant, emptyClade)} partitions, then
     * recursively pick partitions on the variant's subtree.
     *
     * <p><b>TODO — height strategy not yet supported.</b> Heights are always
     * assigned with the "One" semantics (parent = max(child) + 1, SA leaves
     * get parent height for branch length zero), regardless of the
     * {@code heightStrategy} argument. Implementing
     * {@code MeanOccurredHeights} / {@code CommonAncestorHeights} correctly
     * under SJ requires per-SA-context height stats: a leaf can appear as
     * an SA child in some trees and a non-SA child in others, and
     * {@link Clade#getMeanOccurredHeight()} averages over both, which is
     * the wrong number whichever way you use it.
     */
    static Tree buildTree(AbstractCCD ccd, Clade emptyClade,
                          SamplingStrategy samplingStrategy, HeightSettingStrategy heightStrategy) {
        ccd.tidyUpCacheIfDirty();

        CladePartition rootPart = pickRootPartition(ccd, emptyClade, samplingStrategy);
        Clade variant = variantChildOf(rootPart, emptyClade);

        int[] innerIdx = new int[]{ccd.getSizeOfLeavesArray()};
        Node root = buildSubtree(ccd, variant, samplingStrategy, innerIdx);
        return new Tree(root);
    }

    private static Clade variantChildOf(CladePartition partition, Clade emptyClade) {
        Clade[] kids = partition.getChildClades();
        return (kids[0] == emptyClade) ? kids[1] : kids[0];
    }

    private static CladePartition pickRootPartition(AbstractCCD ccd, Clade emptyClade, SamplingStrategy strategy) {
        ArrayList<CladePartition> partitions = ccd.rootClade.getPartitions();
        if (partitions.isEmpty()) {
            throw new AssertionError("SJ rootClade has no variant partitions.");
        }
        if (strategy == SamplingStrategy.MAP) {
            double bestScore = Double.NEGATIVE_INFINITY;
            CladePartition best = partitions.get(0);
            for (CladePartition p : partitions) {
                Clade variant = variantChildOf(p, emptyClade);
                double score = Math.log(p.getCCP()) + variant.getMaxSubtreeLogCCP();
                if (score > bestScore) {
                    bestScore = score;
                    best = p;
                }
            }
            return best;
        }
        return weightedDraw(partitions, ccd);
    }

    private static Node buildSubtree(AbstractCCD ccd, Clade clade, SamplingStrategy strategy, int[] innerIdx) {
        if (clade.isLeaf()) {
            int leafNr = clade.getCladeInBits().nextSetBit(0);
            String taxonName = ccd.getSomeBaseTree().getTaxaNames()[leafNr];
            Node node = new Node(taxonName);
            node.setNr(leafNr);
            return node;
        }

        CladePartition partition = pickInternalPartition(ccd, clade, strategy);
        if (partition == null) {
            throw new AssertionError("No partitions on clade: " + clade.getCladeInBits());
        }

        Node firstChild = buildSubtree(ccd, partition.getChildClades()[0], strategy, innerIdx);
        Node secondChild = buildSubtree(ccd, partition.getChildClades()[1], strategy, innerIdx);

        Node node = new Node();
        node.setNr(innerIdx[0]++);
        node.addChild(firstChild);
        node.addChild(secondChild);

        // "One" height semantics: parent = max(children) + 1; SA leaves get
        // parent's height (branch length 0).
        double parentHeight = Math.max(safeHeight(firstChild), safeHeight(secondChild)) + 1.0;
        node.setHeight(parentHeight);

        for (Node child : node.getChildren()) {
            if (child.isLeaf()) {
                int leafNr = child.getNr();
                boolean isSA = clade.getCladeInBits().get(ccd.leafArraySize + leafNr);
                child.setHeight(isSA ? parentHeight : 0.0);
            }
        }
        return node;
    }

    private static CladePartition pickInternalPartition(AbstractCCD ccd, Clade clade, SamplingStrategy strategy) {
        if (strategy == SamplingStrategy.MAP) {
            return clade.getMaxSubtreeCCPPartition();
        }
        return weightedDraw(clade.getPartitions(), ccd);
    }

    private static CladePartition weightedDraw(ArrayList<CladePartition> partitions, AbstractCCD ccd) {
        double total = 0;
        for (CladePartition p : partitions) total += p.getCCP();
        double r = ccd.random.nextDouble() * total;
        double cum = 0;
        for (CladePartition p : partitions) {
            cum += p.getCCP();
            if (r < cum) return p;
        }
        return partitions.get(partitions.size() - 1);
    }

    private static double safeHeight(Node n) {
        return n.isLeaf() ? 0.0 : n.getHeight();
    }
}
