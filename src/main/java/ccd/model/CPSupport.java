package ccd.model;

import beast.base.evolution.tree.Node;
import ccd.model.bitsets.BitSet;

import java.util.Map;

final class CPSupport {
    private CPSupport() {
    }

    static int extendedBitWidth(int leafArraySize) {
        return leafArraySize + 1;
    }

    static BitSet leafKey(Node vertex, int leafArraySize) {
        BitSet key = BitSet.newBitSet(extendedBitWidth(leafArraySize));
        key.set(vertex.getNr());
        if (vertex.getLength() != 0) {
            key.set(leafArraySize);
        }
        return key;
    }

    static BitSet extendedSelfKey(Clade firstChildClade, Clade secondChildClade, int leafArraySize) {
        BitSet key = BitSet.newBitSet(extendedBitWidth(leafArraySize));
        key.or(firstChildClade.getCladeInBits());
        key.or(secondChildClade.getCladeInBits());
        return key;
    }

    static void initializeRoot(AbstractCCD ccd, int numLeaves) {
        ccd.leafArraySize = numLeaves;
        int width = extendedBitWidth(numLeaves);
        BitSet rootBitSet = BitSet.newBitSet(width);
        rootBitSet.set(0, width);
        ccd.rootClade = new Clade(rootBitSet, ccd);
        ccd.cladeMapping.put(ccd.rootClade.getCladeInBits(), ccd.rootClade);
    }

    static Clade cladifyVertex(AbstractCCD ccd, Node vertex) {
        if (vertex.isLeaf()) {
            BitSet leafKey = leafKey(vertex, ccd.leafArraySize);
            Clade leafClade = ccd.cladeMapping.get(leafKey);
            if (leafClade == null) {
                leafClade = ccd.addNewClade(leafKey);
            }
            leafClade.increaseOccurrenceCount(vertex.getHeight());
            return leafClade;
        }

        Clade firstChildClade = cladifyVertex(ccd, vertex.getChildren().get(0));
        Clade secondChildClade = cladifyVertex(ccd, vertex.getChildren().get(1));

        BitSet selfKey = extendedSelfKey(firstChildClade, secondChildClade, ccd.leafArraySize);
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

        return currentClade;
    }

    static boolean isSampledAncestor(Clade clade, int leafArraySize) {
        // the last bit, i.e. at index = number of taxa is SA flag
        // Note: we use 0 to represent SA, and 1 to represent non-SA
        return !clade.getCladeInBits().get(leafArraySize);
    }

    static double computeParentHeight(CladePartition partition, Node firstChild, Node secondChild, int leafArraySize) {
        if (isSampledAncestorPartition(partition, leafArraySize)) {
            return Math.max(firstChild.getHeight(), secondChild.getHeight());
        } else {
            return Math.max(firstChild.getHeight(), secondChild.getHeight()) + 1.0;
        }
    }

    static boolean isSampledAncestorPartition(CladePartition partition, int leafArraySize) {
        return isSampledAncestor(partition.getChildClades()[0], leafArraySize) ||
                isSampledAncestor(partition.getChildClades()[1], leafArraySize);
    }

    static Clade reduceCladeCount(AbstractCCD ccd, Node vertex, int leafArraySize) {
        BitSet cladeInBits = BitSet.newBitSet(extendedBitWidth(leafArraySize));

        if (vertex.isLeaf()) {
            BitSet leafKey = leafKey(vertex, leafArraySize);
            Clade leafClade = ccd.cladeMapping.get(leafKey);
            leafClade.decreaseOccurrenceCount(vertex.getHeight());
            return leafClade;
        }

        Clade firstChildClade = reduceCladeCount(ccd, vertex.getChildren().get(0), leafArraySize);
        Clade secondChildClade = reduceCladeCount(ccd, vertex.getChildren().get(1), leafArraySize);

        BitSet selfKey = extendedSelfKey(firstChildClade, secondChildClade, leafArraySize);

        // retrieve clade and reduce count
        Clade currentClade = ccd.cladeMapping.get(selfKey);
        currentClade.decreaseOccurrenceCount(vertex.getHeight());

        // reduce counts for its clade partitions
        CladePartition currentPartition = currentClade.getCladePartition(firstChildClade, secondChildClade);
        currentPartition.decreaseOccurrenceCount(vertex.getHeight());

        ccd.removeCladePartitionIfNecessary(currentClade, currentPartition);

        return currentClade;
    }

    static BitSet computeCladeToNodeMapping(Node vertex, Map<Clade, Node> map, AbstractCCD ccd, int leafArraySize) {
        BitSet bits;
        if (vertex.isLeaf()) {
            bits = leafKey(vertex, leafArraySize);
        } else {
            bits = computeCladeToNodeMapping(vertex.getLeft(), map, ccd, leafArraySize);
            BitSet otherBits = computeCladeToNodeMapping(vertex.getRight(), map, ccd, leafArraySize);
            bits.or(otherBits);
        }

        Clade clade = ccd.getClade(bits);
        if (clade == null) {
            System.err.println("No clade found in CCD for this vertex (" + bits + ").");
        }
        map.put(clade, vertex);

        return bits;
    }

    static Clade computeProbCPVertex(Node vertex, double[] runningProbability, boolean computeLog, AbstractCCD ccd, int leafArraySize) {
        if (vertex.isLeaf()) {
            // set bitSet
            BitSet leafKey = leafKey(vertex, leafArraySize);
            Clade leafClade = ccd.cladeMapping.get(leafKey);

            // probability of the leaf
            if (computeLog) {
                runningProbability[0] += leafClade.getLogProbability();
            } else {
                runningProbability[0] *= leafClade.getProbability();
            }
            return leafClade;
        }

        Clade firstChildClade = computeProbCPVertex(vertex.getChildren().get(0), runningProbability, computeLog, ccd, leafArraySize);
        Clade secondChildClade = computeProbCPVertex(vertex.getChildren().get(1), runningProbability, computeLog, ccd, leafArraySize);

        if (computeLog && runningProbability[0] > 0) {
            return null;
        }

        if ((firstChildClade == null) || (secondChildClade == null)) {
            ccd.setComputedNoProbability(runningProbability, computeLog);
            return null;
        }

        BitSet selfKey = extendedSelfKey(firstChildClade, secondChildClade, leafArraySize);
        Clade currentClade = ccd.cladeMapping.get(selfKey);
        if (currentClade != null) {
            CladePartition partition = currentClade.getCladePartition(firstChildClade, secondChildClade);
            if (partition != null) {
                if (computeLog) {
                    runningProbability[0] += partition.getLogCCP();
                } else {
                    runningProbability[0] *= partition.getCCP();
                }
            } else {
                ccd.setComputedNoProbability(runningProbability, computeLog);
            }
        } else {
            ccd.setComputedNoProbability(runningProbability, computeLog);
        }
        return currentClade;
    }

    /* -- HELD-OUT (LEAVE-ONE-TREE-OUT) TREE PROBABILITY -- */

    static Clade computeTreeProbabilityHeldOut(AbstractCCD ccd, Node root, double[] runningProbability, double alpha, boolean computeLog) {
        Clade variant = computeProbCPVertexHeldOut(ccd, root, runningProbability, alpha, computeLog);
        return variant;
    }

    private static Clade computeProbCPVertexHeldOut(AbstractCCD ccd, Node vertex, double[] runningProbability, double alpha, boolean computeLog) {
        if (vertex.isLeaf()) {
            BitSet leafKey = leafKey(vertex, ccd.leafArraySize);
            Clade leafClade = ccd.cladeMapping.get(leafKey);
            if (leafClade == null) {
                AbstractCCD.setComputedNoProbability(runningProbability, computeLog);
                return null;
            }
            return leafClade;
        }

        Clade firstChildClade = computeProbCPVertexHeldOut(ccd, vertex.getChildren().get(0), runningProbability, alpha, computeLog);
        Clade secondChildClade = computeProbCPVertexHeldOut(ccd, vertex.getChildren().get(1), runningProbability, alpha, computeLog);

        if (firstChildClade == null || secondChildClade == null) {
            AbstractCCD.setComputedNoProbability(runningProbability, computeLog);
            return null;
        }

        BitSet selfKey = extendedSelfKey(firstChildClade, secondChildClade, ccd.leafArraySize);
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

    static String saInfoString(Clade clade, int leafArraySize) {
        BitSet bits = clade.getCladeInBits();
        boolean isNonSA = bits.get(leafArraySize);
        return isNonSA ? "sampled ancestor = false" : "sampled ancestor = true";
    }

}

