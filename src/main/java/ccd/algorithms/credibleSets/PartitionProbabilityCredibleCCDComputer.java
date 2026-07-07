package ccd.algorithms.credibleSets;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import ccd.model.AbstractCCD;
import ccd.model.CCD2;
import ccd.model.Clade;
import ccd.model.CladePartition;
import ccd.model.ExtendedClade;
import ccd.model.bitsets.BitSet;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

/**
 * A credible CCD based on clade partitions, intended for CCD1 and CCD2 models.
 * The credible set information is constructed by continuously removing the clade partitions with the lowest probability.
 *
 * @author Jonathan Klawitter
 */
@SuppressWarnings("unchecked")
public class PartitionProbabilityCredibleCCDComputer extends CredibleCCDComputer {

    /**
     * Default constructor not to be used directly.
     *
     * @param baseCCD the CCD this credible CCD is based on
     */
    protected PartitionProbabilityCredibleCCDComputer(AbstractCCD baseCCD) {
        super(baseCCD);
    }

    @Override
    public double getCredibleLevel(Tree tree) {
        double[] credLevel = new double[]{Double.NEGATIVE_INFINITY};
        if (fullCCD instanceof CCD2) {
            Node root = tree.getRoot();
            getCredibleLevelOfPartitionCCD2(root.getChild(0), root.getChild(1), credLevel);
        } else {
            getCredibleLevelOfPartitionGeneral(tree.getRoot(), credLevel);
        }
        if (credLevel[0] <= 0) {
            return -1;
        }
        return credLevel[0];
    }

    /* Recursive helper method */
    private BitSet getCredibleLevelOfPartitionGeneral(Node node, double[] credLevel) {
        BitSet bits;
        if (node.isLeaf()) {
            bits = BitSet.newBitSet(fullCCD.getNumberOfLeaves());
            bits.set(node.getNr());
        } else {
            bits = getCredibleLevelOfPartitionGeneral(node.getLeft(), credLevel);
            if (credLevel[0] == 0) {
                return bits;
            }
            BitSet other = getCredibleLevelOfPartitionGeneral(node.getRight(), credLevel);
            if (credLevel[0] == 0) {
                return bits;
            }

            Clade leftChild = fullCCD.getClade(bits);
            Clade rightChild = fullCCD.getClade(other);

            bits.or(other);
            Clade parent = fullCCD.getClade(bits);
            if (parent == null) {
                credLevel[0] = 0;
                // out.println("clade not found, bits = " + bits);
                // out.println("clades in ccd:");
                // for (Clade clade : this.getClades()) {
                //     out.println("clade = " + clade);
                // }
                return bits;
            }

            CladePartition partition = parent.getCladePartition(leftChild, rightChild);
            if (partition == null) {
                credLevel[0] = 0;
            } else {
                credLevel[0] = Math.max(credLevel[0], partitionMinCredibility.get(partition));
            }
        }

        return bits;
    }

    /* Recursive helper method */
    private ExtendedClade[] getCredibleLevelOfPartitionCCD2(Node leftVertex, Node rightVertex, double[] credLevel) {
        BitSet leftInBits = BitSet.newBitSet(fullCCD.getSizeOfLeavesArray());
        BitSet rightInBits = BitSet.newBitSet(fullCCD.getSizeOfLeavesArray());

        ExtendedClade[] leftChildren = getCredibleLevelOfPartitionChildrenCCD2(leftVertex, leftInBits, credLevel);
        if (!leftVertex.isLeaf() && (leftChildren == null)) {
            return null;
        }
        ExtendedClade[] rightChildren = getCredibleLevelOfPartitionChildrenCCD2(rightVertex, rightInBits, credLevel);
        if (!rightVertex.isLeaf() && (rightChildren == null)) {
            return null;
        }

        ExtendedClade leftClade = ((CCD2) fullCCD).getExtendedClade(leftInBits, rightInBits);
        ExtendedClade rightClade = ((CCD2) fullCCD).getExtendedClade(rightInBits, leftInBits);
        if ((leftClade == null) || (rightClade == null)) {
            credLevel[0] = 0;
            return null;
        }

        if (!leftVertex.isLeaf()) {
            CladePartition partition = leftClade.getCladePartition(leftChildren[0], leftChildren[1]);
            if (partition != null) {
                credLevel[0] = Math.max(credLevel[0], partitionMinCredibility.get(partition));
            } else {
                credLevel[0] = 0;
                return null;
            }
        }
        if (!rightVertex.isLeaf()) {
            CladePartition partition = rightClade.getCladePartition(rightChildren[0], rightChildren[1]);
            if (partition != null) {
                credLevel[0] = Math.max(credLevel[0], partitionMinCredibility.get(partition));
            } else {
                credLevel[0] = 0;
                return null;
            }
        }

        return new ExtendedClade[]{leftClade, rightClade};
    }

    /* Recursive helper method */
    private ExtendedClade[] getCredibleLevelOfPartitionChildrenCCD2(Node vertex, BitSet cladeInBits, double[] credLevel) {
        ExtendedClade[] leftChildren = null;
        if (vertex.isLeaf()) {
            int index = vertex.getNr();
            cladeInBits.set(index);
        } else {
            leftChildren = getCredibleLevelOfPartitionCCD2(vertex.getChildren().get(0), vertex.getChildren().get(1), credLevel);

            if (leftChildren == null) {
                return null;
            }

            cladeInBits.or(leftChildren[0].getCladeInBits());
            cladeInBits.or(leftChildren[1].getCladeInBits());
        }

        return leftChildren;
    }

    @Override
    public AbstractCCD getCredibleCCD(double alpha) {
        // TODO
        return null;
    }

    @Override
    protected void computeCredibleSetInformationHelper() {
        double numPartitionsToHandle = partialCCD.getNumberOfCladePartitions() - partialCCD.getNumberOfLeaves() + 1;
        // System.out.println("partialCCD.numPartitions = " + partialCCD.getNumberOfCladePartitions());
        // System.out.println("numPartitionsToHandle    = " + numPartitionsToHandle);

        List<Clade> cladesToHandle = new ArrayList<>(partialCCD.getNumberOfClades() - partialCCD.getNumberOfLeaves());
        for (Clade clade : partialCCD.getClades()) {
            if (!clade.isLeaf()) {
                cladesToHandle.add(clade);
            }
        }
        cladesToHandle.sort(Comparator.comparingInt(Clade::size));

        do {
            CladePartition nextPartition = getNextCladePartition();

            // remove clade and tidy up CCD graph and store results
            List[] removedObjects = reduceCCD(nextPartition);

            // then for the removed clades and partitions, set the resulting credible levels
            setCredibleLevels(removedObjects);

            // double remainingNontrivialPartitions = numPartitionsToHandle - (partialCCD.getNumberOfCladePartitions() - minNumClades + 1);
            // while ((remainingNontrivialPartitions * 100 / numPartitionsToHandle) > donePercentage) {
            //     donePercentage++;
            // System.out.print(".");
            // }

            // recompute CCPs and clade probabilities in partial CCD
            if (removedObjects[0] != null) {
                cladesToHandle.removeAll((List<Clade>) removedObjects[0]);
            }
            double reductionFactor = recomputeCCPs(cladesToHandle);
            partialCCD.computeCladeProbabilities();

            // update all values
            remainingProbability *= reductionFactor;
            partialCCD.setCacheAsDirty();
            partialCCD.computeCladeProbabilities();
            // System.out.println("reductionFactor = " + reductionFactor);
            // System.out.println("handled partitions: " + partitionMinCredibility.size());
            // System.out.println("remaining partitions: " + partialCCD.getNumberOfCladePartitions());
            // System.out.println("remaining clades: " + partialCCD.getNumberOfClades());

            // prevMAPTree = currentMAPTree;
            // currentMAPTree = new WrappedBeastTree(partialCCD.getMAPTree());
            // writeCurrentResult(remainingProbability);

            // we exist when we have reduced the CCD to a single tree
        } while (partialCCD.getNumberOfClades() != minNumClades);

        // finally, for the clades and clade partitions remaining, also set the credible levels
        for (Clade clade : partialCCD.getClades()) {
            setCredibleLevel(clade, remainingProbability);
            for (CladePartition partition : clade.getPartitions()) {
                setCredibleLevel(partition, remainingProbability);
            }
        }
    }

    /**
     * Find the next clade partition to remove based on its probability as per this strategy.
     *
     * @return the next clade partition to remove
     */
    private CladePartition getNextCladePartition() {
        CladePartition nextPartition = null;
        double min = Double.POSITIVE_INFINITY;

        // pick clade partition with the lowest probability
        for (Clade clade : partialCCD.getClades()) {
            if (clade.isLeaf()) {
                continue;
            }
            for (CladePartition partition : clade.getPartitions()) {
                double p = clade.getProbability() * partition.getCCP();

                // for ties in probability, we pick the clade partition with the larger parent clade
                if ((p < min) || ((p == min) && (clade.size() > nextPartition.getParentClade().size()))) {
                    min = p;
                    nextPartition = partition;
                }
            }
        }

        if (nextPartition == null) {
            throw new AssertionError("Couldn't find clade partition to remove.");
        }

        // System.out.println("> partition to remove picked");
        // System.out.println("partition  = " + nextPartition);
        // System.out.println("clade.prob = " + nextPartition.getParentClade().getProbability());
        // System.out.println("min p      = " + min);

        return nextPartition;
    }

    /**
     * Removes the given clade partition from the partial CCD and then tidies up the CCD;
     * any further removed clades and clade partitions are returned in lists.
     *
     * @param partitionToRemove the clade partition removed from the partial CCD
     * @return array of first list of clades removed and then list of clade partitions removed
     */
    private List[] reduceCCD(CladePartition partitionToRemove) {
        // removing the given partition from its parent clade might be all that is to do,
        // though it could lead to further removals
        Clade parent = partitionToRemove.getParentClade();

        // first the cases where other clades need to be removed as well
        // 1. parent clade has only the given clade partition
        if (parent.getPartitions().size() == 1) {
            return reduceCCD(parent);
        }

        // 2.left child clade has only one parent clade (via the given clade partition)
        Clade child = partitionToRemove.getChildClades()[0];
        if (child.getParentClades().size() == 1) {
            return reduceCCD(child);
        }

        // 2.right child clade has only one parent clade (via the given clade partition)
        child = partitionToRemove.getChildClades()[1];
        if (child.getParentClades().size() == 1) {
            return reduceCCD(child);
        }

        // 3. otherwise, only the given clade partition is removed
        parent.removePartition(partitionToRemove);
        List<CladePartition> removedPartitions = new ArrayList<>();
        removedPartitions.add(partitionToRemove);

        return new List[]{null, removedPartitions};
    }

    /**
     * Recompute the CCPs of all clade partitions in the current partial CCD, namely,
     * where clade partitions have been removed:
     * <ul>
     *  <li>renormalize for child clade partitions</li>
     *  <li>parent clade partitions loose probability proportionally</li>
     * </ul>
     *
     * @param cladesToHandle clades still in the partial CCD, sorted by increasing size, updates might be needed
     * @return the reduction in probability
     */
    private double recomputeCCPs(List<Clade> cladesToHandle) {
        // asserts that clades are sorted by increasing size
        for (Clade clade : cladesToHandle) {
            double sum = 0;
            for (CladePartition partition : clade.getPartitions()) {
                sum += partition.getCCP();
            }
            if (sum < 1) {
                // update parent clade partitions
                for (Clade parentClade : clade.getParentClades()) {
                    CladePartition parentPartition = parentClade.getCladePartition(clade);
                    parentPartition.setCCP(parentPartition.getCCP() * sum);
                }

                // renormalize
                for (CladePartition partition : clade.getPartitions()) {
                    partition.setCCP(Math.min(1.0, partition.getCCP() / sum));
                }
            }

            if (clade.isRoot()) {
                return sum;
            }
        }
        return -1;
    }
}
