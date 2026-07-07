package ccd.algorithms.credibleSets;

import beast.base.evolution.tree.Tree;
import ccd.model.AbstractCCD;
import ccd.model.CCD0;
import ccd.model.CCD2;
import ccd.model.Clade;
import ccd.model.CladePartition;
import ccd.model.ExtendedClade;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public abstract class CredibleCCDComputer implements ICredibleSet {

    /**
     * CCD for which this computer computes the credible set information.
     * Structure and parameters should not be altered (but values could be computed and set).
     */
    final AbstractCCD fullCCD;

    /**
     * Partial CCD with clades and clade splits removed from full CCD during the algorithm,
     * so not all method act as expected on a full CCD.
     */
    AbstractCCD partialCCD;

    /**
     * The credible set information computed for clade strategies:
     * For each clade, the minimum credible set it is contained in
     * (all are contained in the 100% one, not all in the 50% one)
     * with the smallest credible set equal to the last remaining tree.<br>
     * Based on clades of the full CCD.
     */
    Map<Clade, Double> cladeMinCredibility;

    /**
     * The credible set information computed for clade partition strategies:
     * For each clade partition, the minimum credible set it is contained in
     * (all are contained in the 100% one, not all in the 50% one)
     * with the smallest credible set equal to the last remaining tree.<br>
     * Based on clades partitions of the full CCD.
     */
    Map<CladePartition, Double> partitionMinCredibility;

    /**
     * Default constructor not to be used directly; use factory method {@link #getCredibleCCDComputer(AbstractCCD, CredibleSetType)}.
     *
     * @param baseCCD the CCD this credible CCD is based on
     */
    protected CredibleCCDComputer(AbstractCCD baseCCD) {
        this.fullCCD = baseCCD;
        this.cladeMinCredibility = new HashMap<>(fullCCD.getNumberOfClades());
        this.partitionMinCredibility = new HashMap<>(fullCCD.getNumberOfCladePartitions());
    }

    /**
     * Get a CredibleCCDComputer based on the given strategy for the given CCD.
     * Does not need to be called to compute credible levels of trees, because CCDs use them internally,
     * but might be helpful for more extensive experiments and direct access not provided by CCDs.
     *
     * @param baseCCD  the CCD the credibleCCD is based on
     * @param strategy the strategy/type of credible CCD use; not all strategies of enum are valid options
     * @return a CredibleCCDComputer based on the given strategy
     */
    public static CredibleCCDComputer getCredibleCCDComputer(AbstractCCD baseCCD, CredibleSetType strategy) {
        CredibleCCDComputer credibleCCD;

        credibleCCD = switch (strategy) {
            case CladeProbability -> new CladeProbabilityCredibleCCDConstructor(baseCCD);
            case PartitionProbability -> new PartitionProbabilityCredibleCCDComputer(baseCCD);
            // case CladeCount -> null;
            // case CladeTrees -> null;
            case TreeSampling, Frequency -> throw new IllegalArgumentException(
                    "The credible set strategy " + strategy.name() + " does not yield a credible CCD.");
            default -> throw new IllegalArgumentException("TODO implement");
        };

        credibleCCD.initializeCredibleSetInformation();

        return credibleCCD;
    }

    /** Minimum number of clades in CCD to remain to contain a full tree; used in construction algorithm. */
    int minNumClades;

    /** Remaining probability in current partial CCD maintained during construction algorithm. */
    double remainingProbability;

    /** Quick access to root clade for construction algorithm. */
    Clade rootClade;

    /**
     * Top level initialization and construction method, with specific work done by different subclasses.
     */
    private void initializeCredibleSetInformation() {
        // create copy from which we can remove clades and partitions
        partialCCD = fullCCD.copy();
        if (partialCCD instanceof CCD0) {
            ((CCD0) partialCCD).forbidReinitializing();
        }

        // ensure both CCDs are fully initialized
        fullCCD.getMaxTreeProbability();
        partialCCD.getMaxTreeProbability();

        fullCCD.computeCladeProbabilities();
        partialCCD.computeCladeProbabilities();

        fullCCD.computeCladeSumCladeCredibilities();
        // if (type == CredibleSetType.CladeTrees) {
        //     fullCCD.getNumberOfTrees();
        //     partialCCD.getNumberOfTrees();
        // }

        // needed during the algorithm
        rootClade = partialCCD.getRootClade();
        minNumClades = 2 * partialCCD.getNumberOfLeaves() - 1;
        remainingProbability = 1.0;

        computeCredibleSetInformationHelper();
    }

    /**
     * Set the credible level of the clades and clade partition removed in this round of the algorithm.
     *
     * @param removedObjects list of clades followed by list of clade partitions
     */
    protected void setCredibleLevels(List[] removedObjects) {
        if (removedObjects[0] != null) {
            for (Clade removedClade : (List<Clade>) removedObjects[0]) {
                setCredibleLevel(removedClade, remainingProbability);
            }
        }
        if (removedObjects[1] != null) {
            for (CladePartition removedPartition : (List<CladePartition>) removedObjects[1]) {
                setCredibleLevel(removedPartition, remainingProbability);
            }
        }
    }

    /**
     * Using sampling with rejection, get a tree from the alpha credible CCD.
     *
     * @param alpha credible level in (0,1]
     * @return tree sampled from alpha credible CCD
     */
    @Override
    public Tree sampleTreeFromCredibleSet(double alpha) {
        return null; // TODO
    }

    /**
     * Construct and return a credible CCD for the given credible level alpha,
     * that is, the base CCD of this computer is reduced to the clades and clade partitions
     * that make up alpha of the probability mass as computed.
     * Useful when you want to sample a lot  or test tree containment a lot for this credible CCD.
     *
     * @param alpha credible level
     * @return credible CCD for the given credible level alpha
     */
    public abstract AbstractCCD getCredibleCCD(double alpha);

    /**
     * @return the credible level information of clades
     */
    public Map<Clade, Double> getCladeMinCredibility() {
        if (cladeMinCredibility == null) {
            initializeCredibleSetInformation();
        }
        return cladeMinCredibility;
    }

    /**
     * @return the credible level information of clade partitions
     */
    public Map<CladePartition, Double> getPartitionMinCredibility() {
        return partitionMinCredibility;
    }

    /** Strategy/type-specific initialization and construction method. */
    protected abstract void computeCredibleSetInformationHelper();

    /**
     * Set the credible level information of the given clade (by the algorithm likely from the partial CCD)
     * for the corresponding clade of the underlying full CCD.
     *
     * @param clade of partial or full CCD; credible level information stored for corresponding clade in full CCD
     * @param alpha the credible level of the given clade
     */
    void setCredibleLevel(Clade clade, double alpha) {
        Clade originalClade;
        if (clade instanceof ExtendedClade && !clade.isLeaf()) { // CCD2
            originalClade = ((CCD2) fullCCD).getExtendedClade(clade.getCladeInBits(), ((ExtendedClade) clade).getSibling().getCladeInBits());
        } else {
            originalClade = fullCCD.getClade(clade.getCladeInBits());
        }
        cladeMinCredibility.put(originalClade, alpha);
    }

    /**
     * Set the credible level information of the given clade partition (by the algorithm likely from the partial CCD)
     * for the corresponding clade partition of the underlying full CCD.
     *
     * @param partition of partial or full CCD; credible level information stored for corresponding clade in full CCD
     * @param alpha     the credible level of the given clade
     */
    void setCredibleLevel(CladePartition partition, double alpha) {
        Clade fullParent;
        if (partition.getParentClade() instanceof ExtendedClade) { // CCD2
            fullParent = ((CCD2) fullCCD).getExtendedClade(partition.getParentClade().getCladeInBits(), ((ExtendedClade) partition.getParentClade()).getSibling().getCladeInBits());
        } else {
            fullParent = fullCCD.getClade(partition.getParentClade().getCladeInBits());
        }
        Clade fullSmallerChild;

        if (partition.getParentClade() instanceof ExtendedClade && !partition.getSmallerChild().isLeaf()) { // CCD2
            fullSmallerChild = ((CCD2) fullCCD).getExtendedClade(partition.getSmallerChild().getCladeInBits(), ((ExtendedClade) partition.getSmallerChild()).getSibling().getCladeInBits());
        } else {
            fullSmallerChild = fullCCD.getClade(partition.getSmallerChild().getCladeInBits());
        }
        for (CladePartition fullPartition : fullParent.getPartitions()) {
            if (fullPartition.getSmallerChild() == fullSmallerChild) {
                partitionMinCredibility.put(fullPartition, alpha);
                return;
            }
        }
    }

    /**
     * Removes the given clade from the underlying CCD and then proceeds to tidy up the CCD,
     * i.e. when the removal results in another clade without parent or without clade partitions,
     * that one gets removed as well, and so on.
     *
     * @param cladeToRemove the initial clade to remove from the underlying CCD
     * @return list of clades removed followed by list of clade partitions removed
     */
    List[] reduceCCD(Clade cladeToRemove) {
        /* NOTE: Known partial code redundancy with AbstractCCD#tidyUpCCDGraph
        for slightly different behavior and to extract more information. */

        double oldP = cladeToRemove.getProbability();

        List<Clade> removedClades = new ArrayList<>();
        List<CladePartition> removedPartitions = new ArrayList<>();

        ArrayList<Clade> cladesToRemove = new ArrayList<>();
        cladesToRemove.add(cladeToRemove);

        while (!cladesToRemove.isEmpty()) {
            Clade nextToRemove = cladesToRemove.remove(cladesToRemove.size() - 1);

            // do not allow removal of trivial clades
            if (nextToRemove.isLeaf() || (nextToRemove.isRoot())) {
                String errorMessage = "Building partial CCD for credible set, illegal request to remove " + (nextToRemove.isLeaf() ? "leaf" : "root") + "!";
                System.err.println("\n" + errorMessage);
                System.err.println("Initial clade requested to remove:   " + cladeToRemove);
                System.err.println("which had the following probability: " + oldP);
                System.err.println("Current clade processed to remove:   " + nextToRemove);
                System.err.println("CCD number of leaves:                " + partialCCD.getNumberOfLeaves());
                System.err.println("CCD number of clades remaining:      " + partialCCD.getNumberOfClades());
                System.err.println("Remaining clades: ");
                for (Clade clade : partialCCD.getClades()) {
                    if (!clade.isLeaf()) {
                        System.err.println(clade.getProbability() + " - " + clade);
                    }
                }
                throw new AssertionError(errorMessage);
            }

            // do not repeatedly process clades
            if (removedClades.contains(nextToRemove)) {
                continue;
            }
            removedClades.add(nextToRemove);

            // remove from mapping
            partialCCD.getCladeMapping().remove(nextToRemove.getCladeInBits());
            if (partialCCD instanceof CCD2) {
                partialCCD.getClades().remove(nextToRemove);
            }

            // remove connection to parent clades ...
            for (Clade parent : nextToRemove.getParentClades()) {
                parent.getChildClades().remove(nextToRemove);

                // ... and update parent clades
                for (CladePartition parentPartition : parent.getPartitions()) {
                    // there can only be one partition that contains
                    // cladeToRemove under the parent clade
                    if (parentPartition.containsChildClade(nextToRemove)) {
                        parent.getPartitions().remove(parentPartition);
                        removedPartitions.add(parentPartition);
                        if (parent.getPartitions().isEmpty()) {
                            cladesToRemove.add(parent);
                        }

                        // also update other child
                        Clade otherChild = parentPartition.getOtherChildClade(nextToRemove);
                        parent.getChildClades().remove(otherChild);
                        otherChild.getParentClades().remove(parent);
                        if (otherChild.getParentClades().isEmpty()) {
                            cladesToRemove.add(otherChild);
                        }

                        // as only one clade partition of the parent can contain the current clade
                        // we can break the loop here after having processed this one
                        break;
                    }
                }
            }

            // remove all connection to children
            for (CladePartition partition : nextToRemove.getPartitions()) {
                removedPartitions.add(partition);
                for (Clade child : partition.getChildClades()) {
                    child.getParentClades().remove(nextToRemove);
                    if (child.getParentClades().isEmpty()) {
                        cladesToRemove.add(child);
                    }
                }
            }

            // do not have to empty partitions, parent or child clades,
            // since we dump the clade anyway
        }

        return new List[]{removedClades, removedPartitions};
    }
}

    /*-
    private WrappedBeastTree prevMAPTree;
    private WrappedBeastTree currentMAPTree;
    private BufferedWriter writer;
    private String separator;
    private void writeHeader() {
        if (writer != null) {
            try {
                writer.write("p" + separator
                        + "#clades" + separator
                        + "#splits" + separator
                        + "log(#trees)" + separator
                        + "entropy" + separator
                        + "mapTreeChange");
                writer.newLine();
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        }
    }

    private void writeCurrentResult(double remainingProbability) {
        if (writer != null) {
            boolean mapTreeChanged = TreeDistances.robinsonsFouldDistance(currentMAPTree, prevMAPTree) != 0;

            try {
                writer.write(remainingProbability + separator
                        + partialCCD.getNumberOfClades() + separator
                        + partialCCD.getNumberOfCladePartitions() + separator
                        + logBigInteger(partialCCD.getNumberOfTrees()) + separator
                        + partialCCD.getEntropy() + separator
                        + mapTreeChanged);
                writer.newLine();
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        }
    }

    public static double logBigInteger(BigInteger val) {
        int precision = Math.max((int) (Math.log(val.bitLength()) / Math.log(2)), 20); // Ensure sufficient precision
        BigDecimal bigDecimalVal = new BigDecimal(val);
        int scale = bigDecimalVal.scale();

        // Scale value for improved precision
        BigDecimal scaledValue = bigDecimalVal.movePointLeft(scale);

        // Compute the logarithm using BigDecimal
        double log2 = Math.log(scaledValue.doubleValue());

        // Adjust the logarithm based on the scale
        return log2 + scale * Math.log(10);
    }
    */
