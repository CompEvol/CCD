package ccd.algorithms.credibleSets;

import beast.base.evolution.tree.Tree;
import ccd.model.AbstractCCD;
import ccd.model.CCD0;
import ccd.model.Clade;
import ccd.model.WrappedBeastTree;
import ccd.model.bitsets.BitSet;

import java.util.List;

/**
 * A credible CCD based on clades, intended for CCD0 models.
 * The credible set information is constructed by continuously removing the clade
 * chosen by specific strategy of a concrete child class.
 *
 * @author Jonathan Klawitter
 */
public abstract class CladeBasedCredibleCCDComputer extends CredibleCCDComputer {

    /**
     * Default constructor not to be used directly.
     *
     * @param baseCCD the CCD this credible CCD is based on
     */
    protected CladeBasedCredibleCCDComputer(AbstractCCD baseCCD) {
        super(baseCCD);
    }

    @Override
    public double getCredibleLevel(Tree tree) {
        WrappedBeastTree wrapped = new WrappedBeastTree(tree);
        double maxCredLevel = 0;

        for (BitSet bits : wrapped.getClades()) {
            double credLevel = getCredibleLevelOfClade(bits);
            if (credLevel == 0) {
                return -1.0;
            }
            if (credLevel > maxCredLevel) {
                maxCredLevel = credLevel;
            }
        }

        return maxCredLevel;
    }

    /* Helper method */
    private double getCredibleLevelOfClade(BitSet cladeInBits) {
        Clade clade = fullCCD.getClade(cladeInBits);
        return (clade == null) ? 0 : cladeMinCredibility.get(clade);
    }

    @Override
    public AbstractCCD getCredibleCCD(double alpha) {
        // TODO
        return null;
    }

    @Override
    protected void computeCredibleSetInformationHelper() {
        // double numCladesToHandle = partialCCD.getNumberOfClades() - minNumClades;
        do {
            Clade nextClade = getNextClade();

            // remove clade and tidy up CCD graph and store results
            List[] removedObjects = reduceCCD(nextClade);

            // then for the removed clades and partitions, set the resulting credible levels
            setCredibleLevels(removedObjects);

            // double remainingNontrivialClades = numCladesToHandle - (partialCCD.getNumberOfClades() - minNumClades);
            // while ((remainingNontrivialClades * 100 / numCladesToHandle) > donePercentage) {
            //     donePercentage++;
            //     // System.out.print(".");
            // }

            // recompute values in partial CCD
            // if (partialCCD instanceof CCD0) {
            partialCCD.resetSumCladeCredibilities();
            CCD0.setPartitionProbabilities(rootClade);
            partialCCD.computeCladeProbabilities();
            // }

            // recomputeValues;
            remainingProbability = rootClade.computeSumCladeCredibilities() / fullCCD.getRootClade().getSumCladeCredibilities();
            partialCCD.setCacheAsDirty();
            partialCCD.computeCladeProbabilities();
            // prevMAPTree = currentMAPTree;
            // currentMAPTree = new WrappedBeastTree(partialCCD.getMAPTree());
            // writeCurrentResult(remainingProbability);

            // we exist when we have reduced the CCD to a CCD on a single tree
        } while (partialCCD.getNumberOfClades() != minNumClades);

        for (Clade clade : partialCCD.getClades()) {
            setCredibleLevel(clade, remainingProbability);
        }

    }

    abstract Clade getNextClade();
}
