package ccd.algorithms.credibleSets;

import ccd.model.AbstractCCD;
import ccd.model.Clade;

/**
 * A credible CCD based on clades, intended for CCD0 models.
 * The credible set information is constructed by continuously removing the clade with the lowest probability.
 *
 * @author Jonathan Klawitter
 */
public class CladeProbabilityCredibleCCDConstructor extends CladeBasedCredibleCCDComputer {

    /**
     * Default constructor not to be used directly.
     *
     * @param baseCCD the CCD this credible CCD is based on
     */
    protected CladeProbabilityCredibleCCDConstructor(AbstractCCD baseCCD) {
        super(baseCCD);
    }

    @Override
    Clade getNextClade() {
        double min = Double.POSITIVE_INFINITY;
        Clade minClade = null;

        for (Clade clade : partialCCD.getClades()) {
            if (clade.isLeaf() || clade.isRoot()) {
                continue;
            }
            if (Math.abs(clade.getProbability() - 1.0) < AbstractCCD.PROBABILITY_ROUNDING_EPSILON) {
                continue;
            }

            double current = clade.getProbability();

            if ((current < min) || ((current == min) && (clade.size() > minClade.size()))) {
                min = current;
                minClade = clade;
            }
        }

        return minClade;
    }
}
