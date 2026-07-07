package ccd.algorithms.regularisation;

import ccd.algorithms.CCDCombiner;
import ccd.model.CCD1;
import ccd.model.Clade;
import ccd.model.CladePartition;

public class CCD1Regularisor {

    CCDRegularisationStrategy strategy;
    Double value;
    /** Second pseudocount (beta) for NNI-derived splits under {@link CCDRegularisationStrategy#AdditiveXY}. */
    Double value2;

    public CCD1Regularisor(CCDRegularisationStrategy strategy, double value) {
        this.strategy = strategy;
        this.value = value;
    }

    /**
     * Two-level additive smoothing: pseudocount {@code value} (alpha) for
     * observed and split-expanded splits, and {@code value2} (beta) for
     * NNI-derived splits (those touching an NNI-expanded clade). Intended with
     * {@link CCDRegularisationStrategy#AdditiveXY}.
     *
     * @param strategy regularisation strategy (use AdditiveXY)
     * @param value    pseudocount alpha for regular splits
     * @param value2   pseudocount beta for NNI-derived splits
     */
    public CCD1Regularisor(CCDRegularisationStrategy strategy, double value, double value2) {
        this.strategy = strategy;
        this.value = value;
        this.value2 = value2;
    }

    public CCD1Regularisor(CCDRegularisationStrategy strategy) {
        this.strategy = strategy;
    }

    /** Whether a split (parent -&gt; {child0, child1}) is NNI-derived, i.e. touches a novel clade. */
    public static boolean isNNISplit(Clade parent, CladePartition partition) {
        return parent.isNNIExpanded()
                || partition.getChildClades()[0].isNNIExpanded()
                || partition.getChildClades()[1].isNNIExpanded();
    }

    /* Two-level additive smoothing for one clade. */
    private void regulariseAdditiveXY(Clade clade) {
        double alpha = this.value;
        double beta = this.value2;

        double denominator = clade.getNumberOfOccurrences();
        for (CladePartition partition : clade.getPartitions()) {
            denominator += isNNISplit(clade, partition) ? beta : alpha;
        }

        for (CladePartition partition : clade.getPartitions()) {
            if (partition.isCCPSet()) {
                throw new AssertionError("Cannot regularize CCD1 with set CCPs, only with sample-based ones.");
            }
            double addend = isNNISplit(clade, partition) ? beta : alpha;
            partition.setCCP((partition.getNumberOfOccurrences() + addend) / denominator);
        }
    }

    public void regularise(CCD1 ccd) {
        boolean priorBasedStrategy = (strategy == CCDRegularisationStrategy.PriorOne) || (strategy == CCDRegularisationStrategy.PriorScaled);
        // we use a table for caching prior probabilities (if needed for strategy)
        double[][] priorProbabilityTable = priorBasedStrategy ? new double[ccd.getNumberOfLeaves()][] : null;

        for (Clade clade : ccd.getClades()) {
            if (clade.isLeaf() || clade.isCherry() || (clade.getNumberOfPartitions() == 1)) {
                continue;
            }

            if (strategy == CCDRegularisationStrategy.AdditiveXY) {
                regulariseAdditiveXY(clade);
                continue;
            }

            double denominator = clade.getNumberOfOccurrences() + switch (strategy) {
                case AdditiveOne -> clade.getNumberOfPartitions();
                case AdditiveX -> clade.getNumberOfPartitions() * value;
                case AdditiveXY -> throw new AssertionError("AdditiveXY handled separately");
                case PriorOne -> 1;
                case PriorScaled -> value;
            };

            double priorSum = 0;
            double[] priorProbs = priorBasedStrategy ? new double[clade.getNumberOfPartitions()] : null;
            int i = 0;
            if (priorBasedStrategy) {
                // compute all prior CCPs (assuming contemporaneous case)
                for (CladePartition partition : clade.getPartitions()) {
                    double priorP = CCDCombiner.getPartitionPriorProbability(partition, priorProbabilityTable);
                    priorProbs[i++] = priorP;
                    priorSum += priorP;
                }

                // normalize prior probabilities
                for (i = 0; i < priorProbs.length; i++) {
                    priorProbs[i] = priorProbs[i] / priorSum;
                }
            }

            i = 0;
            for (CladePartition partition : clade.getPartitions()) {
                if (partition.isCCPSet()) {
                    throw new AssertionError("Cannot regularize CCD1 with set CCPs, only with sample-based ones.");
                }

                double addend = switch (strategy) {
                    case AdditiveOne -> 1;
                    case AdditiveX -> value;
                    case AdditiveXY -> throw new AssertionError("AdditiveXY handled separately");
                    case PriorOne -> priorProbs[i];
                    case PriorScaled -> priorProbs[i] * value;
                };

                double p = (partition.getNumberOfOccurrences() + addend) / denominator;
                partition.setCCP(p);

                i++;
            }
        }

    }

    @Override
    public String toString() {
        return (this.value == null) ? this.strategy.toString() : this.strategy.toString() + " with value " + this.value;
    }

}
