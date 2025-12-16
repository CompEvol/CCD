package ccd.algorithms.regularisation;

import ccd.algorithms.CCDCombiner;
import ccd.model.CCD1;
import ccd.model.Clade;
import ccd.model.CladePartition;

public class CCD1Regularisor {

    CCDRegularisationStrategy strategy;
    Double value;

    public CCD1Regularisor(CCDRegularisationStrategy strategy, double value) {
        this.strategy = strategy;
        this.value = value;
    }

    public CCD1Regularisor(CCDRegularisationStrategy strategy) {
        this.strategy = strategy;
    }

    public void regularize(CCD1 ccd) {
        boolean priorBasedStrategy = (strategy == CCDRegularisationStrategy.PriorOne) || (strategy == CCDRegularisationStrategy.PriorScaled);
        // we use a table for caching prior probabilities (if needed for strategy)
        double[][] priorProbabilityTable = priorBasedStrategy ? new double[ccd.getNumberOfLeaves()][] : null;

        for (Clade clade : ccd.getClades()) {
            if (clade.isLeaf() || clade.isCherry() || (clade.getNumberOfPartitions() == 1)) {
                continue;
            }

            double denominator = clade.getNumberOfOccurrences() + switch (strategy) {
                case AdditiveOne -> clade.getNumberOfPartitions();
                case AdditiveX -> clade.getNumberOfPartitions() * value;
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
