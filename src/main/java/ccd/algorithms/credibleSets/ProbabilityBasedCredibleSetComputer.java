package ccd.algorithms.credibleSets;

import beast.base.evolution.tree.Tree;
import ccd.model.ITreeDistribution;

import java.util.Arrays;

/**
 * A credible set method for a {@link ITreeDistribution}, i.e a CCD or a sample distribution,
 * based on sampling trees from the distribution and then using their probability as credible level thresholds.
 * More precisely, first sample a specified number of trees, compute their probabilities, and then sort them.
 * A tree with a higher probability than a sample tree at index i with alpha = i / numSamples has then credible level alpha.
 *
 * @author Jonathan Klawitter
 */
public class ProbabilityBasedCredibleSetComputer implements ICredibleSet {

    /** The default number of sampled trees used to compute the credible intervals. */
    public static final int DEFAULT_NUM_SAMPLES = 100000;

    /** The default number of thresholds used for credible levels implying the credible level precision. */
    public static final int DEFAULT_PRECISION = 100;

    /** The distribution this credible set is build upon. */
    private final ITreeDistribution treeDistribution;

    /** The number of samples used to build this credible set. */
    private final int numSamples;

    /** The number of thresholds used to determine credible levels ins this credible set. */
    private final int precision;

    /** The stored probability thresholds determining the credible levels. */
    private double[] sampledProbabilities;

    /**
     * A probability-based credible set with default number of sampled trees and default precision (number of thresholds)
     * on the given tree distribution.
     *
     * @param treeDistribution distribution of trees this credible set is based on
     */
    public ProbabilityBasedCredibleSetComputer(final ITreeDistribution treeDistribution) {
        this(treeDistribution, DEFAULT_NUM_SAMPLES, DEFAULT_PRECISION);
    }

    /**
     * A tree based credible set using the given number of sampled trees with default precision (number of thresholds)
     * on the given tree distribution.
     *
     * @param treeDistribution distribution of trees this credible set is based on
     * @param numberOfSamples  number of sampled trees used to determine credible levels
     */
    public ProbabilityBasedCredibleSetComputer(final ITreeDistribution treeDistribution, final int numberOfSamples) {
        this(treeDistribution, numberOfSamples, DEFAULT_PRECISION);
    }

    /**
     * A tree based credible set using the given number of sampled trees and the given precision (number of thresholds)
     * on the given tree distribution.
     *
     * @param treeDistribution distribution of trees this credible set is based on
     * @param numberOfSamples  number of sampled trees used to determine credible levels
     * @param precision        number of thresholds stored giving the credible level precision; preferably a divisor of numberOfSamples
     */
    public ProbabilityBasedCredibleSetComputer(ITreeDistribution treeDistribution, int numberOfSamples, int precision) {
        this.treeDistribution = treeDistribution;
        this.numSamples = numberOfSamples;

        if (numberOfSamples < precision) {
            System.err.printf("Request to create tree-based credible set with less sampled trees (%d) " +
                            "than precision (%d). Thus precision set to specified number of samples.",
                    numberOfSamples, precision);
            precision = numberOfSamples;
        }
        this.precision = precision;

        initializeCredibleSetInformation();
    }

    /* Initialization method. */
    private void initializeCredibleSetInformation() {
        sampledProbabilities = new double[numSamples + 1];
        // sampledLogProbabilities = new double[numSamples + 1];
        sampledProbabilities[0] = 0;
        // sampledLogProbabilities[0] = Double.NEGATIVE_INFINITY;
        for (int i = 1; i <= numSamples; i++) {
            sampledProbabilities[i] = treeDistribution.sampleTreeProbability();
            // System.out.println(sampledProbabilities[i]);
            // sampledLogProbabilities[i] = ccd.sampleTreeLogProbability();
        }

        Arrays.sort(sampledProbabilities);
        // reverse
        for (int i = 0; i < sampledProbabilities.length / 2; i++) {
            int rightIndex = sampledProbabilities.length - 1 - i;
            double tmp = sampledProbabilities[i];
            sampledProbabilities[i] = sampledProbabilities[rightIndex];
            sampledProbabilities[rightIndex] = tmp;
        }

        // filter to precision
        if (precision < numSamples) {
            double[] downsampledProbabilities = new double[precision];
            double stepSize = numSamples / (double) precision;
            for (int i = 1; i <= precision; i++) {
                int sampleIndex = (int) Math.round(i * stepSize);
                downsampledProbabilities[i - 1] = sampledProbabilities[sampleIndex];
            }
            sampledProbabilities = downsampledProbabilities;
        }

        // for (int i = 0; i < sampledProbabilities.length; i++) {
        //     System.out.printf("i = %d, threshold = %f\n", i, sampledProbabilities[i]);
        // }
    }

    @Override
    public double getCredibleLevel(Tree tree) {
        // Double prob = (Double) tree.getRoot().getMetaData(AbstractCCD.PROB_SUBTREE_KEY);
        double prob = treeDistribution.getProbabilityOfTree(tree);
        if (prob == 0) {
            return -1;
        }

        int indexOfNextSmallest = findIndexOfNextSmallest(prob) + 1;
        return (indexOfNextSmallest / (double) sampledProbabilities.length);
    }

    /* Helper method. */
    private int findIndexOfNextSmallest(double target) {
        // binary search
        int left = 0;
        int right = sampledProbabilities.length - 1;
        int result = right;

        // System.out.println("\ntarget = " + target);
        while (left <= right) {
            int mid = left + (right - left) / 2;
            if (target > sampledProbabilities[mid]) {
                result = mid;
                right = mid - 1; // search in left half
            } else {
                left = mid + 1; // search in right half
            }
        }

        // System.out.println("result = " + result);
        return result;
    }

    @Override
    public Tree sampleTreeFromCredibleSet(double alpha) {
        return null; // TODO
    }

    public double[] getNumberOfTreesPerCredibleLevel() {
        return getNumberOfTrees(false);
    }

    public double[] getLogNumberOfTreesPerCredibleLevel() {
        return getNumberOfTrees(true);
    }

    private double[] getNumberOfTrees(boolean inLogUnits) {
        double[] nums = new double[DEFAULT_PRECISION];
        double[] meanProbability = getMeanProbabilityPerCredibleLevel();
        double bucketProbability = 1.0 / DEFAULT_PRECISION;
        double logBucketProbability = Math.log(bucketProbability);
        // System.out.println("bucketProbability = " + bucketProbability);
        for (int i = 0; i < meanProbability.length; i++) {
            if (inLogUnits) {
                if (meanProbability[i] == 0) {
                    nums[i] = Double.POSITIVE_INFINITY; // or any other handling mechanism
                } else {
                    double logMeanProbability = Math.log(meanProbability[i]);
                    nums[i] = logBucketProbability - logMeanProbability;
                }
            } else {
                nums[i] = bucketProbability / meanProbability[i];
            }
        }


        return nums;
    }

    public double[] getMeanProbabilityPerCredibleLevel() {
        double[] meanProbability = new double[DEFAULT_PRECISION];
        int bucketLength = sampledProbabilities.length / DEFAULT_PRECISION;
        // System.out.println("sampledProbabilities: " + Arrays.toString(sampledProbabilities));
        // System.out.println("bucketLength = " + bucketLength);
        for (int i = 0; i < DEFAULT_PRECISION; i++) {
            double sum = 0.0;
            int index = i * bucketLength;
            for (int j = 0; j < bucketLength; j++) {
                sum += sampledProbabilities[index++];
            }
            meanProbability[i] = sum / bucketLength;
        }
        // System.out.println("mean probs: " + Arrays.toString(meanProbability));
        return meanProbability;
    }
}
