package ccd.algorithms.prior;

import org.apache.commons.numbers.combinatorics.BinomialCoefficient;

import java.util.Arrays;
import java.util.Random;
import java.util.stream.DoubleStream;

public class CladeSplitProbabilitySimulator {

    static int binomialCoefficientThreshold = 66; // result can fit in long
    private static long[][] combination;

    public CladeSplitProbabilitySimulator(int n) {
        setUpCache(n);
    }

    private void setUpCache(int n) {
        // combination: Cache table n choose k. For n up to binomialCoefficientThreshold (included)
        combination = new long[binomialCoefficientThreshold + 1][binomialCoefficientThreshold + 1];
        for (int i = 0; i <= binomialCoefficientThreshold; i++) {
            for (int j = 0; j <= i; j++) {
                combination[i][j] = BinomialCoefficient.value(i, j);
            }
        }
    }

    public static long binomialCoeff(int n, int k) {
        if (n <= binomialCoefficientThreshold) {
            return combination[n][k];
        } else {
            long x = 1;
            for (int i = 0; i < (k - 1); i++) {
                x *= (n - i) / (k - i);
            }
            return x;
        }
    }

    public static double getSample(double theta, double[] timesLeft, double[] timesRight) {

        double finalProbability = 1;

        Random random = new Random();

        // samplingTimes: distinct sampling times
        double[] samplingTimes = DoubleStream.concat(Arrays.stream(timesLeft), Arrays.stream(timesRight))
                .distinct().sorted().toArray();

        // sampledLeftAtTime, sampledRightAtTime: how many samples are added at a sampling time
        int[] sampledLeftAtTime = new int[samplingTimes.length];
        int[] sampledRightAtTime = new int[samplingTimes.length];

        int iLeft = 0;
        int iRight = 0;
        for (int i = 0; i < samplingTimes.length; i++) {
            while ((iLeft < timesLeft.length) && (timesLeft[iLeft] == samplingTimes[i])) {
                sampledLeftAtTime[i]++;
                iLeft++;
            }
            while ((iRight < timesRight.length) && (timesRight[iRight] == samplingTimes[i])) {
                sampledRightAtTime[i]++;
                iRight++;
            }
        }

        int currentTimeIndex = 0;
        double currentTime = samplingTimes[currentTimeIndex];
        double lastSamplingTime = samplingTimes[samplingTimes.length - 1];

        // initialise how many are active, how many are to be added
        int activeLeftSize = sampledLeftAtTime[currentTimeIndex];
        int activeRightSize = sampledRightAtTime[currentTimeIndex];

        int activeNodesCount;
        int totalPairCount;
        int validPairCount;
        int validPairCountLeft;
        int validPairCountRight;

        while (currentTime != lastSamplingTime) {

            activeNodesCount = activeLeftSize + activeRightSize;
            totalPairCount = activeNodesCount * (activeNodesCount - 1) / 2; // activeNodesCount choose 2
            validPairCountLeft = activeLeftSize * (activeLeftSize - 1) / 2; // activeLeft choose 2
            validPairCountRight = activeRightSize * (activeRightSize - 1) / 2; // activeRight choose 2
            validPairCount = validPairCountLeft + validPairCountRight;

            // if there's only 1 active node left, cannot coalesce. Go to the next sampling time
            if (activeNodesCount == 1) {
                currentTimeIndex += 1;
                currentTime = samplingTimes[currentTimeIndex];
                activeLeftSize += sampledLeftAtTime[currentTimeIndex];
                activeRightSize += sampledRightAtTime[currentTimeIndex];
            } else { // there are more than 1 active node, can coalesce
                // if there is valid pair, draw time
                if (validPairCount != 0) {
                    double rate = (activeNodesCount * (activeNodesCount - 1.0)) / (theta * 2.0);
                    double x = -Math.log(random.nextDouble()) / rate;
                    currentTime += x;
                    // if time too large, go to next sampling time
                    // if currentTime (the time after adding the coalescent time increment) > the next sampling time
                    if (currentTime > samplingTimes[currentTimeIndex + 1]) {
                        currentTimeIndex += 1;
                        currentTime = samplingTimes[currentTimeIndex];
                        activeLeftSize += sampledLeftAtTime[currentTimeIndex];
                        activeRightSize += sampledRightAtTime[currentTimeIndex];
                    } else { // At least one group has >= 2 active lineages that can coalesce
                        int randomNum = random.nextInt(validPairCount);
                        if (randomNum < validPairCountLeft) { // left pair coalescence
                            activeLeftSize -= 1; // -2 children, + 1 parent
                        } else { // right pair coalescence
                            activeRightSize -= 1;
                        }
                        finalProbability *= (double) validPairCount / totalPairCount;
                    }
                } else { // no valid pari, go to next time, p_nc
                    // the probability that no coalescent before next sample, exp( -tau * (k choose 2) / theta )
                    currentTimeIndex += 1;
                    double nextTime = samplingTimes[currentTimeIndex];
                    double tau = nextTime - currentTime;
                    double noCoalescentProb = Math.exp(-tau * ((double) totalPairCount / theta));
                    finalProbability *= noCoalescentProb;
                    currentTime = samplingTimes[currentTimeIndex];
                    activeLeftSize += sampledLeftAtTime[currentTimeIndex];
                    activeRightSize += sampledRightAtTime[currentTimeIndex];
                }
            }
        }

        finalProbability *= contempFormula((activeLeftSize + activeRightSize), activeLeftSize);

        return finalProbability;
    }

    public static double contempFormula(int parentSize, int leftSize) {
        return (2.0 / (parentSize - 1)) * (1.0 / binomialCoeff(parentSize, leftSize));
    }

    public static void main(String[] args) { // small test case

        double theta = 1.0;
        double[] timesLeft = new double[]{0, 0};
        double[] timesRight = new double[]{1};

        double totalProb = 0;
        int simReps = 10000;

        int numOfTaxa = 3;
        CladeSplitProbabilitySimulator testSim = new CladeSplitProbabilitySimulator(numOfTaxa);

        for (int i = 0; i < simReps; i++) {
            totalProb += testSim.getSample(theta, timesLeft, timesRight);
        }
        double meanProb = totalProb / simReps;
        System.out.println(meanProb);
        // The correct result should approach 0.7547470392190384
    }

}