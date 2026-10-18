package ccd.algorithms.prior;

import org.apache.commons.numbers.core.ArithmeticUtils;
import org.apache.commons.numbers.combinatorics.BinomialCoefficient;

import java.io.IOException;
import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.stream.DoubleStream;

/**
 * TODO
 *
 * @author Sophie Yang, Remco Bouckaert, Jonathan Klawitter
 */
public class PriorCoalescentCCPComputer {

    private enum Approach {Exact, Heuristic}

    public static final int BIG_DECIMAL_PRECISION = 48; // 36
    public static final MathContext MC = new MathContext(BIG_DECIMAL_PRECISION, RoundingMode.HALF_UP);
    public static final BigDecimal ROUNDING_EPSILON = new BigDecimal(-1e-10);

    /**
     * Construct new CladePartitionProbabilityComputer with caches set up for fast computation.
     *
     * @param numTaxa           size of largest parent clade; used to set up caches
     * @param maxLargeCladeSize size of the largest larger child clade; used for appropriate cache sizes
     * @param maxSmallCladeSize size of the largest smaller child clade
     */
    public PriorCoalescentCCPComputer(int numTaxa, int maxLargeCladeSize, int maxSmallCladeSize) {
        if (maxLargeCladeSize >= numTaxa) {
            throw new AssertionError(String.format("Child clade sizes (%d, %d) should be" +
                            "at most number of largest parent clade size (%d) minus one.",
                    maxLargeCladeSize, maxSmallCladeSize, numTaxa));
        }
        if (maxLargeCladeSize < maxSmallCladeSize) {
            throw new AssertionError(String.format("Largest large child clade size (%d) should be" +
                            "at least largest smaller child clade size (%d).",
                    maxLargeCladeSize, maxSmallCladeSize));
        }

        setUpCache(numTaxa, maxLargeCladeSize, maxSmallCladeSize);
    }

    /**
     * Initialize caches of frequently computed values.
     *
     * @param numTaxa        size of the largest parent clade
     * @param largeCladeSize size of the largest larger child clade
     * @param smallCladeSize size of the largest smaller child clade
     */
    private void setUpCache(int numTaxa, int largeCladeSize, int smallCladeSize) {
        initializeBinomialCoefficientCache(numTaxa, largeCladeSize);
        initializePairPickingCache(largeCladeSize, smallCladeSize);
        initializeFallDivRiseTermCache(numTaxa);
        initializeParenthesesDivFactorialsCache(numTaxa);
    }

    /**
     *
     * @param theta: rescaled eﬀective population size
     * @param timesLeft: left child clade taxa sampled time
     * @param timesRight: right child clade taxa sampled time
     */
    public double getPriorCCPExact(double theta, double[] timesLeft, double[] timesRight) throws IOException {
        return getPriorCCP(theta, timesLeft, timesRight, Approach.Exact, 0);
    }

    /**
     *
     * @param theta: rescaled eﬀective population size
     * @param timesLeft: left child clade taxa sampled time
     * @param timesRight: right child clade taxa sampled time
     * @param epsilon: threshold for heuristic
     */
    public double getPriorCCPHeuristically(double theta, double[] timesLeft, double[] timesRight, double epsilon) throws IOException {
        return getPriorCCP(theta, timesLeft, timesRight, Approach.Heuristic, epsilon);
    }

    private double getPriorCCP(double theta, double[] timesLeft, double[] timesRight, Approach approach, double epsilon) throws IOException {
        Arrays.sort(timesLeft);
        Arrays.sort(timesRight);

        // we want to assume that the left clade is larger than the right one
        // if they are equal in length, we assume the clade sampled first is the left one
        if ((timesLeft.length < timesRight.length) || ((timesLeft.length == timesRight.length) && (timesLeft[0] > timesRight[0]))) {
            double[] tmp = timesLeft;
            timesLeft = timesRight;
            timesRight = tmp;
        }

        // for (int i = 0; i < timesLeft.length; i++) {
        //    System.out.println("timesLeft = " + timesLeft[i]);
        //}
        // for (int i = 0; i < timesRight.length; i++) {
        //    System.out.println("timesRight = " + timesRight[i]);
        //}

        // nLeft, nRight: total number of taxa in the left and right clade
        int nLeft = timesLeft.length;
        int nRight = timesRight.length;

        // times: distinct sampling times
        double[] samplingTimes = DoubleStream.concat(Arrays.stream(timesLeft), Arrays.stream(timesRight))
                .distinct().sorted().toArray();

        // sampledLeftAtTime, sampledRightAtTime: how many samples are added at a sampling time
        int[] sampledLeftAtTime = new int[samplingTimes.length];
        int[] sampledRightAtTime = new int[samplingTimes.length];

        // taus: sampling intervals
        double[] taus = new double[samplingTimes.length - 1];

        // write the values of sampledLeftAtTime and sampledRightAtTime
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
            if (i != (samplingTimes.length - 1)) {
                taus[i] = samplingTimes[i + 1] - samplingTimes[i];
            }
        }

        // dimensions of the table of the DP
        // since we can have 0 to n1 active left, we have nLeft + 1 rows
        // since we can have 0 to n2 active right, we have nRight + 1 columns
        final int rows = nLeft + 1;
        final int cols = nRight + 1;
        // use flattened table for speedup
        // thus T[i][j] = T[i * cols + j]
        final BigDecimal[] T = new BigDecimal[rows * cols];
        // furthermore, the DP loops over several tables in theory,
        // but we can use a single table and overwrite previous values
        // since any entry T[i,j] at time i only depends
        // on entries T[k,l] at time i-1 with index k >= i and l >=j

        // initializing the table with initially active lineages
        int leftSampledSoFar = sampledLeftAtTime[0];
        int rightSampledSoFar = sampledRightAtTime[0];
        T[leftSampledSoFar * cols + rightSampledSoFar] = BigDecimal.ONE;

        // delta: number of new taxa added at a sampling time
        int deltaJ = leftSampledSoFar;
        int deltaK = rightSampledSoFar;

        // for heuristic method
        BigDecimal epsilonBig = new BigDecimal(epsilon);
        BigDecimal threshold = BigDecimal.ZERO;

        // main loop of the DP
        for (int i = 0; i < taus.length; i++) { // loop through tau
            double tau = taus[i];

            // System.out.println("\ntau = " + tau);

            BigDecimal maxP = BigDecimal.ZERO; // max cell in the prev table
            if (approach == Approach.Heuristic) {
                for (int j = 0; j <= leftSampledSoFar; j++) {
                    for (int k = 0; k <= rightSampledSoFar; k++) {
                        if ((T[j * cols + k] != null) && (T[j * cols + k].compareTo(maxP) > 0)) {
                            maxP = T[j * cols + k];
                        }
                    }
                }
                // update threshold
                threshold = epsilonBig.multiply(maxP);
            }

            // find rectangle in matrix containing entries > epsilon
            int boundLoJ = leftSampledSoFar;
            int boundLoK = rightSampledSoFar;
            int boundHiJ = 0;
            int boundHiK = 0;
            if (approach == Approach.Heuristic) {
                for (int j = leftSampledSoFar; j >= 0; j--) {
                    for (int k = rightSampledSoFar; k >= 0; k--) {
                        if (T[j * cols + k] == null) {
                            continue;
                        }
                        if (T[j * cols + k].compareTo(threshold) > 0) {
                            boundLoJ = Math.min(boundLoJ, j);
                            boundLoK = Math.min(boundLoK, k);
                            boundHiJ = Math.max(boundHiJ, j);
                            boundHiK = Math.max(boundHiK, k);
                        }
                    }
                }
            }

            // rows
            for (int j = 0; j <= leftSampledSoFar; j++) {
                // if at least one left lineage has been sampled,
                // we can skip row 0, since we cannot go back to zero left lineages
                if (leftSampledSoFar != 0 && j == 0) continue;

                // columns
                for (int k = 0; k <= rightSampledSoFar; k++) {
                    // if at least one right lineage has been sampled,
                    // we can skip column 0, since we cannot go back to zero right lineages
                    if (rightSampledSoFar != 0 && k == 0) continue;

                    // System.out.println("row j = " + j + ", column k = " + k);

                    BigDecimal prob = BigDecimal.ZERO;

                    // get the overlap area
                    int minJ, maxJ, minK, maxK;
                    if (approach == Approach.Exact) {
                        minJ = j;
                        maxJ = leftSampledSoFar;
                        minK = k;
                        maxK = rightSampledSoFar;
                    } else {
                        minJ = Math.max(j, Math.max(deltaJ, boundLoJ));
                        maxJ = Math.min(boundHiJ, leftSampledSoFar);
                        minK = Math.max(k, Math.max(deltaK, boundLoK));
                        maxK = Math.min(boundHiK, rightSampledSoFar);
                    }

                    // System.out.println("minJ = " + minJ + ", maxJ = " + maxJ + ", minK = " + minK + ", maxK = " + maxK);

                    // recall that we can overwrite old entries,
                    // because we only use values to the right and below
                    for (int jPrev = minJ; jPrev <= maxJ; jPrev++) {
                        for (int kPrev = minK; kPrev <= maxK; kPrev++) {
                            BigDecimal prevProb = T[jPrev * cols + kPrev];
                            // System.out.println("prevProb = " + prevProb);

                            if (prevProb.compareTo(BigDecimal.valueOf(0)) == 0) continue;

                            // System.out.println("jPrev = " + jPrev + ", kPrev = " + kPrev);

                            int x = jPrev + kPrev;
                            int y = j + k;
                            BigDecimal coalesProb1 = coalescentFormula(x, y, tau, theta);
                            BigDecimal coalesProb2 = getPairPickingValue(jPrev, kPrev, j, k);
                            BigDecimal term = coalesProb1.multiply(coalesProb2).multiply(prevProb);

                            prob = prob.add(term);
                            prob = checkProbability(prob);
                        }
                    }
                    T[j * cols + k] = prob;
                }
            }

            // update the number of lineages at sampling time
            deltaJ = sampledLeftAtTime[i + 1];
            deltaK = sampledRightAtTime[i + 1];
            leftSampledSoFar += deltaJ;
            rightSampledSoFar += deltaK;

            // since we sample new lineages,
            // we shift the whole table correspondingly from bottom right to top left
            // System.out.println("shift");
            for (int j = (leftSampledSoFar - deltaJ); j >= 0; j--) {
                if (((leftSampledSoFar - deltaJ) != 0) && (j == 0)) continue;
                for (int k = (rightSampledSoFar - deltaK); k >= 0; k--) {
                    if (((rightSampledSoFar - deltaK) != 0) && (k == 0)) continue;
                    BigDecimal value = T[j * cols + k];
                    T[j * cols + k] = BigDecimal.ZERO;
                    T[(j + deltaJ) * cols + k + deltaK] = value;
                }
            }
        }

        // we have now reached the contemporaneous case,
        // so we combine all entries  with the probability of combining into the right clades now
        BigDecimal finalProbability = BigDecimal.ZERO;
        for (int prev_j = 1; prev_j <= leftSampledSoFar; prev_j++) {
            for (int prev_k = 1; prev_k <= rightSampledSoFar; prev_k++) {
                BigDecimal prev_prob = T[prev_j * cols + prev_k];
                if (!prev_prob.equals(BigDecimal.ZERO)) {
                    BigDecimal contemporaneousProb = getContemporaneousProbability(prev_j + prev_k, prev_j);
                    finalProbability = finalProbability.add(prev_prob.multiply(contemporaneousProb));
                    // checkProbability(finalProbability);
                }
            }
        }

        return finalProbability.doubleValue();
    }

    /**
     * Compute the probability that the given number of (previously) active lineages, x,
     * coalesce to the given number of now active lineages, y,
     * for a given time interval tau and given.
     * The formula used can be found in <a href="https://doi.org/10.1016/0040-5809(84)90027-3">Tavaré, 1984
     * (Line-of-Descent and Genealogical Processes, and Their Applications in Population Genetics Models)</a> Eq (6.1).
     *
     * @param prevActiveLineages number of lineages to start with
     * @param nowActiveLineages  number of lineages after coalescent events
     * @param tau                time interval
     * @param theta              population size
     * @return probability of coalescence of start lineages to target lineages for tau and theta
     */
    private BigDecimal coalescentFormula(int prevActiveLineages, int nowActiveLineages, double tau, double theta) {
        // System.out.println("prev = " + prevActiveLineages + ", now = " + nowActiveLineages + ", tau = " + tau + ", theta = " + theta);
        /*- For
            x = # prev lineages
            y = # remaining lineages
            T = time = tau / theta
         we are computing the probability of x lineages coalescing to y lineages in T
         with the following formula:

            g_{x,y}(T) = \sum_{k=y}^{x} e^{-k(k-1)T/2} (-1)^{k-y} (2k-1) \frac{ y_(k-1) x_[k] }{ y!(k-y)!x_(k) }

         where
            x_[k] is the falling factorial: x_[k] = x * (x-1) * (x-2) * ... * (x - k + 1)
            x_(k) is the rising factorial: x_(k) = x * (x+1) * (x+2) * ... * (x + k - 1)

         Note that in the paper, the notations is p(k,t) = e^{-k(k-1)T/2}.

         To implement this formula efficiently and trying to avoid over/underflows, we are breaking it apart:
         1. p(k,t) * (-1)^{k-y}
         2. fallDivRiseTerm(x, k)       = (2k-1) \frac{         x_[k] }{          x_(k) }
         3. riseDivFactorialsTerm(y, k) =        \frac{ y_(k-1)       }{ y!(k-y)!       }

         */

        BigDecimal probBig = new BigDecimal(0);
        final double time = tau / (2 * theta);
        final BigDecimal timeBig = new BigDecimal(tau / (2 * theta), MC); // term T/2, where T = tau / theta
        // System.out.println("time = " + time);

        boolean positiveSign = true;
        final BigDecimal minusOne = BigDecimal.valueOf(-1);
        if (nowActiveLineages != 1) {
            // System.out.println("if nowActiveLineages != 1");
            for (int k = nowActiveLineages; k <= prevActiveLineages; k++) {
                // System.out.println("coalesce formula k = " + k);
                BigDecimal exponent = new BigDecimal(-k * (k - 1)).multiply(timeBig);
                BigDecimal expPartBig = exp(exponent);
                // System.out.println("*k = " + k);
                // System.out.println("exponent = " + exponent);
                // System.out.println("expPart = " + expPartBig);

                BigDecimal fractionBig = getFallDivRiseTerm(prevActiveLineages, k)
                        .multiply(getRiseDivFactorialsTerm(nowActiveLineages, k));
                // System.out.println("fraction = " + fractionBig.doubleValue());

                // System.out.println("positiveSign = " + positiveSign);
                BigDecimal term = positiveSign ? expPartBig : expPartBig.multiply(minusOne);
                // System.out.println("term before times fraction = " + term);
                term = term.multiply(fractionBig, MC);
                // System.out.println("term = " + term.doubleValue());

                probBig = probBig.add(term, MC);
                // System.out.println("term = " + term.doubleValue() + ", prob = " + probBig.doubleValue());

                positiveSign = !positiveSign;
            }
        } else {
            // System.out.println("else nowActiveLineages == 1");
            // the formula is simpler if only one lineage remains; see also Tavaré Eq 6.2
            for (int k = 2; k <= prevActiveLineages; k++) {
                double expPart = Math.exp(-k * (k - 1) * time);
                BigDecimal expPartBig = BigDecimal.valueOf(expPart);

                BigDecimal fractionBig = getFallDivRiseTerm(prevActiveLineages, k);

                BigDecimal term = positiveSign ? expPartBig : expPartBig.multiply(minusOne);
                positiveSign = !positiveSign;
                term = term.multiply(fractionBig);

                probBig = probBig.add(term);
                // System.out.println("'k = " + k + ", term = " + term.doubleValue() + ", prob = " + probBig.doubleValue());
            }

            probBig = BigDecimal.ONE.subtract(probBig);
        }

        // System.out.println("before check prob in coalescent formula");
        probBig = checkProbability(probBig);
        return probBig;
    }

    /**
     * Get the probability of a clade split of the parent clade of given size such that one child clade has the given size
     * under the coalescent model with all leaves at the same time.
     *
     * @param parentSize size of parent clade C
     * @param childSize  size of child clade A in a clade split {A, B} of C
     * @return CCP of clade split {A, B} of parent clade C under contemporaneous coalescent case
     */
    public static BigDecimal getContemporaneousProbability(int parentSize, int childSize) {
        // System.out.println("parentSize = " + parentSize + ", childSize = " + childSize);
        BigDecimal r;
        if (parentSize > BINOMIAL_COEFFICIENT_LONG_THRESHOLD) {
            r = BigDecimal.ONE.divide(new BigDecimal(getBinomialCoefficientDouble(parentSize, childSize)), MC);
        } else {
            r = BigDecimal.ONE.divide(new BigDecimal(getBinomialCoefficientLong(parentSize, childSize)), MC);
        }
        BigDecimal term = new BigDecimal(2).divide(new BigDecimal(parentSize - 1), MC);
        BigDecimal prob = term.multiply(r);
        prob = checkProbability(prob);
        return prob;
    }


    /* -- FALL-DIV-RISE CACHE -- */

    /**
     * Cache for (2k-1) * x_[k] / x_(k); use method to access values.<br>
     * x_[k] is the falling factorial: x_[k] = x * (x-1) * (x-2) * ... * (x - (k + 1))<br>
     * x_(k) is the rising factorial: x_(k) = x * (x+1) * (x+2) * ... * (x + (k - 1))<br>
     * First dimension: x, i.e., number of lineages to coalesce, (or prevActiveLineages)<br>
     * Second dimension: k, i.e., number of factors in the product
     */
    private BigDecimal[][] fallDivRiseCache;

    /**
     * Initialize cache for falling/rising factorials part of coalescent formula.
     *
     * @param maxLineages maximum number of active lineages that may coalesce; corresponds to x in the Eq 6.2 of Tavaré
     */
    private void initializeFallDivRiseTermCache(int maxLineages) {

        // there is always at least one now active lineage
        fallDivRiseCache = new BigDecimal[maxLineages][];
        for (int x = 1; x <= fallDivRiseCache.length; x++) {
            int xIndex = x - 1;

            // define the length of the 2nd dimension of the cache
            // always at least one remaining lineage
            fallDivRiseCache[xIndex] = new BigDecimal[x];

            for (int k = 1; k <= x; k++) {
                int kIndex = k - 1;

                /*-
                 For x lineages and variable k, we are computing
                    fallDivRiseTerm(x, k) = (2k-1)  x_[k] / x_(k)
                 where
                    x_[k] is the falling factorial: x_[k] = x * (x-1) * (x-2) * ... * (x - (k - 1))
                    x_(k) is the rising factorial:  x_(k) = x * (x+1) * (x+2) * ... * (x + (k - 1))

                 To avoid huge numbers, we directly compute fractions (x-i)/(x+i).
                 */

                BigDecimal result = new BigDecimal(2 * k - 1); // result is initialised to be the poly term

                for (int i = 1; i < k; i++) { // for loop through terms in the falling / rising factorial
                    BigDecimal term = new BigDecimal(x - i).divide(new BigDecimal(x + i), MC);
                    result = result.multiply(term);
                }

                fallDivRiseCache[xIndex][kIndex] = result;
            }
        }
    }

    /**
     * @param prevActiveLineages previously active lineages
     * @param k                  number of steps
     * @return fall-div-rise term of coalescent formula
     */
    private BigDecimal getFallDivRiseTerm(int prevActiveLineages, int k) {
        return fallDivRiseCache[prevActiveLineages - 1][k - 1];
    }


    /* -- RISE-DIV-FACTORIALS CACHE -- */

    /**
     * Cache for y_(k-1) / [ y!(k-y)! ]; use method to access values.<br>
     * x_(k) is the rising factorial: x_(k) = x * (x+1) * (x+2) * ... * (x + (k - 1))<br>
     * First dimension: y, i.e., number of lineages to remain, (or nowActiveLineages)<br>
     * Second dimension: k, i.e., number of factors in the product
     */
    private BigDecimal[][] riseDivFactorialsCache;

    /**
     * Initialize cache for rising factorial / pair of factorials part of coalescent formula.
     *
     * @param maxLineages maximum number of lineages that may coalesce; corresponds to x in the Eq 6.2 of Tavaré
     */
    private void initializeParenthesesDivFactorialsCache(int maxLineages) {
        // there is always at least one now active lineage y and formula equals 1 for y=1,
        riseDivFactorialsCache = new BigDecimal[maxLineages - 1][];
        for (int y = 2; y <= maxLineages; y++) {
            int yIndex = y - 2;

            // note that k can go up to max number of lineages, but starts always at y
            riseDivFactorialsCache[yIndex] = new BigDecimal[maxLineages - y + 1];

            for (int k = y; k <= maxLineages; k++) {
                int kIndex = k - y;

                /*-
                 For x lineages and variable k, with k >= y, we are computing
                    riseDivFactorialsTerm(y, k) = y_(k-1) / [ y!(k-y)! ]
                        = y * (y+1) * ... * (y + (k-2))   // = y_(k-1); k-1 terms is always longest
                        /[y * (y-1) * ... * 2             // = y!;      y-1 terms
                          1 * 2 * ...  * (k - y)]         // = (k-y)!;  k-y terms
                 since
                    x_(k) is the rising factorial:  x_(k) = x * (x+1) * (x+2) * ... * (x + (k-1))

                 To avoid huge numbers, we directly compute fractions per "column".
                 Note that we can skip the first "column" y/(y*1).
                 We have to be smart about going through the columns, because they have different length.
                 */

                BigDecimal result = BigDecimal.ONE;
                int diff = k - y;
                if (diff < (y - 1)) {
                    // so third row is shortest
                    for (int i = 1; i < diff; i++) {
                        BigDecimal numerator = new BigDecimal(y + i);
                        BigDecimal denominator = new BigDecimal((y - i) * (i + 1));
                        result = result.multiply(numerator).divide(denominator, MC);
                    }

                    // now we have already done diff-1 "columns" (and also first "column"),
                    // so there remain y-1-(diff-1) = y - diff terms of second line (i.e., y!)
                    for (int i = diff; i < y - 1; i++) {
                        BigDecimal numerator = new BigDecimal(y + i);
                        BigDecimal denominator = new BigDecimal(y - i);
                        result = result.multiply(numerator).divide(denominator, MC);
                    }

                    // now we have done y-1 "columns" and only first row remains
                    for (int i = y - 1; i < k - 1; i++) {
                        BigDecimal nominator = new BigDecimal(y + i);
                        result = result.multiply(nominator);
                    }
                } else {
                    // so second row is shortest
                    for (int i = 1; i < y - 1; i++) {
                        BigDecimal numerator = new BigDecimal(y + i);
                        BigDecimal denominator = new BigDecimal((y - i) * (i + 1));
                        result = result.multiply(numerator).divide(denominator, MC);
                    }

                    // now we have already done y-2 "columns" (and also first "column"),
                    for (int i = y - 1; i < diff; i++) {
                        BigDecimal numerator = new BigDecimal(y + i);
                        BigDecimal denominator = new BigDecimal(i + 1);
                        result = result.multiply(numerator).divide(denominator, MC);
                    }

                    // now we have done diff "columns" and only first row remains
                    for (int i = diff; i < k - 1; i++) {
                        BigDecimal numerator = new BigDecimal(y + i);
                        result = result.multiply(numerator);
                    }
                }

                riseDivFactorialsCache[yIndex][kIndex] = result;
            }
        }
    }

    /**
     * @param nowActiveLineages now active lineages
     * @param k                 number of steps
     * @return rise-div-factorials term of coalescent formula
     */
    private BigDecimal getRiseDivFactorialsTerm(int nowActiveLineages, int k) {
        return riseDivFactorialsCache[nowActiveLineages - 2][k - nowActiveLineages];
    }


    /* -- PAIR PICKING FORMULA -- */

    /** Cache for pair picking formula; use method to access values. */
    private BigDecimal[][][][] pairPickingCache;

    /**
     * Initialize the cache for pair picking values with appropriate sizes.
     *
     * @param largeCladeSize size of largest large child clade
     * @param smallCladeSize size of largest small child clade
     */
    private void initializePairPickingCache(int largeCladeSize, int smallCladeSize) {
        if ((largeCladeSize <= 0) || (smallCladeSize <= 0)) {
            throw new IllegalArgumentException(String.format("Cannot initialize for non-positive clade size %d.",
                    Math.min(largeCladeSize, smallCladeSize)));
        }

        // indexing of cache is [jPrev][kPrev][dj][dk]
        // observations:
        // - jPrev in [1, ..., largeCladeSize] => jPrev values, so index for jPrev shifted by -1
        // - kPrev in [0, ..., smallCladeSize] (note that even after swap)
        // -- but also kPrev <= jPrev, which we can assume this by pair picking formula symmetry
        // -- so kPrev in [0, ..., jPrev]
        // - dj in [0, ..., jPrev - 1] because at least one remains
        // - dk in [0, ..., kPrev - 1] because at least one remains

        // System.out.println("largeCladeSize = " + largeCladeSize);
        pairPickingCache = new BigDecimal[largeCladeSize][][][];
        for (int jIndex = 0; jIndex < pairPickingCache.length; jIndex++) {
            int jPrev = jIndex + 1;
            int kSize = Math.min(jPrev + 1, smallCladeSize + 1);
            pairPickingCache[jIndex] = new BigDecimal[kSize][][];

            for (int kPrev = 0; kPrev < kSize; kPrev++) {
                pairPickingCache[jIndex][kPrev] = new BigDecimal[jPrev][Math.max(kPrev, 1)];
            }
        }
    }

    /**
     * For the given parameters, get the number of orders valid coalescent events can occur
     * multiplied by their probability, which together we here call the PairPicking formula.
     *
     * @param jPrev number of left active lineages at last time step
     * @param kPrev number of right active lineages at last time step
     * @param j     number of left active lineages at current time step
     * @param k     number of right active lineages at current time step
     * @return result of pair picking formula for given parameters
     */
    private BigDecimal getPairPickingValue(int jPrev, int kPrev, int j, int k) {
        if (jPrev < kPrev) {
            int tmp = jPrev;
            jPrev = kPrev;
            kPrev = tmp;

            tmp = j;
            j = k;
            k = tmp;
        }

        int dj = jPrev - j;
        int dk = kPrev - k;

        // get value from cache, lazily. Only initialise if we need it.
        BigDecimal r = pairPickingCache[jPrev - 1][kPrev][dj][dk];
        if (r == null) {
            final int prev = jPrev + kPrev;
            r = BigDecimal.ONE;
            for (int m = 0; m < dj; m++) {
                BigDecimal numerator = new BigDecimal((jPrev - m) * (jPrev - m - 1));
                BigDecimal denonimator = new BigDecimal((prev - m) * (prev - m - 1));
                BigDecimal term = numerator.divide(denonimator, MC);
                r = r.multiply(term, MC);
                // r *= (dj + dk - m)/(m + 1.0);
            }
            for (int n = 0; n < dk; n++) {
                BigDecimal numerator = new BigDecimal((kPrev - n) * (kPrev - n - 1.0));
                BigDecimal denonimator = new BigDecimal((prev - dj - n) * (prev - dj - n - 1));
                BigDecimal term = numerator.divide(denonimator, MC);
                r = r.multiply(term, MC);
                // r *= (dk - n)/(dj + n + 1.0);
            }
            final int n = dj + dk;
            if (n > BINOMIAL_COEFFICIENT_LONG_THRESHOLD) {
                r = r.multiply(new BigDecimal(getBinomialCoefficientDouble(n, dj)), MC);
            } else {
                r = r.multiply(new BigDecimal(getBinomialCoefficientLong(n, dj)), MC);
            }
            pairPickingCache[jPrev - 1][kPrev][dj][dk] = r;
        }

        // checkProbability(r);
        return r;
    }


    /* -- BINOMIAL COEFFICIENT CACHE -- */

    /** Cache for binomial coefficients up to a maximum size; use method to access values. */
    private static long[][] binomialCoefficientCache;

    /** Threshold for caching table and overflow concerns for binomial coefficient computation using long. */
    private static final int BINOMIAL_COEFFICIENT_LONG_THRESHOLD = 63;

    /** Threshold for binomial coefficient computation using double. */
    private static final int BINOMIAL_COEFFICIENT_DOUBLE_THRESHOLD = 1028;

    /** Shift between n and index in cache for binomial coefficient (n choose k). */
    private static final int N_SHIFT = 4;
    /** Shift between k and index in cache for binomial coefficient (n choose k). */
    private static final int K_SHIFT = 2;

    /**
     * Initialize the cache for binomial coefficients with given parameters limiting the required size.
     *
     * @param numTaxa        max n in (n choose k)
     * @param largeCladeSize max k in (n choose k)
     */
    private void initializeBinomialCoefficientCache(int numTaxa, int largeCladeSize) {
        // combination: Cache table n choose k. For n up to binomialCoefficientThreshold (included)
        // the following cases are handled explicitly
        //      (k choose k) = (k choose 0) = 1
        //      (k choose k-1) = (k choose 1) = k
        // so neither (2 choose 1) nor (3 choose 2) = (3 choose 1) are ever called.
        // hence index for n is shifted by 3+1 and index for k is shifted by 1+1
        // (+1 because indexing starts at 0)
        int maxBinomialCoefficient = Math.min(BINOMIAL_COEFFICIENT_LONG_THRESHOLD, numTaxa);
        binomialCoefficientCache = new long[maxBinomialCoefficient - N_SHIFT + 1][];

        for (int n = 4; n < binomialCoefficientCache.length + N_SHIFT; n++) {
            int nIndex = n - N_SHIFT;
            binomialCoefficientCache[nIndex] = new long[n / 2 + n % 2 + 1 - K_SHIFT];

            // start with k = 2
            binomialCoefficientCache[nIndex][0] = BinomialCoefficient.value(n, 2);

            for (int k = 3; (k - K_SHIFT < binomialCoefficientCache[nIndex].length) && (k <= largeCladeSize); k++) {
                int kIndex = k - K_SHIFT;

                // (n choose k) = (n choose k-1) * (n - k + 1) / k
                if (n <= 61) {
                    binomialCoefficientCache[nIndex][kIndex] =
                            binomialCoefficientCache[nIndex][kIndex - 1] * (n - k + 1) / k;
                } else {
                    // following apache BinomialCoefficient.value(n, k)
                    long upper = (n - k + 1);
                    long d = (long) ArithmeticUtils.gcd(upper, k);
                    binomialCoefficientCache[nIndex][kIndex] = binomialCoefficientCache[nIndex][kIndex - 1] / ((long) k / d) * (upper / d);
                }
            }
        }
    }

    /**
     * Get binomial coefficient (n choose k) within cache bounds.
     *
     * @param n of (n choose k)
     * @param k of (n choose k)
     * @return binomial coefficient (n choose k) if in cache bounds.
     */
    public static long getBinomialCoefficientLong(int n, int k) {
        if (n > BINOMIAL_COEFFICIENT_LONG_THRESHOLD) {
            throw new AssertionError(String.format("Should not use `long` binomial coefficient (n choose k) for n > %d;" +
                    "here called with n = %d", BINOMIAL_COEFFICIENT_LONG_THRESHOLD, n));
        }

        // by symmetry of binomial coefficient, can assume that k is at most half of n
        if (k > n / 2) {
            k = n - k;
        }

        // explicitly handled cases
        if (k == 0) { // also (k == n), but k <= n/2
            return 1;
        }
        if (k == 1) { // also (k == n - 1), but k <= n/2
            return n;
        }

        if (binomialCoefficientCache == null) {
            // in case method is called statically
            return BinomialCoefficient.value(n, k);
        }

        return binomialCoefficientCache[n - N_SHIFT][k - K_SHIFT];
    }

    /**
     * Get binomial coefficient (n choose k) within cache bounds.
     *
     * @param n of (n choose k)
     * @param k of (n choose k)
     * @return binomial coefficient (n choose k) if in cache bounds.
     */
    public static double getBinomialCoefficientDouble(int n, int k) {
        if (n > BINOMIAL_COEFFICIENT_DOUBLE_THRESHOLD) {
            throw new AssertionError(String.format("n = %d is too large for binomial coefficient computation, threshold is %d.",
                    n, BINOMIAL_COEFFICIENT_DOUBLE_THRESHOLD));
        }

        // by symmetry of binomial coefficient, can assume that k is at most half of n
        if (k > n / 2) {
            k = n - k;
        }

        double coefficient = 1.0;
        double upper = n - k + 1;
        for (int lower = 1; lower <= k; lower++) {
            coefficient *= upper / lower;
            upper++;
        }
        return coefficient;
    }


    /* -- SANITY CHECKS  & UTIL -- */

    /**
     * Checks whether given probability is non-negative
     * and throws an exception otherwise.
     *
     * @param p probability to check being non-negative
     */
    private static BigDecimal checkProbability(BigDecimal p) {
        if (p.compareTo(BigDecimal.ZERO) < 0) {
            if (p.compareTo(ROUNDING_EPSILON) >= 0) {
                return BigDecimal.ZERO;
            } else {
                System.err.println("Negative p: " + p);
                // return BigDecimal.ZERO;
                throw new IllegalArgumentException(String.format("Probability cannot be negative, but is %f!", p));
            }
        }

        if (p.compareTo(BigDecimal.ONE) > 0) {
            // System.err.println("too large p: " + p);
            // return BigDecimal.ONE;
            throw new IllegalArgumentException(String.format("too large p: %f!", p));
        }

        return p;
    }

    /**
     * Taylor series expansion of e^x
     *
     * @param x
     * @return
     */
    public static BigDecimal exp(BigDecimal x) {
        BigDecimal result = BigDecimal.ONE;
        BigDecimal term = BigDecimal.ONE;
        BigDecimal i = BigDecimal.ONE;
        BigDecimal tolerance = BigDecimal.ONE.scaleByPowerOfTen(-MC.getPrecision());

        while (term.abs().compareTo(tolerance) > 0) {
            term = term.multiply(x, MC).divide(i, MC);
            result = result.add(term, MC);
            i = i.add(BigDecimal.ONE);
        }

        return result;
    }

    // small test case
    public static void main(String[] args) throws IOException {

        Random random = new Random(123);
        // Random random = new Random();

        int numOfTaxa = 10;
        int largeCladeSize = (numOfTaxa / 2);
        if (numOfTaxa % 2 == 1) {
            largeCladeSize = (numOfTaxa / 2) + 1;
        }
        PriorCoalescentCCPComputer testCase = new PriorCoalescentCCPComputer(numOfTaxa, largeCladeSize, (numOfTaxa - largeCladeSize));

        double epsilon = 1e-6;
        double theta = 1.0;

        List<Double> timeList = new ArrayList<>(); // all taxa times
        double[] timesLeft = new double[largeCladeSize];
        double[] timesRight = new double[numOfTaxa - largeCladeSize];

        for (int i = 1; i <= numOfTaxa; i++) { // create all taxa
            double time = (double) (i - 1) / (numOfTaxa - 1); // total time = 1
            timeList.add(time);
        }

        for (int i = 0; i < largeCladeSize; i++) { // randomly pick some times to be in the left clade
            int randomIndex = random.nextInt(timeList.size());
            double randomTime = timeList.remove(randomIndex);
            timesLeft[i] = randomTime;
        }

        for (int i = 0; i < timeList.size(); i++) { // put the rest of the times in the right clade
            double rightTime = timeList.get(i);
            timesRight[i] = rightTime;
        }

        long startingExact = System.currentTimeMillis();
        double resultExact = testCase.getPriorCCPExact(theta, timesLeft, timesRight);
        long endingExact = System.currentTimeMillis();
        long timeExact = (endingExact - startingExact);
        System.out.println("result Exact = " + resultExact);
        System.out.println("time Exact = " + timeExact);

        long startingHeuristic = System.currentTimeMillis();
        double resultHeuristic = testCase.getPriorCCPHeuristically(theta, timesLeft, timesRight, epsilon);
        long endingHeuristic = System.currentTimeMillis();
        long timeHeuristic = (endingHeuristic - startingHeuristic);
        System.out.println("result Heuristic = " + resultHeuristic);
        System.out.println("time Heuristic = " + timeHeuristic);
    }

}
