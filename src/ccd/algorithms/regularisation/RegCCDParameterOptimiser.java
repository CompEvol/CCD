package ccd.algorithms.regularisation;

import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator;
import ccd.model.RegCCD;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.SearchInterval;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;
import org.apache.commons.math3.optim.univariate.UnivariatePointValuePair;

import java.io.IOException;

public class RegCCDParameterOptimiser {

    public static UnivariateFunction defineFunction(RegCCD regCCD, TreeAnnotator.TreeSet treeSet) {

        UnivariateFunction f = alpha -> {

            double totalLogProb = 0.0;
            int nonNegInfLogProbTreeCount = 0;

            try {
                treeSet.reset();
                while (treeSet.hasNext()) {
                    Tree currentTree = treeSet.next();
                    double prob = regCCD.getProbOfHeldOutTree(currentTree, alpha);
                    double logProb = regCCD.getLogProbOfHeldOutTree(currentTree, alpha);
                    if (logProb != Double.NEGATIVE_INFINITY) {
                        totalLogProb += logProb;
                        nonNegInfLogProbTreeCount++;
                    }
                    if ((logProb != Double.NEGATIVE_INFINITY && prob == 0) ||
                            (logProb == Double.NEGATIVE_INFINITY && prob != 0)) {
                        throw new IllegalStateException(String.format("Inconsistent probability detected: prob = %g, logProb = %g", prob, logProb));
                    }
                }
                double avgLogProb = totalLogProb / nonNegInfLogProbTreeCount;
            } catch (IOException e) {
                throw new RuntimeException(e);
            }

            double avgLogProb = totalLogProb / nonNegInfLogProbTreeCount;
            System.out.println(String.format("testing alpha = %.5f, avgLogProb = %.5f", alpha, avgLogProb));
            return avgLogProb;
        };
        return f;
    }

    public static double optimise(UnivariateFunction f, double relTol, double absTol, int maxEval, double startingPoint) {
        final double lowerBound = 0.0, upperBound = 2;
        BrentOptimizer opt = new BrentOptimizer(relTol, absTol);
        UnivariatePointValuePair sol = opt.optimize(
                new MaxEval(maxEval),
                new UnivariateObjectiveFunction(f),
                GoalType.MAXIMIZE,
                new SearchInterval(lowerBound, upperBound, startingPoint)
        );
        return sol.getPoint();
    }
}
