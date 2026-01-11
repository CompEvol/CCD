package ccd.model;

import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator;
import ccd.algorithms.regularisation.CCD1Regularisor;
import ccd.algorithms.regularisation.CCDRegularisationStrategy;
import ccd.algorithms.regularisation.RegCCDParameterOptimiser;
import org.apache.commons.math3.analysis.UnivariateFunction;

import java.util.List;

public class OptRegCCD extends RegCCD {

    public OptRegCCD(int numLeaves) {
        super(numLeaves);
    }

    public OptRegCCD(List<Tree> trees, double burnin) {
        super(trees, burnin);
    }

    public OptRegCCD(TreeAnnotator.TreeSet treeSet) {
        super(treeSet);
    }

    @Override
    protected void postConstruction() {
        expandRegCCD();
        double optiAlpha = getOptimalAlpha();
        regulariseOptRegCCD(optiAlpha);
        System.out.println("Regularisation on CCD with optimised parameter = " + optiAlpha + " is done");
    }

    @Override
    public String toString() {
        return "Optimised " + super.toString();
    }

    public double getOptimalAlpha() {
        UnivariateFunction f = RegCCDParameterOptimiser.defineFunction(this, this.baseTrees);
        double threshold = 1e-6;
        double startingPoint;
        do { startingPoint = Math.random(); } while (startingPoint == 0.0);
        double optiAlpha = RegCCDParameterOptimiser.optimise(f, threshold, threshold, 200, startingPoint);
        return optiAlpha;
    }

    public void regulariseOptRegCCD(double alpha) {
        CCD1Regularisor regulariser = new CCD1Regularisor(CCDRegularisationStrategy.AdditiveX, alpha);
        regulariser.regularise(this);
    }
}
