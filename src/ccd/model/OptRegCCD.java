package ccd.model;

import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator;
import ccd.algorithms.regularisation.CCD1Regularisor;
import ccd.algorithms.regularisation.CCDRegularisationStrategy;
import ccd.algorithms.regularisation.RegCCDParameterOptimiser;
import org.apache.commons.math3.analysis.UnivariateFunction;

import java.util.List;

public class OptRegCCD extends RegCCD {

    /**
     * Constructor for an empty CCD. Trees can then be processed one by one.
     *
     * @param numLeaves      number of leaves of the trees that this CCD will be based on
     */
    public OptRegCCD(int numLeaves) {
        super(numLeaves);
    }

    /**
     * Constructor for a {@link OptRegCCD} based on the given collection of trees
     * with specified burn-in.
     *
     * @param trees  the trees whose distribution is approximated by the resulting
     *               {@link OptRegCCD}
     * @param burnin value between 0 and 1 of what percentage of the given trees
     *               should be discarded as burn-in
     */
    public OptRegCCD(List<Tree> trees, double burnin) {
        super(trees, burnin);
    }

    /**
     * Constructor for a {@link OptRegCCD} based on the given collection of trees
     * (not containing any burnin trees).
     *
     * @param treeSet        an iterable set of trees, which contains no burnin trees,
     *                       whose distribution is approximated by the resulting
     *                       {@link OptRegCCD}; all of its trees are used
     */
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
