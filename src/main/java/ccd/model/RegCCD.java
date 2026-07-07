package ccd.model;

import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator;
import ccd.algorithms.regularisation.CCDExpansion;
import ccd.algorithms.regularisation.CCD1Regularisor;
import ccd.algorithms.regularisation.CCDRegularisationStrategy;
import java.util.List;

public class RegCCD extends CCD1 {

    /**
     * Constructor for an empty CCD. Trees can then be processed one by one.
     *
     * @param numLeaves      number of leaves of the trees that this CCD will be based on
     */
    public RegCCD(int numLeaves) {
        super(numLeaves, true);
    }

    /**
     * Constructor for a {@link RegCCD} based on the given collection of trees
     * with specified burn-in.
     *
     * @param trees  the trees whose distribution is approximated by the resulting
     *               {@link RegCCD}
     * @param burnin value between 0 and 1 of what percentage of the given trees
     *               should be discarded as burn-in
     */
    public RegCCD(List<Tree> trees, double burnin) {
        super(trees, burnin);
        postConstruction();
    }

    /**
     * Constructor for a {@link RegCCD} based on the given collection of trees
     * (not containing any burnin trees).
     *
     * @param treeSet        an iterable set of trees, which contains no burnin trees,
     *                       whose distribution is approximated by the resulting
     *                       {@link RegCCD}; all of its trees are used
     */
    public RegCCD(TreeAnnotator.TreeSet treeSet) {
        super(treeSet, true);
        postConstruction();
    }

    protected void postConstruction() {
        expandRegCCD();
        regulariseRegCCD();
        System.out.println("Regularisation on CCD with default parameter = 0.4 is done");
    }

    @Override
    public String toString() {
        return "Regularised " + super.toString();
    }

    public void expandRegCCD() {
        CCDExpansion expansionStrategy = new CCDExpansion();
        expansionStrategy.expandCCD(this);
    }

    public void regulariseRegCCD() {
        CCD1Regularisor regulariser = new CCD1Regularisor(CCDRegularisationStrategy.AdditiveX, 0.4);
        regulariser.regularise(this);
    }
}
