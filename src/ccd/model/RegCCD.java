package ccd.model;

import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator;
import ccd.algorithms.regularisation.CCDExpansion;
import ccd.algorithms.regularisation.CCD1Regularisor;
import ccd.algorithms.regularisation.CCDRegularisationStrategy;

import java.util.List;

public class RegCCD extends CCD1 {

    public RegCCD(List<Tree> trees, double burnin) {
        super(trees, burnin);
        expandRegCCD();
        regulariseRegCCD();
    }

    public RegCCD(TreeAnnotator.TreeSet treeSet) {
        super(treeSet, true);
        expandRegCCD();
        regulariseRegCCD();
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
        CCD1Regularisor regulariser = new CCD1Regularisor(CCDRegularisationStrategy.AdditiveX, 0.5);
        regulariser.regularise(this);
    }
}
