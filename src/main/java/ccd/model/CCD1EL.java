package ccd.model;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator;
import ccd.model.bitsets.BitSet;

import java.util.List;
import java.util.Map;

public class CCD1EL extends CCD1 {

    public CCD1EL(List<Tree> trees, double burnin) {
        super(trees, burnin);
    }

    public CCD1EL(TreeAnnotator.TreeSet treeSet) {
        super(treeSet);
    }

    public CCD1EL(TreeAnnotator.TreeSet treeSet, int numTreesToUse) {
        super(treeSet, numTreesToUse);
    }

    public CCD1EL(TreeAnnotator.TreeSet treeSet, boolean storeBaseTrees) {
        super(treeSet, storeBaseTrees);
    }

    public CCD1EL(int numLeaves, boolean storeBaseTrees) {
        super(numLeaves, storeBaseTrees);
    }

    /* -- ROOT INITIALIZATION -- */
    @Override
    protected void initializeRootClade(int numLeaves) {
        ELSupport.initializeRoot(this, numLeaves);
    }

    /* -- TREE INSERTION -- */
    @Override
    protected Clade cladifyVertex(Node vertex) {
        return ELSupport.cladifyVertex(this, vertex);
    }

    /**
     * @return whether this clade is marked as a sampled ancestor
     */
    @Override
    public boolean isSampledAncestor(Clade clade) {
        return ELSupport.isSampledAncestor(clade, leafArraySize);
    }

    /* -- SAMPLING & MAP -- */

    @Override
    protected double computeParentHeight(CladePartition partition, Node firstChild, Node secondChild) {
        return ELSupport.computeParentHeight(partition, firstChild, secondChild, leafArraySize);
    }

    /* -- HELPER METHODS -- */

    @Override
    protected Clade reduceCladeCount(Node vertex) {
        return ELSupport.reduceCladeCount(this, vertex, leafArraySize);
    }

    @Override
    protected BitSet computeCladeToNodeMapping(Node vertex, Map<Clade, Node> map) {
        return ELSupport.computeCladeToNodeMapping(vertex, map, this, leafArraySize);
    }

    /* -- TREE PROBABILITY -- */

    @Override
    protected Clade computeProbabilityOfVertex(Node vertex, double[] runningProbability, boolean computeLog) {
        return ELSupport.computeProbELVertex(vertex, runningProbability, computeLog, this, leafArraySize);
    }

    /* -- HELD-OUT (LEAVE-ONE-TREE-OUT) PROBABILITY -- */

    @Override
    public double getProbOfHeldOutTree(Tree tree, double alpha) {
        resetCacheIfProbabilitiesDirty();
        double[] runningProbability = new double[]{1};
        ELSupport.computeTreeProbabilityHeldOut(this, tree.getRoot(), runningProbability, alpha, false);
        return runningProbability[0];
    }

    @Override
    public double getLogProbOfHeldOutTree(Tree tree, double alpha) {
        resetCacheIfProbabilitiesDirty();
        double[] runningProbability = new double[]{0};
        ELSupport.computeTreeProbabilityHeldOut(this, tree.getRoot(), runningProbability, alpha, true);
        return runningProbability[0];
    }

    /* -- MISC -- */
    @Override
    public String toString() {
        return "CCD1-EL " + super.toString().replaceFirst("CCD1 ", "");
    }

    @Override
    protected String getSampledAncestorInfoString(Clade clade) {
        return ELSupport.saInfoString(clade, leafArraySize);
    }

    @Override
    public AbstractCCD copy() {
        CCD1EL copy = new CCD1EL(this.getSizeOfLeavesArray(), false);
        copy.baseTrees.add(this.getSomeBaseTree());
        copy.numBaseTrees = this.getNumberOfBaseTrees();
        AbstractCCD.buildCopy(this, copy);
        return copy;
    }
}