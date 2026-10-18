package ccd.model;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator.TreeSet;

import java.util.List;

/**
 * CCD1 with the extended-split (EI) clade identity for sampled-ancestor trees.
 *
 * <p>EI is the tightest of the four CCD1 SA variants formalised in
 * {@code CCD_sampled_ancestor/algorithms/alexei_c_four_sa_models.pdf}
 * (Table 1, row 4): each internal clade carries SA flags for its direct leaf
 * children only (no propagation), and tree probability is a product over
 * conditional split probabilities under the full extended-clade identity.
 * EI requires every extended split to have been observed as a whole, so
 * its support is the sampled topologies and SA patterns only — no expand
 * step. This is the model {@code ccd-js} implements.
 *
 * <p>The placeholder {@code rootClade} carries one
 * {@code (variant, emptyClade)} partition per observed root variant, with
 * marginal probability {@code variantCount / numBaseTrees}. Identity bitset
 * layout, key construction, sampling, and tree-probability traversal are all
 * shared with {@link CCD0EI} via {@link EISupport}.
 *
 * @see CCD0EI for the CCD0-normalised companion.
 */
public class CCD1EI extends CCD1 {

    /**
     * Shared empty sister used in rootClade's (variant, emptyClade) partitions.
     * Populated in {@link #initializeRootClade} since the super constructor
     * calls cladifyTree before subclass field initialisers run.
     */
    protected Clade emptyClade;

    public CCD1EI(List<Tree> trees, double burnin) {
        super(trees, burnin);
    }

    public CCD1EI(TreeSet treeSet) {
        super(treeSet);
    }

    public CCD1EI(TreeSet treeSet, int numTreesToUse) {
        super(treeSet, numTreesToUse);
    }

    public CCD1EI(TreeSet treeSet, boolean storeBaseTrees) {
        super(treeSet, storeBaseTrees);
    }

    public CCD1EI(int numLeaves, boolean storeBaseTrees) {
        super(numLeaves, storeBaseTrees);
    }

    /* -- ROOT INITIALISATION -- */

    @Override
    protected void initializeRootClade(int numLeaves) {
        this.emptyClade = EISupport.initializeEIRoot(this, numLeaves);
    }

    /* -- TREE INSERTION -- */

    @Override
    protected Clade cladifyVertex(Node vertex) {
        return EISupport.cladifyTreeRoot(this, emptyClade, vertex);
    }

    /* -- SA DETECTION & DISPLAY -- */

    @Override
    public boolean isSampledAncestor(Clade clade) {
        return EISupport.isSampledAncestor(clade, leafArraySize);
    }

    @Override
    protected String getSampledAncestorInfoString(Clade clade) {
        return EISupport.saInfoString(clade, leafArraySize);
    }

    /* -- TREE PROBABILITY -- */

    @Override
    protected Clade computeProbabilityOfVertex(Node vertex, double[] runningProbability, boolean computeLog) {
        return EISupport.computeTreeProbability(this, emptyClade, vertex, runningProbability, computeLog);
    }

    /* -- HELD-OUT (LEAVE-ONE-TREE-OUT) PROBABILITY -- */
    // CCD1's inherited held-out traversal assumes the CCD1/CCD1EL leaf-SA-bit
    // encoding and fails on EI clade identities; route through EISupport instead.

    @Override
    public double getProbOfHeldOutTree(Tree tree, double alpha) {
        resetCacheIfProbabilitiesDirty();
        double[] runningProbability = new double[]{1};
        EISupport.computeTreeProbabilityHeldOut(this, emptyClade, tree.getRoot(), runningProbability, alpha, false);
        return runningProbability[0];
    }

    @Override
    public double getLogProbOfHeldOutTree(Tree tree, double alpha) {
        resetCacheIfProbabilitiesDirty();
        double[] runningProbability = new double[]{0};
        EISupport.computeTreeProbabilityHeldOut(this, emptyClade, tree.getRoot(), runningProbability, alpha, true);
        return runningProbability[0];
    }

    /* -- SAMPLING & MAP -- */

    @Override
    protected Tree getTreeBasedOnStrategy(SamplingStrategy samplingStrategy,
                                          HeightSettingStrategy heightStrategy) {
        return EISupport.buildTree(this, emptyClade, samplingStrategy, heightStrategy);
    }

    /* -- MISC -- */

    @Override
    public String toString() {
        return "CCD1-EI " + super.toString().replaceFirst("CCD1 ", "");
    }

    @Override
    public AbstractCCD copy() {
        CCD1EI copy = new CCD1EI(this.getSizeOfLeavesArray(), false);
        copy.baseTrees.add(this.getSomeBaseTree());
        copy.numBaseTrees = this.getNumberOfBaseTrees();
        AbstractCCD.buildCopy(this, copy);
        // buildCopy overwrites cladeMapping entries for the placeholder root
        // and emptyClade with fresh copies, so re-resolve the field reference.
        copy.emptyClade = copy.cladeMapping.get(this.emptyClade.getCladeInBits());
        return copy;
    }
}
