package ccd.model;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator.TreeSet;
import ccd.model.bitsets.BitSet;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * CCD0 normalisation on the extended-split (SJ) clade identity for
 * sampled-ancestor trees. The CCD0-normalised companion to {@link CCD1SJ}:
 * same identity, but partition probabilities are derived from clade
 * frequencies via the expand-and-normalise step (CCD0), instead of from
 * observed split counts alone (CCD1).
 *
 * <p>Identity construction, key layout, sampling, MAP selection, and tree
 * probability traversal are shared with {@link CCD1SJ} via
 * {@link SJSupport}; see that class for the SJ identity definition. What's
 * specific to CCD0SJ:
 *
 * <h3>Tree probability</h3>
 *
 * <p>{@code P(T)} is the product of normalised CCPs along the tree path
 * (under the SJ identity), folding in the placeholder root's
 * {@code (variant, emptyClade)} CCP at the top. CCPs come from CCD0's
 * expand-and-normalise step run on the SJ extended-clade graph: for each
 * clade {@code C}, the CCP of partition {@code (A, B)} is proportional to
 * the product of subtree clade-credibility masses of {@code A} and
 * {@code B}, normalised across all partitions on {@code C}. With expansion
 * (see {@link #expand()}), the support extends beyond the observed sample
 * to trees the model considers consistent with the observed clade
 * frequencies and the SA-flag constraints.
 *
 * <h3>Empty sister as a neutral element</h3>
 *
 * <p>Under CCD0 normalisation each clade contributes a subtree mass during
 * recursion. The shared empty sister has no taxa and no partitions, so we
 * pre-set its sumCladeCredibilities to 1 in {@link #initializeRootClade}
 * and short-circuit its handling in
 * {@link #setPartitionProbabilities(Clade, boolean)}. This makes a
 * placeholder-root partition {@code (variant, empty)} score equal to the
 * variant's own subtree mass, so root variants compete cleanly under one
 * normalisation.
 *
 * <h3>Expand step</h3>
 *
 * <p>The CCD0 expand step is overridden ({@link #expand()}) for two
 * reasons: (1) it must operate on taxa-only bitmasks for disjointness/union
 * checks, since the SA flags otherwise spill out of the size-n bucket
 * arrays in {@link CCD0#expand()}, and (2) it enforces SJ-flag
 * compatibility on registered partitions
 * ({@link SJSupport#isSJCompatible}), so the model never assigns
 * probability to trees that violate the parent's SA-flag identity.
 */
public class CCD0SJ extends CCD0 {

    /**
     * Shared empty sister used in rootClade's (variant, emptyClade) partitions.
     * Populated in {@link #initializeRootClade} since the super constructor
     * calls cladifyTree before subclass field initialisers run.
     */
    protected Clade emptyClade;

    public CCD0SJ(List<Tree> trees, double burnin) {
        super(trees, burnin);
    }

    public CCD0SJ(TreeSet treeSet) {
        super(treeSet);
    }

    public CCD0SJ(TreeSet treeSet, int numTreesToUse) {
        super(treeSet, numTreesToUse);
    }

    public CCD0SJ(TreeSet treeSet, boolean storeBaseTrees) {
        super(treeSet, storeBaseTrees);
    }

    public CCD0SJ(int numLeaves, boolean storeBaseTrees) {
        super(numLeaves, storeBaseTrees);
    }

    /* -- ROOT INITIALISATION -- */

    @Override
    protected void initializeRootClade(int numLeaves) {
        this.emptyClade = SJSupport.initializeSJRoot(this, numLeaves);
        // CCD0's recursive setPartitionProbabilities short-circuits on a
        // positive sumCladeCredibilities, so seed the empty sister with the
        // neutral value 1.0 before the first call.
        this.emptyClade.setSumCladeCredibilities(1.0);
    }

    /* -- TREE INSERTION -- */

    @Override
    protected Clade cladifyVertex(Node vertex) {
        return SJSupport.cladifyTreeRoot(this, emptyClade, vertex);
    }

    /* -- SA DETECTION & DISPLAY -- */

    @Override
    public boolean isSampledAncestor(Clade clade) {
        return SJSupport.isSampledAncestor(clade, leafArraySize);
    }

    @Override
    protected String getSampledAncestorInfoString(Clade clade) {
        return SJSupport.saInfoString(clade, leafArraySize);
    }

    /* -- EXPAND -- */

    /**
     * Adds clade partitions where parent and children clades were observed
     * but the partition itself was not. Operates on taxa-only bitmasks for
     * disjointness/union checks (so the SA-augmented bitset width does not
     * break the size-n bucket arrays in {@link CCD0#expand()}). A discovered
     * partition (A, B) is registered on every SA-variant of the parent
     * clade sharing that base mask whose SA flags are compatible with the
     * partition (see {@link SJSupport#isSJCompatible}).
     *
     * <p>The placeholder {@code rootClade} (size 0, sentinel bit set) and
     * the shared {@code emptyClade} are skipped: the placeholder root only
     * ever carries (variant, emptyClade) partitions added during
     * cladification, and normalisation across variants is handled by
     * {@link CCD0#setPartitionProbabilities(Clade, boolean)}.
     */
    @Override
    protected void expand() {
        // Collect non-trivial real clades (skip placeholder root, empty,
        // leaves, cherries — the latter two have no expandable partitions).
        List<Clade> realClades = new ArrayList<>();
        Map<BitSet, List<Clade>> taxaMaskToClades = new HashMap<>();
        for (Clade c : cladeMapping.values()) {
            if (c == rootClade || c == emptyClade) continue;
            realClades.add(c);
            BitSet taxaMask = c.getCladeInBitsTaxaOnly();
            taxaMaskToClades.computeIfAbsent(taxaMask, k -> new ArrayList<>()).add(c);
        }

        int n = realClades.size();
        for (int i = 0; i < n; i++) {
            Clade A = realClades.get(i);
            BitSet taxaA = A.getCladeInBitsTaxaOnly();
            for (int j = i + 1; j < n; j++) {
                Clade B = realClades.get(j);
                BitSet taxaB = B.getCladeInBitsTaxaOnly();

                // Taxa-disjointness check via base mask.
                if (taxaA.intersects(taxaB)) continue;

                // Compute taxa union and look up all parent clades sharing
                // that taxa-mask (any number of SA variants).
                BitSet unionMask = (BitSet) taxaA.clone();
                unionMask.or(taxaB);

                List<Clade> parents = taxaMaskToClades.get(unionMask);
                if (parents == null) continue;

                for (Clade parent : parents) {
                    if (parent == A || parent == B) continue;
                    if (parent.size() <= 2) continue; // cherries already complete
                    if (parent.getCladePartition(A, B) != null) continue;
                    if (!SJSupport.isSJCompatible(parent, A, B, leafArraySize)) continue;
                    parent.createCladePartition(A, B);
                }
            }
        }
    }

    /* -- PROBABILITY (CCD0 normalisation) -- */

    /**
     * Short-circuits the empty sister to a neutral 1.0 so CCD0's recursion
     * does not stumble over a clade with size 0 and no partitions; otherwise
     * delegates to CCD0's frequency-based normalisation.
     */
    @Override
    public double setPartitionProbabilities(Clade clade, boolean useCladeParameters) {
        if (clade == emptyClade) {
            clade.setSumCladeCredibilities(1.0);
            return 1.0;
        }
        return super.setPartitionProbabilities(clade, useCladeParameters);
    }

    /* -- TREE PROBABILITY -- */

    @Override
    protected Clade computeProbabilityOfVertex(Node vertex, double[] runningProbability, boolean computeLog) {
        return SJSupport.computeTreeProbability(this, emptyClade, vertex, runningProbability, computeLog);
    }

    /* -- SAMPLING & MAP -- */

    /**
     * Builds a tree under the SJ identity. The {@code heightStrategy}
     * argument is currently ignored — see
     * {@link SJSupport#buildTree(AbstractCCD, Clade, SamplingStrategy, HeightSettingStrategy)}
     * for the TODO and the reason.
     */
    @Override
    protected Tree getTreeBasedOnStrategy(SamplingStrategy samplingStrategy,
                                          HeightSettingStrategy heightStrategy) {
        return SJSupport.buildTree(this, emptyClade, samplingStrategy, heightStrategy);
    }
    
    /* -- MISC -- */

    @Override
    public String toString() {
        return "CCD0-SJ " + super.toString().replaceFirst("CCD0 ", "");
    }

    @Override
    public AbstractCCD copy() {
        CCD0SJ copy = new CCD0SJ(this.getSizeOfLeavesArray(), false);
        copy.baseTrees.add(this.getSomeBaseTree());
        copy.numBaseTrees = this.getNumberOfBaseTrees();
        AbstractCCD.buildCopy(this, copy);
        // buildCopy overwrites cladeMapping entries for the placeholder root
        // and emptyClade with fresh copies, so re-resolve the field reference.
        copy.emptyClade = copy.cladeMapping.get(this.emptyClade.getCladeInBits());
        copy.emptyClade.setSumCladeCredibilities(1.0);
        return copy;
    }
}
