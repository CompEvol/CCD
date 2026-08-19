package ccd.model;

import beast.base.evolution.tree.Tree;
import ccd.algorithms.credibleSets.CredibleSetType;
import ccd.model.bitsets.BitSet;

import java.math.BigInteger;
import java.util.Collection;

public interface ITreeDistribution {

    /**
     * Randomly samples a tree (without heights set) from this distribution.
     *
     * @return a tree sampled from this distribution
     */
    public Tree sampleTree();

    /**
     * Randomly samples a tree from this distribution with heights set with the
     * given strategy.
     *
     * @param heightStrategy the strategy used to set the heights of the tree vertices
     * @return a tree sampled from this distribution
     */
    public Tree sampleTree(HeightSettingStrategy heightStrategy);

    /**
     * Return the probability of a randomly sampled tree from this distribution.
     *
     * @return probability of randomly sampled tree
     */
    public double sampleTreeProbability();

    public double sampleTreeLogProbability();

    /**
     * Returns the tree (without heights set) with maximum probability in this
     * distribution.
     *
     * @return the tree with maximum probability
     */
    public Tree getMAPTree();

    /**
     * Returns the tree with maximum probability in this distribution with
     * heights set with the given strategy.
     *
     * @param heightStrategy the strategy used to set the heights of the tree vertices
     * @return the tree with maximum probability
     */
    public Tree getMAPTree(HeightSettingStrategy heightStrategy);

    /** @return the maximum probability of any tree in this distribution */
    public double getMaxTreeProbability();

    /**
     * Return the probability of the given tree in this distribution.
     *
     * @param tree whose probability is requested
     * @return the probability of the given tree
     */
    public double getProbabilityOfTree(Tree tree);

    public double getLogProbabilityOfTree(Tree tree);

    /**
     * Returns whether this distribution contains the given tree.
     *
     * @param tree tested whether contained
     * @return whether this distribution contains the given tree
     */
    public boolean containsTree(Tree tree);

    /**
     * Returns the size of this distribution's <em>support</em>, that is, the number of distinct
     * tree topologies to which it assigns non-zero probability.
     *
     * <p>For a CCD whose support is exactly the set of topologies its graph represents (CCD0,
     * CCD1, regCCD) this equals the number of topologies of the graph. For a full-support model
     * (KRegCCD, MRegCCD, {@link UniformEscapeCCD}) it is the number of rooted topologies on the
     * taxon set, which is strictly larger: those models place probability outside their graph.
     * Implementations must report the support, not the graph, so that
     * {@code getNumberOfTrees()} and {@link #containsTree(Tree)} agree.
     *
     * <p>The result is a {@link BigInteger} because it overflows {@code long} at a handful of
     * taxa and {@code double} not long after: the number of rooted topologies needs 840 bits on
     * 129 taxa and 1093 bits on 160, and a {@code double} overflows past about 1024 bits. Take
     * logarithms with {@link AbstractCCD#logBigInteger(BigInteger)}, not via
     * {@code doubleValue()}, which is infinite from roughly 155 taxa upwards.
     *
     * @return the number of topologies with non-zero probability under this distribution
     */
    public BigInteger getNumberOfTrees();

    /**
     * Return the probability of the given clade.
     *
     * @param cladeInBits whose probability id requested
     * @return the probability of the given clade
     */
    public double getCladeProbability(BitSet cladeInBits);

    /** @return number of leaves/taxa of the trees in this distribution */
    public int getNumberOfLeaves();

    /** @return the number of distinct clades in this distribution */
    public int getNumberOfClades();

    /** @return all clades of this distribution */
    public Collection<Clade> getClades();

    /**
     * Return the min credible level of the given tree, that is,
     * the smallest credible set that still contains this tree.
     *
     * @param tree whose cred level is requested
     * @param type what type of information to use
     * @return min credible level of given tree
     */
    public double getCredibleLevel(Tree tree, CredibleSetType type);

}
