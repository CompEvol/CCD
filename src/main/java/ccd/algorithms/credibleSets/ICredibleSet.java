package ccd.algorithms.credibleSets;

import beast.base.evolution.tree.Tree;

/**
 * A credible set of an {@link ccd.model.ITreeDistribution} supporting (at least) the most important methods.
 *
 * @author Jonathan Klawitter
 */
public interface ICredibleSet {

    /**
     * Returns whether the given tree lies in the alpha credible set.
     *
     * @param tree  contained in the alpha credible set?
     * @param alpha the credible level in (0,1]
     * @return whether the given tree lies in the alpha credible set
     */
    default boolean isTreeInCredibleSet(Tree tree, double alpha) {
        return getCredibleLevel(tree) <= alpha;
    }

    /**
     * Compute and return the credible level of the given tree in the underlying {@link ccd.model.ITreeDistribution}.
     *
     * @param tree whose credible level is requested
     * @return credible level of the given tree in (0,1] or -1 if tree is not in distribution
     */
    double getCredibleLevel(Tree tree);

    /**
     * Sample a tree from the alpha credible set of the underlying  {@link ccd.model.ITreeDistribution}.
     *
     * @param alpha credible level in (0,1]
     * @return a tree sampled from the alpha
     */
    Tree sampleTreeFromCredibleSet(double alpha);

    /**
     * Converts and rounds the given credible level to an integer percentage between 0 and 100.
     * So a value x in [1, 100] represents being (already) in the x%-credible level
     * while a value 0, represents not being in the distribution at all.
     *
     * @param credibleLevel assumed in [0,1]
     * @return credible level converted to int as percentage
     */
    default int convertCredibleLevel(double credibleLevel) {
        if (credibleLevel <= 0) {
            return 0;
        } else {
            return (int) Math.ceil(credibleLevel * 100);
        }
    }
}
