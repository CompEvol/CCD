package ccd.tools;

import beast.base.core.Log;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeUtils;
import ccd.model.AbstractCCD;
import ccd.model.BitSet;
import ccd.model.CCD2;
import ccd.model.Clade;

import java.text.DecimalFormat;

/**
 * This is an abstract base class for CCD0 and CCD1 point estimates, sharing sanity checks.
 */
public abstract class PointEstimate {

    public static final String WARNING_LINE = "- - - - - WARNING WARNING - - - - - WARNING WARNING - - - - - ";

    /**
     * Check whether the given computed MAP tree of a CCD passes some simple sanity check.
     *
     * @param mapTree   tree that gets sanity checked here
     * @param firstTree the first tree of the tree set used to construct the CCD
     * @param ccd       from which the mapTree is from
     * @return whether the MAP tree passes the sanity check
     */
    abstract protected boolean sanityCheck(Tree mapTree, Tree firstTree, AbstractCCD ccd);

    /**
     * A sanity check producing warnings if the first tree equals the MAP tree.
     *
     * @param mapTree   tree that gets sanity checked here
     * @param firstTree the first tree of the tree set used to construct the CCD
     * @param ccd       from which the mapTree is from
     * @return whether the MAP tree passes the sanity check
     */
    public static boolean equalsFirstTreeCheck(Tree mapTree, Tree firstTree, AbstractCCD ccd, boolean verbose) {
        if (ccd.getProbabilityOfTree(mapTree) == 1) {
            // there is only one topology, so of course MAP-tree would equal the first tree
            // in fact, a sanity check is not even required
            return true;
        }

        String newick1 = TreeUtils.sortedNewickTopology(mapTree.getRoot(), false).trim();
        String newick2 = TreeUtils.sortedNewickTopology(firstTree.getRoot(), false).trim();

        if (newick1.equals(newick2)) {
            if (verbose) {
                Log.warning(WARNING_LINE + "\n");
                Log.warning("The summary tree equals the first tree in the tree set!");
                Log.warning("Unless you expect to have very informative data,");
                Log.warning("this strongly suggests burn-in was not removed, and the summary tree is not valid.");
                Log.warning("In this case, we recommend analyzing the necessary burn-in and rerunning TreeAnnotator.");
                Log.warning("If you removed enough burn-in and this warning keeps showing,");
                Log.warning("then this might be just a coincidence or simply reflect your data.");
                Log.warning("\n" + WARNING_LINE);
            }
            return false;
        }
        return true;
    }

    /**
     * A sanity check producing warnings if there are any clades in the MAP tree
     * that are only supported by a single tree.
     *
     * @param mapTree tree that gets sanity checked here
     * @param ccd     from which the mapTree is from
     * @return whether the MAP tree passes this sanity check
     */
    public static boolean singlySupportedCheck(Tree mapTree, AbstractCCD ccd, boolean verbose) {
        /*- Sanity check: We count how many of the clades in the mapTree
        only have a single observation in the tree set.
        The tree may have been under monophyly constraints, so we also count how
        many clades are under 100% support. */
        int[] counts = new int[2];
        if (ccd instanceof CCD2) {
            CCD2PointEstimate.traverse(mapTree.getRoot().getLeft(), mapTree.getRoot().getRight(), counts, (CCD2) ccd);
        } else {
            traverse(mapTree.getRoot(), counts, ccd);
        }
        int singlySupported = counts[0];
        int monophyletics = counts[1];
        int nodeCount = mapTree.getInternalNodeCount();

        if (singlySupported >= 1) {
            double percentage = singlySupported / (double) (nodeCount - monophyletics) * 100;
            if (verbose) {
                Log.warning(WARNING_LINE + "\n");
                Log.warning("There are " + singlySupported + " clades supported by only a single tree in the tree set;");
                DecimalFormat df = new DecimalFormat("0.##");
                Log.warning("more precisely, " + df.format(percentage) + "% of non-monophyletic non-trivial clades sare supported by only a single tree.");
                Log.warning("Unless you expect to have very noisy data, this could be due to");
                Log.warning("not having removed enough burn-in as CCD1 and CCD2 are sensitive to that.");
                Log.warning("In this case, we recommend analyzing the necessary burn-in and rerunning TreeAnnotator.");
                Log.warning("If you removed enough burn-in and this warning keeps showing,");
                Log.warning("then this might be just a coincidence or simply reflect your data.");
                Log.warning("\n" + WARNING_LINE);
            }
            return false;
        }
        return true;
    }

    /* Recursive helper method to count clades observed only once and monophyletic clades */
    protected static BitSet traverse(Node node, int[] counts, AbstractCCD ccd) {
        if (node.isLeaf()) {
            BitSet bitset = BitSet.newBitSet(ccd.getSizeOfLeavesArray());
            bitset.set(node.getNr());
            return bitset;
        } else {
            BitSet bitset = traverse(node.getLeft(), counts, ccd);
            BitSet bitset2 = traverse(node.getRight(), counts, ccd);
            bitset.or(bitset2);

            Clade clade = ccd.getClade(bitset);
            int count = clade.getNumberOfOccurrences();
            if (count == 1) {
                counts[0]++;
            }
            if (clade.isMonophyletic()) {
                counts[1]++;
            }
            return bitset;
        }

    }
}
