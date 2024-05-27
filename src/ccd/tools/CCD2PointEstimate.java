package ccd.tools;

import beast.base.core.Description;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator;
import beastfx.app.treeannotator.TreeAnnotator.TreeSet;
import beastfx.app.treeannotator.services.TopologySettingService;
import ccd.model.AbstractCCD;
import ccd.model.BitSet;
import ccd.model.CCD2;
import ccd.model.ExtendedClade;
import ccd.model.HeightSettingStrategy;

import java.io.IOException;
import java.io.PrintStream;


@Description("TreeAnnotator plugin for setting the tree topology as CCD1 MAP tree")
public class CCD2PointEstimate extends PointEstimate implements TopologySettingService {

    @Override
    public Tree setTopology(TreeSet treeSet, PrintStream progressStream, TreeAnnotator annotator)
            throws IOException {

        progressStream.println("CCD2 MAP tree computation");
        progressStream.println("0              25             50             75            100");
        progressStream.println("|--------------|--------------|--------------|--------------|");

        treeSet.reset();
        Tree tree = treeSet.next();
        Tree firstTree = tree;
        CCD2 ccd = new CCD2(tree.getLeafNodeCount(), false);

        int k = treeSet.totalTrees - treeSet.burninCount;
        int percentageDone = 0;
        int i = 0;
        while (tree != null) {
            // progress bar reporting
            while ((62 * i) / k > percentageDone) {
                progressStream.print("*");
                progressStream.flush();
                percentageDone++;
            }

            ccd.addTree(tree);
            tree = treeSet.hasNext() ? treeSet.next() : null;
            i++;
        }
        progressStream.println();

        Tree mapTree = ccd.getMAPTree(HeightSettingStrategy.One);
        sanityCheck(mapTree, firstTree, ccd);

        return mapTree;
    }

    @Override
    public String getServiceName() {
        return "CCD2";
    }

    @Override
    public String getDescription() {
        return "MAP (CCD2)";
    }

    @Override
    protected boolean sanityCheck(Tree mapTree, Tree firstTree, AbstractCCD ccd) {
        try {
            return singlySupportedCheck(mapTree, ccd, true)
                    && equalsFirstTreeCheck(mapTree, firstTree, ccd, true);
        } catch (Throwable e) {
            // we do not want sanity checks to ruin the job,
            // so report any issues but otherwise ignore it.
            e.printStackTrace();
            return false;
        }
    }

    /* Recursive helper method to count clades observed only once and monophyletic clades */
    protected static BitSet[] traverse(Node left, Node right, int[] counts, CCD2 ccd) {
        BitSet leftInBits = BitSet.newBitSet(ccd.getSizeOfLeavesArray());
        BitSet rightInBits = BitSet.newBitSet(ccd.getSizeOfLeavesArray());

        BitSet[] leftChildren = null;
        if (left.isLeaf()) {
            leftInBits.set(left.getNr());
        } else {
            leftChildren = traverse(left.getLeft(), left.getRight(), counts, ccd);
            leftInBits = leftChildren[0];
            leftInBits.or(leftChildren[1]);
        }
        BitSet[] rightChildren = null;
        if (right.isLeaf()) {
            rightInBits.set(right.getNr());
        } else {
            rightChildren = traverse(right.getLeft(), right.getRight(), counts, ccd);
            rightInBits = rightChildren[0];
            rightInBits.or(rightChildren[1]);
        }

        if (!left.isLeaf()) {
            ExtendedClade leftClade = ccd.getExtendedClade(leftInBits, rightInBits);
            int count = leftClade.getNumberOfOccurrences();
            if (count == 1) {
                counts[0]++;
            }
            if (leftClade.isMonophyletic()) {
                counts[1]++;
            }
        }
        if (!right.isLeaf()) {
            ExtendedClade rightClade = ccd.getExtendedClade(rightInBits, leftInBits);
            int count = rightClade.getNumberOfOccurrences();
            if (count == 1) {
                counts[0]++;
            }
            if (rightClade.isMonophyletic()) {
                counts[1]++;
            }
        }

        return new BitSet[]{leftInBits, rightInBits};
    }

}
