package ccd.tools;

import beast.base.core.Description;
import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator;
import beastfx.app.treeannotator.TreeAnnotator.TreeSet;
import beastfx.app.treeannotator.services.TopologySettingService;
import ccd.model.AbstractCCD;
import ccd.model.CCD2;
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
        return "CCD1";
    }

    @Override
    public String getDescription() {
        return "MAP (CCD2)";
    }

    @Override
    protected boolean sanityCheck(Tree mapTree, Tree firstTree, AbstractCCD ccd) {
        try {
            return singlySupportedCheck(mapTree, ccd)
                    && equalsFirstTreeCheck(mapTree, firstTree, ccd);
        } catch (Throwable e) {
            // we do not want sanity checks to ruin the job,
            // so report any issues but otherwise ignore it.
            e.printStackTrace();
            return false;
        }
    }

}
