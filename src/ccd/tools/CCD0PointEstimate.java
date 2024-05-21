package ccd.tools;

import beast.base.core.Description;
import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator;
import beastfx.app.treeannotator.TreeAnnotator.TreeSet;
import beastfx.app.treeannotator.services.TopologySettingService;
import ccd.model.AbstractCCD;
import ccd.model.CCD0;
import ccd.model.HeightSettingStrategy;

import java.io.IOException;
import java.io.PrintStream;

@Description("TreeAnnotator plugin for setting the tree topology as CCD1 MAP tree")
public class CCD0PointEstimate extends PointEstimate implements TopologySettingService {

    @Override
    public Tree setTopology(TreeSet treeSet, PrintStream progressStream, TreeAnnotator annotator)
            throws IOException {

        progressStream.println("CCD0 MAP tree computation");
        progressStream.println("0              25             50             75            100");
        progressStream.println("|--------------|--------------|--------------|--------------|");

        treeSet.reset();
        Tree tree = treeSet.next();
        Tree firstTree = tree;
        CCD0 ccd = new CCD0(tree.getLeafNodeCount(), false);
        ccd.setProgressStream(progressStream);

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
        ccd.initialize();

        Tree mapTree = ccd.getMAPTree(HeightSettingStrategy.One);
        sanityCheck(mapTree, firstTree, ccd);

        return mapTree;
    }

    @Override
    public String getServiceName() {
        return "CCD0";
    }

    @Override
    public String getDescription() {
        return "MAP (CCD0)";
    }

    @Override
    protected boolean sanityCheck(Tree mapTree, Tree firstTree, AbstractCCD ccd) {
        // nothing to do for CCD0 MAP tree
        // could analyse singly supported clades,
        // but that would mean we do post-analysis work of the users
        return true;
    }
}
