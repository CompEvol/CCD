package ccd.tools;

import java.io.IOException;
import java.io.PrintStream;

import beast.base.core.Description;
import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator;
import beastfx.app.treeannotator.TreeAnnotator.TreeSet;
import beastfx.app.treeannotator.services.TopologySettingService;
import ccd.model.AbstractCCD;
import ccd.model.CCD0;
import ccd.model.HeightSettingStrategy;

@Description("TreeAnnotator plugin for setting the tree topology as CCD1 MAP tree")
public class CCD0ApproxPointEstimate extends PointEstimate implements TopologySettingService {

    @Override
    public Tree setTopology(TreeSet treeSet, PrintStream progressStream, TreeAnnotator annotator)
            throws IOException {

        progressStream.println("CCD0Approx MAP tree computation");
        progressStream.println("0              25             50             75            100");
        progressStream.println("|--------------|--------------|--------------|--------------|");

        treeSet.reset();
        Tree tree = treeSet.next();
        Tree firstTree = tree;
        int CCD0ApproxMultiplier = 50;
        if (System.getProperty("CCD0ApproxMultiplier") != null) {
        	try {
        		CCD0ApproxMultiplier = Integer.valueOf(System.getProperty("CCD0ApproxMultiplier"));
        	} catch (NumberFormatException e) {
        		progressStream.println("Could not parse CCD0ApproxMultiplier property (" + System.getProperty("CCD0ApproxMultiplier") + ").");
        		CCD0ApproxMultiplier = 50;
        		progressStream.println("Using default value of " + CCD0ApproxMultiplier + ".");
        	}
        }
    
        CCD0 ccd = new CCD0(tree.getLeafNodeCount(), false, CCD0ApproxMultiplier * tree.getLeafNodeCount());
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
        return "CCD0Approx";
    }

    @Override
    public String getDescription() {
        return "MAP (Approximate CCD0)";
    }

    @Override
    protected boolean sanityCheck(Tree mapTree, Tree firstTree, AbstractCCD ccd) {
        // nothing to do for CCD0 MAP tree
        // could analyse singly supported clades,
        // but that would mean we do post-analysis work of the users
        return true;
    }

}
