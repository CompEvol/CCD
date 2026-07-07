package ccd.tools;

import java.io.IOException;
import java.io.PrintStream;

import beast.base.core.Citation;
import beast.base.core.Description;
import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator;
import beastfx.app.treeannotator.TreeAnnotator.TreeSet;
import beastfx.app.treeannotator.services.TopologySettingService;
import ccd.model.AbstractCCD;
import ccd.model.CCD0;
import ccd.model.HeightSettingStrategy;


@Citation(value = "Baele et al. (2024). bioRxiv.\n" +
        "HIPSTR: highest independent posterior subtree reconstruction in TreeAnnotator X.",
        DOI = "https://doi.org/10.1101/2024.12.08.627395")
@Description("TreeAnnotator plugin for setting the tree topology as HIPSTR tree = CCD0 MAP tree without expand ")
public class HIPSTR extends PointEstimate implements TopologySettingService {
	
    @Override
    public Tree setTopology(TreeSet treeSet, PrintStream progressStream, TreeAnnotator annotator)
            throws IOException {

        progressStream.println("HIPSTR tree computation");
        progressStream.println("0              25             50             75            100");
        progressStream.println("|--------------|--------------|--------------|--------------|");

        treeSet.reset();
        Tree tree = treeSet.next();
        Tree firstTree = tree;
    
        CCD0 ccd = new CCD0(tree.getLeafNodeCount(), false, 0);
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
        return "HIPSTR";
    }

    @Override
    public String getDescription() {
        return "HIPSTR (CCD0 without expansion)";
    }

    @Override
    protected boolean sanityCheck(Tree mapTree, Tree firstTree, AbstractCCD ccd) {
        // nothing to do for HIPSTR tree
        // could analyse singly supported clades,
        // but that would mean we do post-analysis work of the users
        return true;
    }

}
