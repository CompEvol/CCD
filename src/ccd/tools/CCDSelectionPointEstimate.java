package ccd.tools;

import beast.base.core.Description;
import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator;
import beastfx.app.treeannotator.TreeAnnotator.TreeSet;
import beastfx.app.treeannotator.services.TopologySettingService;
import ccd.model.AbstractCCD;
import ccd.model.CCD0;
import ccd.model.CCD1;
import ccd.model.CCD2;
import ccd.model.HeightSettingStrategy;

import java.io.IOException;
import java.io.PrintStream;

/**
 * This point estimate class computes a topology given a tree set by construct
 * a CCD0, a CCD1, and a CCD2, and then picking the best model based on the Akaike information criterion (AIC).
 */

@Description("TreeAnnotator plugin for setting the tree topology as CCD1 MAP tree")
public class CCDSelectionPointEstimate extends PointEstimate implements TopologySettingService {

    @Override
    public Tree setTopology(TreeSet treeSet, PrintStream progressStream, TreeAnnotator annotator)
            throws IOException {

        progressStream.println("CCDx MAP tree computation with x in {0, 1, 2} and AIC-based selection");
        progressStream.println("\nConstruct CCDs");
        progressStream.println("0              25             50             75            100");
        progressStream.println("|--------------|--------------|--------------|--------------|");


        treeSet.reset();
        Tree tree = treeSet.next();
        Tree firstTree = tree;
        int n = tree.getLeafNodeCount();
        CCD0 ccd0 = new CCD0(n, false);
        CCD1 ccd1 = new CCD1(n, false);
        CCD2 ccd2 = new CCD2(n, false);

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

            ccd0.addTree(tree);
            ccd1.addTree(tree);
            ccd2.addTree(tree);
            tree = treeSet.hasNext() ? treeSet.next() : null;
            i++;
        }
        progressStream.println();
        ccd0.setProgressStream(progressStream);
        ccd0.initialize();

        // since we constructed CCDs tree by tree, we need to set the TreeSet manually
        ccd0.setBaseTreeSet(treeSet);
        ccd1.setBaseTreeSet(treeSet);
        ccd2.setBaseTreeSet(treeSet);

        double aic0 = ccd0.getAICScore();
        double aic1 = ccd1.getAICScore();
        double aic2 = ccd2.getAICScore();

        progressStream.print("\nComputing AIC scores... ");
        progressStream.println("which are:");
        // progressStream.println("CCD0: " + aic0);
        // progressStream.println("CCD1: " + aic1);
        // progressStream.println("CCD2: " + aic2);
        progressStream.println("CCD0: " + aic0 + " (" + ccd0.getNumberOfClades() + " clades)");
        progressStream.println("CCD1: " + aic1 + " (" + ccd1.getNumberOfCladePartitions() + " clade partitions)");
        progressStream.println("CCD2: " + aic2 + " (" + ccd2.getNumberOfCladePartitions() + " clade partitions)");
        String winner = "";
        Tree mapTree = null;
        if ((aic0 < aic1) && (aic0 < aic2)) {
            progressStream.println("Hence, the CCD0 MAP tree is computed.\n");
            mapTree = ccd0.getMAPTree(HeightSettingStrategy.One);
        } else if ((aic1 <= aic0) && (aic1 < aic2)) {
            progressStream.println("Hence, the CCD1 MAP tree is computed.\n");
            mapTree = ccd1.getMAPTree(HeightSettingStrategy.One);
            sanityCheck(mapTree, firstTree, ccd1);
        } else {
            progressStream.println("Hence, the CCD2 MAP tree is computed.\n");
            mapTree = ccd2.getMAPTree(HeightSettingStrategy.One);
            sanityCheck(mapTree, firstTree, ccd2);
        }

        return mapTree;
    }

    @Override
    public String getServiceName() {
        return "CCDx-AIC";
    }

    @Override
    public String getDescription() {
        return "MAP (CCD, AIC selected)";
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

}
