package ccd.experiments;

import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator;
import ccd.algorithms.LoadOrStoreTrees;
import ccd.algorithms.TreeDistances;
import ccd.model.KRegCCD;
import ccd.model.WrappedBeastTree;


import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;

public class KRegPointEstimateTime {

    public static String dataName;
    public static int subsampleSize;
    public static int startRep;
    public static int endRep;

    private static final double MU = KRegCCD.DEFAULT_MU;
    private static final double ALPHA = KRegCCD.DEFAULT_ALPHA;
    private static final int K = KRegCCD.DEFAULT_RESERVE_DEPTH;

    public static void main(String[] args) throws IOException {

        if (args.length > 1) {
            dataName = args[0];
            subsampleSize = Integer.parseInt(args[1]);
            startRep = Integer.parseInt(args[2]);
            endRep = Integer.parseInt(args[3]);
        } else {
            dataName = "Coal320";
            subsampleSize = 1000;
            startRep = 37;
            endRep = 37;
        }

        long experimentStart = System.nanoTime();

        // output file
        String outputPathName = "/Users/sophieyang/Desktop/local_test/output_" + dataName + "_sub" + subsampleSize + "_pointEst_" + startRep + ".csv";
        File outputFile = new File(outputPathName);
        FileWriter fileWriter = new FileWriter(outputFile);
        PrintWriter writer = new PrintWriter(fileWriter);
        // header
        String separator = ",";
        StringBuilder sb = new StringBuilder();
        sb.append("rep").append(separator);
        sb.append("RF_kreg");
        writer.println(sb.toString());
        writer.flush();

        for (int rep = startRep; rep <= endRep; rep++) {
            long repStart = System.nanoTime();
            System.out.println();
            System.out.println("==================================================");
            System.out.println("Starting replicate " + rep);
            System.out.println("==================================================");

            String treesDir = "/Users/sophieyang/Desktop/wcss_full/" + dataName + "/rep" + rep + "/run1/"
                    + dataName.substring(0, 4).toLowerCase() + "-n" + dataName.substring(4, dataName.length()) + "-" + rep + ".trees";
            String trueTreeDir = "/Users/sophieyang/Desktop/wcss_full/" + dataName + "/rep" + rep + "/"
                    + dataName.substring(0, 4).toLowerCase() + "-n" + dataName.substring(4, dataName.length()) + "-" + rep + "_true_phi.trees";
            // ----------------------------
            // Load trees
            // ----------------------------
            long t = System.nanoTime();
            File treesFile = new File(treesDir);
            List<Tree> treesList = LoadOrStoreTrees.loadTrees(treesFile, 0, subsampleSize);
            TreeAnnotator.MemoryFriendlyTreeSet trueTreeSet = new TreeAnnotator().new MemoryFriendlyTreeSet(trueTreeDir, 0);
            trueTreeSet.reset();
            WrappedBeastTree trueTree = new WrappedBeastTree(trueTreeSet.next());
            System.out.printf("Loaded trees in %.2f s%n", (System.nanoTime() - t) / 1e9);
            // ----------------------------
            // Build KReg
            // ----------------------------
            t = System.nanoTime();
            System.out.println("Building KReg...");
            // KRegCCD kreg = KRegCCD.withOptimisedMu(treesList);
            KRegCCD kreg = new KRegCCD(treesList, 0, MU, ALPHA, K);
            System.out.printf("Built KReg in %.2f s%n", (System.nanoTime() - t) / 1e9);
            // ----------------------------
            // MAP tree
            // ----------------------------
            t = System.nanoTime();
            WrappedBeastTree mapTreeKreg = new WrappedBeastTree(kreg.getMAPTree());
            System.out.printf("Computed MAP tree in %.2f s%n", (System.nanoTime() - t) / 1e9);
            // ----------------------------
            // RF distance
            // ----------------------------
            t = System.nanoTime();
            int rf = TreeDistances.robinsonsFouldDistance(trueTree, mapTreeKreg);
            System.out.printf("Computed RF distance in %.2f s%n", (System.nanoTime() - t) / 1e9);
            // ----------------------------
            // Write output
            // ----------------------------
            t = System.nanoTime();
            sb = new StringBuilder();
            sb.append(rep).append(separator);
            sb.append(rf);
            writer.println(sb.toString());
            writer.flush();
            System.out.printf("Wrote output in %.2f s%n", (System.nanoTime() - t) / 1e9);
            // ----------------------------
            // Rep summary
            // ----------------------------
            double repTime = (System.nanoTime() - repStart) / 1e9;
            System.out.println("------------------------------------------");
            System.out.printf("Replicate %d completed in %.2f seconds (%.2f minutes)%n", rep, repTime, repTime / 60.0);
            System.out.println("------------------------------------------");
        }
        writer.close();
        double totalTime = (System.nanoTime() - experimentStart) / 1e9;
        System.out.println();
        System.out.println("==================================================");
        System.out.printf("Entire experiment completed in %.2f seconds (%.2f minutes, %.2f hours)%n",
                totalTime,
                totalTime / 60.0,
                totalTime / 3600.0);
        System.out.println("==================================================");
    }
}

