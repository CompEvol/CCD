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

public class KregPointEstimateExperiment {

    public static String dataName;
    public static int subsampleSize;
    public static int startRep; // the starting index of reps, used in the case when nesi killed job half way
    public static int endRep; // the ending index of reps, used in the case when nesi killed job half way

    private static final double MU = KRegCCD.DEFAULT_MU;
    private static final double ALPHA = KRegCCD.DEFAULT_ALPHA;
    private static final int K = KRegCCD.DEFAULT_RESERVE_DEPTH;

    public static void main(String[] args) throws IOException {

        if (args.length > 1) { // if bash file pass arguments
            dataName = args[0];
            subsampleSize = Integer.parseInt(args[1]);
            startRep = Integer.parseInt(args[2]);
            endRep = Integer.parseInt(args[3]);
        } else { // default local test
            dataName = "Yule50";
            subsampleSize = 10;
            startRep = 1;
            endRep = 1;
        }

        // output file
        // String outputPathName = "/nesi/nobackup/uoa04397/sophie/smoothing/KReg/pointEst/output_" + dataName + "_sub" + subsampleSize + "_pointEst_" + startRep + ".csv";
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

            // String treesDir = "/nesi/nobackup/uoa04397/sophie/wcss_full/" + dataName + "/rep" + rep + "/run1/"
            //         + dataName.substring(0, 4).toLowerCase() + "-n" + dataName.substring(4, dataName.length()) + "-" + rep + ".trees";
            // String trueTreeDir = "/nesi/nobackup/uoa04397/sophie/wcss_full/" + dataName + "/rep" + rep + "/"
            //         + dataName.substring(0, 4).toLowerCase() + "-n" + dataName.substring(4, dataName.length()) + "-" + rep + "_true_phi.trees";

            String treesDir = "/Volumes/DYNABOOK/wcss_full/" + dataName + "/rep" + rep + "/run1/"
                    + dataName.substring(0, 4).toLowerCase() + "-n" + dataName.substring(4, dataName.length()) + "-" + rep + ".trees";
            String trueTreeDir = "/Volumes/DYNABOOK/wcss_full/" + dataName + "/rep" + rep + "/"
                    + dataName.substring(0, 4).toLowerCase() + "-n" + dataName.substring(4, dataName.length()) + "-" + rep + "_true_phi.trees";

            File treesFile = new File(treesDir);
            List<Tree> treesList = LoadOrStoreTrees.loadTrees(treesFile, 0, subsampleSize);
            TreeAnnotator.MemoryFriendlyTreeSet trueTreeSet = new TreeAnnotator().new MemoryFriendlyTreeSet(trueTreeDir, 0);
            trueTreeSet.reset();
            WrappedBeastTree trueTree = new WrappedBeastTree(trueTreeSet.next());

            // KRegCCD kreg = KRegCCD.withOptimisedMu(treesList);
            KRegCCD kreg = new KRegCCD(treesList, 0, MU, ALPHA, K);
            WrappedBeastTree mapTreeKreg = new WrappedBeastTree(kreg.getMAPTree());
            int rf = TreeDistances.robinsonsFouldDistance(trueTree, mapTreeKreg);

            sb = new StringBuilder();
            sb.append(rep).append(separator);
            sb.append(rf);
            writer.println(sb.toString());
            writer.flush();
        }
    }
}
