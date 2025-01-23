package ccd.model;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator;
import org.apache.commons.math3.stat.descriptive.rank.Median;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

/**
 * Calculates MAP trees for different expansion strategies.
 * <p>
 * Usage: java ... input_tree_file.java output_map_tree_file.java
 */
public class Experiment2 {
    static File TREES_DIR = new File("/Users/tobiaochsner/Documents/Thesis/Test Hipstr/trees");
    static File OUTPUT_DIR = new File("/Users/tobiaochsner/Documents/Thesis/Test Hipstr/CCD0");

    public static void main(String[] args) throws Exception {
        runBenchmark();
    }

    private static void runBenchmark() throws IOException {
        for (File treeFile : TREES_DIR.listFiles()) {
            TreeAnnotator.MemoryFriendlyTreeSet treeSet = new TreeAnnotator().new MemoryFriendlyTreeSet(treeFile.getAbsolutePath(), 0);

            try {
                AbstractCCD ccd0 = new CCD0(treeSet, false, 0);
                System.out.println(ccd0.getNumberOfCladePartitions());
                Tree mapTree = ccd0.getMAPTree(HeightSettingStrategy.MeanOccurredHeights);

                storeTrees(List.of(mapTree), new File(OUTPUT_DIR, treeFile.getName()));
            } catch (Exception e) {}
        }
    }

    private static void storeTrees(List<Tree> trees, File treeFile) throws FileNotFoundException {
        try (PrintStream handle = new PrintStream(treeFile)) {
            trees.get(0).init(handle);

            for (int i = 0; i < trees.size(); i++) {
                Tree tree = trees.get(i);

                int id = (tree.getID() == null) ? i : Integer.parseInt(tree.getID().replace("_", ""));
                tree.log(id, handle);
                handle.println();
            }

            trees.get(trees.size() - 1).close(handle);
        }
    }

}
