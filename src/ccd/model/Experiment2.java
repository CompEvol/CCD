package ccd.model;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.parser.NexusParser;
import beastfx.app.treeannotator.TreeAnnotator;
import ccd.algorithms.TreeDistances;
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
        for (File treeFile : OUTPUT_DIR.listFiles()) {
            if (!treeFile.getName().endsWith("trees")) continue;

            File file1 = new File(new File("/Users/tobiaochsner/Documents/Thesis/Test Hipstr/CCD0"), treeFile.getName());
            File file2 = new File(new File("/Users/tobiaochsner/Documents/Thesis/Test Hipstr/Hipstr"), treeFile.getName());

            if (!file1.exists() || !file2.exists()) continue;

            Tree tree1 = loadTrees(file1).get(0);
            Tree tree2 = loadTrees(file2).get(0);

            System.out.println(treeFile.getName() + ": " + TreeDistances.robinsonsFouldDistance(
                    new WrappedBeastTree(tree1),
                    new WrappedBeastTree(tree2)
            ));
        }

//        runBenchmark();
    }

    private static void runBenchmark() throws IOException {
        for (File treeFile : TREES_DIR.listFiles()) {
            if (!treeFile.getName().endsWith("trees")) continue;

            try {
                List<Tree> trees = loadTrees(treeFile);
                int numLeaves = trees.get(0).getLeafNodeCount();

                AbstractCCD ccd0 = new CCD0(numLeaves, 0);

                for (Tree tree : trees) {
                    ccd0.addTree(tree);
                }
                ccd0.initialize();

                System.out.println(ccd0.getNumberOfCladePartitions());
                Tree mapTree = ccd0.getMAPTree(HeightSettingStrategy.MeanOccurredHeights);

                storeTrees(List.of(mapTree), new File(OUTPUT_DIR, treeFile.getName()));
            } catch (Exception e) {
                System.out.println(e);
            }
        }
    }

    public static List<Tree> loadTrees(File treeFile) throws IOException {
        NexusParser parser = new NexusParser();
        parser.parseFile(treeFile);
        return parser.trees;
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
