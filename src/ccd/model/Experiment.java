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
 *
 *  Usage: java ... input_tree_file.java output_map_tree_file.java
 */
public class Experiment {
    static File TREE_FILE = new File("/Users/tobiaochsner/Documents/Thesis/Validation/data/mcmc_runs/yule-50_61.trees");
    static File OUTPUT_DIR = new File("/Users/tobiaochsner/Documents/Thesis/Test Hipstr");

    static int[] Ks = new int[]{0};
    static boolean RUN_CCD_1 = false;

    public static void main(String[] args) throws Exception {
        if (args.length == 2) {
            TREE_FILE = new File(args[0]);
            OUTPUT_DIR = new File(args[1]);
        }

        if (args.length == 3) {
            TREE_FILE = new File(args[0]);
            OUTPUT_DIR = new File(args[1]);
            Ks = new int[] {Integer.parseInt(args[2])};
        }

        runBenchmark();
    }

    private static void runBenchmark() throws IOException {
        TreeAnnotator.MemoryFriendlyTreeSet treeSet = new TreeAnnotator().new MemoryFriendlyTreeSet(TREE_FILE.getAbsolutePath(), 0);

        long timeStartNs;
        long timeEndNs;
        long durationS;
        Tree mapTree;

        // run CCD 1

        if (RUN_CCD_1) {
            timeStartNs = System.nanoTime();

            AbstractCCD ccd1 = new CCD1(treeSet, false);
            mapTree = ccd1.getMAPTree(HeightSettingStrategy.MeanOccurredHeights);

            timeEndNs = System.nanoTime();
            durationS = (int) Math.ceil((timeEndNs - timeStartNs) / 1e9);

            storeTrees(List.of(mapTree), new File(OUTPUT_DIR, TREE_FILE.getName() + ".map.ccd1." + durationS + "s.trees"));
            System.out.println("CCD1 took " + (timeEndNs - timeStartNs) / 1e9 + "s");

            System.out.println("-".repeat(20));
            System.out.println("CCD1 results:");
            evaluateTree(mapTree, ccd1);
        }

        // run CCD with different expansions

        for (int K : Ks) {
            timeStartNs = System.nanoTime();

            AbstractCCD ccd0 = new CCD0(treeSet, false, K);
            System.out.println(ccd0.getNumberOfCladePartitions());
            mapTree = ccd0.getMAPTree(HeightSettingStrategy.MeanOccurredHeights);

            timeEndNs = System.nanoTime();
            durationS = (int) Math.ceil((timeEndNs - timeStartNs) / 1e9);

            storeTrees(List.of(mapTree), new File(OUTPUT_DIR, TREE_FILE.getName() + ".map.ccd0-" + K + "." + durationS + "s.trees"));
            System.out.println("CCD0 with " + K + " expansion took " + (timeEndNs - timeStartNs) / 1e9 + "s");

            System.out.println("We have " + treeSet.totalTrees + " trees with " + ccd0.getNumberOfClades() + " clades");

            System.out.println("-".repeat(20));
            System.out.println("CCD0 results (expansion of " + K + "):");
            evaluateTree(mapTree, ccd0);
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

    private static void evaluateTree(Tree tree, AbstractCCD ccd) {
        List<Clade> clades = new ArrayList<>(ccd.getNumberOfLeaves() - 1);
        getClades(ccd, tree.getRoot(), clades);

        double sum = -clades.stream().mapToDouble(Clade::getCladeCredibility).map(Math::log).sum();

        int count99 = (int) clades.stream().filter(x -> (x.getCladeCredibility() >= 0.99)).count();
        int count95 = (int) clades.stream().filter(x -> (x.getCladeCredibility() >= 0.95)).count();
        int count50 = (int) clades.stream().filter(x -> (x.getCladeCredibility() >= 0.5)).count();

        double[] cladeCreds = clades.stream().mapToDouble(Clade::getCladeCredibility).sorted().toArray();
        Median medianCalculator = new Median();
        double median = medianCalculator.evaluate(cladeCreds);

        System.out.println("median = " + median);
        System.out.println("count99 = " + count99);
        System.out.println("count95 = " + count95);
        System.out.println("count50 = " + count50);
        System.out.println("sum = " + sum);
    }

    private static BitSet getClades(AbstractCCD ccd, Node vertex, List<Clade> clades) {
        BitSet bits = BitSet.newBitSet(ccd.getSizeOfLeavesArray());

        if (vertex.isLeaf()) {
            bits.set(vertex.getNr());
        } else {
            bits = getClades(ccd, vertex.getLeft(), clades);
            BitSet otherBits = getClades(ccd, vertex.getRight(), clades);
            bits.or(otherBits);

            clades.add(ccd.getClade(bits));
        }

        return bits;
    }

}
