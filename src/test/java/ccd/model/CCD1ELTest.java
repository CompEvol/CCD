package ccd.model;

import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beastfx.app.treeannotator.TreeAnnotator;
import ccd.algorithms.TreeDistances;
import org.junit.jupiter.api.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Reproduces the four-taxon worked example from doc/sa-models.tex.
 * <p>
 * Three sampled trees on taxa {A,B,C,D}, each with count 1:
 * Tree 1: ((A,B:0),(C,D))   — B is SA
 * Tree 2: ((A:0,B),(C,D))   — A is SA
 * Tree 3: (((A,B),C:0),D)   — C is SA
 * <p>
 * Expected probabilities (in eighty-firsts):
 * EL:  18, 18, 9
 * EI:  27, 27, 27
 */

public class CCD1ELTest {

    private static final List<String> TAXA = Arrays.asList("A", "B", "C", "D");

    private static Tree parse(String newick) {
        return new TreeParser(TAXA, newick, 1, false);
    }

    private static List<Tree> sampleTreeList() {
        List<Tree> trees = new ArrayList<>();
        trees.add(parse("((A:1,B:0):1,(C:1,D:1):1):0;"));
        trees.add(parse("((A:0,B:1):1,(C:1,D:1):1):0;"));
        trees.add(parse("(((A:1,B:1):1,C:0):1,D:2):0;"));
        return trees;
    }

    @Test
    public void graphSizeTest() {
        List<Tree> trees = sampleTreeList();
        CCD1EL ccd1el = new CCD1EL(trees, 0.0);
        CCD1EI ccd1ei = new CCD1EI(trees, 0.0);
        CCD0EL ccd0el = new CCD0EL(trees, 0.0);
        CCD0EI ccd0ei = new CCD0EI(trees, 0.0);
        System.out.println("ccd1 EL # clades = " + ccd1el.getNumberOfClades());
        System.out.println("ccd1 EL # splits = " + ccd1el.getNumberOfCladePartitions());
        System.out.println("ccd1 EI # clades = " + ccd1ei.getNumberOfClades());
        for (Clade clade : ccd1ei.getClades()) {
            System.out.println(clade);
        }
        System.out.println("ccd1 EI # splits = " + ccd1ei.getNumberOfCladePartitions());
        System.out.println("ccd0 EL # clades = " + ccd0el.getNumberOfClades());
        System.out.println("ccd0 EL # splits = " + ccd0el.getNumberOfCladePartitions());
        System.out.println("ccd0 EI # clades = " + ccd0ei.getNumberOfClades());
        System.out.println("ccd0 EI # splits = " + ccd0ei.getNumberOfCladePartitions());
    }

    @Test
    public void distanceTest() {
        List<Tree> trees = sampleTreeList();
        WrappedBeastTreeWithSampledAncestor tree1 = new WrappedBeastTreeWithSampledAncestor(trees.get(0));
        WrappedBeastTreeWithSampledAncestor tree2 = new WrappedBeastTreeWithSampledAncestor(trees.get(1));
        System.out.println("RF distance (SA WrappedBeastTree) = " + TreeDistances.robinsonsFouldDistance(tree1, tree2));
        WrappedBeastTree tree1NonSA = new WrappedBeastTree(trees.get(0));
        WrappedBeastTree tree2NonSA = new WrappedBeastTree(trees.get(1));
        System.out.println("RF distance (non-SA WrappedBeastTree) = " + TreeDistances.robinsonsFouldDistance(tree1NonSA, tree2NonSA));
    }

    private static TreeAnnotator.MemoryFriendlyTreeSet sampleTreeSet() throws IOException {
        String dataPath = "/Users/zyan598/Documents/GitHub/CCD_sampled_ancestor/example_trees/4taxa_alexei.trees";
        TreeAnnotator.MemoryFriendlyTreeSet treeSet = new TreeAnnotator().new MemoryFriendlyTreeSet(dataPath, 0);
        return treeSet;
    }

    private static double rounded(double x) {
        return Math.round(x * 81 * 1e6) / 1e6 / 81.0;
    }

    @Test
    public void treeSetTest() throws IOException {
        org.junit.jupiter.api.Assumptions.assumeTrue(
                new java.io.File("/Users/zyan598/Documents/GitHub/CCD_sampled_ancestor/example_trees/4taxa_alexei.trees").exists(),
                "skipped: SA example tree file not present on this machine");
        TreeAnnotator.MemoryFriendlyTreeSet treeSet = sampleTreeSet();
        // CCD1 ccd = new CCD1EL(treeSet, false);
        // WrappedBeastTree mapTreeCCD1 = new WrappedBeastTree(ccd.getMAPTree());
        // System.out.println("CCD MAP tree log prob = " + ccd.getMaxLogTreeProbability());
    }

    @Test
    public void testCCD1ELProbabilities() {
        List<Tree> trees = sampleTreeList();
        CCD1EL el = new CCD1EL(trees, 0.0);

        for (Clade clade : el.getClades()) {
            System.out.println(clade);
        }

        double p1 = el.getProbabilityOfTree(trees.get(0));
        double p2 = el.getProbabilityOfTree(trees.get(1));
        double p3 = el.getProbabilityOfTree(trees.get(2));

        assertEquals(18.0 / 81.0, rounded(p1), 1e-9, "EL P(Tree 1)");
        assertEquals(18.0 / 81.0, rounded(p2), 1e-9, "EL P(Tree 2)");
        assertEquals(9.0 / 81.0, rounded(p3), 1e-9, "EL P(Tree 3)");

    }
}
