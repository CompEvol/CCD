package ccd.model;

import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import org.junit.jupiter.api.Test;

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
public class CCD1EITest {

    private static final List<String> TAXA = Arrays.asList("A", "B", "C", "D");

    private static Tree parse(String newick) {
        return new TreeParser(TAXA, newick, 1, false);
    }

    private static List<Tree> sampleTrees() {
        List<Tree> trees = new ArrayList<>();
        trees.add(parse("((A:1,B:0):1,(C:1,D:1):1):0;"));
        trees.add(parse("((A:0,B:1):1,(C:1,D:1):1):0;"));
        trees.add(parse("(((A:1,B:1):1,C:0):1,D:2):0;"));
        return trees;
    }

    private static double rounded(double x) {
        return Math.round(x * 81 * 1e6) / 1e6 / 81.0;
    }

    @Test
    public void testCCD1EIProbabilities() {
        List<Tree> trees = sampleTrees();
        CCD1EI ei = new CCD1EI(trees, 0.0);

        for (Clade clade : ei.getClades()) {
            System.out.println(clade);
            for (CladePartition p : clade.getPartitions()) {
                System.out.println(p);
            }
        }

        double p1 = ei.getProbabilityOfTree(trees.get(0));
        double p2 = ei.getProbabilityOfTree(trees.get(1));
        double p3 = ei.getProbabilityOfTree(trees.get(2));

        assertEquals(27.0 / 81.0, rounded(p1), 1e-9, "EI P(Tree 1)");
        assertEquals(27.0 / 81.0, rounded(p2), 1e-9, "EI P(Tree 2)");
        assertEquals(27.0 / 81.0, rounded(p3), 1e-9, "EI P(Tree 3)");

        assertEquals(1.0, p1 + p2 + p3, 1e-9, "EI total mass on sampled trees");
    }

    @Test
    public void testCCD1ELProbabilities() {
        List<Tree> trees = sampleTrees();
        CCD1EL el = new CCD1EL(trees, 0.0);

        double p1 = el.getProbabilityOfTree(trees.get(0));
        double p2 = el.getProbabilityOfTree(trees.get(1));
        double p3 = el.getProbabilityOfTree(trees.get(2));

        assertEquals(18.0 / 81.0, rounded(p1), 1e-9, "EL P(Tree 1)");
        assertEquals(18.0 / 81.0, rounded(p2), 1e-9, "EL P(Tree 2)");
        assertEquals(9.0 / 81.0, rounded(p3), 1e-9, "EL P(Tree 3)");
    }

    @Test
    public void testEILogProbabilityMatches() {
        List<Tree> trees = sampleTrees();
        CCD1EI ei = new CCD1EI(trees, 0.0);

        for (Tree t : trees) {
            double p = ei.getProbabilityOfTree(t);
            double logp = ei.getLogProbabilityOfTree(t);
            assertEquals(Math.log(p), logp, 1e-9, "log/prob consistency");
        }
    }

    /**
     * Canonical signature for a sampled tree: Newick over taxon names with
     * SA leaves (branch length 0) marked. Used to bucket sampled trees by
     * topology+SA-pattern for the empirical-distribution check.
     */
    private static String signature(beast.base.evolution.tree.Node node) {
        if (node.isLeaf()) {
            String name = node.getID();
            boolean sa = node.getLength() == 0.0;
            return name + (sa ? ":SA" : "");
        }
        java.util.List<String> kids = new java.util.ArrayList<>();
        for (beast.base.evolution.tree.Node c : node.getChildren()) {
            kids.add(signature(c));
        }
        java.util.Collections.sort(kids);
        return "(" + kids.get(0) + "," + kids.get(1) + ")";
    }

    @Test
    public void testEISamplingMatchesAnalytical() {
        List<Tree> trees = sampleTrees();
        CCD1EI ei = new CCD1EI(trees, 0.0);

        String sig1 = signature(trees.get(0).getRoot());
        String sig2 = signature(trees.get(1).getRoot());
        String sig3 = signature(trees.get(2).getRoot());

        java.util.Map<String, Integer> counts = new java.util.HashMap<>();
        int N = 30000;
        for (int i = 0; i < N; i++) {
            Tree sampled = ei.sampleTree();
            String sig = signature(sampled.getRoot());
            counts.merge(sig, 1, Integer::sum);
        }

        int c1 = counts.getOrDefault(sig1, 0);
        int c2 = counts.getOrDefault(sig2, 0);
        int c3 = counts.getOrDefault(sig3, 0);

        // Sampled-only support: nothing outside the three known trees.
        assertEquals(N, c1 + c2 + c3, "all samples land on one of the three sampled trees");

        // Each tree gets ~1/3. With N=30000, std dev ~sqrt(N*p*(1-p))=~86.
        // Tolerance 5σ ~ 430 → use 500 absolute, plenty of headroom.
        double p1 = c1 / (double) N;
        double p2 = c2 / (double) N;
        double p3 = c3 / (double) N;
        assertEquals(1.0 / 3.0, p1, 0.02, "empirical P(Tree 1)");
        assertEquals(1.0 / 3.0, p2, 0.02, "empirical P(Tree 2)");
        assertEquals(1.0 / 3.0, p3, 0.02, "empirical P(Tree 3)");
    }

    @Test
    public void testMAPTreePicksDominantVariant() {
        // Tree 1 appears twice → marginal 2/4; Tree 2 and Tree 3 once each → 1/4.
        // With distinct counts there's no root-level tie and MAP must return
        // the Tree-1 variant (including its SA pattern: B is the SA).
        List<Tree> trees = new ArrayList<>();
        trees.add(parse("((A:1,B:0):1,(C:1,D:1):1):0;"));
        trees.add(parse("((A:1,B:0):1,(C:1,D:1):1):0;"));
        trees.add(parse("((A:0,B:1):1,(C:1,D:1):1):0;"));
        trees.add(parse("(((A:1,B:1):1,C:0):1,D:2):0;"));

        CCD1EI ei = new CCD1EI(trees, 0.0);
        Tree map = ei.getMAPTree();

        assertEquals(signature(trees.get(0).getRoot()), signature(map.getRoot()),
                "MAP tree should be the duplicated Tree-1 variant");
    }
}
