package test.ccd.model;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beastfx.app.inputeditor.BeautiPanelConfig;
import ccd.model.*;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * CCD0SJ — clade-frequency parameterisation (CCD0) on the SJ extended-clade
 * identity for sampled-ancestor trees.
 *
 * <p>Reuses the four-taxon worked example from {@code doc/sa-models.tex}:
 * three sampled trees on taxa {A,B,C,D}, each with count 1:
 * <pre>
 *   Tree 1: ((A,B:0),(C,D))   — B is SA
 *   Tree 2: ((A:0,B),(C,D))   — A is SA
 *   Tree 3: (((A,B),C:0),D)   — C is SA
 * </pre>
 *
 * <p>All three trees share the same root clade (taxa = ABCD, no SA at root),
 * so the placeholder root has a single (variant, empty) partition with
 * marginal 1. Under CCD0 normalisation, the root variant's partitions get
 * scored by the product of child clade-frequencies; with the JS-style expand
 * step ({@link CCD0SJ#expand()}) one new root-level partition
 * ({A,B}-noSA, {C,D}-noSA) is added, plus two new partitions on
 * {A,B,C}-SA on C pairing C with the SA-flagged {A,B} cherries.
 *
 * <p>Resulting tree probabilities:
 * <pre>
 *   P(Tree 1) = 2/9, P(Tree 2) = 2/9, P(Tree 3) = 1/9
 * </pre>
 * The remaining 4/9 mass is on three CCD0-expanded trees not in the sample.
 */
public class CCD0SJTest {

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

    private static List<Tree> unsampleTrees() {
        List<Tree> trees = new ArrayList<>();
        trees.add(parse("((A:1,B:1):1,(C:1,D:1):1):0;"));
        trees.add(parse("(((A:0,B:1):1,C:0):1,D:1):0;"));
        trees.add(parse("(((A:1,B:0):1,C:0):1,D:1):0;"));
        trees.add(parse("(((A:1,B:0):1,C:1):1,D:1):0;"));
        trees.add(parse("(((A:0,B:1):1,C:1):1,D:1):0;"));
        trees.add(parse("(((A:1,B:1):1,C:1):1,D:1):0;"));
        return trees;
    }

    @Test
    public void testCCD0SJProbabilities() {
        List<Tree> trees = sampleTrees();
        CCD0SJ ccd = new CCD0SJ(trees, 0.0);

        for (Clade clade : ccd.getClades()) {
            System.out.println(clade);
            for (CladePartition p : clade.getPartitions()) {
                System.out.println(p);
            }
        }

        double p1 = ccd.getProbabilityOfTree(trees.get(0));
        double p2 = ccd.getProbabilityOfTree(trees.get(1));
        double p3 = ccd.getProbabilityOfTree(trees.get(2));

        System.out.println("p1 = " + p1);

        assertEquals(2.0 / 9.0, p1, 1e-9, "P(Tree 1)");
        assertEquals(2.0 / 9.0, p2, 1e-9, "P(Tree 2)");
        assertEquals(1.0 / 9.0, p3, 1e-9, "P(Tree 3)");
    }

    @Test
    public void testCCD0CPProbabilities() {
        List<Tree> trees = sampleTrees();
        List<Tree> unsampledTrees = unsampleTrees();
        CCD0CP ccd = new CCD0CP(trees, 0.0);
        for (Clade clade : ccd.getClades()) {
            System.out.println(clade);
            for (CladePartition p : clade.getPartitions()) {
                System.out.println(p);
            }
        }
        double p1 = ccd.getProbabilityOfTree(trees.get(0));
        double p2 = ccd.getProbabilityOfTree(trees.get(1));
        double p3 = ccd.getProbabilityOfTree(trees.get(2));
        double p4 = ccd.getProbabilityOfTree(unsampledTrees.get(0));
        double p5 = ccd.getProbabilityOfTree(unsampledTrees.get(1));
        double p6 = ccd.getProbabilityOfTree(unsampledTrees.get(2));
        double p7 = ccd.getProbabilityOfTree(unsampledTrees.get(3));
        double p8 = ccd.getProbabilityOfTree(unsampledTrees.get(4));
        double p9 = ccd.getProbabilityOfTree(unsampledTrees.get(5));
        System.out.println("sampled tree1 = " + p1);
        System.out.println("sampled tree2 = " + p2);
        System.out.println("sampled tree3 = " + p3);
        System.out.println("unsampled 1 = " + p4);
        System.out.println("unsampled 2 = " + p5);
        System.out.println("unsampled 3 = " + p6);
        System.out.println("unsampled 4 = " + p7);
        System.out.println("unsampled 5 = " + p8);
        System.out.println("unsampled 6 = " + p9);
        double sum = p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9;
        System.out.println("all trees prob sum = " + sum);
    }

    @Test
    public void testRootCCPsSumToOne() {
        List<Tree> trees = sampleTrees();
        CCD0SJ ccd = new CCD0SJ(trees, 0.0);

        // Single root variant → placeholder root partition CCP = 1.0.
        assertEquals(1, ccd.getRootClade().getPartitions().size(),
                "placeholder root has one (variant, empty) partition (only one root variant)");
        assertEquals(1.0, ccd.getRootClade().getPartitions().get(0).getCCP(), 1e-9);

        // The variant's own partitions (the actual root-split distribution)
        // must sum to 1 under CCD0 normalisation. The variant is the non-empty
        // child of the placeholder root partition.
        CladePartition rootVariantPart = ccd.getRootClade().getPartitions().get(0);
        ccd.model.Clade variant = rootVariantPart.getChildClades()[0].size() == 0
                ? rootVariantPart.getChildClades()[1]
                : rootVariantPart.getChildClades()[0];
        double ccpSum = 0;
        for (CladePartition p : variant.getPartitions()) {
            ccpSum += p.getCCP();
        }
        assertEquals(1.0, ccpSum, 1e-9, "variant's partition CCPs sum to 1");
    }

    @Test
    public void testLogProbabilityMatches() {
        List<Tree> trees = sampleTrees();
        CCD0SJ ccd = new CCD0SJ(trees, 0.0);

        for (Tree t : trees) {
            double p = ccd.getProbabilityOfTree(t);
            double logp = ccd.getLogProbabilityOfTree(t);
            assertEquals(Math.log(p), logp, 1e-9, "log/prob consistency");
        }
    }

    /**
     * With no SA leaves, CCD0SJ should produce the same tree probabilities as
     * a standard CCD0 — the SA-flag bits are never set and the placeholder
     * root collapses to a single (variant, empty) partition with CCP 1.
     */
    @Test
    public void testNoSAReducesToCCD0() {
        List<Tree> trees = new ArrayList<>();
        trees.add(parse("((A:1,B:1):1,(C:1,D:1):1):0;"));
        trees.add(parse("((A:1,C:1):1,(B:1,D:1):1):0;"));
        trees.add(parse("((A:1,B:1):1,(C:1,D:1):1):0;"));

        CCD0 baseline = new CCD0(trees, 0.0);
        CCD0SJ sj = new CCD0SJ(trees, 0.0);

        for (Tree t : trees) {
            double pBaseline = baseline.getProbabilityOfTree(t);
            double pSj = sj.getProbabilityOfTree(t);
            assertEquals(pBaseline, pSj, 1e-9,
                    "CCD0SJ matches CCD0 when there are no SA leaves");
        }
    }

    @Test
    public void testSADetection() {
        List<Tree> trees = sampleTrees();
        CCD0SJ ccd = new CCD0SJ(trees, 0.0);

        boolean foundSAClade = false;
        for (ccd.model.Clade c : ccd.getClades()) {
            if (ccd.isSampledAncestor(c)) {
                foundSAClade = true;
                break;
            }
        }
        assertTrue(foundSAClade, "at least one clade should be flagged as carrying an SA");
    }

    /**
     * Canonical signature for a sampled tree: Newick over taxon names with
     * SA leaves (branch length 0) marked.
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
    public void testMAPTreePicksDominantVariant() {
        // Two trees with one root-split variant, one with the other; under
        // CCD0SJ scoring (subtree clade-credibility mass), the duplicated
        // variant has strictly greater MAP score and must win.
        List<Tree> trees = new ArrayList<>();
        trees.add(parse("((A:1,B:0):1,(C:1,D:1):1):0;"));
        trees.add(parse("((A:1,B:0):1,(C:1,D:1):1):0;"));
        trees.add(parse("(((A:1,B:1):1,C:0):1,D:2):0;"));

        CCD0SJ ccd = new CCD0SJ(trees, 0.0);
        Tree map = ccd.getMAPTree();

        assertNotNull(map);
        assertEquals(signature(trees.get(0).getRoot()), signature(map.getRoot()),
                "MAP tree should be the duplicated SA-on-B variant");
    }
}
