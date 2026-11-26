package ccd.model;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;

import org.junit.jupiter.api.Test;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;

/**
 * Trees with unnamed single-child nodes, such as SIEVE trees with a root above a trunk
 * (CompEvol/CCD#22), give the same CCD as the trees without them.
 */
public class SingleChildNodesTest {

    private static final List<String> TAXA = Arrays.asList("A", "B", "C", "D", "E");

    /** Binary trees without the root edge and ';'. The last one has a sampled ancestor (C:0). */
    private static final String[] BINARY = {
            "(((A:1,B:1):1,C:2):1,(D:1,E:1):2)",
            "(((A:1,C:1):1,B:2):1,(D:1,E:1):2)",
            "((A:2,(B:1,C:1):1):1,(D:1,E:1):2)",
            "(((A:1,B:1):1,C:2):1,(D:1,E:1):2)",
            "(((A:1,B:1):1,C:0):1,(D:1,E:1):1)",
    };

    private static Tree parse(String newick) {
        return new TreeParser(TAXA, newick, 1, false);
    }

    private static List<Tree> binaryTrees(boolean withSampledAncestor) {
        List<Tree> trees = new ArrayList<>();
        for (int i = 0; i < BINARY.length - (withSampledAncestor ? 0 : 1); i++) {
            trees.add(parse(BINARY[i] + ";"));
        }
        return trees;
    }

    /** The same trees, each under a root with a single child, as in SIEVE; one also has one inside. */
    private static List<Tree> singleChildTrees(boolean withSampledAncestor) {
        List<Tree> trees = new ArrayList<>();
        for (int i = 0; i < BINARY.length - (withSampledAncestor ? 0 : 1); i++) {
            String newick = BINARY[i];
            if (i == 1) {
                newick = newick.replace("(A:1,C:1):1", "((A:1,C:1):0.5):0.5");
            }
            trees.add(parse("(" + newick + ":0.5);"));
        }
        return trees;
    }

    private static void assertSameDistribution(String name, Function<List<Tree>, AbstractCCD> make,
                                               boolean withSampledAncestor) {
        List<Tree> binary = binaryTrees(withSampledAncestor);
        AbstractCCD expected = make.apply(binary);
        AbstractCCD actual = make.apply(singleChildTrees(withSampledAncestor));

        assertEquals(expected.getNumberOfClades(), actual.getNumberOfClades(), name + ": clades");
        for (Tree tree : binary) {
            assertEquals(expected.getProbabilityOfTree(tree), actual.getProbabilityOfTree(tree), 1e-12,
                    name + ": probability of " + tree.getRoot().toNewick());
        }
        // a tree with single-child nodes is scored as the tree without them
        for (Tree tree : singleChildTrees(withSampledAncestor)) {
            assertEquals(expected.getProbabilityOfTree(collapsed(tree)), actual.getProbabilityOfTree(tree), 1e-12,
                    name + ": probability of single-child tree");
        }
        assertEquals(clades(expected.getMAPTree()), clades(actual.getMAPTree()), name + ": MAP tree");
    }

    private static Tree collapsed(Tree tree) {
        return SingleChildNodes.collapse(tree, null);
    }

    /** The clades of a tree, as sets of taxon names. */
    private static Set<Set<String>> clades(Tree tree) {
        Set<Set<String>> clades = new HashSet<>();
        for (Node node : tree.getNodesAsArray()) {
            Set<String> clade = new HashSet<>();
            for (Node leaf : node.getAllLeafNodes()) {
                clade.add(leaf.getID());
            }
            if (node.isLeaf()) {
                clade.add(node.getID());
            }
            clades.add(clade);
        }
        return clades;
    }

    @Test
    void ccdsIgnoreUnnamedSingleChildNodes() {
        assertSameDistribution("CCD0", t -> new CCD0(t, 0), false);
        assertSameDistribution("CCD1", t -> new CCD1(t, 0), false);
        assertSameDistribution("CCD2", t -> new CCD2(t, 0), false);
    }

    @Test
    void sampledAncestorCCDsIgnoreUnnamedSingleChildNodes() {
        assertSameDistribution("CCD0EL", t -> new CCD0EL(t, 0), true);
        assertSameDistribution("CCD1EL", t -> new CCD1EL(t, 0), true);
        assertSameDistribution("CCD0EI", t -> new CCD0EI(t, 0), true);
        assertSameDistribution("CCD1EI", t -> new CCD1EI(t, 0), true);
    }

    @Test
    void collapsedCopyKeepsLeavesHeightsAndMetadata() {
        Tree tree = parse("((((A:1,B:1)[&x=1]:1,C:2):1,(D:1,E:1):2)[&posterior=0.5]:0.5);");
        Tree copy = collapsed(tree);
        assertEquals(tree.getLeafNodeCount(), copy.getLeafNodeCount());
        assertEquals(tree.getNodeCount() - 1, copy.getNodeCount());
        for (int i = 0; i < tree.getLeafNodeCount(); i++) {
            assertEquals(tree.getNode(i).getID(), copy.getNode(i).getID());
            assertEquals(tree.getNode(i).getHeight(), copy.getNode(i).getHeight());
        }
        assertEquals(3.0, copy.getRoot().getHeight(), 1e-12);
        assertEquals("posterior=0.5", copy.getRoot().metaDataString);
        assertEquals(copy.getNodeCount() - 1, copy.getRoot().getNr());
        for (int i = 0; i < copy.getNodeCount(); i++) {
            assertEquals(i, copy.getNode(i).getNr());
        }
    }

    @Test
    void namedSingleChildNodeIsRefused() {
        // a sampled ancestor written as a labelled internal node
        Tree tree = parse("(((A:1,B:1):1)F:1,(C:1,(D:0.5,E:0.5):0.5):1);");
        IllegalArgumentException e = assertThrows(IllegalArgumentException.class, () -> new CCD0(List.of(tree), 0));
        assertTrue(e.getMessage().contains("'F'"), e.getMessage());
        assertTrue(e.getMessage().contains("sampled ancestor"), e.getMessage());
    }

    @Test
    void collapsingCanBeTurnedOff() {
        CCD1 ccd = new CCD1(5, false);
        ccd.setCollapseSingleChildNodes(false);
        Tree tree = singleChildTrees(false).get(0);
        IllegalArgumentException e = assertThrows(IllegalArgumentException.class, () -> ccd.addTree(tree));
        assertTrue(e.getMessage().contains("ccd.collapseSingleChildNodes"), e.getMessage());

        ccd.addTree(binaryTrees(false).get(0)); // binary trees are still fine
        assertEquals(1, ccd.getNumberOfBaseTrees());
    }

    @Test
    void removeTreeTakesOutTheCollapsedCopy() {
        List<Tree> trees = singleChildTrees(false);
        CCD1 ccd = new CCD1(5, true);
        for (Tree tree : trees) {
            ccd.addTree(tree);
        }
        ccd.removeTree(trees.get(1), true);
        assertEquals(trees.size() - 1, ccd.getNumberOfBaseTrees());
        assertEquals(trees.size() - 1, ccd.getBaseTrees().size());
        assertEquals(0, ccd.getProbabilityOfTree(trees.get(1)), 1e-12, "the only tree with clade {A,C} is gone");
    }

    @Test
    void cladeToNodeMapUsesTheGivenTreesNodes() {
        List<Tree> trees = singleChildTrees(false);
        CCD1 ccd = new CCD1(trees, 0);
        Tree tree = trees.get(0);
        Map<Clade, Node> map = ccd.getCladeToNodeMap(tree);
        assertEquals(tree.getNodeCount() - 1, map.size(), "every node but the single-child root");
        Set<Node> nodes = new HashSet<>(Arrays.asList(tree.getNodesAsArray()));
        for (Node node : map.values()) {
            assertTrue(nodes.contains(node));
        }
        assertSame(tree.getRoot().getChild(0), map.get(ccd.getRootClade()));
    }
}
