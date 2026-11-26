package ccd.model;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.junit.jupiter.api.Test;

import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.util.Randomizer;
import ccd.algorithms.credibleSets.ProbabilityBasedCredibleSetComputer;

/**
 * CCD2.getLogProbabilityOfTree is the log of getProbabilityOfTree (#24). It used to be
 * inherited from AbstractCCD, which walks plain clades and gave -Infinity for every tree.
 */
public class CCD2LogProbabilityTest {

    private static final List<String> TAXA = Arrays.asList("A", "B", "C", "D", "E");

    private static Tree parse(String newick) {
        return new TreeParser(TAXA, newick, 1, false);
    }

    private static List<Tree> trees() {
        List<Tree> trees = new ArrayList<>();
        for (String newick : new String[]{
                "(((A:1,B:1):1,C:2):1,(D:1,E:1):2);",
                "(((A:1,B:1):1,C:2):1,(D:1,E:1):2);",
                "(((A:1,B:1):1,C:2):1,(D:1,E:1):2);",
                "(((A:1,C:1):1,B:2):1,(D:1,E:1):2);",
                "((A:2,(B:1,C:1):1):1,(D:1,E:1):2);",
                "(((A:1,B:1):1,(C:1,D:1):1):1,E:3);",
                "((((A:1,B:1):1,C:2):1,D:3):1,E:4);"}) {
            trees.add(parse(newick));
        }
        return trees;
    }

    @Test
    void logProbabilityIsLogOfProbability() {
        List<Tree> trees = trees();
        CCD2 ccd = new CCD2(trees, 0);
        for (Tree tree : trees) {
            double p = ccd.getProbabilityOfTree(tree);
            assertTrue(p > 0, "a base tree has positive probability");
            assertEquals(Math.log(p), ccd.getLogProbabilityOfTree(tree), 1e-12, tree.getRoot().toNewick());
        }
        Tree map = ccd.getMAPTree();
        assertEquals(Math.log(ccd.getProbabilityOfTree(map)), ccd.getLogProbabilityOfTree(map), 1e-12);
    }

    @Test
    void treeOutsideTheSupportHasLogProbabilityMinusInfinity() {
        CCD2 ccd = new CCD2(trees(), 0);
        Tree unseen = parse("(((A:1,E:1):1,C:2):1,(D:1,B:1):2);");
        assertEquals(0, ccd.getProbabilityOfTree(unseen));
        assertEquals(Double.NEGATIVE_INFINITY, ccd.getLogProbabilityOfTree(unseen));
        // every clade and split is in CCD1, but {B,C}|{A} under {A,B,C} only next to {D,E}, not next to {D}
        Tree newPairing = parse("((((B:1,C:1):1,A:2):1,D:3):1,E:4);");
        assertTrue(new CCD1(trees(), 0).getProbabilityOfTree(newPairing) > 0);
        assertEquals(0, ccd.getProbabilityOfTree(newPairing));
        assertEquals(Double.NEGATIVE_INFINITY, ccd.getLogProbabilityOfTree(newPairing));
    }

    @Test
    void credibleLevelOfTheMapTreeIsFound() {
        Randomizer.setSeed(24);
        CCD2 ccd = new CCD2(trees(), 0);
        ProbabilityBasedCredibleSetComputer cred = new ProbabilityBasedCredibleSetComputer(ccd, 1000);
        double level = cred.getCredibleLevel(ccd.getMAPTree());
        assertTrue(level > 0 && level <= 1, "credible level of the MAP tree was " + level);
    }
}
