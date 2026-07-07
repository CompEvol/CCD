package test.ccd.model;

import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import ccd.model.*;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * * Tree 1: ((A,B),C)   — no SA
 * * Tree 2: ((A,B),C:0)   — C is SA
 * * Tree 3: ((A,B:0),C:0)   — B,C are SA
 * <p>
 * Expected probabilities under CCD1:
 * CP:  2/9, 4/9, 2/9
 * SJ:  1/3, 1/3, 1/3
 */

public class ThreeTaxaRootSATest {

    private static final List<String> TAXA = Arrays.asList("A", "B", "C");

    private static Tree parse(String newick) {
        return new TreeParser(TAXA, newick, 1, false);
    }

    private static List<Tree> sampleTreeList() {
        List<Tree> trees = new ArrayList<>();
        trees.add(parse("((A:1,B:1):1,C:1):0;"));
        trees.add(parse("((A:1,B:1):1,C:0):0;"));
        trees.add(parse("((A:1,B:0):1,C:0):0;"));
        return trees;
    }

    @Test
    public void testProbabilities() {
        List<Tree> trees = sampleTreeList();
        CCD0CP ccd = new CCD0CP(trees, 0.0);
        // CCD0SJ ccd = new CCD0SJ(trees, 0.0);

        for (Clade clade : ccd.getClades()) {
            System.out.println(clade);
            for (CladePartition p : clade.getPartitions()) {
                System.out.println(p);
            }
        }

        System.out.println("p(tree1) = " + ccd.getProbabilityOfTree(trees.get(0)));
        System.out.println("p(tree2) = " + ccd.getProbabilityOfTree(trees.get(1)));
        System.out.println("p(tree3) = " + ccd.getProbabilityOfTree(trees.get(2)));
    }
}
