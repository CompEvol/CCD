package test.ccd.model;

import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import ccd.model.CCD0CP;
import ccd.model.CCD1CP;
import ccd.model.Clade;
import ccd.model.CladePartition;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * * Tree 1: ((A:0,B),(C,D))   — A is SA
 * * Tree 2: ((A,B:0),(C,D))   — B is SA
 * * Tree 3: (((A,B),C:0),D)   — C is SA
 * <p>
 * Expected probabilities under CCD1 (in eighty-firsts):
 * CP:  18, 18, 9
 * SJ:  27, 27, 27
 */

public class FourTaxaTest2 {
    private static final List<String> TAXA = Arrays.asList("A", "B", "C", "D");

    private static Tree parse(String newick) {
        return new TreeParser(TAXA, newick, 1, false);
    }

    private static List<Tree> sampleTreeList() {
        List<Tree> trees = new ArrayList<>();
        trees.add(parse("((A:0,B:1):1,(C:1,D:1):1):0;"));
        trees.add(parse("((A:1,B:0):1,(C:1,D:1):1):0;"));
        trees.add(parse("(((A:1,B:1):1,C:0):1,D:1):0;"));
        return trees;
    }

    @Test
    public void testProbabilities() {
        List<Tree> trees = sampleTreeList();
        CCD0CP ccd = new CCD0CP(trees, 0.0);

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
