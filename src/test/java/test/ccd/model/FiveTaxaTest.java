package test.ccd.model;

import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import ccd.model.CCD0CP;
import ccd.model.Clade;
import ccd.model.CladePartition;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class FiveTaxaTest {
    private static final List<String> TAXA = Arrays.asList("A", "B", "C", "D", "E");

    private static Tree parse(String newick) {
        return new TreeParser(TAXA, newick, 1, false);
    }

    private static List<Tree> sampleTreeList() {
        List<Tree> trees = new ArrayList<>();
        trees.add(parse("((A:0,B:1):1,(C:1,(D:1,E:1):1):1):0;"));
        trees.add(parse("(((A:1,B:1):1,C:0):1,(D:1,E:1):1):0;"));
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
    }
}
