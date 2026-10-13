package ccd.model;

import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class FourTaxaRootSATest2 {
    private static final List<String> TAXA = Arrays.asList("A", "B", "C", "D");

    private static Tree parse(String newick) {
        return new TreeParser(TAXA, newick, 1, false);
    }

    private static List<Tree> sampleTreeList1() {
        List<Tree> trees = new ArrayList<>();
        trees.add(parse("(((A:1,B:0):1,C:1):1,D:0):0"));
        trees.add(parse("(((A:0,B:1):1,C:1):1,D:1):0"));
        trees.add(parse("((A:1,B:1):1,(C:1,D:0):1):0"));
        return trees;
    }

    private static List<Tree> sampleTreeList() {
        List<Tree> trees = new ArrayList<>();
        trees.add(parse("((A:0,B:1):1,(C:1,D:1):1):0;"));
        trees.add(parse("((A:1,B:0):1,(C:1,D:1):1):0;"));
        trees.add(parse("(((A:1,B:0):1,C:1):1,D:0):0;"));
        return trees;
    }

    @Test
    public void testProbabilities() {
        List<Tree> trees = sampleTreeList();
        // CCD1EL ccd = new CCD1EL(trees, 0.0);
        CCD0EL ccd = new CCD0EL(trees, 0.0);
        // CCD1EI ccd = new CCD1EI(trees, 0.0);
        // CCD0EI ccd = new CCD0EI(trees, 0.0);

        System.out.println("num of partitions = " + ccd.getNumberOfCladePartitions());

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
