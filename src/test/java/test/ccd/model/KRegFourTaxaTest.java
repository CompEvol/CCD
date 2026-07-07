package test.ccd.model;

import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import ccd.model.Clade;
import ccd.model.CladePartition;
import ccd.model.KRegCCD;
import ccd.model.RegCCD;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class KRegFourTaxaTest {

    private static final List<String> TAXA = Arrays.asList("A", "B", "C", "D");

    private static Tree parse(String newick) {
        return new TreeParser(TAXA, newick, 1, false);
    }

    private static List<Tree> sampleTreeList() {
        List<Tree> trees = new ArrayList<>();
        trees.add(parse("(((A:1,B:1):1,C:1):1,D:1):0;"));
        trees.add(parse("(((D:1,C:1):1,B:1):1,A:1):0;"));
        return trees;
    }

    private static List<Tree> expandedTrees() {
        List<Tree> trees = new ArrayList<>();
        trees.add(parse("((A:1,B:1):1,(C:1,D:1):1):0;"));
        return trees;
    }

    private static List<Tree> kregTrees() {
        List<Tree> trees = new ArrayList<>();
        trees.add(parse("(((A:1,B:1):1,D:1):1,C:1):0;"));
        trees.add(parse("(((C:1,D:1):1,A:1):1,B:1):0;"));
        trees.add(parse("(((A:1,C:1):1,B:1):1,D:1):0;"));
        trees.add(parse("(((B:1,C:1):1,A:1):1,D:1):0;"));
        trees.add(parse("(((B:1,C:1):1,D:1):1,A:1):0;"));
        trees.add(parse("(((B:1,D:1):1,C:1):1,A:1):0;"));
        trees.add(parse("(((A:1,D:1):1,B:1):1,C:1):0;"));
        trees.add(parse("(((B:1,D:1):1,A:1):1,C:1):0;"));
        trees.add(parse("(((A:1,C:1):1,D:1):1,B:1):0;"));
        trees.add(parse("(((A:1,D:1):1,C:1):1,B:1):0;"));
        trees.add(parse("((A:1,C:1):1,(B:1,D:1):1):0;"));
        trees.add(parse("((A:1,D:1):1,(B:1,C:1):1):0;"));
        return trees;
    }

    @Test
    public void testProbabilities() {
        List<Tree> trees = sampleTreeList();
        List<Tree> expandedTrees = expandedTrees();
        List<Tree> kregTrees = kregTrees();
        KRegCCD ccd = KRegCCD.withOptimisedParameters(trees);
        // RegCCD ccd = new RegCCD(trees,0);

        for (Clade clade : ccd.getClades()) {
            System.out.println(clade);
            for (CladePartition p : clade.getPartitions()) {
                System.out.println(p);
            }
        }

        System.out.println("p(tree1) = " + ccd.getProbabilityOfTree(trees.get(0)));
        System.out.println("p(tree2) = " + ccd.getProbabilityOfTree(trees.get(1)));
        System.out.println("p(tree3) = " + ccd.getProbabilityOfTree(expandedTrees.get(0)));

        for (int i = 0; i < kregTrees.size(); i++) {
            System.out.println(ccd.getProbabilityOfTree(kregTrees.get(i)));
        }
    }
}
