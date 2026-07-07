package test.ccd.model;

import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beastfx.app.treeannotator.TreeAnnotator;
import ccd.model.*;
import org.junit.jupiter.api.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class CCD0Test {

    private static final List<String> TAXA = Arrays.asList("A", "B", "C", "D");

    private static Tree parse(String newick) {
        return new TreeParser(TAXA, newick, 1, false);
    }

    private static List<Tree> sampledTreeListNonSA() {
        List<Tree> trees = new ArrayList<>();
        trees.add(parse("(((A:1,B:1):1,C:2):1,D:3):0;"));
        trees.add(parse("(((D:1,C:1):1,B:2):1,A:3):0;"));
        return trees;
    }

    private static List<Tree> sampledTreeListSA() {
        List<Tree> trees = new ArrayList<>();
        trees.add(parse("(((A:1,B:1):1,C:2):1,D:0):0;"));
        trees.add(parse("(((D:1,C:1):1,B:2):1,A:0):0;"));
        // trees.add(parse("((A:1,C:1):1,(B:1,D:1):1):0;"));
        // trees.add(parse("(((D:1,C:1):1,B:2):1,A:1):0;"));
        return trees;
    }

    private static TreeAnnotator.MemoryFriendlyTreeSet sampledTreeSetNonSA() throws IOException {
        String dataPath = "/Users/zyan598/Documents/GitHub/CCD_sampled_ancestor/example_trees/ccd0_nonSA.trees";
        TreeAnnotator.MemoryFriendlyTreeSet treeSet = new TreeAnnotator().new MemoryFriendlyTreeSet(dataPath, 0);
        return treeSet;
    }

    private static List<Tree> unsampledTreeList() {
        List<Tree> trees = new ArrayList<>();
        trees.add(parse("((A:1,B:1):1,(C:1,D:1):1):0;"));
        return trees;
    }

    @Test
    public void testNonSA() throws IOException {
        List<Tree> treeList = sampledTreeListNonSA();
        List<Tree> unsampledTreeList = unsampledTreeList();
        CCD0 ccd = new CCD0(treeList, 0);
        for (Clade clade : ccd.getClades()) {
            System.out.println(clade);
            for (CladePartition p : clade.getPartitions()) {
                System.out.println(p);
            }
        }
        double p1 = ccd.getProbabilityOfTree(unsampledTreeList.get(0));
        System.out.println("unsampled tree prob = " + p1);
    }

    @Test
    public void testCCD0CP() throws IOException {
        List<Tree> treeList = sampledTreeListSA();
        List<Tree> unsampledTreeList = unsampledTreeList();
        CCD0CP ccd = new CCD0CP(treeList, 0);
        for (Clade clade : ccd.getClades()) {
            System.out.println(clade);
            for (CladePartition p : clade.getPartitions()) {
                System.out.println(p);
            }
        }
        double p1 = ccd.getProbabilityOfTree(unsampledTreeList.get(0));
        System.out.println("p1 = " + p1);
    }

    @Test
    public void testCCD0SJ() throws IOException {
        List<Tree> treeList = sampledTreeListSA();
        List<Tree> unsampledTreeList = unsampledTreeList();
        CCD0SJ ccd = new CCD0SJ(treeList, 0);
        for (Clade clade : ccd.getClades()) {
            System.out.println(clade);
            for (CladePartition p : clade.getPartitions()) {
                System.out.println(p);
            }
        }
        double p1 = ccd.getProbabilityOfTree(unsampledTreeList.get(0));
        System.out.println("p1 = " + p1);
    }

    // @Test
    // public void testCCD0CP2() throws IOException {
    //     List<Tree> treeList = sampledTreeList2();
    //     List<Tree> unsampledTreeList = unsampledTreeList();
    //     TreeAnnotator.MemoryFriendlyTreeSet treeSet = sampledTreeSet();
    //     CCD0CP ccd = new CCD0CP(treeSet);
    //     for (Clade clade : ccd.getClades()) {
    //         System.out.println(clade);
    //         for (CladePartition p : clade.getPartitions()) {
    //             System.out.println(p);
    //         }
    //     }
    //     double p1 = ccd.getProbabilityOfTree(unsampledTreeList.get(0));
    //     double p2 = ccd.getProbabilityOfTree(treeList.get(0));
    //     double p3 = ccd.getProbabilityOfTree(treeList.get(1));
    //     System.out.println("sampled tree 1 = " + p2);
    //     System.out.println("sampled tree 2 = " + p3);
    //     System.out.println("unsampled tree 1 = " + p1);
    // }
    //
    // @Test
    // public void testCCD0SJ() throws IOException {
    //     List<Tree> treeList = sampledTreeList();
    //     TreeAnnotator.MemoryFriendlyTreeSet treeSet = sampledTreeSet();
    //     CCD0SJ ccd = new CCD0SJ(treeSet);
    //     double p1 = ccd.getProbabilityOfTree(treeList.get(0));
    //     System.out.println("p1 = " + p1);
    // }

    /**
     * The {@link ccd.model.AbstractCCD#AbstractCCD(List, double)} constructor
     * has two branches; the {@code burnin == 0} branch must still set
     * {@code numBaseTrees} so {@link Clade#getCladeCredibility()} is finite.
     * Before the fix, {@code burnin == 0} left {@code numBaseTrees = 0}, so
     * {@code numOccurrences / 0 = Infinity} for every clade and the CCD0
     * constructor blew up with NaN inside {@code setPartitionProbabilities}.
     */
    @Test
    public void testZeroBurninSetsBaseTreeCount() {
        List<Tree> trees = new ArrayList<>();
        trees.add(parse("((A:1,B:1):1,(C:1,D:1):1):0;"));
        trees.add(parse("((A:1,C:1):1,(B:1,D:1):1):0;"));

        CCD0 ccd = new CCD0(trees, 0.0);

        assertEquals(trees.size(), ccd.getNumberOfBaseTrees());
        for (Clade c : ccd.getClades()) {
            double cred = c.getCladeCredibility();
            assertTrue(Double.isFinite(cred), "clade credibility must be finite, got " + cred + " for " + c);
        }
    }

    /**
     * Drives CCD0.setPartitionProbabilities into its log-sum-exp underflow
     * branch (CCD0.java:860 "if (sumSubtreeProbabilities == 0)") and checks
     * that the resulting partition CCPs are normalised. With the sign error
     * at CCD0.java:894 (logSum = logMax - log(intermSum), should be +) the
     * CCPs of the root's two partitions come out as 2.0 each instead of 0.5,
     * so they sum to 4.0 rather than 1.0.
     * <p>
     * Setup: two 4-taxon trees that disagree on the root split, giving the
     * root clade two partitions: {A,B}|{C,D} and {A,C}|{B,D}. Setting every
     * non-leaf clade's parameter to 1e-200 makes each partition's product
     * (1e-200 * 1e-200 = 1e-400) underflow to exactly 0, so the sum is 0 and
     * the underflow branch fires.
     */
    @Test
    public void testRootPartitionCCPsNormaliseUnderUnderflow() {
        List<Tree> trees = new ArrayList<>();
        trees.add(parse("((A:1,B:1):1,(C:1,D:1):1):0;"));
        trees.add(parse("((A:1,C:1):1,(B:1,D:1):1):0;"));

        CCD0 ccd = new CCD0(trees, 0.0);

        for (Clade clade : ccd.getClades()) {
            if (!clade.isLeaf()) {
                clade.setCladeParameter(1e-200);
            }
        }
        ccd.resetSumCladeCredibilities();
        ccd.setPartitionProbabilities(ccd.getRootClade(), true);

        double ccpSum = 0.0;
        for (CladePartition p : ccd.getRootClade().getPartitions()) {
            ccpSum += p.getCCP();
        }

        assertEquals(1.0, ccpSum, 1e-9);
    }
}
