package ccd.model;

import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import org.junit.jupiter.api.Test;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * {@link ITreeDistribution#getNumberOfTrees()} must report the size of the <em>support</em>, the
 * number of topologies with non-zero probability, and so must agree with
 * {@link ITreeDistribution#containsTree(Tree)}.
 *
 * <p>Uses the same four-taxon worked example as {@link UniformEscapeCCDTest}: two sampled trees,
 * a backbone covering only those two, and 15 rooted topologies in total. A full-support model
 * must therefore report 15, not 2 -- the failure mode this test exists to catch is a full-support
 * model inheriting the graph-based count from {@link AbstractCCD}.
 */
public class NumberOfTreesContractTest {

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

    @Test
    public void ccd1ReportsItsGraph() {
        CCD1 ccd1 = new CCD1(sampleTreeList(), 0.0);
        assertEquals(BigInteger.TWO, ccd1.getNumberOfTrees(),
                "CCD1's support is exactly the topologies its graph represents");
    }

    @Test
    public void kRegCCDReportsFullSupport() {
        KRegCCD kreg = new KRegCCD(sampleTreeList(), 0.0, 0.05, 0.4, 2);
        assertEquals(BigInteger.valueOf(15), kreg.getNumberOfTrees(),
                "KRegCCD is full support, so all 15 rooted topologies on four taxa");
    }

    @Test
    public void mRegCCDReportsFullSupport() {
        MRegCCD mreg = new MRegCCD(sampleTreeList(), 0.0, 0.05);
        assertEquals(BigInteger.valueOf(15), mreg.getNumberOfTrees(),
                "MRegCCD is full support, so all 15 rooted topologies on four taxa");
    }

    /** The count and containsTree describe the same set, so they must not disagree. */
    @Test
    public void fullSupportCountAgreesWithContainsTree() {
        for (ITreeDistribution d : List.of(
                new KRegCCD(sampleTreeList(), 0.0, 0.05, 0.4, 2),
                new MRegCCD(sampleTreeList(), 0.0, 0.05),
                new UniformEscapeCCD(sampleTreeList(), 0.0, 0.05))) {
            boolean full = d.getNumberOfTrees()
                    .equals(AbstractCCD.numberOfRootedTopologies(d.getNumberOfLeaves()));
            assertTrue(full, d.getClass().getSimpleName()
                    + " claims every tree via containsTree, so must count every tree");
            // A topology absent from the two sampled trees, hence outside the backbone graph.
            assertTrue(d.containsTree(parse("((A:1,B:1):1,(C:1,D:1):1):0;")),
                    d.getClass().getSimpleName() + " should contain an off-graph topology");
        }
    }

    @Test
    public void numberOfRootedTopologiesMatchesDoubleFactorial() {
        assertEquals(BigInteger.ONE, AbstractCCD.numberOfRootedTopologies(1));
        assertEquals(BigInteger.ONE, AbstractCCD.numberOfRootedTopologies(2));
        assertEquals(BigInteger.valueOf(3), AbstractCCD.numberOfRootedTopologies(3));
        assertEquals(BigInteger.valueOf(15), AbstractCCD.numberOfRootedTopologies(4));
        assertEquals(BigInteger.valueOf(105), AbstractCCD.numberOfRootedTopologies(5));
    }

    /**
     * The reason {@link AbstractCCD#logBigInteger(BigInteger)} exists: the counts here overflow
     * {@code double}, so the naive {@code Math.log(v.doubleValue())} returns infinity.
     */
    @Test
    public void logBigIntegerSurvivesValuesBeyondDoubleRange() {
        // The naive conversion survives 129 taxa (840 bits) but not 160 (1093 bits): a double
        // overflows past ~1024 bits. Several data sets in routine use are on the wrong side of
        // that boundary, which is why the helper exists.
        BigInteger ok = AbstractCCD.numberOfRootedTopologies(129);
        assertEquals(840, ok.bitLength(), "129 taxa needs 840 bits");
        assertTrue(Double.isFinite(ok.doubleValue()), "129 taxa still fits a double");
        assertEquals(582.13, AbstractCCD.logBigInteger(ok), 0.01);

        BigInteger big = AbstractCCD.numberOfRootedTopologies(160);
        assertEquals(1093, big.bitLength(), "160 taxa needs 1093 bits");
        assertTrue(Double.isInfinite(big.doubleValue()), "the naive conversion overflows here");
        assertTrue(Double.isFinite(AbstractCCD.logBigInteger(big)), "the helper does not");

        // Exact on small values, and consistent with the identity log(a*b) = log a + log b.
        assertEquals(Math.log(15), AbstractCCD.logBigInteger(BigInteger.valueOf(15)), 1e-12);
        BigInteger a = AbstractCCD.numberOfRootedTopologies(60);
        assertEquals(AbstractCCD.logBigInteger(a) + AbstractCCD.logBigInteger(big),
                AbstractCCD.logBigInteger(a.multiply(big)), 1e-6);

        assertThrows(IllegalArgumentException.class,
                () -> AbstractCCD.logBigInteger(BigInteger.ZERO));
    }
}
