package test.ccd.model;

import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import ccd.model.UniformEscapeCCD;
import org.junit.jupiter.api.Test;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * The talk's worked example: two trees on four taxa, all 15 rooted topologies enumerated.
 *
 * <p>The plain CCD1 backbone built from {@code (((A,B),C),D)} and {@code (((D,C),B),A)} covers only
 * those two topologies (its root only ever splits as {@code ABC|D} or {@code BCD|A}), each with
 * probability 1/2. {@link UniformEscapeCCD} adds a lump of escape mass {@code epsilon} spread evenly
 * over the remaining {@code M = 15 - 2 = 13} topologies, so:
 * <pre>
 *   the 2 in-support trees:  (1 - epsilon) * 1/2
 *   the 13 off-support trees: epsilon / 13
 *   total:                    (1 - epsilon) + epsilon = 1.
 * </pre>
 */
public class UniformEscapeCCDTest {

    private static final List<String> TAXA = Arrays.asList("A", "B", "C", "D");
    private static final double EPS = 0.05;

    private static Tree parse(String newick) {
        return new TreeParser(TAXA, newick, 1, false);
    }

    /** The two sampled trees; the CCD1 backbone covers exactly these. */
    private static List<Tree> sampleTreeList() {
        List<Tree> trees = new ArrayList<>();
        trees.add(parse("(((A:1,B:1):1,C:1):1,D:1):0;"));
        trees.add(parse("(((D:1,C:1):1,B:1):1,A:1):0;"));
        return trees;
    }

    /** The other 13 rooted topologies on {A,B,C,D} (1 balanced split of the root + 12 others). */
    private static List<Tree> offSupportTrees() {
        List<Tree> trees = new ArrayList<>();
        trees.add(parse("((A:1,B:1):1,(C:1,D:1):1):0;"));     // root split AB|CD (unobserved)
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
    public void testUniformEscapeDistribution() {
        List<Tree> sampled = sampleTreeList();
        List<Tree> offSupport = offSupportTrees();
        UniformEscapeCCD ccd = new UniformEscapeCCD(sampled, 0.0, EPS);

        // The backbone covers exactly the 2 sampled topologies; the full model covers all 15.
        assertEquals(BigInteger.valueOf(15), ccd.getNumberOfTrees(), "full support over all 15 rooted topologies");
        assertEquals(BigInteger.valueOf(13), ccd.getNumberOfOffSupportTrees(), "M = 15 - 2 off-support topologies");

        // In-support trees: (1 - eps) * P_CCD1 = (1 - eps) * 1/2 each.
        double inSupport = (1.0 - EPS) * 0.5;
        for (Tree t : sampled) {
            assertEquals(inSupport, ccd.getProbabilityOfTree(t), 1e-12, "in-support probability");
            assertEquals(Math.log(inSupport), ccd.getLogProbabilityOfTree(t), 1e-9, "in-support log-probability");
            assertTrue(ccd.containsTree(t), "in-support tree is contained");
        }

        // Off-support trees: each gets eps / 13.
        double offProb = EPS / 13.0;
        double escapeLogProb = Math.log(EPS) - Math.log(13.0);
        assertEquals(escapeLogProb, ccd.getEscapeLogProbability(), 1e-12, "escape log-probability = log(eps/13)");
        for (Tree t : offSupport) {
            assertEquals(offProb, ccd.getProbabilityOfTree(t), 1e-12, "off-support probability = eps/13");
            assertEquals(escapeLogProb, ccd.getLogProbabilityOfTree(t), 1e-9, "off-support log-probability");
            assertTrue(ccd.containsTree(t), "off-support tree is still contained (full support)");
        }

        // The whole distribution is exactly normalised over all 15 topologies.
        double total = 0.0;
        for (Tree t : sampled) total += ccd.getProbabilityOfTree(t);
        for (Tree t : offSupport) total += ccd.getProbabilityOfTree(t);
        assertEquals(1.0, total, 1e-12, "the 15 topology probabilities sum to 1");
    }

    @Test
    public void testZeroEscapeRecoversCCD1() {
        // With epsilon = 0 the model is exactly the CCD1 backbone: off-support trees have probability 0.
        UniformEscapeCCD ccd = new UniformEscapeCCD(sampleTreeList(), 0.0, 0.0);
        for (Tree t : sampleTreeList()) {
            assertEquals(0.5, ccd.getProbabilityOfTree(t), 1e-12, "backbone probability unchanged");
        }
        for (Tree t : offSupportTrees()) {
            assertEquals(0.0, ccd.getProbabilityOfTree(t), 1e-12, "no escape mass => off-support probability 0");
            assertTrue(Double.isInfinite(ccd.getLogProbabilityOfTree(t)), "off-support log-probability is -inf");
        }
    }

    @Test
    public void testNumberOfRootedTopologies() {
        assertEquals(BigInteger.ONE, UniformEscapeCCD.numberOfRootedTopologies(2), "(2*2-3)!! = 1");
        assertEquals(BigInteger.valueOf(3), UniformEscapeCCD.numberOfRootedTopologies(3), "3!! = 3");
        assertEquals(BigInteger.valueOf(15), UniformEscapeCCD.numberOfRootedTopologies(4), "5!! = 15");
        assertEquals(BigInteger.valueOf(105), UniformEscapeCCD.numberOfRootedTopologies(5), "7!! = 105");
        assertEquals(BigInteger.valueOf(945), UniformEscapeCCD.numberOfRootedTopologies(6), "9!! = 945");
    }
}
