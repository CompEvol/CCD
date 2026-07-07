package test.ccd.algorithms;

import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import ccd.algorithms.NNICladeExpansion;
import ccd.algorithms.NNICladeExpansion.NovelClade;
import ccd.algorithms.NNICladeExpansion.PairingMode;
import ccd.algorithms.NNICladeExpansion.Provenance;
import ccd.algorithms.regularisation.CCDExpansion;
import ccd.algorithms.regularisation.NNICladeExpander;
import ccd.model.CCD0;
import ccd.model.Clade;
import ccd.model.CladePartition;
import ccd.model.NNIRegCCD;
import ccd.model.bitsets.BitSet;
import org.junit.jupiter.api.Test;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Tests for {@link NNICladeExpansion} against hand-computed NNI neighbourhoods.
 */
public class NNICladeExpansionTest {

    private static final List<String> TAXA5 = Arrays.asList("A", "B", "C", "D", "E");

    private static Tree parse(String newick) {
        return new TreeParser(TAXA5, newick, 1, false);
    }

    /** Taxa of a novel clade as a sorted letter string, e.g. "ABD". */
    private static String name(BitSet taxa, Tree tree) {
        StringBuilder sb = new StringBuilder();
        String[] names = tree.getTaxaNames();
        for (int i = taxa.nextSetBit(0); i >= 0; i = taxa.nextSetBit(i + 1)) {
            sb.append(names[i]);
        }
        return sb.toString();
    }

    private static Set<String> novelNames(NNICladeExpansion exp, Tree tree) {
        Set<String> names = new TreeSet<>();
        for (NovelClade novel : exp.getNovelClades()) {
            names.add(name(novel.taxa, tree));
        }
        return names;
    }

    /**
     * Single tree (((A,B),C),D). The NNI neighbourhood of its three internal
     * edges produces exactly the novel clades {ABD, CD, AC, BC}; with one tree
     * the graph-wide and co-occurring modes must coincide.
     */
    @Test
    public void testSingleTreeExactNeighbourhood() {
        List<Tree> trees = new ArrayList<>();
        Tree tree = parse("(((A:1,B:1):1,C:2):1,D:3):0;");
        trees.add(tree);
        CCD0 ccd = new CCD0(trees, 0.0);

        NNICladeExpansion broad = new NNICladeExpansion(ccd, PairingMode.GRAPH_WIDE).compute();
        NNICladeExpansion strict = new NNICladeExpansion(ccd, PairingMode.CO_OCCURRING).compute();

        Set<String> expected = new TreeSet<>(Arrays.asList("ABD", "CD", "AC", "BC"));
        assertEquals(expected, novelNames(broad, tree), "graph-wide novel clades");
        assertEquals(expected, novelNames(strict, tree), "co-occurring novel clades");
    }

    /**
     * Every novel clade's NNI seed split is {kept, sibling}, and both of those
     * are clades that already exist in C - which is what makes the Phase 2
     * extension feasible. Check this on the single-tree case.
     */
    @Test
    public void testSeedSplitConstituentsExist() {
        List<Tree> trees = new ArrayList<>();
        trees.add(parse("(((A:1,B:1):1,C:2):1,D:3):0;"));
        CCD0 ccd = new CCD0(trees, 0.0);

        NNICladeExpansion exp = new NNICladeExpansion(ccd, PairingMode.GRAPH_WIDE).compute();
        for (NovelClade novel : exp.getNovelClades()) {
            assertTrue(!novel.provenances.isEmpty(), "novel clade has a provenance");
            for (Provenance prov : novel.provenances) {
                assertTrue(ccd.getClades().contains(prov.kept()), "kept piece is in C");
                assertTrue(ccd.getClades().contains(prov.sibling()), "sibling is in C");
                // novel taxa == kept ∪ sibling
                BitSet union = (BitSet) prov.kept().getCladeInBitsTaxaOnly().clone();
                union.or(prov.sibling().getCladeInBitsTaxaOnly());
                assertEquals(novel.taxa, union, "novel taxa equal kept ∪ sibling");
            }
        }
    }

    /**
     * Two 5-taxon trees that share the clade {A,B,C} but split it differently
     * ({AB,C} under {ABCD} in tree 1; {AC,B} under {ABCE} in tree 2). Graph-wide
     * pairing recombines split contexts that never co-occurred, so it strictly
     * exceeds the co-occurring (true sample-NNI) neighbourhood.
     */
    @Test
    public void testGraphWideExceedsCoOccurring() {
        List<Tree> trees = new ArrayList<>();
        Tree t1 = parse("((((A:1,B:1):1,C:1):1,D:1):1,E:1):0;");
        Tree t2 = parse("((((A:1,C:1):1,B:1):1,E:1):1,D:1):0;");
        trees.add(t1);
        trees.add(t2);
        CCD0 ccd = new CCD0(trees, 0.0);

        NNICladeExpansion broad = new NNICladeExpansion(ccd, PairingMode.GRAPH_WIDE).compute();
        NNICladeExpansion strict = new NNICladeExpansion(ccd, PairingMode.CO_OCCURRING).compute();

        Set<String> broadNames = novelNames(broad, t1);
        Set<String> strictNames = novelNames(strict, t1);

        Set<String> expectedBroad = new TreeSet<>(Arrays.asList(
                "ABD", "CD", "ACD", "BD", "ABE", "CE", "ACE", "BE", "DE", "BC"));
        Set<String> expectedStrict = new TreeSet<>(Arrays.asList(
                "DE", "ABD", "CD", "BC", "ACE", "BE"));
        assertEquals(expectedBroad, broadNames, "graph-wide novel clades");
        assertEquals(expectedStrict, strictNames, "co-occurring novel clades");

        // strict neighbourhood is contained in the graph-wide one
        assertTrue(broadNames.containsAll(strictNames), "strict ⊆ broad");

        // the cross-tree recombinations are exactly the graph-wide extras
        Set<String> onlyBroad = new TreeSet<>(broadNames);
        onlyBroad.removeAll(strictNames);
        assertEquals(new TreeSet<>(Arrays.asList("ACD", "BD", "ABE", "CE")), onlyBroad,
                "graph-wide-only (cross-tree) clades");
    }

    /**
     * The read-only novel-tree count must equal what the actual graph mutation
     * (NNICladeExpander) produces for the same expansion. Checked on the
     * single-tree neighbourhood.
     */
    @Test
    public void testNovelTreeCountMatchesMutationSingleTree() {
        Tree t = parse("(((A:1,B:1):1,C:2):1,D:3):0;");

        List<Tree> trees1 = new ArrayList<>();
        trees1.add(t);
        CCD0 readOnly = new CCD0(trees1, 0.0);
        BigInteger treesBefore = readOnly.getNumberOfTrees();
        NNICladeExpansion exp = new NNICladeExpansion(readOnly, PairingMode.GRAPH_WIDE);
        BigInteger predictedAfter = exp.getNumberOfTreesAfterExpansion();
        BigInteger predictedNovel = exp.getNumberOfNovelTrees();

        List<Tree> trees2 = new ArrayList<>();
        trees2.add(parse("(((A:1,B:1):1,C:2):1,D:3):0;"));
        CCD0 mutated = new CCD0(trees2, 0.0);
        new NNICladeExpander(PairingMode.GRAPH_WIDE).expandClades(mutated);
        BigInteger actualAfter = mutated.getNumberOfTrees();

        assertEquals(actualAfter, predictedAfter, "predicted post-expansion tree count");
        assertEquals(actualAfter.subtract(treesBefore), predictedNovel, "predicted novel tree count");
        // CCD itself must remain unchanged by the read-only count
        assertEquals(treesBefore, readOnly.getNumberOfTrees(), "read-only count must not mutate CCD");
    }

    /**
     * Same check for the two-tree case in both pairing modes. Graph-wide and
     * co-occurring expansions add different sets of clades, so the predicted
     * tree counts must follow.
     */
    @Test
    public void testNovelTreeCountMatchesMutationTwoTrees() {
        String n1 = "((((A:1,B:1):1,C:1):1,D:1):1,E:1):0;";
        String n2 = "((((A:1,C:1):1,B:1):1,E:1):1,D:1):0;";

        for (PairingMode mode : PairingMode.values()) {
            List<Tree> trees1 = new ArrayList<>();
            trees1.add(parse(n1));
            trees1.add(parse(n2));
            CCD0 readOnly = new CCD0(trees1, 0.0);
            BigInteger treesBefore = readOnly.getNumberOfTrees();
            NNICladeExpansion exp = new NNICladeExpansion(readOnly, mode);
            BigInteger predictedAfter = exp.getNumberOfTreesAfterExpansion();
            BigInteger predictedNovel = exp.getNumberOfNovelTrees();

            List<Tree> trees2 = new ArrayList<>();
            trees2.add(parse(n1));
            trees2.add(parse(n2));
            CCD0 mutated = new CCD0(trees2, 0.0);
            new NNICladeExpander(mode).expandClades(mutated);
            BigInteger actualAfter = mutated.getNumberOfTrees();

            assertEquals(actualAfter, predictedAfter,
                    "predicted post-expansion tree count [" + mode + "]");
            assertEquals(actualAfter.subtract(treesBefore), predictedNovel,
                    "predicted novel tree count [" + mode + "]");
            assertEquals(treesBefore, readOnly.getNumberOfTrees(),
                    "read-only count must not mutate CCD [" + mode + "]");
        }
    }

    /**
     * Post-(NNI + split-expansion) count must equal what running both passes
     * for real (NNICladeExpander + CCDExpansion) produces. Two-tree case in
     * both pairing modes.
     */
    @Test
    public void testTreeCountAfterFullExpansionMatchesMutation() {
        String n1 = "((((A:1,B:1):1,C:1):1,D:1):1,E:1):0;";
        String n2 = "((((A:1,C:1):1,B:1):1,E:1):1,D:1):0;";

        for (PairingMode mode : PairingMode.values()) {
            List<Tree> trees1 = new ArrayList<>();
            trees1.add(parse(n1));
            trees1.add(parse(n2));
            CCD0 readOnly = new CCD0(trees1, 0.0);
            BigInteger treesBefore = readOnly.getNumberOfTrees();
            NNICladeExpansion exp = new NNICladeExpansion(readOnly, mode);
            BigInteger predictedFull = exp.getNumberOfTreesAfterFullExpansion();
            BigInteger predictedNovelFull = exp.getNumberOfNovelTreesAfterFullExpansion();

            List<Tree> trees2 = new ArrayList<>();
            trees2.add(parse(n1));
            trees2.add(parse(n2));
            CCD0 mutated = new CCD0(trees2, 0.0);
            new NNICladeExpander(mode).expandClades(mutated);
            new CCDExpansion().expandCCD(mutated);
            BigInteger actualFull = mutated.getNumberOfTrees();

            assertEquals(actualFull, predictedFull,
                    "predicted post-full-expansion tree count [" + mode + "]");
            assertEquals(actualFull.subtract(treesBefore), predictedNovelFull,
                    "predicted novel-after-full tree count [" + mode + "]");
            // NNI-only count must be <= full count, and full must be >= before
            assertTrue(predictedFull.compareTo(exp.getNumberOfTreesAfterExpansion()) >= 0,
                    "full expansion >= NNI-only [" + mode + "]");
            assertTrue(predictedFull.compareTo(treesBefore) >= 0,
                    "full expansion >= before [" + mode + "]");
            // CCD itself must remain unchanged by the read-only count
            assertEquals(treesBefore, readOnly.getNumberOfTrees(),
                    "read-only count must not mutate CCD [" + mode + "]");
        }
    }

    /**
     * Probability mass on NNI-expanded trees must equal what a graph walk over
     * the fully-built {@link NNIRegCCD} (NNI + CCDExpansion + regularisation)
     * reports. Two-tree case in both pairing modes.
     */
    @Test
    public void testProbabilityOfNNIExpandedTreesMatchesGraphWalk() {
        String n1 = "((((A:1,B:1):1,C:1):1,D:1):1,E:1):0;";
        String n2 = "((((A:1,C:1):1,B:1):1,E:1):1,D:1):0;";
        double alpha = 0.4;

        for (PairingMode mode : PairingMode.values()) {
            List<Tree> trees1 = new ArrayList<>();
            trees1.add(parse(n1));
            trees1.add(parse(n2));
            CCD0 readOnly = new CCD0(trees1, 0.0);
            double predicted = new NNICladeExpansion(readOnly, mode)
                    .getProbabilityOfNNIExpandedTrees(alpha);

            List<Tree> trees2 = new ArrayList<>();
            trees2.add(parse(n1));
            trees2.add(parse(n2));
            NNIRegCCD reg = new NNIRegCCD(trees2, 0.0, mode, alpha);
            double actual = 1.0 - existingOnlyMass(reg.getRootClade(), new HashMap<>());

            assertEquals(actual, predicted, 1e-12,
                    "P(NNI-expanded trees) [" + mode + "]");
            // CCD itself must remain unchanged
            assertEquals(2, readOnly.getNumberOfBaseTrees(),
                    "read-only must not mutate CCD [" + mode + "]");
        }
    }

    /** Total CCP-mass of trees rooted at c that use only existing clades. */
    private static double existingOnlyMass(Clade c, Map<Clade, Double> memo) {
        if (c.isNNIExpanded()) return 0.0;
        if (c.isLeaf()) return 1.0;
        Double cached = memo.get(c);
        if (cached != null) return cached;
        double m = 0.0;
        for (CladePartition p : c.getPartitions()) {
            Clade l = p.getChildClades()[0];
            Clade r = p.getChildClades()[1];
            if (l.isNNIExpanded() || r.isNNIExpanded()) continue;
            m += p.getCCP() * existingOnlyMass(l, memo) * existingOnlyMass(r, memo);
        }
        memo.put(c, m);
        return m;
    }

    /**
     * Co-occurring mode needs all base trees; a CCD built without storing them
     * should fail loudly rather than silently producing wrong results.
     */
    @Test
    public void testCoOccurringRequiresBaseTrees() {
        // empty CCD with storeBaseTrees = false: only the first tree is kept,
        // so getBaseTrees() returns null and co-occurring mode cannot run
        CCD0 ccd = new CCD0(5, false);
        ccd.addTree(parse("((((A:1,B:1):1,C:1):1,D:1):1,E:1):0;"));
        ccd.addTree(parse("((((A:1,C:1):1,B:1):1,E:1):1,D:1):0;"));
        ccd.initialize();

        boolean threw = false;
        try {
            new NNICladeExpansion(ccd, PairingMode.CO_OCCURRING).compute();
        } catch (IllegalStateException expected) {
            threw = true;
        }
        // graph-wide still works without stored base trees
        NNICladeExpansion broad = new NNICladeExpansion(ccd, PairingMode.GRAPH_WIDE).compute();
        assertTrue(broad.getNumberOfNovelClades() >= 0);
        assertTrue(threw, "co-occurring mode should require stored base trees");
    }
}
