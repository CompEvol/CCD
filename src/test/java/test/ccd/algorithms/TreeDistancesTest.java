package test.ccd.algorithms;

import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import ccd.algorithms.TreeDistances;
import ccd.model.WrappedBeastTreeWithSampledAncestor;
import org.junit.jupiter.api.Test;

import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Tests for {@link TreeDistances#sampledAncestorRobinsonFoulds}.
 *
 * <p>All tests use 4-taxon trees on {A, B, C, D}. Branch length 0 marks a
 * sampled ancestor (SA). The cuttable clade set of an SA tree is the set of
 * taxa-only internal-node clades plus the singleton {i} for every non-SA
 * leaf i; SA-RF is the size of the symmetric difference between the two
 * trees' cuttable clade sets.
 */
public class TreeDistancesTest {

    private static final List<String> TAXA = Arrays.asList("A", "B", "C", "D");

    private static WrappedBeastTreeWithSampledAncestor wrap(String newick) {
        Tree t = new TreeParser(TAXA, newick, 1, false);
        return new WrappedBeastTreeWithSampledAncestor(t);
    }

    /* ---------- Identity and symmetry ---------- */

    @Test
    public void testIdentityNoSA() {
        WrappedBeastTreeWithSampledAncestor t = wrap("((A:1,B:1):1,(C:1,D:1):1):0;");
        assertEquals(0, TreeDistances.sampledAncestorRobinsonFoulds(t, t));
    }

    @Test
    public void testIdentityWithSA() {
        WrappedBeastTreeWithSampledAncestor t = wrap("((A:1,B:0):1,(C:1,D:1):1):0;");
        assertEquals(0, TreeDistances.sampledAncestorRobinsonFoulds(t, t));
    }

    @Test
    public void testIdenticalTreesDistinctInstances() {
        WrappedBeastTreeWithSampledAncestor t1 = wrap("((A:1,B:0):1,(C:1,D:1):1):0;");
        WrappedBeastTreeWithSampledAncestor t2 = wrap("((A:1,B:0):1,(C:1,D:1):1):0;");
        assertEquals(0, TreeDistances.sampledAncestorRobinsonFoulds(t1, t2));
    }

    @Test
    public void testSymmetry() {
        WrappedBeastTreeWithSampledAncestor t1 = wrap("((A:1,B:0):1,(C:1,D:1):1):0;");
        WrappedBeastTreeWithSampledAncestor t2 = wrap("((A:1,C:1):1,(B:1,D:0):1):0;");
        int d12 = TreeDistances.sampledAncestorRobinsonFoulds(t1, t2);
        int d21 = TreeDistances.sampledAncestorRobinsonFoulds(t2, t1);
        assertEquals(d12, d21, "SA-RF must be symmetric");
    }

    /* ---------- SA-only differences (same topology) ---------- */

    @Test
    public void testSingleSAFlipContributesOne() {
        // Same topology, A is SA only in T1.
        // T1 cuttable: internals {A,B},{C,D},{ABCD}; singletons {B},{C},{D}.
        // T2 cuttable: same internals; singletons {A},{B},{C},{D}.
        // Symm diff: {{A}} (only in T2). Distance = 1.
        WrappedBeastTreeWithSampledAncestor t1 = wrap("((A:0,B:1):1,(C:1,D:1):1):0;");
        WrappedBeastTreeWithSampledAncestor t2 = wrap("((A:1,B:1):1,(C:1,D:1):1):0;");
        assertEquals(1, TreeDistances.sampledAncestorRobinsonFoulds(t1, t2));
    }

    @Test
    public void testSwappedSAOnSameTopologyContributesTwo() {
        // Same topology; A is SA in T1, B is SA in T2.
        // T1 singletons {B,C,D}; T2 singletons {A,C,D}. Symm diff = {{A},{B}}.
        // Internals shared. Distance = 2.
        WrappedBeastTreeWithSampledAncestor t1 = wrap("((A:0,B:1):1,(C:1,D:1):1):0;");
        WrappedBeastTreeWithSampledAncestor t2 = wrap("((A:1,B:0):1,(C:1,D:1):1):0;");
        assertEquals(2, TreeDistances.sampledAncestorRobinsonFoulds(t1, t2));
    }

    /* ---------- Topology-only differences (no SAs) ---------- */

    @Test
    public void testTopologyOnlyDifference() {
        // No SAs in either tree, so singletons coincide and only internal
        // clades contribute. T1 internals {A,B},{C,D},{ABCD};
        // T2 internals {A,C},{B,D},{ABCD}. Symm diff size = 4.
        WrappedBeastTreeWithSampledAncestor t1 = wrap("((A:1,B:1):1,(C:1,D:1):1):0;");
        WrappedBeastTreeWithSampledAncestor t2 = wrap("((A:1,C:1):1,(B:1,D:1):1):0;");
        assertEquals(4, TreeDistances.sampledAncestorRobinsonFoulds(t1, t2));
    }

    /* ---------- Combined topology + SA differences ---------- */

    @Test
    public void testCombinedTopologyAndSADifference() {
        // T1 = ((A,B),(C,D)) with A,C SA  (one SA per cherry, valid).
        //   Internals {A,B},{C,D},{ABCD}; singletons {B},{D}.
        //   Cuttable = {{A,B},{C,D},{ABCD},{B},{D}}.
        // T2 = ((A,D),(B,C)) with B,D SA  (one SA per cherry, valid).
        //   Internals {A,D},{B,C},{ABCD}; singletons {A},{C}.
        //   Cuttable = {{A,D},{B,C},{ABCD},{A},{C}}.
        // Intersection = {{ABCD}}.
        // Symm diff size = 4 + 4 = 8.
        WrappedBeastTreeWithSampledAncestor t1 = wrap("((A:0,B:1):1,(C:0,D:1):1):0;");
        WrappedBeastTreeWithSampledAncestor t2 = wrap("((A:1,D:0):1,(B:0,C:1):1):0;");
        assertEquals(8, TreeDistances.sampledAncestorRobinsonFoulds(t1, t2));
    }

    /* ---------- Reduces to standard RF on SA-free trees ---------- */

    @Test
    public void testReducesToFullRFWhenNoSAs() {
        // For two SA-free trees, SA-RF should equal 2 * (the existing
        // halved robinsonsFouldDistance), because singletons are shared and
        // internal clades contribute symmetrically.
        WrappedBeastTreeWithSampledAncestor t1 = wrap("((A:1,B:1):1,(C:1,D:1):1):0;");
        WrappedBeastTreeWithSampledAncestor t2 = wrap("((A:1,C:1):1,(B:1,D:1):1):0;");
        int saRF = TreeDistances.sampledAncestorRobinsonFoulds(t1, t2);
        int halvedRF = TreeDistances.robinsonsFouldDistance(t1, t2);
        assertEquals(2 * halvedRF, saRF,
                "SA-RF on SA-free trees must equal full (unhalved) standard RF");
    }

    /* ---------- Pure topology RF (taxa-only, halved convention) ---------- */

    @Test
    public void testTopologyRFIgnoresSADifferences() {
        // Same topology, different SA assignments → topology RF = 0.
        WrappedBeastTreeWithSampledAncestor t1 = wrap("((A:0,B:1):1,(C:1,D:0):1):0;");
        WrappedBeastTreeWithSampledAncestor t2 = wrap("((A:1,B:0):1,(C:0,D:1):1):0;");
        assertEquals(0, TreeDistances.topologyRobinsonFoulds(t1, t2),
                "topology RF must ignore SA flag differences when topology is identical");
    }

    @Test
    public void testTopologyRFOnDifferentTopologies() {
        // Different topologies, no SAs in either tree.
        // T1 internals (taxa-only): {A,B}, {C,D}, {ABCD}
        // T2 internals (taxa-only): {A,C}, {B,D}, {ABCD}
        // |T1 \ T2| = 2 (the two cherries differ).
        WrappedBeastTreeWithSampledAncestor t1 = wrap("((A:1,B:1):1,(C:1,D:1):1):0;");
        WrappedBeastTreeWithSampledAncestor t2 = wrap("((A:1,C:1):1,(B:1,D:1):1):0;");
        assertEquals(2, TreeDistances.topologyRobinsonFoulds(t1, t2));
    }

    @Test
    public void testTopologyRFOnDifferentTopologiesWithSAs() {
        // Different topologies AND different SAs; the topology RF should be
        // the same as if no SAs were present (it strips them).
        WrappedBeastTreeWithSampledAncestor t1 = wrap("((A:0,B:1):1,(C:0,D:1):1):0;");
        WrappedBeastTreeWithSampledAncestor t2 = wrap("((A:1,D:0):1,(B:0,C:1):1):0;");
        // T1 internals (taxa-only): {A,B}, {C,D}, {ABCD}
        // T2 internals (taxa-only): {A,D}, {B,C}, {ABCD}
        // |T1 \ T2| = 2.
        assertEquals(2, TreeDistances.topologyRobinsonFoulds(t1, t2));
    }

    @Test
    public void testTopologyRFAgreesWithStandardRFOnSAFreeTrees() {
        // When no SAs are present the wide bitsets of WrappedBeastTreeWithSampledAncestor
        // collapse to taxa-only, so robinsonsFouldDistance and topologyRobinsonFoulds
        // must agree.
        WrappedBeastTreeWithSampledAncestor t1 = wrap("((A:1,B:1):1,(C:1,D:1):1):0;");
        WrappedBeastTreeWithSampledAncestor t2 = wrap("((A:1,C:1):1,(B:1,D:1):1):0;");
        assertEquals(TreeDistances.robinsonsFouldDistance(t1, t2),
                TreeDistances.topologyRobinsonFoulds(t1, t2),
                "topologyRF must equal standard RF on SA-free trees");
    }

    @Test
    public void testTopologyRFSymmetryOnBinaryTrees() {
        // For binary trees on the same taxa, |T1 \ T2| = |T2 \ T1|.
        WrappedBeastTreeWithSampledAncestor t1 = wrap("((A:0,B:1):1,(C:1,D:1):1):0;");
        WrappedBeastTreeWithSampledAncestor t2 = wrap("((A:1,D:0):1,(B:1,C:1):1):0;");
        assertEquals(TreeDistances.topologyRobinsonFoulds(t1, t2),
                TreeDistances.topologyRobinsonFoulds(t2, t1),
                "topology RF must be symmetric on binary trees");
    }
}
