package ccd.model;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beastfx.app.treeannotator.TreeAnnotator;
import ccd.model.bitsets.BitSet;
import ccd.tools.CCDToolUtil;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Where does a real held-out tree's probability actually go under {@link CRegCCD}?
 *
 * <p>Walks every internal node of every held-out tree, classifies its split into the four classes,
 * and tallies how many nodes fall in each class and how much log probability each class contributes.
 * This measures the concern that a class-4 split (neither child observed) reconnects to the CCD only
 * by chance: if class 4 is both rare and responsible for a large share of the total loss, the
 * uniform-within-class-4 prior is the binding weakness.
 *
 * <p>Also reports, for each class-4 node encountered, how many observed clades survive intact inside
 * the two novel children -- i.e. how much backbone a class-4 split destroys.
 */
public class ClassUsageAnalysis {

    private static final String PATH = System.getProperty("ccd.trees", "");
    private static final int N = Integer.parseInt(System.getProperty("ccd.n", "1000"));

    private static List<String> newickCache;
    private static List<String> taxaCache;

    private static void load() throws Exception {
        if (newickCache != null) {
            return;
        }
        TreeAnnotator.TreeSet ts = CCDToolUtil.getTreeSet(PATH, 10);
        ts.reset();
        List<String> nwk = new ArrayList<>();
        List<String> taxa = null;
        while (ts.hasNext()) {
            Tree t = ts.next();
            if (taxa == null) {
                String[] byNr = new String[t.getLeafNodeCount()];
                for (Node leaf : t.getExternalNodes()) {
                    byNr[leaf.getNr()] = leaf.getID();
                }
                taxa = new ArrayList<>(List.of(byNr));
            }
            nwk.add(t.getRoot().toNewick() + ";");
        }
        newickCache = nwk;
        taxaCache = taxa;
    }

    private static List<Tree> read(int count, double from, double to) throws Exception {
        load();
        List<String> pool = newickCache.subList((int) (from * newickCache.size()),
                (int) (to * newickCache.size()));
        List<Tree> out = new ArrayList<>();
        double step = Math.max(1.0, pool.size() / (double) count);
        for (int i = 0; i < count && (int) (i * step) < pool.size(); i++) {
            out.add(new TreeParser(taxaCache, pool.get((int) (i * step)), 0, false));
        }
        return out;
    }

    @Test
    public void classUsageOnHeldOutTrees() throws Exception {
        Assumptions.assumeTrue(!PATH.isEmpty() && new File(PATH).exists(),
                "set -Dccd.trees=/path/to/x.trees");

        CRegCCD ccd = new CRegCCD(read(N, 0.0, 0.5), 0.0,
                Double.parseDouble(System.getProperty("ccd.alpha", "0.4")),
                Double.parseDouble(System.getProperty("ccd.alpha1", "0.4")),
                Double.parseDouble(System.getProperty("ccd.alpha2", "0.05")));
        List<Tree> test = read(N, 0.5, 1.0);

        long[] nodes = new long[4];
        double[] logp = new double[4];
        long class4Nodes = 0;
        long survivingObserved = 0;
        long shatteredObserved = 0;

        for (Tree t : test) {
            Map<Node, BitSet> bits = new HashMap<>();
            computeBits(t.getRoot(), bits, ccd.getSizeOfLeavesArray());
            for (Node v : t.getNodesAsArray()) {
                if (v.isLeaf()) {
                    continue;
                }
                BitSet cb = bits.get(v);
                BitSet ab = bits.get(v.getChildren().get(0));
                BitSet bb = bits.get(v.getChildren().get(1));
                int cls = ccd.splitClass(cb, ab, bb);
                nodes[cls]++;
                logp[cls] += ccd.logSplitProbability(cb, ab, bb,
                        ccd.getAlpha(), ccd.getAlpha1(), ccd.getAlpha2());
                if (cls == 3) {
                    class4Nodes++;
                    // how much observed structure did this split preserve vs destroy?
                    for (Clade obs : ccd.getClades()) {
                        if (obs.size() < 2 || obs.size() >= cb.cardinality()) {
                            continue;
                        }
                        BitSet o = obs.getCladeInBits();
                        BitSet tmp = BitSet.newBitSet(o);
                        tmp.andNot(cb);
                        if (!tmp.isEmpty()) {
                            continue; // not inside this clade at all
                        }
                        if (subset(o, ab) || subset(o, bb)) {
                            survivingObserved++;
                        } else {
                            shatteredObserved++;
                        }
                    }
                }
            }
        }

        long totalNodes = nodes[0] + nodes[1] + nodes[2] + nodes[3];
        double totalLogp = logp[0] + logp[1] + logp[2] + logp[3];
        System.out.printf("%n=== %s: class usage over %d held-out trees (%s) ===%n",
                new File(PATH).getName(), test.size(), ccd);
        System.out.printf("%-28s %10s %8s %14s %10s %12s%n",
                "class", "nodes", "% nodes", "total logP", "% logP", "mean logP");
        String[] names = {"1 observed split", "2 both children obs.",
                "3 one child observed", "4 neither observed"};
        for (int j = 0; j < 4; j++) {
            System.out.printf("%-28s %10d %7.2f%% %14.1f %9.2f%% %12.3f%n",
                    names[j], nodes[j], 100.0 * nodes[j] / totalNodes, logp[j],
                    100.0 * logp[j] / totalLogp, nodes[j] == 0 ? 0 : logp[j] / nodes[j]);
        }
        System.out.printf("total mean logP per tree = %.2f%n", totalLogp / test.size());
        if (class4Nodes > 0) {
            long tot = survivingObserved + shatteredObserved;
            System.out.printf("class-4 splits: %d; observed clades below them: %d intact (%.1f%%), "
                            + "%d shattered (%.1f%%)%n",
                    class4Nodes, survivingObserved, 100.0 * survivingObserved / tot,
                    shatteredObserved, 100.0 * shatteredObserved / tot);
        } else {
            System.out.println("no class-4 splits occurred in any held-out tree");
        }
    }

    private static boolean subset(BitSet a, BitSet c) {
        BitSet tmp = BitSet.newBitSet(a);
        tmp.andNot(c);
        return tmp.isEmpty();
    }

    private static BitSet computeBits(Node v, Map<Node, BitSet> bits, int leafArraySize) {
        BitSet b = BitSet.newBitSet(leafArraySize);
        if (v.isLeaf()) {
            b.set(v.getNr());
        } else {
            b.or(computeBits(v.getChildren().get(0), bits, leafArraySize));
            b.or(computeBits(v.getChildren().get(1), bits, leafArraySize));
        }
        bits.put(v, b);
        return b;
    }


    /**
     * Size profile of the splits that introduce two novel clades: parent size and the sizes of the
     * two novel children. If the smaller side is consistently small, grading class 2 by the size of
     * the smaller side would concentrate its mass where the real novelty is, instead of on the
     * balanced splits that dominate a uniform draw.
     */
    @Test
    public void novelSplitSizeProfile() throws Exception {
        Assumptions.assumeTrue(!PATH.isEmpty() && new File(PATH).exists(),
                "set -Dccd.trees=/path/to/x.trees");
        CRegCCD ccd = new CRegCCD(read(N, 0.0, 0.5), 0.0,
                Double.parseDouble(System.getProperty("ccd.alpha", "0.4")),
                Double.parseDouble(System.getProperty("ccd.alpha1", "0.4")),
                Double.parseDouble(System.getProperty("ccd.alpha2", "0.05")));
        List<Tree> test = read(N, 0.5, 1.0);

        System.out.printf("%n=== %s: size profile of two-novel-clade splits ===%n",
                new File(PATH).getName());
        System.out.printf("%-8s %-10s %-10s %-14s %-16s%n",
                "parent m", "small side", "large side", "small/parent", "uniform E[small]");
        int count = 0;
        double sumFrac = 0.0;
        java.util.Map<Integer, Integer> smallSizes = new java.util.TreeMap<>();
        for (Tree t : test) {
            Map<Node, BitSet> bits = new HashMap<>();
            computeBits(t.getRoot(), bits, ccd.getSizeOfLeavesArray());
            for (Node v : t.getNodesAsArray()) {
                if (v.isLeaf()) {
                    continue;
                }
                BitSet cb = bits.get(v);
                BitSet ab = bits.get(v.getChildren().get(0));
                BitSet bb = bits.get(v.getChildren().get(1));
                if (ccd.splitClass(cb, ab, bb) != 3) {
                    continue;
                }
                int m = cb.cardinality();
                int small = Math.min(ab.cardinality(), bb.cardinality());
                int large = Math.max(ab.cardinality(), bb.cardinality());
                count++;
                sumFrac += small / (double) m;
                smallSizes.merge(small, 1, Integer::sum);
                if (count <= 40) {
                    System.out.printf("%-8d %-10d %-10d %-14.3f %-16.1f%n",
                            m, small, large, small / (double) m, m / 2.0);
                }
            }
        }
        if (count == 0) {
            System.out.println("no two-novel-clade splits in the held-out set");
            return;
        }
        System.out.printf("%d such splits; mean smaller-side fraction = %.3f "
                        + "(a uniform bipartition would give ~0.5)%n", count, sumFrac / count);
        System.out.println("distribution of smaller-side size: " + smallSizes);
    }
}
