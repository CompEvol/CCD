package ccd.model;

import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import ccd.model.bitsets.BitSet;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * Exploratory: sizes of the four split classes at a clade, for the class-based smoothing proposal.
 *
 * <p>At a clade {@code C} of {@code m} taxa every one of the {@code 2^(m-1) - 1} bipartitions falls
 * into exactly one of: (1) observed split; (2) unobserved, both children observed clades
 * (the CCD0 expansion); (3) unobserved, exactly one child observed; (4) unobserved, neither child
 * observed. Classes 1-3 are at most polynomial in the number of observed clades; class 4 is
 * essentially all of {@code 2^(m-1)}. This prints the four sizes and the probability the smoothed
 * model would put on class 1, under a per-split constant pseudocount versus a per-class total.
 */
public class SplitClassSizeAnalysis {

    private static List<Tree> randomTrees(int nTaxa, int nTrees, long seed) {
        List<String> taxa = new ArrayList<>();
        for (int i = 0; i < nTaxa; i++) {
            taxa.add("T" + i);
        }
        Random rng = new Random(seed);
        List<Tree> out = new ArrayList<>();
        for (int t = 0; t < nTrees; t++) {
            List<String> pool = new ArrayList<>(taxa);
            while (pool.size() > 1) {
                int i = rng.nextInt(pool.size());
                String a = pool.remove(i);
                int j = rng.nextInt(pool.size());
                String b = pool.remove(j);
                pool.add("(" + a + "," + b + ")");
            }
            out.add(new TreeParser(taxa, pool.get(0) + ";", 1, false));
        }
        return out;
    }

    @Test
    public void classSizesAtRoot() {
        System.out.printf("%-6s %-7s %-8s %-8s %-8s %-14s %-14s %-14s%n",
                "taxa", "trees", "|A1|", "|A2|", "|A3|", "|A4|", "P(A1) per-split", "P(A1) per-class");
        for (int nTaxa : new int[]{12, 20, 30, 40}) {
            int nTrees = 1000;
            List<Tree> trees = randomTrees(nTaxa, nTrees, 42L);
            CCD0 ccd0 = new CCD0(trees, 0);
            Clade root = null;
            for (Clade c : ccd0.getClades()) {
                if (c.size() == nTaxa) {
                    root = c;
                }
            }
            if (root == null) {
                continue;
            }

            int a1 = 0;
            int a2 = 0;
            for (CladePartition p : root.getPartitions()) {
                if (p.getNumberOfOccurrences() > 0) {
                    a1++;
                } else {
                    a2++;
                }
            }

            // |A3|: observed proper subclades whose complement within root is NOT an observed clade
            BitSet rootBits = root.getCladeInBits();
            int a3 = 0;
            for (Clade d : ccd0.getClades()) {
                if (d.size() >= root.size()) {
                    continue;
                }
                BitSet db = d.getCladeInBits();
                BitSet tmp = (BitSet) db.clone();
                tmp.and(rootBits);
                if (!tmp.equals(db)) {
                    continue; // not a subclade of root
                }
                BitSet comp = (BitSet) rootBits.clone();
                comp.andNot(db);
                if (ccd0.getClade(comp) == null) {
                    a3++;
                }
            }

            double total = Math.pow(2, nTaxa - 1) - 1;
            double a4 = total - a1 - a2 - a3;

            // per-split constant pseudocount alpha on every class
            double alpha = 0.4;
            double fC = nTrees;
            double denomSplit = fC + alpha * (a1 + a2 + a3 + a4);
            double pA1Split = (fC + alpha * a1) / denomSplit;

            // per-class total pseudocount alpha (spread within each class)
            double denomClass = fC + alpha * 4;
            double pA1Class = (fC + alpha) / denomClass;

            System.out.printf("%-6d %-7d %-8d %-8d %-8d %-14.4g %-14.6g %-14.6g%n",
                    nTaxa, nTrees, a1, a2, a3, a4, pA1Split, pA1Class);
        }
    }
}
