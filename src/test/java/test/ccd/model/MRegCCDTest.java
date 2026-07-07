package test.ccd.model;

import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import ccd.model.CCD0;
import ccd.model.Clade;
import ccd.model.CladePartition;
import ccd.model.MRegCCD;
import ccd.model.bitsets.BitSet;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Validates the one-parameter per-new-split {@link MRegCCD}: that with the full reserve depth it is an
 * exactly normalised distribution on enumerable taxon sets, and that truncating the reserve without a
 * tail correction super-normalises it (the artefact that made order-2 falsely appear to close the gap
 * to KRegCCD in the RSV2 experiment).
 */
public class MRegCCDTest {

    /* -- labelled rooted-binary-topology enumeration (as in GRegCCDTest) -- */
    private sealed interface T permits Leaf, Node {}
    private record Leaf(String name) implements T {}
    private record Node(T l, T r) implements T {}

    private static String topo(T t) {
        if (t instanceof Leaf lf) return lf.name();
        Node n = (Node) t;
        return "(" + topo(n.l()) + "," + topo(n.r()) + ")";
    }

    private static List<T> insertAll(T t, String x) {
        List<T> out = new ArrayList<>();
        out.add(new Node(new Leaf(x), t));
        if (t instanceof Node n) {
            for (T l2 : insertAll(n.l(), x)) out.add(new Node(l2, n.r()));
            for (T r2 : insertAll(n.r(), x)) out.add(new Node(n.l(), r2));
        }
        return out;
    }

    private static List<T> allTopologies(List<String> taxa) {
        List<T> trees = new ArrayList<>();
        trees.add(new Leaf(taxa.get(0)));
        for (int i = 1; i < taxa.size(); i++) {
            List<T> next = new ArrayList<>();
            for (T t : trees) next.addAll(insertAll(t, taxa.get(i)));
            trees = next;
        }
        return trees;
    }

    private static T cat(String... names) {
        T t = new Node(new Leaf(names[0]), new Leaf(names[1]));
        for (int i = 2; i < names.length; i++) t = new Node(t, new Leaf(names[i]));
        return t;
    }

    private static List<Tree> trees(List<String> taxa, List<T> shapes) {
        List<Tree> out = new ArrayList<>();
        for (T s : shapes) out.add(new TreeParser(taxa, topo(s) + ";", 1, false));
        return out;
    }

    private static double totalMass(MRegCCD m, List<String> taxa) {
        double sum = 0.0;
        for (T t : allTopologies(taxa)) {
            Tree tree = new TreeParser(taxa, topo(t) + ";", 1, false);
            sum += Math.exp(m.getLogProbabilityOfTree(tree));
        }
        return sum;
    }

    @Test
    public void fullDepthIsExactlyNormalised() {
        for (double mu : new double[]{0.02, 0.1, 0.25}) {
            check5(mu);
            check6(mu);
        }
    }

    private void check5(double mu) {
        List<String> taxa = Arrays.asList("A", "B", "C", "D", "E");
        List<Tree> train = trees(taxa,
                List.of(cat("A", "B", "C", "D", "E"), cat("D", "C", "B", "A", "E")));
        MRegCCD m = new MRegCCD(train, 0.0, mu, taxa.size(), false); // full depth -> no omitted tail
        double sum = totalMass(m, taxa);
        System.out.printf("MRegCCD 5 taxa mu=%.2f full-depth SUM = %.12f%n", mu, sum);
        assertEquals(1.0, sum, 1e-9, "MRegCCD at full reserve depth must be exactly normalised");
    }

    private void check6(double mu) {
        List<String> taxa = Arrays.asList("A", "B", "C", "D", "E", "F");
        List<Tree> train = trees(taxa,
                List.of(new Node(cat("A", "B", "C", "D"), new Node(new Leaf("E"), new Leaf("F"))),
                        new Node(cat("D", "C", "B", "A"), new Node(new Leaf("E"), new Leaf("F")))));
        MRegCCD m = new MRegCCD(train, 0.0, mu, taxa.size(), false);
        double sum = totalMass(m, taxa);
        System.out.printf("MRegCCD 6 taxa mu=%.2f full-depth SUM = %.12f%n", mu, sum);
        assertEquals(1.0, sum, 1e-9, "MRegCCD at full reserve depth must be exactly normalised");
    }

    @Test
    public void m2EqualsCCD0ExpandedSplits() {
        // A diverse but incomplete training set so many clades are observed without all their splits,
        // creating recombinations (CCD0-expanded splits = M_2).
        List<String> taxa = Arrays.asList("A", "B", "C", "D", "E", "F");
        List<T> all = allTopologies(taxa);
        List<T> picks = new ArrayList<>();
        for (int i = 0; i < all.size(); i += 47) picks.add(all.get(i)); // ~20 trees spread across the space
        List<Tree> train = trees(taxa, picks);

        MRegCCD mreg = new MRegCCD(train, 0.0, 0.05);
        CCD0 ccd0 = new CCD0(train, 0);

        int checked = 0, withRecomb = 0;
        for (Clade c : mreg.getClades()) {
            if (c.size() < 2) continue;
            BitSet bits = c.getCladeInBits();
            int m2 = mreg.reserveCounts(bits)[2];

            Clade c0 = ccd0.getClade(bits);
            int expanded = 0;
            for (CladePartition p : c0.getPartitions()) {
                if (p.getNumberOfOccurrences() == 0) expanded++; // CCD0 split that was never observed
            }
            assertEquals(expanded, m2,
                    "M_2 must equal the CCD0-expanded split count for clade " + bits);
            checked++;
            if (m2 > 0) withRecomb++;
        }
        System.out.printf("M_2 == CCD0-expanded splits verified on %d clades (%d with recombinations)%n",
                checked, withRecomb);
        assertTrue(withRecomb > 0, "test should exercise clades that actually have recombinations");
    }

    @Test
    public void samplerMatchesScorer() {
        // E_q[-log q(S)] = H(q): if the simulator draws from q AND reports the correct log q, the mean
        // sampled -logp equals the entropy computed by enumeration with the scorer. Use the tail-OFF
        // model so the self-consistent sampler and the scorer are the same distribution.
        List<String> taxa = Arrays.asList("A", "B", "C", "D", "E", "F");
        List<T> all = allTopologies(taxa);
        List<T> picks = new ArrayList<>();
        for (int i = 0; i < all.size(); i += 31) picks.add(all.get(i));
        List<Tree> train = trees(taxa, picks);
        MRegCCD m = new MRegCCD(train, 0.0, 0.1, taxa.size(), false); // full depth, tail off

        // true entropy and normalisation by enumeration (scorer)
        double sum = 0.0, H = 0.0;
        for (T t : all) {
            Tree tree = new TreeParser(taxa, topo(t) + ";", 1, false);
            double logq = m.getLogProbabilityOfTree(tree);
            double q = Math.exp(logq);
            sum += q;
            H -= q * logq;
        }
        assertEquals(1.0, sum, 1e-9, "tail-off full-depth model must be normalised");

        // Monte-Carlo entropy from the sampler
        int N = 300_000;
        double s1 = 0.0, s2 = 0.0;
        for (int i = 0; i < N; i++) {
            double logp = m.sampleTreeLogProbability();
            s1 += -logp;
            s2 += logp * logp;
        }
        double hHat = s1 / N;
        double se = Math.sqrt(Math.max(0, s2 / N - hHat * hHat) / N);
        System.out.printf("MRegCCD sampler: H_enum=%.5f  H_MC=%.5f +/- %.5f (%.1f SE off)%n",
                H, hHat, se, Math.abs(hHat - H) / se);
        assertEquals(H, hHat, Math.max(5 * se, 0.01),
                "sampler entropy must match the scorer's enumerated entropy");
    }

    @Test
    public void truncatedReserveSuperNormalises() {
        List<String> taxa = Arrays.asList("A", "B", "C", "D", "E", "F");
        List<Tree> train = trees(taxa,
                List.of(new Node(cat("A", "B", "C", "D"), new Node(new Leaf("E"), new Leaf("F"))),
                        new Node(cat("D", "C", "B", "A"), new Node(new Leaf("E"), new Leaf("F")))));
        double mu = 0.2;
        double full = totalMass(new MRegCCD(train, 0.0, mu, taxa.size(), false), taxa);
        double order2 = totalMass(new MRegCCD(train, 0.0, mu, 2, false), taxa); // M2 only, no tail
        System.out.printf("MRegCCD 6 taxa mu=%.2f: full-depth SUM=%.9f  order-2 SUM=%.9f%n", mu, full, order2);
        assertEquals(1.0, full, 1e-9, "full depth normalised");
        assertTrue(order2 > 1.0 + 1e-4,
                "order-2 reserve (no tail) must super-normalise (sum > 1), got " + order2);
    }
}
