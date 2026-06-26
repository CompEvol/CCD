package test.ccd.model;

import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import ccd.model.GRegCCD;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.TreeMap;

import static org.junit.jupiter.api.Assertions.assertEquals;

/** Validates the GRegCCD model class: exact normalisation and P(T) = eps^(s(T))/Z, on enumerable sets. */
public class GRegCCDTest {

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

    private void check(String label, List<String> taxa, List<T> sample, double eps) {
        List<Tree> sampleTrees = new ArrayList<>();
        for (T s : sample) sampleTrees.add(new TreeParser(taxa, topo(s) + ";", 1, false));
        GRegCCD m = new GRegCCD(sampleTrees, eps);
        double Z = Math.exp(m.getLogPartitionFunction());

        TreeMap<Integer, double[]> byNovel = new TreeMap<>(); // s -> {count, sumP}
        double sum = 0;
        for (T t : allTopologies(taxa)) {
            Tree tree = new TreeParser(taxa, topo(t) + ";", 1, false);
            double p = m.getProbabilityOfTree(tree);
            int s = m.novelSplitCount(tree);
            // P(T) = (product of observed split counts) * eps^s / Z, so P*Z/eps^s is a positive integer.
            double ratio = p * Z / Math.pow(eps, s);
            assertEquals(Math.rint(ratio), ratio, 1e-7 * Math.max(ratio, 1.0),
                    "P(T) not of the form (int) * eps^s / Z for " + topo(t));
            sum += p;
            byNovel.computeIfAbsent(s, k -> new double[2]);
            byNovel.get(s)[0]++; byNovel.get(s)[1] += p;
        }
        System.out.printf("%n%s: %d taxa, eps=%.2f, Z=%.6f -> tiers (new splits: count):%n", label, taxa.size(), eps, Z);
        for (var e : byNovel.entrySet())
            System.out.printf("  %d: %.0f%n", e.getKey(), e.getValue()[0]);
        System.out.printf("  SUM = %.12f%n", sum);
        assertEquals(1.0, sum, 1e-10, "GRegCCD must be exactly normalised");
    }

    @Test
    public void exactlyNormalised() {
        check("4 taxa", Arrays.asList("A", "B", "C", "D"),
                List.of(cat("A", "B", "C", "D"), cat("D", "C", "B", "A")), 0.1);
        check("5 taxa +E", Arrays.asList("A", "B", "C", "D", "E"),
                List.of(cat("A", "B", "C", "D", "E"), cat("D", "C", "B", "A", "E")), 0.1);
        check("6 taxa +(E,F)", Arrays.asList("A", "B", "C", "D", "E", "F"),
                List.of(new Node(cat("A", "B", "C", "D"), new Node(new Leaf("E"), new Leaf("F"))),
                        new Node(cat("D", "C", "B", "A"), new Node(new Leaf("E"), new Leaf("F")))), 0.1);
        // a different eps to confirm normalisation holds generally
        check("5 taxa +E (eps=0.03)", Arrays.asList("A", "B", "C", "D", "E"),
                List.of(cat("A", "B", "C", "D", "E"), cat("D", "C", "B", "A", "E")), 0.03);
    }
}
