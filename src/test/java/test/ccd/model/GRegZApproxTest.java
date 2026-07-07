package test.ccd.model;

import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import ccd.model.GRegCCD;
import ccd.model.GRegZApprox;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/** Compares the large-n GRegZApprox partition function to GRegCCD's exact Z on small sets. */
public class GRegZApproxTest {

    private sealed interface T permits Leaf, Node {}
    private record Leaf(String name) implements T {}
    private record Node(T l, T r) implements T {}

    private static String topo(T t) {
        if (t instanceof Leaf lf) return lf.name();
        Node n = (Node) t;
        return "(" + topo(n.l()) + "," + topo(n.r()) + ")";
    }

    private static T cat(String... names) {
        T t = new Node(new Leaf(names[0]), new Leaf(names[1]));
        for (int i = 2; i < names.length; i++) t = new Node(t, new Leaf(names[i]));
        return t;
    }

    private void compare(String label, List<String> taxa, List<T> sample, double eps) {
        List<Tree> trees = new ArrayList<>();
        for (T s : sample) trees.add(new TreeParser(taxa, topo(s) + ";", 1, false));
        double exact = new GRegCCD(trees, eps).getLogPartitionFunction();
        double approx = GRegZApprox.fromTrees(trees).logZ(eps);
        double relErr = Math.abs(Math.exp(approx - exact) - 1.0);
        System.out.printf("%-22s eps=%.3f  logZ exact=%.6f  approx=%.6f  relErr=%.3e%n",
                label, eps, exact, approx, relErr);
    }

    @Test
    public void approxMatchesExact() {
        var t4 = Arrays.asList("A", "B", "C", "D");
        var s4 = List.<T>of(cat("A", "B", "C", "D"), cat("D", "C", "B", "A"));
        compare("4 taxa", t4, s4, 0.1);
        compare("4 taxa", t4, s4, 0.01);

        var t5 = Arrays.asList("A", "B", "C", "D", "E");
        var s5 = List.<T>of(cat("A", "B", "C", "D", "E"), cat("D", "C", "B", "A", "E"));
        compare("5 taxa +E", t5, s5, 0.1);
        compare("5 taxa +E", t5, s5, 0.01);

        var t6 = Arrays.asList("A", "B", "C", "D", "E", "F");
        var s6 = List.<T>of(
                new Node(cat("A", "B", "C", "D"), new Node(new Leaf("E"), new Leaf("F"))),
                new Node(cat("D", "C", "B", "A"), new Node(new Leaf("E"), new Leaf("F"))));
        compare("6 taxa +(E,F)", t6, s6, 0.1);
        compare("6 taxa +(E,F)", t6, s6, 0.01);
    }
}
