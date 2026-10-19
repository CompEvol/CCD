package ccd.model;

import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator;
import ccd.tools.CCDToolUtil;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.Test;
import java.io.File;
import java.util.*;

/** Checks the newick round-trip used by the comparison harness against direct tree objects. */
public class MunroParseCheck {
    private static final String PATH = System.getProperty("ccd.trees", "");

    @Test
    public void roundTripPreservesTopologies() throws Exception {
        Assumptions.assumeTrue(!PATH.isEmpty() && new File(PATH).exists());
        TreeAnnotator.TreeSet ts = CCDToolUtil.getTreeSet(PATH, 10);
        ts.reset();
        List<Tree> direct = new ArrayList<>();
        while (ts.hasNext() && direct.size() < 400) direct.add(ts.next());
        Tree t0 = direct.get(0);
        List<String> taxa = new ArrayList<>(List.of(t0.getTaxaNames()));
        while (taxa.remove(null)) { }
        System.out.println("taxa[0..3] = " + taxa.subList(0, Math.min(4, taxa.size())));
        String nwk = t0.getRoot().toShortNewick(false);
        System.out.println("round-tripped newick head: " + nwk.substring(0, Math.min(150, nwk.length())));

        // compare clade sets: direct vs reparsed
        int mismatches = 0;
        for (int i = 0; i < direct.size(); i++) {
            String s = direct.get(i).getRoot().toShortNewick(false) + ";";
            Tree re = new beast.base.evolution.tree.TreeParser(taxa, s, 0, false);
            Set<String> a = clades(direct.get(i)), b = clades(re);
            if (!a.equals(b)) mismatches++;
        }
        System.out.printf("clade-set mismatches after round trip: %d of %d trees%n",
                mismatches, direct.size());
    }

    private static Set<String> clades(Tree t) {
        Set<String> out = new TreeSet<>();
        collect(t.getRoot(), out);
        return out;
    }

    private static SortedSet<String> collect(beast.base.evolution.tree.Node v, Set<String> out) {
        SortedSet<String> s = new TreeSet<>();
        if (v.isLeaf()) { s.add(v.getID() != null ? v.getID() : String.valueOf(v.getNr())); }
        else for (beast.base.evolution.tree.Node c : v.getChildren()) s.addAll(collect(c, out));
        out.add(String.join(",", s));
        return s;
    }
}
