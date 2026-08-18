package ccd.model;

import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator;
import ccd.tools.CCDToolUtil;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/** How much of MRegCCD's scoring cost is its reserve depth, versus the boundary enumeration itself? */
public class MRegDepthTimingTest {

    private static final String PATH = System.getProperty("ccd.trees", "");

    private static List<Tree> read(int count, int skip) throws Exception {
        TreeAnnotator.TreeSet ts = CCDToolUtil.getTreeSet(PATH, 10);
        ts.reset();
        List<Tree> all = new ArrayList<>();
        while (ts.hasNext()) {
            all.add(ts.next());
        }
        List<Tree> pool = all.subList(skip * all.size() / 2, (skip + 1) * all.size() / 2);
        List<Tree> out = new ArrayList<>();
        double step = Math.max(1.0, pool.size() / (double) count);
        for (int i = 0; i < count && (int) (i * step) < pool.size(); i++) {
            out.add(pool.get((int) (i * step)));
        }
        return out;
    }

    @Test
    public void depthVersusCost() throws Exception {
        Assumptions.assumeTrue(!PATH.isEmpty() && new File(PATH).exists());
        List<Tree> test = read(200, 1);
        System.out.printf("%n=== %s: MRegCCD cost by reserve depth ===%n", new File(PATH).getName());
        System.out.printf("%-7s %-9s %-12s %-12s %-14s%n", "depth", "impl", "construct", "score/200", "mean logP");
        for (int depth : new int[]{2, 3, 4}) {
            for (String which : new String[]{"MRegCCDSlow", "MRegCCD"}) {
                long t0 = System.nanoTime();
                MRegCCDSlow m = which.equals("MRegCCDSlow")
                        ? new MRegCCDSlow(read(500, 0), 0.0, MRegCCDSlow.DEFAULT_MU, depth, true)
                        : new MRegCCD(read(500, 0), 0.0, MRegCCDSlow.DEFAULT_MU, depth, true);
                long build = (System.nanoTime() - t0) / 1_000_000L;
                t0 = System.nanoTime();
                double sum = 0;
                for (Tree t : test) {
                    sum += m.getLogProbabilityOfTree(t);
                }
                long score = (System.nanoTime() - t0) / 1_000_000L;
                System.out.printf("%-7d %-9s %9dms %9dms %14.6f%n",
                        depth, which, build, score, sum / test.size());
            }
        }
    }
}
