package ccd.tools;

import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator;
import ccd.model.CCD1;
import ccd.model.CRegCCD;
import ccd.model.KRegCCD;
import ccd.model.MRegCCD;
import ccd.model.RegCCD;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.Set;
import java.util.List;
import java.util.Map;

/**
 * Held-out predictive comparison of the CCD variants across many prepared posteriors, in one JVM.
 *
 * <pre>
 *   java -cp ... ccd.tools.HeadToHeadBenchmark &lt;prepared-dir&gt; &lt;out.csv&gt; &lt;pertree-dir&gt; [models]
 * </pre>
 *
 * <p>Each dataset must already have been prepared into {@code base.fit/val/train/test.trees}: the
 * model is fitted on {@code train}, its hyperparameters chosen on the disjoint {@code fit}/{@code
 * val} pair, and scored on {@code test}. Preparing separately means each source chain is parsed
 * once rather than once per model, and running every dataset in one process means the work is
 * modelling rather than process startup.
 *
 * <p>Per-tree log probabilities are written per dataset and model, so paired statistics can be
 * computed afterwards between any two models without re-running anything.
 */
public class HeadToHeadBenchmark {

    private static final List<String> ALL =
            List.of("CCD1", "RegCCD", "KRegCCD", "MRegCCD", "CRegCCD");

    public static void main(String[] args) throws Exception {
        File prepared = new File(args[0]);
        File out = new File(args[1]);
        File perTree = new File(args[2]);
        List<String> models = args.length > 3 ? Arrays.asList(args[3].split(",")) : ALL;
        perTree.mkdirs();

        List<String> bases = new ArrayList<>();
        for (File f : prepared.listFiles((d, n) -> n.endsWith(".train.trees"))) {
            bases.add(f.getName().replaceAll("\\.train\\.trees$", ""));
        }
        bases.sort(String::compareTo);
        System.out.printf("%d datasets, models %s%n", bases.size(), models);

        // Resume support. A sweep over this corpus takes about an hour, and a JVM that dies partway
        // (heap exhaustion, SIGBUS from an exhausted swap, an external kill) used to cost every
        // dataset already scored. So read back whatever rows the CSV already holds and skip those
        // datasets, appending rather than truncating. Delete the CSV to force a clean run.
        Set<String> done = new LinkedHashSet<>();
        boolean header = false;
        if (out.exists() && out.length() > 0) {
            List<String> lines = java.nio.file.Files.readAllLines(out.toPath());
            for (int i = 0; i < lines.size(); i++) {
                String line = lines.get(i);
                if (line.isBlank()) continue;
                if (i == 0 && line.startsWith("dataset,")) { header = true; continue; }
                done.add(line.substring(0, line.indexOf(',')));
            }
            System.out.printf("resuming: %d datasets already in %s, skipping them%n",
                    done.size(), out.getName());
        }
        PrintWriter csv = new PrintWriter(new FileWriter(out, true));
        PrintWriter failed = new PrintWriter(new FileWriter(new File(out.getPath() + ".failed"), true));
        long t0 = System.nanoTime();

        for (String base : bases) {
            if (done.contains(base)) continue;
            Map<String, String> row = new LinkedHashMap<>();
            try {
                List<Tree> train = read(prepared, base, "train");
                List<Tree> test = read(prepared, base, "test");
                row.put("dataset", base);
                row.put("taxa", String.valueOf(train.get(0).getLeafNodeCount()));
                row.put("nTrain", String.valueOf(train.size()));
                row.put("nTest", String.valueOf(test.size()));
                row.put("testProtocol",
                        new File(prepared, base + ".siblingtest").exists() ? "sibling-run" : "chain-half");
                row.put("CCD1_entropy", fmt(new CCD1(read(prepared, base, "train"), 0.0).getEntropy()));

                for (String m : models) {
                    runModel(m, prepared, base, test, row, perTree);
                }
                if (!header) {
                    csv.println(String.join(",", row.keySet()));
                    header = true;
                }
                csv.println(String.join(",", row.values()));
                csv.flush();
                System.out.printf("  %-70s %s taxa%n", base.substring(0, Math.min(70, base.length())),
                        row.get("taxa"));
            } catch (Throwable e) {
                failed.printf("%s,%s: %s%n", base, e.getClass().getSimpleName(), e.getMessage());
                failed.flush();
                System.out.printf("  FAILED %s (%s)%n", base, e.getClass().getSimpleName());
            }
        }
        csv.close();
        failed.close();
        System.out.printf("done in %.1f s%n", (System.nanoTime() - t0) / 1e9);
    }

    private static void runModel(String model, File dir, String base, List<Tree> test,
                                 Map<String, String> row, File perTree) {
        try {
            List<Tree> val = read(dir, base, "val");
            String params;
            long build;
            Scorer scorer;

            switch (model) {
                case "CCD1" -> {
                    long t = System.nanoTime();
                    CCD1 m = new CCD1(read(dir, base, "train"), 0.0);
                    build = ms(t);
                    params = "-";
                    scorer = m::getLogProbabilityOfTree;
                }
                case "RegCCD" -> {
                    double best = pick(new double[]{0.01, 0.05, 0.1, 0.2, 0.4, 0.8, 1.0}, a -> {
                        RegCCD f = new RegCCD(read(dir, base, "fit"), 0.0, a);
                        return mean(f::getLogProbabilityOfTree, val);
                    });
                    long t = System.nanoTime();
                    RegCCD m = new RegCCD(read(dir, base, "train"), 0.0, best);
                    build = ms(t);
                    params = "alpha=" + fmt(best);
                    scorer = m::getLogProbabilityOfTree;
                }
                case "KRegCCD" -> {
                    KRegCCD f = new KRegCCD(read(dir, base, "fit"), 0.0, 0.005, 0.4);
                    double best = pick(new double[]{0.00002, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05},
                            mu -> mean(t -> f.getLogProbabilityOfTree(t, mu), val));
                    long t = System.nanoTime();
                    KRegCCD m = new KRegCCD(read(dir, base, "train"), 0.0, best, 0.4);
                    m.precomputeReserves();
                    build = ms(t);
                    params = "alpha=0.4;mu=" + fmt(best);
                    scorer = m::getLogProbabilityOfTree;
                }
                case "MRegCCD" -> {
                    MRegCCD f = new MRegCCD(read(dir, base, "fit"), 0.0, MRegCCD.DEFAULT_MU);
                    double best = pick(new double[]{0.0002, 0.001, 0.002, 0.008, 0.0159, 0.05, 0.1},
                            mu -> mean(t -> f.getLogProbabilityOfTree(t, mu), val));
                    long t = System.nanoTime();
                    MRegCCD m = new MRegCCD(read(dir, base, "train"), 0.0, best);
                    build = ms(t);
                    params = "mu=" + fmt(best);
                    scorer = m::getLogProbabilityOfTree;
                }
                case "CRegCCD" -> {
                    CRegCCD f = new CRegCCD(read(dir, base, "fit"), 0.0);
                    // alpha is a per-split pseudocount, so alpha > 1 lets the prior outweigh a real
                    // observation on every split at once: at alpha = 12 a split seen once is only
                    // 1.08x as probable as one never seen, and the counts stop carrying information.
                    // The one reason to allow alpha > 1 is that MCMC samples are autocorrelated, so
                    // f(S) overstates the evidence by roughly N/ESS; capping at 5 admits that
                    // correction down to ESS = 200 out of the 1000 training trees and no further.
                    double[] gridA = {0.002, 0.01, 0.05, 0.2, 0.4, 1.0, 2.0, 5.0};
                    // alpha1 and alpha2 are class totals rather than per-split pseudocounts, so the
                    // argument above does not bound them; the bands are set by what is ever chosen.
                    // No dataset selected alpha1 above 2.0 or alpha2 above 1.0, so the 12.0 band is
                    // dead weight for both; each keeps one band of headroom above its observed max
                    // so that the top selected value is interior rather than pinned to a ceiling.
                    double[] gridA1 = {0.002, 0.01, 0.05, 0.2, 0.4, 1.0, 2.0, 5.0};
                    // alpha2 was selected at the old floor of 0.002 on 43% of datasets, which is a
                    // boundary rather than an optimum, so the floor drops four bands to locate it.
                    double[] gridA2 = {0.000016, 0.00008, 0.0004, 0.002, 0.01, 0.05, 0.2, 0.4, 1.0, 2.0};
                    double bA = 0.4, b1 = 0.4, b2 = 0.05, bestScore = Double.NEGATIVE_INFINITY;
                    for (double a : gridA) {
                        for (double a1 : gridA1) {
                            for (double a2 : gridA2) {
                                double sc = mean(t -> f.getLogProbabilityOfTree(t, a, a1, a2), val);
                                if (sc > bestScore) {
                                    bestScore = sc; bA = a; b1 = a1; b2 = a2;
                                }
                            }
                        }
                    }
                    long t = System.nanoTime();
                    CRegCCD m = new CRegCCD(read(dir, base, "train"), 0.0, bA, b1, b2);
                    build = ms(t);
                    params = "alpha=" + fmt(bA) + ";alpha1=" + fmt(b1) + ";alpha2=" + fmt(b2);
                    scorer = m::getLogProbabilityOfTree;
                }
                default -> throw new IllegalArgumentException("unknown model " + model);
            }

            long t = System.nanoTime();
            int covered = 0;
            double sum = 0;
            try (PrintWriter w = new PrintWriter(new File(perTree, base + "__" + model + ".txt"))) {
                for (Tree x : test) {
                    double lp = scorer.logP(x);
                    w.println(lp);
                    if (Double.isFinite(lp)) { covered++; sum += lp; }
                }
            }
            long score = ms(t);
            row.put(model + "_params", params);
            row.put(model + "_coverage", fmt(covered / (double) test.size()));
            row.put(model + "_meanLogP", covered == test.size() ? fmt(sum / test.size()) : "-inf");
            row.put(model + "_constructMs", String.valueOf(build));
            row.put(model + "_scoreMs", String.valueOf(score));
        } catch (Throwable e) {
            row.put(model + "_params", "FAILED:" + e.getClass().getSimpleName());
            row.put(model + "_coverage", "");
            row.put(model + "_meanLogP", "");
            row.put(model + "_constructMs", "");
            row.put(model + "_scoreMs", "");
        }
    }

    interface Scorer { double logP(Tree t); }
    interface Obj { double at(double x); }

    private static double pick(double[] grid, Obj f) {
        double best = grid[0], bestScore = Double.NEGATIVE_INFINITY;
        for (double x : grid) {
            double s = f.at(x);
            if (s > bestScore) { bestScore = s; best = x; }
        }
        return best;
    }

    private static double mean(Scorer s, List<Tree> trees) {
        double sum = 0; int n = 0;
        for (Tree t : trees) {
            double lp = s.logP(t);
            if (Double.isFinite(lp)) { sum += lp; n++; }
        }
        return n == 0 ? Double.NEGATIVE_INFINITY : sum / n;
    }

    private static List<Tree> read(File dir, String base, String part) {
        try {
            TreeAnnotator.TreeSet ts = CCDToolUtil.getTreeSet(
                    new File(dir, base + "." + part + ".trees").getAbsolutePath(), 0);
            ts.reset();
            List<Tree> out = new ArrayList<>();
            while (ts.hasNext()) out.add(ts.next());
            return out;
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    private static long ms(long t0) { return (System.nanoTime() - t0) / 1_000_000L; }
    /**
     * Formats a hyperparameter without losing it. A fixed "%.4f" collapsed every grid value below
     * 5e-5 to "0.0000" -- the alpha2 bands at 1.6e-5, 8e-5 and 4e-4 and the KRegCCD mu band at
     * 2e-5 all became the same unreadable string -- so the selected value could not be recovered
     * from the CSV. Double.toString round-trips exactly and Python's float() parses it.
     */
    private static String fmt(double d) { return Double.toString(d); }
}
