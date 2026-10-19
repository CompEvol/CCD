package ccd.model;

import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator;
import ccd.tools.CCDToolUtil;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Head-to-head held-out predictive comparison of the regularised CCD variants on a real posterior
 * tree sample.
 *
 * <p>Protocol follows the manuscript's RSV2 comparison: the model is trained on trees drawn from the
 * first half of the chain and scored on trees from the second half, so the test trees are genuinely
 * out of training. Each model's hyperparameters are selected on an inner fit/validation split of the
 * training half alone, then the model is rebuilt on the whole training half and scored on the test
 * half. Reported per model: support coverage (fraction of test trees with positive probability), mean
 * log probability over all test trees, and -- for the full-support models -- the paired per-tree
 * comparison against KRegCCD.
 *
 * <p>Point the test at a tree file with {@code -Dccd.trees=/path/to/x.trees}; it is skipped when no
 * file is given. {@code -Dccd.n=1000} sets the number of training and test trees.
 */
public class RealDataHeadToHeadTest {

    private static String PATH = System.getProperty("ccd.trees", "");
    private static final int N = Integer.parseInt(System.getProperty("ccd.n", "1000"));
    private static final double BURNIN_PERCENT = 10;

    /** Which models to run, comma separated; default all. Running one model at a time isolates
     *  cost and stops a model that fails on a dataset from losing the others' results for it. */
    private static final String MODELS = System.getProperty("ccd.models", "");
    /** Directory for per-tree log probabilities, one file per dataset and model. Paired statistics
     *  are computed from these afterwards, so models need not run together. */
    private static final String PERTREE = System.getProperty("ccd.pertree", "");
    /**
     * A second chain to draw the test set from, for analyses that were run more than once. Splitting
     * one chain in half conflates model quality with convergence: the halves differ partly because
     * the sampler was still moving. An independent run removes that, so where replicate runs exist
     * the test set is taken from a sibling run and the whole of the primary chain trains.
     */
    private static final String TEST_TREES = System.getProperty("ccd.testtrees", "");

    /** Directory of prepared subsets; when set, the source chain is never touched. */
    private static final String PREPARED = System.getProperty("ccd.prepared", "");

    private static List<Tree> readPrepared(String part) throws Exception {
        String base = new File(PATH).getName().replaceAll("\\.trees$", "");
        File f = new File(PREPARED, base + "." + part + ".trees");
        TreeAnnotator.TreeSet ts = CCDToolUtil.getTreeSet(f.getAbsolutePath(), 0);
        ts.reset();
        List<Tree> out = new ArrayList<>();
        while (ts.hasNext()) {
            out.add(ts.next());
        }
        if (chainTrees == null) {
            chainTrees = out.size();
            System.out.printf("prepared %s: %d trees, %d taxa%n",
                    base, out.size(), out.get(0).getLeafNodeCount());
        }
        return out;
    }

    private static boolean runs(String model) {
        return MODELS.isEmpty() || ("," + MODELS + ",").contains("," + model + ",");
    }

    private static void dumpPerTree(String model, Scorer s, List<Tree> test) throws Exception {
        if (PERTREE.isEmpty()) {
            return;
        }
        new File(PERTREE).mkdirs();
        String name = new File(PATH).getName().replaceAll("[^A-Za-z0-9._-]", "_");
        try (java.io.PrintWriter w = new java.io.PrintWriter(
                new File(PERTREE, name + "__" + model + ".txt"))) {
            for (Tree t : test) {
                w.println(s.logP(t));
            }
        }
    }

    private interface Scorer {
        double logP(Tree t);
    }

    private static Integer chainTrees;

    /** One cheap pass to size the stride before caching anything. */
    private static int countTrees() throws Exception {
        TreeAnnotator.TreeSet ts = CCDToolUtil.getTreeSet(PATH, (int) BURNIN_PERCENT);
        ts.reset();
        int n = 0;
        while (ts.hasNext()) {
            ts.next();
            n++;
        }
        return n;
    }

    private static List<String> newickCache;
    private static List<String> taxaCache;

    /**
     * Parses the tree file once, keeping each topology as a newick string.
     *
     * <p>The topologies are written with node numbers rather than taxon labels, and the taxon list
     * supplies the mapping. Labels cannot be round-tripped safely: bracketed genus names are standard
     * for unvalidated nomenclature, so real files contain labels like
     * {@code GluRS-B_AF_Bact_[Eubacterium]_eligens_...}, and a newick parser reads {@code [} as a
     * comment opener. Stripping the brackets corrupts the name; keeping them fails to parse.
     * MunroParseCheck verifies that the round trip preserves clade sets exactly.
     */
    /** Upper bound on cached topologies; the experiment never needs more than a few thousand. */
    private static final int MAX_CACHED = 8000;

    private static void load() throws Exception {
        if (newickCache != null) {
            return;
        }
        int total = countTrees();
        int stride = Math.max(1, (total + MAX_CACHED - 1) / MAX_CACHED);
        TreeAnnotator.TreeSet ts = CCDToolUtil.getTreeSet(PATH, (int) BURNIN_PERCENT);
        ts.reset();
        List<String> nwk = new ArrayList<>();
        List<String> taxa = null;
        // Keep at most MAX_CACHED topologies, evenly spaced through the chain. The experiment draws
        // evenly spaced subsets anyway, so this changes nothing statistically, but it bounds memory:
        // holding every tree of a 36k-tree chain on several hundred taxa exhausts the heap.
        int chain = 0;
        while (ts.hasNext()) {
            Tree t = ts.next();
            chain++;
            if (taxa == null) {
                taxa = new ArrayList<>(List.of(t.getTaxaNames()));
                while (taxa.remove(null)) {
                    // getTaxaNames() is sized for all nodes on some inputs; drop the internal slots
                }
            }
            if (stride > 1 && (chain - 1) % stride != 0) {
                continue;
            }
            nwk.add(t.getRoot().toShortNewick(false) + ";");
        }
        newickCache = nwk;
        taxaCache = taxa;
        chainTrees = chain;
        System.out.printf("loaded %d of %d trees (stride %d), %d taxa from %s%n",
                nwk.size(), chain, stride, taxa.size(), new File(PATH).getName());
    }

    /**
     * {@code count} evenly spaced trees from the chain segment {@code [from, to)} (as fractions of
     * the post-burn-in chain), freshly parsed on every call because the CCD constructors take
     * ownership of the trees they are given.
     *
     * <p>Segments must be disjoint for hyperparameter selection to be honest: scoring trees that are
     * also in the fitted set drives every regularisation parameter to zero, because the backbone
     * already fits its own training trees and any reserved mass is then pure loss.
     */
    private static List<Tree> read(int count, double from, double to) throws Exception {
        load();
        List<String> pool = newickCache.subList((int) (from * newickCache.size()),
                (int) (to * newickCache.size()));
        List<Tree> out = new ArrayList<>();
        double step = Math.max(1.0, pool.size() / (double) count);
        for (int i = 0; i < count && (int) (i * step) < pool.size(); i++) {
            out.add(new beast.base.evolution.tree.TreeParser(taxaCache, pool.get((int) (i * step)),
                    0, false));
        }
        return out;
    }

    /** Trees for the final models and for scoring: first half of the chain trains, second half tests. */
    private static List<Tree> train(int count) throws Exception {
        return PREPARED.isEmpty() ? read(count, 0.0, 0.5) : readPrepared("train");
    }

    /** Disjoint inner split of the training half: [0, 0.25) fits, [0.25, 0.5) validates. */
    private static List<Tree> fitSet(int count) throws Exception {
        return PREPARED.isEmpty() ? read(count, 0.0, 0.25) : readPrepared("fit");
    }

    private static List<Tree> valSet(int count) throws Exception {
        return PREPARED.isEmpty() ? read(count, 0.25, 0.5) : readPrepared("val");
    }

    /** Machine-readable output: one row per dataset, appended to -Dccd.csv if set. */
    private static final String CSV = System.getProperty("ccd.csv", "");
    private static final Map<String, String> CSV_CELLS = new java.util.LinkedHashMap<>();

    /** Milliseconds since a System.nanoTime() mark. */
    private static long ms(long t0) {
        return (System.nanoTime() - t0) / 1_000_000L;
    }

    private static void cell(String key, Object value) {
        CSV_CELLS.put(key, String.valueOf(value));
    }

    private static void writeCsv() throws Exception {
        if (CSV.isEmpty()) {
            return;
        }
        File f = new File(CSV);
        boolean header = !f.exists() || f.length() == 0;
        try (java.io.PrintWriter w = new java.io.PrintWriter(new java.io.FileWriter(f, true))) {
            if (header) {
                w.println(String.join(",", CSV_CELLS.keySet()));
            }
            w.println(String.join(",", CSV_CELLS.values()));
        }
    }

    private static double[] score(Scorer s, List<Tree> test) {
        int covered = 0;
        double sum = 0.0;
        for (Tree t : test) {
            double lp = s.logP(t);
            if (Double.isFinite(lp)) {
                covered++;
                sum += lp;
            }
        }
        return new double[]{covered, covered == 0 ? Double.NEGATIVE_INFINITY : sum / covered};
    }

    /** Mean log probability over every test tree, treating unsupported trees as -infinity. */
    private static double meanAll(Scorer s, List<Tree> test) {
        double sum = 0.0;
        for (Tree t : test) {
            double lp = s.logP(t);
            if (!Double.isFinite(lp)) {
                return Double.NEGATIVE_INFINITY;
            }
            sum += lp;
        }
        return sum / test.size();
    }

    /**
     * Writes the four tree subsets this experiment uses -- fit, validation, train, test -- as small
     * NEXUS files, so the source chain is parsed once rather than once per model.
     *
     * <p>Trees are written with a Translate block and numeric newick: taxon labels appear only in the
     * Translate block, which keeps labels containing brackets (standard for unvalidated genus names)
     * out of the newick, where a parser would read them as comments.
     */
    @Test
    public void prepareSubsets() throws Exception {
        String outDir = System.getProperty("ccd.prepare", "");
        Assumptions.assumeTrue(!outDir.isEmpty() && !PATH.isEmpty() && new File(PATH).exists());
        String base = new File(PATH).getName().replaceAll("\\.trees$", "");
        File dir = new File(outDir);
        dir.mkdirs();
        load();
        // With a sibling run the whole primary chain trains and the sibling supplies the test set;
        // otherwise the chain is split in half, as the manuscript's RSV2 comparison does.
        String[][] parts = TEST_TREES.isEmpty()
                ? new String[][]{{"fit", "0.0", "0.25"}, {"val", "0.25", "0.5"},
                                 {"train", "0.0", "0.5"}, {"test", "0.5", "1.0"}}
                : new String[][]{{"fit", "0.0", "0.5"}, {"val", "0.5", "1.0"},
                                 {"train", "0.0", "1.0"}};
        for (String[] part : parts) {
            int count = part[0].equals("fit") || part[0].equals("val") ? N / 2 : N;
            List<Tree> sub = read(count, Double.parseDouble(part[1]), Double.parseDouble(part[2]));
            File out = writeSubset(dir, base, part[0], sub);
            System.out.printf("wrote %s (%d trees)%n", out.getName(), sub.size());
        }
        if (!TEST_TREES.isEmpty()) {
            writeSubset(dir, base, "test", readFrom(TEST_TREES, N));
            // marker so the run phase can record which protocol produced the test set
            try (java.io.PrintWriter w = new java.io.PrintWriter(
                    new File(dir, base + ".siblingtest"))) {
                w.println(new File(TEST_TREES).getName());
            }
            System.out.printf("test set taken from sibling run %s%n", new File(TEST_TREES).getName());
        }
    }

    /** Evenly spaced trees from the whole post-burn-in chain of another file. */
    private static List<Tree> readFrom(String path, int count) throws Exception {
        TreeAnnotator.TreeSet ts = CCDToolUtil.getTreeSet(path, (int) BURNIN_PERCENT);
        ts.reset();
        List<Tree> all = new ArrayList<>();
        while (ts.hasNext()) {
            all.add(ts.next());
        }
        List<Tree> out = new ArrayList<>();
        double step = Math.max(1.0, all.size() / (double) count);
        for (int i = 0; i < count && (int) (i * step) < all.size(); i++) {
            out.add(all.get((int) (i * step)));
        }
        return out;
    }

    /** Writes one subset as NEXUS with a Translate block and numeric newick. */
    private static File writeSubset(File dir, String base, String part, List<Tree> sub)
            throws Exception {
        File out = new File(dir, base + "." + part + ".trees");
        try (java.io.PrintWriter w = new java.io.PrintWriter(out)) {
            w.println("#NEXUS");
            w.println("Begin taxa;");
            w.println("\tDimensions ntax=" + taxaCache.size() + ";");
            w.println("\tTaxlabels");
            for (String t : taxaCache) {
                w.println("\t\t'" + t + "'");
            }
            w.println("\t\t;");
            w.println("End;");
            w.println("Begin trees;");
            w.println("\tTranslate");
            for (int i = 0; i < taxaCache.size(); i++) {
                w.println("\t\t" + (i + 1) + " '" + taxaCache.get(i) + "'"
                        + (i + 1 < taxaCache.size() ? "," : ""));
            }
            w.println("\t\t;");
            for (int i = 0; i < sub.size(); i++) {
                w.println("tree STATE_" + i + " = " + numericNewick(sub.get(i).getRoot()) + ";");
            }
            w.println("End;");
        }
        return out;
    }

    /** Newick using 1-based taxon numbers, matching the Translate block written above. */
    private static String numericNewick(beast.base.evolution.tree.Node v) {
        if (v.isLeaf()) {
            return String.valueOf(v.getNr() + 1);
        }
        StringBuilder sb = new StringBuilder("(");
        for (int i = 0; i < v.getChildCount(); i++) {
            sb.append(i > 0 ? "," : "").append(numericNewick(v.getChild(i)));
        }
        return sb.append(")").toString();
    }

    /**
     * Runs every prepared dataset in one JVM. A separate Maven invocation per dataset spends far
     * more time starting up than modelling: sweeping CRegCCD over 90 datasets took 11 minutes of
     * wall clock for 2.4 seconds of actual work.
     */
    @Test
    public void batchOverPreparedDatasets() throws Exception {
        String batch = System.getProperty("ccd.batch", "");
        Assumptions.assumeTrue(!batch.isEmpty());
        File dir = new File(batch);
        List<String> bases = new ArrayList<>();
        for (File f : dir.listFiles((d, n) -> n.endsWith(".train.trees"))) {
            bases.add(f.getName().replaceAll("\\.train\\.trees$", ""));
        }
        java.util.Collections.sort(bases);
        System.out.printf("batch: %d prepared datasets%n", bases.size());
        int ok = 0;
        for (String base : bases) {
            PATH = new File(dir, base + ".trees").getAbsolutePath();
            chainTrees = null;
            newickCache = null;
            taxaCache = null;
            CSV_CELLS.clear();
            try {
                headToHeadOnRealData();
                ok++;
            } catch (Throwable e) {
                System.out.printf("FAILED %s (%s: %s)%n", base,
                        e.getClass().getSimpleName(), String.valueOf(e.getMessage()));
                try (java.io.PrintWriter w = new java.io.PrintWriter(
                        new java.io.FileWriter(CSV + ".failed", true))) {
                    w.println(base + "," + e.getClass().getSimpleName());
                }
            }
        }
        System.out.printf("batch complete: %d/%d%n", ok, bases.size());
    }

    @Test
    public void headToHeadOnRealData() throws Exception {
        // in prepared mode the path names the dataset and locates its subsets; the original
        // chain file need not be present
        Assumptions.assumeTrue(!PATH.isEmpty()
                        && (new File(PATH).exists() || !PREPARED.isEmpty()),
                "set -Dccd.trees=/path/to/x.trees to run this comparison");

        List<Tree> probe = train(N);
        int nTaxa = probe.get(0).getLeafNodeCount();
        cell("dataset", new File(PATH).getName().replace(",", ";"));
        // which protocol produced the test set: a sibling MCMC run, or the second half of this chain
        cell("testProtocol", new File(PREPARED, new File(PATH).getName()
                .replaceAll("\\.trees$", "") + ".siblingtest").exists()
                ? "sibling-run" : "chain-half");
        cell("taxa", nTaxa);
        cell("chainTrees", chainTrees);
        cell("nTrain", probe.size());
        cell("nTest", N);
        // CCD1 entropy of the training sample, recorded for every dataset regardless of which
        // models run: it is the manuscript's measure of how much topological uncertainty a
        // posterior actually holds, and datasets below a few nats cannot discriminate the models.
        CCD1 entropyProbe = new CCD1(train(N), 0.0);
        cell("CCD1_entropy", String.format("%.3f", entropyProbe.getEntropy()));
        System.out.printf("%n=== %s: %d taxa, %d train / %d test trees ===%n",
                new File(PATH).getName(), nTaxa, probe.size(), N);

        List<Tree> test = PREPARED.isEmpty() ? read(N, 0.5, 1.0) : readPrepared("test");
        List<Tree> val = valSet(N / 2);
        int nFit = N / 2;

        // ---- CCD1 ----
        if (runs("CCD1")) {
        List<Tree> ccd1Trees = train(N);
        long tCcd1 = System.nanoTime();
        CCD1 ccd1 = new CCD1(ccd1Trees, 0.0);
        cell("CCD1_constructMs", ms(tCcd1));
        report("CCD1", "-", ccd1::getLogProbabilityOfTree, test);
        }

        // ---- RegCCD: alpha on validation ----
        if (runs("RegCCD")) {
        double bestAlpha = 0.4;
        double bestAlphaScore = Double.NEGATIVE_INFINITY;
        for (double alpha : new double[]{0.01, 0.05, 0.1, 0.2, 0.4, 0.8, 1.0}) {
            RegCCD m = new RegCCD(fitSet(nFit), 0.0, alpha);
            double sc = score(m::getLogProbabilityOfTree, val)[1];
            if (sc > bestAlphaScore) {
                bestAlphaScore = sc;
                bestAlpha = alpha;
            }
        }
        List<Tree> regTrees = train(N);
        long tReg = System.nanoTime();
        RegCCD reg = new RegCCD(regTrees, 0.0, bestAlpha);
        cell("RegCCD_constructMs", ms(tReg));
        report("RegCCD", String.format("alpha=%.2f", bestAlpha), reg::getLogProbabilityOfTree, test);
        }

        // ---- KRegCCD: mu on validation at alpha = 0.4 (the manuscript's setting) ----
        KRegCCD kreg = null;
        if (runs("KRegCCD")) {
        KRegCCD kFit = new KRegCCD(fitSet(nFit), 0.0, 0.005, 0.4);
        double bestMu = 0.005;
        double bestMuScore = Double.NEGATIVE_INFINITY;
        for (double mu : new double[]{0.00002, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05}) {
            final double m = mu;
            double sc = score(t -> kFit.getLogProbabilityOfTree(t, m), val)[1];
            if (sc > bestMuScore) {
                bestMuScore = sc;
                bestMu = mu;
            }
        }
        List<Tree> kregTrees = train(N);
        long tKreg = System.nanoTime();
        kreg = new KRegCCD(kregTrees, 0.0, bestMu, 0.4);
        kreg.precomputeReserves();   // KRegCCD defers its reserve solve; charge it to construction
        cell("KRegCCD_constructMs", ms(tKreg));
        report("KRegCCD", String.format("alpha=0.4, mu=%.5f", bestMu),
                kreg::getLogProbabilityOfTree, test);
        }

        // ---- MRegCCD: mu on validation (isolated: a baseline crash must not lose the row) ----
        MRegCCD mreg = null;
        try {
            if (!runs("MRegCCD")) {
                throw new IllegalStateException("skipped by ccd.models");
            }
            MRegCCD mFit = new MRegCCD(fitSet(nFit), 0.0, MRegCCD.DEFAULT_MU);
            double bestMMu = MRegCCD.DEFAULT_MU;
            double bestMMuScore = Double.NEGATIVE_INFINITY;
            for (double mu : new double[]{0.00005, 0.0002, 0.001, 0.002, 0.008, 0.0159, 0.05, 0.1}) {
                final double m = mu;
                double sc = score(t -> mFit.getLogProbabilityOfTree(t, m), val)[1];
                if (sc > bestMMuScore) {
                    bestMMuScore = sc;
                    bestMMu = mu;
                }
            }
            List<Tree> mregTrees = train(N);
            long tMreg = System.nanoTime();
            mreg = new MRegCCD(mregTrees, 0.0, bestMMu);
            cell("MRegCCD_constructMs", ms(tMreg));
            report("MRegCCD", String.format("mu=%.5f", bestMMu), mreg::getLogProbabilityOfTree, test);
        } catch (Throwable e) {
            mreg = null;
            if (runs("MRegCCD")) System.out.printf("MRegCCD   SKIPPED (%s: %s)%n",
                    e.getClass().getSimpleName(), String.valueOf(e.getMessage()));
        }

        // ---- CRegCCD: (alpha, alpha1, alpha2) on validation ----
        CRegCCD creg = null;
        if (runs("CRegCCD")) {
        CRegCCD cFit = new CRegCCD(fitSet(nFit), 0.0);
        double[] grid = {0.002, 0.01, 0.05, 0.2, 0.4, 1.0, 2.0, 5.0, 12.0};
        double[] bestC = {0.0, 0.4, 0.4, 0.4};
        double bestCScore = Double.NEGATIVE_INFINITY;
        for (double b2 : grid) {
            for (double b3 : grid) {
                for (double b4 : grid) {
                    double sc = score(t -> cFit.getLogProbabilityOfTree(t, b2, b3, b4), val)[1];
                    if (sc > bestCScore) {
                        bestCScore = sc;
                        bestC = new double[]{0.0, b2, b3, b4};
                    }
                }
            }
        }
        List<Tree> cregTrees = train(N);
        long tCreg = System.nanoTime();
        creg = new CRegCCD(cregTrees, 0.0, bestC[1], bestC[2], bestC[3]);
        cell("CRegCCD_constructMs", ms(tCreg));
        report("CRegCCD", String.format("alpha=%.3f, alpha1=%.3f, alpha2=%.3f", bestC[1], bestC[2], bestC[3]),
                creg::getLogProbabilityOfTree, test);
        }

        // ---- paired comparison of the full-support models against KRegCCD ----
        if (kreg == null) {
            writeCsv();
            return;   // paired statistics are computed from the per-tree dumps instead
        }
        System.out.printf("%npaired per-tree comparison against KRegCCD (n = %d test trees):%n", test.size());
        if (mreg != null) {
            paired("MRegCCD", mreg::getLogProbabilityOfTree, kreg::getLogProbabilityOfTree, test);
        }
        if (creg != null) {
            paired("CRegCCD", creg::getLogProbabilityOfTree, kreg::getLogProbabilityOfTree, test);
        }

        writeCsv();
    }

    private static void report(String name, String params, Scorer s, List<Tree> test)
            throws Exception {
        dumpPerTree(name, s, test);
        long t0 = System.nanoTime();
        double[] sc = score(s, test);
        long scoreMs = ms(t0);            // one pass over the test set
        double all = meanAll(s, test);
        cell(name + "_scoreMs", scoreMs);
        cell(name + "_params", params.replace(",", ";"));
        cell(name + "_coverage", String.format("%.4f", sc[0] / test.size()));
        cell(name + "_meanLogP", Double.isFinite(all) ? String.format("%.4f", all) : "-inf");
        System.out.printf("%-9s %-27s coverage %6.1f%%   mean logP(covered) %10.2f   mean logP(all) %s%n",
                name, params, 100.0 * sc[0] / test.size(), sc[1],
                Double.isFinite(all) ? String.format("%10.2f", all) : "      -inf");
    }

    private static void paired(String name, Scorer a, Scorer baseline, List<Tree> test) {
        // records name_vs_KRegCCD_{diff,se,wins,n} in addition to printing
        int wins = 0;
        double sumDiff = 0.0;
        double sumSq = 0.0;
        for (Tree t : test) {
            double d = a.logP(t) - baseline.logP(t);
            if (d > 0) {
                wins++;
            }
            sumDiff += d;
            sumSq += d * d;
        }
        int n = test.size();
        double mean = sumDiff / n;
        double se = Math.sqrt(Math.max(0, sumSq / n - mean * mean) / n);
        cell(name + "_vs_KRegCCD_diff", String.format("%.4f", mean));
        cell(name + "_vs_KRegCCD_se", String.format("%.4f", se));
        cell(name + "_vs_KRegCCD_wins", wins);
        cell(name + "_vs_KRegCCD_n", n);
        System.out.printf("  %-9s mean log-ratio %+8.2f +/- %.2f nats/tree, better on %d/%d trees%n",
                name, mean, se, wins, n);
    }

    /** Scale check: sampling from a real 129-taxon posterior must be fast, valid and
     *  self-consistent with the scorer. */
    @Test
    public void samplingOnRealData() throws Exception {
        // in prepared mode the path names the dataset and locates its subsets; the original
        // chain file need not be present
        Assumptions.assumeTrue(!PATH.isEmpty()
                        && (new File(PATH).exists() || !PREPARED.isEmpty()),
                "set -Dccd.trees=/path/to/x.trees to run this comparison");
        CRegCCD ccd = new CRegCCD(train(N), 0.0);
        int nTaxa = ccd.getSizeOfLeavesArray();
        long t0 = System.nanoTime();
        int draws = 200;
        for (int i = 0; i < draws; i++) {
            Tree t = ccd.sampleTree();
            if (t.getLeafNodeCount() != nTaxa || t.getNodeCount() != 2 * nTaxa - 1) {
                throw new AssertionError("invalid sampled tree");
            }
            double stamped = (Double) t.getRoot().getMetaData(CCD1.LOG_PROB_SUBTREE_KEY);
            double scored = ccd.getLogProbabilityOfTree(t);
            if (Math.abs(stamped - scored) > 1e-9) {
                throw new AssertionError("stamped " + stamped + " != scored " + scored);
            }
        }
        double secs = (System.nanoTime() - t0) / 1e9;
        System.out.printf("CRegCCD sampling on %s: %d taxa, %d trees in %.2f s (%.1f ms/tree), "
                + "all valid and self-consistent%n",
                new File(PATH).getName(), nTaxa, draws, secs, 1000 * secs / draws);
    }

    /** MAP and entropy on a real posterior: correctness certificate and wall-clock cost. */
    @Test
    public void mapAndEntropyOnRealData() throws Exception {
        // in prepared mode the path names the dataset and locates its subsets; the original
        // chain file need not be present
        Assumptions.assumeTrue(!PATH.isEmpty()
                        && (new File(PATH).exists() || !PREPARED.isEmpty()),
                "set -Dccd.trees=/path/to/x.trees to run this comparison");
        CRegCCD creg = new CRegCCD(train(N), 0.0, 2.0, 0.4, 0.05);

        long t0 = System.nanoTime();
        double maxLog = creg.getMaxLogTreeProbability();
        boolean certified = creg.isMAPCertifiedGlobal();
        double bound = creg.getOffBackboneBound();
        double mapSecs = (System.nanoTime() - t0) / 1e9;

        t0 = System.nanoTime();
        Tree map = creg.getMAPTree();
        double treeSecs = (System.nanoTime() - t0) / 1e9;
        double scored = creg.getLogProbabilityOfTree(map);

        t0 = System.nanoTime();
        double[] h = creg.getEntropyMonteCarlo(20_000);
        double entSecs = (System.nanoTime() - t0) / 1e9;

        System.out.printf("%n=== %s: CRegCCD MAP and entropy ===%n", new File(PATH).getName());
        System.out.printf("MAP DP + certificate : %.2f s, max logP = %.4f, "
                        + "off-backbone bound = %.4f, certified global = %s%n",
                mapSecs, maxLog, bound, certified);
        System.out.printf("MAP tree traceback   : %.2f s, scored logP = %.4f (matches: %s)%n",
                treeSecs, scored, Math.abs(scored - maxLog) < 1e-9);
        System.out.printf("entropy (20k draws)  : %.2f s, H = %.3f +/- %.3f nats%n",
                entSecs, h[0], h[1]);
    }

    /** Deterministic entropy recursion vs the unbiased Monte-Carlo estimator on a real posterior. */
    @Test
    public void entropyRecursionVersusMonteCarloOnRealData() throws Exception {
        // in prepared mode the path names the dataset and locates its subsets; the original
        // chain file need not be present
        Assumptions.assumeTrue(!PATH.isEmpty()
                        && (new File(PATH).exists() || !PREPARED.isEmpty()),
                "set -Dccd.trees=/path/to/x.trees to run this comparison");
        System.out.printf("%n=== %s: CRegCCD entropy, recursion vs Monte Carlo ===%n",
                new File(PATH).getName());
        System.out.printf("%-22s %-11s %-9s %-22s %-9s %-10s%n",
                "pseudocounts", "recursion", "rec (s)", "Monte Carlo", "MC (s)", "difference");
        for (double[] p : new double[][]{{2.0, 0.4, 0.05}, {5.0, 2.0, 0.4}, {0.4, 0.4, 0.4}}) {
            CRegCCD ccd = new CRegCCD(train(N), 0.0, p[0], p[1], p[2]);
            long t0 = System.nanoTime();
            double rec = ccd.getEntropyRecursive();
            double recSecs = (System.nanoTime() - t0) / 1e9;
            t0 = System.nanoTime();
            double[] mc = ccd.getEntropyMonteCarlo(200_000);
            double mcSecs = (System.nanoTime() - t0) / 1e9;
            double diff = rec - mc[0];
            System.out.printf("a=(%.2f,%.2f,%.2f)%-6s %-11.4f %-9.2f %8.4f +/-%.4f %-9.2f %+.4f (%+.3f%%)%n",
                    p[0], p[1], p[2], "", rec, recSecs, mc[0], mc[1], mcSecs, diff, 100 * diff / mc[0]);
        }
    }

    /** Does the exact A_0/A_1 MAP search stay tractable on a real posterior? */
    @Test
    public void exactMapOnRealData() throws Exception {
        // in prepared mode the path names the dataset and locates its subsets; the original
        // chain file need not be present
        Assumptions.assumeTrue(!PATH.isEmpty()
                        && (new File(PATH).exists() || !PREPARED.isEmpty()),
                "set -Dccd.trees=/path/to/x.trees to run this comparison");
        CRegCCD creg = new CRegCCD(train(N), 0.0, 0.4, 0.4, 0.05);
        double backbone = creg.getMaxLogTreeProbability();
        System.out.printf("%n=== %s: MAP search by allowed A_1 depth ===%n", new File(PATH).getName());
        System.out.printf("backbone (A_0 only) max logP = %.4f%n", backbone);
        System.out.printf("%-6s %-12s %-12s %-12s %-10s %-8s%n",
                "maxA1", "max logP", "improvement", "bound", "certified", "states");
        for (int k = 0; k <= 3; k++) {
            long t0 = System.nanoTime();
            CRegCCD.MapResult r = creg.solveMAP(k);
            double secs = (System.nanoTime() - t0) / 1e9;
            if (!r.complete()) {
                System.out.printf("%-6d exceeded the state budget after %d states (%.1f s)%n",
                        k, r.statesExplored(), secs);
                break;
            }
            System.out.printf("%-6d %-12.4f %-12.4f %-12.4f %-10s %-8d (%.1f s)%n",
                    k, r.maxLogProbability(), r.maxLogProbability() - backbone,
                    r.offBackboneBound(), r.a2Excluded(), r.statesExplored(), secs);
        }
    }
}
