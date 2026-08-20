package ccd.tools;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator;
import ccd.model.*;

import java.io.File;
import java.util.*;

/**
 * Measures what CRegCCD, KRegCCD and MRegCCD charge for each of CRegCCD's four split classes.
 *
 * <pre>
 *   java -cp ... ccd.tools.SplitClassCost &lt;prepared-dir&gt; &lt;base&gt; alpha alpha1 alpha2 muK muM
 * </pre>
 *
 *
 * The classes are a property of the split and the training set, not of any model, so they are
 * computed here directly from the training trees rather than from any model's internals:
 *   1 observed  -- the split of C was seen in training
 *   2 expanded  -- C and both parts are observed clades, but never as this split of C
 *   3 one novel -- exactly one part is an observed clade
 *   4 two novel -- neither part is an observed clade
 *
 * Each held-out tree is then summarised by how many splits of each class it contains, and the
 * per-tree log probability under each model is regressed on those counts. KRegCCD and MRegCCD
 * price maximal novel *regions* rather than individual splits, so the fit is not exact by
 * construction; the residual scatter is reported so that the extent of the approximation is
 * visible rather than assumed away.
 */
public class SplitClassCost {

    static List<Tree> read(File dir, String base, String part) throws Exception {
        File f = new File(dir, base + "." + part + ".trees");
        TreeAnnotator.MemoryFriendlyTreeSet s =
                new TreeAnnotator().new MemoryFriendlyTreeSet(f.getAbsolutePath(), 0);
        s.reset();
        List<Tree> out = new ArrayList<>();
        Tree t;
        while ((t = s.next()) != null) out.add(t);
        return out;
    }

    /** Leaf-index bitset of the clade below a node, as a sorted string key. */
    static String clade(Node n, Map<Node, java.util.BitSet> memo, int nTaxa) {
        return bits(n, memo, nTaxa).toString();
    }

    static java.util.BitSet bits(Node n, Map<Node, java.util.BitSet> memo, int nTaxa) {
        java.util.BitSet b = memo.get(n);
        if (b != null) return b;
        b = new java.util.BitSet(nTaxa);
        if (n.isLeaf()) {
            b.set(n.getNr());
        } else {
            for (Node c : n.getChildren()) b.or(bits(c, memo, nTaxa));
        }
        memo.put(n, b);
        return b;
    }

    static void collect(Tree t, Set<String> clades, Set<String> splits, int nTaxa) {
        Map<Node, java.util.BitSet> memo = new HashMap<>();
        for (Node n : t.getNodesAsArray()) {
            if (n.isLeaf()) continue;
            clades.add(clade(n, memo, nTaxa));
            List<Node> ch = n.getChildren();
            if (ch.size() == 2) {
                String l = clade(ch.get(0), memo, nTaxa), r = clade(ch.get(1), memo, nTaxa);
                splits.add(clade(n, memo, nTaxa) + "|" + (l.compareTo(r) < 0 ? l + "," + r : r + "," + l));
            }
        }
    }

    /** Size bins for the parent clade: the analytic claim is about how the penalty scales with
     *  the number of taxa under the clade being split, so the bins are geometric. */
    static final int[] BIN_HI = {4, 8, 16, 32, 64, Integer.MAX_VALUE};
    static final String[] BIN_NAME = {"2-4", "5-8", "9-16", "17-32", "33-64", "65+"};
    /** mean parent-clade size actually seen in each (class, bin) cell, so the fitted penalty can be
     *  checked against an analytic prediction in m rather than against a bin label. */
    static double[] sumM = new double[4 * 6];

    static int bin(int m) {
        for (int i = 0; i < BIN_HI.length; i++) if (m <= BIN_HI[i]) return i;
        return BIN_HI.length - 1;
    }

    /** counts of the four classes crossed with parent-clade size bin, for one held-out tree. */
    static int[] classify(Tree t, Set<String> clades, Set<String> splits, int nTaxa) {
        Map<Node, java.util.BitSet> memo = new HashMap<>();
        int[] c = new int[4 * BIN_HI.length];
        for (Node n : t.getNodesAsArray()) {
            if (n.isLeaf() || n.getChildren().size() != 2) continue;
            String C = clade(n, memo, nTaxa);
            String l = clade(n.getChildren().get(0), memo, nTaxa);
            String r = clade(n.getChildren().get(1), memo, nTaxa);
            boolean lo = clades.contains(l) || n.getChildren().get(0).isLeaf();
            boolean ro = clades.contains(r) || n.getChildren().get(1).isLeaf();
            String key = C + "|" + (l.compareTo(r) < 0 ? l + "," + r : r + "," + l);
            int cls = splits.contains(key) ? 0 : (lo && ro) ? 1 : (lo || ro) ? 2 : 3;
            int m = bits(n, memo, nTaxa).cardinality();
            int idx = cls * BIN_HI.length + bin(m);
            c[idx]++;
            sumM[idx] += m;
        }
        return c;
    }

    /** A node is blue for KRegCCD/MRegCCD purposes when its split introduces a novel clade. */
    static boolean isBlue(Node n, Set<String> clades, Map<Node, java.util.BitSet> memo, int nTaxa) {
        if (n.isLeaf() || n.getChildren().size() != 2) return false;
        for (Node ch : n.getChildren()) {
            if (!ch.isLeaf() && !clades.contains(clade(ch, memo, nTaxa))) return true;
        }
        return false;
    }

    /**
     * Least squares WITHOUT an intercept; returns {b1..bk, R2}.
     *
     * The four class counts sum to the number of internal nodes, which is the same for every tree
     * on a fixed taxon set, so a design matrix with an intercept is singular and the split between
     * intercept and coefficients is arbitrary. Without the intercept each coefficient is directly
     * the mean log probability contributed by one split of that class, which is the quantity of
     * interest. R2 is taken about zero accordingly.
     */
    static double[] ols(double[][] X, double[] y) {
        int n = y.length, p = X[0].length;
        double[][] A = new double[p][p];
        double[] b = new double[p];
        for (int i = 0; i < n; i++) {
            double[] xi = X[i];
            for (int j = 0; j < p; j++) {
                b[j] += xi[j] * y[i];
                for (int k = 0; k < p; k++) A[j][k] += xi[j] * xi[k];
            }
        }
        for (int i = 0; i < p; i++) A[i][i] += 1e-9;
        double[] beta = solve(A, b);
        double ssTot = 0, ssRes = 0;
        for (int i = 0; i < n; i++) {
            double pred = 0;
            for (int j = 0; j < p; j++) pred += beta[j] * X[i][j];
            ssRes += (y[i] - pred) * (y[i] - pred);
            ssTot += y[i] * y[i];
        }
        double[] out = new double[p + 1];
        System.arraycopy(beta, 0, out, 0, p);
        out[p] = 1 - ssRes / ssTot;
        return out;
    }

    static double[] solve(double[][] A, double[] b) {
        int n = b.length;
        double[][] M = new double[n][n + 1];
        for (int i = 0; i < n; i++) {
            System.arraycopy(A[i], 0, M[i], 0, n);
            M[i][n] = b[i];
        }
        for (int col = 0; col < n; col++) {
            int piv = col;
            for (int r = col + 1; r < n; r++) if (Math.abs(M[r][col]) > Math.abs(M[piv][col])) piv = r;
            double[] tmp = M[col]; M[col] = M[piv]; M[piv] = tmp;
            for (int r = 0; r < n; r++) {
                if (r == col || M[col][col] == 0) continue;
                double f = M[r][col] / M[col][col];
                for (int c2 = col; c2 <= n; c2++) M[r][c2] -= f * M[col][c2];
            }
        }
        double[] x = new double[n];
        for (int i = 0; i < n; i++) x[i] = M[i][n] / M[i][i];
        return x;
    }

    public static void main(String[] args) throws Exception {
        File dir = new File(args[0]);
        String base = args[1];
        double alpha = Double.parseDouble(args[2]), a1 = Double.parseDouble(args[3]);
        double a2 = Double.parseDouble(args[4]);
        double muK = Double.parseDouble(args[5]), muM = Double.parseDouble(args[6]);

        List<Tree> test = read(dir, base, "test");
        int nTaxa = test.get(0).getLeafNodeCount();
        Set<String> clades = new HashSet<>(), splits = new HashSet<>();
        for (Tree t : read(dir, base, "train")) collect(t, clades, splits, nTaxa);

        int n = test.size();
        int NB = BIN_HI.length, P = 4 * NB;
        double[][] X = new double[n][P];
        int[] tot = new int[P];
        for (int i = 0; i < n; i++) {
            int[] c = classify(test.get(i), clades, splits, nTaxa);
            for (int j = 0; j < P; j++) { X[i][j] = c[j]; tot[j] += c[j]; }
        }
        System.out.printf("%s: %d taxa, %d held-out trees, %d internal splits each%n",
                base, nTaxa, n, nTaxa - 1);
        String[] CLS = {"observed", "expanded", "one-novel", "two-novel"};
        int tt = n * (nTaxa - 1);
        System.out.print("split classes:");
        for (int k = 0; k < 4; k++) {
            int sum = 0;
            for (int bnd = 0; bnd < NB; bnd++) sum += tot[k * NB + bnd];
            System.out.printf("  %s %d (%.2f%%)", CLS[k], sum, 100.0 * sum / tt);
        }
        System.out.println();

        Map<String, double[]> ys = new LinkedHashMap<>();
        CRegCCD cre = new CRegCCD(read(dir, base, "train"), 0.0, alpha, a1, a2);
        KRegCCD kre = new KRegCCD(read(dir, base, "train"), 0.0, muK, 0.4);
        MRegCCD mre = new MRegCCD(read(dir, base, "train"), 0.0, muM);
        double[] yc = new double[n], yk = new double[n], ym = new double[n];
        for (int i = 0; i < n; i++) {
            yc[i] = cre.getLogProbabilityOfTree(test.get(i));
            yk[i] = kre.getLogProbabilityOfTree(test.get(i));
            ym[i] = mre.getLogProbabilityOfTree(test.get(i));
        }
        ys.put("CRegCCD", yc); ys.put("KRegCCD", yk); ys.put("MRegCCD", ym);

        System.out.printf("%nmean log probability contributed per split, by class and parent clade size%n");
        for (Map.Entry<String, double[]> e : ys.entrySet()) {
            double[] r = ols(X, e.getValue());
            // OLS makes residuals orthogonal to every column, and the columns sum to a constant,
            // so the fitted mean must equal the actual mean. Printing both is the check that the
            // decomposition is real rather than a plausible-looking artefact of the solver.
            double actual = 0, fitted = 0;
            for (int i = 0; i < n; i++) {
                actual += e.getValue()[i] / n;
                for (int j = 0; j < P; j++) fitted += r[j] * X[i][j] / n;
            }
            System.out.printf("%n%s  (mean logP %.2f, fitted %.2f, R2 %.4f)%n",
                    e.getKey(), actual, fitted, r[P]);
            System.out.printf("  %-10s", "class\\size");
            for (String bn : BIN_NAME) System.out.printf("%9s", bn);
            System.out.println();
            for (int k = 0; k < 4; k++) {
                System.out.printf("  %-10s", CLS[k]);
                for (int bnd = 0; bnd < NB; bnd++) {
                    if (tot[k * NB + bnd] < 20) System.out.printf("%9s", ".");
                    else System.out.printf("%9.2f", r[k * NB + bnd]);
                }
                System.out.println();
            }
            if (e.getKey().equals("CRegCCD")) {
                System.out.printf("  %-10s", "mean m");
                for (int bnd = 0; bnd < NB; bnd++) {
                    int cnt = tot[3 * NB + bnd];
                    if (cnt < 20) System.out.printf("%9s", ".");
                    else System.out.printf("%9.1f", sumM[3 * NB + bnd] / cnt);
                }
                System.out.println("   <- for the two-novel row");
            }
        }
        // How many terms does each model actually score? CRegCCD is a chain-rule CCD over every
        // clade in the tree, so it contributes one factor per internal node -- including the novel
        // clades themselves. KRegCCD and MRegCCD cluster novel nodes into maximal regions and
        // score each region once at its top, so the nodes interior to a region never receive a
        // conditional split distribution at all. Counting the nodes is exact, unlike regressing
        // the score on class counts: the interior count is largely determined by the novel splits
        // above it, so those columns are collinear and their coefficients are not identifiable.
        long internal = 0, novelClade = 0, regionTop = 0, blue = 0;
        for (Tree t : test) {
            Map<Node, java.util.BitSet> memo = new HashMap<>();
            for (Node nd : t.getNodesAsArray()) {
                if (nd.isLeaf() || nd.getChildren().size() != 2) continue;
                internal++;
                boolean isB = isBlue(nd, clades, memo, nTaxa);
                if (isB) blue++;
                if (!clades.contains(clade(nd, memo, nTaxa)) && !nd.isRoot()) novelClade++;
                if (isB) {
                    Node par = nd.getParent();
                    if (par == null || !isBlue(par, clades, memo, nTaxa)) regionTop++;
                }
            }
        }
        double perTree = 1.0 / test.size();
        System.out.printf("%nterms scored per held-out tree (%d internal nodes each)%n", nTaxa - 1);
        System.out.printf("  internal nodes                      %8.2f%n", internal * perTree);
        System.out.printf("  novel clades (no observed counterpart) %5.2f%n", novelClade * perTree);
        System.out.printf("  blue nodes (split introduces a novel clade) %2.2f%n", blue * perTree);
        System.out.printf("  maximal novel regions               %8.2f%n", regionTop * perTree);
        System.out.printf("  CRegCCD factors                     %8.2f  (one per internal node)%n",
                internal * perTree);
        System.out.printf("  KRegCCD / MRegCCD factors           %8.2f  (red nodes + one per region)%n",
                (internal - blue + regionTop) * perTree);
        System.out.printf("  difference                          %8.2f  factors per tree that CRegCCD%n",
                (blue - regionTop) * perTree);
        System.out.printf("%36s  pays and the others do not%n", "");

        System.out.println();
        System.out.println("Coefficients are nats per split of that class; more negative is a heavier");
        System.out.println("penalty. CRegCCD spreads the class total alpha2 uniformly over |A_4| ~ 2^(m-1)");
        System.out.println("splits, so its two-novel penalty should fall linearly in m with slope -log 2 =");
        System.out.println("-0.693 nats per taxon; KRegCCD and MRegCCD price the same split as a power of");
        System.out.println("an escape rate eps ~ mu / O(s^2), which is logarithmic in m. Compare the");
        System.out.println("two-novel row against the mean-m row to check that.");
    }
}
