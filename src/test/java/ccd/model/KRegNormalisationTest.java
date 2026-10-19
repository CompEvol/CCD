package ccd.model;

import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Total probability mass of {@link KRegCCD}, by brute-force enumeration of every rooted topology.
 *
 * <p>Documents the model's normalisation exactly as the manuscript states it:
 * <ul>
 *   <li>on four taxa the model is exactly normalised (no blue-region boundary part is itself a
 *       reserving clade, so the region decomposition is tight);</li>
 *   <li>from six taxa on it is <em>sub</em>-normalised -- total mass is below 1, never above --
 *       and the deficit scales as {@code mu^2}, which is the maximality deficit rather than the
 *       {@code O(mu^(k+1))} reserve truncation (it persists at full reserve depth).</li>
 * </ul>
 */
public class KRegNormalisationTest {

    private static final List<String> TAXA4 = Arrays.asList("A", "B", "C", "D");
    private static final List<String> TAXA5 = Arrays.asList("A", "B", "C", "D", "E");
    private static final List<String> TAXA6 = Arrays.asList("A", "B", "C", "D", "E", "F");
    private static final List<String> TAXA7 = Arrays.asList("A", "B", "C", "D", "E", "F", "G");

    private static List<Tree> trees(List<String> taxa, String... newicks) {
        List<Tree> out = new ArrayList<>();
        for (String nwk : newicks) {
            out.add(new TreeParser(taxa, nwk, 1, false));
        }
        return out;
    }

    private static List<Tree> training4() {
        return trees(TAXA4, "(((A:1,B:1):1,C:1):1,D:1):0;", "((A:1,B:1):1,(C:1,D:1):1):0;");
    }

    private static List<Tree> training6() {
        return trees(TAXA6,
                "(((((A:1,B:1):1,C:1):1,D:1):1,E:1):1,F:1):0;",
                "((((D:1,C:1):1,B:1):1,A:1):1,(E:1,F:1):1):0;");
    }

    private static List<Tree> training7() {
        return trees(TAXA7,
                "((((((A:1,B:1):1,C:1):1,D:1):1,E:1):1,F:1):1,G:1):0;",
                "(((((D:1,C:1):1,B:1):1,A:1):1,(E:1,F:1):1):1,G:1):0;");
    }

    /** Total mass under the full-support score. */
    private static double totalMass(KRegCCD ccd, List<String> taxa) {
        double mass = 0.0;
        for (Tree t : allRootedTopologies(taxa)) {
            mass += Math.exp(ccd.getLogProbabilityOfTree(t));
        }
        return mass;
    }

    @Test
    public void exactlyNormalisedOnFourTaxa() {
        for (double mu : new double[]{0.001, 0.005, 0.05}) {
            KRegCCD ccd = new KRegCCD(training4(), 0.0, mu, 0.4, 2, KRegCCD.TailMode.NONE,
                    KRegCCD.NovelMode.FLAT);
            double mass = totalMass(ccd, TAXA4);
            System.out.printf("KRegCCD 4 taxa mu=%-6.3f totalMass = %.12f%n", mu, mass);
            assertEquals(1.0, mass, 1e-9, "four-taxon model must normalise exactly");
        }
    }

    /** Five taxa: exactness is not a taxon-count property but depends on whether any escape region
     *  has a reserving clade on its boundary, which the training set determines. */
    @Test
    public void fiveTaxaNormalisationDependsOnTrainingSet() {
        List<Tree> caterpillar = trees(TAXA5, "((((A:1,B:1):1,C:1):1,D:1):1,E:1):0;");
        List<Tree> mixed = trees(TAXA5,
                "((((A:1,B:1):1,C:1):1,D:1):1,E:1):0;",
                "(((A:1,B:1):1,(C:1,D:1):1):1,E:1):0;");
        for (double mu : new double[]{0.005, 0.05}) {
            for (String name : new String[]{"caterpillar", "mixed"}) {
                List<Tree> training = name.equals("caterpillar") ? caterpillar : mixed;
                KRegCCD ccd = new KRegCCD(training, 0.0, mu, 0.4, 2, KRegCCD.TailMode.NONE,
                        KRegCCD.NovelMode.FLAT);
                double mass = totalMass(ccd, TAXA5);
                System.out.printf("KRegCCD 5 taxa mu=%-6.3f %-12s totalMass = %.9f  (1-mass = %+.3e)%n",
                        mu, name, mass, 1.0 - mass);
                assertTrue(mass <= 1.0 + 1e-12, "must never super-normalise, got " + mass);
            }
        }
    }

    @Test
    public void subNormalisedFromSixTaxaWithMuSquaredDeficit() {
        for (List<String> taxa : List.of(TAXA6, TAXA7)) {
            List<Tree> training = (taxa == TAXA6) ? training6() : training7();
            double prevDeficit = Double.NaN;
            double prevMu = Double.NaN;
            for (double mu : new double[]{0.05, 0.005, 0.001}) {
                for (KRegCCD.TailMode tm : KRegCCD.TailMode.values()) {
                    KRegCCD ccd = new KRegCCD(training, 0.0, mu, 0.4, 2, tm,
                            KRegCCD.NovelMode.FLAT);
                    double mass = totalMass(ccd, taxa);
                    double deficit = 1.0 - mass;
                    System.out.printf("KRegCCD %d taxa mu=%-6.3f %-8s totalMass = %.9f  (1-mass = %+.3e)%n",
                            taxa.size(), mu, tm, mass, deficit);
                    assertTrue(mass <= 1.0 + 1e-12,
                            "model must never super-normalise, got " + mass);
                    if (tm == KRegCCD.TailMode.NONE) {
                        if (!Double.isNaN(prevDeficit)) {
                            // a mu^2 deficit shrinks by the square of the mu ratio
                            double ratio = prevDeficit / deficit;
                            double expected = (prevMu / mu) * (prevMu / mu);
                            System.out.printf("    deficit shrank %.1fx as mu fell %.0fx "
                                            + "(mu^2 predicts %.0fx)%n",
                                    ratio, prevMu / mu, expected);
                            assertEquals(expected, ratio, 0.25 * expected,
                                    "deficit must scale as mu^2");
                        }
                        prevDeficit = deficit;
                        prevMu = mu;
                    }
                }
            }
        }
    }

    private static List<Tree> allRootedTopologies(List<String> taxa) {
        List<Tree> trees = new ArrayList<>();
        for (String shape : shapes(taxa)) {
            trees.add(new TreeParser(taxa, shape + ";", 1, false));
        }
        return trees;
    }

    private static List<String> shapes(List<String> taxa) {
        List<String> out = new ArrayList<>();
        if (taxa.size() == 1) {
            out.add(taxa.get(0) + ":1");
            return out;
        }
        String first = taxa.get(0);
        List<String> rest = taxa.subList(1, taxa.size());
        int n = rest.size();
        for (int mask = 0; mask < (1 << n); mask++) {
            List<String> left = new ArrayList<>();
            left.add(first);
            List<String> right = new ArrayList<>();
            for (int i = 0; i < n; i++) {
                if ((mask & (1 << i)) != 0) {
                    left.add(rest.get(i));
                } else {
                    right.add(rest.get(i));
                }
            }
            if (right.isEmpty()) {
                continue;
            }
            for (String l : shapes(left)) {
                for (String r : shapes(right)) {
                    out.add("(" + l + "," + r + "):1");
                }
            }
        }
        return out;
    }
}
