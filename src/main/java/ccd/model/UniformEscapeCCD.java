package ccd.model;

import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator.TreeSet;

import java.math.BigInteger;
import java.util.List;

/**
 * The <em>trivial</em> full-support baseline: a plain {@link CCD1} backbone with a single lump of
 * escape mass {@code epsilon} spread <em>uniformly</em> over every rooted topology the backbone does
 * not already cover.
 *
 * <p>
 * It exists to make one point precise. A {@link CCD1} (sampled splits only) assigns probability
 * <em>exactly zero</em> to any tree containing an unseen clade or split, so on held-out data it has a
 * coverage hole and an infinite cross-entropy. Full coverage <em>by itself</em> is trivial to buy:
 * reserve a total mass {@code epsilon} and divide it equally among the
 * {@code M = (2n-3)!! - K} topologies outside the backbone's support, where {@code K} is the number
 * of topologies the backbone covers and {@code (2n-3)!!} is the number of rooted topologies on
 * {@code n} taxa. Concretely, for a tree {@code T},
 * <pre>
 *   P(T) = (1 - epsilon) * P_CCD1(T)   if T is in the backbone's support,
 *        = epsilon / M                  otherwise.
 * </pre>
 * The mass is exactly normalised: {@code (1 - epsilon) * 1 + M * (epsilon / M) = 1}.
 *
 * <p>
 * Because every off-support tree shares the same minimal probability {@code epsilon / M}, this model
 * is a deliberate <em>strawman</em>. Under the PIT / credible-level statistic
 * {@code u(T) = P_{S~q}[ q(S) >= q(T) ]}, every held-out tree with novel structure scores the global
 * minimum probability and so lands at {@code u ≈ 1} — they all pile into the top PIT bin (for
 * {@code epsilon} smaller than the bin width), producing a calibration spike. And on held-out
 * likelihood each novel tree costs {@code log(epsilon) - log(M) ≈ -log((2n-3)!!)} nats — finite, but
 * catastrophically negative because the escape budget is smeared uniformly over an astronomically
 * large set instead of being concentrated near the observed trees.
 *
 * <p>
 * That contrast is the point: {@link KRegCCD} reaches full support with the <em>same</em> kind of
 * escape budget but places it <em>structurally</em> (geometric in tree distance, concentrated near
 * the sample), which both flattens the PIT and keeps the held-out likelihood high. This class is the
 * uniform-placement control that isolates "structural placement" as the source of KRegCCD's
 * advantage.
 *
 * <h3>Scope and assumptions</h3>
 * The topology count {@code (2n-3)!!} is for fully resolved rooted bifurcating trees on {@code n}
 * labelled taxa (no sampled ancestors, no polytomies) — the setting of the PIT / RSV2 held-out
 * experiments. The escape probability is reported via the log, since {@code M} overflows
 * {@code double} for even modest {@code n}; the linear {@link #getProbabilityOfTree(Tree)} of an
 * off-support tree therefore underflows to {@code 0} for large {@code n} while
 * {@link #getLogProbabilityOfTree(Tree)} stays finite (that is what the held-out / PIT drivers use).
 *
 * <h3>Sampling</h3>
 * {@link #sampleTreeLogProbability()} and {@link #sampleTreeProbability()} correctly sample the full
 * mixture (an {@code epsilon} fraction returns the escape value), which is all the PIT reference
 * distribution needs. {@link #sampleTree()} is inherited and returns a concrete topology from the
 * CCD1 backbone only — it does not synthesise the astronomically many escape topologies. This baseline
 * is a held-out <em>scoring</em> device, not a generative model.
 *
 * @author Alexei Drummond
 */
public class UniformEscapeCCD extends CCD1 {

    /** Total escape mass reserved for (and spread uniformly over) the off-support topologies. */
    private final double escapeMass;

    /** Precomputed {@code log(1 - escapeMass)}, the discount applied to every in-support tree. */
    private final double log1mEscapeMass;

    /**
     * Cached escape log-probability {@code log(epsilon) - log(M)} and the backbone topology count
     * {@code K} it was computed for; recomputed lazily if {@code K} changes (e.g. trees added later).
     * {@link Double#NEGATIVE_INFINITY} means there are no off-support topologies, so the model
     * collapses to the plain CCD1 backbone.
     */
    private double escapeLogProbCache = Double.NaN;
    private BigInteger cachedSupportSize = null;

    /* -- CONSTRUCTORS -- */

    /**
     * Constructor from a collection of trees with burn-in.
     *
     * @param trees      the trees whose CCD1 backbone is escape-smoothed
     * @param burnin     fraction in {@code [0,1]} of leading trees discarded as burn-in
     * @param escapeMass total escape mass {@code epsilon} spread uniformly over the off-support
     *                   topologies; must be in {@code [0,1)}
     */
    public UniformEscapeCCD(List<Tree> trees, double burnin, double escapeMass) {
        super(trees, burnin);
        this.escapeMass = checkEscapeMass(escapeMass);
        this.log1mEscapeMass = Math.log(1.0 - this.escapeMass);
    }

    /**
     * Constructor from a tree set.
     *
     * @param treeSet    an iterable set of trees with no burn-in
     * @param escapeMass total escape mass {@code epsilon} in {@code [0,1)}
     */
    public UniformEscapeCCD(TreeSet treeSet, double escapeMass) {
        super(treeSet, true);
        this.escapeMass = checkEscapeMass(escapeMass);
        this.log1mEscapeMass = Math.log(1.0 - this.escapeMass);
    }

    /**
     * Constructor for an empty CCD (trees added later or via {@link #copy()}).
     *
     * @param numLeaves      number of leaves of the trees this CCD will be based on
     * @param storeBaseTrees whether to store the trees used to build this CCD
     * @param escapeMass     total escape mass {@code epsilon} in {@code [0,1)}
     */
    public UniformEscapeCCD(int numLeaves, boolean storeBaseTrees, double escapeMass) {
        super(numLeaves, storeBaseTrees);
        this.escapeMass = checkEscapeMass(escapeMass);
        this.log1mEscapeMass = Math.log(1.0 - this.escapeMass);
    }

    private static double checkEscapeMass(double escapeMass) {
        if (!(escapeMass >= 0.0 && escapeMass < 1.0)) {
            throw new IllegalArgumentException("escapeMass must be in [0, 1), but was " + escapeMass);
        }
        return escapeMass;
    }

    /* -- ESCAPE BOOKKEEPING -- */

    /** @return the total escape mass {@code epsilon}. */
    public double getEscapeMass() {
        return escapeMass;
    }

    /** @return the number of rooted topologies outside the backbone's support, {@code (2n-3)!! - K}. */
    public BigInteger getNumberOfOffSupportTrees() {
        return numberOfRootedTopologies(getNumberOfLeaves()).subtract(super.getNumberOfTrees());
    }

    /**
     * @return the (natural) log-probability of any single off-support tree, {@code log(epsilon/M)};
     * {@link Double#NEGATIVE_INFINITY} if the backbone already covers every topology.
     */
    public double getEscapeLogProbability() {
        BigInteger k = super.getNumberOfTrees();
        if (cachedSupportSize == null || !cachedSupportSize.equals(k)) {
            cachedSupportSize = k;
            BigInteger offSupport = numberOfRootedTopologies(getNumberOfLeaves()).subtract(k);
            escapeLogProbCache = (escapeMass <= 0.0 || offSupport.signum() <= 0)
                    ? Double.NEGATIVE_INFINITY
                    : Math.log(escapeMass) - logBigInteger(offSupport);
        }
        return escapeLogProbCache;
    }

    /* -- PROBABILITY -- */

    @Override
    public double getProbabilityOfTree(Tree tree) {
        double p = super.getProbabilityOfTree(tree);
        if (p > 0) {
            return (1.0 - escapeMass) * p;
        }
        double escapeLogProb = getEscapeLogProbability();
        return Double.isFinite(escapeLogProb) ? Math.exp(escapeLogProb) : 0.0;
    }

    @Override
    public double getLogProbabilityOfTree(Tree tree) {
        double logP = super.getLogProbabilityOfTree(tree);
        if (Double.isFinite(logP)) {
            return log1mEscapeMass + logP;
        }
        return getEscapeLogProbability();
    }

    @Override
    public boolean containsTree(Tree tree) {
        // Full support: with positive escape mass every rooted topology has positive probability.
        // With no escape mass (or a saturated backbone) fall back to the backbone's support.
        return Double.isFinite(getEscapeLogProbability()) || super.containsTree(tree);
    }

    @Override
    public BigInteger getNumberOfTrees() {
        BigInteger total = numberOfRootedTopologies(getNumberOfLeaves());
        BigInteger k = super.getNumberOfTrees();
        return total.compareTo(k) > 0 ? total : k;
    }

    /* -- SAMPLING (log-prob mixture; see class doc) -- */

    @Override
    public double sampleTreeLogProbability() {
        double escapeLogProb = getEscapeLogProbability();
        if (Double.isFinite(escapeLogProb) && random.nextDouble() < escapeMass) {
            return escapeLogProb;
        }
        return log1mEscapeMass + super.sampleTreeLogProbability();
    }

    @Override
    public double sampleTreeProbability() {
        double escapeLogProb = getEscapeLogProbability();
        if (Double.isFinite(escapeLogProb) && random.nextDouble() < escapeMass) {
            return Math.exp(escapeLogProb);
        }
        return (1.0 - escapeMass) * super.sampleTreeProbability();
    }

    /* -- OTHER -- */

    @Override
    protected double getNumberOfParameters() {
        // the CCD1 backbone's parameters plus the single escape parameter epsilon
        return super.getNumberOfParameters() + 1;
    }

    @Override
    public AbstractCCD copy() {
        UniformEscapeCCD copy = new UniformEscapeCCD(this.getSizeOfLeavesArray(), false, this.escapeMass);
        copy.baseTrees.add(this.getSomeBaseTree());
        copy.numBaseTrees = this.getNumberOfBaseTrees();

        AbstractCCD.buildCopy(this, copy);

        return copy;
    }

    @Override
    public String toString() {
        return "UniformEscapeCCD [epsilon=" + escapeMass + "] " + super.toString();
    }

    /* -- STATIC HELPERS -- */

    /**
     * The number of distinct rooted bifurcating topologies on {@code n} labelled taxa, {@code (2n-3)!!}.
     * Returns {@code 1} for {@code n <= 2}.
     *
     * @param n number of taxa
     * @return {@code (2n-3)!!} as a {@link BigInteger}
     */
    public static BigInteger numberOfRootedTopologies(int n) {
        return AbstractCCD.numberOfRootedTopologies(n);
    }

}
