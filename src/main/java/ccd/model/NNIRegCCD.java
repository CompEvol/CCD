package ccd.model;

import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator;
import ccd.algorithms.NNICladeExpansion.PairingMode;
import ccd.algorithms.regularisation.CCD1Regularisor;
import ccd.algorithms.regularisation.CCDExpansion;
import ccd.algorithms.regularisation.CCDRegularisationStrategy;
import ccd.algorithms.regularisation.NNICladeExpander;

import java.util.List;

/**
 * A {@link RegCCD} whose clade set is additionally expanded with the clades
 * reachable by a single nearest-neighbour-interchange (NNI) move (see
 * {@link NNICladeExpander}).
 *
 * <p>
 * Plain regCCD ({@link RegCCD}) assigns nonzero probability to trees with
 * unsampled <em>splits</em> of sampled clades, but a tree containing an
 * <em>unsampled clade</em> still gets zero probability. This class closes that
 * gap at the clade level: it proposes the NNI-reachable novel clades, wires them
 * into the graph with {@link CCDExpansion} exactly as for observed clades, and
 * prices every split (observed, split-expanded, and NNI) with the same
 * additive-{@code alpha} smoothing. No new free parameter is introduced --- the
 * novel clades carry no marginal parameter, and their probability is induced by
 * the CCD1 product of conditional clade probabilities.
 * </p>
 *
 * <p>
 * Pipeline:
 * <ol>
 * <li>build the CCD1 graph from the trees;</li>
 * <li>split expansion ({@link CCDExpansion}) over observed clades, so NNI sees
 * the full set of compatible observed-clade splits as recombination contexts;</li>
 * <li>NNI clade expansion ({@link PairingMode}) --- add novel clades and their
 * NNI-implied splits;</li>
 * <li>split expansion ({@link CCDExpansion}) over all clades, now including
 * compatible splits involving the novel clades;</li>
 * <li>additive-{@code alpha} regularisation;</li>
 * <li>finalise deterministic single-split novel clades.</li>
 * </ol>
 * </p>
 *
 * @author Claude (CCD-Sophie)
 */
public class NNIRegCCD extends RegCCD {

    /** Default regularisation pseudocount, matching {@link RegCCD}. */
    public static final double DEFAULT_ALPHA = 0.4;

    private final PairingMode pairingMode;
    private final double alpha;
    /** Second pseudocount for NNI splits; NaN means single-alpha smoothing. */
    private final double beta;
    private NNICladeExpander cladeExpander;

    /**
     * @param trees  the trees whose distribution is approximated
     * @param burnin fraction of given trees discarded as burn-in
     * @param mode   how splits are paired into NNI recombinations
     */
    public NNIRegCCD(List<Tree> trees, double burnin, PairingMode mode) {
        this(trees, burnin, mode, DEFAULT_ALPHA);
    }

    /**
     * @param trees  the trees whose distribution is approximated
     * @param burnin fraction of given trees discarded as burn-in
     * @param mode   how splits are paired into NNI recombinations
     * @param alpha  additive-smoothing pseudocount
     */
    public NNIRegCCD(List<Tree> trees, double burnin, PairingMode mode, double alpha) {
        this(trees, burnin, mode, alpha, Double.NaN);
    }

    /**
     * Two-level smoothing: NNI-derived splits get pseudocount {@code beta} while
     * observed and split-expanded splits keep {@code alpha}. With {@code beta <
     * alpha} the speculative NNI recombinations dilute observed splits less.
     *
     * @param trees  the trees whose distribution is approximated
     * @param burnin fraction of given trees discarded as burn-in
     * @param mode   how splits are paired into NNI recombinations
     * @param alpha  pseudocount for observed and split-expanded splits
     * @param beta   pseudocount for NNI-derived splits (NaN for single-alpha)
     */
    public NNIRegCCD(List<Tree> trees, double burnin, PairingMode mode, double alpha, double beta) {
        super(trees, burnin, false);
        this.pairingMode = mode;
        this.alpha = alpha;
        this.beta = beta;
        nniPostConstruction();
    }

    /**
     * @param treeSet trees without burn-in used to build the CCD
     * @param mode    how splits are paired into NNI recombinations
     * @param alpha   additive-smoothing pseudocount
     */
    public NNIRegCCD(TreeAnnotator.TreeSet treeSet, PairingMode mode, double alpha) {
        this(treeSet, mode, alpha, Double.NaN);
    }

    /**
     * @param treeSet trees without burn-in used to build the CCD
     * @param mode    how splits are paired into NNI recombinations
     * @param alpha   pseudocount for observed and split-expanded splits
     * @param beta    pseudocount for NNI-derived splits (NaN for single-alpha)
     */
    public NNIRegCCD(TreeAnnotator.TreeSet treeSet, PairingMode mode, double alpha, double beta) {
        super(treeSet, false);
        this.pairingMode = mode;
        this.alpha = alpha;
        this.beta = beta;
        nniPostConstruction();
    }

    private void nniPostConstruction() {
        // 1. split expansion among observed clades, so NNI has the full set of
        //    compatible observed-clade splits to recombine (GRAPH_WIDE picks up
        //    novel clades from CCD0-expanded contexts that CCD1 would miss; for
        //    CO_OCCURRING the enumeration walks base trees and is unaffected)
        new CCDExpansion().expandCCD(this);

        // 2. NNI clade expansion on the split-expanded observed graph
        cladeExpander = new NNICladeExpander(pairingMode);
        cladeExpander.expandClades(this);

        // 3. split expansion over all clades (observed + novel), as in plain regCCD
        new CCDExpansion().expandCCD(this);

        // 3. regularisation over the full graph: single-alpha, or two-level alpha/beta
        if (Double.isNaN(beta)) {
            regulariseRegCCD(alpha);
        } else {
            new CCD1Regularisor(CCDRegularisationStrategy.AdditiveXY, alpha, beta).regularise(this);
        }

        // 4. finalise deterministic single-split novel clades (0/0 otherwise)
        NNICladeExpander.fixZeroOccurrenceSingletonCCPs(this);

        System.out.println("NNI clade expansion [" + pairingMode + "]: added "
                + cladeExpander.getCladesAdded() + " clades; regularised with alpha = " + alpha
                + (Double.isNaN(beta) ? "" : ", beta = " + beta));
    }

    /** @return the pairing mode used for the NNI clade expansion */
    public PairingMode getPairingMode() {
        return pairingMode;
    }

    /** @return the number of novel clades added by the NNI expansion */
    public int getNumberOfNNIClades() {
        return (cladeExpander == null) ? 0 : cladeExpander.getCladesAdded();
    }

    @Override
    public String toString() {
        return "NNI-expanded " + super.toString();
    }
}
