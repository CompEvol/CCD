package ccd.algorithms.regularisation;

import ccd.algorithms.NNICladeExpansion;
import ccd.algorithms.NNICladeExpansion.NovelClade;
import ccd.algorithms.NNICladeExpansion.PairingMode;
import ccd.algorithms.NNICladeExpansion.Provenance;
import ccd.model.AbstractCCD;
import ccd.model.Clade;
import ccd.model.CladePartition;

/**
 * Applies the clade-level NNI expansion to a CCD graph by inserting the novel
 * clades it generates, together with the splits NNI implies.
 *
 * <p>
 * This is the clade-level analogue of {@link CCDExpansion}: where
 * {@code CCDExpansion} adds unobserved <em>splits</em> among already-present
 * clades, this adds unobserved <em>clades</em> (and the splits the NNI move
 * implies) so that a CCD1/regCCD can assign nonzero probability to trees that
 * contain a clade no sampled tree had. The enumeration is done read-only by
 * {@link NNICladeExpansion}; here we mutate the graph.
 * </p>
 *
 * <p>
 * For each novel clade {@code K = kept ∪ sibling} with provenance parent
 * {@code P} and dropped grandchild {@code B}, we add the seed split
 * {@code K -> {kept, sibling}} and the parent split {@code P -> {K, B}}. No
 * conditional probabilities are assigned here; a subsequent
 * {@link CCD1Regularisor} prices the new (zero-count) splits with the same
 * additive-{@code alpha} smoothing as the observed and {@code CCDExpansion}
 * splits. A novel clade carries no marginal parameter --- under the CCD1
 * parameterisation its probability is induced by the product of conditional
 * clade probabilities along the path to it.
 * </p>
 *
 * @author Claude (CCD-Sophie)
 */
public class NNICladeExpander {

    private final PairingMode mode;

    private int cladesAdded = 0;
    private int seedSplitsAdded = 0;
    private int parentSplitsAdded = 0;

    /**
     * @param mode how splits are paired into NNI recombinations
     *             ({@link PairingMode#GRAPH_WIDE} or
     *             {@link PairingMode#CO_OCCURRING})
     */
    public NNICladeExpander(PairingMode mode) {
        this.mode = mode;
    }

    /**
     * Enumerates the NNI-reachable novel clades of the given CCD and inserts
     * them, with their NNI-implied splits, into the CCD graph.
     *
     * @param ccd the CCD to expand in place (its observed splits are used for
     *            the enumeration, so call before any split-level expansion)
     */
    public void expandClades(AbstractCCD ccd) {
        NNICladeExpansion enumeration = new NNICladeExpansion(ccd, mode).compute();

        for (NovelClade novel : enumeration.getNovelClades()) {
            Clade k = ccd.getOrCreateClade(novel.taxa);
            k.setNNIExpanded(true);
            cladesAdded++;

            for (Provenance prov : novel.provenances) {
                // NNI seed split of the novel clade: K -> {kept, sibling}
                if (k.getCladePartition(prov.kept(), prov.sibling()) == null) {
                    k.createCladePartition(prov.kept(), prov.sibling());
                    seedSplitsAdded++;
                }
                // recombined parent split: P -> {K, dropped}
                Clade parent = prov.parent();
                if (parent.getCladePartition(k, prov.dropped()) == null) {
                    parent.createCladePartition(k, prov.dropped());
                    parentSplitsAdded++;
                }
            }
        }

        ccd.setCacheAsDirty();
    }

    /** @return number of novel clades inserted */
    public int getCladesAdded() {
        return cladesAdded;
    }

    /** @return number of NNI seed splits inserted */
    public int getSeedSplitsAdded() {
        return seedSplitsAdded;
    }

    /** @return number of recombined parent splits inserted */
    public int getParentSplitsAdded() {
        return parentSplitsAdded;
    }

    /**
     * Sets the conditional clade probability to 1 for every clade that has a
     * single split and no observed occurrences. Such clades (e.g. novel cherries
     * or novel clades with only their NNI seed split) are skipped by
     * {@link CCD1Regularisor}, and their count-based CCP would be 0/0; this makes
     * them deterministic, as there is only one way to resolve them.
     *
     * <p>
     * Call after regularisation. Existing observed clades with a single split
     * have positive occurrences, so their CCP is well defined and they are left
     * untouched.
     * </p>
     *
     * @param ccd the CCD whose novel single-split clades are finalised
     */
    public static void fixZeroOccurrenceSingletonCCPs(AbstractCCD ccd) {
        for (Clade clade : ccd.getClades()) {
            if (clade.getNumberOfPartitions() == 1 && clade.getNumberOfOccurrences() == 0) {
                CladePartition only = clade.getPartitions().get(0);
                if (!only.isCCPSet()) {
                    only.setCCP(1.0);
                }
            }
        }
    }
}
