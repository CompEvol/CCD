package ccd.algorithms.credibleSets;

import ccd.model.AbstractCCD;
import ccd.model.CCD0;
import ccd.model.CCD1;
import ccd.model.SampleDistribution;

/**
 * Strategies for credible sets and credible intervals for {@link ccd.model.ITreeDistribution}s.
 * Also used by algorithms to pre-construct the credible set information.
 * Note all strategies make sense for every type of distributions, so not all are necessarily supported.
 */
public enum CredibleSetType {
    /**
     * Main strategy - Credible set and credible levels based on sampling trees from the distribution;
     * applicable for all distributions.
     */
    TreeSampling("Main strategy - Build the credible set information by sampling many trees."),
    /**
     * Credible set and credible levels based on frequencies of trees;
     * only applicable for {@link SampleDistribution}.
     */
    Frequency("Build the credible set information based on frequencies of trees (for SampleDistribution)"),
    /**
     * Credible set and credible levels based on credible CCD, build by greedy removal of clades based on their <em>probability</em> ;
     * applicable for all {@link AbstractCCD} but intended for {@link CCD0} models.
     */
    CladeProbability("Build the credible set by greedily taking out clades based on their probability."),
    /**
     * Credible set and credible levels based on credible CCD, build by greedy removal of clade partitions based on their <em>probability</em>;
     * applicable for all {@link AbstractCCD} but intended for {@link CCD1} models.
     */
    PartitionProbability("Build the credible set by greedily taking out clade partitions based on their probability (for CCD1 and CCD2)."),
    /**
     * (EXPERIMENTAL) Credible set and credible levels based on credible CCD, build by greedy removal of clades based on their <em>per tree ratio</em>.;
     * applicable for all {@link AbstractCCD} but intended for {@link CCD0} models.
     */
    CladeTrees("(Experimental) Build the credible set by greedily taking out clades based on their probability per tree ratio."),
    /**
     * (EXPERIMENTAL) Credible set and credible levels based on credible CCD, build by greedy removal of clades based on their <em>initial frequency</em>.;
     * applicable for all {@link AbstractCCD} but intended for {@link CCD0} models. Note that would equal probability for CCD1.
     */
    CladeCount("(Experimental) Build the credible set by greedily taking out clades based on their initial clade count.");

    /** Strategy description for UI tools. */
    private final String description;

    /**
     * Default constructor.
     *
     * @param description required information
     */
    CredibleSetType(String description) {
        this.description = description;
    }

    @Override
    public String toString() {
        return description;
    }
}
