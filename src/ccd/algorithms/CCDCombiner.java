package ccd.algorithms;

import ccd.model.AbstractCCD;
import ccd.model.BitSet;
import ccd.model.CCD0;
import ccd.model.Clade;
import ccd.model.CladePartition;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

/**
 * This class provides methods to combine CCDs into a single CCD.
 *
 * @author Jonathan Klawitter
 */
public class CCDCombiner {

    /**
     * Creates a new CCD that contains all clades and clade partitions of the given CCDs;
     * their probabilities and other attributes are not set to something specific.
     *
     * @param baseCCDs the CCDs to be combined
     * @return combined CCD with all clades and clade partitions of the given CCDs
     */
    public static AbstractCCD union(AbstractCCD[] baseCCDs) {
        AbstractCCD unionCCD = baseCCDs[0].copy();

        for (int i = 1; i < baseCCDs.length; i++) {
            AbstractCCD nextCCD = baseCCDs[i];
            expandCCDWithOther(unionCCD, nextCCD);
        }

        return unionCCD;
    }

    /**
     * Expands the first CCD with the clades and clade partitions of the second CCD.
     * No probabilities are set.
     *
     * @param expandingCCD the CCD that gets expanded
     * @param referenceCCD the CCD whose clades and clade partitions the other CCD is
     *                     expanded with
     */
    private static void expandCCDWithOther(AbstractCCD expandingCCD, AbstractCCD referenceCCD) {
        List<Clade> nextClades = referenceCCD.getClades().stream()
                .sorted(Comparator.comparingInt(x -> x.size())).toList();

        for (Clade referenceClade : nextClades) {
            if (referenceClade.isLeaf()) {
                continue;
            }

            // if a clade is not in the unionCCD yet, we add it
            Clade expandingClade = expandingCCD.getClade(referenceClade.getCladeInBits());
            if (expandingClade == null) {
                expandingClade = referenceClade.copy(expandingCCD);
                expandingCCD.getCladeMapping().put(expandingClade.getCladeInBits(), expandingClade);
            }

            // if a clade partition is not in the unionCCD yet, we add it
            partitionLoop:
            for (CladePartition referencePartition : referenceClade.getPartitions()) {
                for (CladePartition expandingPartition : expandingClade.getPartitions()) {
                    if (expandingPartition
                            .containsClade(referencePartition.getChildClades()[1].getCladeInBits())) {
                        // unionCCD already contains this partition
                        continue partitionLoop;
                    }
                }

                expandingClade.createCladePartition(expandingCCD.getClade(referencePartition.getChildClades()[0].getCladeInBits()),
                        expandingCCD.getClade(referencePartition.getChildClades()[1].getCladeInBits()));
            }
        }
    }

    /**
     * Combines the gives CCDs into a single CCD containing all their clades and clade partitions
     * and with clade partition probabilities set to the (weighted) average.
     *
     * @param baseCCDs the CCDs to be combined
     * @return combined CCD with clade partition probabilities set to the (weighted) average
     */
    public static AbstractCCD average(AbstractCCD[] baseCCDs) {
        // 1. create CCD that contains all clade and clade splits of all base CCDs
        AbstractCCD combinedCCD = union(baseCCDs);

        // 2. compute probability for each clade split
        // 2.1 make sure everything is initialized
        for (AbstractCCD ccd : baseCCDs) {
            ccd.computeCladeProbabilities();
        }

        // 2.2 loop over clades and then clade partitions
        int numCCDs = baseCCDs.length;
        for (Clade clade : combinedCCD.getClades()) {
            if (clade.isLeaf()) {
                continue;
            }

            BitSet cladeInBits = clade.getCladeInBits();
            Clade[] baseClades = new Clade[numCCDs];

            // compute clade weights per CCD
            double sumProbabilities = 0;
            for (int i = 0; i < numCCDs; i++) {
                baseClades[i] = baseCCDs[i].getClade(cladeInBits);
                if (baseClades[i] != null) {
                    sumProbabilities += baseClades[i].getProbability();
                }
            }
            double[] cladeWeights = new double[numCCDs];
            for (int i = 0; i < numCCDs; i++) {
                cladeWeights[i] = (baseClades[i] != null) ?
                        (baseClades[i].getProbability() / sumProbabilities) : 0;
            }

            // set probability for each clade partition
            for (CladePartition partition : clade.getPartitions()) {
                Clade childOne = partition.getChildClades()[0];
                Clade childTwo = partition.getChildClades()[1];
                double avgProb = 0;
                for (int i = 0; i < numCCDs; i++) {
                    if (baseClades[i] != null) {
                        // get corresponding partition in baseCCD, if it exists
                        Clade baseChildOne = baseCCDs[i].getClade(childOne.getCladeInBits());
                        Clade baseChildTwo = baseCCDs[i].getClade(childTwo.getCladeInBits());
                        CladePartition basePartition = baseClades[i].getCladePartition(baseChildOne, baseChildTwo);

                        if (basePartition != null) {
                            avgProb += cladeWeights[i] * basePartition.getCCP();
                        }
                    }
                }

                partition.setCCP(avgProb);
            }
        }

        return combinedCCD;
    }

    /**
     * Combines the gives CCD0s into a single CCD0 containing all their clades and clade partitions
     * but not more (the expand step is not executed); clade partition probabilities are set based
     * as normal in a CCD0. Asserts that the base CCD0s nearly have the same number of base trees.
     *
     * @param baseCCDs the CCDs to be combined
     * @return combined CCD with clade partition probabilities set to the (weighted) average
     */
    public static CCD0 combineCCD0Unexpanded(AbstractCCD[] baseCCDs) {
        // 1. input assertions
        // 1.1 check/cast base CCDs to CCD0
        CCD0[] baseCCD0s = new CCD0[baseCCDs.length];
        for (int i = 0; i < baseCCDs.length; i++) {
            baseCCD0s[i] = (CCD0) baseCCDs[i];
        }
        // 1.2 check that base CCDs have roughly same amount of trees
        int numberOfFirstBaseTrees = baseCCD0s[0].getNumberOfBaseTrees();
        int totalNumberOfBaseTrees = numberOfFirstBaseTrees;
        int maxNumTreesDifference = 2;
        for (int i = 1; i < baseCCD0s.length; i++) {
            int numberOfBaseTrees = baseCCD0s[i].getNumberOfBaseTrees();
            totalNumberOfBaseTrees += numberOfBaseTrees;
            if (Math.abs(numberOfBaseTrees - numberOfFirstBaseTrees) > maxNumTreesDifference) {
                throw new AssertionError("Combining CCD0s with too different amount of trees.");
            }
        }

        // 2. create CCD that contains all clade and clade splits of all base CCDs
        CCD0 combinedCCD = (CCD0) union(baseCCD0s);

        // 2. set all counters as sum of the base CCDs
        // 2.1 number of base trees
        combinedCCD.setNumBaseTrees(totalNumberOfBaseTrees);
        // 2.2 number of clade occurrences
        int numCCDs = baseCCD0s.length;
        for (Clade clade : combinedCCD.getClades()) {
            BitSet cladeInBits = clade.getCladeInBits();
            int totalOccurrences = 0;
            for (CCD0 ccd0 : baseCCD0s) {
                Clade baseClade = ccd0.getClade(cladeInBits);
                if (baseClade != null) {
                    totalOccurrences += baseClade.getNumberOfOccurrences();
                }
            }
            clade.setNumberOfOccurrences(totalOccurrences);
        }

        combinedCCD.setPartitionProbabilities(combinedCCD.getRootClade());

        System.out.println("num base trees = " + combinedCCD.getNumberOfBaseTrees());
        System.out.println("combinedCCD.getRootClade().getNumberOfOccurrences() = " + combinedCCD.getRootClade().getNumberOfOccurrences());

        return combinedCCD;
    }


    /**
     * Combines the gives CCDs into a single CCD with the probability of clade partitions
     * set to the product of the respective partitions.
     * Hence, if a clade partition (or even clade) does not appear in one of the base CCDs,
     * then the probability is zero.
     * The returned CCD is tidied up and so does not contain any such clade partitions
     * or clades without partitions (besides leaves and the root clade); cf. {@link AbstractCCD#tidyUpCCDGraph()}.
     *
     * <p>
     * <i>Combining:</i> Roughly speaking, for k base CCDs, the combined CCP of a
     * clade partition is the product of its clade partition in the base CCDs times
     * divided by k-1 its prior probability times an unknown normalization constant.
     * The prior probability has to be taken out k-1 times, because otherwise the
     * prior would contribute too many times in the CCP of the combined inference.
     * Note that the prior probability depends on the model used in the Bayesian
     * MCMC algorithm. <br>
     * <i>Remark:</i> Currently only the prior probabilities for clade partitions
     * under standard birth-death and coalescent models with contemporaneous leaves
     * is supported.
     * </p>
     *
     * @param baseCCDs the CCDs to be combined
     * @return combined CCD with clade partition probabilities set to the product
     */
    public static AbstractCCD product(AbstractCCD[] baseCCDs) {
        /* Caching table for prior probabilities of clade partitions
        with respect to parent clade size and smaller clade size. */
        double[][] priorProbabilityTable = new double[baseCCDs[0].getNumberOfLeaves()][];

        // combine all CCDs into one and multiply probabilities
        AbstractCCD combinedCCD = baseCCDs[0].copy();
        for (int i = 1; i < baseCCDs.length; i++) {
            combine(combinedCCD, baseCCDs[i], priorProbabilityTable);
        }
        // ... and then normalize them
        normalize(combinedCCD);

        return combinedCCD;
    }

    /**
     * Extends {@link CCDCombiner#product(AbstractCCD[])} by also using spiking to combine the given CCDs
     * with a product of clade partition probabilities.
     *
     * <p>
     * <i>Spiking:</i> Under the basic combining method, when a clade or clade
     * partition does not appear in all base CCDs, then its probability is zero in
     * the combined CCD. This can have cascading effects as then any clade partition
     * containing this clade has to be removed as well. In the worst case, all
     * clades get removed. This happens exactly when the intersection of the trees
     * represented by the base CCDs is empty. To address this problem, when
     * <i>spiking</i> is used, all base CCDs are <i>expanded</i> to include all
     * clades and clade partitions in the union of the base CCDs. The CCPs of the
     * clade partitions of a particular clade are then spiked as if this clade was
     * observed one more time with each clade partition observed proportional to its
     * prior probability.
     * </p>
     *
     * @param baseCCDs the CCDs to be combined
     * @return combined CCD with clade partition probabilities
     * set to the product after using expansion and spiking
     */
    public static AbstractCCD spikedProduct(AbstractCCD[] baseCCDs) {
        /* Caching table for prior probabilities of clade partitions
        with respect to parent clade size and smaller clade size. */
        double[][] priorProbabilityTable = new double[baseCCDs[0].getNumberOfLeaves()][];

        // for spiking, we use copies of the base CCDs
        System.out.println("- initialize");
        AbstractCCD[] copiedCCDs = new AbstractCCD[baseCCDs.length];
        for (int i = 1; i < baseCCDs.length; i++) {
            copiedCCDs[i] = baseCCDs[i].copy();
        }
        baseCCDs = copiedCCDs;

        System.out.println("- expand CCDs");
        // as combinedCCD we use the union of all base CCDs
        baseCCDs[0] = union(baseCCDs);
        AbstractCCD combinedCCD = baseCCDs[0];

        // we then use the combinedCCD to expand all other base CCDs
        for (int i = 1; i < baseCCDs.length; i++) {
            expandCCDWithOther(baseCCDs[i], combinedCCD);
        }
        System.out.println(" expanded, num clades: " + combinedCCD.getNumberOfClades());

        // we can then spike all CCDs
        spike(baseCCDs, priorProbabilityTable);

        for (int i = 1; i < baseCCDs.length; i++) {
            combine(combinedCCD, baseCCDs[i], priorProbabilityTable);
        }

        normalize(combinedCCD);

        return combinedCCD;
    }

    /* Helper method - integrates nextCCD into combinedCCD. */
    private static void combine(AbstractCCD combinedCCD, AbstractCCD nextCCD, double[][] priorProbabilityTable) {
        // for every clade, we go through all clade partitions and multiply their CCPs;
        // if no expanded, spiked CCDs are used, then it might happen
        // that clade or clade partition does not exist in both CCDs
        for (Clade clade : combinedCCD.getClades()) {
            if (clade.isLeaf()) {
                continue;
            }

            ArrayList<CladePartition> partitionsToRemove = new ArrayList<>();
            Clade otherClade = nextCCD.getClade(clade.getCladeInBits());
            if (otherClade == null) {
                // the other CCD does not have this clade,
                // so we remove it from this one
                // (indirectly by dropping all its partitions);
                // note that this does not happen with expanded, spiked CCDs
                partitionsToRemove.addAll(clade.getPartitions());
            } else {
                ArrayList<CladePartition> otherPartitions = otherClade.getPartitions();
                partitionLoop:
                for (CladePartition partition : clade.getPartitions()) {
                    for (CladePartition otherPartition : otherPartitions) {
                        if (otherPartition
                                .containsClade(partition.getChildClades()[1].getCladeInBits())) {
                            // so otherClade has the same partition as the
                            // current clade and we can thus update their probabilities
                            double prior = getPartitionPriorProbability(partition, priorProbabilityTable);
                            double unnormalizedCCP = partition.getCCP() * otherPartition.getCCP()
                                    / prior;
                            partition.setCCP(unnormalizedCCP);
                            continue partitionLoop;
                        }
                    }

                    // if this point reached, then otherClades does not have
                    // this partition, so we have to remove it;
                    // note that this does not happen with expanded, spiked CCDs
                    partitionsToRemove.add(partition);
                }

            }

            // remove partitions only occurring in the combined CCD
            for (CladePartition partition : partitionsToRemove) {
                Clade firstChild = partition.getChildClades()[0];
                Clade secondChild = partition.getChildClades()[1];

                clade.getChildClades().remove(firstChild);
                clade.getChildClades().remove(secondChild);

                firstChild.getParentClades().remove(clade);
                secondChild.getParentClades().remove(clade);
                clade.getPartitions().remove(partition);
            }
        } // handling clades

        if (!combinedCCD.tidyUpCCDGraph()) {
            System.err.println("After combining CCDs, we got an empty tree distribution.");
        }
    }

    /**
     * Helper method - computes the prior probability of the given partition under birth-death
     * or coalescent models with <i>contemporaneous</i> leaves. Uses caching for efficiency.
     *
     * @param partition             whose prior probability is computed
     * @param priorProbabilityTable caching table
     * @return the prior probability of the given partition under birth-death or
     * coalescent models with contemporaneous leaves
     */
    private static double getPartitionPriorProbability(CladePartition partition, double[][] priorProbabilityTable) {
        int nParent = partition.getParentClade().size();
        int nChild = partition.getChildClades()[0].size();
        if (nChild > (nParent / 2)) {
            nChild = nParent - nChild;
        }

        // use caching
        if (priorProbabilityTable[nParent - 1] == null) {
            priorProbabilityTable[nParent - 1] = new double[nParent - 1];
        }
        double prob = priorProbabilityTable[nParent - 1][nChild - 1];

        // if not cached, compute value
        if (prob == 0) {
            double binom = 1;
            for (int i = 0; i < nChild; i++) {
                binom *= (nParent - i) / (double) (nChild - i);
            }
            prob = (2 / (double) (nParent - 1)) * (1.0 / binom);
            priorProbabilityTable[nParent - 1][nChild - 1] = prob;
        }

        return prob;
    }

    /**
     * Helper method to normalize the clade partition probabilities of the given CCD.
     *
     * @param ccd whose clade partition probabilities get e normalized
     */
    private static void normalize(AbstractCCD ccd) {
        for (Clade clade : ccd.getClades()) {
            if (clade.isLeaf()) {
                continue;
            }

            // normalize the CCPs of the partitions
            double sumProbabilities = clade.getPartitions().stream().mapToDouble(x -> x.getCCP())
                    .sum();
            clade.getPartitions().forEach(x -> x.setCCP(x.getCCP() / sumProbabilities));
        }
    }

    /* Helper method - spikes the expanded CCDs. */
    private static void spike(AbstractCCD[] ccds, double[][] priorProbabilityTable) {
        System.out.println("- spike CCDs");
        for (AbstractCCD ccd : ccds) {
            // for each clade, spike the CCPs of its partitions as if we had one
            // more observation distributed among the observed partitions
            // proportional to their prior probability
            for (Clade clade : ccd.getClades()) {
                double sumPriors = 0;
                for (CladePartition partition : clade.getPartitions()) {
                    sumPriors += getPartitionPriorProbability(partition, priorProbabilityTable);
                }
                double normalizer = 1 / sumPriors;
                // note that normalizer is the same for each CCD, so could be
                // cached; since the values it is based on get cached, probably not
                // worth it

                for (CladePartition partition : clade.getPartitions()) {
                    double newCCP = (partition.getNumberOfOccurrences()
                            + normalizer * getPartitionPriorProbability(partition, priorProbabilityTable))
                            / (partition.getParentClade().getNumberOfOccurrences() + 1);
                    partition.setCCP(newCCP);
                }
            }
        }
    }
}
