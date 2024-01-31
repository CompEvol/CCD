package ccd.model;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator.TreeSet;

/**
 * This class represents a tree distribution via the conditional clade
 * distribution (CCD) as defined by
 * <a href="https://doi.org/10.1093/sysbio/syt014">B. Larget, 2013</a>. Clade
 * partition probabilities (CCP) are set based on the frequency f(C1, C2) of the
 * clade partition {C1, C2} divided by the frequency f(C) of the parent clade C,
 * i.e. f(C1,C2)/f(C). The CCD graph and computations are then handled by the
 * parent class {@link AbstractCCD}.
 * 
 * <p>
 * It can be constructed in one go by providing a set of trees or maintained
 * (e.g. during an MCMC run) by adding and removing trees.
 * </p>
 * 
 * <p>
 * It can be used to sample trees with different {@link SamplingStrategy} and
 * {@link HeightSettingStrategy}, compute their probabilities, compute the
 * entropy of the distribution, ...
 * </p>
 * 
 * <p>
 * The MAP tree of this distribution is the tree with highest CCP.
 * </p>
 * 
 * @author Jonathan Klawitter
 */
public class CCD1 extends AbstractCCD {

	/* -- CONSTRUCTORS & CONSTRUCTION METHODS -- */
	/**
	 * Constructor for a {@link CCD1} based on the given collection of trees
	 * with specified burn-in.
	 * 
	 * @param trees
	 *            the trees whose distribution is approximated by the resulting
	 *            {@link CCD1}
	 * @param burnin
	 *            value between 0 and 1 of what percentage of the given trees
	 *            should be discarded as as burn-in
	 */
	public CCD1(List<Tree> trees, double burnin) {
		super(trees, burnin);
	}

	/**
	 * Constructor for a {@link CCD1} based on the given collection of trees
	 * (not containing any burnin trees).
	 * 
	 * @param treeSet
	 *            an iterable set of trees, which contains no burnin trees,
	 *            whose distribution is approximated by the resulting
	 *            {@link CCD1}
	 * @param storeBaseTrees
	 *            whether to store the trees used to create this CCD
	 */
	public CCD1(TreeSet treeSet, boolean storeBaseTrees) {
		super(treeSet, storeBaseTrees);
	}

	/**
	 * Constructor for an empty CDD. Trees can then be processed one by one.
	 * 
	 * @param numLeaves
	 *            number of leaves of the trees that this CCD will be based on
	 * @param storeBaseTrees
	 *            whether to store the trees used to create this CCD;
	 *            recommended not to when huge set of trees is used
	 */
	public CCD1(int numLeaves, boolean storeBaseTrees) {
		super(numLeaves, storeBaseTrees);
	}


	/* -- GENERAL GETTERS & SETTERS -- */
	/**
	 * @return number of leaves/taxa remaining in the trees represented by this
	 *         CCD, which only differs from {@link CCD1#getNumberOfLeaves()} if
	 *         this is a {@link FilteredConditionalCladeDistribution}
	 */
	public int getNumberOfActualLeaves() {
		return this.getNumberOfLeaves();
	}

	/**
	 * @return the taxa of this CCD as BitSet with length based on original
	 *         number of leaves (only different if this a
	 *         {@link FilteredConditionalCladeDistribution})
	 */
	@Override
	public BitSet getTaxaAsBitSet() {
		return this.rootClade.getCladeInBits();
	}


	/* -- STATE MANAGEMENT - STATE MANAGEMENT -- */
	@Override
	public void setAsDirty() {
		super.setAsDirty();
	}

	@Override
	protected void tidyUpIfDirty() {
		resetCachedValuesIfProbabilitesDirty();
	}


	/* -- PROBABILITY, POINT ESTIMATE & SAMPLING METHODS -- */
	// all handled by parent class AbstractCCD as long as clade partition
	// probabilities are set correctly


	/* -- OTHER METHODS -- */
	public List<CladePartition> findAttachmentPointsOfClade(Clade attachingClade) {
		// TODO WORK IN PRORGESS
		ArrayList<CladePartition> parentPartitions = new ArrayList<CladePartition>();
		for (Clade clade : cladeMapping.values()) {
			if (clade == attachingClade) {
				continue;
			} else if (clade.containsClade(attachingClade)) {
				for (CladePartition partition : clade.getPartitions()) {
					if (partition.containsChildClade(attachingClade)) {
						parentPartitions.add(partition);
					}
				}
			}
		}

		return parentPartitions;
	}

	@Override
	public String toString() {
		return "CCD0 " + super.toString();
	}

	/**
	 * Create a (deep) copy of this CCD, so with copies of the Clades and
	 * CladePartitions; copies at most one stored tree.
	 * 
	 * @return a (deep) copy of this CCD
	 */
	public CCD1 copy() {
		CCD1 copy = new CCD1(this.getNumberOfLeaves(), false);
		copy.baseTrees.add(this.getSomeBaseTree());

		for (Clade originalClade : this.getClades()) {
			Clade copiedClade = originalClade.copy(copy);
			copy.cladeMapping.put(originalClade.getCladeInBits(), copiedClade);
		}
		copy.rootClade = cladeMapping.get(this.getRootClade().getCladeInBits());

		for (Clade originalClade : this.getClades()) {
			for (CladePartition originalPartition : originalClade.getPartitions()) {
				Clade copiedParent = copy.cladeMapping.get(originalClade.getCladeInBits());
				Clade copiedChildFirst = copy.cladeMapping
						.get(originalPartition.getChildClades()[0].getCladeInBits());
				Clade copiedChildSecond = copy.cladeMapping
						.get(originalPartition.getChildClades()[1].getCladeInBits());

				CladePartition copiedPartition = copiedParent.createCladePartition(copiedChildFirst,
						copiedChildSecond, true);
				if (originalPartition.getNumberOfOccurrences() <= 0) {
					copiedPartition.setCCP(originalPartition.getCCP());
				} else {
					copiedPartition.increaseOccurrenceCountBy(
							originalPartition.getNumberOfOccurrences(),
							originalPartition.getMeanOccurredHeight());
				}
			}
		}

		return copy;
	}

}
