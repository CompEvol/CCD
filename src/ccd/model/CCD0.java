package ccd.model;

import java.util.ArrayList;
import java.util.List;

import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator.TreeSet;

/**
 * This class represents a tree distribution using the CCD graph (via the parent
 * class {@link AbstractCCD}). Tree probabilities are set proportional to the
 * product of Monte Carlo probabilities of the samples used to construct this
 * CCD0 but normalized to have a tree distribution.
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
 * The MAP tree of this distribution is the tree with highest product of clade
 * credibility. Furthermore, with {@link CCD0#getMSCCTree()} the tree with
 * highest sum of clade credibility can be retrieved.
 * </p>
 * 
 * @author Jonathan Klawitter
 */
public class CCD0 extends AbstractCCD {

	/* -- CONSTRUCTORS & CONSTRUCTION METHODS -- */
	/**
	 * Constructor for a {@link CCD0} based on the given collection of trees
	 * with specified burn-in.
	 * 
	 * @param trees
	 *            the trees whose distribution is approximated by the resulting
	 *            {@link CCD0}
	 * @param burnin
	 *            value between 0 and 1 of what percentage of the given trees
	 *            should be discarded as as burn-in
	 */
	public CCD0(List<Tree> trees, double burnin) {
		super(trees, burnin);
		initialize();
	}

	/**
	 * Constructor for a {@link CCD0} based on the given collection of trees
	 * (not containing any burnin trees).
	 * 
	 * @param treeSet
	 *            an iterable set of trees, which contains no burnin trees,
	 *            whose distribution is approximated by the resulting
	 *            {@link CCD0}
	 * @param storeBaseTrees
	 *            whether to store the trees used to create this CCD
	 */
	public CCD0(TreeSet treeSet, boolean storeTrees) {
		super(treeSet, storeTrees);
		initialize();
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
	public CCD0(int numLeaves, boolean storeTrees) {
		super(numLeaves, storeTrees);
	}

	/* -- STATE MANAGEMENT & INITIALIZATION -- */
	/** Whether the CCD graph is expanded and CCPs renormalized. */
	private boolean dirtyStructure = true;

	@Override
	public void setAsDirty() {
		super.setAsDirty();
		this.dirtyStructure = true;
	}

	@Override
	protected void tidyUpIfDirty() {
		if (dirtyStructure) {
			super.resetCachedValues();
			this.initialize();
		}
	}

	/**
	 * Set up this CCD0 after adding/removing trees. Assumes reset/not-set
	 * values in CCD graph.
	 */
	public void initialize() {
		// need to find all clade partitions that could exist but were not
		// observed in base trees
		// System.out.print("\nExpanding CCD graph ... ");
		expand();

		// then need to set clade partition probabilities
		// which normalizes the product of clade probabilities
		// System.out.print("setting probabilities ... ");
		setPartitionProbabilities(this.rootClade);

		// System.out.println("done.\n");
		this.dirtyStructure = false;
	}

	/*
	 * Expand CCD graph with clade partitions where parent and children were
	 * observed, but not that clade partition.
	 */
	private void expand() {
		int n = this.getNumberOfLeaves();
		int[] expandedClades = new int[n];

		// 1. clade buckets
		// for easier matching of child clades, we want to group them by size
		// 1.i sort clades
		Clade[] clades = cladeMapping.values().stream().sorted((x, y) -> Integer
				.compare(x.getCladeInBits().cardinality(), y.getCladeInBits().cardinality()))
				.toArray(Clade[]::new);

		/*- System.out.println("num clades: " + clades.length);
		 System.out.println("num partitions start: " + cladeMapping.values()
		 .stream().mapToInt(x -> x.partitions.size()).sum());*/

		// 1.ii init clade buckets
		List<List<Clade>> cladeBuckets = new ArrayList<List<Clade>>(n);
		for (int i = 0; i < n; i++) {
			if (i < 3) {
				cladeBuckets.add(new ArrayList<Clade>(n));
			} else {
				cladeBuckets.add(new ArrayList<Clade>());
			}
		}
		// 1.iii fill clade buckets
		for (int i = 0; i < clades.length; i++) {
			cladeBuckets.get(clades[i].size() - 1).add(clades[i]);
		}
		// int x = 0;
		// for (List<Clade> list : cladeBuckets) {
		// System.out.println("clade bucket " + (x++ + 1) + " has size: " +
		// list.size());
		// for (Clade clade : list) {
		// System.out.print(clade.getCladeInBits() + ", ");
		// }
		// System.out.println("");
		// }

		// 2. find missing clade partitions
		// we go through clades in increasing size, then check for each clade of
		// at most half potential parent's size, whether we can find a partner
		for (int i = 0; i < clades.length; i++) {
			Clade parent = clades[i];

			// we skip leaves and cherries as they have no/only one partition
			if (parent.isLeaf() || parent.isCherry()) {
				continue;
			}

			// otherwise we check if we find a larger partner clade for any
			// smaller clade that together partition the parent clade
			for (int j = 1; j <= parent.size() / 2; j++) {
				for (Clade smallChild : cladeBuckets.get(j - 1)) {
					if (parent.contains(smallChild.getCladeInBits())
							&& !smallChild.parentClades.contains(parent)) {
						// the small child is a subclade of the parent clade
						// but is not yet in a partition of the parent clade;
						// hence we check all potential partner clades
						for (Clade largeChild : cladeBuckets.get(parent.size() - j - 1)) {
							if (parent.contains(largeChild.getCladeInBits())
									&& !largeChild.intersects(smallChild)) {
								// ... whether they are subclades of the parent,
								// but the children's clade do not intersect
								parent.createCladePartition(smallChild, largeChild);
								expandedClades[parent.size() - 1]++;
								// System.err.println("+");
							}
						}
					}
				}
			}
		}

		/*-double[] avgExpandedClades = new double[n];
		int i = 0;
		System.out.println("k, #c, avg1, #exp, avgNew, ratio");
		DecimalFormatSymbols symbols = new DecimalFormatSymbols(Locale.US);
		DecimalFormat df = new DecimalFormat("0.##", symbols);
		for (List<Clade> list : cladeBuckets) {
			avgExpandedClades[i] = (list.size() != 0) ? expandedClades[i] / (double) list.size()
					: 0.0;
			int sum = list.stream().map(x -> x.getPartitions().size()).mapToInt(Integer::intValue)
					.sum();
			double avg = (list.size() != 0) ? sum / (double) list.size() : 0;
		
			String next = (i + 1) + ", ";
			next += list.size() + ", ";
			next += df.format(avg) + ", ";
			next += expandedClades[i] + ", ";
			next += df.format(avgExpandedClades[i]) + ", ";
			next += (avg != 0) ? df.format(avgExpandedClades[i] / avg) : 0;
			System.out.println(next);
			i++;
		}
		
		System.out.println("num partitions end: "
				+ cladeMapping.values().stream().mapToInt(x -> x.partitions.size()).sum());*/
	}

	/* Recursive helper method */
	private double setPartitionProbabilities(Clade clade) {
		if (clade.getSumCladeCredibilities() > 0) {
			return clade.getSumCladeCredibilities();
		}

		if (clade.isLeaf()) {
			// a leaf has no partition, sum of probabilities is 1
			return 1.0;
		} else if (clade.isCherry()) {
			// a cherry has only one partition
			clade.partitions.get(0).setCCP(1);
			clade.setSumCladeCredibilities(1);
			return clade.getCladeCredibility();
		} else {
			// other might have more partitions
			double sumSubtreeProbabilities = 0.0;
			double[] sumPartitionSubtreeProbabilities = new double[clade.getPartitions().size()];

			// compute sum of probabilities over all partitions ...
			int i = 0;
			for (CladePartition partition : clade.getPartitions()) {
				double sumPartitionSubtreeCladeCredibility = setPartitionProbabilities(
						partition.getChildClades()[0])
						* setPartitionProbabilities(partition.getChildClades()[1]);
				sumPartitionSubtreeProbabilities[i] = sumPartitionSubtreeCladeCredibility;
				sumSubtreeProbabilities += sumPartitionSubtreeProbabilities[i];
				i++;
			}
			// ... and then normalize
			i = 0;
			for (CladePartition partition : clade.getPartitions()) {
				double probability = sumPartitionSubtreeProbabilities[i]
						/ sumSubtreeProbabilities;
				partition.setCCP(probability);
				i++;
			}

			// combined with probability of clade, we get sum of all subtree
			// probabilities
			double sumCladeCredibilities = sumSubtreeProbabilities * clade.getCladeCredibility();
			clade.setSumCladeCredibilities(sumCladeCredibilities);
			return sumCladeCredibilities;
		}

	}

	/* -- POINT ESTIMATE METHODS -- */
	/**
	 * Returns the tree with maximum sum of clade probabilities (#occurrences /
	 * #num trees) of all trees in this CCD. Note that the maximum is with
	 * respect to the sum of clade credibilities, unlike the MAP tree or an MCC
	 * tree, which uses the product.
	 * 
	 * @return the tree with maximum sum of clade probabilities
	 */
	public Tree getMSCCTree() {
		return this.getMSCCTree(HeightSettingStrategy.None);
	}

	/**
	 * Returns the tree with maximum sum of clade probabilities (#occurrences /
	 * #num trees) of all trees in this CCD and heights set with the given
	 * strategy. Note that the maximum is with respect to the sum of clade
	 * credibilities, unlike the MAP tree or an MCC tree, which uses the
	 * product. Returns the tree with maximum sum of clade credibilities
	 * (#occurrences / #num trees) of all trees in this CCD and heights set with
	 * the given strategy. Note that the maximum is with respect to the sum of
	 * clade credibilities, unlike the MCC tree, which uses the product.
	 * 
	 * @param heightStrategy
	 *            the strategy used to set the heights of the tree vertices
	 * @return the tree with maximum sum of clade probabilities
	 */
	public Tree getMSCCTree(HeightSettingStrategy heightStrategy) {
		return getTreeBasedOnStrategy(SamplingStrategy.MaxSumCladeCredibility, heightStrategy);
	}

	/* -- OTHER METHODS -- */
	@Override
	public String toString() {
		return "CCD1 " + super.toString();
	}

}
