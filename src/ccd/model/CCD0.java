package ccd.model;

import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator.TreeSet;

import java.io.PrintStream;
import java.util.*;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.RejectedExecutionException;

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
     * @param trees  the trees whose distribution is approximated by the resulting
     *               {@link CCD0}
     * @param burnin value between 0 and 1 of what percentage of the given trees
     *               should be discarded as burn-in
     */
    public CCD0(List<Tree> trees, double burnin) {
        super(trees, burnin);
        initialize();
    }

    /**
     * Constructor for a {@link CCD0} based on the given collection of trees
     * (not containing any burnin trees).
     *
     * @param treeSet        an iterable set of trees, which contains no burnin trees,
     *                       whose distribution is approximated by the resulting
     *                       {@link CCD0}
     * @param storeBaseTrees whether to store the trees used to create this CCD
     */
    public CCD0(TreeSet treeSet, boolean storeBaseTrees) {
        super(treeSet, storeBaseTrees);
//        System.out.println("xx");
        initialize();
    }

    /**
     * Constructor for an empty CDD. Trees can then be processed one by one.
     *
     * @param numLeaves      number of leaves of the trees that this CCD will be based on
     * @param storeBaseTrees whether to store the trees used to create this CCD;
     *                       recommended not to when huge set of trees is used
     */
    public CCD0(int numLeaves, boolean storeBaseTrees) {
        super(numLeaves, storeBaseTrees);
    }

    /* -- STATE MANAGEMENT & INITIALIZATION -- */
    /**
     * Whether the CCD graph is expanded and CCPs renormalized.
     */
    private boolean dirtyStructure = true;
    private PrintStream progressStream;
    
    public void setProgressStream(PrintStream progressStream) {
    	this.progressStream = progressStream;
    }
    
    @Override
    public void setCacheAsDirty() {
        super.setCacheAsDirty();
        this.dirtyStructure = true;
    }

    @Override
    protected void tidyUpCacheIfDirty() {
        if (dirtyStructure) {
            super.resetCache();
            this.initialize();
        }
    }

    @Override
    protected void checkCladePartitionRemoval(Clade clade, CladePartition partition) {
        // do nothing here
    }

    /**
     * Set up this CCD0 after adding/removing trees. Assumes reset/not-set
     * values in CCD graph.
     */
    @Override
    public void initialize() {
        // need to find all clade partitions that could exist but were not
        // observed in base trees
        //System.out.print("\nExpanding CCD graph ... ");
        expand();

        // then need to set clade partition probabilities
        // which normalizes the product of clade probabilities
        //System.out.print("setting probabilities ... ");
        setPartitionProbabilities(this.rootClade);

        //System.out.println("done.\n");
        this.dirtyStructure = false;
    }

    /*
     * Expand CCD graph with clade partitions where parent and children were
     * observed, but not that clade partition.
     */
    private void expand() {
        // System.out.print("Expand ... ");
        int n = this.getNumberOfLeaves();

        // 1. clade buckets
        // for easier matching of child clades, we want to group them by size
        // 1.i sort clades
        Clade [] clades = cladeMapping.values().stream().sorted(Comparator.comparingInt(x -> x.getCladeInBits().cardinality()))
                .toArray(Clade[]::new);

		expand(clades, n);
        
//        for (int i = 0; i < clades.length; i++) {
//        	final Clade parent = clades[i];
//        	if (parent.isMonophyletic()) {
//        		List<Clade> clades0 = new ArrayList<>();
//        		for (int j = 0; j < i; j++) {
//        			if (clades[j] != null && parent.containsClade(clades[j])) {
//        				clades0.add(clades[j]);
//        				clades[j] = null;
//        			}
//        		}
//        		clades0.add(parent);
//        		expand(clades0, n);
//        	}
//        }
    }
        
    private void expand(Clade [] clades, int n) {
    	progressStream.println("processing " + clades.length + " clades");

        // 1.ii init clade buckets
        List<Set<Clade>> cladeBuckets = new ArrayList<Set<Clade>>(n);
        for (int i = 0; i < n; i++) {
            if (i < 3) {
                cladeBuckets.add(new HashSet<Clade>(n));
            } else {
                cladeBuckets.add(new HashSet<Clade>());
            }
        }
        // 1.iii fill clade buckets
        for (Clade clade : clades) {
            cladeBuckets.get(clade.size() - 1).add(clade);
        }

        // 2. find missing clade partitions
        threadCount = Runtime.getRuntime().availableProcessors();
    	done = new HashSet<>();
        if (threadCount <= 1 || clades.length < 20000) {
        	findPartitions(clades, cladeBuckets);
        } else {
            try {
            	System.out.println("Running with " + threadCount + " threads");
                countDown = new CountDownLatch(threadCount);
                ExecutorService exec = Executors.newFixedThreadPool(threadCount);
                // kick off the threads
                int end = clades.length;
                for (int i = 0; i < threadCount; i++) {
                    CoreRunnable coreRunnable = new CoreRunnable(clades, cladeBuckets, i, end);
                    exec.execute(coreRunnable);
                }
                countDown.await();
                exec. shutdownNow();
            } catch (RejectedExecutionException | InterruptedException e) {
            }
        	
        }

        // System.out.println("done.");
    }
    
    
    private CountDownLatch countDown = null;
    private int progressed = 0;
    private int threadCount = 4;
    private Set<Clade> done;

    class CoreRunnable implements java.lang.Runnable {
    	Clade[] clades;
    	List<Set<Clade>> cladeBuckets;
    	int start;
    	int end;
    	
        CoreRunnable(Clade[] clades, List<Set<Clade>> cladeBuckets, int start, int end) {
        	this.clades = clades;
        	this.cladeBuckets = cladeBuckets;
        	this.start = start;
        	this.end = end;
        }

        @Override
		public void run() {
        	int i = start;
        	try {
	            BitSet helperBits = new BitSet(clades[0].getCladeInBits().length());
	        	while (i < end) {
	            	findPartitions(clades, cladeBuckets, helperBits, i);
	        		i += threadCount;
	        		
	                if (progressStream != null) {
	                	while (progressed < i * 61 / clades.length) {
	                    	progressStream.print("*");
	                    	progressed++;
	                	}
	                }
	        	}
        	} catch (Throwable t) {
        		t.printStackTrace();
        	}
        	// progressStream.println("finished " + start + " " + i);
            countDown.countDown();
        }

    } // CoreRunnable


    private void findPartitions(Clade[] clades, List<Set<Clade>> cladeBuckets, BitSet helperBits, int i) {

        Clade parent = clades[i];

        // we skip leaves and cherries as they have no/only one partition
        if (parent.isLeaf() || parent.isCherry()) {
            return;
        }

        BitSet parentBits = parent.getCladeInBits();

        // otherwise we check if we find a larger partner clade for any
        // smaller clade that together partition the parent clade;
        for (int j = 1; j <= parent.size() / 2; j++) {
        	processSubClades(helperBits, cladeBuckets, j, parentBits, parent);
        }

            // remove clades below monophyletic clades
            if (parent.isMonophyletic()) {
                Set<Clade> descendants = parent.getDescendantClades(true);
                for (Clade descendant : descendants) {
                	synchronized(this) {
                		done.add(descendant);
                		//cladeBuckets.get(descendant.size() - 1).remove(descendant);
                	}
                }
            }
    }
    
    private void findPartitions(Clade [] clades, List<Set<Clade>> cladeBuckets) {
        BitSet helperBits = new BitSet(clades[0].getCladeInBits().length());

        // we go through clades in increasing size, then check for each clade of
        // at most half potential parent's size, whether we can find a partner
        int progressed = 0;
        for (int i = 0; i < clades.length; i++) {
            Clade parent = clades[i];
            if (progressStream != null) {
            	while (progressed < i * 61 / clades.length) {
                	progressStream.print("*");
                	progressed++;
            	}
            }

            // we skip leaves and cherries as they have no/only one partition
            if (parent.isLeaf() || parent.isCherry()) {
                continue;
            }

            BitSet parentBits = parent.getCladeInBits();

            // otherwise we check if we find a larger partner clade for any
            // smaller clade that together partition the parent clade;
            for (int j = 1; j <= parent.size() / 2; j++) {
            	processSubClades(helperBits, cladeBuckets, j, parentBits, parent);
            }

            // remove clades below monophyletic clades
            if (parent.isMonophyletic()) {
                Set<Clade> descendants = parent.getDescendantClades(true);
                for (Clade descendant : descendants) {
                    cladeBuckets.get(descendant.size() - 1).remove(descendant);
                }
            }
        }
        if (progressStream != null) {
        	progressStream.println();
        }
    }

    private void processSubClades(BitSet helperBits, List<Set<Clade>> cladeBuckets, int j, BitSet parentBits, Clade parent) {
    	
//    	Clade[] children = null;
//    	
//    	synchronized(this) {
//    		children = cladeBuckets.get(j - 1).toArray(new Clade[] {});
//    	}
    	
        for (Clade smallChild : cladeBuckets.get(j - 1)) {
        	if (!done.contains(smallChild)) {
	            BitSet smallChildBits = smallChild.getCladeInBits();
	
	            helperBits.clear();
	            helperBits.or(parentBits);
	            helperBits.and(smallChildBits);
	
	            if (helperBits.equals(smallChildBits)
	                    && !smallChild.parentClades.contains(parent)) {
	                // the small child is a subclade of the parent clade
	                // but is not yet in a partition of the parent clade;
	                // hence we check all potential partner clades
	
	                // here helperBits equal smallChildBits, so with an XOR
	                // with the parentBits we get the bits of the potential partner clade
	                helperBits.xor(parentBits);
	
	                if (cladeMapping.containsKey(helperBits)) {
	                    Clade largeChild = cladeMapping.get(helperBits);
	                    if (threadCount > 1) {
	                    	synchronized(this) {
	                    		parent.createCladePartition(smallChild, largeChild);
	                    	}
	                    } else {
	                    	parent.createCladePartition(smallChild, largeChild);
	                    }
	                }
	            }
        	}
        }
	}

	private static void findPartitionsOld(Clade[] clades, List<List<Clade>> cladeBuckets, int[] expandedClades) {
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
            if (clade.partitions.isEmpty()) {
                throw new AssertionError("Cherry should contain a clade split.");
            }
            clade.partitions.get(0).setCCP(1);
            clade.setSumCladeCredibilities(clade.getCladeCredibility());
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
                double probability;
                if (sumSubtreeProbabilities == 0) {
                    throw new AssertionError("Sum of subtree probabilities should not be 0.");
                }
                probability = sumPartitionSubtreeProbabilities[i] / sumSubtreeProbabilities;
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
     * @param heightStrategy the strategy used to set the heights of the tree vertices
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
