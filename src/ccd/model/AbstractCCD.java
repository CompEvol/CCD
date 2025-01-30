package ccd.model;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator.TreeSet;

import java.io.IOException;
import java.io.PrintStream;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Random;
import java.util.Set;

/**
 * <p>
 * This abstract class represents a tree distribution based on CCD graph. It
 * implements the tree distribution methods of {@link ITreeDistribution} like
 * sampling and getting the MAP tree. For this it uses the probabilities set in
 * the CCD graph by the implementing classes.
 * </p>
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
 * @author Jonathan Klawitter
 */
public abstract class AbstractCCD implements ITreeDistribution {

    /** Key for posterior support (clade probability) of vertex based on CCD. */
    public static final String CLADE_SUPPORT_KEY = "posterior";

    /** Whether to print information during construction, etc. */
    public static boolean verbose = true;

    /** Stream to print information to. */
    public static PrintStream out = System.out;

    /** Threshold used for rounding errors when adding up probabilities, to round back to 0 or 1. */
    public final static double PROBABILITY_ROUNDING_EPSILON = 1e-10;

    /** Threshold used for throwing error when probability that much out of bounds (mostly above 1). */
    public final static double PROBABILITY_ERROR = 1e-5;

    /**
     * The trees this CCD is based on (burnin trees removed).
     */
    protected List<Tree> baseTrees;

    /** The trees this CCD is based on (burnin trees removed). */
    protected TreeSet baseTreeSet;

    /**
     * Whether to store the trees used to create this CCD.
     */
    private final boolean storeBaseTrees;

    /**
     * The number of trees this CCD is based on.
     */
    protected int numBaseTrees = 0;

    /**
     * Burnin for trees (percentage of ignored first trees).
     */
    protected double burnin = 0;

    /**
     * Clade on all taxa, representing the root of the CCD graph.
     */
    protected Clade rootClade;

    /**
     * Number of leaves/taxa the trees this CCD is based on.
     */
    protected int leafArraySize;

    /**
     * Mapping from a BitSet representation to a clade; used to ensure
     * uniqueness of clades. Assumes that all BitSets have length equal to the
     * number of leaves.
     */
    protected Map<BitSet, Clade> cladeMapping;

    /**
     * Random used for sampling and tie breaking
     */
    protected Random random = new Random();


    /* -- CONSTRUCTORS & CONSTRUCTION METHODS -- */

    /**
     * Constructor for a {@link AbstractCCD} based on the given collection of
     * trees with specified burn-in.
     *
     * @param trees  the trees used to build and populate the CCD graph of this
     *               {@link AbstractCCD}
     * @param burnin value between 0 and 1 of what percentage of the given trees
     *               should be discarded as burn-in
     */
    public AbstractCCD(List<Tree> trees, double burnin) {
        this(trees.get(0).getLeafNodeCount(), true);

        this.burnin = burnin;
        List<Tree> treesToUse;
        if (burnin == 0) {
            treesToUse = trees;
        } else {
            int numDiscardedTrees = (int) (trees.size() * burnin);
            int numUsedTrees = trees.size() - numDiscardedTrees;
            this.numBaseTrees = numUsedTrees;
            treesToUse = new ArrayList<Tree>(numUsedTrees);
            treesToUse.addAll(trees.subList(numDiscardedTrees, trees.size()));
        }

        for (Tree tree : treesToUse) {
            cladifyTree(tree);
        }
    }

    /**
     * Constructor for a {@link AbstractCCD} based on the given collection of
     * trees (not containing any burnin trees).
     *
     * @param treeSet        an iterable set of trees, which contains no burnin trees and
     *                       that are used to build and populate the CCD graph of this
     *                       {@link AbstractCCD}
     * @param storeBaseTrees whether to store the trees used to create this CCD
     */
    public AbstractCCD(TreeSet treeSet, boolean storeBaseTrees) {
        this(storeBaseTrees);
        this.baseTreeSet = treeSet;
        this.burnin = 0;
        try {
            treeSet.reset();
            Tree tree = treeSet.next();
            initializeRootClade(tree.getLeafNodeCount());

            if (verbose) {
                out.println("Constructing CCD with " + (treeSet.totalTrees - treeSet.burninCount) + " trees...");
            }

            while (tree != null) {
                this.numBaseTrees++;
                cladifyTree(tree);

                // report progress
                if (verbose) {
                    if (numBaseTrees % 10 == 0) {

                        System.out.print(".");
                        out.flush();
                    }
                    if (numBaseTrees % 1000 == 0) {
                        out.println(" (" + numBaseTrees + ")");
                    }
                }

                tree = treeSet.hasNext() ? treeSet.next() : null;
            }

            if (verbose) {
                out.println(" ...done.");
            }

        } catch (IOException e) {
            System.err.println("Error reading in trees to create CCD.");
        }
    }

    /**
     * Constructor to start with an empty CDD graph. Trees can then be processed
     * one by one.
     *
     * @param numLeaves  number of leaves of the trees that this CCD will be based on
     * @param storeTrees whether to store the trees used to create this CCD;
     *                   recommended not to when huge set of trees is used
     */
    public AbstractCCD(int numLeaves, boolean storeTrees) {
        this(storeTrees);
        initializeRootClade(numLeaves);
    }

    /* Base constructor */
    protected AbstractCCD(boolean storeTrees) {
        this.storeBaseTrees = storeTrees;
        this.cladeMapping = new HashMap<BitSet, Clade>();
        this.baseTrees = new ArrayList<Tree>(storeTrees ? 1000 : 1);
    }

    /* Initialization helper method */
    protected void initializeRootClade(int numLeaves) {
        this.leafArraySize = numLeaves;

        BitSet rootBitSet = BitSet.newBitSet(leafArraySize);
        rootBitSet.set(0, numLeaves);

        this.rootClade = new Clade(rootBitSet, this);
        cladeMapping.put(rootClade.getCladeInBits(), rootClade);
    }

    /**
     * Add and process the given tree into the CCD graph.
     *
     * @param tree to be added and processed into the CCD graph
     */
    public void addTree(Tree tree) {
        this.numBaseTrees++;
        this.cladifyTree(tree);
        this.setCacheAsDirty();
    }

    /* Helper method; process one tree into this CCD */
    protected void cladifyTree(Tree tree) {
        if (storeBaseTrees) {
            this.baseTrees.add(tree);
        } else if (this.baseTrees.isEmpty()) {
            this.baseTrees.add(tree);
        }
        cladifyVertex(tree.getRoot());
    }

    /* Recursive helper method */
    private Clade cladifyVertex(Node vertex) {
        BitSet cladeInBits = BitSet.newBitSet(leafArraySize);
        Clade firstChildClade = null;
        Clade secondChildClade = null;

        if (vertex.isLeaf()) {
            int index = vertex.getNr();
            cladeInBits.set(index);
        } else {
            firstChildClade = cladifyVertex(vertex.getChildren().get(0));
            secondChildClade = cladifyVertex(vertex.getChildren().get(1));

            cladeInBits.or(firstChildClade.getCladeInBits());
            cladeInBits.or(secondChildClade.getCladeInBits());
        }

        Clade currentClade = cladeMapping.get(cladeInBits);
        if (currentClade == null) {
            currentClade = addNewClade(cladeInBits);
        }
        currentClade.increaseOccurrenceCount(vertex.getHeight());

        if (!vertex.isLeaf()) {
            CladePartition currentPartition = currentClade.getCladePartition(firstChildClade,
                    secondChildClade);
            if (currentPartition == null) {
                currentPartition = currentClade.createCladePartition(firstChildClade,
                        secondChildClade);
            }
            currentPartition.increaseOccurrenceCount(vertex.getHeight());
        }

        return currentClade;
    }

    /**
     * Adds and returns a new clade to this CCD based on the given BitSet.
     * Assumes that the clade does not exist yet; otherwise future behaviour is undefined.
     *
     * @param cladeInBits BitSet describing new clade; should not exist yet
     * @return newly added clade
     */
    protected Clade addNewClade(BitSet cladeInBits) {
        Clade clade = new Clade(cladeInBits, this);
        cladeMapping.put(cladeInBits, clade);
        return clade;
    }

    /**
     * Removes the given tree from set of trees the CCD graph is based on
     * (reduces the number of occurrences for each of the tree's clades and
     * partitions by one). The behavior is unspecified if tree wasn't used
     * to build CCD graph or has been removed (to a count of 0) before.
     *
     * @param tree           to be taken out from set of trees this CCD is based on
     * @param tidyUpCCDGraph whether to keep the CCD graph tidy, that is, not keeping any
     *                       nontrivial clades without child clade partitions and parent
     *                       clades
     */
    public void removeTree(Tree tree, boolean tidyUpCCDGraph) {
        if (this.storeBaseTrees && !this.baseTrees.remove(tree)) {
            System.err.println("WARNING: Removing tree from CCD that was not part of it.");
        }

        this.reduceCladeCount(tree.getRoot());
        this.numBaseTrees--;

        if (tidyUpCCDGraph) {
            this.tidyUpCCDGraph(false);
        }

        this.setCacheAsDirty();
    }

    /* Recursive helper method */
    private Clade reduceCladeCount(Node vertex) {
        // 1. build BitSet to retrieve clade and call recursion
        BitSet cladeInBits = BitSet.newBitSet(leafArraySize);
        Clade firstChildClade = null;
        Clade secondChildClade = null;

        if (vertex.isLeaf()) {
            int index = vertex.getNr();
            cladeInBits.set(index);
        } else {
            firstChildClade = reduceCladeCount(vertex.getChildren().get(0));
            secondChildClade = reduceCladeCount(vertex.getChildren().get(1));

            cladeInBits.or(firstChildClade.getCladeInBits());
            cladeInBits.or(secondChildClade.getCladeInBits());
        }

        // 2. retrieve clade and reduce count
        Clade currentClade = this.cladeMapping.get(cladeInBits);
        currentClade.decreaseOccurrenceCount(vertex.getHeight());

        // 3. reduce counts for its clade partitions
        if (!vertex.isLeaf()) {
            CladePartition currentPartition = currentClade.getCladePartition(firstChildClade, secondChildClade);
            currentPartition.decreaseOccurrenceCount(vertex.getHeight());

            removeCladePartitionIfNecessary(currentClade, currentPartition);
        }

        return currentClade;
    }

    /**
     * Removes a clade partition if it shouldn't be kept in a CCD,
     * e.g. after its count has been decreased.
     *
     * @param clade     parent of the clade partition
     * @param partition to be checked and potentially removed
     * @return whether the clade partition has been removed
     */
    protected abstract boolean removeCladePartitionIfNecessary(Clade clade, CladePartition partition);

    /**
     * Tidy up this CCD by removing nontrivial clades that have no clade
     * partitions and clade partitions whose parent or child clades have been
     * removed.
     *
     * @param renormalize whether to renormalize CCPs of clade partitions under a clade where ones was removed
     * @return whether this CCD is still complete, i.e. returns {@code false}
     * when any leaf is not reachable anymore from the root
     */
    public boolean tidyUpCCDGraph(boolean renormalize) {
        boolean complete = true;
        ArrayList<Clade> cladesToRemove = new ArrayList<Clade>();

        int numCladesRemoved = 0;
        int numPartitionsRemoved = 0;

        // check if any (nontrivial) clade should be removed, because
        // - has no parents anymore
        // - has no children (clades/clade partitions) anymore
        for (Clade clade : this.getCladeMapping().values()) {
            if (clade.isLeaf()) {
                if (clade.getParentClades().isEmpty()) {
                    // System.out.println("- leaf "
                    //         + clade.getCladeInBits() + " (" + this.baseTrees.get(0)
                    //         .getNode(clade.getCladeInBits().nextSetBit(0)).getID()
                    //         + ") has no parent clade!");
                    complete = false;
                }
            } else if (clade == this.rootClade) {
                if (clade.getPartitions().isEmpty()) {
                    out.println("- root clade has no partitions!");
                    // complete = false;
                }
            } else if (clade.getParentClades().isEmpty() || clade.getPartitions().isEmpty()
                    || (clade.getNumberOfOccurrences() == 0)) {
                // out.println("- clade found to remove: " + clade.getCladeInBits());
                cladesToRemove.add(clade);
            }
        }

        // when removing a clades, have to check their (old) parents/children
        while (!cladesToRemove.isEmpty()) {
            Clade cladeToRemove = cladesToRemove.remove(cladesToRemove.size() - 1);
            // out.println("clade to remove: " + cladeToRemove.getCladeInBits());

            // keep trivial clades
            if (cladeToRemove.isLeaf() || (cladeToRemove == this.rootClade)) {
                out.println("- request to remove " + (cladeToRemove.isLeaf() ? "leaf" : "root") + "!");
                complete = false;
                continue;
            }

            // a) take clade out of mapping,
            // which we might have already done before
            // (clade could be added to removal list multiple times)
            if (this.cladeMapping.remove(cladeToRemove.getCladeInBits()) == null) {
                continue;
            }
            setCacheAsDirty();
            numCladesRemoved++;

            // b) remove connection to parent clades ...
            for (Clade parent : cladeToRemove.getParentClades()) {
                if (!parent.getChildClades().remove(cladeToRemove)) {
                    continue;
                }

                // ... and update parent clades
                for (CladePartition parentPartition : parent.getPartitions()) {
                    // there can only be one partition that contains
                    // cladeToRemove under the parent clade
                    if (parentPartition.containsChildClade(cladeToRemove)) {
                        parent.getPartitions().remove(parentPartition);
                        numPartitionsRemoved++;
                        if (parent.getPartitions().isEmpty()) {
                            cladesToRemove.add(parent);
                        } else {
                            if (renormalize) {
                                parent.normalizeCCPs();
                            }
                        }

                        // also update other child
                        Clade otherChild = parentPartition.getOtherChildClade(cladeToRemove);
                        parent.getChildClades().remove(otherChild);
                        otherChild.getParentClades().remove(parent);
                        if (otherChild.getParentClades().isEmpty()) {
                            cladesToRemove.add(otherChild);
                        }

                        break;
                    }
                }

            }

            // c) remove all connection of children to parent
            for (CladePartition partition : cladeToRemove.getPartitions()) {
                for (Clade child : partition.getChildClades()) {
                    child.getParentClades().remove(cladeToRemove);
                    if (child.getParentClades().isEmpty()) {
                        cladesToRemove.add(child);
                    }
                }
            }
            // do not have to empty partitions of clade, since it is dumped anyway
        }

        /*- if ((numCladesRemoved > 0) || (numPartitionsRemoved > 0)) {
            out.println("Tidying up CCD - removed " + numCladesRemoved + " clades and " + numPartitionsRemoved + " clade partitions.");
            if (!complete) {
                out.println("=> CCD is not complete anymore (contains no tree or only incomplete ones)!");
            }
        }*/

        return complete;
    }

    /**
     * @param random used in this CCD for sampling
     */
    public void setRandom(Random random) {
        this.random = random;
    }


    /* -- GENERAL & CCD GRAPH GETTERS -- */

    /**
     * @return number of leaves/taxa of the trees this CCD is build on (which
     * might be less than the taxa existing in this CCD, namely, if it
     * is filtered)
     */
    @Override
    public int getNumberOfLeaves() {
        return leafArraySize;
    }

    /**
     * Returns the size of the arrays working with the leaves.
     * This should differ only from the number of leaves (cf. {@link AbstractCCD#getNumberOfLeaves()})
     * if leaves have been removed, for example, when filtering in a {@link FilteredCCD}.
     *
     * @return the size of the arrays working with the leaves
     */
    public int getSizeOfLeavesArray() {
        return leafArraySize;
    }

    /**
     * @return the root clade of the CCD graph
     */
    public Clade getRootClade() {
        return rootClade;
    }

    @Override
    public int getNumberOfClades() {
        return cladeMapping.size();
    }

    @Override
    public Collection<Clade> getClades() {
        return cladeMapping.values();
    }

    /**
     * @return mapping from BitSet to clade
     */
    public Map<BitSet, Clade> getCladeMapping() {
        return cladeMapping;
    }

    /**
     * @return the number of distinct clade partitions of this CCD
     */
    public int getNumberOfCladePartitions() {
        int count = 0;
        for (Clade clade : cladeMapping.values()) {
            count += clade.getNumberOfPartitions();
        }
        return count;
    }

    /**
     * @return the taxa of this CCD as BitSet with length based on original
     * number of leaves (only different if this a
     * {@link FilteredCCD})
     */
    public BitSet getTaxaAsBitSet() {
        return this.rootClade.getCladeInBits();
    }

    /**
     * @return the clade represented by the given BitSet if it is in this CCD
     * and {@code null} otherwise
     */
    public Clade getClade(BitSet cladeAsBitSet) {
        return this.cladeMapping.get(cladeAsBitSet);
    }

    /**
     * @return number of trees this CCD is based on
     */
    public int getNumberOfBaseTrees() {
        return this.numBaseTrees;
    }

    /**
     * <i>Use only with extra care and if you know what you are doing!</i>
     * Manually set the number of base trees of this CCD.
     *
     * @param newNumberOfBaseTrees new value
     */
    public void setNumBaseTrees(int newNumberOfBaseTrees) {
        this.numBaseTrees = newNumberOfBaseTrees;
    }

    /** @return trees this CCD is based on if stored, otherwise {@code null} */
    public List<Tree> getBaseTrees() {
        if (this.baseTrees == null) {
            return null;
        } else {
            return (baseTrees.size() == this.numBaseTrees) ? this.baseTrees : null;
        }
    }

    /** @return tree set used to construct this CCD if stored, otherwise {@code null} */
    public TreeSet getBaseTreeSet() {
        return baseTreeSet;
    }

    /**
     * Set the base TreeSet of this CCD.
     * Should be the actual TreeSet used to construct,
     * which might have been done tree by tree for some reason.
     *
     * @param baseTreeSet should be the one used to construct the CCD
     */
    public void setBaseTreeSet(TreeSet baseTreeSet) {
        this.baseTreeSet = baseTreeSet;
    }

    /**
     * @return some tree this CCD is based on, which might be needed to access
     * tree and taxa information
     */
    public Tree getSomeBaseTree() {
        return baseTrees.get(0);
    }

    /** @return whether this CCD stores all trees used to construct it */
    protected boolean storesBaseTrees() {
        return storeBaseTrees;
    }


    /* -- STATE MANAGEMENT - STATE MANAGEMENT -- */

    /** Whether cached probability values are out of date. */
    protected boolean probabilitiesDirty = false;

    /** Whether cached entropy values are out of date. */
    protected boolean entropyDirty = false;

    /** Whether cached numbers of topologies are out of date. */
    protected boolean numberOfTopologiesDirty = false;

    /** Whether cached numbers of topologies are out of date. */
    protected boolean commonAncestorHeightsDirty = true;

    /**
     * Sets CCD as dirty lazily, meaning cached values become out of date and
     * (lazily) clades and clade partitions are not told about being dirty!
     */
    public void setCacheAsDirty() {
        this.probabilitiesDirty = true;
        this.entropyDirty = true;
        this.numberOfTopologiesDirty = true;
        this.commonAncestorHeightsDirty = true;
        this.cladeMinCredibility = null;
    }

    /**
     * Tidy up CCD graph, probabilities, ... before computing a values.
     */
    protected abstract void tidyUpCacheIfDirty();

    /* Helper method. */
    protected void resetCache() {
        for (Clade clade : cladeMapping.values()) {
            clade.resetCachedValues();
        }
        probabilitiesDirty = false;
        entropyDirty = false;
        numberOfTopologiesDirty = false;
        commonAncestorHeightsDirty = true;
    }

    /* Helper method. */
    protected void resetCacheIfProbabilitiesDirty() {
        if (probabilitiesDirty) {
            resetCache();
        }
    }


    /* -- DISTRIBUTION VALUES -- */

    /**
     * Compute the entropy of the tree distribution modeled by this CCD. Based
     * on the algorithm by <a href="https://dx.doi.org/10.1093/sysbio/syw042">Lewis et al.
     * 2016</a>, Appendix 1 # Algorithm 1
     *
     * @return the entropy of the tree distribution modeled by this CCD
     */
    public double getEntropyLewis() {
        if (entropyDirty) {
            resetCache();
        }

        return this.rootClade.getEntropy();
    }

    /**
     * Compute the entropy of the tree distribution modeled by this CCD
     * with the simple clade partition based formula.
     *
     * @return the entropy of the tree distribution modeled by this CCD
     */
    public double getEntropy() {
        computeCladeProbabilitiesIfDirty();
        double testro = 0;
        for (Clade clade : getClades()) {
            for (CladePartition partition : clade.getPartitions()) {
                double logS = partition.getLogCCP();
                testro += partition.getProbability() * logS;
            }
        }

        return -testro;
    }

    @Override
    public BigInteger getNumberOfTrees() {
        if (numberOfTopologiesDirty) {
            resetCache();
        }

        return this.rootClade.getNumberOfTopologies();
    }

    /**
     * Returns the AIC score of this CCD.
     * The number of parameters depends on the specific CCD.
     * The log likelihood is computed as the log of the product of the probability of each tree used to construct this CCD;
     * thus the method requires that this CCD knows the TreeSet it was constructed with.
     *
     * @return the AIC score of this CCD
     * @throws IOException
     */
    public double getAICScore() throws IOException {
        if (this.baseTreeSet == null) {
            System.err.println("Cannot compute AIC score as CCD not constructed from TreeSet.");
            return -1;
        }

        double twoK = 2 * this.getNumberOfParameters();
        double logL = 0.0;

        baseTreeSet.reset();
        while (baseTreeSet.hasNext()) {
            logL += Math.log(getProbabilityOfTree(baseTreeSet.next()));
        }

        return twoK - 2 * logL;
    }

    /** @return the number of parameters this CCD model has */
    abstract protected double getNumberOfParameters();

    /** @return the log likelihood of this CCD based on the trees it was constructed with */
    public double getLogLikelihood() {
        double logP = 0;
        for (Clade c : this.getClades()) {
            for (CladePartition p : c.getPartitions()) {
                int x = p.getNumberOfOccurrences();
                if (x > 0) {
                    logP += x * p.getLogCCP();
                }
            }
        }
        return logP;
    }

    /**
     * Computes and returns the Fair Proportion Diversity Index of the taxa in this CCD
     * by using branch length derived from heights set by the given strategy
     *
     * @param heightStrategy used to set clade heights
     * @return Fair Proportion Diversity Index of the taxa in this CCD
     */
    public double[] getFairProportionIndex(HeightSettingStrategy heightStrategy) {
        this.computeCladeProbabilitiesIfDirty();
        if (heightStrategy == HeightSettingStrategy.CommonAncestorHeights) {
            setupCommonAncestorHeightsIfDirty();
        }

        double max = 0;

        double[] index = new double[this.getSizeOfLeavesArray()];

        for (Clade parent : this.getClades()) {
            double pClade = parent.getProbability();

            for (CladePartition partition : parent.getPartitions()) {
                double pPartition = partition.getCCP();

                for (Clade child : partition.getChildClades()) {
                    double branchLength = 0;
                    if (heightStrategy == HeightSettingStrategy.CommonAncestorHeights) {
                        branchLength = parent.getCommonAncestorHeight() - child.getCommonAncestorHeight();
                    } else if (heightStrategy == HeightSettingStrategy.MeanOccurredHeights) {
                        branchLength = parent.getMeanOccurredHeight() - child.getMeanOccurredHeight();
                    }

                    if (branchLength < 0) {
                        throw new AssertionError("Negative branch length.");
                    }

                    int size = child.size();
                    double diversity = pClade * pPartition * branchLength / size;

                    BitSet bitset = child.getCladeInBits();
                    int i = bitset.nextSetBit(0);
                    while (i != -1) {
                        index[i] += diversity;
                        i = bitset.nextSetBit(i + 1);
                    }
                }
            }
        }

        return index;
    }


    /* -- POINT ESTIMATE / SAMPLING METHODS -- */

    @Override
    public Tree sampleTree() {
        return sampleTree(HeightSettingStrategy.None);
    }

    @Override
    public Tree sampleTree(HeightSettingStrategy heightStrategy) {
        return getTreeBasedOnStrategy(SamplingStrategy.Sampling, heightStrategy);
    }

    @Override
    public Tree getMAPTree() {
        return this.getMAPTree(HeightSettingStrategy.One);
    }

    @Override
    public Tree getMAPTree(HeightSettingStrategy heightStrategy) {
        return getTreeBasedOnStrategy(SamplingStrategy.MAP, heightStrategy);
    }

    /* Helper for methods to assign indices to inner vertices */
    private int runningInnerIndex;

    /* Helper for methods to assign indices to leaves */
    // private int runningLeafIndex;

    /* Strategy based tree sampling method */
    protected Tree getTreeBasedOnStrategy(SamplingStrategy samplingStrategy, HeightSettingStrategy heightStrategy) {
        tidyUpCacheIfDirty();
        computeCladeProbabilitiesIfDirty();

        if (heightStrategy == HeightSettingStrategy.CommonAncestorHeights) {
            setupCommonAncestorHeightsIfDirty();
        }

        runningInnerIndex = this.getSizeOfLeavesArray();
        Node root = getVertexBasedOnStrategy(this.rootClade, samplingStrategy, heightStrategy);

        if (this instanceof FilteredCCD) {
            return new FilteredTree(root);
        } else {
            return new Tree(root);
        }
    }

    /* Recursive helper method */
    private Node getVertexBasedOnStrategy(Clade clade, SamplingStrategy samplingStrategy, HeightSettingStrategy heightStrategy) {
        // computeCladeProbabilitiesIfDirty();

        Node vertex = null;
        if (clade.isLeaf()) {
            int leafNr = clade.getCladeInBits().nextSetBit(0);
            String taxonName = this.getSomeBaseTree().getTaxaNames()[leafNr];

            vertex = new Node(taxonName);
            vertex.setNr(leafNr);
            vertex.setMetaData(CLADE_SUPPORT_KEY, 1.0);
            // vertex.setNr(runningLeafIndex++);
            if (heightStrategy != null) {
                // common ancestor, mean height and ONE all use the same height for leaves
                vertex.setHeight(clade.getMeanOccurredHeight());
            }
        } else {
            CladePartition partition = getPartitionBasedOnStrategy(clade, samplingStrategy);
            if (partition == null) {
                throw new AssertionError("Unsuccessful to find clade partition of clade: " + clade.getCladeInBits());
            }

            Node firstChild = getVertexBasedOnStrategy(partition.getChildClades()[0],
                    samplingStrategy, heightStrategy);
            Node secondChild = getVertexBasedOnStrategy(partition.getChildClades()[1],
                    samplingStrategy, heightStrategy);

            // These are not needed and only make the output newick longer
            vertex = new Node();
            vertex.setNr(runningInnerIndex++);
            double cladeProbability = clade.getProbability();
            vertex.setMetaData(CLADE_SUPPORT_KEY, cladeProbability);
            String posteriorSupport = CLADE_SUPPORT_KEY + "=" + cladeProbability;
            if (vertex.metaDataString != null) {
                vertex.metaDataString += "," + posteriorSupport;
            } else {
                vertex.metaDataString = posteriorSupport;
            }
            vertex.addChild(firstChild);
            vertex.addChild(secondChild);

            if (heightStrategy == HeightSettingStrategy.MeanOccurredHeights) {
                vertex.setHeight(clade.getMeanOccurredHeight());
            } else if (heightStrategy == HeightSettingStrategy.One) {
                double height = Math.max(firstChild.getHeight(), secondChild.getHeight()) + 1;
                vertex.setHeight(height);
            } else if (heightStrategy == HeightSettingStrategy.CommonAncestorHeights) {
                // out.println("\nvertex = " + vertex);
                // out.println("vertex.getHeight() = " + vertex.getHeight());
                // out.println("clade.getCommonAncestorHeight() = " + clade.getCommonAncestorHeight());
                vertex.setHeight(clade.getCommonAncestorHeight());
                // out.println("vertex.getHeight() = " + vertex.getHeight());

                if (Double.isNaN(clade.getCommonAncestorHeight())) {
                    System.err.println("\nNaN height!");
                    System.err.println("clade = " + clade);
                    System.err.println("clade.getCommonAncestorHeight() = " + clade.getCommonAncestorHeight());
                }

                if (vertex.getHeight() < 0) {
                    System.err.println("\nVertex with negative height");
                    System.err.println("vertex.getHeight() = " + vertex.getHeight());
                    System.err.println("clade.getCommonAncestorHeight =  " + clade.getCommonAncestorHeight());
                    System.err.println("clade.getMeanOccurredHeight =  " + clade.getMeanOccurredHeight());
                }
                if ((vertex.getHeight() - vertex.getLeft().getHeight()) < 0) {
                    System.err.println("\nNegative branch length, L");
                    System.err.println("branchLength = " + (vertex.getHeight() - vertex.getLeft().getHeight()));
                    System.err.println("parent = " + vertex);
                    System.err.println("childL = " + vertex.getLeft());
                }
                if ((vertex.getHeight() - vertex.getRight().getHeight()) < 0) {
                    System.err.println("\nNegative branch length, R");
                    System.err.println("branchLength = " + (vertex.getHeight() - vertex.getLeft().getHeight()));
                    System.err.println("parent = " + vertex);
                    System.err.println("childR = " + vertex.getRight());
                }
            }
        }

        return vertex;
    }

    /**
     * Returns the probability of the most likely tree. Note that this can
     * underflow for large trees. It is recommended to use {@link #getMaxLogTreeProbability()}
     * instead.
     *
     * @return probability of the most likely tree.
     */
    @Override
    public double getMaxTreeProbability() {
        return Math.exp(this.getMaxLogTreeProbability());
    }

    /**
     * @return the log probability of the most likely tree.
     */
    public double getMaxLogTreeProbability() {
        tidyUpCacheIfDirty();
        resetCacheIfProbabilitiesDirty();

        return this.rootClade.getMaxSubtreeLogCCP();
    }

    /* Helper method */
    private CladePartition getPartitionBasedOnStrategy(Clade clade, SamplingStrategy samplingStrategy) {
        CladePartition partition = null;
        switch (samplingStrategy) {
            case MAP: {
                partition = clade.getMaxSubtreeCCPPartition();
                break;
            }
            case Sampling: {
                ArrayList<CladePartition> partitions = clade.getPartitions();

                // the sum of probabilities over all partitions of a clade
                // should be 1, so we can sample with a random value
                double sampleWithMe = random.nextDouble();

                double probabilitySum = 0;
                for (CladePartition nextPartition : partitions) {
                    probabilitySum += nextPartition.getCCP();
                    if (sampleWithMe < probabilitySum) {
                        partition = nextPartition;
                        break;
                    }
                }

                // sum might not exactly add up to 1.0, so for robustness then
                // pick the last partition
                if (partition == null) {
                    partition = partitions.get(partitions.size() - 1);
                }

                break;
            }
            case MaxSumCladeCredibility:
                partition = clade.getMaxSubtreeSumCladeCredibilityPartition();
                break;
            default:
                throw new IllegalArgumentException("Unexpected value: " + samplingStrategy);
        }
        return partition;
    }

    /* Helper method. */
    protected void setupCommonAncestorHeightsIfDirty() {
        if (this.commonAncestorHeightsDirty) {
            setupCommonAncestorHeights();
            this.commonAncestorHeightsDirty = false;
        }
    }

    /**
     * Set the common ancestor heights in all clades based on the treeset or list of trees used to construct this CCD.
     * Note that if single trees were used, then this only works if the CCD was constructed with the parameter to store the base trees.
     */
    protected void setupCommonAncestorHeights() {
        if ((this.baseTreeSet == null) && (this.baseTrees == null || this.numBaseTrees != this.baseTrees.size())) {
            throw new AssertionError("Method to set common ancestor heights called, " +
                    "but neither are the base trees nor a base treeset stored.");
        }

        // - make sure heights are properly reset
        // - for monophyletic clades (which includes leaves and root),
        // the common ancestor height equals the mean observed height
        for (Clade clade : this.getClades()) {
            clade.setCommonAncestorHeight(0);
            if (clade.isLeaf() || clade.isRoot() || clade.isMonophyletic()) {
                clade.setCommonAncestorHeight(clade.getMeanOccurredHeight());
            }
        }

        if (storeBaseTrees) {
            for (Tree tree : baseTrees) {
                extractCommonAncestorHeightsFromTree(tree);
            }
        } else {
            try {
                baseTreeSet.reset();
                while (baseTreeSet.hasNext()) {
                    Tree tree = baseTreeSet.next();
                    extractCommonAncestorHeightsFromTree(tree);
                }

            } catch (IOException e) {
                throw new RuntimeException("Error opening/using trees file used to construct CCD.");
            }
        }

        // validation
        for (Clade clade : this.getClades()) {
            // out.println("clade: " + clade);
            // out.println("clade.h: " + clade.getCommonAncestorHeight());

            // if (clade.isLeaf() || clade.isRoot() || clade.isMonophyletic()) {
            //     if (Math.abs(clade.getMeanOccurredHeight() - clade.getCommonAncestorHeight()) != 0.0) {
            //         out.println("\nMean and CA heights differ where they shouldn't!");
            //         out.println("clade = " + clade);
            //         out.println("clade.getCommonAncestorHeight() = " + clade.getCommonAncestorHeight());
            //         out.println("clade.getMeanOccurredHeight() = " + clade.getMeanOccurredHeight());
            //     }
            // }

            if (clade.getCommonAncestorHeight() < 0) {
                out.println("\nNegative height!");
                out.println("clade = " + clade);
                out.println("clade.getCommonAncestorHeight() = " + clade.getCommonAncestorHeight());
            }

            if (clade.isRoot()) {
                continue;
            }

            for (Clade parent : clade.getParentClades()) {
                if (parent.getCommonAncestorHeight() - clade.getCommonAncestorHeight() < 0) {
                    out.println("\nNegative branch length!");
                    out.println("parent:   " + parent);
                    out.println("parent.h: " + parent.getCommonAncestorHeight());
                    out.println("clade:    " + clade);
                    out.println("clade.h:  " + clade.getCommonAncestorHeight());
                }
            }
        }
    }

    /* Helper method */
    private void extractCommonAncestorHeightsFromTree(Tree tree) {
        WrappedBeastTree wrappedTree = new WrappedBeastTree(tree);

        for (Clade clade : this.getClades()) {
            if (clade.isLeaf() || clade.isRoot() || clade.isMonophyletic()) {
                continue;
            }

            double height = wrappedTree.getCommonAncestorHeightOfClade(clade.getCladeInBits());
            double additive = height / this.getNumberOfBaseTrees();

            clade.setCommonAncestorHeight(clade.getCommonAncestorHeight() + additive);
        }
    }

    /**
     * Sets the heights of the given tree based on the given strategy.
     *
     * @param tree                  to receive heights
     * @param heightSettingStrategy (not all implemented yet TODO)
     */
    public Tree setHeights(Tree tree, HeightSettingStrategy heightSettingStrategy) {
        switch (heightSettingStrategy) {
            case MeanOccurredHeights -> setMeanOccurredHeights(tree);
            case CommonAncestorHeights -> setCommonAncestorHeights(tree);
            default -> System.err.println("This height setting strategy is not implemented yet."); // TODO
        }

        return tree;
    }

    /* Helper method */
    private void setCommonAncestorHeights(Tree tree) {
        setupCommonAncestorHeightsIfDirty();

        // initialize BitSet representation of trees
        WrappedBeastTree wrappedTree = new WrappedBeastTree(tree);
        WrappedBeastTree[] wrappedUsedTrees = new WrappedBeastTree[numBaseTrees];
        int i = 0;
        for (Tree inputTree : baseTrees) {
            wrappedUsedTrees[i++] = new WrappedBeastTree(inputTree);
        }

        for (Node vertex : tree.getNodesAsArray()) {
            if (vertex == null) {
                // if filtered, there can be null entries in vertex array
                continue;
            }

            double meanHeight = 0;
            if (vertex.isLeaf()) {
                for (Tree inputTree : baseTrees) {
                    double height = inputTree.getNode(vertex.getNr()).getHeight();
                    meanHeight += height / baseTrees.size();
                }
            } else {
                BitSet cladeInBits = wrappedTree.getCladeInBits(vertex.getNr());
                for (WrappedBeastTree wrappedUsedTree : wrappedUsedTrees) {
                    double height = wrappedUsedTree.getCommonAncestorHeightOfClade(cladeInBits);
                    meanHeight += height / baseTrees.size();
                }
            }

            vertex.setHeight(meanHeight);
        }
    }

    /* Helper method */
    private void setMeanOccurredHeights(Tree tree) {
        WrappedBeastTree wrappedTree = new WrappedBeastTree(tree);

        for (Node vertex : tree.getNodesAsArray()) {
            if (vertex == null) {
                // if filtered, there can be null entries in vertex array
                continue;
            }

            BitSet cladeInBits = wrappedTree.getCladeInBits(vertex.getNr());
            Clade clade = cladeMapping.get(cladeInBits);

            if (clade != null) {
                vertex.setHeight(clade.getMeanOccurredHeight());
            } else {
                System.err.println("Could not set height of vertex as clade " + cladeInBits
                        + " does not exist in CCD.");
            }
        }
    }

    /**
     * Returns a map from each clade with a respective vertex in the given tree to that vertex.
     *
     * @param tree to whose vertices we map from clades of this CCD
     * @return a map from each clade with a respective vertex in the given tree to that vertex
     */
    public Map<Clade, Node> getCladeToNodeMap(Tree tree) {
        Map<Clade, Node> map = new HashMap<>(tree.getNodeCount());
        computeCladeToNodeMapping(tree.getRoot(), map);
        return map;
    }

    /* Helper method */
    private BitSet computeCladeToNodeMapping(Node vertex, Map<Clade, Node> map) {
        BitSet bits;
        if (vertex.isLeaf()) {
            bits = BitSet.newBitSet(getSizeOfLeavesArray());
            bits.set(vertex.getNr());
        } else {
            bits = computeCladeToNodeMapping(vertex.getLeft(), map);
            BitSet otherBits = computeCladeToNodeMapping(vertex.getRight(), map);
            bits.or(otherBits);
        }

        Clade clade = getClade(bits);
        if (clade == null) {
            System.err.println("No clade found in CCD for this vertex (" + bits + ").");
        }
        map.put(clade, vertex);

        return bits;
    }


    /* -- PROBABILITY - PROBABILITY -- */

    @Override
    public double getProbabilityOfTree(Tree tree) {
        resetCacheIfProbabilitiesDirty();

        double[] runningProbability = new double[]{1};
        computeProbabilityOfVertex(tree.getRoot(), runningProbability);

        return runningProbability[0];
    }

    /* Recursive helper method */
    private Clade computeProbabilityOfVertex(Node vertex, double[] runningProbability) {
        BitSet cladeInBits = BitSet.newBitSet(leafArraySize);

        if (vertex.isLeaf()) {
            int index = vertex.getNr();
            cladeInBits.set(index);

            // leaf has probability 1, so no changes to runningProbability

            return cladeMapping.get(cladeInBits);
        } else {
            Clade firstChildClade = computeProbabilityOfVertex(vertex.getChildren().get(0),
                    runningProbability);
            Clade secondChildClade = computeProbabilityOfVertex(vertex.getChildren().get(1),
                    runningProbability);

            if ((firstChildClade == null) || (secondChildClade == null)) {
                runningProbability[0] = 0;
                return null;
            }

            cladeInBits.or(firstChildClade.getCladeInBits());
            cladeInBits.or(secondChildClade.getCladeInBits());

            Clade currentClade = cladeMapping.get(cladeInBits);
            if (currentClade != null) {

                CladePartition partition = currentClade.getCladePartition(firstChildClade,
                        secondChildClade);
                if (partition != null) {
                    runningProbability[0] *= partition.getCCP();
                } else {
                    runningProbability[0] = 0;
                }
            } else {
                runningProbability[0] = 0;
            }

            return currentClade;
        }
    }

    @Override
    public boolean containsTree(Tree tree) {
        return (getProbabilityOfTree(tree) > 0);
    }

    @Override
    public double getCladeProbability(BitSet cladeInBits) {
        resetCacheIfProbabilitiesDirty();

        Clade clade = cladeMapping.get(cladeInBits);
        if (clade == null) {
            return 0;
        } else {
            if (clade.getProbability() < 0) {
                computeCladeProbabilities();
            }
            return clade.getProbability();
        }
    }

    /** Compute the probabilities of all clades in this distribution unless there are cached values. */
    public void computeCladeProbabilitiesIfDirty() {
        resetCacheIfProbabilitiesDirty();
        if (rootClade.getProbability() < 0) {
            computeCladeProbabilities();
        }
    }

    /** Compute the probabilities of all clades in this distribution. */
    public void computeCladeProbabilities() {
        resetCacheIfProbabilitiesDirty();

        for (Clade clade : this.getClades()) {
            if (!clade.isLeaf()) {
                clade.setProbability(-1);
            }
        }

        // the probability of a clade in a CCD is given by the sum of products
        // of probabilities along any path from the root to that clade;
        // hence, to compute it for each clade, we use a BFS-like traversal of
        // the CCD graph where a clade is handled only if all edges from parent
        // clades have been used (we keep track of this with counters)
        HashMap<Clade, Integer> visitCountMap = new HashMap<>(cladeMapping.size());
        Queue<Clade> queue = new LinkedList<>();

        rootClade.setProbability(1);
        queue.add(rootClade);
        visitCountMap.put(rootClade, 0);

        while (!queue.isEmpty()) {
            Clade clade = queue.poll();
            int count = visitCountMap.get(clade);
            if (count != clade.getNumberOfParentClades()) {
                // clade was not visited often enough or was already handled
                continue;
            } else if (clade.isLeaf()) {
                continue;
            }

            double parentProbability = clade.getProbability();

            for (CladePartition partition : clade.getPartitions()) {
                for (Clade childClade : partition.getChildClades()) {
                    if (childClade.isLeaf()) {
                        continue;
                    }

                    // if reset, value is -1, so have to start with 0
                    double childProbability = Math.max(0, childClade.getProbability());
                    // probability of child clade is sum
                    childProbability += parentProbability * partition.getCCP();

                    if (childProbability > 1 + PROBABILITY_ERROR) {
                        System.err.println("\nComputed invalid probability.");
                        System.err.println("parent clade = " + clade);
                        System.err.println("    with probability: " + parentProbability);
                        System.err.println("child clade = " + childClade);
                        System.err.println("    w old probability:   " + childClade.getProbability());
                        System.err.println("    and new probability: " + childProbability);
                        System.err.println("Partition = " + partition);
                        System.err.println("    with CCP: " + partition.getCCP());
                        System.err.println("List all partitions of parent clade: ");
                        for (CladePartition partition1 : clade.getPartitions()) {
                            System.err.println("\t" + partition1);
                        }
                        throw new AssertionError("Computation of invalid probability.");
                    }

                    if (childProbability > (1 + PROBABILITY_ROUNDING_EPSILON)) {
                        System.err.println("\nComputed probability for a clade above 1+epsilon.");
                        System.err.println("childClade = " + childClade);
                        System.err.println("oldChildProbability = " + childClade.getProbability());
                        System.err.println("newChildProbability = " + childProbability);
                        System.err.println("parentProbability = " + parentProbability);
                        System.err.println("partition.ccp = " + partition.getCCP());
                        System.err.println("Could just come from  rounding errors.");
                    }

                    // rounding correction
                    if (childProbability > 1 && (childProbability < (1 + PROBABILITY_ROUNDING_EPSILON))) {
                        childProbability = 1.0;
                    }

                    childClade.setProbability(childProbability);

                    if (visitCountMap.containsKey(childClade)) {
                        int childCount = visitCountMap.get(childClade) + 1;
                        visitCountMap.put(childClade, childCount);
                    } else {
                        visitCountMap.put(childClade, 1);
                    }
                    queue.add(childClade);
                }
            }

            visitCountMap.put(clade, -1);
        }

    }

    /** Compute the sum of clade credibilities of all clades in this distribution. */
    public void computeCladeSumCladeCredibilities() {
        resetCacheIfProbabilitiesDirty();
        resetSumCladeCredibilities();
        rootClade.computeSumCladeCredibilities();
    }

    /**
     * Resets the sum of clade credibilities of all clades in this distribution,
     * i.e. to -1 for all non-leaf, non-cherry clades.
     */
    public void resetSumCladeCredibilities() {
        for (Clade clade : this.getClades()) {
            clade.resetSumCladeCredibilities();
        }
    }


    /* -- CREDIBLE SET - CREDIBLE SET -- */

    /** Whether credible level is computed based on clade or clade partition credible levels. */
    public static enum CredibleLevelType {
        Clade, CladePartition;
    }

    /**
     * Credible level information based on clades.
     * Stored value for a clade represents its credible level,
     * i.e., the smallest threshold for which the clade is in the credible set
     * and for the next smaller one it is not.
     */
    private Map<Clade, Double> cladeMinCredibility;

    /**
     * Credible level information based on clade partitions.
     * Stored value for a partition represents its credible level,
     * i.e., the smallest threshold for which the partition is in the credible set
     * and for the next smaller one it is not.
     */
    private Map<CladePartition, Double> partitionMinCredibility;

    public Map<Clade, Double> getCladeMinCredibility() {
        return cladeMinCredibility;
    }

    public void setCladeMinCredibility(Map<Clade, Double> cladeMinCredibility) {
        this.cladeMinCredibility = cladeMinCredibility;
    }

    public Map<CladePartition, Double> getPartitionMinCredibility() {
        return partitionMinCredibility;
    }

    /**
     * @param type of credible set information
     * @return whether the credible set information has been computed and set
     */
    public boolean isCredibleSetInformationInitialized(CredibleLevelType type) {
        return switch (CredibleLevelType.Clade) {
            case Clade -> cladeMinCredibility != null;
            case CladePartition -> partitionMinCredibility != null;
            default -> false;
        };
    }

    public void setPartitionMinCredibility(Map<CladePartition, Double> partitionMinCredibility) {
        this.partitionMinCredibility = partitionMinCredibility;
    }

    /**
     * Return the min credible level of the given tree, that is,
     * the smallest credible set that still contains this tree.
     *
     * @param tree whose cred level is requested
     * @param type what type of information to use
     * @return min credible level of given tree
     */
    public double getCredibleLevelOfTree(Tree tree, CredibleLevelType type) {
        return switch (type) {
            case Clade -> getCredibleLevelOfTreeBasedOnClades(tree);
            case CladePartition -> getCredibleLevelOfTreeBasedOnPartitions(tree);
            default -> 0;
        };
    }

    /* Helper method */
    private double getCredibleLevelOfTreeBasedOnClades(Tree tree) {
        WrappedBeastTree wrapped = new WrappedBeastTree(tree);
        double maxCredLevel = 0;
        for (BitSet bits : wrapped.getClades()) {
            double credLevel = getCredibleLevelOfClade(bits);
            if (credLevel == 0) {
                return credLevel;
            }
            if (credLevel > maxCredLevel) {
                maxCredLevel = credLevel;
            }
        }

        return maxCredLevel;
    }

    /* Helper method */
    private double getCredibleLevelOfClade(BitSet cladeInBits) {
        Clade clade = this.getClade(cladeInBits);
        return (clade == null) ? 0 : cladeMinCredibility.get(clade);
    }

    /* Helper method */
    private double getCredibleLevelOfTreeBasedOnPartitions(Tree tree) {
        double[] credLevel = new double[]{Double.NEGATIVE_INFINITY};
        getCredibleLevelOfPartition(tree.getRoot(), credLevel);
        if (credLevel[0] <= 0) {
            return 0;
        }
        return credLevel[0];
    }

    /* Recursive helper method */
    private BitSet getCredibleLevelOfPartition(Node node, double[] credLevel) {
        BitSet bits;
        if (node.isLeaf()) {
            bits = BitSet.newBitSet(this.getNumberOfLeaves());
            bits.set(node.getNr());
        } else {
            bits = getCredibleLevelOfPartition(node.getLeft(), credLevel);
            if (credLevel[0] == 0) {
                return bits;
            }
            BitSet other = getCredibleLevelOfPartition(node.getRight(), credLevel);
            if (credLevel[0] == 0) {
                return bits;
            }

            Clade leftChild = this.getClade(bits);
            Clade rightChild = this.getClade(other);

            bits.or(other);
            Clade parent = this.getClade(bits);
            if (parent == null) {
                credLevel[0] = 0;
                out.println("clade not found, bits = " + bits);
                out.println("clades in ccd:");
                for (Clade clade : this.getClades()) {
                    out.println("clade = " + clade);
                }
                return bits;
            }

            CladePartition partition = parent.getCladePartition(leftChild, rightChild);
            if (partition == null) {
                credLevel[0] = 0;
            } else {
                credLevel[0] = Math.max(credLevel[0], this.getPartitionMinCredibility().get(partition));
            }
        }

        return bits;
    }

    /**
     * Converts and rounds the given credible level to an integer percentage between 0 and 100.
     * So a value x in [1, 100] represents being (already) in the x%-credible level
     * while a value 0, represents not being in the distribution at all.
     *
     * @param credibleLevel assumed in [0,1]
     * @return credible level converted to int as percentage
     */
    public int convertCredibleLevel(double credibleLevel) {
        if (credibleLevel <= 0) {
            return 0;
        } else {
            return (int) Math.ceil(credibleLevel * 100);
        }
    }


    /*-- DISTANCES - DISTANCES -- */

    /**
     * Returns the average RF distance of the given tree to the trees of this
     * CCD weighted by their probability.
     *
     * @param tree whose average RF distance we compute
     * @return average RF distance of the given tree to this CCD
     */
    public double averageRFDistances(Tree tree) {
        WrappedBeastTree wrappedTree = new WrappedBeastTree(tree);
        HashMap<Clade, Double> cladeRFs = new HashMap<Clade, Double>(this.cladeMapping.size());
        return averageRFDistance(this.rootClade, wrappedTree, cladeRFs);
    }

    /* Recursive helper method */
    private double averageRFDistance(Clade clade, WrappedBeastTree wrappedTree,
                                     HashMap<Clade, Double> cladeRFs) {
        double cladeRF = 0.0;
        if (clade.isLeaf()) {
            cladeRF = 0.0;
        } else {
            for (CladePartition partition : clade.getPartitions()) {
                Clade firstChild = partition.getChildClades()[0];
                Clade secondChild = partition.getChildClades()[1];

                double firstRF = cladeRFs.containsKey(firstChild) ? cladeRFs.get(firstChild)
                        : averageRFDistance(firstChild, wrappedTree, cladeRFs);
                double secondRF = cladeRFs.containsKey(secondChild) ? cladeRFs.get(secondChild)
                        : averageRFDistance(secondChild, wrappedTree, cladeRFs);

                cladeRF += partition.getCCP() * (firstRF + secondRF);
            }

            cladeRF += wrappedTree.containsClade(clade.getCladeInBits()) ? 0 : 1;
        }

        cladeRFs.put(clade, cladeRF);
        return cladeRF;
    }

    public double lostProbability(Set<Clade> excludedClades) {
        HashSet<Clade> handledClades = new HashSet<>(this.getNumberOfClades());

        return lostProbability(this.getRootClade(), excludedClades, handledClades);
    }

    private double lostProbability(Clade clade, Set<Clade> excludedClades, Set<Clade> handledClades) {
        if (clade.isLeaf() || clade.isCherry() || handledClades.contains(clade)) {
            return 0.0;
        }

        double lostProbability = 0.0;
        handledClades.add(clade);
//        out.println("o");

        for (CladePartition partition : clade.getPartitions()) {
            Clade firstChild = partition.getChildClades()[0];
            Clade secondChild = partition.getChildClades()[1];

            if (excludedClades.contains(firstChild) || excludedClades.contains(secondChild)) {
//                out.print("x");
                lostProbability += partition.getCCP();
            } else {
//                out.print("z");
                lostProbability += lostProbability(firstChild, excludedClades, handledClades)
                        + lostProbability(secondChild, excludedClades, handledClades);
            }
        }
        return lostProbability;
    }


    /* -- OTHER METHODS -- */

    /**
     * Create a (deep) copy of this CCD, so with copies of the Clades and
     * CladePartitions; copies at most one stored tree.
     *
     * @return a (deep) copy of this CCD
     */
    public abstract AbstractCCD copy();

    /* General copy helper method for all types of CCD */
    protected static void buildCopy(AbstractCCD original, AbstractCCD copy) {
        for (Clade originalClade : original.getClades()) {
            Clade copiedClade = originalClade.copy(copy);
            copy.cladeMapping.put(originalClade.getCladeInBits(), copiedClade);
        }
        copy.rootClade = copy.cladeMapping.get(original.getRootClade().getCladeInBits());

        for (Clade originalClade : original.getClades()) {
            for (CladePartition originalPartition : originalClade.getPartitions()) {
                Clade copiedParent = copy.cladeMapping.get(originalClade.getCladeInBits());
                Clade copiedChildFirst = copy.cladeMapping
                        .get(originalPartition.getChildClades()[0].getCladeInBits());
                Clade copiedChildSecond = copy.cladeMapping
                        .get(originalPartition.getChildClades()[1].getCladeInBits());

                CladePartition copiedPartition = copiedParent.createCladePartition(copiedChildFirst,
                        copiedChildSecond, true);
                if ((originalPartition.getNumberOfOccurrences() <= 0) || originalPartition.isCCPSet()) {
                    copiedPartition.setCCP(originalPartition.getCCP());
                } else {
                    copiedPartition.increaseOccurrenceCountBy(
                            originalPartition.getNumberOfOccurrences(),
                            originalPartition.getMeanOccurredHeight());
                }
            }
        }
    }

    /**
     * Returns the names of the taxa specified by the mask as "{name1, name2,
     * ...}".
     *
     * @param mask specifies which taxa the names are requested for
     * @return names of taxa as "{name1, name2, ...}"
     */
    public String getTaxaNames(BitSet mask) {
        return this.getTaxaNames(mask, ", ");
    }

    /**
     * Returns the names of the taxa specified by the mask
     * concatenated with the given separator.
     *
     * @param mask      specifies which taxa the names are requested for
     * @param separator used between two taxa
     * @return names of taxa concatenated with given separator
     */
    public String getTaxaNames(BitSet mask, String separator) {
        Tree tree = this.getSomeBaseTree();
        StringBuilder taxa = new StringBuilder("{");
        for (int j = mask.nextSetBit(0); j >= 0; j = mask.nextSetBit(j + 1)) {
            taxa.append(tree.getNode(j).getID()).append(separator);
        }
        return taxa.substring(0, taxa.length() - separator.length()) + "}";
    }

    /**
     * Returns a set of the names of the taxa specified by the mask.
     *
     * @param mask specifies which taxa the names are requested for
     * @return set of names of taxa
     */
    public Set<String> getTaxaNamesList(BitSet mask) {
        Tree tree = this.getSomeBaseTree();
        Set<String> names = new HashSet<>(mask.cardinality());
        int index = 0;
        for (int j = mask.nextSetBit(0); j >= 0; j = mask.nextSetBit(j + 1)) {
            names.add(tree.getNode(j).getID());
        }
        return names;
    }

    @Override
    public String toString() {
        return "[number of leaves: " + this.getNumberOfLeaves() + ", number of clades: "
                + this.getNumberOfClades() + ", max probability: " + this.getMaxTreeProbability()
                + ", entropy: " + this.getEntropy() + ", taxa: " + this.getTaxaAsBitSet() + "]";
    }

    public abstract void initialize();

}
