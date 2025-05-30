package ccd.model;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeUtils;
import beastfx.app.treeannotator.TreeAnnotator.TreeSet;
import ccd.algorithms.TreeDistances;
import ccd.algorithms.credibleSets.CredibleSetType;
import ccd.algorithms.credibleSets.ProbabilityBasedCredibleSetComputer;
import ccd.model.bitsets.BitSet;

import java.io.IOException;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

/**
 * This class represent a most simple distribution of trees, namely, it is
 * simply based on a set of given trees. The probability is determined by the
 * number of occurrences.
 *
 * @author Jonathan Klawitter
 */
public class SampleDistribution implements ITreeDistribution {

    /** Set of trees this distribution is based on. */
    private TreeSet[] treeSets;

    /** Different trees (topologies) with counts of this distribution. */
    private List<WrappedBeastTree> trees = new ArrayList<>(100);

    private HashMap<String, WrappedBeastTree> treeMap = new HashMap<>(100);

    /** The number of trees this CCD is based on. */
    protected int numBaseTrees = 0;

    /** Number of leaves/taxa the trees this CCD is based on. */
    private int numLeaves;

    /**
     * Clades appearing in the trees of this distribution with mapping from a
     * BitSet representation to a clade. Assumes that all BitSets have length
     * equal to the number of leaves.
     */
    private Map<BitSet, Clade> cladeMapping;

    /** Random used by this distribution to sample trees. */
    private Random random;

    /** Whether distribution is tidy (trees sorted, clades extracted). */
    private boolean isDirty = false;

    /**
     * Constructor using one TreeSet that is also used as storage of base trees.
     *
     * @param treeSet burnin assumed to be set; all of its trees are used
     */
    public SampleDistribution(TreeSet treeSet) {
        this(treeSet, treeSet.totalTrees - treeSet.burninCount);
    }

    /**
     * Constructor using one TreeSet that is also used as storage of base trees.
     *
     * @param treeSet       burnin assumed to be set
     * @param numTreesToUse the number of trees to use from the treeSet
     */
    public SampleDistribution(TreeSet treeSet, int numTreesToUse) {
        this.treeSets = new TreeSet[1];
        this.treeSets[0] = treeSet;

        // get all topologies with counts
        initializeTrees(treeSet, numTreesToUse);

        // get all clades with counts
        initializeClades();

        // sort them in descending order
        sortTrees();
    }

    /**
     * Constructor using multiple TreeSets that are also used as storage of base
     * trees using all trees.
     *
     * @param treeSets burnins assumed to be set
     */
    public SampleDistribution(TreeSet[] treeSets) {
        this.treeSets = treeSets;

        // get all topologies with counts
        int sum = 0;
        for (TreeSet treeSet : treeSets) {
            int numTrees = treeSet.totalTrees - treeSet.burninCount;
            initializeTrees(treeSet, numTrees);
            sum += numTrees;
        }
        this.numBaseTrees = sum;

        initialize();
    }

    /**
     * Constructor using multiple TreeSets that are also used as storage of base
     * trees.
     *
     * @param treeSets burnins assumed to be set
     */
    public SampleDistribution(TreeSet[] treeSets, int[] numTreesToUse) {
        this.treeSets = treeSets;

        // get all topologies with counts
        int sum = 0;
        int i = 0;
        for (TreeSet treeSet : treeSets) {
            int numTrees = Math.min(treeSet.totalTrees - treeSet.burninCount, numTreesToUse[i++]);
            initializeTrees(treeSet, numTrees);
            sum += numTrees;
        }
        this.numBaseTrees = sum;

        initialize();
    }

    /**
     * Constructor for empty {@link SampleDistribution} where trees can then be
     * added one by one; these base trees are not stored.
     *
     * @param numberOfLeaves the number of leaves the trees of this distribution have
     */
    public SampleDistribution(int numberOfLeaves) {
        this.numLeaves = numberOfLeaves;
    }

    /* Constructor helper method */
    private void initializeTrees(TreeSet treeSet, int numTreesToUse) {
        int numUsedTrees = 0;
        try {
            treeSet.reset();
            Tree tree = treeSet.next();
            this.numLeaves = tree.getLeafNodeCount();

            while ((tree != null) && (numUsedTrees < numTreesToUse)) {
                numUsedTrees++;
                this.numBaseTrees++;
                processTree(tree);
                tree = treeSet.next();
            }
        } catch (IOException e) {
            System.err.println("Error reading in trees to create distribution.");
        }
    }

    /**
     * Add and process the given tree into this distribution.
     *
     * @param tree to be added and processed into this distribution
     */
    public void addTree(Tree tree) {
        this.numBaseTrees++;
        this.processTree(tree);
        this.isDirty = true;
    }

    /* Constructor helper method */
    private void processTree(Tree tree) {
        String newick = TreeUtils.sortedNewickTopology(tree.getRoot(), true).trim();

        if (treeMap.containsKey(newick)) {
            treeMap.get(newick).increaseCount();
        } else {
            // tree was not equal to any existing tree
            WrappedBeastTree wrappedTree = new WrappedBeastTree(tree);
            treeMap.put(newick, wrappedTree);
            trees.add(wrappedTree);
        }
    }

    /** Set up this SampleDistribution after adding trees. */
    public void initialize() {
        initializeClades();
        sortTrees();
        isDirty = false;
    }

    /* Constructor helper method */
    private void initializeClades() {
        this.cladeMapping = new HashMap<BitSet, Clade>();
        for (WrappedBeastTree wrappedTree : trees) {
            for (BitSet cladeInBits : wrappedTree.getClades()) {
                Clade clade = cladeMapping.get(cladeInBits);
                if (clade == null) {
                    clade = new Clade(cladeInBits, null);
                    cladeMapping.put(cladeInBits, clade);
                }
                clade.increaseOccurrenceCountBy(wrappedTree.getCount());
            }
        }
    }

    /* Constructor helper method */
    private void sortTrees() {
        // decreasing order
        trees.sort((x, y) -> Integer.compare(y.getCount(), x.getCount()));
    }

    /** Tidy up distribution after trees have been added one by one. */
    public void tidyUpIfDirty() {
        if (isDirty) {
            initialize();
        }
    }

    /**
     * Computes and returns the credible set of this distribution for the given credible level.
     *
     * @param credibleLevel how much probability to keep (e.g. 0.95)
     * @return a SampleDistribution that represents the credible set of this one
     */
    public SampleDistribution getCredibilitySet(final double credibleLevel) {
        double sum = 0;
        double numTrees = this.numBaseTrees;
        int totalCount = 0;
        SampleDistribution credi = new SampleDistribution(this.numLeaves);

        for (WrappedBeastTree tree : trees) {
            credi.trees.add(tree);
            totalCount += tree.getCount();

            sum += tree.getCount() / numTrees;
            if (sum >= credibleLevel) {
                break;
            }
        }

        credi.initializeClades();
        credi.numBaseTrees = totalCount;

        return credi;
    }

    /** @return trees of this distribution, sorted by decreasing probability */
    public List<WrappedBeastTree> getTrees() {
        return trees;
    }

    @Override
    public Tree sampleTree() {
        return sampleTree(HeightSettingStrategy.None);
    }

    @Override
    public Tree sampleTree(HeightSettingStrategy heightStrategy) {
        this.tidyUpIfDirty();

        int pick = this.getRandom().nextInt(numBaseTrees);
        int next = 0;
        for (WrappedBeastTree tree : trees) {
            next += tree.getCount();

            if (pick <= next) {
                return tree.getWrappedTree();
            }
        }

        System.err.println("No tree sampled. Suspect number of tree counts does not add up to number of base trees.");
        return null;
    }

    @Override
    public double sampleTreeProbability() {
        this.tidyUpIfDirty();

        int pick = this.getRandom().nextInt(numBaseTrees);
        int next = 0;
        for (WrappedBeastTree tree : trees) {
            next += tree.getCount();

            if (pick <= next) {
                return tree.getCount() / (double) this.numBaseTrees;
            }
        }

        System.err.println("No tree probability sampled. Suspect number of tree counts does not add up to number of base trees.");
        return -1;
    }

    @Override
    public int getNumberOfLeaves() {
        return numLeaves;
    }

    @Override
    public int getNumberOfClades() {
        return cladeMapping.size();
    }

    @Override
    public Collection<Clade> getClades() {
        return cladeMapping.values();
    }

    @Override
    public BigInteger getNumberOfTrees() {
        return BigInteger.valueOf(trees.size());
    }

    /** @return the number of trees in this distribution (as int) */
    public int getNumberOfTreesInt() {
        return trees.size();
    }

    /** @return number of trees this distribution is based on */
    public int getNumberOfBaseTrees() {
        return numBaseTrees;
    }

    @Override
    public Tree getMAPTree() {
        return getMAPTree(HeightSettingStrategy.None);
    }

    @Override
    public Tree getMAPTree(HeightSettingStrategy heightStrategy) {
        this.tidyUpIfDirty();

        Tree tree = this.trees.get(0).getWrappedTree();

        setHeights(tree, heightStrategy);

        return tree;
    }

    /* Helper method */
    private void setHeights(Tree tree, HeightSettingStrategy heightStrategy) {
        switch (heightStrategy) {
            case None:
                return;
            case One:
                setHeightsToOne(tree.getRoot());
                return;
            default:
                // TODO
                System.err.println(
                        "Height setting strategy (" + heightStrategy + ") not implemented yet.");
                break;
        }
    }

    /* Recursive helper method */
    private void setHeightsToOne(Node vertex) {
        if (vertex.isLeaf()) {
            vertex.setHeight(0);
        } else {
            Node firstChild = vertex.getChild(0);
            setHeightsToOne(firstChild);

            Node secondChild = vertex.getChild(1);
            setHeightsToOne(secondChild);

            vertex.setHeight(Math.max(firstChild.getHeight(), secondChild.getHeight()) + 1);
        }
    }

    @Override
    public double getMaxTreeProbability() {
        this.tidyUpIfDirty();

        return this.trees.get(0).getCount() / (double) this.numBaseTrees;
    }

    /** @return the MCC (maximum clade credibility tree) of this distribution */
    public Tree getMCCTree() {
        this.tidyUpIfDirty();

        Tree bestTree = null;
        double maxCC = 0;
        for (WrappedBeastTree tree : trees) {
            double cladeCredibility = 1;
            for (BitSet cladeInBits : tree.getNontrivialClades()) {
                Clade clade = cladeMapping.get(cladeInBits);
                if (!clade.isLeaf()) {
                    cladeCredibility *= clade.getNumberOfOccurrences() / (double) this.numBaseTrees;
                }
            }

            if (cladeCredibility > maxCC) {
                maxCC = cladeCredibility;
                bestTree = tree.getWrappedTree();
            }
        }

        return bestTree;
    }

    @Override
    public double getProbabilityOfTree(Tree tree) {
        this.tidyUpIfDirty();

        WrappedBeastTree wrappedTree = new WrappedBeastTree(tree);
        for (WrappedBeastTree wrappedOther : trees) {
            if (TreeDistances.robinsonsFouldDistance(wrappedOther, wrappedTree) == 0) {
                return wrappedOther.getCount() / (double) this.numBaseTrees;
            }
        }
        return 0;
    }

    @Override
    public boolean containsTree(Tree tree) {
        WrappedBeastTree wrappedTree = new WrappedBeastTree(tree);
        for (WrappedBeastTree wrappedOther : trees) {
            if (TreeDistances.robinsonsFouldDistance(wrappedOther, wrappedTree) == 0) {
                return true;
            }
        }
        return false;
    }

    @Override
    public double getCladeProbability(BitSet cladeInBits) {
        Clade clade = cladeMapping.get(cladeInBits);
        if (clade == null) {
            return 0;
        }

        return clade.getNumberOfOccurrences() / (double) this.getNumberOfBaseTrees();
    }

    /** @return the entropy of this tree distribution */
    public double getEntropy() {
        double entropy = 0;

        for (WrappedBeastTree tree : trees) {
            double count = tree.getCount();
            double p = count / numBaseTrees;
            entropy -= p * Math.log(p);
        }

        return entropy;
    }


    /* -- CREDIBLE SET - CREDIBLE SET -- */

    private ProbabilityBasedCredibleSetComputer treeBasedCredibleSetComputer;

    @Override
    public double getCredibleLevel(Tree tree, CredibleSetType type) {
        if (type == CredibleSetType.Frequency) {
            WrappedBeastTree wrappedInputTree = new WrappedBeastTree(tree);

            double sum = 0;
            for (WrappedBeastTree wrappedStoredTree : trees) {
                sum += wrappedStoredTree.getCount();

                if (TreeDistances.robinsonsFouldDistance(wrappedStoredTree, wrappedInputTree) == 0) {
                    return sum / this.numBaseTrees;
                }
            }

            return -1;
        } else if (type == CredibleSetType.TreeSampling) {
            if (treeBasedCredibleSetComputer == null) {
                treeBasedCredibleSetComputer = new ProbabilityBasedCredibleSetComputer(this);
            }
            return treeBasedCredibleSetComputer.getCredibleLevel(tree);
        } else {
            throw new IllegalArgumentException("Credible set type not supported by sample distributions.");
        }
    }


    /* -- OTHER - OTHER -- */

    /**
     * Computes and returns the Fair Proportion Diversity Index of the taxa in this sample distribution.
     *
     * @return Fair Proportion Diversity Index of the taxa in this sample distribution
     */
    public double[] getFairProportionIndex() {
        int numLeaves = trees.get(0).getWrappedTree().getLeafNodeCount();
        double[] index = new double[numLeaves];
        double numTrees = (double) this.numBaseTrees;

        for (WrappedBeastTree wrappedTree : trees) {
            Tree tree = wrappedTree.getWrappedTree();
            double[] currentIndex = getFairProportionIndex(tree);
            double treeWeight = wrappedTree.getCount() / numTrees;

            for (int i = 0; i < index.length; i++) {
                index[i] += currentIndex[i] * treeWeight;
            }
        }

        return index;
    }

    /**
     * Computes and returns the Fair Proportion Diversity Index of the taxa for the given tree.
     *
     * @param tree for which FP is computed
     * @return the FP values of the taxa for the given tree
     */
    public static double[] getFairProportionIndex(Tree tree) {
        int numLeaves = tree.getLeafNodeCount();
        double[] index = new double[numLeaves];

        getFairProportionIndex(tree.getRoot(), numLeaves, index);

        return index;
    }

    /* Recursive helper method */
    private static BitSet getFairProportionIndex(Node vertex, int numLeaves, double[] index) {
        if (vertex.isLeaf()) {
            BitSet bitset = BitSet.newBitSet(numLeaves);
            bitset.set(vertex.getNr());
            return bitset;
        }

        BitSet leftBits = getFairProportionIndex(vertex.getLeft(), numLeaves, index);
        double branchLength = vertex.getHeight() - vertex.getLeft().getHeight();
        int size = leftBits.cardinality();
        double diversity = branchLength / size;

        if (branchLength < 0) {
            System.out.println("branchLength = " + branchLength);
            System.out.println("size = " + size);
            System.out.println("diversity = " + diversity);
            System.out.println("parent = " + vertex);
            System.out.println("childL = " + vertex.getLeft());
            throw new AssertionError("Negative branch length.");
        }

        int i = leftBits.nextSetBit(0);
        while (i != -1) {
            index[i] += diversity;
            i = leftBits.nextSetBit(i + 1);
        }

        BitSet rightBits = getFairProportionIndex(vertex.getRight(), numLeaves, index);
        branchLength = vertex.getHeight() - vertex.getRight().getHeight();
        size = rightBits.cardinality();
        diversity = branchLength / size;

        i = rightBits.nextSetBit(0);
        while (i != -1) {
            index[i] += diversity;
            i = rightBits.nextSetBit(i + 1);
        }

        if (branchLength < 0) {
            System.out.println("branchLength = " + branchLength);
            System.out.println("size = " + size);
            System.out.println("diversity = " + diversity);
            System.out.println("parent = " + vertex);
            System.out.println("childR = " + vertex.getRight());
        }


        leftBits.or(rightBits);
        return leftBits;
    }

    /** Set the Random of this distribution. */
    public void setRandom(Random random) {
        this.random = random;
    }

    /** @return the Random of this distribution */
    public Random getRandom() {
        if (this.random == null) {
            this.random = new Random();
        }

        return random;
    }
}
