package ccd.algorithms;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import ccd.model.AbstractCCD;
import ccd.model.BitSet;
import ccd.model.Clade;
import ccd.model.HeightSettingStrategy;

import java.util.ArrayList;
import java.util.List;

/**
 * This class provides a method to compute the greedy majority rule consensus
 * (MRC) tree.
 * <p>
 * Currently uses an algorithm that is not asymptotically optimal.
 *
 * @author Jonathan Klawitter
 */
public class GreedyMRCTree {

    /**
     * Computes the greedy majority rule consensus (MRC) tree based on the
     * clades in the given CCD; uses registered occurrences as clade counts.
     * Tides are broken arbitrarily but deterministically.
     *
     * @param ccd            based on whose clades and clade counts the MRC tree is
     *                       constructed
     * @param heightStrategy
     * @return greedy MRC tree based on clades of given CCD
     */
    public static Tree constructTree(AbstractCCD ccd, HeightSettingStrategy heightStrategy) {
        // we need all clades with counts; here counts are stored in the clades
        // as number of occurrences
        ArrayList<Clade> clades = new ArrayList<>(ccd.getNumberOfClades());
        clades.addAll(ccd.getClades());
        // System.out.println("num base clades: " + clades.size());

        // we then first find a set of compatible clades, that will form the tree
        ArrayList<Clade> compatibleClades = computeCompatibleClades(clades);
        // System.out.println("num comp. clades: " + compatibleClades.size());

        // sort clades by increasing size, to build tree bottom-up
        List<Clade> tempClades = compatibleClades.stream()
                .sorted((x, y) -> Integer.compare(x.size(), y.size())).toList();
        ArrayList<Clade> vertexClades = new ArrayList<>(tempClades.size() + 1);
        vertexClades.addAll(tempClades);
        vertexClades.add(ccd.getRootClade());

        // construct tree bottom-up using an array where we only store the
        // current vertices of constructed subtrees at lowest index of their leaves;
        // all other vertices are taken out
        int n = ccd.getNumberOfLeaves();
        Node[] vertices = new Node[n];

        // create leaves
        Tree baseTree = ccd.getSomeBaseTree();
        for (int i = 0; i < n; i++) {
            String taxonName = baseTree.getNode(i).getID();
            Node vertex = new Node(taxonName);
            vertex.setNr(i);
            if (heightStrategy == HeightSettingStrategy.One) {
                vertex.setHeight(0);
            }

            // each leaf is stored at its index
            vertices[i] = vertex;
        }

        // create inner vertices
        int runningIndex = n;
        for (Clade clade : vertexClades) {
            BitSet cladeInBits = clade.getCladeInBits();
            int firstIndex = cladeInBits.nextSetBit(0);
            Node firstChild = vertices[firstIndex];

            int secondIndex = firstIndex;
            Node secondChild = null;
            do {
                secondIndex = cladeInBits.nextSetBit(secondIndex + 1);
                secondChild = vertices[secondIndex];
            } while (secondChild == null);

            String id = runningIndex + "";
            Node vertex = new Node(id);
            vertex.setNr(runningIndex++);
            vertex.addChild(firstChild);
            vertex.addChild(secondChild);

            // update vertex array
            vertices[firstIndex] = vertex;
            vertices[secondIndex] = null;

            // if the sizes of the children do not add up to clade size, then
            // the tree is not fully resolved and we have to continue searching
            int childSum = firstChild.getLeafNodeCount() + secondChild.getLeafNodeCount();
			/*-if (childSum < clade.size()) {
				System.out.println("parent clade size: " + clade.size() + ", " + clade);
				System.out.println("fc, index: " + firstIndex + ", clade size: "
						+ firstChild.getLeafNodeCount());
				System.out.println("sc, index: " + secondIndex + ", clade size: "
						+ secondChild.getLeafNodeCount());
			}*/
            while (childSum < clade.size()) {
                int nextIndex = secondIndex;
                Node nextChild = null;
                do {
                    nextIndex = cladeInBits.nextSetBit(nextIndex + 1);
                    nextChild = vertices[nextIndex];
                } while (nextChild == null);
                vertex.addChild(nextChild);
                childSum += nextChild.getLeafNodeCount();

                vertices[nextIndex] = null;
            }

            if (heightStrategy == HeightSettingStrategy.One) {
                double height = Math.max(firstChild.getHeight(), secondChild.getHeight()) + 1;
                vertex.setHeight(height);
            }

            // System.out.println("");
        }

        Tree tree = new Tree(vertices[0]);

        return tree;
    }

    /* Helper method - uses pairwise comparison to build compatible set */
    private static ArrayList<Clade> computeCompatibleClades(List<Clade> clades) {
        // 1. we sort them in descending order of counts
        List<Clade> sortedClades = clades.stream().filter(x -> !x.isLeaf()).filter(x -> !x.isRoot())
                .sorted((x, y) -> Integer.compare(y.getNumberOfOccurrences(),
                        x.getNumberOfOccurrences())) // sorted decreasingly
                .toList();

        // 2. build up compatible clade set by pairwise compatibility check
        ArrayList<Clade> compatibleClades = new ArrayList<>();
        for (Clade testedClade : sortedClades) {
            boolean compatible = true;
            for (Clade clade : compatibleClades) {
                // two clades are compatible if one is contained in the other or they are disjoint;
                // so they are not compatible, if they intersect properly
                if (BitSetUtil.intersectProperly(clade.getCladeInBits(),
                        testedClade.getCladeInBits())) {
                    compatible = false;
                    break;
                }
            }
            if (compatible) {
                compatibleClades.add(testedClade);
            }
        }
        return compatibleClades;
    }

}
