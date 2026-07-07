package ccd.algorithms;

import beast.base.evolution.tree.Node;
import ccd.model.WrappedBeastTree;
import ccd.model.WrappedBeastTreeWithSampledAncestor;
import ccd.model.bitsets.BitSet;

import java.util.*;

/**
 * This class provides distance computation methods for two beast trees.
 * {@link ccd.model.WrappedBeastTree} are used to cache values, in particular, BitSets representing the clades.
 */
public class TreeDistances {

    /**
     * Compute the Robinsons-Fould (RF) distance (divided by 2) of the two given
     * trees, that is, the number of clades in the first tree that are not in
     * the second tree. Assumption is that the two trees are binary otherwise
     * the current implementation is not symmetric.
     *
     * @param first  tree, assumed to be binary
     * @param second tree, assumed to be binary
     * @return the RF distance (divided by 2) for the two given trees
     */
    public static int robinsonsFouldDistance(WrappedBeastTree first, WrappedBeastTree second) {
        // if input trees are of type WrappedBeastTreeWithSampledAncestor,
        // clades with the same taxa but different SA info would be identified as different clades,
        // therefore increase the distance
        ArrayList<BitSet> firstClades = first.getNontrivialClades();
        ArrayList<BitSet> secondClades = second.getNontrivialClades();
        firstClades.removeAll(secondClades);
        return firstClades.size();
    }

    /**
     * Symmetric difference of the two trees' sampled-ancestor leaf sets:
     * {@code |SA(first) △ SA(second)|}. This counts only SA-status
     * disagreements per leaf and ignores topology entirely; two trees with
     * completely different topologies but identical SA sets get distance 0.
     *
     * <p>This is <em>not</em> the sampled-ancestor Robinson-Foulds distance;
     * it is one of its two components (the asymmetric singleton penalty).
     * For the full SA-RF metric that combines topology and SA-status
     * differences, see
     * {@link #sampledAncestorRobinsonFoulds(WrappedBeastTreeWithSampledAncestor, WrappedBeastTreeWithSampledAncestor)}.
     */
    public static int sampledAncestorSymmetricDistance(WrappedBeastTreeWithSampledAncestor first, WrappedBeastTreeWithSampledAncestor second) {
        ArrayList<Node> firstTreeSANodes = new ArrayList<>(first.getSampledAncestors());
        ArrayList<Node> secondTreeSANodes = new ArrayList<>(second.getSampledAncestors());

        ArrayList<Node> onlyFirst = new ArrayList<>(firstTreeSANodes);
        onlyFirst.removeAll(secondTreeSANodes);

        ArrayList<Node> onlySecond = new ArrayList<>(secondTreeSANodes);
        onlySecond.removeAll(firstTreeSANodes);

        return onlyFirst.size() + onlySecond.size();
    }

    /**
     * Sampled-ancestor Robinson-Foulds distance between two SA trees on the
     * same taxa. The "cuttable" clade set of a tree {@code T = (tau, sigma)}
     * is
     * <pre>
     *   C_cut(T) = { C(v) : v internal in tau } u { {i} : sigma_i = 0 }
     * </pre>
     * the taxa-only internal-node clades, plus singleton clades for the
     * non-SA leaves only (an SA leaf has a zero-length parent branch and
     * cannot be separated from its parent's clade by an edge cut, so no
     * {@code {i}} singleton is contributed). The distance is the size of
     * the symmetric difference,
     * <pre>
     *   d_SA-RF(T, T') = |C_cut(T) △ C_cut(T')|.
     * </pre>
     *
     * <p>Properties:
     * <ul>
     *   <li>Reduces to standard rooted RF (full symmetric difference) when
     *       neither tree contains SAs.
     *   <li>Asymmetric singleton penalty: leaf {@code i} being SA in only
     *       one tree contributes 1 (not 2) to the distance, since the SA
     *       tree contributes no {@code {i}} clade.
     *   <li>Symmetric and non-negative; zero iff the two trees have the
     *       same internal clades and the same SA pattern.
     * </ul>
     *
     * <p>This method returns the full symmetric-difference cardinality,
     * unlike {@link #robinsonsFouldDistance(WrappedBeastTree, WrappedBeastTree)}
     * which returns the halved count and relies on the binary-tree symmetry
     * {@code |T \\ T'| = |T' \\ T|} that the asymmetric singleton penalty
     * here breaks.
     */
    public static int sampledAncestorRobinsonFoulds(WrappedBeastTreeWithSampledAncestor first,
                                                    WrappedBeastTreeWithSampledAncestor second) {
        int n = first.getWrappedTree().getLeafNodeCount();
        if (second.getWrappedTree().getLeafNodeCount() != n) {
            throw new IllegalArgumentException(
                    "SA-RF requires both trees to have the same number of taxa; got "
                            + n + " and " + second.getWrappedTree().getLeafNodeCount());
        }

        Set<BitSet> cuttableFirst = cuttableCladeSet(first, n);
        Set<BitSet> cuttableSecond = cuttableCladeSet(second, n);

        int onlyFirst = 0;
        for (BitSet c : cuttableFirst) if (!cuttableSecond.contains(c)) onlyFirst++;
        int onlySecond = 0;
        for (BitSet c : cuttableSecond) if (!cuttableFirst.contains(c)) onlySecond++;

        return onlyFirst + onlySecond;
    }

    /**
     * Pure topology Robinson-Foulds distance on SA trees: counts internal-node
     * clades (taxa-only, ignoring SA flags) that appear in {@code first} but
     * not {@code second}.
     *
     * <p>Useful as a complement to
     * {@link #robinsonsFouldDistance(WrappedBeastTree, WrappedBeastTree)} when
     * called on {@link WrappedBeastTreeWithSampledAncestor} instances: that
     * method uses the wide bitsets ({@code [0,n)} taxa, {@code [n,2n)} SA
     * flags) so two trees with the same taxa-only clades but different
     * direct-SA-leaf-children patterns at a clade are counted as differing.
     * This method strips the SA flags and considers only base topology.
     *
     * <p>Follows the halved convention of
     * {@link #robinsonsFouldDistance(WrappedBeastTree, WrappedBeastTree)}: for
     * binary trees on the same taxa,
     * {@code |topology(T1) \\ topology(T2)| = |topology(T2) \\ topology(T1)|},
     * so this is half the symmetric difference.
     */
    public static int topologyRobinsonFoulds(WrappedBeastTreeWithSampledAncestor first,
                                             WrappedBeastTreeWithSampledAncestor second) {
        int n = first.getWrappedTree().getLeafNodeCount();
        if (second.getWrappedTree().getLeafNodeCount() != n) {
            throw new IllegalArgumentException(
                    "topologyRF requires both trees to have the same number of taxa; got "
                            + n + " and " + second.getWrappedTree().getLeafNodeCount());
        }
        Set<BitSet> firstClades = internalCladesTaxaOnly(first, n);
        Set<BitSet> secondClades = internalCladesTaxaOnly(second, n);
        int onlyFirst = 0;
        for (BitSet c : firstClades) if (!secondClades.contains(c)) onlyFirst++;
        return onlyFirst;
    }

    private static Set<BitSet> cuttableCladeSet(WrappedBeastTreeWithSampledAncestor t, int n) {
        Set<BitSet> set = internalCladesTaxaOnly(t, n);
        for (Node leaf : t.getWrappedTree().getExternalNodes()) {
            if (leaf.getLength() != 0.0) {
                BitSet singleton = BitSet.newBitSet(n);
                singleton.set(leaf.getNr());
                set.add(singleton);
            }
        }
        return set;
    }

    private static Set<BitSet> internalCladesTaxaOnly(WrappedBeastTreeWithSampledAncestor t, int n) {
        Set<BitSet> set = new HashSet<>();
        collectInternalCladesTaxaOnly(t.getWrappedTree().getRoot(), n, set);
        return set;
    }

    private static BitSet collectInternalCladesTaxaOnly(Node v, int n, Set<BitSet> internalClades) {
        if (v.isLeaf()) {
            BitSet leafBits = BitSet.newBitSet(n);
            leafBits.set(v.getNr());
            return leafBits;
        }
        BitSet cladeBits = BitSet.newBitSet(n);
        for (Node child : v.getChildren()) {
            cladeBits.or(collectInternalCladesTaxaOnly(child, n, internalClades));
        }
        internalClades.add(cladeBits);
        return cladeBits;
    }

    /**
     * Compute the path distance of the two given trees (by
     * <a href="https://doi.org/10.1093/sysbio/42.2.126">Steel and Penny,
     * 1993</a>), defines as follows. Define matrix with d[i][j] the length of
     * path between leaf i and leaf j. Then path distance is square-root of sum
     * of the squares of the differences d[i][j](T) - d[i][j](T').
     *
     * @param first
     * @param second
     * @return the path distance of the two given trees
     */
    public static double pathDistance(WrappedBeastTree first, WrappedBeastTree second) {
        double[][] firstMatrix = first.getPathDistanceMatrix();
        double[][] secondMatrix = second.getPathDistanceMatrix();

        double squaredSum = 0;
        for (int i = 0; i < firstMatrix.length; i++) {
            for (int j = i + 1; j < firstMatrix[0].length; j++) {
                double dif = firstMatrix[i][j] - secondMatrix[i][j];
                squaredSum += dif * dif;
            }
        }

        return Math.sqrt(squaredSum);
    }

    /**
     * Compute the matching distance of the two given trees, a refinement of the
     * RF distance; see
     * <a href="https://doi.org/10.2478/amcs-2013-0050">Bogdanowcz and Giaro,
     * 2013</a>: we compute a min-weight pereft matching between clades of the
     * first tree and the second tree where the weight for one pair of clades is
     * given by the size of their symmetric difference.
     *
     * @param first
     * @param second
     * @return the matching distance for the two given trees
     */
    public static int matchingDistance(WrappedBeastTree first, WrappedBeastTree second) {
        ArrayList<BitSet> firstClades = first.getNontrivialClades();
        ArrayList<BitSet> secondClades = second.getNontrivialClades();
        firstClades.removeAll(secondClades);
        secondClades.removeAll(first.getNontrivialClades());

        // remark:
        // why can we take out the clades that appear in both trees?
        // suppose we used all clades when finding a min-weight perfect
        // matching. then it might indeed happen that a clade C appearing in T
        // and T' is not matched up, but instead paired with D in T and E in T'.
        // by design of the cost function, there exists then an alternative
        // solution
        // where C is paired with C and D with E at the same cost.

        // System.out.println("Computing matchingDistance of pair of trees");
        // System.out.println("First tree: " +
        // first.getWrappedTree().toString());
        // System.out.println("Second tree: " +
        // second.getWrappedTree().toString());
        //
        // System.out.println("first clades size: " + firstClades.size());
        // System.out.println("second clades size: " + secondClades.size());

        int n = firstClades.size();
        int[] costMatrix = new int[n * n]; // flattened
        // System.out.println("cost matrix: ");
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                BitSet copy = (BitSet) firstClades.get(i).clone();
                copy.xor(secondClades.get(j));
                costMatrix[i * n + j] = copy.cardinality();
                // System.out.print(costMatrix[i * n + j] + ",");
            }
            // System.out.println("");
        }

        return minCostMatchingInt(costMatrix, n);
    }

    /* WIP */
    public static double phylogeneticInformationDistance(WrappedBeastTree first,
                                                         WrappedBeastTree second) {
        // double spi = sharedPhylogeneticInformation(first, second);
        // TODO
        return 0;
    }

    /* WIP */
    public static double sharedPhylogeneticInformation(WrappedBeastTree first,
                                                       WrappedBeastTree second) {
        ArrayList<BitSet> firstClades = first.getNontrivialClades();
        ArrayList<BitSet> secondClades = second.getNontrivialClades();

        // TODO figure out if we can take out matching clades

        // 1. set up cost matrix
        int n = firstClades.get(0).size();
        int m = firstClades.size();
        double[] costMatrix = new double[m * m]; // flattened

        // System.out.println("cost matrix: ");
        for (int i = 0; i < m; i++) {
            BitSet firstClade = firstClades.get(i);
            // System.out.println("firstClade: " + firstClade);
            for (int j = 0; j < m; j++) {
                double hShared;
                BitSet secondClade = secondClades.get(j);
                // System.out.println(" secondClade: " + secondClade);

                // determine if same, subset/disjoint, or conflicting
                if (firstClade.equals(secondClade)) {
                    // case 1: firstClade and secondClade equal
                    // System.out.println(" case 1: equal");
                    hShared = -Math.log(priorMonophylyProbability(firstClade.cardinality(), n));
                } else if (BitSetUtil.intersectProperly(firstClade, secondClade)) {
                    // case 2: firstClade and secondClade intersect properly
                    // System.out.println(" case 2: conflicting");
                    // so nothing to do here
                    hShared = 0;
                } else {
                    // otherwise
                    double pFirst = priorMonophylyProbability(firstClade.cardinality(), n);
                    double hFirst = -Math.log(pFirst);
                    double pSecond = priorMonophylyProbability(secondClade.cardinality(), n);
                    double hSecond = -Math.log(pSecond);

                    double pBoth; // p(C1, C2) = p(C1 | C2) * p(C2) = p(C2 | C1)
                    // * p(C1)
                    double pOneConditionalOther;
                    if (BitSetUtil.disjoint(firstClade, secondClade)) {
                        // case 3: firstClade and secondClade are disjoint
                        // System.out.println(" case 3: disjoint");
                        pOneConditionalOther = priorMonophylyProbability(firstClade.cardinality(),
                                n - secondClade.cardinality());
                        pBoth = pOneConditionalOther * pSecond;
                    } else {
                        // case 4: one contains the other
                        if (BitSetUtil.contains(secondClade, firstClade)) {
                            // case 4a: firstClade contained in secondClade
                            // System.out.println(" case 4a: first in second");
                            pOneConditionalOther = priorMonophylyProbability(
                                    firstClade.cardinality(), secondClade.cardinality());
                            // System.out.println(" pConditional: " +
                            // pOneConditionalOther);
                            // System.out.println(" pSecond: " + pSecond);
                            pBoth = pOneConditionalOther * pSecond;
                            // System.out.println(" pBoth: " + pBoth);
                        } else {
                            // case 4b: secondClade contained in firstClade
                            // System.out.println(" case 4a: second in first");
                            pOneConditionalOther = priorMonophylyProbability(
                                    secondClade.cardinality(), firstClade.cardinality());
                            pBoth = pOneConditionalOther * pFirst;
                        }
                    }

                    double hBoth = -Math.log(pBoth);
                    // System.out.println(" hFirst: " + hFirst);
                    // System.out.println(" hSecond: " + hSecond);
                    // System.out.println(" hBoth: " + hBoth);
                    hShared = hFirst + hSecond - hBoth;
                }

                costMatrix[i * m + j] = hShared;
                // System.out.print(costMatrix[i * m + j] + ",");
            }
            // System.out.println("");
        }

        // 2. solve max cost matching
        return optCostMatchingDouble(costMatrix, m, false);
    }

    double[][] priorMonophylyProbabilities;

    public static double priorMonophylyProbability(int a, int n) {
        // TODO add caching
        double p = 2;
        p = (a + n) / ((double) a * (a + 1));
        double binom = 1;
        for (int i = 0; i < a; i++) {
            binom *= (n - i) / (double) (a - i);
        }
        return p / binom;
    }

    private static int minCostMatchingInt(int[] costMatrix, int n) {
        // # vertex indices:
        // indices on left side: 0 to n-1
        // (adjusted) indices on right side: n to 2n-1
        int sinkIndex = 2 * n;

        // whether a node on the left has been matched
        // then no augmenting path starts from it
        boolean[] matched = new boolean[2 * n];

        // # edge indices (flattened):
        // index of edge from left i to right j is i*n + j
        boolean[] edgeUsed = new boolean[n * n];

        int totalCost = 0;

        PriorityQueue<Vertistance> queue = new PriorityQueue<Vertistance>(new Vertistance());
        int[] prev = null;
        // int i = 1;
        do {
            // System.out.println("## dijkstra run " + i++ + " ##");

            // settled maintains whether node has been handled yet
            boolean[] settled = new boolean[2 * n];
            // prev points backwards along the augmenting path
            prev = new int[sinkIndex + 1];
            Arrays.fill(prev, -1);
            // distance to that node from source
            int distance[] = new int[sinkIndex + 1];
            Arrays.fill(distance, Integer.MAX_VALUE);

            // # Dijkstra algorithm to find min-cost augmenting path
            // from s to all on left side

            // initialize left side
            for (int j = 0; j < n; j++) {
                // no augmenting path starts at an already matched left node
                if (!matched[j]) {
                    queue.add(new Vertistance(j, 0));
                    distance[j] = 0;
                }
            }
            // back and forth between left and right
            while (!queue.isEmpty()) {
                int currentIndexAdjusted = queue.remove().index;

                if (settled[currentIndexAdjusted]) {
                    continue;
                } else {
                    settled[currentIndexAdjusted] = true;
                }

                boolean onLeftSide = (currentIndexAdjusted < n);
                int currentIndex = onLeftSide ? currentIndexAdjusted : currentIndexAdjusted - n;
                int currentDistance = distance[currentIndexAdjusted];

                // System.out.println("\n# currentVertex: " +
                // currentIndexAdjusted);
                // process all neighbours of current vertex
                for (int neighbourIndex = 0; neighbourIndex < n; neighbourIndex++) {
                    int neighbourIndexAdjusted = onLeftSide ? n + neighbourIndex : neighbourIndex;
                    // System.out.println("neighbour: " +
                    // neighbourIndexAdjusted);

                    // if current neighbor node hasn't already been processed
                    // and from left to right [right to left],
                    // the edge must not [must] be in the matching,
                    // then update distances
                    if ((!settled[neighbourIndexAdjusted]) //
                            && ((onLeftSide && !edgeUsed[currentIndex * n + neighbourIndex])
                            || (!onLeftSide
                            && edgeUsed[neighbourIndex * n + currentIndex]))) {

                        // left-to-right we add the edge cost,
                        // right-to-left we subtract it
                        int edgeDistance = onLeftSide
                                ? costMatrix[currentIndex * n + neighbourIndex]
                                : -costMatrix[neighbourIndex * n + currentIndex];
                        int newDistance = currentDistance + edgeDistance;
                        // System.out.println(
                        // currentDistance + " + " + edgeDistance + " = " +
                        // newDistance);

                        // if we lower the distance for the neighbour,
                        // then update values and add it to queue
                        if (newDistance < distance[neighbourIndexAdjusted]) {
                            distance[neighbourIndexAdjusted] = newDistance;
                            prev[neighbourIndexAdjusted] = currentIndex;
                            queue.add(new Vertistance(neighbourIndexAdjusted, newDistance));

                            // System.out.println("updated neighbour distance: "
                            // + distance[neighbourIndexAdjusted]);
                        }
                    }
                }

                // but don't forget distance to sink
                // (if on right side and not matched yet)
                if (!onLeftSide && !matched[currentIndexAdjusted]
                        && (currentDistance < distance[sinkIndex])) {
                    distance[sinkIndex] = currentDistance;
                    prev[sinkIndex] = currentIndex;
                }
            }

            // int oldTotalCost = totalCost;
            if (prev[sinkIndex] != -1) {
                boolean onLeftSide = false;
                int leftIndex = -1;
                int rightIndex = prev[sinkIndex];
                matched[rightIndex + n] = true;
                // System.out.println("new augmenting path");
                // System.out.print("t - " + rightIndex + "r - ");
                do {
                    if (onLeftSide) {
                        // augmenting path backwards from left to right was over
                        // used edge
                        rightIndex = prev[leftIndex];
                        matched[rightIndex + n] = true;
                        edgeUsed[leftIndex * n + rightIndex] = false;
                        totalCost -= costMatrix[leftIndex * n + rightIndex];

                        // System.out.print(rightIndex + "r - ");
                    } else {
                        // augmenting path backwards from right to left was over
                        // unused edge, which now becomes used
                        leftIndex = prev[rightIndex + n];
                        matched[leftIndex] = true;
                        edgeUsed[leftIndex * n + rightIndex] = true;
                        totalCost += costMatrix[leftIndex * n + rightIndex];

                        // System.out
                        // .println(leftIndex + "l - " + ((prev[leftIndex] !=
                        // -1) ? "" : "s"));
                    }
                    onLeftSide = !onLeftSide;
                } while (prev[leftIndex] != -1);
            }
            // System.out
            // .println("totalCost: " + totalCost + " (" + (totalCost -
            // oldTotalCost) + ")\n");

        } while (prev[sinkIndex] != -1);

        return totalCost;
    }

    private static double optCostMatchingDouble(double[] costMatrix, int n, boolean minimize) {
        // # vertex indices:
        // indices on left side: 0 to n-1
        // (adjusted) indices on right side: n to 2n-1
        int sinkIndex = 2 * n;

        // whether a node on the left has been matched
        // then no augmenting path starts from it
        boolean[] matched = new boolean[2 * n];

        // # edge indices (flattened):
        // index of edge from left i to right j is i*n + j
        boolean[] edgeUsed = new boolean[n * n];

        double totalCost = 0;

        PriorityQueue<Vertistance> queue = new PriorityQueue<Vertistance>(new Vertistance());
        int[] prev = null;
        // int i = 1;
        do {
            // System.out.println("## dijkstra run " + i++ + " ##");

            // settled maintains whether node has been handled yet
            boolean[] settled = new boolean[2 * n];
            // prev points backwards along the augmenting path
            prev = new int[sinkIndex + 1];
            Arrays.fill(prev, -1);
            // distance to that node from source
            double distance[] = new double[sinkIndex + 1];
            if (minimize) {
                Arrays.fill(distance, Double.MAX_VALUE);
            } else {
                Arrays.fill(distance, -1);
            }

            // # Dijkstra algorithm to find min-cost augmenting path
            // from s to all on left side

            // initialize left side
            for (int j = 0; j < n; j++) {
                // no augmenting path starts at an already matched left node
                if (!matched[j]) {
                    queue.add(new Vertistance(j, 0.0));
                    distance[j] = 0.0;
                }
            }
            // back and forth between left and right
            while (!queue.isEmpty()) {
                int currentIndexAdjusted = queue.remove().index;

                if (settled[currentIndexAdjusted]) {
                    continue;
                } else {
                    settled[currentIndexAdjusted] = true;
                }

                boolean onLeftSide = (currentIndexAdjusted < n);
                int currentIndex = onLeftSide ? currentIndexAdjusted : currentIndexAdjusted - n;
                double currentDistance = distance[currentIndexAdjusted];

                // System.out.println("\n# currentVertex: " +
                // currentIndexAdjusted);
                // process all neighbours of current vertex
                for (int neighbourIndex = 0; neighbourIndex < n; neighbourIndex++) {
                    int neighbourIndexAdjusted = onLeftSide ? n + neighbourIndex : neighbourIndex;
                    // System.out.println("neighbour: " +
                    // neighbourIndexAdjusted);

                    // if current neighbor node hasn't already been processed
                    // and from left to right [right to left],
                    // the edge must not [must] be in the matching,
                    // then update distances
                    if ((!settled[neighbourIndexAdjusted]) //
                            && ((onLeftSide && !edgeUsed[currentIndex * n + neighbourIndex])
                            || (!onLeftSide
                            && edgeUsed[neighbourIndex * n + currentIndex]))) {

                        // left-to-right we add the edge cost,
                        // right-to-left we subtract it
                        double edgeDistance = onLeftSide
                                ? costMatrix[currentIndex * n + neighbourIndex]
                                : -costMatrix[neighbourIndex * n + currentIndex];
                        double newDistance = currentDistance + edgeDistance;

                        // if we lower the distance for the neighbour,
                        // then update values and add it to queue
                        if ((minimize && (newDistance < distance[neighbourIndexAdjusted]))
                                || (!minimize
                                && (newDistance > distance[neighbourIndexAdjusted]))) {

                            // System.out.println(
                            // currentDistance + " + " + edgeDistance + " = " +
                            // newDistance);

                            distance[neighbourIndexAdjusted] = newDistance;
                            prev[neighbourIndexAdjusted] = currentIndex;
                            queue.add(new Vertistance(neighbourIndexAdjusted, newDistance));
                        }
                    }
                }

                // but don't forget distance to sink
                // (if on right side and not matched yet)
                if (!onLeftSide && !matched[currentIndexAdjusted] && //
                        ((minimize && (currentDistance < distance[sinkIndex]))
                                || (!minimize && (currentDistance > distance[sinkIndex])))) {
                    distance[sinkIndex] = currentDistance;
                    prev[sinkIndex] = currentIndex;
                }
            }

            // double oldTotalCost = totalCost;
            if (prev[sinkIndex] != -1) {
                boolean onLeftSide = false;
                int leftIndex = -1;
                int rightIndex = prev[sinkIndex];
                matched[rightIndex + n] = true;

                // System.out.println("new augmenting path");
                // System.out.print("t - " + rightIndex + "r - ");

                do {
                    if (onLeftSide) {
                        // augmenting path backwards from left to right was over
                        // used edge
                        rightIndex = prev[leftIndex];
                        matched[rightIndex + n] = true;
                        edgeUsed[leftIndex * n + rightIndex] = false;
                        totalCost -= costMatrix[leftIndex * n + rightIndex];

                        // System.out.print(rightIndex + "r - ");

                    } else {
                        // augmenting path backwards from right to left was over
                        // unused edge, which now becomes used
                        leftIndex = prev[rightIndex + n];
                        matched[leftIndex] = true;
                        edgeUsed[leftIndex * n + rightIndex] = true;
                        totalCost += costMatrix[leftIndex * n + rightIndex];

                        // System.out
                        // .println(leftIndex + "l - " + ((prev[leftIndex] !=
                        // -1) ? "" : "s"));

                    }
                    onLeftSide = !onLeftSide;
                } while (prev[leftIndex] != -1);
            }
            // System.out
            // .println("totalCost: " + totalCost + " (" + (totalCost -
            // oldTotalCost) + ")\n");

        } while (prev[sinkIndex] != -1);

        return totalCost;
    }

}

class Vertistance implements Comparator<Vertistance> {
    public int index;
    public int distanceInt;
    public double distanceDouble;

    public Vertistance() {
    }

    public Vertistance(int node, int cost) {
        this.index = node;
        this.distanceInt = cost;
        this.distanceDouble = -1;
    }

    public Vertistance(int node, double cost) {
        this.index = node;
        this.distanceDouble = cost;
        this.distanceInt = -1;
    }

    @Override
    public int compare(Vertistance node1, Vertistance node2) {
        if (node1.distanceInt >= 0) {
            if (node1.distanceInt < node2.distanceInt) {
                return -1;
            }
            if (node1.distanceInt > node2.distanceInt) {
                return 1;
            }
            return 0;
        } else {
            if (node1.distanceDouble < node2.distanceDouble) {
                return -1;
            }
            if (node1.distanceDouble > node2.distanceDouble) {
                return 1;
            }
            return 0;

        }
    }
}
