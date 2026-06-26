package ccd.model;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator;
import ccd.model.bitsets.BitSet;

import java.util.*;

public class CCD0CP extends CCD0 {

    public CCD0CP(List<Tree> trees, double burnin) {
        super(trees, burnin);
    }

    public CCD0CP(TreeAnnotator.TreeSet treeSet) {
        super(treeSet);
    }

    public CCD0CP(TreeAnnotator.TreeSet treeSet, int numTreesToUse) {
        super(treeSet, numTreesToUse);
    }

    public CCD0CP(TreeAnnotator.TreeSet treeSet, boolean storeBaseTrees) {
        super(treeSet, storeBaseTrees);
    }

    public CCD0CP(int numLeaves, boolean storeBaseTrees) {
        super(numLeaves, storeBaseTrees);
    }

    /* -- ROOT INITIALIZATION -- */
    @Override
    protected void initializeRootClade(int numLeaves) {
        CPSupport.initializeRoot(this, numLeaves);
    }

    @Override
    protected Clade cladifyVertex(Node vertex) {
        return CPSupport.cladifyVertex(this, vertex);
    }

    @Override
    public boolean isSampledAncestor(Clade clade) {
        return CPSupport.isSampledAncestor(clade, leafArraySize);
    }

    @Override
    protected String getSampledAncestorInfoString(Clade clade) {
        return CPSupport.saInfoString(clade, leafArraySize);
    }

    @Override
    protected double computeParentHeight(CladePartition partition, Node firstChild, Node secondChild) {
        return CPSupport.computeParentHeight(partition, firstChild, secondChild, leafArraySize);
    }

    @Override
    protected void expand() {
        List<Clade> realClades = new ArrayList<>();
        Map<BitSet, List<Clade>> taxaMaskToClades = new HashMap<>();
        for (Clade c : cladeMapping.values()) {
            realClades.add(c);
            BitSet taxaMask = c.getCladeInBitsTaxaOnly();
            taxaMaskToClades.computeIfAbsent(taxaMask, k -> new ArrayList<>()).add(c);
        }

        int n = realClades.size();
        for (int i = 0; i < n; i++) {
            Clade A = realClades.get(i);
            BitSet taxaA = A.getCladeInBitsTaxaOnly();
            for (int j = i + 1; j < n; j++) {
                Clade B = realClades.get(j);
                BitSet taxaB = B.getCladeInBitsTaxaOnly();

                // Taxa-disjointness check via base mask.
                if (taxaA.intersects(taxaB)) continue;

                // Compute taxa union and look up all parent clades sharing
                // that taxa-mask (any number of SA variants).
                BitSet unionMask = (BitSet) taxaA.clone();
                unionMask.or(taxaB);

                List<Clade> parents = taxaMaskToClades.get(unionMask);
                if (parents == null) continue;

                for (Clade parent : parents) {
                    if (parent == A || parent == B) continue;
                    if (parent.size() == 2) {
                        // cannot have a cherry with both children being SA
                        if (isSampledAncestor(A) && isSampledAncestor(B)) continue;
                    }
                    if (parent.getCladePartition(A, B) != null) continue;
                    parent.createCladePartition(A, B);
                }
            }
        }
    }

    @Override
    protected Clade computeProbabilityOfVertex(Node vertex, double[] runningProbability, boolean computeLog) {
        return CPSupport.computeProbCPVertex(vertex, runningProbability, computeLog, this, leafArraySize);
    }

    @Override
    public double setPartitionProbabilities(Clade clade, boolean useCladeParameters) {

        // usually initialise CCD0 uses getCladeCredibility
        double cladeValue = useCladeParameters ? clade.getCladeParameter() : clade.getCladeCredibility();

        if (clade.isLeaf()) {
            // Leaves carry cladeValue (their extended-clade credibility), NOT 1, so the
            // leaf branch must run BEFORE the memo guard: resetSumCladeCredibilities sets
            // leaves to 1, which the guard would otherwise return as a stale CCP.
            clade.setSumCladeCredibilities(cladeValue);
            return cladeValue;
        } else {
            // Memoise internal clades: the clade DAG has heavy shared substructure, so
            // without this each clade is recomputed once per path that reaches it
            // (exponential). Sentinel -1 = not computed this pass
            // (Clade.resetSumCladeCredibilities sets internal clades to -1).
            if (clade.getSumCladeCredibilities() >= 0) {
                return clade.getSumCladeCredibilities();
            }
            double sumSubtreeProbabilities = 0.0;
            double[] sumPartitionSubtreeProbabilities = new double[clade.getPartitions().size()];

            // compute sum of probabilities over all partitions ...
            int i = 0;
            for (CladePartition partition : clade.getPartitions()) {
                sumPartitionSubtreeProbabilities[i] =
                        setPartitionProbabilities(partition.getChildClades()[0], useCladeParameters)
                                * setPartitionProbabilities(partition.getChildClades()[1], useCladeParameters);
                sumSubtreeProbabilities += sumPartitionSubtreeProbabilities[i];
                i++;
            }

            // ... and then normalize
            if (sumSubtreeProbabilities == 0) {
                // probability of this subtree is so small, that we have an underflow problem;
                // we can try to use log transformed probabilities to set CCPs
                // but the probability of the clade will still become zero

                // try log-sum-exp trick
                double logMax = Double.NEGATIVE_INFINITY;
                double[] logProbs = new double[clade.getPartitions().size()];
                i = 0;
                for (CladePartition partition : clade.getPartitions()) {
                    double left = partition.getChildClades()[0].getSumCladeCredibilities();
                    double right = partition.getChildClades()[1].getSumCladeCredibilities();
                    logProbs[i] = Math.log(left) + Math.log(right);
                    logMax = Math.max(logProbs[i], logMax);
                    i++;
                }

                // log(Σ exp(xᵢ)) = max + log(Σ exp(xᵢ - max))
                double intermSum = 0;
                for (int j = 0; j < logProbs.length; j++) {
                    if (Double.isFinite(logProbs[j])) {
                        intermSum += Math.exp(logProbs[j] - logMax);
                    }
                }
                double logSum = logMax + Math.log(intermSum);

                // normalizing with normalized log pi = log pi - log sum
                i = 0;
                double pSum = 0;
                for (CladePartition partition : clade.getPartitions()) {
                    if (!Double.isFinite(logProbs[i])) {
                        // we are still encountering underflows
                        throw new UnderflowException("An underflow has occurred.");
                    } else {
                        double logProbability = logProbs[i] - logSum;
                        double probability = Math.exp(logProbability);
                        pSum += probability;

                        if (Double.isNaN(probability)) {
                            out.println("NaN probability = " + probability);
                            out.println("logProbability = " + logProbability);
                            out.println("logProbs[i] = " + logProbs[i]);
                            out.println("logSum = " + logSum);
                        }

                        partition.setCCP(probability);
                    }
                    i++;
                }
            } else {
                i = 0;
                for (CladePartition partition : clade.getPartitions()) {
                    double probability = sumPartitionSubtreeProbabilities[i] / sumSubtreeProbabilities;
                    if (Double.isNaN(probability)) {
                        out.println("clade = " + clade);
                        out.println("partition = " + partition);
                        out.println("sumPartitionSubtreeProbabilities = " + sumPartitionSubtreeProbabilities[i]);
                        out.println("sumSubtreeProbabilities = " + sumSubtreeProbabilities);
                    }
                    partition.setCCP(probability);
                    i++;
                }
            }

            // combined with probability of clade, we get sum of all subtree probabilities
            double sumCladeCredibilities = sumSubtreeProbabilities * cladeValue;
            clade.setSumCladeCredibilities(sumCladeCredibilities);
            return sumCladeCredibilities;
        }
    }

    @Override
    public String toString() {
        return "CCD0-CP " + super.toString().replaceFirst("CCD0 ", "");
    }

    @Override
    public AbstractCCD copy() {
        CCD0CP copy = new CCD0CP(this.getSizeOfLeavesArray(), false);
        copy.baseTrees.add(this.getSomeBaseTree());
        copy.numBaseTrees = this.getNumberOfBaseTrees();
        AbstractCCD.buildCopy(this, copy);
        return copy;
    }

}
