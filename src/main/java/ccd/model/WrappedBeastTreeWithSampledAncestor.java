package ccd.model;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import ccd.model.bitsets.BitSet;

import java.util.ArrayList;
import java.util.List;

public class WrappedBeastTreeWithSampledAncestor extends WrappedBeastTree {

    public WrappedBeastTreeWithSampledAncestor(Tree wrappedTree) {
        super(wrappedTree);
    }

    @Override
    protected BitSet initCladeBitSet(Node vertex) {
        BitSet cladeAsBitSet;
        if (vertex.isLeaf()) {
            cladeAsBitSet = leafKey(vertex.getNr());
        } else {
            for (Node child : vertex.getChildren()) {
                initCladeBitSet(child);
            }
            cladeAsBitSet = extendedSelfKey(vertex);
        }
        cladeOfVertex[vertex.getNr()] = cladeAsBitSet;
        return cladeAsBitSet;
    }

    private BitSet leafKey(int taxonIndex) {
        BitSet key = BitSet.newBitSet(2 * wrappedTree.getLeafNodeCount());
        key.set(taxonIndex);
        return key;
    }

    private BitSet extendedSelfKey(Node vertex) {
        BitSet key = BitSet.newBitSet(2 * wrappedTree.getLeafNodeCount());
        collectTaxaBits(vertex, key);
        for (Node child : vertex.getChildren()) {
            if (child.isLeaf() && child.getLength() == 0) {
                key.set(wrappedTree.getLeafNodeCount() + child.getNr());
            }
        }
        return key;
    }

    private void collectTaxaBits(Node vertex, BitSet key) {
        if (vertex.isLeaf()) {
            key.set(vertex.getNr());
        } else {
            for (Node child : vertex.getChildren()) {
                collectTaxaBits(child, key);
            }
        }
    }

    public int getNumberOfSampledAncestors() {
        int numberOfSampledAncestors = 0;
        List<Node> leaves = wrappedTree.getExternalNodes();
        for (Node leaf : leaves) {
            if (leaf.getLength() == 0) {
                numberOfSampledAncestors++;
            }
        }
        return numberOfSampledAncestors;
    }

    public ArrayList<Node> getSampledAncestors() {
        ArrayList<Node> sampledAncestorNodes = new ArrayList<>();
        List<Node> leaves = wrappedTree.getExternalNodes();
        for (Node leaf : leaves) {
            if (leaf.getLength() == 0) {
                sampledAncestorNodes.add(leaf);
            }
        }
        return sampledAncestorNodes;
    }
}
