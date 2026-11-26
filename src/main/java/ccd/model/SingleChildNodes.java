package ccd.model;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;

/**
 * Internal nodes with a single child, such as the root above the trunk in SIEVE trees
 * (CompEvol/CCD#22). The CCDs need every internal node to have two children, so
 * {@link AbstractCCD} collapses unnamed single-child nodes before it processes a tree.
 * A named single-child node is refused, as it may be a sampled ancestor written as a
 * labelled internal node rather than as a zero-length leaf.
 */
final class SingleChildNodes {

    private SingleChildNodes() {
    }

    /**
     * @return the number of unnamed single-child nodes in the tree
     * @throws IllegalArgumentException if a single-child node is named
     */
    static int countUnnamed(Tree tree) {
        int count = 0;
        for (Node node : tree.getNodesAsArray()) {
            if (node != null && node.getChildCount() == 1) {
                if (node.getID() != null) {
                    throw new IllegalArgumentException("Node '" + node.getID() + "' has a single child, "
                            + "but CCD needs every internal node to have two children. "
                            + "If it is a sampled ancestor, write it as a leaf with a zero-length branch, "
                            + "the form the SA package uses.");
                }
                count++;
            }
        }
        return count;
    }

    /**
     * Returns a copy of the tree in which every single-child node is replaced by its child.
     * Leaves keep their numbers; internal nodes are renumbered after them, with the root last.
     * Heights, IDs and metadata are kept.
     *
     * @param copyToOriginal if not null, receives the original node of every node in the copy
     */
    static Tree collapse(Tree tree, Map<Node, Node> copyToOriginal) {
        List<Node> internal = new ArrayList<>();
        Node root = copy(tree.getRoot(), internal, copyToOriginal);
        int nr = tree.getLeafNodeCount();
        for (Node node : internal) { // post-order, so the root comes last
            node.setNr(nr++);
        }
        return new Tree(root);
    }

    private static Node copy(Node node, List<Node> internal, Map<Node, Node> copyToOriginal) {
        while (node.getChildCount() == 1) {
            node = node.getChild(0);
        }
        Node copy = new Node();
        copy.setID(node.getID());
        copy.metaDataString = node.metaDataString;
        copy.lengthMetaDataString = node.lengthMetaDataString;
        for (String name : node.getMetaDataNames()) {
            copy.setMetaData(name, node.getMetaData(name));
        }
        for (String name : node.getLengthMetaDataNames()) {
            copy.setLengthMetaData(name, node.getLengthMetaData(name));
        }
        copy.setHeight(node.getHeight()); // after the metadata: a "date" key would set the height
        if (node.isLeaf()) {
            copy.setNr(node.getNr());
        } else {
            for (Node child : node.getChildren()) {
                copy.addChild(copy(child, internal, copyToOriginal));
            }
            internal.add(copy);
        }
        if (copyToOriginal != null) {
            copyToOriginal.put(copy, node);
        }
        return copy;
    }
}
