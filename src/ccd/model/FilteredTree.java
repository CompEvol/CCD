package ccd.model;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;

public class FilteredTree extends Tree {

    public FilteredTree(final Node rootNode) {
        setRoot(rootNode);
        initArrays();
    }

    @Override
    public void setRoot(Node root) {
        this.root = root;
        nodeCount = this.root.getNr() + 1;
    }
}
