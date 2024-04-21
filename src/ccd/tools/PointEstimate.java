package ccd.tools;

import beast.base.core.Log;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import ccd.model.AbstractCCD;
import ccd.model.BitSet;
import ccd.model.Clade;

/** base class for CCD0 and CCD1 point estimates, sharing sanity checks **/
public class PointEstimate {

	protected Tree doSanityCheck(Tree maxCCDTree, Tree firstTree, AbstractCCD ccd) {
		try {
			if (!sanityCheckPassed(maxCCDTree, ccd)) {
				
				// check if maxCCDTree equals the first tree
				String newick1 = toSortedNewick(maxCCDTree.getRoot(), new int[1]);
				String newick2 = toSortedNewick(firstTree.getRoot(), new int[1]);
				if (newick1.equals(newick2)) {
					Log.warning("The summary tree equals the first tree in the tree set!");
					Log.warning("This strongly suggests burn-in was not removed, and the summary tree is not valid.");

				
					Log.warning("The first tree is removed and the summary tree based on all following trees.");
					ccd.removeTree(firstTree, true);
					maxCCDTree = ccd.getMAPTree();
				}
			}
		} catch (Throwable e) {
			// we do not want sanity checks to ruin the job, 
			// so report any issues but otherwise ignore it.
			e.printStackTrace();
		}
		return maxCCDTree;
	}


	private String toSortedNewick(Node node, int[] maxNodeInClade) {
        StringBuilder buf = new StringBuilder();

        if (!node.isLeaf()) {
            buf.append("(");
            String child1 = toSortedNewick(node.getChild(0), maxNodeInClade);
            int child1Index = maxNodeInClade[0];
            String child2 = toSortedNewick(node.getChild(1), maxNodeInClade);
            int child2Index = maxNodeInClade[0];
            if (child1Index > child2Index) {
                buf.append(child2);
                buf.append(",");
                buf.append(child1);
            } else {
                buf.append(child1);
                buf.append(",");
                buf.append(child2);
                maxNodeInClade[0] = child1Index;
            }
            buf.append(")");
      } else {
            maxNodeInClade[0] = node.getNr();
            buf.append(node.getNr());
        }
        return buf.toString();
    }

	private boolean sanityCheckPassed(Tree maxCCDTree, AbstractCCD ccd) {
		// Sanity check: if burn-in was not properly removed, the starting tree has a 
		// good chance of becoming the MAP tree. To check this, we count how many of
		// the clades in maxCCDTree only have a single observation in the tree set.
		// The tree may have been under monophyly constraints, so we also count how
		// many clades are under 100% support.
		int [] counts = new int[2];
		traverse(maxCCDTree.getRoot(), counts, ccd);
		int singletons = counts[0];
		int monophyletics = counts[1];
		int nodeCount = maxCCDTree.getNodeCount();
		
		// A warning is produced if there is
		// o at least one singleton
		// o at most 10 clades that are not singleton or monophyletics
		if (singletons >= 1 && nodeCount - singletons - monophyletics <= 10) {
			Log.warning("WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING\n");
			Log.warning("There are "+ singletons +" clades supported by only a single tree in the tree set");
			Log.warning("This is likeley to happen when trees from burn-in ended up in the tree set.");
			Log.warning("Make sure brun-in is removed properly.");
			Log.warning("\nWARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING\n");
			return false;
		} else if (singletons >= 1) {
			Log.warning("Note: there are " +singletons + " clades in the summary tree supported by a single tree only.");
			Log.warning("The tree may not be well supported because there is a lot of topological uncertainty or ");
			Log.warning("burn-in was not removed from the tree set properly.");
			return false;
		}
		return true;
	}

	private BitSet traverse(Node node, int[] counts, AbstractCCD ccd) {
		if (node.isLeaf()) {
			BitSet bitset = BitSet.newBitSet(ccd.getNumberOfLeaves());
			bitset.set(node.getNr());
			return bitset;
		} else {
			BitSet bitset = traverse(node.getLeft(), counts, ccd);
			BitSet bitset2 = traverse(node.getRight(), counts, ccd);
			bitset.or(bitset2);
		
			Clade clade = ccd.getCladeMapping().get(bitset);
			int count = clade.getNumberOfOccurrences();
			if (count == 1) {
				 counts[0]++;
			}
			if (clade.isMonophyletic()) {
				counts[1]++;
			}
			return bitset;
		}
		
	}
}
