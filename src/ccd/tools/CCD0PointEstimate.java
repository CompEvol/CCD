package ccd.tools;

import java.io.IOException;
import java.io.PrintStream;

import beast.base.core.Description;
import beast.base.core.Log;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator;
import beastfx.app.treeannotator.TreeAnnotator.TreeSet;
import beastfx.app.treeannotator.services.TopologySettingService;
import ccd.model.AbstractCCD;
import ccd.model.BitSet;
import ccd.model.CCD0;
import ccd.model.Clade;


@Description("TreeAnnotator plugin for setting tree topology by maximum CCD0")
public class CCD0PointEstimate extends PointEstimate implements TopologySettingService {

	@Override
	public Tree setTopology(TreeSet treeSet, PrintStream progressStream, TreeAnnotator annotator)
			throws IOException {

		progressStream.println("Maximum CCD0 Point Esitmate");
		progressStream.println("0              25             50             75            100");
		progressStream.println("|--------------|--------------|--------------|--------------|");
		
		treeSet.reset();
		Tree tree = treeSet.next();
		Tree firstTree = tree;
		CCD0 ccd = new CCD0(tree.getLeafNodeCount(),
				false);
		ccd.setProgressStream(progressStream);

		int k = treeSet.totalTrees - treeSet.burninCount;
		int percentageDone = 0;
		int i = 0;
		while (tree != null) {
			while ((62*i) /k > percentageDone) {
				progressStream.print("*");
				progressStream.flush();
				percentageDone++;
			}

			ccd.addTree(tree);
			tree = treeSet.hasNext() ? treeSet.next() : null;
			i++;
		}
		progressStream.println();

		Tree maxCCDTree = ccd.getMAPTree();
		
		
		// sanity checks
		doSanityCheck(maxCCDTree, firstTree, ccd);
		
		// set non-zero branch lenghts
		// otherwise all internal nodes are considered to be ancestral
		// which messes up the height, length and posterior annotations
		for (Node n : maxCCDTree.getRoot().getAllLeafNodes()) {
			double h = 0;
			do {
				n.setHeight(h);
				n = n.getParent();
				if (n != null) {
					h = Math.max(n.getLeft().getHeight() + 1, h);
					h = Math.max(n.getRight().getHeight() + 1, h);
				}
			} while (n != null);
		}
		return maxCCDTree;
	}

	@Override
	public String getServiceName() {
		return "CCD0";
	}

	@Override
	public String getDescription() {
		return "Conditional Clade Distribution 0";
	}

}
