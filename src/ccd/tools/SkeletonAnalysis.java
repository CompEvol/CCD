package ccd.tools;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Log;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Runnable;
import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.tools.Application;
import beastfx.app.treeannotator.TreeAnnotator.MemoryFriendlyTreeSet;
import beastfx.app.util.OutFile;
import beastfx.app.util.TreeFile;
import ccd.algorithms.RogueDetection;
import ccd.model.AbstractCCD;
import ccd.model.BitSet;
import ccd.model.CCDType;
import ccd.model.FilteredCCD;
import ccd.model.HeightSettingStrategy;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

import static ccd.algorithms.RogueDetection.TerminationStrategy.*;

@Description("Analyses the skeleton of a CCD for given trees based on the total rogue scores (based on entropy)")
public class SkeletonAnalysis extends Runnable {
    // input
    final public Input<TreeFile> treeInput = new Input<>("trees", "trees file to construct CCD with and analyse", Validate.REQUIRED);
    final public Input<Integer> burnInPercentageInput = new Input<>("burnin", "percentage of trees to used as burn-in (and will be ignored)", 10);

    // configuration
    // - for skeleton computation
    final public Input<String> ccdTypeInput = new Input<>("ccdType", "CCD0 or CCD1", "CCD0");
    final public Input<Integer> maxCladeSizeInput = new Input<>("maxCladeSize", "maximum size for clade to be considered for removing", 10);
    final public Input<Double> minProbabilityInput = new Input<>("minProbability", "minimum probability for clade to analyse", 0.5);
    final public Input<RogueDetection.RogueDetectionStrategy> strategyInput = new Input<>("detectionStrategy",
            "rogue detection strategy (default recommended); one of " + Arrays.toString(RogueDetection.RogueDetectionStrategy.values()),
            RogueDetection.RogueDetectionStrategy.Entropy,
            RogueDetection.RogueDetectionStrategy.values());
    final public Input<RogueDetection.TerminationStrategy> terminationInput = new Input<>("terminationStrategy",
            "termination strategy (default recommended); one of " + Arrays.toString(RogueDetection.TerminationStrategy.values()),
            RogueDetection.TerminationStrategy.Entropy,
            RogueDetection.TerminationStrategy.values());
    final public Input<Double> terminationThresholdInput = new Input<>("terminationThreshold",
            "threshold for termination strategy (if not default or exhaustive strategy)");

    // - for MAP tree
    enum hss {CA, MH, ONE}

    final public Input<hss> heightSettingStrategyInput = new Input<>("heightSettingStrategy",
            "heights used in MAP tree output, can be CA (Mean of Least Common Ancestor heights), MH (mean (sampled) height), or ONE",
            hss.CA,
            hss.values());
    final public Input<OutFile> excludeInput = new Input<>("exclude", "file name of text file containing taxa to exclude from filtering" +
            " -- can be comma, tab or newline delimited.");

    // output
    final public Input<OutFile> outputInput = new Input<>("out", "reduced tree output file; the given tree set will not be filtered if not specified");

    @Override
    public void initAndValidate() {
    }

    @Override
    public void run() throws Exception {
        System.out.println("# Skeleton Analysis");
        System.out.println("> with the following parameters...");
        System.out.println("    trees file: " + treeInput.get().getPath());
        System.out.println("    burnin:     " + burnInPercentageInput.get());
        System.out.println("    CCD type:   " + CCDType.fromName(ccdTypeInput.get()));

        // create base CCD
        MemoryFriendlyTreeSet treeSet = CCDToolUtil.getTreeSet(treeInput, burnInPercentageInput.get());
        AbstractCCD ccd = CCDToolUtil.getCCDTypeByName(treeSet, ccdTypeInput);

        // parse parameters
        RogueDetection.TerminationStrategy terminationStrategy = terminationInput.get();
        if (terminationThresholdInput.get() == null) {
            switch (terminationStrategy) {
                case Entropy -> terminationStrategy.setThreshold(ENTROPY_THRESHOLD_DEFAULT);
                case MaxProbability -> terminationStrategy.setThreshold(MAX_PROB_THRESHOLD_DEFAULT);
                case NumRogues -> terminationStrategy.setThreshold(Math.min(NUM_ROGUES_THRESHOLD_DEFAULT, (double) ccd.getNumberOfLeaves() / 2));
                case Support -> terminationStrategy.setThreshold(SUPPORT_THRESHOLD_DEFAULT);
            }
        } else {
            terminationStrategy.setThreshold(terminationThresholdInput.get());
        }

        // run analysis
        int maxCladeSize = Math.max(maxCladeSizeInput.get(), 1);
        maxCladeSize = Math.min(maxCladeSize, ccd.getNumberOfLeaves() - 1);
        double minProbability = minProbabilityInput.get();
        ArrayList<AbstractCCD> ccds = RogueDetection.detectRoguesWhileImproving(ccd,
                maxCladeSize, strategyInput.get(), terminationStrategy, minProbability, true);

        // report removed clades
        System.out.println("Final rogue removal sequence to obtain skeleton");
        System.out.println("n - entropy - num clades - removed taxa");
        for (AbstractCCD ccdi : ccds) {
            System.out.println(ccdi.getRootClade().getCladeInBits().cardinality() + " - " //
                    + ccdi.getEntropy() + " - " //
                    + ccdi.getNumberOfClades() + " - "//
                    + ((ccdi instanceof FilteredCCD) ? ccd.getTaxaNames(((FilteredCCD) ccdi).getRemovedTaxaMask()) : ""));
        }

        // filter tree set if outputInput is specified
        filterTrees(treeSet, ccds);

        // annotate MAP tree
        System.out.println("\n> annotate CCD MAP tree");
        FilteredCCD lastCCD = (FilteredCCD) ccds.get(ccds.size() - 1);
        HeightSettingStrategy hss = HeightSettingStrategy.fromName(heightSettingStrategyInput.get().toString());
        Tree tree = lastCCD.getMAPTree(hss);
        Set<BitSet> rogues = RogueDetection.extractRogues(ccds);
        RogueDetection.annotateRoguePlacements(ccd, lastCCD, rogues, tree);
        System.out.println("The resulting annotated tree is: ");
        System.out.println(tree.getRoot().toNewick());
    }

    /* Extracted code that filters source treeset. */
    private void filterTrees(MemoryFriendlyTreeSet treeSet, ArrayList<AbstractCCD> ccds) throws IOException {
        if (outputInput.get() != null && !outputInput.get().getName().equals("[[none]]")) {

            System.out.println("\n> Output filtered trees to " + outputInput.get().getPath());
            treeSet.reset();

            // compute which taxa to include and which to exclude
            Set<String> taxaToInclude = new HashSet<>();
            Set<String> taxaToExclude = new HashSet<>();
            Tree tree = treeSet.next();
            for (int i = 0; i < tree.getLeafNodeCount(); i++) {
                taxaToInclude.add(tree.getNode(i).getID());
            }

            // load set of taxa to exclude from deletion
            Set<String> excludeFromDeletion = new HashSet<>();
            if (excludeInput.get() != null && !excludeInput.get().getName().equals("[[none]]")) {
                String str = BeautiDoc.load(excludeInput.get());
                String[] strs = str.split("[,\t\n]");
                for (String str2 : strs) {
                    excludeFromDeletion.add(str2);
                    if (!taxaToInclude.contains(str2)) {
                        Log.warning("Taxon name >>" + str2 + "<< could not be found in tree");
                    }
                }
            }

            FilteredCCD lastCCD = (FilteredCCD) ccds.get(ccds.size() - 1);
            AbstractCCD ccd = ccds.get(0);
            BitSet taxaToKeepBits = lastCCD.getTaxaAsBitSet();
            Set<String> taxaNamesToKeep = ccd.getTaxaNamesList(taxaToKeepBits);
            taxaNamesToKeep.addAll(excludeFromDeletion);

            // process file
            PrintStream out = new PrintStream(outputInput.get());
            treeSet.reset();
            while (treeSet.hasNext()) {
                Node root = treeSet.next().getRoot();
                root = filterTree(root, taxaToInclude);
                out.println(root.toNewick());
            }
            out.close();
        }

        System.out.println("... done.");
    }

    /**
     * Recursive helper method that kicks out all taxa not specified to keep in subtree below given vertex.
     *
     * @param node          on which to recurse
     * @param taxaToInclude taxon in this set is kept, all others are excluded
     * @return this vertex, a new one for some changes, or null if subtree below is completely lost
     */
    public static Node filterTree(Node node, Set<String> taxaToInclude) {
        /* Copied from Babel package TaxonFiler */
        if (node.isLeaf()) {
            if (taxaToInclude.contains(node.getID())) {
                return node;
            } else {
                return null;
            }
        } else {
            Node left_ = node.getLeft();
            Node right_ = node.getRight();
            left_ = filterTree(left_, taxaToInclude);
            right_ = filterTree(right_, taxaToInclude);
            if (left_ == null && right_ == null) {
                return null;
            }
            if (left_ == null) {
                return right_;
            }
            if (right_ == null) {
                return left_;
            }
            node.removeAllChildren(false);
            node.addChild(left_);
            node.addChild(right_);
            return node;
        }
    }

    public static void main(String[] args) throws Exception {
        new Application(new SkeletonAnalysis(), "Rogue Entropy Sweeper", args);
    }

}
