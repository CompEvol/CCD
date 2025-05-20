package ccd.tools;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.tree.Tree;
import beastfx.app.tools.Application;
import beastfx.app.treeannotator.TreeAnnotator.MemoryFriendlyTreeSet;
import beastfx.app.util.OutFile;
import beastfx.app.util.TreeFile;
import ccd.algorithms.credibleSets.CredibleCCDComputer;
import ccd.algorithms.credibleSets.CredibleSetType;
import ccd.algorithms.credibleSets.ICredibleSet;
import ccd.algorithms.credibleSets.ProbabilityBasedCredibleSetComputer;
import ccd.model.AbstractCCD;
import ccd.model.CCD0;

import java.io.FileWriter;
import java.io.IOException;

import static ccd.algorithms.credibleSets.ProbabilityBasedCredibleSetComputer.DEFAULT_NUM_SAMPLES;

@Description("Compute probability and credible level of given tree in CCD of given tree set")
public class TreeCredibleLevel extends beast.base.inference.Runnable {
    // input
    final public Input<TreeFile> treeInput = new Input<>("trees", "trees file to construct CCD", Input.Validate.REQUIRED);
    final public Input<Integer> burnInPercentageInput = new Input<>("burnin", "percentage of trees that is burnin (and will be ignored)", 10);
    final public Input<TreeFile> testTreeInput = new Input<>("tree", "tree file of tree to compute probability in CCD", Input.Validate.REQUIRED);

    // configuration
    final public Input<Boolean> quietInput = new Input<>("quiet", "'true' to only output entropy value and nothing else, 'false' otherwise", false);
    final public Input<String> ccdTypeInput = new Input<>("ccdType", "CCD0, CCD1, or CCD2", "CCD1");
    final public Input<String> methodInput = new Input<>("method", "'probability' for a probability-based credible set (default), " +
            "'credibleCCD' for a credible CCD", "probability");
    final public Input<Integer> numSamplesInput = new Input<>("numsamples",
            "if using probability-based credible set, you may set the number of trees sampled " +
                    "(default: " + DEFAULT_NUM_SAMPLES + ")", DEFAULT_NUM_SAMPLES);

    // output
    final public Input<OutFile> outputInput = new Input<>("out", "output file", Input.Validate.REQUIRED);

    @Override
    public void initAndValidate() {
        // nothing to do
    }

    @Override
    public void run() throws Exception {
        boolean quiet = quietInput.get();

        if (!quiet) {
            Log.info("# Compute the credible level of a tree in a CCD (type: " + ccdTypeInput.get() + ")");
        }
        String method = methodInput.get();
        if (!(method.equalsIgnoreCase("probability") || method.equalsIgnoreCase("credibleCCD"))) {
            Log.info("Unknown method, fall back to default 'probability'.");
            method = "probability";
        }

        // init CCD
        MemoryFriendlyTreeSet treeSetCCD = CCDToolUtil.getTreeSet(treeInput, burnInPercentageInput.get());
        AbstractCCD ccd = CCDToolUtil.getCCDTypeByName(treeSetCCD, ccdTypeInput);

        // init test tree
        MemoryFriendlyTreeSet treeSetTest = CCDToolUtil.getTreeSet(testTreeInput, 0);
        treeSetTest.reset();
        Tree testTree = treeSetTest.next();

        // credible set init
        if (!quiet) {
            Log.info("Build credible set.");
        }
        ICredibleSet cred = null;
        if (method.equalsIgnoreCase("probability")) {
            int numSamples = numSamplesInput.get();
            cred = new ProbabilityBasedCredibleSetComputer(ccd, numSamples);
        } else {
            cred = CredibleCCDComputer.getCredibleCCDComputer(ccd, (ccd instanceof CCD0) ? CredibleSetType.CladeProbability : CredibleSetType.PartitionProbability);
        }

        // calculations
        double p = ccd.getProbabilityOfTree(testTree);
        int credLevel = (int) Math.ceil(cred.getCredibleLevel(testTree));

        // output
        if (!quiet) {
            Log.info("Write results .");
        }
        try (FileWriter writer = new FileWriter(outputInput.get())) {
            writer.write(p + "\n");
            writer.write(credLevel + "\n");
        } catch (IOException e) {
            e.printStackTrace();
        }

        if (!quiet) {
            Log.info("Done.");
        }
    }

    public static void main(String[] args) throws Exception {
        new Application(new TreeCredibleLevel(), "CCD-based Tree Credible Level and Probability Calculator", args);
    }

}
