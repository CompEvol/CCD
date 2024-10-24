package ccd.tools;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beastfx.app.tools.Application;
import beastfx.app.treeannotator.TreeAnnotator.MemoryFriendlyTreeSet;
import beastfx.app.util.TreeFile;
import ccd.model.AbstractCCD;
import ccd.model.CCDType;

@Description("Calculates the phylogenetic entropy of the posterior tree distribution")
public class EntropyCalculator extends beast.base.inference.Runnable {
    // input
    final public Input<TreeFile> treeInput = new Input<>("trees", "trees file for which to compute entropy", Input.Validate.REQUIRED);
    final public Input<Integer> burnInPercentageInput = new Input<>("burnin", "percentage of trees to used as burn-in (and will be ignored)", 10);

    // configuration
    final public Input<String> ccdTypeInput = new Input<>("ccdType", "CCD0, CCD1, or CCD2", "CCD0");

    @Override
    public void initAndValidate() {
        // nothing to do
    }

    @Override
    public void run() throws Exception {
        Log.info("# Entropy Calculator");
        Log.info("> with the following parameters...");
        Log.info("    trees file: " + treeInput.get().getPath());
        Log.info("    burnin:     " + burnInPercentageInput.get());
        Log.info("    CCD type:   " + CCDType.fromName(ccdTypeInput.get()));

        MemoryFriendlyTreeSet treeSet = CCDToolUtil.getTreeSet(treeInput, burnInPercentageInput.get());
        AbstractCCD ccd = CCDToolUtil.getCCDTypeByName(treeSet, ccdTypeInput);

        Log.info("\nPhylogenetic Entropy: " + ccd.getEntropy());
    }

    public static void main(String[] args) throws Exception {
        new Application(new EntropyCalculator(), "Entropy Calculator", args);
    }

}
