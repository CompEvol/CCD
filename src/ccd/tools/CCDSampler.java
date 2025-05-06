package ccd.tools;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Runnable;
import beastfx.app.tools.Application;
import beastfx.app.treeannotator.TreeAnnotator;
import beastfx.app.util.TreeFile;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.Random;

import ccd.model.AbstractCCD;
import ccd.model.CCDType;
import ccd.model.HeightSettingStrategy;

@Description("Allows to sample from a CCD{0,1} based on a input set of trees")
public class CCDSampler extends Runnable {
    final public Input<TreeFile> treeInput = new Input<>("trees", "trees file to construct CCD with and analyse", Input.Validate.REQUIRED);
    final public Input<String> outputInput = new Input<>("out", "file name for output newick trees", Input.Validate.REQUIRED);
    final public Input<Integer> burnInPercentageInput = new Input<>("burnin", "percentage of trees to used as burn-in (and will be ignored)", 10);
    final public Input<String> ccdTypeInput = new Input<>("ccdType", "CCD0 or CCD1", "CCD0");
    final public Input<Integer> sampleSizeInput = new Input<>("length", "number of trees sampled from CCD", 1000);
    final public Input<Long> seedInput = new Input<>("seed", "seed for random for chain generation");

    @Override
    public void initAndValidate() {
        // nothing to do
    }

    @Override
    public void run() throws Exception {
        Log.info.println("# Sample Trees from CCD");
        Log.info.println("> with the following parameters...");
        Log.info.println("    trees file:  " + treeInput.get().getPath());
        Log.info.println("    burnin:      " + burnInPercentageInput.get());
        Log.info.println("    CCD type:    " + CCDType.fromName(ccdTypeInput.get()));
        Log.info.println("    #samples:    " + sampleSizeInput.get());
        Log.info.println("    output file: " + outputInput.get());

        TreeAnnotator.MemoryFriendlyTreeSet treeSet = CCDToolUtil.getTreeSet(treeInput, burnInPercentageInput.get());
        AbstractCCD ccd = CCDToolUtil.getCCDTypeByName(treeSet, ccdTypeInput);

        long seed = (seedInput.get() != null) ? seedInput.get() : System.currentTimeMillis();
        ccd.setRandom(new Random(seed));

        try (BufferedWriter bufferedOutputWriter = new BufferedWriter(new FileWriter(outputInput.get()))) {
            for (int i = 0; i < sampleSizeInput.get(); i++) {
                Tree tree = ccd.sampleTree(HeightSettingStrategy.CommonAncestorHeights);
                bufferedOutputWriter.write(tree.getRoot().toNewick());
                bufferedOutputWriter.newLine();
            }
        }
    }

    public static void main(String[] args) throws Exception {
        new Application(new CCDSampler(), "Sample Assistant", args);
    }
}
