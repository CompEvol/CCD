package ccd.tools;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beastfx.app.tools.Application;
import beastfx.app.treeannotator.TreeAnnotator.MemoryFriendlyTreeSet;
import beastfx.app.util.OutFile;
import beastfx.app.util.TreeFile;
import ccd.model.AbstractCCD;
import ccd.model.CCDType;

@Description("Calculates the phylogenetic entropy of the posterior tree distribution")
public class EntropyCalculator extends beast.base.inference.Runnable {
    // input
    final public Input<List<TreeFile>> treeInput = new Input<>("trees", "trees file for which to compute entropy", new ArrayList<>(), Input.Validate.REQUIRED);
    final public Input<Integer> burnInPercentageInput = new Input<>("burnin", "percentage of trees to used as burn-in (and will be ignored)", 10);
    final public Input<Boolean> verboseInput = new Input<>("verbose", "print info to stdout -- set to false to minimise output", true);
    final public Input<OutFile> outputInput = new Input<>("output", "name of the output file to store results. If not specified, use stdout");

    // configuration
    final public Input<String> ccdTypeInput = new Input<>("ccdType", "CCD0, CCD1, or CCD2", "CCD0");

    @Override
    public void initAndValidate() {
        // nothing to do
    }

    @Override
    public void run() throws Exception {
    	
    	
    	PrintStream trace = null;
    	if (outputInput.get() != null) {
    		trace = new PrintStream(outputInput.get());
    		trace.println("sample\tfile\tentropy");
    	}
    		
    	int k = 0;
    	
    	for (TreeFile treefile : treeInput.get()) {
    		if (verboseInput.get()) {
		        Log.info("# Entropy Calculator");
		        Log.info("> with the following parameters...");
		        Log.info("    trees file: " + treefile.getPath());
		        Log.info("    burnin:     " + burnInPercentageInput.get());
		        Log.info("    CCD type:   " + CCDType.fromName(ccdTypeInput.get()));
    		}
	
	        MemoryFriendlyTreeSet treeSet = CCDToolUtil.getTreeSet(treefile.getPath(), burnInPercentageInput.get());
	        AbstractCCD ccd = CCDToolUtil.getCCDTypeByName(treeSet, ccdTypeInput);
	        
	        double entropy = ccd.getEntropy();
	
	        if (verboseInput.get()) {
	        	Log.info("\nPhylogenetic Entropy: " + entropy);
	        }
	        if (trace != null) {
	        	trace.println(k + "\t" + treefile.getName() + "\t" + entropy);
	        }
	        k++;
    	}
    }

    public static void main(String[] args) throws Exception {
        new Application(new EntropyCalculator(), "Entropy Calculator", args);
    }

}
