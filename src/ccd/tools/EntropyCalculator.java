package ccd.tools;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.tree.Tree;
import beastfx.app.tools.Application;
import beastfx.app.treeannotator.TreeAnnotator;
import beastfx.app.treeannotator.TreeAnnotator.TreeSet;
import beastfx.app.util.TreeFile;
import ccd.model.CCD0;

import java.util.ArrayList;
import java.util.List;

@Description("Calculates entropy of a tree set")
public class EntropyCalculator extends beast.base.inference.Runnable {
    final public Input<List<TreeFile>> treeInput = new Input<>("trees", "trees to include in dissonance calculation", new ArrayList<>());
    final public Input<Integer> burnInPercentageInput = new Input<>("burnin", "percentage of trees to used as burn-in (and will be ignored)", 10);
    final public Input<Boolean> quietInput = new Input<>("quiet", "'true' to only output entropy value and nothing else, 'false' otherwise", false);
    final public Input<Boolean> dissonanceInput = new Input<>("dissonance", "calculate dissonance based on halving the tree set into two tree sets", false);
    final public Input<Boolean> summariseInput = new Input<>("summarise", "print summary of mean and variance of the entropies for all trees", false);

    @Override
    public void initAndValidate() {
        // nothing to do
    }

    @Override
    public void run() throws Exception {

        long start = System.currentTimeMillis();
        if (!quietInput.get()) {
            Log.info("# starting Entropy Calculator");
        }


        List<Double> entropies = new ArrayList<>();
        for (TreeFile t : treeInput.get()) {
            // init treeSets
            TreeSet treeSets = new TreeAnnotator().new MemoryFriendlyTreeSet(t.getPath(),
                    burnInPercentageInput.get());
            treeSets.reset();
            int numLeaves = treeSets.next().getLeafNodeCount();
            treeSets.reset();

            // init counts
            int numTrees = treeSets.totalTrees - treeSets.burninCount;

            // init CCDS
            CCD0 ccds = new CCD0(numLeaves, false);

            // process trees
            if (!quietInput.get()) {
                System.out.println("- processing trees");
            }
            int percentSize = numTrees / 50;
            if (!quietInput.get()) {
                Log.warning("#trees #clades #partitions dissonance time");
            }
            boolean dissonance = dissonanceInput.get();

            CCD0 ccdsFirstHalf = null;
            CCD0 ccdsSecondHalf = null;
            if (dissonance) {
                ccdsFirstHalf = new CCD0(numLeaves, false);
                ccdsSecondHalf = new CCD0(numLeaves, false);
            }

            for (int i = 0; i < numTrees; i++) {
                if (i % percentSize == 0 && i > 0) {
                    // System.out.print("+");
                    long end = System.currentTimeMillis();
                    if (!quietInput.get()) {
                        Log.warning.println(i + " " + (end - start) / 1000.0 + " seconds");
                    }
                }

                Tree tree = treeSets.next();
                ccds.addTree(tree);
                if (dissonance) {
                    if (i < numTrees / 2) {
                        ccdsFirstHalf.addTree(tree);
                    } else {
                        ccdsSecondHalf.addTree(tree);
                    }
                }
            }
            if (!quietInput.get()) {
                Log.warning("");
            }

            if (quietInput.get()) {
                Log.info(ccds.getEntropy() + (dissonance ? " " + (ccds.getEntropy() - (ccdsFirstHalf.getEntropy() + ccdsSecondHalf.getEntropy()) / 2.0) : ""));
            } else {
                Log.info("Entropy " + ccds.getEntropy() + (dissonance ? " " + (ccds.getEntropy() - (ccdsFirstHalf.getEntropy() + ccdsSecondHalf.getEntropy()) / 2.0) : ""));
            }
            entropies.add(ccds.getEntropy());
        }

        if (summariseInput.get()) {
            double sum = 0;
            for (double d : entropies) {
                sum += d;
            }
            double mean = sum / entropies.size();
            double sum2 = 0;
            for (double d : entropies) {
                sum2 += (d - mean) * (d - mean);
            }
            sum2 /= entropies.size();
            double stdev = Math.sqrt(sum2);
            Log.info("Mean entropy: " + mean + " (" + stdev + ")");
        }

        long end = System.currentTimeMillis();
        if (!quietInput.get()) {
            Log.info("Done in  " + (end - start) / 1000.0 + " seconds");
        }
    }

    public static void main(String[] args) throws Exception {
        new Application(new EntropyCalculator(), "Entropy Calculator", args);
    }

}
