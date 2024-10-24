package ccd.tools;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Log;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Runnable;
import beastfx.app.tools.Application;
import beastfx.app.treeannotator.TreeAnnotator.MemoryFriendlyTreeSet;
import beastfx.app.util.OutFile;
import beastfx.app.util.TreeFile;
import ccd.algorithms.RogueDetection;
import ccd.model.AbstractCCD;
import ccd.model.CCDType;
import ccd.model.Clade;
import ccd.model.HeightSettingStrategy;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Map;

@Description("Analyse the clades in a CCD for given trees based on their clade rogue score (based on entropy)")
public class RogueAnalysis extends Runnable {

    public static final String ROGUE_SCORE_KEY = "rogueScore";
    public static final String SEPARATOR = "\t";

    // input
    final public Input<TreeFile> treeInput = new Input<>("trees", "trees file to construct CCD with and analyse", Validate.REQUIRED);
    final public Input<Integer> burnInPercentageInput = new Input<>("burnin", "percentage of trees that is burnin (and will be ignored)", 10);

    // configuration
    // - for rogue score
    final public Input<String> ccdTypeInput = new Input<>("ccdType", "CCD0 or CCD1", "CCD0");
    final public Input<Integer> maxCladeSizeInput = new Input<>("maxCladeSize", "maximum size for clade to be analysed", 10);
    final public Input<Double> minProbabilityInput = new Input<>("minProbability", "minimum probability for clade to be analysed", 0.1);

    // - for MAP tree
    enum hss {CA, MH, ONE}

    final public Input<hss> heightSettingStrategyInput = new Input<>("heightSettingStrategy",
            "heights used in MAP tree output, can be CA (Mean of Least Common Ancestor heights), MH (mean (sampled) height), or ONE",
            hss.CA,
            hss.values());
    //    final public Input<HeightSettingStrategy> heightSettingStrategyInput = new Input<>("heightSettingStrategy",
    //            "heights used in MAP tree output, can be CA (Mean of Least Common Ancestor heights), MH (mean (sampled) height), or ONE",
    //            HeightSettingStrategy.CommonAncestorHeights,
    //            HeightSettingStrategy.values());

    // output
    final public Input<OutFile> outputInput = new Input<>("out", "file name for output (without file ending)," +
            "will be used with '.csv' for rogue score and '.trees' for annotated MAP trees", Validate.REQUIRED);
    final public Input<String> separatorInput = new Input<>("separator",
            "separator used in csv file; default is tab", SEPARATOR);

    @Override
    public void initAndValidate() {
    }

    @Override
    public void run() throws Exception {
        System.out.println("# Rogue Analysis");
        System.out.println("> with the following parameters...");
        System.out.println("    trees file: " + treeInput.get().getPath());
        System.out.println("    burnin:     " + burnInPercentageInput.get());
        System.out.println("    CCD type:   " + CCDType.fromName(ccdTypeInput.get()));

        // create base CCD
        MemoryFriendlyTreeSet treeSet = CCDToolUtil.getTreeSet(treeInput, burnInPercentageInput.get());
        AbstractCCD ccd = CCDToolUtil.getCCDTypeByName(treeSet, ccdTypeInput);
        ccd.computeCladeProbabilities();

        // computation
        System.out.println("> computing rogue scores with...");
        int maxCladeSize = Math.max(maxCladeSizeInput.get(), 1);
        maxCladeSize = Math.min(maxCladeSize, ccd.getNumberOfLeaves() - 1);
        double minProbability = minProbabilityInput.get();
        System.out.println("    max clade size: " + maxCladeSize);
        System.out.println("    min clade prob: " + minProbability);
        for (Clade clade : ccd.getClades()) {
            if (cladeConditionSatisfied(clade, maxCladeSize, minProbability)) {
                double rogueScore = RogueDetection.computeCladeRogueScore(ccd, clade, null);
                clade.addData(ROGUE_SCORE_KEY, rogueScore);
            }
        }

        // output scores
        System.out.print("> writing clade rogue scores to file ");
        String separator = separatorInput.get();
        String infoOutputFileName = outputInput.get() + ".csv";
        System.out.println(infoOutputFileName);
        FileWriter fw = new FileWriter(infoOutputFileName);
        BufferedWriter bw = new BufferedWriter(fw);
        String header = "size" + separator + "rogueScore" + separator + "clade";
        bw.write(header);
        bw.newLine();
        for (Clade clade : ccd.getClades()) {
            StringBuilder sb = new StringBuilder();
            if (cladeConditionSatisfied(clade, maxCladeSize, minProbability)) {
                sb.append(clade.size()).append(separator).append(clade.data.get(ROGUE_SCORE_KEY)).append(separator);
                sb.append('"').append(ccd.getTaxaNames(clade.getCladeInBits())).append('"');
                bw.write(sb.toString());
                bw.newLine();
            }
        }
        bw.close();
        fw.close();

        // output MAP tree
        System.out.print("> writing annotated CCD MAP tree to file ");
        infoOutputFileName = outputInput.get() + ".trees";
        System.out.println(infoOutputFileName);
        HeightSettingStrategy hss = HeightSettingStrategy.fromName(heightSettingStrategyInput.get().toString());
        Tree tree = ccd.getMAPTree(hss);
        Map<Clade, Node> map = ccd.getCladeToNodeMap(tree);
        for (Clade clade : ccd.getClades()) {
            if (clade.size() <= maxCladeSize) {
                Node vertex = map.get(clade);
                if (vertex != null) {
                    double rogueScore;
                    if ((clade.data == null) || (clade.data.get(ROGUE_SCORE_KEY) == null)) {
                        rogueScore = RogueDetection.computeCladeRogueScore(ccd, clade, null);
                        clade.addData(ROGUE_SCORE_KEY, rogueScore);
                    } else {
                        rogueScore = (Double) clade.data.get(ROGUE_SCORE_KEY);
                    }

                    String placementInfo = ROGUE_SCORE_KEY + "=" + rogueScore;
                    if (vertex.metaDataString != null) {
                        vertex.metaDataString += "," + placementInfo;
                    } else {
                        vertex.metaDataString = placementInfo;
                    }
                }
            }
        }
        createNexusTreeFile(infoOutputFileName, tree);

        Log.warning("Done");
    }

    private boolean cladeConditionSatisfied(Clade clade, int maxCladeSize, double minProbability) {
        return (clade.size() <= maxCladeSize) && (clade.getProbability() >= minProbability);
    }

    /**
     * Write tree to nexus file.
     *
     * @param outputFileName filename
     * @param tree           to write
     * @throws IOException ...
     */
    public static void createNexusTreeFile(String outputFileName, Tree tree) throws IOException {
        FileWriter fw = new FileWriter(outputFileName);
        BufferedWriter bw = new BufferedWriter(fw);

        bw.write("#NEXUS");
        bw.newLine();
        bw.write("begin trees;");
        bw.newLine();
        bw.write("\ttree CCDMAPtree = [&R] ");
        bw.write(tree.getRoot().toNewick() + ";");
        bw.newLine();
        bw.write("end trees;");
        bw.newLine();

        bw.close();
        fw.close();
    }

    public static void main(String[] args) throws Exception {
        System.out.println(Arrays.toString(args));
        new Application(new RogueAnalysis(), "Rogue Analysis", args);
    }

}
