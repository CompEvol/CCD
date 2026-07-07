package ccd.tools;

import beast.base.core.Input;
import beast.base.core.Log;
import beastfx.app.treeannotator.TreeAnnotator;
import beastfx.app.util.TreeFile;
import ccd.model.AbstractCCD;
import ccd.model.CCD0;
import ccd.model.CCD1;
import ccd.model.CCD2;
import ccd.model.CCDType;
import ccd.model.KRegCCD;
import ccd.model.OptRegCCD;
import ccd.model.RegCCD;

import java.io.IOException;

/**
 * Static methods provided for tools {@link beastfx.app.tools.Application} in the CCD package
 * to reduce redundant code.
 */
public class CCDToolUtil {

    /**
     * Get MemoryFriendlyTreeSet from given filename input and given burnin inputs.
     *
     * @param treeInput        trees file input
     * @param burnInPercentage burnin input checked to be in bounds, altered otherwise
     * @return MemoryFriendlyTreeSet for given inputs
     * @throws IOException ...
     */
    public static TreeAnnotator.MemoryFriendlyTreeSet getTreeSet(Input<TreeFile> treeInput, int burnInPercentage) throws IOException {
        return getTreeSet(treeInput.get().getPath(), burnInPercentage);
    }

    /**
     * Get MemoryFriendlyTreeSet from given filename and given burnin inputs.
     *
     * @param treeFilePath     filename of trees
     * @param burnInPercentage burnin input checked to be in bounds, altered otherwise
     * @return MemoryFriendlyTreeSet for given inputs
     * @throws IOException ...
     */
    public static TreeAnnotator.MemoryFriendlyTreeSet getTreeSet(String treeFilePath, int burnInPercentage) throws IOException {
        int burnin = Math.max(burnInPercentage, 0);
        if (burnin >= 100) {
            Log.warning("Specified burnin input too high - (" + burnin + " >= 100%); set to default of 10%.");
            burnin = 10;
        }
        return new TreeAnnotator().new MemoryFriendlyTreeSet(treeFilePath, burnin);
    }

    /**
     * Get CCD with type specified by input on given treeset.
     *
     * @param treeSet      used to construct CCD
     * @param ccdTypeInput input to specify type
     * @return CCD of given type input on given treeset
     */
    public static AbstractCCD getCCDTypeByName(TreeAnnotator.MemoryFriendlyTreeSet treeSet, Input<String> ccdTypeInput) {
        AbstractCCD ccd;
        CCDType ccdType = CCDType.fromName(ccdTypeInput.get());
        if (ccdType == CCDType.CCD0) {
            ccd = new CCD0(treeSet);
        } else if (ccdType == CCDType.CCD1) {
            ccd = new CCD1(treeSet);
        } else if (ccdType == CCDType.CCD2) {
            ccd = new CCD2(treeSet);
        } else {
            throw new IllegalArgumentException("Illegal CCD type.");
        }
        return ccd;
    }

    /**
     * Materialises all trees of a (burn-in-free) tree set into a list, e.g. for cross-validation
     * that needs to partition the trees.
     *
     * @param treeSet the tree set to read
     * @return all of its trees, in order
     * @throws IOException if the underlying tree file cannot be read
     */
    public static java.util.List<beast.base.evolution.tree.Tree> treesFromSet(
            TreeAnnotator.MemoryFriendlyTreeSet treeSet) throws IOException {
        java.util.List<beast.base.evolution.tree.Tree> trees = new java.util.ArrayList<>();
        treeSet.reset();
        while (treeSet.hasNext()) {
            trees.add(treeSet.next());
        }
        return trees;
    }

    public static AbstractCCD getCCDTypeByName(TreeAnnotator.MemoryFriendlyTreeSet treeSet, CCDType ccdType) {
        AbstractCCD ccd;
        if (ccdType == CCDType.CCD0) {
            ccd = new CCD0(treeSet);
        } else if (ccdType == CCDType.CCD1) {
            ccd = new CCD1(treeSet);
        } else if (ccdType == CCDType.CCD2) {
            ccd = new CCD2(treeSet);
        } else if (ccdType == CCDType.RegCCD) {
            ccd = new RegCCD(treeSet);
        } else if (ccdType == CCDType.OptRegCCD) {
            ccd = new OptRegCCD(treeSet);
        } else if (ccdType == CCDType.KRegCCD) {
            ccd = new KRegCCD(treeSet, KRegCCD.DEFAULT_MU, KRegCCD.DEFAULT_ALPHA,
                    KRegCCD.DEFAULT_RESERVE_DEPTH, KRegCCD.TailMode.BOUND);
        } else {
            throw new IllegalArgumentException("Illegal CCD type.");
        }
        return ccd;
    }
}
