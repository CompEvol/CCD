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
}
