open module ccd {
    requires beast.pkgmgmt;
    requires beast.base;
    requires static beast.fx;
    requires static javafx.controls;
    requires org.apache.commons.math4.legacy;

    exports ccd.algorithms;
    exports ccd.algorithms.credibleSets;
    exports ccd.algorithms.regularisation;
    exports ccd.model;
    exports ccd.model.bitsets;
    exports ccd.tools;

    provides beastfx.app.treeannotator.services.TopologySettingService with
        ccd.tools.CCD0PointEstimate,
        ccd.tools.CCD1PointEstimate,
        ccd.tools.CCD2PointEstimate,
        ccd.tools.CCD0ApproxPointEstimate,
        ccd.tools.HIPSTR;

    provides beast.base.core.BEASTInterface with
        ccd.tools.CCDSampler,
        ccd.tools.EntropyCalculator,
        ccd.tools.SkeletonAnalysis,
        ccd.tools.RogueAnalysis,
        ccd.tools.DissonanceCalculator,
        ccd.tools.TreeCredibleLevel;
}
