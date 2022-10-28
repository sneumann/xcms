# On the long run it would be nice to have all generics in here.
## Alphabetically ordered.

## A
setGeneric("absent", function(object, class, minfrac) standardGeneric("absent"))
setGeneric("absMz", function(object, ...) standardGeneric("absMz"))
setGeneric("absMz<-", function(object, value) standardGeneric("absMz<-"))
setGeneric("absRt", function(object, ...) standardGeneric("absRt"))
setGeneric("absRt<-", function(object, value) standardGeneric("absRt<-"))
setGeneric("addParams", function(object, ...) standardGeneric("addParams"))
setGeneric("addParams<-", function(object, value) standardGeneric("addParams<-"))
setGeneric("addProcessHistory", function(object, ...)
    standardGeneric("addProcessHistory"))

#' @aliases adjustRtime ObiwarpParam-class
#'
#' @title Alignment: Retention time correction methods.
#'
#' @description
#'
#' The `adjustRtime` method(s) perform retention time correction (alignment)
#' between chromatograms of different samples. Alignment is performed by defaul
#' on MS level 1 data. Retention times of spectra from other MS levels, if
#' present, are subsequently adjusted based on the adjusted retention times
#' of the MS1 spectra. Note that calling `adjustRtime` on a *xcms* result object
#' will remove any eventually present previous alignment results as well as
#' any correspondence analysis results.
#'
#' The alignment method can be specified (and configured) using a dedicated
#' `param` argument.
#'
#' Supported `param` objects are:
#'
#' - `ObiwarpParam`: performs retention time adjustment based on the full m/z -
#'   rt data using the *obiwarp* method (Prince (2006)). It is based on the
#'   [original code](http://obi-warp.sourceforge.net) but supports in addition
#'   alignment of multiple samples by aligning each against a *center* sample.
#'   The alignment is performed directly on the [profile-matrix] and can hence
#'   be performed independently of the peak detection or peak grouping.
#'
#' @section Subset-based alignment:
#'
#' All alignment methods allow to perform the retention time correction on a
#' user-selected subset of samples (e.g. QC samples) after which all samples
#' not part of that subset will be adjusted based on the adjusted retention
#' times of the *closest* subset sample (close in terms of index within `object`
#' and hence possibly injection index). It is thus suggested to load MS data
#' files in the order in which their samples were injected in the measurement
#' run(s).
#'
#' How the non-subset samples are adjusted depends also on the parameter
#' `subsetAdjust`: with `subsetAdjust = "previous"`, each non-subset
#' sample is adjusted based on the closest *previous* subset sample which
#' results in most cases with adjusted retention times of the non-subset
#' sample being identical to the subset sample on which the adjustment bases.
#' The second, default, option is `subsetAdjust = "average"` in which case
#' each non subset sample is adjusted based on the average retention time
#' adjustment from the previous and following subset sample. For the average,
#' a weighted mean is used with weights being the inverse of the distance of
#' the non-subset sample to the subset samples used for alignment.
#'
#' See also section *Alignment of experiments including blanks* in the
#' *xcms* vignette for more details.
#'
#' @param binSize \code{numeric(1)} defining the bin size (in mz dimension)
#'     to be used for the \emph{profile matrix} generation. See \code{step}
#'     parameter in \code{\link{profile-matrix}} documentation for more details.
#'
#' @param BPPARAM parallel processing setup. Defaults to `BPPARAM = bpparam()`.
#'     See [bpparam()] for details.
#'
#' @param centerSample \code{integer(1)} defining the index of the center sample
#'     in the experiment. It defaults to
#'     \code{floor(median(1:length(fileNames(object))))}. Note that if
#'     \code{subset} is used, the index passed with \code{centerSample} is
#'     within these subset samples.
#'
#' @param chunkSize For `adjustRtime` if `object` is either an `MsExperiment` or
#'     `XcmsExperiment`: `integer(1)` defining the number of files (samples)
#'     that should be loaded into memory and processed at the same time.
#'     Alignment is then performed in parallel (per sample) on this subset of
#'     loaded data. This setting thus allows to balance between memory
#'     demand and speed (due to parallel processing). Because parallel
#'     processing can only performed on the subset of data currently loaded
#'     into memory in each iteration, the value for `chunkSize` should match
#'     the defined  parallel setting setup. Using a parallel processing setup
#'     using 4 CPUs (separate processes) but using `chunkSize = `1` will not
#'     perform any parallel processing, as only the data from one sample is
#'     loaded in memory at a time. On the other hand, setting `chunkSize` to
#'     the total number of samples in an experiment will load the full MS data
#'     into memory and will thus in most settings cause an out-of-memory error.
#'
#' @param distFun For `ObiwarpParam`: `character(1)` defining the distance
#'     function to be used. Allowed values are `"cor"` (Pearson's correlation),
#'     `"cor_opt"` (calculate only 10% diagonal band of distance matrix;
#'     better runtime), `"cov"` (covariance), `"prd"` (product) and `"euc"`
#'     (Euclidian distance). The default value is `distFun = "cor_opt"`.
#'
#' @param factorDiag For `ObiwarpParam`: `numeric(1)` defining the local weight
#'     applied to diagonal moves in the alignment.
#'
#' @param factorGap For `ObiwarpParam`: `numeric(1)` defining the local weight
#'     for gap moves in the alignment.
#'
#' @param gapExtend For `ObiwarpParam`: `numeric(1)` defining the penalty for
#'     gap enlargement. The default value for `gapExtend` depends on the value
#'     of `distFun`: for `distFun = "cor"` and `distFun = "cor_opt"` it is
#'     `2.4`, `distFun = "cov"` `11.7`, for `distFun = "euc"` `1.8` and for
#'     `distFun = "prd"` `7.8`.
#'
#' @param gapInit For `ObiwarpParam`: `numeric(1)` defining the penalty for gap
#'     opening. The default value for depends on the value of `distFun`:
#'     `distFun = "cor"` and `distFun = "cor_opt"` it is `0.3`, for
#'     `distFun = "cov"` and `distFun = "prd"` `0.0` and for `distFun = "euc"`
#'     `0.9`.
#'
#' @param initPenalty For `ObiwarpParam`: `numeric(1)` defining the penalty for
#'     initiating an alignment (for local alignment only).
#'
#' @param localAlignment For `ObiwarpParam`: `logical(1)` whether a local
#'     alignment should be performed instead of the default global alignment.
#'
#' @param msLevel For `adjustRtime`: `integer(1)` defining the MS level on
#'     which the alignment should be performed.
#'
#' @param object For `adjustRtime`: an [OnDiskMSnExp()], [XCMSnExp()],
#'     [MsExperiment()] or [XcmsExperiment()] object.
#'
#' @param param The parameter object defining the alignment method (and its
#'     setting).
#'
#' @param response For `ObiwarpParam`: `numeric(1)` defining the
#'     *responsiveness* of warping with `response = 0` giving linear warping on
#'     start and end points and `response = 100` warping using all bijective
#'     anchors.
#'
#' @param subset For `ObiwarpParam`: `integer` with the indices of samples
#'     within the experiment on which the alignment models should be estimated.
#'     Samples not part of the subset are adjusted based on the closest subset
#'     sample. See *Subset-based alignment* section for details.
#'
#' @param subsetAdjust For `ObiwarpParam`: `character(1)` specifying the method
#'     with which non-subset samples should be adjusted. Supported options are
#'     `"previous"` and `"average"` (default). See *Subset-based alignment*
#'     section for details.
#'
#' @param value For all assignment methods: the value to set/replace.
#'
#' @param x An `ObiwarpParam`.
#'
#' @param ... ignored.
#'
#' @return
#'
#' `adjustRtime` on an `OnDiskMSnExp` or `XCMSnExp` object will return an
#' `XCMSnExp` object with the alignment results.
#'
#' `adjustRtime` on an `MsExperiment` or `XcmsExperiment` will return an
#' `XcmsExperiment` with the adjusted retention times stored in an new
#' *spectra variable* `rtime_adjusted` in the object's `spectra`.
#'
#' @name adjustRtime
#'
#' @family retention time correction methods
#'
#' @author Colin Smith, Johannes Rainer
#'
#' @references
#'
#' Prince, J. T., and Marcotte, E. M. (2006) "Chromatographic Alignment of
#' ESI-LC-MS Proteomic Data Sets by Ordered Bijective Interpolated Warping"
#' *Anal. Chem.*, 78 (17), 6140-6152.
#'
#' @md
setGeneric("adjustRtime", function(object, param, ...)
    standardGeneric("adjustRtime"))

setGeneric("adjustedRtime", function(object, ...) standardGeneric("adjustedRtime"))
setGeneric("adjustedRtime<-", function(object, value)
    standardGeneric("adjustedRtime<-"))
setGeneric("ampTh", function(object, ...) standardGeneric("ampTh"))
setGeneric("ampTh<-", function(object, value) standardGeneric("ampTh<-"))
setGeneric("AutoLockMass", function(object) standardGeneric("AutoLockMass"))

## B
setGeneric("baseValue", function(object, ...) standardGeneric("baseValue"))
setGeneric("baseValue<-", function(object, value) standardGeneric("baseValue<-"))
setGeneric("binSize", function(object, ...) standardGeneric("binSize"))
setGeneric("binSize<-", function(object, value) standardGeneric("binSize<-"))
setGeneric("bw", function(object) standardGeneric("bw"))
setGeneric("bw<-", function(object, value) standardGeneric("bw<-"))

## C
setGeneric("calibrate", function(object, ...) standardGeneric("calibrate"))
setGeneric("checkBack", function(object, ...) standardGeneric("checkBack"))
setGeneric("centerSample", function(object) standardGeneric("centerSample"))
setGeneric("centerSample<-", function(object, value)
    standardGeneric("centerSample<-"))
setGeneric("checkBack<-", function(object, value) standardGeneric("checkBack<-"))
setGeneric("chromPeaks", function(object, ...) standardGeneric("chromPeaks"))
setGeneric("chromPeaks<-", function(object, value)
    standardGeneric("chromPeaks<-"))
setGeneric("chromPeakData", function(object, ...) standardGeneric("chromPeakData"))
setGeneric("chromPeakData<-", function(object, value)
    standardGeneric("chromPeakData<-"))
setGeneric("collect", function(object, ...) standardGeneric("collect"))
setGeneric("consecMissedLimit", function(object, ...)
    standardGeneric("consecMissedLimit"))
setGeneric("consecMissedLimit<-", function(object, value)
    standardGeneric("consecMissedLimit<-"))
setGeneric("correlate", function(x, y, ...) standardGeneric("correlate"))
setGeneric("criticalValue", function(object, ...)
    standardGeneric("criticalValue"))
setGeneric("criticalValue<-", function(object, value)
    standardGeneric("criticalValue<-"))

## D
setGeneric("deepCopy", function(object) standardGeneric("deepCopy"))
setGeneric("diffreport", function(object, ...) standardGeneric("diffreport"))
setGeneric("distance", function(object, ...) standardGeneric("distance"))
setGeneric("distance<-", function(object, value) standardGeneric("distance<-"))
setGeneric("distFun", function(object) standardGeneric("distFun"))
setGeneric("distFun<-", function(object, value) standardGeneric("distFun<-"))
setGeneric("dropAdjustedRtime", function(object, ...)
    standardGeneric("dropAdjustedRtime"))
setGeneric("dropChromPeaks", function(object, ...)
    standardGeneric("dropChromPeaks"))
setGeneric("dropFeatureDefinitions", function(object, ...)
    standardGeneric("dropFeatureDefinitions"))
setGeneric("dropFilledChromPeaks", function(object, ...)
    standardGeneric("dropFilledChromPeaks"))

## E
setGeneric("expandMz", function(object, ...)
    standardGeneric("expandMz"))
setGeneric("expandMz<-", function(object, value)
    standardGeneric("expandMz<-"))
setGeneric("expandRt", function(object, ...)
    standardGeneric("expandRt"))
setGeneric("expandRt<-", function(object, value)
    standardGeneric("expandRt<-"))
setGeneric("extraPeaks", function(object, ...)
    standardGeneric("extraPeaks"))
setGeneric("extraPeaks<-", function(object, value)
    standardGeneric("extraPeaks<-"))
setGeneric("extractMsData", function(object, ...)
    standardGeneric("extractMsData"))

## F
setGeneric("factorDiag", function(object) standardGeneric("factorDiag"))
setGeneric("factorDiag<-", function(object, value) standardGeneric("factorDiag<-"))
setGeneric("factorGap", function(object) standardGeneric("factorGap"))
setGeneric("factorGap<-", function(object, value) standardGeneric("factorGap<-"))
setGeneric("family", function(object, ...) standardGeneric("family"))
setGeneric("family<-", function(object, value) standardGeneric("family<-"))
setGeneric("featureDefinitions", function(object, ...)
    standardGeneric("featureDefinitions"))
setGeneric("featureDefinitions<-", function(object, value)
    standardGeneric("featureDefinitions<-"))
setGeneric("featureValues", function(object, ...)
    standardGeneric("featureValues"))
setGeneric("fileIndex", function(object) standardGeneric("fileIndex"))
setGeneric("fileIndex<-", function(object, value) standardGeneric("fileIndex<-"))
setGeneric("filepaths", function(object) standardGeneric("filepaths"))
setGeneric("filepaths<-", function(object, value) standardGeneric("filepaths<-"))
setGeneric("fillChromPeaks", function(object, param, ...)
    standardGeneric("fillChromPeaks"))
setGeneric("fillPeaks.chrom", function(object, ...)
    standardGeneric("fillPeaks.chrom"))
setGeneric("fillPeaks.MSW", function(object, ...)
    standardGeneric("fillPeaks.MSW"))
setGeneric("fillPeaks", function(object, ...) standardGeneric("fillPeaks"))
setGeneric("filterChromPeaks", function(object, ...)
    standardGeneric("filterChromPeaks"))
setGeneric("filterColumnsIntensityAbove", function(object, ...)
    standardGeneric("filterColumnsIntensityAbove"))
setGeneric("filterColumnsKeepTop", function(object, ...)
    standardGeneric("filterColumnsKeepTop"))

#' @aliases findChromPeaks
#'
#' @title Chromatographic Peak Detection
#'
#' @description
#'
#' The `findChromPeaks` method performs chromatographic peak detection on
#' LC/GC-MS data. The peak detection algorithm can be selected, and configured,
#' using the `param` argument.
#'
#' Supported `param` objects are:
#'
#' - [CentWaveParam()]: chromatographic peak detection using the *centWave*
#'   algorithm.
#'
#' - [CentWavePredIsoParam()]: *centWave* with predicted isotopes. Peak
#'   detection uses a two-step centWave-based approach considering also feature
#'   isotopes.
#'
#' - [MatchedFilterParam()]: peak detection using the *matched filter*
#'   algorithm.
#'
#' - [MassifquantParam()]: peak detection using the Kalman filter-based
#'   *massifquant* method.
#'
#' - [MSWParam()]: single-spectrum non-chromatography MS data peak detection.
#'
#' For specific examples see the help pages of the individual parameter classes
#' listed above.
#'
#' @param add `logical(1)` (if `object` contains already chromatographic peaks,
#'     i.e. is either an `XCMSnExp` or `XcmsExperiment`) whether chromatographic
#'     peak detection results should be **added** to existing results. By
#'     default (`add = FALSE`) any additional `findChromPeaks` call on a result
#'     object will remove previous results.
#'
#' @param BPPARAM Parallel processing setup. Uses by default the system-wide
#'     default setup. See [bpparam()] for more details.
#'
#' @param chunkSize `integer(1)` for `object` being an `MsExperiment` or
#'     [XcmsExperiment()]: defines the number of files (samples) for which the
#'     full peaks data (m/z and intensity values) should be loaded into memory
#'     at the same time. Peak detection is then performed in parallel (per
#'     sample) on this subset of loaded data. This setting thus allows to
#'     balance between memory demand and speed (due to parallel processing) of
#'     the peak detection. Because parallel processing can only performed on
#'     the subset of data loaded currently into memory (in each iteration), the
#'     value for `chunkSize` should be match the defined  parallel setting
#'     setup. Using a parallel processing setup using 4 CPUs (separate
#'     processes) but using `chunkSize = `1` will not perform any parallel
#'     processing, as only the data from one sample is loaded in memory at a
#'     time. On the other hand, setting `chunkSize` to the total number of
#'     samples in an experiment will load the full MS data into memory and
#'     will thus in most settings cause an out-of-memory error.
#'     By setting `chunkSize = -1` the peak detection will be performed
#'     separately, and in parallel, for each sample. This will however not work
#'     for all `Spectra` *backends* (see eventually [Spectra()] for details).
#'
#' @param msLevel `integer(1)` defining the MS level on which the
#'     chromatographic peak detection should be performed.
#'
#' @param object The data object on which to perform the peak detection. Can be
#'     an [OnDiskMSnExp()], [XCMSnExp()], [MChromatograms()] or [MsExperiment()]
#'     object.
#'
#' @param param The parameter object selecting and configuring the algorithm.
#'
#' @param ... Optional parameters.
#'
#' @name findChromPeaks
#'
#' @family peak detection methods
#'
#' @seealso
#'
#' [plotChromPeaks()] to plot identified chromatographic peaks for one file.
#'
#' [refineChromPeaks()] for methods to *refine* or clean identified
#' chromatographic peaks.
#'
#' [manualChromPeaks()] to manually add/define chromatographic peaks.
#'
#' @author Johannes Rainer
#'
#' @md
setGeneric("findChromPeaks", function(object, param, ...)
           standardGeneric("findChromPeaks"))

setGeneric("findMZ", function(object, find, ppmE=25, print=TRUE)
    standardGeneric("findMZ"))
setGeneric("findmzROI", function(object, ...) standardGeneric("findmzROI"))
setGeneric("findneutral", function(object, find, ppmE=25, print=TRUE)
    standardGeneric("findneutral"))
setGeneric("findKalmanROI", function(object, ...)
    standardGeneric("findKalmanROI"))
setGeneric("findPeaks", function(object, ...) standardGeneric("findPeaks"))
setGeneric("findPeaks.centWave", function(object, ...)
    standardGeneric("findPeaks.centWave"))
setGeneric("findPeaks.addPredictedIsotopeFeatures", function(object, ...)
    standardGeneric("findPeaks.addPredictedIsotopeFeatures"))
setGeneric("findPeaks.centWaveWithPredictedIsotopeROIs", function(object, ...)
    standardGeneric("findPeaks.centWaveWithPredictedIsotopeROIs"))
setGeneric("findPeaks.massifquant", function(object, ...)
    standardGeneric("findPeaks.massifquant"))
setGeneric("findPeaks.matchedFilter", function(object, ...)
    standardGeneric("findPeaks.matchedFilter"))
setGeneric("findPeaks.MSW", function(object, ...)
    standardGeneric("findPeaks.MSW"))
setGeneric("findPeaks.MS1", function(object, ...)
    standardGeneric("findPeaks.MS1"))
setGeneric("firstBaselineCheck", function(object, ...)
    standardGeneric("firstBaselineCheck"))
setGeneric("firstBaselineCheck<-", function(object, value)
    standardGeneric("firstBaselineCheck<-"))
setGeneric("fitgauss", function(object, ...) standardGeneric("fitgauss"))
setGeneric("fitgauss<-", function(object, value) standardGeneric("fitgauss<-"))
setGeneric("fwhm", function(object, ...) standardGeneric("fwhm"))
setGeneric("fwhm<-", function(object, value) standardGeneric("fwhm<-"))

## G
setGeneric("gapExtend", function(object) standardGeneric("gapExtend"))
setGeneric("gapExtend<-", function(object, value) standardGeneric("gapExtend<-"))
setGeneric("gapInit", function(object) standardGeneric("gapInit"))
setGeneric("gapInit<-", function(object, value) standardGeneric("gapInit<-"))
setGeneric("getEIC", function(object, ...) standardGeneric("getEIC"))
setGeneric("getMsnScan", function(object, ...) standardGeneric("getMsnScan"))
setGeneric("getPeaks", function(object, ...) standardGeneric("getPeaks"))
setGeneric("getScan", function(object, ...) standardGeneric("getScan"))
setGeneric("getSpec", function(object, ...) standardGeneric("getSpec"))
setGeneric("getXcmsRaw", function(object, ...) standardGeneric("getXcmsRaw"))
## There's too many "group" methods here...
setGeneric("group.density", function(object, ...) standardGeneric("group.density"))
setGeneric("group.mzClust", function(object, ...) standardGeneric("group.mzClust"))
setGeneric("group.nearest", function(object, ...) standardGeneric("group.nearest"))
setGeneric("group", function(object, ...) standardGeneric("group"))


#' @aliases groupChromPeaks
#'
#' @title Correspondence: group chromatographic peaks across samples
#'
#' @description
#'
#' The `groupChromPeaks` method performs a correspondence analysis i.e., it
#' groups chromatographic peaks across samples to define the LC-MS *features*.
#' The correspondence algorithm can be selected, and configured, using the
#' `param` argument.
#'
#' Supported `param` objects are:
#'
#' - [PeakDensityParam()]: correspondence using the *peak density* method.
#'
#' For specific examples and description of the method and settings see the
#' help pages of the individual parameter classes listed above.
#'
#' @param add `logical(1)` (if `object` contains already chromatographic peaks,
#'     i.e. is either an `XCMSnExp` or `XcmsExperiment`) whether chromatographic
#'     peak detection results should be **added** to existing results. By
#'     default (`add = FALSE`) any additional `findChromPeaks` call on a result
#'     object will remove previous results.
#'
#' @param msLevel `integer(1)` defining the MS level on which the
#'     chromatographic peak detection should be performed.
#'
#' @param object The data object on which the correspondence analysis should be
#'     performed. Can be an [XCMSnExp()], [XcmsExperiment()] object.
#'
#' @param param The parameter object selecting and configuring the algorithm.
#'
#' @param ... Optional parameters.
#'
#' @name groupChromPeaks
#'
#' @family peak grouping methods
#'
#' @author Johannes Rainer
#'
#' @md
setGeneric("groupChromPeaks", function(object, param, ...)
           standardGeneric("groupChromPeaks"))

setGeneric("groupidx", function(object) standardGeneric("groupidx"))
setGeneric("groupidx<-", function(object, value) standardGeneric("groupidx<-"))
setGeneric("groupnames", function(object, ...) standardGeneric("groupnames"))
setGeneric("groups", function(object) standardGeneric("groups"))
setGeneric("groups<-", function(object, value) standardGeneric("groups<-"))
setGeneric("groupval", function(object, ...) standardGeneric("groupval"))

## H
setGeneric("hasMSn", function(object, ...) standardGeneric("hasMSn"))
setGeneric("hasAdjustedRtime", function(object, ...)
    standardGeneric("hasAdjustedRtime"))
setGeneric("hasFeatures", function(object, ...)
    standardGeneric("hasFeatures"))
setGeneric("hasFilledChromPeaks", function(object, ...)
    standardGeneric("hasFilledChromPeaks"))
setGeneric("hasChromPeaks", function(object, ...)
    standardGeneric("hasChromPeaks"))
setGeneric("hasFilledChromPeaks", function(object, ...)
    standardGeneric("hasFilledChromPeaks"))


## I
setGeneric("image", function(x, ...) standardGeneric("image"))
setGeneric("impute<-", function(object, value) standardGeneric("impute<-"))
setGeneric("index", function(object, ...) standardGeneric("index"))
setGeneric("index<-", function(object, value) standardGeneric("index<-"))
setGeneric("integrate")
##setGeneric("integrate", function(object, ...) standardGeneric("integrate"))
setGeneric("integrate<-", function(object, value) standardGeneric("integrate<-"))
setGeneric("initPenalty", function(object) standardGeneric("initPenalty"))
setGeneric("initPenalty<-", function(object, value) standardGeneric("initPenalty<-"))

## K
setGeneric("kNN", function(object, ...) standardGeneric("kNN"))
setGeneric("kNN<-", function(object, value) standardGeneric("kNN<-"))


## L
setGeneric("levelplot", function(x, data, ...) standardGeneric("levelplot"))
setGeneric("localAlignment", function(object) standardGeneric("localAlignment"))
setGeneric("localAlignment<-", function(object, value) standardGeneric("localAlignment<-"))
setGeneric("loadRaw", function(object, ...) standardGeneric("loadRaw"))

## M
##setGeneric("max", function(x, ...) standardGeneric("max"))
setGeneric("max")
setGeneric("max<-", function(object, value) standardGeneric("max<-"))
setGeneric("maxCharge", function(object) standardGeneric("maxCharge"))
setGeneric("maxCharge<-", function(object, value) standardGeneric("maxCharge<-"))
setGeneric("maxFeatures", function(object) standardGeneric("maxFeatures"))
setGeneric("maxFeatures<-", function(object, value) standardGeneric("maxFeatures<-"))
setGeneric("maxIso", function(object) standardGeneric("maxIso"))
setGeneric("maxIso<-", function(object, value) standardGeneric("maxIso<-"))
setGeneric("makeacqNum", function(object, freq, start=1) standardGeneric("makeacqNum"))
setGeneric("minFraction", function(object) standardGeneric("minFraction"))
setGeneric("minFraction<-", function(object, value) standardGeneric("minFraction<-"))
setGeneric("minNoiseLevel", function(object, ...) standardGeneric("minNoiseLevel"))
setGeneric("minNoiseLevel<-", function(object, value)
    standardGeneric("minNoiseLevel<-"))
setGeneric("minSamples", function(object) standardGeneric("minSamples"))
setGeneric("minSamples<-", function(object, value) standardGeneric("minSamples<-"))
setGeneric("mslevel", function(object, ...) standardGeneric("mslevel"))
setGeneric("mslevel<-", function(object, value) standardGeneric("mslevel<-"))
setGeneric("msnparent2ms", function(object, ...) standardGeneric("msnparent2ms"))
setGeneric("msn2ms", function(object, ...) standardGeneric("msn2ms"))
setGeneric("mzdiff", function(object, ...) standardGeneric("mzdiff"))
setGeneric("mzdiff<-", function(object, value) standardGeneric("mzdiff<-"))
setGeneric("mzrange", function(object, ...) standardGeneric("mzrange"))
setGeneric("mzCenterFun", function(object, ...) standardGeneric("mzCenterFun"))
setGeneric("mzCenterFun<-", function(object, value)
    standardGeneric("mzCenterFun<-"))
setGeneric("mzIntervalExtension", function(object, ...)
    standardGeneric("mzIntervalExtension"))
setGeneric("mzIntervalExtension<-", function(object, value)
    standardGeneric("mzIntervalExtension<-"))
setGeneric("mzVsRtBalance", function(object, ...)
    standardGeneric("mzVsRtBalance"))
setGeneric("mzVsRtBalance<-", function(object, value)
    standardGeneric("mzVsRtBalance<-"))

## N
setGeneric("nearbyPeak", function(object, ...) standardGeneric("nearbyPeak"))
setGeneric("nearbyPeak<-", function(object, value) standardGeneric("nearbyPeak<-"))
setGeneric("noise", function(object, ...) standardGeneric("noise"))
setGeneric("noise<-", function(object, value) standardGeneric("noise<-"))

## P
setGeneric("peakGroupsMatrix", function(object, ...)
    standardGeneric("peakGroupsMatrix"))
setGeneric("peakGroupsMatrix<-", function(object, value)
    standardGeneric("peakGroupsMatrix<-"))
setGeneric("peakScaleRange", function(object, ...)
    standardGeneric("peakScaleRange"))
setGeneric("peakScaleRange<-", function(object, value)
    standardGeneric("peakScaleRange<-"))
setGeneric("peakTable", function(object, ...) standardGeneric("peakTable"))
setGeneric("peakThr", function(object, ...) standardGeneric("peakThr"))
setGeneric("peakThr<-", function(object, value) standardGeneric("peakThr<-"))
setGeneric("peakwidth", function(object, ...) standardGeneric("peakwidth"))
setGeneric("peakwidth<-", function(object, value) standardGeneric("peakwidth<-"))
setGeneric("plotChrom", function(object, ...) standardGeneric("plotChrom"))
setGeneric("plotChromPeakDensity", function(object, ...)
    standardGeneric("plotChromPeakDensity"))
setGeneric("plotChromatogramsOverlay", function(object, ...)
    standardGeneric("plotChromatogramsOverlay"))
setGeneric("plotEIC", function(object, ...) standardGeneric("plotEIC"))
setGeneric("plotPeaks", function(object, ...) standardGeneric("plotPeaks"))
setGeneric("plotRaw", function(object, ...) standardGeneric("plotRaw"))
setGeneric("plotrt", function(object, ...) standardGeneric("plotrt"))
setGeneric("plotScan", function(object, ...) standardGeneric("plotScan"))
setGeneric("plotSpec", function(object, ...) standardGeneric("plotSpec"))
setGeneric("plotSurf", function(object, ...) standardGeneric("plotSurf"))
setGeneric("plotTIC", function(object, ...) standardGeneric("plotTIC"))
setGeneric("plotTree", function(object, ...) standardGeneric("plotTree"))
setGeneric("ppm", function(object, ...) standardGeneric("ppm"))
setGeneric("ppm<-", function(object, value) standardGeneric("ppm<-"))
setGeneric("prefilter", function(object, ...) standardGeneric("prefilter"))
setGeneric("prefilter<-", function(object, value) standardGeneric("prefilter<-"))
setGeneric("present", function(object, class, minfrac) standardGeneric("present"))
setGeneric("processDate", function(object, ...) standardGeneric("processDate"))
setGeneric("processDate<-", function(object, value) standardGeneric("processDate<-"))
setGeneric("processInfo", function(object, ...) standardGeneric("processInfo"))
setGeneric("processInfo<-", function(object, value) standardGeneric("processInfo<-"))
setGeneric("processParam", function(object, ...) standardGeneric("processParam"))
setGeneric("processParam<-", function(object, value)
    standardGeneric("processParam<-"))
setGeneric("processType", function(object, ...) standardGeneric("processType"))
setGeneric("processType<-", function(object, value) standardGeneric("processType<-"))
setGeneric("processHistory", function(object, ...) standardGeneric("processHistory"))
setGeneric("profinfo", function(object) standardGeneric("profinfo"))
setGeneric("profinfo<-", function(object, value) standardGeneric("profinfo<-"))
setGeneric("profMat", function(object, ...) standardGeneric("profMat"))
setGeneric("profMedFilt", function(object, ...) standardGeneric("profMedFilt"))
setGeneric("profMethod", function(object) standardGeneric("profMethod"))
setGeneric("profMethod<-", function(object, value) standardGeneric("profMethod<-"))
setGeneric("profRange", function(object, ...) standardGeneric("profRange"))
setGeneric("profStep", function(object) standardGeneric("profStep"))
setGeneric("profStep<-", function(object, value) standardGeneric("profStep<-"))
setGeneric("profStepPad<-", function(object, value) standardGeneric("profStepPad<-"))
setGeneric("profMz", function(object) standardGeneric("profMz"))
setGeneric("progressCallback", function(object) standardGeneric("progressCallback"))
setGeneric("progressCallback<-", function(object, value) standardGeneric("progressCallback<-"))
setGeneric("progressInfoUpdate", function(object) standardGeneric("progressInfoUpdate"))

## R
setGeneric("rawEIC", function(object, ...) standardGeneric("rawEIC"))
setGeneric("rawMat", function(object, ...) standardGeneric("rawMat"))
setGeneric("rawMZ", function(object, ...) standardGeneric("rawMZ"))
setGeneric("refineChromPeaks", function(object, param, ...)
    standardGeneric("refineChromPeaks"))
setGeneric("removeIntensity", function(object, ...) standardGeneric("removeIntensity"))
setGeneric("response", function(object) standardGeneric("response"))
setGeneric("response<-", function(object, value) standardGeneric("response<-"))
setGeneric("retcor", function(object, ...) standardGeneric("retcor"))
setGeneric("retcor.peakgroups", function(object, ...) standardGeneric("retcor.peakgroups"))
setGeneric("retcor.obiwarp", function(object, ...) standardGeneric("retcor.obiwarp"))
setGeneric("revMz", function(object, ...) standardGeneric("revMz"))
setGeneric("ridgeLength", function(object, ...) standardGeneric("ridgeLength"))
setGeneric("ridgeLength<-", function(object, value) standardGeneric("ridgeLength<-"))
setGeneric("roiList", function(object, ...) standardGeneric("roiList"))
setGeneric("roiList<-", function(object, value) standardGeneric("roiList<-"))
setGeneric("roiScales", function(object, ...) standardGeneric("roiScales"))
setGeneric("roiScales<-", function(object, value) standardGeneric("roiScales<-"))
setGeneric("rtrange", function(object) standardGeneric("rtrange"))

## S
setGeneric("sampclass", function(object) standardGeneric("sampclass"))
setGeneric("sampclass<-", function(object, value) standardGeneric("sampclass<-"))
setGeneric("sampleGroups", function(object) standardGeneric("sampleGroups"))
setGeneric("sampleGroups<-", function(object, value)
    standardGeneric("sampleGroups<-"))
setGeneric("sampnames", function(object) standardGeneric("sampnames"))
setGeneric("sampnames<-", function(object, value) standardGeneric("sampnames<-"))
setGeneric("scales", function(object, ...) standardGeneric("scales"))
setGeneric("scales<-", function(object, value) standardGeneric("scales<-"))
setGeneric("scanrange", function(object, ...) standardGeneric("scanrange"))
setGeneric("scanrange<-", function(object, value) standardGeneric("scanrange<-"))
setGeneric("sigma", function(object, value) standardGeneric("sigma"))
setGeneric("sigma<-", function(object, value) standardGeneric("sigma<-"))
setGeneric("showError", function(object, ...) standardGeneric("showError"))
setGeneric("smooth<-", function(object, value) standardGeneric("smooth<-"))
setGeneric("snthresh", function(object, ...) standardGeneric("snthresh"))
setGeneric("snthresh<-", function(object, value) standardGeneric("snthresh<-"))
setGeneric("snthreshIsoROIs", function(object, ...)
    standardGeneric("snthreshIsoROIs"))
setGeneric("snthreshIsoROIs<-", function(object, value)
    standardGeneric("snthreshIsoROIs<-"))
setGeneric("sortMz", function(object, ...) standardGeneric("sortMz"))
setGeneric("span", function(object, ...) standardGeneric("span"))
setGeneric("span<-", function(object, value) standardGeneric("span<-"))
setGeneric("specDist", function(object, ...) standardGeneric("specDist"))
setGeneric("specDist.meanMZmatch",
           function(peakTable1, peakTable2, matchdist=1, matchrate=1,
                    mzabs=0.001, mzppm=10, symmetric=TRUE)
               standardGeneric("specDist.meanMZmatch"))
setGeneric("specDist.cosine",
           function(peakTable1, peakTable2, mzabs = 0.001, mzppm = 10,
                    mzExp = 0.6, intExp = 3, nPdiff = 2, nPmin = 8,
                    symmetric = FALSE)
               standardGeneric("specDist.cosine"))
setGeneric("specDist.peakCount",
           function(peakTable1, peakTable2, mzabs=0.001, mzppm=10,symmetric=FALSE)
               standardGeneric("specDist.peakCount"))
setGeneric("steps", function(object, ...) standardGeneric("steps"))
setGeneric("steps<-", function(object, value) standardGeneric("steps<-"))
setGeneric("stitch", function(object, lockMass, ...) standardGeneric("stitch"))
setGeneric("stitch.xml", function(object, lockMass) standardGeneric("stitch.xml"))
setGeneric("stitch.netCDF", function(object, lockMass) standardGeneric("stitch.netCDF"))
setGeneric("stitch.netCDF.new", function(object, lockMass) standardGeneric("stitch.netCDF.new"))
setGeneric("subset<-", function(object, value) standardGeneric("subset<-"))
setGeneric("subsetAdjust", function(object, ...) standardGeneric("subsetAdjust"))
setGeneric("subsetAdjust<-", function(object, value) standardGeneric("subsetAdjust<-"))


## T
setGeneric("tuneIn", function(object, ...) standardGeneric("tuneIn"))
setGeneric("tuneIn<-", function(object, value) standardGeneric("tuneIn<-"))

## U
setGeneric("unions", function(object, ...) standardGeneric("unions"))
setGeneric("unions<-", function(object, value) standardGeneric("unions<-"))

## V
setGeneric("verboseColumns", function(object, ...) standardGeneric("verboseColumns"))
setGeneric("verboseColumns<-", function(object, value)
    standardGeneric("verboseColumns<-"))

## W
setGeneric("withWave", function(object, ...) standardGeneric("withWave"))
setGeneric("withWave<-", function(object, value) standardGeneric("withWave<-"))
setGeneric("write.cdf", function(object, ...) standardGeneric("write.cdf"))
setGeneric("write.mzdata", function(object, ...) standardGeneric("write.mzdata"))
setGeneric("write.mzQuantML", function(object, ...) standardGeneric("write.mzQuantML"))

## X
setGeneric("xcmsSource", function(object, ...) standardGeneric("xcmsSource"))
