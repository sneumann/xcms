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

#' @aliases adjustRtime ObiwarpParam-class PeakGroupsParam-class
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
#' - `PeakGroupsParam`: performs retention time correctoin based on the
#'   alignment of features defined in all/most samples (corresponding to
#'   *house keeping compounds* or marker compounds) (Smith 2006). First the
#'   retention time deviation of these features is described by fitting either a
#'   polynomial (`smooth = "loess"`) or a linear (`smooth = "linear"`) function
#'   to the data points. These are then subsequently used to adjust the
#'   retention time of each spectrum in each sample (even from spectra of
#'   MS levels different than MS 1). Since the function is
#'   based on features (i.e. chromatographic peaks grouped across samples) a
#'   initial correspondence analysis has to be performed **before** using the
#'   [groupChromPeaks()] function. Alternatively, it is also possible to
#'   manually define a `numeric` matrix with retention times of markers in each
#'   samples that should be used for alignment. Such a `matrix` can be passed
#'   to the alignment function using the `peakGroupsMatrix` parameter of the
#'   `PeakGroupsParam` parameter object. By default the `adjustRtimePeakGroups`
#'   function is used to define this `matrix`. This function identifies peak
#'   groups (features) for alignment in `object` based on the parameters defined
#'   in `param`. See also [do_adjustRtime_peakGroups()] for the core API
#'   function.
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
#' @param extraPeaks For `PeakGroupsParam`: `numeric(1)` defining the maximal
#'     number of additional peaks for all samples to be assigned to a peak
#'     group (feature) for retention time correction. For a data set with 6
#'     samples, `extraPeaks = 1` uses all peak groups with a total peak count
#'     `<= 6 + 1`. The total peak count is the total number of peaks being
#'     assigned to a peak group and considers also multiple peaks within a
#'     sample that are assigned to the group.
#'
#' @param factorDiag For `ObiwarpParam`: `numeric(1)` defining the local weight
#'     applied to diagonal moves in the alignment.
#'
#' @param factorGap For `ObiwarpParam`: `numeric(1)` defining the local weight
#'     for gap moves in the alignment.
#'
#' @param family For `PeakGroupsParam`: `character(1)` defining the method for
#'     loess smoothing. Allowed values are `"gaussian"` and `"symmetric"`. See
#'     [loess()] for more information.
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
#' @param minFraction For `PeakGroupsParam`: `numeric(1)` between 0 and 1
#'     defining the minimum required fraction of samples in which peaks for
#'     the peak group were identified. Peak groups passing this criteria will
#'     be aligned across samples and retention times of individual spectra will
#'     be adjusted based on this alignment. For `minFraction = 1` the peak
#'     group has to contain peaks in all samples of the experiment. Note that if
#'     `subset` is provided, the specified fraction is relative to the
#'     defined subset of samples and not to the total number of samples within
#'     the experiment (i.e. a peak has to be present in the specified
#'     proportion of subset samples).
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
#' @param peakGroupsMatrix For `PeakGroupsParam`: optional `matrix` of (raw)
#'     retention times for the (marker) peak groups on which the alignment
#'     should be performed. Each column represents a sample, each row a
#'     feature/peak group. The `adjustRtimePeakGroups` method is used by
#'     default to determine this matrix on the provided `object`.
#'
#' @param response For `ObiwarpParam`: `numeric(1)` defining the
#'     *responsiveness* of warping with `response = 0` giving linear warping on
#'     start and end points and `response = 100` warping using all bijective
#'     anchors.
#'
#' @param smooth For `PeakGroupsParam`: `character(1)` defining the function to
#'     be used to interpolate corrected retention times for all peak groups.
#'     Can be either `"loess"` or `"linear"`.
#'
#' @param span For `PeakGroupsParam`: `numeric(1)` defining the degree of
#'     smoothing (if `smooth = "loess"`). This parameter is passed to the
#'     internal call to [loess()].
#'
#' @param subset For `ObiwarpParam` and `PeakGroupsParam`: `integer` with the
#'     indices of samples within the experiment on which the alignment models
#'     should be estimated.
#'     Samples not part of the subset are adjusted based on the closest subset
#'     sample. See *Subset-based alignment* section for details.
#'
#' @param subsetAdjust For `ObiwarpParam` and `PeakGroupsParam`: `character(1)`
#'     specifying the method with which non-subset samples should be adjusted.
#'     Supported options are `"previous"` and `"average"` (default).
#'     See *Subset-based alignment* section for details.
#'
#' @param value For all assignment methods: the value to set/replace.
#'
#' @param x An `ObiwarpParam` or `PeakGroupsParam` object.
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
#' `ObiwarpParam` and `PeakGroupsParam` return the respective parameter object.
#'
#' `adjustRtimeGroups` returns a `matrix` with the retention times of *marker*
#' features in each sample (each row one feature, each row one sample).
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
#' Smith, C.A., Want, E.J., O'Maille, G., Abagyan, R. and Siuzdak, G. (2006).
#' "XCMS: Processing Mass Spectrometry Data for Metabolite Profiling Using
#' Nonlinear Peak Alignment, Matching, and Identification" *Anal. Chem.*
#' 78:779-787.
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


#' @title Gap Filling
#'
#' @aliases fillChromPeaks
#'
#' @description
#'
#' Gap filling integrate signal in the m/z-rt area of a feature (i.e., a
#' chromatographic peak group) for samples in which no chromatographic
#' peak for this feature was identified and add it to the [chromPeaks()]
#' matrix. Such *filled-in* peaks are indicated with a `TRUE` in column
#' `"is_filled"` in the result object's [chromPeakData()] data frame.
#'
#' The method for gap filling along with its settings can be defined with
#' the `param` argument. Two different approaches are available:
#'
#' - `param = FillChromPeaksParam()`: the default of the original `xcms`
#'   code. Signal is integrated from the m/z and retention time range as
#'   defined in the [featureDefinitions()] data frame, i.e. from the
#'   `"rtmin"`, `"rtmax"`, `"mzmin"` and `"mzmax"`. This method is not
#'   suggested as it underestimates the actual peak area and it is also
#'   not available for `object` being an [XcmsExperiment] object. See
#'   details below for more information and settings for this method.
#'
#' - `param = ChromPeakAreaParam()`: the area from which the signal for a
#'   feature is integrated is defined based on the feature's chromatographic
#'   peak areas. The m/z range is by default defined as the the lower quartile
#'   of chromatographic peaks' `"mzmin"` value to the upper quartile of the
#'   chromatographic peaks' `"mzmax"` values. The retention time range for the
#'   area is defined analogously. Alternatively, by setting `mzmin = median`,
#'   `mzmax = median`, `rtmin = median` and `rtmax = median` in
#'   `ChromPeakAreaParam`, the median `"mzmin"`, `"mzmax"`, `"rtmin"` and
#'   `"rtmax"` values from all detected chromatographic peaks of a feature
#'   would be used instead.
#'   In contrast to the  `FillChromPeaksParam` approach this method uses the
#'   actual identified chromatographic peaks of a feature to define the area
#'   from which the signal should be integrated.
#'
#' @details
#'
#' After correspondence (i.e. grouping of chromatographic peaks across
#' samples) there will always be features (peak groups) that do not include
#' peaks from every sample. The `fillChromPeaks` method defines
#' intensity values for such features in the missing samples by integrating
#' the signal in the m/z-rt region of the feature. Two different approaches
#' to define this region are available: with `ChromPeakAreaParam` the region
#' is defined based on the detected **chromatographic peaks** of a feature,
#' while with `FillChromPeaksParam` the region is defined based on the m/z and
#' retention times of the **feature** (which represent the m/z and retentention
#' times of the apex position of the associated chromatographic peaks). For the
#' latter approach various parameters are available to increase the area from
#' which signal is to be integrated, either by a constant value (`fixedMz` and
#' `fixedRt`) or by a feature-relative amount (`expandMz` and `expandRt`).
#'
#' Adjusted retention times will be used if available.
#'
#' Based on the peak finding algorithm that was used to identify the
#' (chromatographic) peaks, different internal functions are used to
#' guarantee that the integrated peak signal matches as much as possible
#' the peak signal integration used during the peak detection. For peaks
#' identified with the [matchedFilter()] method, signal
#' integration is performed on the *profile matrix* generated with
#' the same settings used also during peak finding (using the same
#' `bin` size for example). For direct injection data and peaks
#' identified with the `MSW` algorithm signal is integrated
#' only along the mz dimension. For all other methods the complete (raw)
#' signal within the area is used.
#'
#' @note
#'
#' The reported `"mzmin"`, `"mzmax"`, `"rtmin"` and
#' `"rtmax"` for the filled peaks represents the actual MS area from
#' which the signal was integrated.
#'
#' No peak is filled in if no signal was present in a file/sample
#' in the respective mz-rt area. These samples will still show a `NA`
#' in the matrix returned by the [featureValues()] method.
#'
#' @param chunkSize For `fillChromPeaks` if `object` is an `XcmsExperiment`:
#'     `integer(1)` defining the number of files (samples)
#'     that should be loaded into memory and processed at the same time.
#'     This setting thus allows to balance between memory
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
#' @param expandMz for `FillChromPeaksParam`: `numeric(1)` defining the
#'     value by which the mz width of peaks should be expanded. Each peak is
#'     expanded in mz direction by `expandMz *` their original m/z width.
#'     A value of `0` means no expansion, a value of `1` grows each peak
#'     by `1 *` the m/z width of the peak resulting in peaks with twice
#'     their original size in m/z direction (expansion by half m/z width
#'     to both sides).
#'
#' @param expandRt for `FillChromPeaksParam`: `numeric(1)`, same as
#'     `expandMz` but for the retention time width.
#'
#' @param fixedMz for `FillChromPeaksParam`: `numeric(1)` defining a constant
#'     factor by which the m/z width of each feature is to be expanded.
#'     The m/z width is expanded on both sides by `fixedMz` (i.e. `fixedMz`
#'     is subtracted from the lower m/z and added to the upper m/z). This
#'     expansion is applied *after* `expandMz` and `ppm`.
#'
#' @param fixedRt for `FillChromPeaksParam`: `numeric(1)` defining a constant
#'     factor by which the retention time width of each factor is to be
#'     expanded. The rt width is expanded on both sides by `fixedRt` (i.e.
#'     `fixedRt` is subtracted from the lower rt and added to the upper rt).
#'     This expansion is applied *after* `expandRt`.
#'
#' @param msLevel `integer(1)` defining the MS level on which peak filling
#'     should be performed (defaults to `msLevel = 1L`). Only peak filling
#'     on one MS level at a time is supported, to fill in peaks for MS
#'     level 1 and 2 run first using `msLevel = 1` and then (on the returned
#'     result object) again with `msLevel = 2`.
#'
#' @param mzmax `function` to be applied to values in the `"mzmax"` column
#'     of all chromatographic peaks of a feature to define the upper m/z
#'     value of the area from which signal for the feature should be
#'    integrated. Defaults to  `mzmax = function(z) quantile(z, probs = 0.75)`
#'    hence using the 75% quantile of all values.
#'
#' @param mzmin `function` to be applied to values in the `"mzmin"` column
#'     of all chromatographic peaks of a feature to define the lower m/z
#'     value of the area from which signal for the feature should be
#'     integrated. Defaults to `mzmin = function(z) quantile(z, probs = 0.25)`
#'     hence using the 25% quantile of all values.
#'
#' @param object `XcmsExperiment` or `XCMSnExp` object with identified and
#'     grouped chromatographic peaks.
#'
#' @param param `ChromPeakAreaParam` or `FillChromPeaksParam` object
#'     defining which approach should be used (see details section).
#'
#' @param ppm for `FillChromPeaksParam`: `numeric(1)` optionally specifying
#'     a *ppm* by which the m/z width of the peak region should be expanded.
#'     For peaks with an m/z width smaller than
#'     `mean(c(mzmin, mzmax)) * ppm / 1e6`, the `mzmin` will be replaced by
#'     `mean(c(mzmin, mzmax)) - (mean(c(mzmin, mzmax)) * ppm / 2 / 1e6)`
#'     `mzmax` by
#'     `mean(c(mzmin, mzmax)) + (mean(c(mzmin, mzmax)) * ppm / 2 / 1e6)`.
#'     This is applied before eventually expanding the m/z width using the
#'     `expandMz` parameter.
#'
#' @param rtmax `function` to be applied to values in the `"rtmax"` column
#'     of all chromatographic peaks of a feature to define the upper rt
#'     value of the area from which signal for the feature should be
#'     integrated. Defaults to `rtmax = function(z) quantile(z, probs = 0.75)`
#'     hence using the 75% quantile of all values.
#'
#' @param rtmin `function` to be applied to values in the `"rtmin"` column
#'     of all chromatographic peaks of a feature to define the lower rt
#'     value of the area from which signal for the feature should be
#'     integrated. Defaults to `rtmin = function(z) quantile(z, probs = 0.25)`
#'     hence using the 25% quantile of all values.
#'
#' @param BPPARAM Parallel processing settings.
#'
#' @param ... currently ignored.
#'
#' @return
#'
#' An [XcmsExperiment] or  `XCMSnExp` object with previously missing
#' chromatographic peaks for features filled into its [chromPeaks()] matrix.
#'
#' @rdname fillChromPeaks
#'
#' @author Johannes Rainer
#'
#' @seealso [groupChromPeaks()] for methods to perform the correspondence.
#'
#' @seealso [featureArea] for the function to define the m/z-retention time
#'     region for each feature.
#'
#' @md
#'
#' @examples
#'
#' ## Load a test data set with identified chromatographic peaks
#' data(faahko_sub)
#' ## Update the path to the files for the local system
#' dirname(faahko_sub) <- system.file("cdf/KO", package = "faahKO")
#' res <- faahko_sub
#'
#' ## Disable parallel processing for this example
#' register(SerialParam())
#'
#' ## Perform the correspondence. We assign all samples to the same group.
#' res <- groupChromPeaks(res,
#'     param = PeakDensityParam(sampleGroups = rep(1, length(fileNames(res)))))
#'
#' ## For how many features do we lack an integrated peak signal?
#' sum(is.na(featureValues(res)))
#'
#' ## Filling missing peak data using the peak area from identified
#' ## chromatographic peaks.
#' res <- fillChromPeaks(res, param = ChromPeakAreaParam())
#'
#' ## How many missing values do we have after peak filling?
#' sum(is.na(featureValues(res)))
#'
#' ## Get the peaks that have been filled in:
#' fp <- chromPeaks(res)[chromPeakData(res)$is_filled, ]
#' head(fp)
#'
#' ## Get the process history step along with the parameters used to perform
#' ## The peak filling:
#' ph <- processHistory(res, type = "Missing peak filling")[[1]]
#' ph
#'
#' ## The parameter class:
#' ph@param
#'
#' ## It is also possible to remove filled-in peaks:
#' res <- dropFilledChromPeaks(res)
#'
#' sum(is.na(featureValues(res)))
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


#' @aliases groupChromPeaks PeakDensityParam-class
#'
#' @aliases NearestPeaksParam-class MzClustParam-class
#'
#' @title Correspondence: group chromatographic peaks across samples
#'
#' @description
#'
#' The `groupChromPeaks` method performs a correspondence analysis i.e., it
#' groups chromatographic peaks across samples to define the LC-MS *features*.
#' The correspondence algorithm can be selected, and configured, using the
#' `param` argument. See documentation of [XcmsExperiment()] and [XCMSnExp()]
#' for information on how to access and extract correspondence results.
#'
#' The correspondence analysis can be performed on chromatographic peaks of
#' any MS level (if present and if chromatographic peak detection has been
#' performed for that MS level) defining features combining these peaks. The
#' MS level can be selected with the parameter `msLevel`. By default, calling
#' `groupChromPeaks` will remove any previous correspondence results. This can
#' be disabled with `add = TRUE`, which will add newly defined features to
#' already present feature definitions.
#'
#' Supported `param` objects are:
#'
#' - `PeakDensityParam`: correspondence using the *peak density* method
#'   (Smith 2006) that groups chromatographic peaks along the retention time
#'   axis within slices of (partially overlapping) m/z ranges. All peaks (from
#'   the same or from different samples) with their apex position being close
#'   on the retention time axis are grouped into a LC-MS feature. See in
#'   addition [do_groupChromPeaks_density()] for the core API function.
#'
#' - `NearestPeaksParam`: performs peak grouping based on the proximity of
#'   chromatographic peaks from different samples in the m/z - rt space similar
#'   to the correspondence method of *mzMine* (Katajamaa 2006). The method
#'   creates first a *master peak list* consisting of all chromatographic peaks
#'   from the sample with the most detected peaks and iteratively calculates
#'   distances to peaks from the sample with the next most number of peaks
#'   grouping peaks together if their *distance* is smaller than the provided
#'   thresholds.
#'   See in addition [do_groupChromPeaks_nearest()] for the core API function.
#'
#' - `MzClustParam`: performs high resolution peak grouping for
#'   **single spectrum** metabolomics data (Kazmi 2006). This method should
#'   **only** be used for such data as the retention time is not considered
#'   in the correspondence analysis.
#'   See in addition [do_groupPeaks_mzClust()] for the core API function.
#'
#' For specific examples and description of the method and settings see the
#' help pages of the individual parameter classes listed above.
#'
#' @param absMz For `NearestPeaksParam` and `MzClustParam`: `numeric(1)`
#'     maximum tolerated distance for m/z values.
#'
#' @param absRt For `NearestPeaksParam`: `numeric(1)` maximum tolerated
#'     distance for retention times.
#'
#' @param add `logical(1)` (if `object` contains already chromatographic peaks,
#'     i.e. is either an `XCMSnExp` or `XcmsExperiment`) whether chromatographic
#'     peak detection results should be **added** to existing results. By
#'     default (`add = FALSE`) any additional `findChromPeaks` call on a result
#'     object will remove previous results.
#'
#' @param binSize For `PeakDensityParam`: `numeric(1)` defining the size of the
#'     overlapping slices in m/z dimension.
#'
#' @param bw For `PeakDensityParam`: `numeric(1)` defining the bandwidth
#'     (standard deviation ot the smoothing kernel) to be used. This argument
#'     is passed to the [density() method.
#'
#' @param kNN For `NearestPeaksParam`: `integer(1)` representing the number of
#'     nearest neighbors to check.
#'
#' @param maxFeatures For `PeakDensityParam`: `numeric(1)` with the maximum
#'     number of peak groups to be identified in a single mz slice.
#'
#' @param minFraction For `PeakDensityParam`: `numeric(1)` defining the minimum
#'     fraction of samples in at least one sample group in which the peaks
#'     have to be present to be considered as a peak group (feature).
#'
#' @param minSamples For `PeakDensityParam`: `numeric(1)` with the minimum
#'     number of samples in at least one sample group in which the peaks have
#'     to be detected to be considered a peak group (feature).
#'
#' @param msLevel `integer(1)` defining the MS level on which the
#'     chromatographic peak detection should be performed.
#'
#' @param mzVsRtBalance For `NearestPeaksParam`: `numeric(1)` representing the
#'     factor by which m/z values are multiplied before calculating the
#'     (euclician) distance between two peaks.
#'
#' @param object The data object on which the correspondence analysis should be
#'     performed. Can be an [XCMSnExp()], [XcmsExperiment()] object.
#'
#' @param ppm For `MzClustParam`: `numeric(1)` representing the relative m/z
#'     error for the clustering/grouping (in parts per million).
#'
#' @param param The parameter object selecting and configuring the algorithm.
#'
#' @param sampleGroups For `PeakDensityParam`: A vector of the same length than
#'     samples defining the sample group assignments (i.e. which samples
#'     belong to which sample
#'     group). This parameter is mandatory for the `PeakDensityParam`
#'     and has to be provided also if there is no sample grouping in the
#'     experiment (in which case all samples should be assigned to the
#'     same group).
#'
#' @param value Replacement value for `<-` methods.
#'
#' @param ... Optional parameters.
#'
#' @return For `groupChromPeaks`: either an [XcmsExperiment()] or [XCMSnExp()]
#'     object with the correspondence result.
#'
#' @name groupChromPeaks
#'
#' @family peak grouping methods
#'
#' @author Colin Smith, Johannes Rainer
#'
#' @references
#'
#' Smith, C.A., Want E.J., O'Maille G., Abagyan R., and Siuzdak G. (2006)
#' "XCMS: Processing Mass Spectrometry Data for Metabolite Profiling Using
#' Nonlinear Peak Alignment, Matching, and Identification" *Anal. Chem.*
#' 78:779-787.
#'
#' Katajamaa, M., Miettinen, J., Oresic, M. (2006) "MZmine: Toolbox for
#' processing and visualization of mass spectrometry based molecular profile
#' data". *Bioinformatics*, 22:634-636.
#'
#' Kazmi S. A., Ghosh, S., Shin, D., Hill, D.W., and Grant, D.F. (2006)
#' "Alignment of high resolution mass spectra: development of a
#' heuristic approach for metabolomics. *Metabolomics* Vol. 2, No. 2, 75-83.
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

#' @title Refine Identified Chromatographic Peaks
#'
#' @aliases FilterIntensityParam-class show,FilterIntensityParam-method
#'
#' @aliases CleanPeaksParam-class show,CleanPeaksParam-method
#'
#' @aliases MergeNeighboringPeaksParam-class show,MergeNeighboringPeaksParam-method
#'
#' @description
#'
#' The `refineChromPeaks` method performs a post-processing of the
#' chromatographic peak detection step to eventually clean and improve the
#' results. The function can be applied to a [XcmsExperiment()] or [XCMSnExp()]
#' object **after** peak detection with [findChromPeaks()]. The type of peak
#' refinement and cleaning can be defined, along with all its settings, using
#' one of the following parameter objects:
#'
#' - `CleanPeaksParam`: remove chromatographic peaks with a retention time
#'   range larger than the provided maximal acceptable width (`maxPeakwidth`).
#'
#' - `FilterIntensityParam`: remove chromatographic peaks with intensities
#'   below the specified threshold. By default (with `nValues = 1`) values in
#'   the `chromPeaks` matrix are evaluated: all peaks with a value in the
#'   column defined with parameter `value` that are `>=` a threshold (defined
#'   with parameter `threshold`) are retained. If `nValues` is larger than 1,
#'   the individual peak intensities from the raw MS files are evaluated:
#'   chromatographic peaks with at least `nValues` mass peaks `>= threshold`
#'   are retained.
#'
#' - `MergeNeighboringPeaksParam`: peak detection sometimes fails to identify a
#'   chromatographic peak correctly, especially for broad peaks and if the peak
#'   shape is irregular (mostly for HILIC data). In such cases several smaller
#'   peaks are reported. Also, peak detection with *centWave* can result in
#'   partially or completely overlapping peaks. This method aims to reduce
#'   such peak detection artifacts by merging chromatographic peaks that are
#'   overlapping or close in RT and m/z dimension (considering also the measured
#'   signal between them). See section *Details for MergeNeighboringPeaksParam*
#'   for details and a comprehensive description of the approach.
#'
#' `refineChromPeaks` methods will always remove feature definitions, because
#' a call to this method can change or remove identified chromatographic peaks,
#' which may be part of features.
#'
#' @section Details for MergeNeighboringPeaksParam:
#'
#' For peak refinement using the `MergeNeighboringPeaksParam`, chromatographic
#' peaks are first expanded in m/z and retention time dimension (based on
#' parameters `expandMz`, `ppm` and `expandRt`) and subsequently grouped into
#' sets of merge candidates if they are (after expansion) overlapping in both
#' m/z and rt (within the **same** sample). Note that **each** peak gets
#' expanded by `expandRt` and `expandMz`, thus peaks differing by less than
#' `2 * expandMz` (or `2 * expandRt`) will be evaluated for merging.
#' Peak merging is performed along the retention time axis, i.e., the peaks are
#' first ordered by their `"rtmin"` and merge candidates are defined iteratively
#' starting with the first peak.
#' Candidate peaks are merged if the
#' average intensity of the 3 data points in the middle position between them
#' (i.e., at half the distance between `"rtmax"` of the first and `"rtmin"` of
#' the second peak) is larger than a certain proportion (`minProp`) of the
#' smaller (`"maxo"`) intensity of both peaks. In cases in which this calculated
#' mid point is not located between the apexes of the two peaks (e.g., if the
#' peaks are largely overlapping) the average signal intensity at half way
#' between the apexes is used instead. Candidate peaks are not merged if all 3
#' data points between them have `NA` intensities.
#'
#' Merged peaks get the `"mz"`, `"rt"`, `"sn"` and `"maxo"` values from the
#' peak with the largest signal (`"maxo"`) as well as its row in the metadata
#' of the peak (`chromPeakData`). The `"rtmin"` and `"rtmax"` of the merged
#' peaks are updated and `"into"` is recalculated based on all signal between
#' `"rtmin"` and `"rtmax"` and the newly defined `"mzmin"` and `"mzmax"` (which
#' is the range of `"mzmin"` and `"mzmax"` of the merged peaks after expanding
#' by `expandMz` and `ppm`). The reported `"mzmin"` and `"mzmax"` for the
#' merged peak represents the m/z range of all non-NA intensities used for the
#' calculation of the peak signal (`"into"`).
#'
#' @param BPPARAM parameter object to set up parallel processing. Uses the
#'     default parallel processing setup returned by `bpparam()`. See
#'     [bpparam()] for details and examples.
#'
#' @param chunkSize For `refineChromPeaks` if `object` is either an
#'     `XcmsExperiment`: `integer(1)` defining the number of files (samples)
#'     that should be loaded into memory and processed at the same time.
#'     Peak refinement is then performed in parallel (per sample) on this subset
#'     data. This setting thus allows to balance between memory
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
#' @param expandRt For `MergeNeighboringPeaksParam`: `numeric(1)` defining by
#'     how many seconds the retention time window is expanded on both sides to
#'     check for overlapping peaks.
#'
#' @param expandMz For `MergeNeighboringPeaksParam`: `numeric(1)` constant
#'     value by which the m/z range of each chromatographic peak is expanded
#'     (on both sides!) to check for overlapping peaks.
#'
#' @param maxPeakwidth For `CleanPeaksParam`: `numeric(1)` defining the maximal
#'     allowed peak width (in retention time).
#'
#' @param minProp For `MergeNeighboringPeaksParam`: `numeric(1)` between `0`
#'     and `1` representing the proporion of intensity required for peaks to be
#'     joined. See description for more details. With default (`minProp = 0.75`)
#'     only peaks are joined if the signal half way between them is larger than
#'     75% of the smallest of the two peak's `"maxo"` (maximal intensity at
#'     peak apex).
#'
#' @param msLevel `integer` defining for which MS level(s) the chromatographic
#'     peaks should be cleaned.
#'
#' @param nValues For `FilterIntensityParam`: `integer(1)` defining the number
#'     of data points (for each chromatographic peak) that have to be
#'     `>= threshold`. Defaults to `nValues = 1`.
#'
#' @param object [XCMSnExp] or [XcmsExperiment] object with identified
#'     chromatographic peaks.
#'
#' @param param Object defining the refinement method and its settings.
#'
#' @param ppm For `MergeNeighboringPeaksParam`: `numeric(1)` defining a m/z
#'     relative value (in parts per million) by which the m/z range of each
#'     chromatographic peak is expanded to check for overlapping peaks.
#'
#' @param threshold For `FilterIntensityParam`: `numeric(1)` defining the
#'     threshold below which peaks are removed.
#'
#' @param value For `FilterIntensityParam`: `character(1)` defining the name
#'     of the column in `chromPeaks` that contains the values to be used for
#'     the filtering.
#'
#' @param ... ignored.
#'
#' @return `XCMSnExp` or [XcmsExperiment] object with the refined
#'     chomatographic peaks.
#'
#' @author Johannes Rainer, Mar Garcia-Aloy
#'
#' @md
#'
#' @family chromatographic peak refinement methods
#'
#' @rdname refineChromPeaks
#'
#' @examples
#'
#' ## Load a test data set with detected peaks
#' data(faahko_sub)
#' ## Update the path to the files for the local system
#' dirname(faahko_sub) <- system.file("cdf/KO", package = "faahKO")
#'
#' ## Disable parallel processing for this example
#' register(SerialParam())
#'
#' ####
#' ## CleanPeaksParam:
#'
#' ## Distribution of chromatographic peak widths
#' quantile(chromPeaks(faahko_sub)[, "rtmax"] - chromPeaks(faahko_sub)[, "rtmin"])
#'
#' ## Remove all chromatographic peaks with a width larger 60 seconds
#' data <- refineChromPeaks(faahko_sub, param = CleanPeaksParam(60))
#'
#' quantile(chromPeaks(data)[, "rtmax"] - chromPeaks(data)[, "rtmin"])
#'
#' ####
#' ## FilterIntensityParam:
#'
#' ## Remove all peaks with a maximal intensity below 50000
#' res <- refineChromPeaks(faahko_sub,
#'     param = FilterIntensityParam(threshold = 50000))
#'
#' nrow(chromPeaks(faahko_sub))
#' nrow(chromPeaks(res))
#'
#' ####
#' ## MergeNeighboringPeaksParam:
#'
#' ## Subset to a single file
#' xd <- filterFile(faahko_sub, file = 1)
#'
#' ## Example of a split peak that will be merged
#' mzr <- 305.1 + c(-0.01, 0.01)
#' chr <- chromatogram(xd, mz = mzr, rt = c(2700, 3700))
#' plot(chr)
#'
#' ## Combine the peaks
#' res <- refineChromPeaks(xd, param = MergeNeighboringPeaksParam(expandRt = 4))
#' chr_res <- chromatogram(res, mz = mzr, rt = c(2700, 3700))
#' plot(chr_res)
#'
#' ## Example of a peak that was not merged, because the signal between them
#' ## is lower than the cut-off minProp
#' mzr <- 496.2 + c(-0.01, 0.01)
#' chr <- chromatogram(xd, mz = mzr, rt = c(3200, 3500))
#' plot(chr)
#' chr_res <- chromatogram(res, mz = mzr, rt = c(3200, 3500))
#' plot(chr_res)
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
