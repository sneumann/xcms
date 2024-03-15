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

#' @aliases adjustRtime ObiwarpParam-class PeakGroupsParam-class LamaParama-class
#'
#' @title Alignment: Retention time correction methods.
#'
#' @description
#'
#' The `adjustRtime` method(s) perform retention time correction (alignment)
#' between chromatograms of different samples/dataset. Alignment is performed by default
#' on MS level 1 data. Retention times of spectra from other MS levels, if
#' present, are subsequently adjusted based on the adjusted retention times
#' of the MS1 spectra. Note that calling `adjustRtime` on a *xcms* result object
#' will remove any eventually present previous alignment results as well as
#' any correspondence analysis results. To run a second round of alignment,
#' raw retention times need to be replaced with adjusted ones using the
#' [applyAdjustedRtime()] function.
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
#' - `PeakGroupsParam`: performs retention time correction based on the
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
#' - `LamaParama`: This function performs retention time correction by aligning
#'   chromatographic data to an external reference dataset (concept and initial
#'   implementation by Carl Brunius). The process involves identifying and
#'   aligning peaks within the experimental chromatographic data, represented
#'   as an `XcmsExperiment` object, to a predefined set of landmark features
#'   called "lamas". These landmark features are characterized by their
#'   mass-to-charge ratio (m/z) and retention time. see [LamaParama()] for more
#'   information on the method.
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
#' @param bs For `LamaParama`: `character(1)` defining the GAM moothing method.
#'     (defaults to thin plate; NB: B- and P-splines have been shown to produce
#'     artefacts).
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
#' @param lamas For `LamaParama`: `matrix` or `data.frame` with the m/z and
#'     retention times values of features (as first and second column) from the
#'     external dataset on which the alignment will be based on.
#'
#' @param localAlignment For `ObiwarpParam`: `logical(1)` whether a local
#'     alignment should be performed instead of the default global alignment.
#'
#' @param method For `LamaParama`:`character(1)` with the type of warping.
#'     Either `method = "gam"` or `method = "loess"` (default).
#'
#' @param minFraction For `PeakGroupsParam`: `numeric(1)` between 0 and 1
#'     defining the minimum required proportion of samples in which peaks for
#'     the peak group were identified. Peak groups passing this criteria will
#'     be aligned across samples and retention times of individual spectra will
#'     be adjusted based on this alignment. For `minFraction = 1` the peak
#'     group has to contain peaks in all samples of the experiment. Note that if
#'     `subset` is provided, the specified fraction is relative to the
#'     defined subset of samples and not to the total number of samples within
#'     the experiment (i.e., a peak has to be present in the specified
#'     proportion of subset samples).
#'
#' @param msLevel For `adjustRtime`: `integer(1)` defining the MS level on
#'     which the alignment should be performed.
#'
#' @param object For `adjustRtime`: an [OnDiskMSnExp()], [XCMSnExp()],
#'     [MsExperiment()] or [XcmsExperiment()] object.
#'
#' @param outlierTolerance For `LamaParama`: `numeric(1)` defining the settings
#'     for outlier removal during the fitting. By default
#'     (with `outlierTolerance = 3`), all data points with absolute residuals
#'     larger than 3 times the mean absolute residual of all data points from
#'     the first, initial fit, are removed from the final model fit.
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
#' @param ppm For `LamaParama`: `numeric(1)` defining the m/z-relative maximal
#'     allowed difference in m/z between `lamas` and chromatographic peaks. Used
#'     for the mapping of identified chromatographic peaks and lamas.
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
#' @param span For `PeakGroupsParam` and `LamaParama`: `numeric(1)` defining
#'     the degree of smoothing (if `smooth = "loess"` or `method = "loess"`).
#'     This parameter is passed to the internal call to [loess()].
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
#' @param tolerance For `LamaParama`: `numeric(1)` defining the absolute
#'     acceptable difference in m/z between lamas and chromatographic peaks.
#'     Used for the mapping of identified chromatographic peaks and `lamas`.
#'
#' @param toleranceRt For `LamaParama`: `numeric(1)` defining the absolute
#'     acceptable difference in retention time between lamas and
#'     chromatographic peaks. Used for the mapping of identified chromatographic
#'     peaks and `lamas`.
#'
#' @param value For all assignment methods: the value to set/replace.
#'
#' @param x An `ObiwarpParam`, `PeakGroupsParam` or `LamaParama` object.
#'
#' @param zeroWeight For `LamaParama`: `numeric(1)`: defines the weight of the
#'     first data point (i.e. retention times of the first lama-chromatographic
#'     peak pair). Values larger than 1 reduce warping problems in the early RT
#'     range.
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
#' `ObiwarpParam`, `PeakGroupsParam` and `LamaParama` return the respective
#' parameter object.
#'
#' `adjustRtimeGroups` returns a `matrix` with the retention times of *marker*
#' features in each sample (each row one feature, each row one sample).
#'
#' @name adjustRtime
#'
#' @family retention time correction methods
#'
#' @seealso [plotAdjustedRtime()] for visualization of alignment results.
#'
#' @author Colin Smith, Johannes Rainer, Philippine Louail, Carl Brunius
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

#' @title Extract an ion chromatogram for each chromatographic peak
#'
#' @name chromPeakChromatograms
#'
#' @description
#'
#' Extract an ion chromatogram (EIC) for each chromatographic peak in an
#' [XcmsExperiment()] object. The result is returned as an [XChromatograms()]
#' of length equal to the number of chromatographic peaks (and one column).
#'
#' @param object An [XcmsExperiment()] with identified chromatographic peaks.
#'
#' @param expandRt `numeric(1)` to eventually expand the retention time range
#'     from which the signal should be integrated. The chromatogram will
#'     contain signal from `chromPeaks[, "rtmin"] - expandRt` to
#'     `chromPeaks[, "rtmax"] + expandRt`. The default is `expandRt = 0`.
#'
#' @param expandMz `numeric(1)` to eventually expand the m/z range
#'     from which the signal should be integrated. The chromatogram will
#'     contain signal from `chromPeaks[, "mzmin"] - expandMz` to
#'     `chromPeaks[, "mzmax"] + expandMz`. The default is `expandMz = 0`.
#'
#' @param aggregationFun `character(1)` defining the function how signals
#'     within the m/z range in each spectrum (i.e. for each discrete retention
#'     time) should be aggregated. The default (`aggregationFun = "max"`)
#'     reports the largest signal for each spectrum.
#'
#' @param peaks optional `character` providing the IDs of the chromatographic
#'     peaks (i.e. the row names of the peaks in `chromPeaks(object)`) for
#'     which chromatograms should be returned.
#'
#' @param return.type `character(1)` specifying the type of the returned object.
#'     Can be either `return.type = "XChromatograms"` (the default) or
#'     `return.type = "MChromatograms"` to return either a chromatographic
#'     object with or without the identified chromatographic peaks,
#'     respectively.
#'
#' @param ... currently ignored.
#'
#' @param progressbar `logical(1)` whether the progress of the extraction
#'     process should be displayed.
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @seealso [featureChromatograms()] to extract an EIC for each feature.
#'
#' @examples
#'
#' ## Load a test data set with detected peaks
#' faahko_sub <- loadXcmsData("faahko_sub2")
#'
#' ## Get EICs for every detected chromatographic peak
#' chrs <- chromPeakChromatograms(faahko_sub)
#' chrs
#'
#' ## Order of EICs matches the order in chromPeaks
#' chromPeaks(faahko_sub) |> head()
#'
#' ## variable "sample_index" provides the index of the sample the EIC was
#' ## extracted from
#' fData(chrs)$sample_index
#'
#' ## Get the EIC for selected peaks only.
#' pks <- rownames(chromPeaks(faahko_sub))[c(6, 12)]
#' pks
#'
#' ## Expand the data on retention time dimension by 15 seconds (on each side)
#' res <- chromPeakChromatograms(faahko_sub, peaks = pks, expandRt = 5)
#' plot(res[1, ])
setGeneric("chromPeakChromatograms", function(object, ...)
    standardGeneric("chromPeakChromatograms"))

setGeneric("chromPeaks", function(object, ...) standardGeneric("chromPeaks"))
setGeneric("chromPeaks<-", function(object, value)
    standardGeneric("chromPeaks<-"))
setGeneric("chromPeakData", function(object, ...)
    standardGeneric("chromPeakData"))
setGeneric("chromPeakData<-", function(object, value)
    standardGeneric("chromPeakData<-"))

#' @title Extract spectra associated with chromatographic peaks
#'
#' @name chromPeakSpectra
#'
#' @description
#'
#' Extract (MS1 or MS2) spectra from an [XcmsExperiment] or [XCMSnExp] object
#' for identified chromatographic peaks. To return spectra for selected
#' chromatographic peaks, their *peak ID* (i.e., row name in the `chromPeaks`
#' matrix) can be provided with parameter `peaks`.
#' For `msLevel = 1L` (only supported for `return.type = "Spectra"` or
#' `return.type = "List"`) MS1 spectra within the retention time boundaries
#' (in the file in which the peak was detected) are returned. For
#' `msLevel = 2L` MS2 spectra are returned for a chromatographic
#' peak if their precursor m/z is within the retention time and m/z range of
#' the chromatographic peak. Parameter `method` allows to define whether all
#' or a single spectrum should be returned:
#'
#' - `method = "all"`: (default): return all spectra for each chromatographic
#'   peak.
#' - `method = "closest_rt"`: return the spectrum with the retention time
#'   closest to the peak's retention time (at apex).
#' - `method = "closest_mz"`: return the spectrum with the precursor m/z
#'   closest to the peaks's m/z (at apex); only supported for `msLevel > 1`.
#' - `method = "largest_tic"`: return the spectrum with the largest total
#'   signal (sum of peaks intensities).
#' - `method = "largest_bpi"`: return the spectrum with the largest peak
#'   intensity (maximal peak intensity).
#' - `method = "signal"`: only for `object` being a `XCMSnExp`: return the
#'   spectrum with the sum of intensities most similar to the peak's apex
#'   signal (`"maxo"`); only supported for `msLevel = 2L`.
#'
#' Parameter `return.type` allows to specify the *type* of the result object.
#' With `return.type = "Spectra"` (the default) a [Spectra] object with all
#' matching spectra is returned. The spectra variable `"peak_id"` of the
#' returned `Spectra` contains the ID of the chromatographic peak (i.e., the
#' rowname of the peak in the `chromPeaks` matrix) for each spectrum.
#' With `return.type = "Spectra"` a `List` of `Spectra` is returned. The
#' length of the list is equal to the number of rows of `chromPeaks`. Each
#' element of the list contains thus a `Spectra` with all spectra for one
#' chromatographic peak (or a `Spectra` of length 0 if no spectrum was found
#' for the respective chromatographic peak).
#'
#' See also the *LC-MS/MS data analysis* vignette for more details and examples.
#'
#' @param object [XcmsExperiment] or [XCMSnExp] object with identified
#'     chromatographic peaks for which spectra should be returned.
#'
#' @param msLevel `integer(1)` defining the MS level of the spectra that
#'     should be returned.
#'
#' @param expandRt `numeric(1)` to expand the retention time range of each
#'     peak by a constant value on each side.
#'
#' @param expandMz `numeric(1)` to expand the m/z range of each peak by a
#'     constant value on each side.
#'
#' @param ppm `numeric(1)` to expand the m/z range of each peak (on each side)
#'     by a value dependent on the peak's m/z.
#'
#' @param method `character(1)` specifying which spectra to include in the
#'     result. Defaults to `method = "all"`. See function description for
#'     details.
#'
#' @param peaks `character`, `logical` or `integer` allowing to specify a
#'     subset of chromatographic peaks in `chromPeaks` for which spectra should
#'     be returned (providing either their ID, a logical vector same length
#'     than `nrow(chromPeaks(x))` or their index in `chromPeaks(x)`). This
#'     parameter overrides `skipFilled`.
#'
#' @param skipFilled `logical(1)` whether spectra for filled-in peaks should
#'     be reported or not.
#'
#' @param return.type `character(1)` defining the type of result object that
#'     should be returned.
#'
#' @param BPPARAM parallel processing setup. Defaults to [bpparam()].
#'
#' @param ... ignored.
#'
#' @return
#'
#' parameter `return.type` allow to specify the type of the returned object:
#'
#' - `return.type = "Spectra"` (default): a `Spectra` object (defined in the
#'   `Spectra` package). The result contains all spectra for all peaks.
#'   Metadata column `"peak_id"` provides the ID of the respective peak
#'   (i.e. its rowname in [chromPeaks()].
#' - `return.type = "List"`: `List` of length equal to the number of
#'   chromatographic peaks is returned, each element being a `Spectra` with
#'   the spectra for one chromatographic peak.
#'
#' For backward compatibility options `"MSpectra"` and `"list"` are also
#' supported but are not suggested.
#'
#' - `return.type = "MSpectra"` (deprecated): a [MSpectra] object with elements being
#'   [Spectrum-class] objects. The result objects contains all spectra
#'   for all peaks. Metadata column `"peak_id"` provides the ID of the
#'   respective peak (i.e. its rowname in [chromPeaks()]).
#' - `return.type = "list"`: `list` of `list`s that are either of length
#'   0 or contain [Spectrum2-class] object(s) within the m/z-rt range. The
#'   length of the list matches the number of peaks.
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @examples
#'
#' ## Read a file with DDA LC-MS/MS data
#' library(MsExperiment)
#' fl <- system.file("TripleTOF-SWATH/PestMix1_DDA.mzML", package = "msdata")
#'
#' dda <- readMsExperiment(fl)
#'
#' ## Perform MS1 peak detection
#' dda <- findChromPeaks(dda, CentWaveParam(peakwidth = c(5, 15),
#'     prefilter = c(5, 1000)))
#'
#' ## Return all MS2 spectro for each chromatographic peaks as a Spectra object
#' ms2_sps <- chromPeakSpectra(dda)
#' ms2_sps
#'
#' ## spectra variable *peak_id* contain the row names of the peaks in the
#' ## chromPeak matrix and allow thus to map chromatographic peaks to the
#' ## returned MS2 spectra
#' ms2_sps$peak_id
#' chromPeaks(dda)
#'
#' ## Alternatively, return the result as a List of Spectra objects. This list
#' ## is parallel to chromPeaks hence the mapping between chromatographic peaks
#' ## and MS2 spectra is easier.
#' ms2_sps <- chromPeakSpectra(dda, return.type = "List")
#' names(ms2_sps)
#' rownames(chromPeaks(dda))
#' ms2_sps[[1L]]
#'
#' ## Parameter `msLevel` allows to define from which MS level spectra should
#' ## be returned. By default `msLevel = 2L` but with `msLevel = 1L` all
#' ## MS1 spectra with a retention time within the retention time range of
#' ## a chromatographic peak can be returned. Alternatively, selected
#' ## spectra can be returned by specifying the selection criteria/method
#' ## with the `method` parameter. Below we extract for each chromatographic
#' ## peak the MS1 spectra with a retention time closest to the
#' ## chromatographic peak's apex position. Alternatively it would also be
#' ## possible to select the spectrum with the highest total signal or
#' ## highest (maximal) intensity.
#' ms1_sps <- chromPeakSpectra(dda, msLevel = 1L, method = "closest_rt")
#' ms1_sps
#'
#' ## Parameter peaks would allow to extract spectra for specific peaks only.
#' ## Peaks can be defined with parameter `peaks` which can be either an
#' ## `integer` with the index of the peak in the `chromPeaks` matrix or a
#' ## `character` with its rowname in `chromPeaks`.
#' chromPeakSpectra(dda, msLevel = 1L, method = "closest_rt", peaks = c(3, 5))
setGeneric("chromPeakSpectra", function(object, ...)
    standardGeneric("chromPeakSpectra"))

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

#' @title Extract ion chromatograms for each feature
#'
#' @description
#'
#' Extract ion chromatograms for features in an [XcmsExperiment] or
#' [XCMSnExp-class] object. The function returns for each feature the
#' extracted ion chromatograms (along with all associated chromatographic
#' peaks) in each sample. The chromatogram is extracted from the m/z - rt
#' region that includes **all** chromatographic peaks of a feature. By default,
#' this region is defined using the range of the chromatographic peaks' m/z
#' and retention times (with `mzmin = min`, `mzmax = max`, `rtmin = min` and
#' `rtmax = max`). For some features, and depending on the data, the m/z and
#' rt range can thus be relatively large. The boundaries of the m/z - rt
#' region can also be restricted by changing parameters `mzmin`, `mzmax`,
#' `rtmin` and `rtmax` to a different functions, such as `median`.
#'
#' By default only chromatographic peaks associated with a feature are
#' included in the returned [XChromatograms] object. For `object` being an
#' `XCMSnExp` object parameter `include` allows also to return all
#' chromatographic peaks with their apex position within the selected
#' region (`include = "apex_within"`) or any chromatographic peak overlapping
#' the m/z and retention time range (`include = "any"`).
#'
#' @note
#'
#' The EIC data of a feature is extracted from every sample using the same
#' m/z - rt area. The EIC in a sample does thus not exactly represent the
#' signal of the actually identified chromatographic peak in that sample.
#' The [chromPeakChromatograms()] function would allow to extract the actual
#' EIC of the chromatographic peak in a specific sample. See also examples
#' below.
#'
#' Parameters `include`, `filled`, `n` and `value` are only supported
#' for `object` being an `XCMSnExp`.
#'
#' When extracting EICs from only the top `n` samples it can happen that one
#' or more of the features specified with `features` are dropped because they
#' have no detected peak in the *top n* samples. The chance for this to happen
#' is smaller if `x` contains also filled-in peaks (with `fillChromPeaks`).
#'
#' @param aggregationFun `character(1)` specifying the name that should be
#'     used to aggregate intensity values across the m/z value range for
#'     the same retention time. The default `"max"` returns a base peak
#'     chromatogram.
#'
#' @param BPPARAM For `object` being an `XcmsExperiment`: parallel processing
#'     setup. Defaults to `BPPARAM = bpparam()`. See [bpparam()] for more
#'     information.
#'
#' @param chunkSize For `object` being an `XcmsExperiment`: `integer(1)`
#'     defining the number of files from which the data should be loaded at
#'     a time into memory. Defaults to `chunkSize = 2L`.
#'
#' @param expandMz `numeric(1)` to expand the m/z range for each chromatographic
#'     peak by a constant value on each side. Be aware that by extending the
#'     m/z range the extracted EIC might **no longer** represent the actual
#'     identified chromatographic peak because intensities of potential
#'     additional mass peaks within each spectra would be aggregated into the
#'     final reported intensity value per spectrum (retention time).
#'
#' @param expandRt `numeric(1)` to expand the retention time range for each
#'     chromatographic peak by a constant value on each side.
#'
#' @param features `integer`, `character` or `logical` defining a subset of
#'     features for which chromatograms should be returned. Can be the index
#'     of the features in `featureDefinitions`, feature IDs (row names of
#'     `featureDefinitions`) or a logical vector.
#'
#' @param filled Only for `object` being an `XCMSnExp`: `logical(1)` whether
#'     filled-in peaks should be included in the result object. The default
#'     is `filled = FALSE`, i.e. only detected peaks are reported.
#'
#' @param include Only for `object` being an `XCMSnExp`: `character(1)`
#'     defining which chromatographic peaks (and related feature definitions)
#'     should be included in the returned [XChromatograms()].
#'     Defaults to `"feature_only"`; See description above for options and
#'     details.
#'
#' @param mzmax `function` defining how the upper boundary of the m/z region
#'     from which the EIC is integrated should be defined. Defaults to
#'     `mzmax = max` thus the largest `"mzmax"` value for all chromatographic
#'     peaks of a feature will be used.
#'
#' @param mzmin `function` defining how the lower boundary of the m/z region
#'     from which the EIC is integrated should be defined. Defaults to
#'     `mzmin = min` thus the smallest `"mzmin"` value for all chromatographic
#'     peaks of a feature will be used.
#'
#' @param n Only for `object` being an `XCMSnExp`: `integer(1)` to optionally
#'     specify the number of *top n* samples from which the EIC should be
#'     extracted.
#'
#' @param object `XcmsExperiment` or `XCMSnExp` object with grouped
#'     chromatographic peaks.
#'
#' @param progressbar `logical(1)` defining whether a progress bar is shown.
#'
#' @param return.type `character(1)` defining how the result should be
#'     returned. At present only `return.type = "XChromatograms"` is
#'     supported and the results are thus returned as an [XChromatograms()]
#'     object.
#'
#' @param rtmax `function` defining how the upper boundary of the rt region
#'     from which the EIC is integrated should be defined. Defaults to
#'     `rtmax = max` thus the largest `"rtmax"` value for all chromatographic
#'     peaks of a feature will be used.
#'
#' @param rtmin `function` defining how the lower boundary of the rt region
#'     from which the EIC is integrated should be defined. Defaults to
#'     `rtmin = min` thus the smallest `"rtmin"` value for all chromatographic
#'     peaks of a feature will be used.
#'
#' @param value Only for `object` being an `XCMSnExp`: `character(1)`
#'     specifying the column to be used to sort the samples. Can be either
#'     `"maxo"` (the default) or `"into"` to use the maximal peak intensity
#'     or the integrated peak area, respectively.
#'
#' @param ... optional arguments to be passed along to the [chromatogram()]
#'     function.
#'
#' @return [XChromatograms()] object. In future, depending on parameter
#'     `return.type`, the data might be returned as a different object.
#'
#' @name featureChromatograms
#'
#' @md
#'
#' @seealso [filterColumnsKeepTop()] to filter the extracted EICs keeping only
#'     the *top n* columns (samples) with the highest intensity.
#'     [chromPeakChromatograms()] for a function to extract an EIC for each
#'     chromatographic peak.
#'
#' @author Johannes Rainer
#'
#' @examples
#'
#' ## Load a test data set with detected peaks
#' faahko_sub <- loadXcmsData("faahko_sub2")
#'
#' ## Disable parallel processing for this example
#' register(SerialParam())
#'
#' ## Perform correspondence analysis
#' xdata <- groupChromPeaks(faahko_sub,
#'     param = PeakDensityParam(minFraction = 0.8, sampleGroups = rep(1, 3)))
#'
#' ## Get the feature definitions
#' featureDefinitions(xdata)
#'
#' ## Extract ion chromatograms for the first 3 features. Parameter
#' ## `features` can be either the feature IDs or feature indices.
#' chrs <- featureChromatograms(xdata,
#'     features = rownames(featureDefinitions)[1:3])
#'
#' ## Plot the EIC for the first feature using different colors for each file.
#' plot(chrs[1, ], col = c("red", "green", "blue"))
#'
#' ## The EICs for all 3 samples use the same m/z and retention time range,
#' ## which was defined using the `featureArea` function:
#' featureArea(xdata, features = rownames(featureDefinitions(xdata))[1:3],
#'     mzmin = min, mzmax = max, rtmin = min, rtmax = max)
#'
#' ## To extract the actual (exact) EICs for each chromatographic peak of
#' ## a feature in each sample, the `chromPeakChromatograms` function would
#' ## need to be used instead. Below we extract the EICs for all
#' ## chromatographic peaks of the first feature. We need to first get the
#' ## IDs of all chromatographic peaks assigned to the first feature:
#' peak_ids <- rownames(chromPeaks(xdata))[featureDefinitions(xdata)$peakidx[[1L]]]
#'
#' ## We can now pass these to the `chromPeakChromatograms` function with
#' ## parameter `peaks`:
#' eic_1 <- chromPeakChromatograms(xdata, peaks = peak_ids)
#'
#' ## To plot these into a single plot we need to use the
#' ## `plotChromatogramsOverlay` function:
#' plotChromatogramsOverlay(eic_1)
setGeneric("featureChromatograms", function(object, ...)
    standardGeneric("featureChromatograms"))

setGeneric("featureDefinitions", function(object, ...)
    standardGeneric("featureDefinitions"))
setGeneric("featureDefinitions<-", function(object, value)
    standardGeneric("featureDefinitions<-"))

#' @title Extract spectra associated with features
#'
#' @name featureSpectra
#'
#' @description
#'
#' This function returns spectra associated with the identified features in
#' the input object. By default, spectra are returned for all features (from
#' all MS levels), but parameter `features` allows to specify/select features
#' for which the result should be returned.
#' Parameter `msLevel` allows to define whether MS level 1 or 2 spectra
#' should be returned. For `msLevel = 1L` all MS1 spectra within the
#' retention time range of each chromatographic peak (in that respective
#' data file) associated with a feature are returned. Note that for samples
#' in which no peak was identified (or even filled-in) no spectra are
#' returned. For `msLevel = 2L` all MS2 spectra with a retention time within
#' the retention time range and their precursor m/z within the m/z range of
#' any chromatographic peak of a feature are returned.
#'
#' See also [chromPeakSpectra()] (used internally to extract spectra for
#' each chromatographic peak of a feature) for additional information,
#' specifically also on parameter `method`. By default (`method = "all"`)
#' all spectra associated with any of the chromatographic peaks of a
#' feature are returned. With any other option for `method`, a single
#' spectrum **per chromatographic peak** will be returned (hence multiple
#' spectra per feature).
#'
#' The ID of each chromatographic peak (i.e. its row name in `chromPeaks`)
#' and each feature (i.e., its row name in `featureDefinitions`) are
#' available in the returned [Spectra()] with spectra variables `"peak_id"`
#' and `"feature_id"`, respectively.
#'
#' @param object [XcmsExperiment] or [XCMSnExp] object with feature defitions.
#'
#' @inheritParams chromPeakSpectra
#'
#' @param features `character`, `logical` or `integer` allowing to specify a
#'     subset of features in `featureDefinitions` for which spectra should
#'     be returned (providing either their ID, a logical vector same length
#'     than `nrow(featureDefinitions(x))` or their index in
#'     `featureDefinitions(x)`). This parameter overrides `skipFilled` and is
#'     only supported for `return.type` being either `"Spectra"` or `"List"`.
#'
#' @param ... additional arguments to be passed along to [chromPeakSpectra()],
#'     such as `method`.
#'
#' @return
#'
#' The function returns either a [Spectra()] (for `return.type = "Spectra"`)
#' or a `List` of `Spectra` (for `return.type = "List"`). For the latter,
#' the order and the length matches parameter `features` (or if no `features`
#' is defined the order of the features in `featureDefinitions(object)`).
#'
#' Spectra variables `"peak_id"` and `"feature_id"` define to which
#' chromatographic peak or feature each individual spectrum is associated
#' with.
#'
#' @author Johannes Rainer
#'
#' @md
setGeneric("featureSpectra", function(object, ...)
    standardGeneric("featureSpectra"))

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
#'   In contrast to the  `FillChromPeaksParam` approach this method uses (all)
#'   identified chromatographic peaks of a feature to define the area
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
#' res <- loadXcmsData("faahko_sub2")
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

#' @rdname XcmsExperiment
setGeneric("filterFeatureDefinitions", function(object, ...)
           standardGeneric("filterFeatureDefinitions"))

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

#' @title Data independent acquisition (DIA): peak detection in isolation windows
#'
#' @description
#'
#' The `findChromPeaksIsolationWindow` function allows to perform a
#' chromatographic peak detection in MS level > 1 spectra of certain isolation
#' windows (e.g. SWATH pockets). The function performs a peak detection,
#' separately for all spectra belonging to the same isolation window and adds
#' them to the [chromPeaks()] matrix of the result object. Information about
#' the isolation window in which they were detected is added to
#' [chromPeakData()] data frame.
#'
#' Note that peak detection with this method does not remove previously
#' identified chromatographic peaks (e.g. on MS1 level using the
#' [findChromPeaks()] function but adds newly identified peaks to the existing
#' [chromPeaks()] matrix.
#'
#' Isolation windows can be defined with the `isolationWindow` parameter, that
#' by default uses the definition of [isolationWindowTargetMz()], i.e.
#' chromatographic peak detection is performed for all spectra with the same
#' isolation window target m/z (seprarately for each file). The parameter
#' `param` allows to define and configure the peak detection algorithm (see
#' [findChromPeaks()] for more information).
#'
#' @param object `MsExperiment`, `XcmsExperiment`, `OnDiskMSnExp` or `XCMSnExp`
#'     object with the DIA data.
#'
#' @param param Peak detection parameter object, such as a
#'     [CentWaveParam-class] object defining and configuring the chromographic
#'     peak detection algorithm.
#'     See also [findChromPeaks()] for more details.
#'
#' @param msLevel `integer(1)` specifying the MS level in which the peak
#'     detection should be performed. By default `msLevel = 2L`.
#'
#' @param isolationWindow `factor` or similar defining the isolation windows in
#'     which the peak detection should be performed with length equal to the
#'     number of spectra in `object`.
#'
#' @param chunkSize if `object` is an `MsExperiment` or `XcmsExperiment`:
#'     `integer(1)` defining the number of files (samples) that should be
#'     loaded into memory and processed at a time. See [findChromPeaks()] for
#'     more information.
#'
#' @param BPPARAM if `object` is an `MsExperiment` or `XcmsExperiment`:
#'     parallel processing setup. See [bpparam()] for more information.
#'
#' @param ... currently not used.
#'
#' @return
#'
#' An `XcmsExperiment` or `XCMSnExp` object with the chromatographic peaks
#' identified in spectra of each isolation window from each file added to the
#' `chromPeaks` matrix.
#' Isolation window definition for each identified peak are stored as additional
#' columns in [chromPeakData()].
#'
#' @author Johannes Rainer, Michael Witting
#'
#' @seealso [reconstructChromPeakSpectra()] for the function to reconstruct
#'     MS2 spectra for each MS1 chromatographic peak.
#'
#' @md
#'
#' @aliases findChromPeaksIsolationWindow
setGeneric("findChromPeaksIsolationWindow", function(object, ...)
           standardGeneric("findChromPeaksIsolationWindow"))

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
#'   axis within slices of (partially overlapping) m/z ranges. By default,
#'   these m/z ranges (bins) have a constant size. By setting `ppm` to a value
#'   larger than 0, m/z dependent bin sizes can be used instead (better
#'   representing the m/z dependent measurement error of some MS instruments).
#'   All peaks (from the same or from different samples) with their apex
#'   position being close on the retention time axis are grouped into a LC-MS
#'   feature. See in addition [do_groupChromPeaks_density()] for the core API
#'   function.
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
#'     For `PeakDensityParam`: `numeric(1)` to define m/z-dependent, increasing
#'     m/z bin sizes. If `ppm = 0` (the default) m/z bins are defined by the
#'     sequence of values from the smallest to the larges m/z value with a
#'     constant bin size of `binSize`. For `ppm` > 0 the size of each bin is
#'     increased in addition by the `ppm` of the (upper) m/z boundary of the
#'     bin. The maximal bin size (used for the largest m/z values) would then
#'     be `binSize` plus `ppm` parts-per-million of the largest m/z value of
#'     all peaks in the data set.
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

#' @title Manual peak integration and feature definition
#'
#' @description
#'
#' The `manualChromPeaks` function allows to *manually* define chromatographic
#' peaks, integrate the intensities within the specified peak area and add
#' them to the object's `chromPeaks` matrix. A peak is not added for a sample
#' if no signal was found in the respective data file.
#'
#' Because chromatographic peaks are added to eventually previously identified
#' peaks, it is suggested to run [refineChromPeaks()] with the
#' [MergeNeighboringPeaksParam()] approach to merge potentially overlapping
#' peaks.
#'
#' The `manualFeatures` function allows to manually group identified
#' chromatographic peaks into features by providing their index in the
#' object's `chromPeaks` matrix.
#'
#' @param BPPARAM parallel processing settings (see [bpparam()] for details).
#'
#' @param chromPeaks For `manualChromPeaks`: `matrix` defining the boundaries
#'     of the chromatographic peaks with one row per chromatographic peak and
#'     columns `"mzmin"`, `"mzmax"`, `"rtmin"` and `"rtmax"` defining the
#'     m/z and retention time region of each peak.
#'
#' @param chunkSize `integer(1)` defining the number of files (samples)
#'     that should be loaded into memory and processed at the same time.
#'     Peak integration is then performed in parallel (per sample) on this
#'     subset data. This setting thus allows to balance between memory
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
#' @param msLevel `integer(1)` defining the MS level in which peak integration
#'     should be performed. Only a single MS level at a time is supported.
#'     Defaults to `msLevel = 1L`.
#'
#' @param object [XcmsExperiment], [XCMSnExp] or [OnDiskMSnExp] object.
#'
#' @param peakIdx For `manualFeatures`: `list` of `integer` vectors with the
#'     indices of chromatographic peaks in the object's `chromPeaks` matrix
#'     that should be grouped into features.
#'
#' @param samples For `manualChromPeaks`: optional `integer` defining
#'     individual samples in which the peak integration should be performed.
#'     Defaults to all samples.
#'
#'
#' @param ... ignored.
#'
#' @return `XcmsExperiment` or `XCMSnExp` with the manually added
#'     chromatographic peaks or features.
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @name manualChromPeaks
#'
#' @examples
#'
#' ## Read a test dataset.
#' fls <- c(system.file("microtofq/MM14.mzML", package = "msdata"),
#'          system.file("microtofq/MM8.mzML", package = "msdata"))
#'
#' ## Define a data frame with some sample annotations
#' ann <- data.frame(
#'     injection_index = 1:2,
#'     sample_id = c("MM14", "MM8"))
#'
#' ## Import the data
#' library(MsExperiment)
#' mse <- readMsExperiment(fls)
#'
#' ## Define some arbitrary peak areas
#' pks <- cbind(
#'     mzmin = c(512, 234.3), mzmax = c(513, 235),
#'     rtmin = c(10, 33), rtmax = c(19, 50)
#' )
#' pks
#'
#' res <- manualChromPeaks(mse, pks)
#' chromPeaks(res)
#'
#' ## Peaks were only found in the second file.
setGeneric("manualChromPeaks", function(object, ...)
    standardGeneric("manualChromPeaks"))

#' @rdname manualChromPeaks
setGeneric("manualFeatures", function(object, ...)
    standardGeneric("manualFeatures"))

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

#' @title Data independent acquisition (DIA): reconstruct MS2 spectra
#'
#' @description
#'
#' *Reconstructs* MS2 spectra for each MS1 chromatographic peak (if possible)
#' for data independent acquisition (DIA) data (such as SWATH). See the
#' *LC-MS/MS analysis* vignette for more details and examples.
#'
#' @details
#'
#' In detail, the function performs for each MS1 chromatographic peak:
#'
#' - Identify all MS2 chromatographic peaks from the isolation window
#'   containing the m/z of the ion (i.e. the MS1 chromatographic peak) with
#'   approximately the same retention time than the MS1 peak (accepted rt shift
#'   can be specified with the `diffRt` parameter).
#' - Correlate the peak shapes of the candidate MS2 chromatographic peaks with
#'   the peak shape of the MS1 peak retaining only MS2 chromatographic peaks
#'   for which the correlation is `> minCor`.
#' - Reconstruct the MS2 spectrum using the m/z of all above selected MS2
#'   chromatographic peaks and their intensity (either `"maxo"` or `"into"`).
#'   Each MS2 chromatographic peak selected for an MS1 peak will thus represent
#'   one **mass peak** in the reconstructed spectrum.
#'
#' The resulting [Spectra()] object provides also the peak IDs of the MS2
#' chromatographic peaks for each spectrum as well as their correlation value
#' with spectra variables *ms2_peak_id* and *ms2_peak_cor*.
#'
#' @param object `XCMSnExp` with identified chromatographic peaks.
#'
#' @param expandRt `numeric(1)` allowing to expand the retention time range
#'     for extracted ion chromatograms by a constant value (for the peak
#'     shape correlation). Defaults to `expandRt = 0` hence correlates only
#'     the signal included in the identified chromatographic peaks.
#'
#' @param diffRt `numeric(1)` defining the maximal allowed difference between
#'     the retention time of the chromatographic peak (apex) and the retention
#'     times of MS2 chromatographic peaks (apex) to consider them as
#'     representing candidate fragments of the original ion.
#'
#' @param minCor `numeric(1)` defining the minimal required correlation
#'     coefficient for MS2 chromatographic peaks to be considered for MS2
#'     spectrum reconstruction.
#'
#' @param intensity `character(1)` defining the column in the `chromPeaks`
#'     matrix that should be used for the intensities of the reconstructed
#'     spectra's peaks. The same value from the MS1 chromatographic peaks will
#'     be used as `precursorIntensity` of the resulting spectra.
#'
#' @param peakId optional `character` vector with peak IDs (i.e. rownames of
#'     `chromPeaks`) of MS1 peaks for which MS2 spectra should be reconstructed.
#'     By default they are reconstructed for all MS1 chromatographic peaks.
#'
#' @param BPPARAM parallel processing setup. See [bpparam()] for more
#'     information.
#'
#' @param return.type `character(1)` defining the type of the returned object.
#'     Only `return.type = "Spectra"` is supported, `return.type = "MSpectra"`
#'     is deprecated.
#'
#' @param ... ignored.
#'
#' @return
#'
#' - [Spectra()] object (defined in the `Spectra` package) with the
#'   reconstructed MS2 spectra for all MS1 peaks in `object`. Contains
#'   empty spectra (i.e. without m/z and intensity values) for MS1 peaks for
#'   which reconstruction was not possible (either no MS2 signal was recorded
#'   or the correlation of the MS2 chromatographic peaks with the MS1
#'   chromatographic peak was below threshold `minCor`. Spectra variables
#'   `"ms2_peak_id"` and `"ms2_peak_cor"` (of type [CharacterList()]
#'   and [NumericList()] with length equal to the number of peaks per
#'   reconstructed MS2 spectrum) providing the IDs and the correlation of the
#'   MS2 chromatographic peaks from which the MS2 spectrum was reconstructed.
#'   As retention time the median retention times of all MS2 chromatographic
#'   peaks used for the spectrum reconstruction is reported. The MS1
#'   chromatographic peak intensity is reported as the reconstructed
#'   spectrum's `precursorIntensity` value (see parameter `intensity` above).
#'
#' @author Johannes Rainer, Michael Witting
#'
#' @md
#'
#' @seealso [findChromPeaksIsolationWindow()] for the function to perform MS2
#'     peak detection in DIA isolation windows and for examples.
#'
#' @aliases reconstructChromPeakSpectra
setGeneric("reconstructChromPeakSpectra", function(object, ...)
    standardGeneric("reconstructChromPeakSpectra"))

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
#'     chromatographic peak is expanded (on each side) to check for overlapping
#'     peaks.
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
#' faahko_sub <- loadXcmsData("faahko_sub2")
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

#' @title Save xcms result objects in a specified format
#'
#' @description
#'
#' The `storeResults` function saves an `object` resulting from processing with
#' the `xcms` package (mainly `XcmsExperiment`). Multiple formats for storing
#' and exporting are available and can be defined by the `param` argument.
#'
#' Supported `param` objects are:
#'
#' - [`RDataParam`]: Save in an .RData format file. The name of the file can be
#'  specified in the `fileName` argument.
#'
#' - [`PlainTextParam`]: Store `MsExperiment` and `XcmsExperiment` objects as a
#' folder of plain text files, folder path defined in the `path` argument.
#'
#' - `MzTabMParam`: Save in MzTab format (to be defined).
#'
#' For specific examples, see the help pages of the individual parameter classes
#' listed above.
#'
#' @param object `MsExperiment` or `XcmsExperiment` The data object that needs
#' to be saved.
#'
#' @param param The parameter object selecting and configuring the format for
#' saving. It can be one of the following classes: [`RDataParam`],
#' [`PlainTextParam`], or `MzTabMParam`.
#'
#' @param ... Optional parameters.
#'
#' @name storeResults
#'
#' @author Philippine Louail
#'
#' @examples
#'
#' ## Load a test data set with detected peaks
#' faahko_sub <- loadXcmsData("faahko_sub2")
#'
#' ## Set up parameter to save as .RData file
#' param <- RDataParam(fileName = "example_xcms_results")
#'
#' ## save as .RData
#' storeResults(object = faahko_sub, param = param)
#'
#' ## Set up parameter to save as a collection of plain text file
#' param <- PlainTextParam(path = "test/path/")
#'
#' ## Save as a collection of plain text files
#' storeResults(object = faahko_sub, param = param)
#'
#' @md
setGeneric("storeResults", function(object, param, ...) standardGeneric("storeResults"))
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
