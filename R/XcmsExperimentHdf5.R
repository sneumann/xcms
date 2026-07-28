#' @include hidden_aliases.R

#' @title xcms result object for very large data sets
#'
#' @aliases XcmsExperimentHdf5-class
#'
#' @name XcmsExperimentHdf5
#'
#' @description
#'
#' The *xcms* result objects [XcmsExperiment()] and [XCMSnExp()] keep all
#' preprocessing results in memory and can thus (depending on the size of the
#' data set) require a large amount of memory. In contrast, the
#' `XcmsExperimentHdf5` class, by using an on-disk data storage mechanism,
#' has a much lower memory footprint allowing also to analyze very large data
#' sets on regular computer systems such as desktop or laptop computers. With
#' some exceptions, including additional parameters, the functionality and
#' usability of this object is identical to the default `XcmsExperiment`
#' object.
#'
#' This help page lists functions that have additional or different parameters
#' or properties than the respective methods for [XcmsExperiment()] objects.
#' For all other functions not listed here the usability is identical to those
#' for the [XcmsExperiment()] object (see the respective help page for
#' information).
#'
#' @details
#'
#' The `XcmsExperimentHdf5` object stores all preprocessing results (except
#' adjusted retention times, which are stored as an additional spectra variable
#' in the object's [Spectra::Spectra()] object), in a file in HDF5 format.
#'
#' `XcmsExperimentHdf5` uses a different naming scheme for chromatographic
#' peaks: for efficiency reasons, chromatographic peak data is organized by
#' sample and MS level. The chrom peak IDs are hence in the format
#' `"CP"<MS level>"S"<sample id><chrom peak index>` with `<MS level>` being
#' the MS level in which the chromatographic peaks were detected and
#' `<sample id>` the ID of the sample (usually related to the index in the
#' original `MsExperiment` object) and the `<chrom peak index>` the index of the
#' chromatographic peak in the chrom peak matrix **of that sample** and
#' MS level.
#'
#' HDF5 files do not support parallel processing, thus preprocessing results
#' need to be stored or loaded sequentially.
#'
#' All functionality for `XcmsExperimentHdf5` objects is optimized to reduce
#' memory demand at the cost of eventually lower performance.
#'
#' @section Conversion between `XcmsExperiment` and `XcmsExperimentHdf5`:
#'
#' To use the `XcmsExperimentHdf5` class for preprocessing results, the
#' `hdf5File` parameter of the [findChromPeaks()] function needs to be defined,
#' specifying the path and name of the HDF5 file to store the results. In
#' addition it is possible to convert a `XcmsExperiment` object to a
#' `XcmsExperimentHdf5` object with the `toXcmsExperimentHdf5()` function. All
#' present preprocessing results will be stored to the specified HDF5 file.
#' To load all preprocessing results into memory and hence change from a
#' `XcmsExperimentHdf5` to a `XcmsExperiment` object, the `toXcmsExperiment()`
#' function can be used.
#'
#' @section Using the HDF5 file-based on-disk data storage:
#'
#' Calling [findChromPeaks()] on an `MsExperiment` using the parameter
#' `hdf5File` will return an instance of the `XcmsExperimentHdf5` class and
#' hence use the on-disk data storage mode described on this page. The results
#' are stored in the file specified with parameter `hdf5File`.
#'
#' @section Subset:
#'
#' - `[`: subset the `XcmsExperimentHdf5` object to the specified samples.
#'   Parameters `keepChromPeaks` (default `TRUE`), `keepAdjustedRtime`
#'   (default `TRUE`) and `keepFeatures` (default `FALSE`) allow to configure
#'   whether present chromatographic peaks, alignment or correspondence results
#'   should be retained. This will only change information in the object (i.e.,
#'   the reference to the respective entries in the HDF5 file), but will
#'   **not** change the content of the HDF5 file. Thus, *reverting* the
#'   retention times of detected chromatographic peaks is **not** supported and
#'   `keepChromPeaks = TRUE` with `keepAdjustedRtime = FALSE` will throw an
#'   error. Note that with `keepChromPeaks = FALSE` also `keepFeatures` is set
#'   to `FALSE`.
#'
#' - `filterChromPeaks()` and `filterFeatureDefinitions()` to filter the
#'   chromatographic peak and correspondence results, respectively. See
#'   documentation below for details. Subset using unsorted or duplicated
#'   indices is not supported.
#'
#' @section Functionality related to chromatographic peaks:
#'
#' - `chromPeaks()` gains parameter `bySample = FALSE` that, if set to `TRUE`
#'   returns a `list` of `chromPeaks` matrices, one for each sample. Due to
#'   the way data is organized in `XcmsExperimentHdf5` objects this is more
#'   efficient than `bySample = FALSE`. Thus, in cases where chrom peak data
#'   is subsequently evaluated or processed by sample, it is suggested to
#'   use `bySample = TRUE`.
#'
#' - `chromPeakData()` gains a new parameter `peaks = character()` which allows
#'   to specify from which chromatographic peaks data should be returned.
#'   For these chromatographic peaks the ID (row name in `chromPeaks()`)
#'   should be provided with the `peaks` parameter. This can reduce the memory
#'   requirement for cases in which only data of some selected chromatographic
#'   peaks needs to be extracted. Also, `chromPeakData()` supports the
#'   `bySample` parameter described for `chromPeaks()` above. All other
#'   parameters present also for `chromPeakData()` of `XcmsExperiment` objects,
#'   such as `columns` are supported.
#'
#' - `filterChromPeaks()` allows to filter the chromatographic peaks specifying
#'   which should be retainend using the `keep` parameter. This can be either
#'   a `logical`, `character` or `integer` vector. Duplicated or unsorted
#'   indices are **not** supported. Eventually present feature definitions
#'   will be updated as well. The function returns the object with the
#'   filtered chromatographic peaks.
#'
#' @section Retention time alignment:
#'
#' - `adjustRtimePeakGroups()` and `adjustRtime()` with `PeakGroupsParam`:
#'   parameter `extraPeaks` of `PeakGroupsParam` is **ignored**. Anchor peaks
#'   are thus only defined using the `minFraction` and the optional `subset`
#'   parameter.
#'
#' @section Correspondence analysis results:
#'
#' - `featureChromPeaks()`: get the mapping between features and
#'   chromatographic peaks. See [featureChromPeaks()] for details.
#'
#' - `featureDefinitions()`: similarly to `featureDefinitions()` for
#'   [XcmsExperiment] objects, this method returns a `data.frame` with the
#'   characteristics for the defined LC-MS features. The function for
#'   `XcmsExperimentHdf5` does however **not** return the `"peakidx"` column
#'   with the indices of the chromatographic peaks per feature. This
#'   information can be extracted with the [featurePeakidx()] function.
#'   Note: the columns of the data frame returned by `featureDefinitions()` are
#'   in alphabetic order.
#'
#' - `featureValues()`: for parameter `value`, the option `value = "index"`
#'   (i.e. returning the index of the chromatographic peaks within the
#'   `chromPeaks()` matrix per feature) is **not** supported.
#'
#' - `filterFeatureDefinitions()`: filter the feature definitions keeping only
#'   the specified features. Parameter `features` can be used to define the
#'   features to retain. It supports a `logical`, `integer` indices or
#'   `character` with the IDs of the features (i.e., their row names in
#'   `featureDefinitions()`). The function returns the input
#'   `XcmsExperimentHdf5` with the filtered content.
#'
#' @param bySample For `chromPeaks()` and `chromPeakData()`: `logical(1)`
#'    whether the data should be returned *by sample*, i.e. as a `list` of
#'    `matrix` or `data.frame` objects, one for each sample.
#'
#' @param columns For `chromPeakData()~: optional `character` allowing to
#'    define a subset of columns that should be included in the returned
#'    data frame. By default (`columns = character()`) the full data is
#'    returned.
#'
#' @param features For `filterFeatureDefinitions()`: defining the features to
#'    keep: either a `logical` with the same length than the number of features,
#'    an `integer` with the indices or a `character` with the ID of the
#'    features to keep.
#'
#' @param hdf5File For `toXcmsExperimentHdf5()`: `character(1)` with the path
#'    and name of the (not yet existing) file where the preprocessing results
#'    should be stored to.
#'
#' @param keep For `filterChromPeaks()`: defining the chromatographic peaks to
#'    keep: either a `logical` with the same length than the number of
#'    chromatographic peaks, an `integer` with the indices or a `character`
#'    with the IDs of the chromatographic peaks to keep.
#'
#' @param method For `filterChromPeaks()`: `character(1)`; currently
#'    only `method = "keep"` is supported.
#'
#' @param msLevel For `chromPeaks()` and `chromPeakData()`: optional `integer`
#'    with the MS level(s) from which the data should be returned. By default
#'    `msLevel = integer()` results from all MS levels are returned (if
#'    present).
#'    For `refineChromPeaks()`: `integer(1)` with the MS level from which
#'    chromatographic peaks should be refined.
#'
#' @param object `XcmsExperimentHdf5` object.
#'
#' @param param *parameter* object defining and configuring the algorithm to
#'    be used.
#'
#' @param peaks For `chromPeakData()`: optional `character` with the ID of
#'    chromatographic peaks (row name in `chromPeaks()`) for which the data
#'    should be returned. By default (`peaks = character()`) the data for all
#'    chromatographic peaks is returned.
#'
#' @param return.type For `chromPeakData()`: `character(1)` specifying the type
#'    of object that should be returned. Can be either
#'    `return.type = "DataFrame"` (the default) to return a `DataFrame`, or
#'    `return.type = "data.frame"` to return the results as a `data.frame`.
#'
#' @param ... additional parameters eventually passed to downstream functions.
#'
#' @return
#'
#' See description of the individual methods for information.
#'
#' @author Johannes Rainerr, Philippine Louail
#'
#' @examples
#'
#' ## Create a MsExperiment object representing the data from an LC-MS
#' ## experiment.
#' library(MsExperiment)
#' library(faahKO)
#'
#' ## Define the raw data files
#' fls <- c(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
#'          system.file('cdf/KO/ko16.CDF', package = "faahKO"),
#'          system.file('cdf/KO/ko18.CDF', package = "faahKO"))
#'
#' ## Define a data frame with the sample characterization
#' df <- data.frame(mzML_file = basename(fls),
#'                 sample = c("ko15", "ko16", "ko18"))
#' ## Importe the data. This will initialize a `Spectra` object representing
#' ## the raw data and assign these to the individual samples.
#' mse <- readMsExperiment(spectraFiles = fls, sampleData = df)
#'
#' ## Perform chromatographic peak detection storing the data in an HDF5 file
#' ## Parameter `hdf5File` has to be provided and needs to be the path and
#' ## name of a (not yet existing) file to which results are going to be
#' ## stored. For the example below we use a temporary file.
#' xmse <- findChromPeaks(mse, param = CentWaveParam(prefilter = c(4, 100000)),
#'     hdf5File = tempfile())
#' xmse
#'
#' ## Extract selected columnds from the chromatographic peak detection
#' ## results
#' chromPeaks(xmse, columns = c("rt", "mz", "into")) |> head()
#'
#' ## Extract the results per sample
#' res <- chromPeaks(xmse, columns = c("rt", "mz", "into"), bySample = TRUE)
#'
#' ## The chromatographic peaks of the second sample:
#' res[[2]] |> head()
#'
#' ## Convert the result object to the in-memory representation:
#' xmse_mem <- toXcmsExperiment(xmse)
#' xmse_mem
#'
NULL

setClass("XcmsExperimentHdf5",
         contains = "XcmsExperiment",
         slots = c(hdf5_file = "character",
                   hdf5_mod_count = "integer",
                   sample_id = "character",
                   chrom_peaks_ms_level = "integer", # use to keep track of the
                                        # MS levels for which chrom peaks are
                                        # available
                   gap_peaks_ms_level = "integer", # gap-filled chrom peaks
                   features_ms_level = "integer"),
         prototype = prototype(
             hdf5_file = character(),
             hdf5_mod_count = 0L,
             sample_id = character(),
             chrom_peaks_ms_level = integer(),
             gap_peaks_ms_level = integer(),
             features_ms_level = integer()
         ))

setValidity("XcmsExperimentHdf5", function(object) {
    .h5_require_rhdf5()
    if (length(object@hdf5_file)) {
        if(!file.exists(object@hdf5_file))
            return(paste0("Data storage file \"", object@hdf5_file,
                          "\" does not exist!"))
        .h5_valid_file(object@hdf5_file, object@hdf5_mod_count)
    }
    if (length(object@sample_id) != nrow(sampleData(object)))
        return(paste0("Corrupt data: number of samples does not match ",
                      "length of sample IDs."))
    TRUE
})

#' @rdname hidden_aliases
setMethod("show", "XcmsExperimentHdf5", function(object) {
    validObject(object)
    getMethod("show", "MsExperiment")(object)
    cat(" xcms results:\n")
    if (length(object@hdf5_file)) {
        if (hasChromPeaks(object)) {
            cat("  - chromatographic peaks in MS level(s):",
                paste0(object@chrom_peaks_ms_level, collapse = ", "), "\n")
        }
    if (hasAdjustedRtime(object))
        cat("  - adjusted retention times\n")
    }
    if (hasFeatures(object))
        cat("  - correspondence results in MS level(s):",
            paste(object@features_ms_level, collapse = ", "), "\n")
    cat(" results storage file:\n")
    cat("   ", object@hdf5_file, "\n", sep = "")
})

################################################################################
##
##        SUBSET AND FILTER
##
################################################################################

#' @rdname hidden_aliases
setMethod("[", "XcmsExperimentHdf5", function(x, i, j, ...) {
    if (!missing(j))
        stop("subsetting by j not supported")
    .h5_subset_xcms_experiment(x, i = i, ...)
})

#' @rdname hidden_aliases
setMethod(
    "filterMsLevel", "XcmsExperimentHdf5",
    function(object, msLevel. = uniqueMsLevels(object)) {
        if (!length(msLevel.))
            return(object)
        if (hasChromPeaks(object)) {
            object@chrom_peaks_ms_level <-
                object@chrom_peaks_ms_level[
                           object@chrom_peaks_ms_level %in% msLevel.]
            object@gap_peaks_ms_level <-
                object@gap_peaks_ms_level[
                       object@gap_peaks_ms_level %in% msLevel.]
        }
        if (hasFeatures(object))
            object@features_ms_level <-
                object@features_ms_level[object@features_ms_level %in% msLevel.]
        getMethod("filterMsLevel", "MsExperiment")(object, msLevel.)
    })

#' @rdname hidden_aliases
setMethod(
    "filterIsolationWindow", "XcmsExperimentHdf5",
    function(object, mz = numeric()) {
        if (length(mz) > 1L)
            mz <- mz[1L]
        object <- filterSpectra(object, filterIsolationWindow, mz = mz)
        if (hasChromPeaks(object) && length(mz)) {
            cn <- .h5_read_data(
                object@hdf5_file, id = object@sample_id[1L], "chrom_peak_data",
                ms_level = object@chrom_peaks_ms_level[1L], read_colnames =TRUE)
            if (all(c("isolationWindowLowerMz", "isolationWindowUpperMz") %in%
                    colnames(cn[[1L]]))) {
                mc <- .h5_filter_chrom_peaks(
                    object, object@chrom_peaks_ms_level,
                    FUN = .which_isolation_window, mz = mz,
                    chrom_peak_data = TRUE)
                object@hdf5_mod_count <- mc
            }
        }
        object
    })

#' @rdname hidden_aliases
setMethod(
    "filterRt", "XcmsExperimentHdf5",
    function(object, rt, msLevel. = uniqueMsLevels(object)) {
        if (missing(rt))
            return(object)
        rt <- range(rt)
        if (hasChromPeaks(object)) {
            msl <- intersect(msLevel., object@chrom_peaks_ms_level)
            mc <- .h5_filter_chrom_peaks(
                object, msl, FUN = .which_in_range,
                range = rt, column = "rt")
            object@hdf5_mod_count <- mc
        }
        getMethod("filterRt", "MsExperiment")(
            object, rt = rt, msLevel. = msLevel.)
    })

#' @rdname hidden_aliases
setMethod(
    "filterMzRange", "XcmsExperimentHdf5",
    function(object, mz = numeric(), msLevel. = uniqueMsLevels(object)) {
        if (missing(mz) || !length(mz))
            return(object)
        mz <- range(mz)
        if (hasChromPeaks(object)) {
            msl <- intersect(msLevel., object@chrom_peaks_ms_level)
            mc <- .h5_filter_chrom_peaks(
                object, msl, FUN = .which_in_range,
                range = mz, column = "mz")
            object@hdf5_mod_count <- mc
        }
        getMethod("filterMzRange", "MsExperiment")(
            object, mz = mz, msLevel. = msLevel.)
    })

################################################################################
##
##        CHROM PEAKS FUNCTIONALITY
##
################################################################################

#' @rdname hidden_aliases
setMethod(
    "findChromPeaks",
    signature(object = "XcmsExperimentHdf5", param = "Param"),
    function(object, param, msLevel = 1L, chunkSize = 2L, add = FALSE,
             ..., BPPARAM = bpparam()) {
        .h5_valid_file(object@hdf5_file, object@hdf5_mod_count)
        if (length(msLevel) > 1)
            stop("Currently only peak detection in a single MS level is ",
                 "supported", call. = FALSE)
        if (!msLevel %in% uniqueMsLevels(object)) {
            warning("No spectra of MS level ", msLevel, " present in 'object'",
                    call. = FALSE)
            return(object)
        }
        if (hasFeatures(object)) {
            message("Remove feature definitions")
            object <- dropFeatureDefinitions(
                object, keepAdjustedRtime = hasAdjustedRtime(object))
        }
        if (hasChromPeaks(object) && !add) {
            message("Remove previously identified chromatographic peaks")
            object <- dropChromPeaks(
                object, keepAdjustedRtime = hasAdjustedRtime(object))
        }
        res <- .mse_spectrapply_chunks(
            object, FUN = .h5_find_chrom_peaks_chunk, msLevel = msLevel,
            param = param, h5_file = object@hdf5_file,
            sample_id = object@sample_id, chunkSize = chunkSize,
            BPPARAM = BPPARAM, add = add & hasChromPeaks(object, msLevel))
        ph <- XProcessHistory(param = param,
                              type. = .PROCSTEP.PEAK.DETECTION,
                              fileIndex. = seq_along(object),
                              msLevel = msLevel)
        object@hdf5_mod_count <- max(unlist(res, use.names = FALSE))
        object@chrom_peaks_ms_level <- unique(
            c(object@chrom_peaks_ms_level, msLevel))
        object <- addProcessHistory(object, ph)
        object
    })

#' @rdname hidden_aliases
setMethod(
    "dropChromPeaks", "XcmsExperimentHdf5",
    function(object, keepAdjustedRtime = FALSE) {
        if (hasChromPeaks(object)) {
            pt <- vapply(object@processHistory, processType, character(1))
            object@processHistory <- dropProcessHistoriesList(
                object@processHistory,
                type = c(.PROCSTEP.PEAK.DETECTION, .PROCSTEP.PEAK.GROUPING,
                         .PROCSTEP.PEAK.FILLING, .PROCSTEP.CALIBRATION,
                         .PROCSTEP.PEAK.REFINEMENT))
            object@chrom_peaks_ms_level <- integer()
            if (hasAdjustedRtime(object) && !keepAdjustedRtime) {
                ## remove if alignment performed AFTER chrom peaks
                nom <- length(pt) + 1L
                idx_cp <- .match_last(.PROCSTEP.PEAK.DETECTION, pt,
                                      nomatch = nom)
                idx_al <- .match_last(.PROCSTEP.RTIME.CORRECTION, pt,
                                      nomatch = nom)
                if (idx_al > idx_cp)
                    object <- dropAdjustedRtime(object)
            }
            if (hasFeatures(object)) {
                object <- dropFeatureDefinitions(
                    object, keepAdjustedRtime = keepAdjustedRtime)
            }
        }
        object
    })

#' @rdname hidden_aliases
setMethod("hasChromPeaks", "XcmsExperimentHdf5",
          function(object, msLevel = integer()) {
              .h5_require_rhdf5()
              if (!length(object)) return(FALSE)
              if (length(object@chrom_peaks_ms_level)) {
                  if (length(msLevel))
                      all(msLevel %in% object@chrom_peaks_ms_level)
                  else TRUE
              } else FALSE
          })

#' @rdname hidden_aliases
setReplaceMethod("chromPeaks", "XcmsExperimentHdf5", function(object, value) {
    stop("Not implemented for ", class(object)[1L])
})

#' @rdname hidden_aliases
setMethod(
    "chromPeaks", "XcmsExperimentHdf5",
    function(object, rt = numeric(), mz = numeric(), ppm = 0,
             msLevel = integer(),  type = c("any", "within", "apex_within"),
             isFilledColumn = FALSE, columns = character(), bySample = FALSE) {
        if (isFilledColumn)
            stop("Parameter 'isFilledColumn = TRUE' is not supported")
        type <- match.arg(type)
        .h5_require_rhdf5()
        if (!length(object))
            return(object@chromPeaks)
        .h5_check_mod_count(object@hdf5_file, object@hdf5_mod_count)
        if (!hasChromPeaks(object, msLevel = msLevel))
            return(object@chromPeaks)
        msl <- object@chrom_peaks_ms_level
        if (length(msLevel))
            msl <- msl[msl %in% msLevel]
        .h5_chrom_peaks(object, msLevel = msl, columns = columns,
                        by_sample = bySample, mz = mz, rt = rt, ppm = ppm,
                        type = type)
    })

#' @rdname hidden_aliases
setReplaceMethod(
    "chromPeakData", "XcmsExperimentHdf5",
    function(object, value) {
    stop("Not implemented for ", class(object)[1L])
})

#' @rdname XcmsExperimentHdf5
setMethod(
    "chromPeakData", "XcmsExperimentHdf5",
    function(object, msLevel = integer(), peaks = character(),
             columns = character(), return.type = c("DataFrame", "data.frame"),
             bySample = FALSE) {
        return.type <- match.arg(return.type)
        .h5_require_rhdf5()
        if (!length(object))
            return(as(object@chromPeakData, return.type))
        .h5_check_mod_count(object@hdf5_file, object@hdf5_mod_count)
        if (!length(msLevel))
            msLevel <- object@chrom_peaks_ms_level
        if (!hasChromPeaks(object, msLevel = msLevel))
            return(as(object@chromPeakData, return.type))
        if (return.type == "DataFrame")
            as(.h5_chrom_peak_data(
                object, msLevel, peaks = peaks, by_sample = FALSE,
                columns = columns), "DataFrame")
        else .h5_chrom_peak_data(
                 object, msLevel, peaks = peaks, by_sample = bySample,
                 columns = columns)
    })

#' @rdname hidden_aliases
setMethod(
    "refineChromPeaks",
    signature(object = "XcmsExperimentHdf5",
              param = "MergeNeighboringPeaksParam"),
    function(object, param, msLevel = 1L, chunkSize = 2L, BPPARAM = bpparam()) {
        .h5_require_rhdf5()
        if (!hasChromPeaks(object, msLevel = msLevel)) {
            warning("No chromatographic peaks for MS level ",
                    msLevel, " present", call. = FALSE)
            return(object)
        }
        if (hasFeatures(object)) {
            message("Removing feature definitions")
            object <- dropFeatureDefinitions(object)
        }
        validObject(param)
        .xmse_apply_chunks(
            object, FUN = .h5_xmse_merge_neighboring_peaks,
            msLevel = msLevel, expandRt = param@expandRt,
            expandMz = param@expandMz, ppm = param@ppm,
            minProp = param@minProp, BPPARAM = BPPARAM,
            chunkSize = chunkSize, SUBSET_FUN = .h5_subset_xcms_experiment,
            keepAdjustedRtime = TRUE, ignoreHistory = TRUE)
        ## Update the @hdf5_mod_count with the one from the file.
        object@hdf5_mod_count <- rhdf5::h5read(object@hdf5_file,
                                               "/header/modcount")[1L]
        xph <- XProcessHistory(param = param, date. = date(),
                               type. = .PROCSTEP.PEAK.REFINEMENT,
                               fileIndex = seq_along(object),
                               msLevel = msLevel)
        object <- addProcessHistory(object, xph)
        validObject(object)
        object
    })

#' @rdname hidden_aliases
setMethod(
    "refineChromPeaks",
    signature(object = "XcmsExperimentHdf5", param = "CleanPeaksParam"),
    function(object, param = CleanPeaksParam(), msLevel = 1L) {
        if (!hasChromPeaks(object, msLevel = msLevel)) {
            warning("No chromatographic peaks for MS level ",
                    msLevel, " present", call. = FALSE)
            return(object)
        }
        if (hasFeatures(object)) {
            message("Removing feature definitions")
            object <- dropFeatureDefinitions(object)
        }
        validObject(param)
        mc <- object@hdf5_mod_count
        for (id in object@sample_id) {
            pks <- .h5_read_data(
                object@hdf5_file, id = id, name = "chrom_peaks",
                ms_level = msLevel, read_colnames = TRUE,
                read_rownames = TRUE)[[1L]]
            keep <- which(pks[, "rtmax"] - pks[, "rtmin"] < param@maxPeakwidth)
            pkd <- .h5_read_data(
                object@hdf5_file, id = id, name = "chrom_peak_data",
                ms_level = msLevel, read_colnames = TRUE)[[1L]]
            pks <- list(pks[keep, , drop = FALSE])
            pkd <- list(pkd[keep, , drop = FALSE])
            names(pks) <- id
            names(pkd) <- id
            .h5_write_data(
                object@hdf5_file, pks, name = "chrom_peaks",
                ms_level = msLevel, replace = TRUE, write_colnames = TRUE,
                write_rownames = TRUE)
            mc <- .h5_write_data(
                object@hdf5_file, pkd, name = "chrom_peak_data",
                ms_level = msLevel, replace = TRUE, write_rownames = FALSE)
        }
        object@hdf5_mod_count <- mc
        xph <- XProcessHistory(
            param = param, date. = date(), type. = .PROCSTEP.PEAK.REFINEMENT,
            fileIndex = seq_along(object), msLevel = msLevel)
        object <- addProcessHistory(object, xph)
        object
    })

#' @rdname refineChromPeaks
setMethod(
    "refineChromPeaks",
    signature(object = "XcmsExperimentHdf5", param = "FilterIntensityParam"),
    function(object, param, msLevel = 1L, ...) {
        if (!hasChromPeaks(object, msLevel = msLevel)) {
            warning("No chromatographic peaks for MS level ",
                    msLevel, " present", call. = FALSE)
            return(object)
        }
        if (hasFeatures(object)) {
            message("Removing feature definitions")
            object <- dropFeatureDefinitions(object)
        }
        validObject(param)
        mc <- object@hdf5_mod_count
        rtc <- c("rtmin", "rtmax")
        mzc <- c("mzmin", "mzmax")
        for (i in seq_along(object@sample_id)) {
            id <- object@sample_id[i]
            pks <- .h5_read_data(
                object@hdf5_file, id = id, name = "chrom_peaks",
                ms_level = msLevel, read_colnames = TRUE,
                read_rownames = TRUE)[[1L]]
            if (param@nValues == 1) {
                keep <- which(pks[, param@value] >= param@threshold)
            } else {
                s <- filterMsLevel(spectra(object[i]), msLevel)
                if (hasAdjustedRtime(object))
                    rt <- s$rtime_adjusted
                else
                    rt <- rtime(s)
                p_data <- peaksData(s, c("mz", "intensity"),
                                    return.type = "list")
                keep <- vapply(seq_len(nrow(pks)), function(z) {
                    rt_idx <- between(rt, pks[z, rtc])
                    vals <- vapply(p_data[rt_idx], .aggregate_intensities,
                                   mzr = pks[z, mzc], numeric(1))
                    sum(vals >= param@threshold, na.rm = TRUE) >= param@nValues
                }, NA)
            }
            pkd <- .h5_read_data(
                object@hdf5_file, id = id, name = "chrom_peak_data",
                ms_level = msLevel, read_colnames = TRUE)[[1L]]
            pks <- list(pks[keep, , drop = FALSE])
            pkd <- list(pkd[keep, , drop = FALSE])
            names(pks) <- id
            names(pkd) <- id
            .h5_write_data(
                object@hdf5_file, pks, name = "chrom_peaks",
                ms_level = msLevel, replace = TRUE, write_colnames = TRUE,
                write_rownames = TRUE)
            mc <- .h5_write_data(
                object@hdf5_file, pkd, name = "chrom_peak_data",
                ms_level = msLevel, replace = TRUE, write_rownames = FALSE)
        }
        object@hdf5_mod_count <- mc
        xph <- XProcessHistory(
            param = param, date. = date(), type. = .PROCSTEP.PEAK.REFINEMENT,
            fileIndex = seq_along(object), msLevel = msLevel)
        object <- addProcessHistory(object, xph)
        object
    })

#' @rdname hidden_aliases
setMethod("hasFilledChromPeaks", "XcmsExperimentHdf5", function(object) {
    length(object@gap_peaks_ms_level) > 0
})

#' @rdname hidden_aliases
setMethod(
    "fillChromPeaks",
    signature(object = "XcmsExperimentHdf5", param = "ChromPeakAreaParam"),
    function(object, param, msLevel = 1L, chunkSize = 2L, BPPARAM = bpparam()) {
        if (length(msLevel) != 1)
            stop("Can only perform peak filling for one MS level at a time.")
        if (!hasFeatures(object, msLevel = msLevel))
            stop("No feature definitions for MS level ", msLevel, " present.")
        .h5_valid_file(object@hdf5_file, object@hdf5_mod_count)
        ## Get integration function and other info.
        ph <- .xmse_process_history(object, .PROCSTEP.PEAK.DETECTION,
                                    msLevel = msLevel)
        fill_fun <- .history2fill_fun(ph)
        mzf <- "wMean"
        if (length(ph) && inherits(ph[[1L]], "XProcessHistory")) {
            prm <- ph[[1L]]@param
            if (any(slotNames(prm) == "mzCenterFun"))
                mzf <- prm@mzCenterFun
        } else
            prm <- MatchedFilterParam()
        mzf <- paste0("mzCenter.", gsub("mzCenter.", "", mzf, fixed = TRUE))
        ## Identify for each feature the samples in which there is a missing
        ## value
        fvals <- is.na(featureValues(object, msLevel = msLevel, method = "sum"))
        fvals <- fvals[which(rowSums(fvals) > 0), , drop = FALSE]
        if (nrow(fvals)) {
            message("Defining MS area to integrate signal from")
            fr <- .h5_features_ms_region(
                object, mzmin = param@mzmin, mzmax = param@mzmax,
                rtmin = param@rtmin, rtmax = param@rtmax,
                features = rownames(fvals), ms_level = msLevel,
                minMzWidthPpm = param@minMzWidthPpm)
            ## define the features to integrate signal for each sample
            frl <- lapply(seq_len(ncol(fvals)), function(i) {
                fr[which(fvals[, i]), , drop = FALSE]
            })
            names(frl) <- seq_along(frl)
            idx <- seq_along(object)
            chunks <- split(idx, ceiling(idx / chunkSize))
            message("Integrating signal from raw data files")
            pb <- progress_bar$new(format = paste0("[:bar] :current/:",
                                                   "total (:percent) in ",
                                                   ":elapsed"),
                                   total = length(chunks), clear = FALSE)
            pb$tick(0)
            res <- lapply(chunks, function(z, ...) {
                mc <- .h5_xmse_integrate_chrom_peaks(
                    .h5_subset_xcms_experiment(
                        object, i = z, keepChromPeaks = TRUE,
                        keepAdjustedRtime = TRUE, keepFeatures = TRUE,
                        ignoreHistory = TRUE),
                    pal = frl[z], intFun = fill_fun, update_features = TRUE,
                    mzCenterFun = mzf, param = prm, msLevel = msLevel,
                    BPPARAM = BPPARAM)
                pb$tick()
                mc
            })
            object@hdf5_mod_count <- max(unlist(res, use.names = FALSE))
            object@gap_peaks_ms_level <- union(object@gap_peaks_ms_level,
                                               msLevel)
            ph <- XProcessHistory(
                param = param, date. = date(), type. = .PROCSTEP.PEAK.FILLING,
                fileIndex = seq_along(object), msLevel = msLevel)
            object <- addProcessHistory(object, ph)
        }
        object
    })

#' @rdname hidden_aliases
setMethod("dropFilledChromPeaks", "XcmsExperimentHdf5", function(object) {
    if (!.hasFilledPeaks(object))
        return(object)
    mc <- object@hdf5_mod_count
    for (msl in object@gap_peaks_ms_level) {
        for (id in object@sample_id) {
            pks <- .h5_read_data(
                object@hdf5_file, id, "chrom_peaks", msl,
                read_colnames = TRUE, read_rownames = TRUE)[[1L]]
            pkd <- .h5_read_data(
                object@hdf5_file, id, "chrom_peak_data", msl,
                read_colnames = TRUE)[[1L]]
            keep <- which(!pkd$is_filled)
            l <- list(pks[keep, , drop = FALSE])
            names(l) <- id
            .h5_write_data(
                object@hdf5_file, l, name = "chrom_peaks", ms_level = msl,
                replace = TRUE, write_colnames = TRUE, write_rownames = TRUE)
            l <- list(extractROWS(pkd, keep))
            names(l) <- id
            .h5_write_data(
                object@hdf5_file, l, name = "chrom_peak_data", ms_level = msl,
                replace = TRUE, write_rownames = FALSE)
            fmap <- .h5_read_data(
                object@hdf5_file, id, "feature_to_chrom_peaks", msl)[[1L]]
            l <- list(fmap[fmap[, 2L] %in% keep, , drop = FALSE])
            names(l) <- id
            mc <- .h5_write_data(
                object@hdf5_file, l, name = "feature_to_chrom_peaks",
                write_colnames = FALSE, write_rownames = FALSE,
                ms_level = msl)
       }
    }
    object@gap_peaks_ms_level <- integer()
    object@hdf5_mod_count <- mc
    object@processHistory <- dropProcessHistoriesList(
        object@processHistory, type = .PROCSTEP.PEAK.FILLING)
                type = c(.PROCSTEP.PEAK.DETECTION, .PROCSTEP.PEAK.GROUPING,
                         .PROCSTEP.PEAK.FILLING, .PROCSTEP.CALIBRATION,
                         .PROCSTEP.PEAK.REFINEMENT)
    object
})

#' @rdname XcmsExperimentHdf5
setMethod(
    "filterChromPeaks", "XcmsExperimentHdf5",
    function(object, keep = rep(TRUE, nrow(chromPeaks(object))),
             method = "keep", ...) {
        method <- match.arg(method)
        object <- switch(
            method,
            keep = {
                p <- unlist(.h5_chrom_peaks_rownames(object), use.names = FALSE)
                idx <- .i2index(keep, ids = p, name = "keep")
                if (anyDuplicated(idx) || is.unsorted(idx))
                    stop("Filtering with duplicated or unsorted indices",
                         " is not supported", call. = FALSE)
                mc <- .h5_filter_chrom_peaks(
                    object, msLevel = object@chrom_peaks_ms_level,
                    FUN = function(x, peak_id) {
                        which(rownames(x) %in% peak_id)
                    }, peak_id = p[idx])
                object@hdf5_mod_count <- mc
                object
            }
        )
        object
    })

#' @rdname hidden_aliases
setMethod(
    "manualChromPeaks", "XcmsExperimentHdf5",
    function(object, chromPeaks = matrix(numeric()),
             samples = seq_along(object), msLevel = 1L,
             chunkSize = 2L, BPPARAM = bpparam()) {
        if (length(msLevel) > 1L)
            stop("Can only add peaks from one MS level at a time.")
        if (is.data.frame(chromPeaks)) chromPeaks <- as.matrix(chromPeaks)
        if (!nrow(chromPeaks)) return(object)
        if (!all(c("mzmin", "mzmax", "rtmin", "rtmax") %in%
                 colnames(chromPeaks)))
            stop("'chromPeaks' lacks one or more of the required colums ",
                 "\"mzmin\", \"mzmax\", \"rtmin\" and \"rtmax\".")
        chromPeaks <- chromPeaks[, c("mzmin", "mzmax", "rtmin", "rtmax"),
                                 drop = FALSE]
        if (!all(samples %in% seq_along(object)))
            stop("'samples' out of bounds")
        if (hasFeatures(object))
            object <- dropFeatureDefinitions(object)
        pal <- lapply(samples, function(z) chromPeaks)
        names(pal) <- samples
        chunks <- split(samples, ceiling(seq_along(samples) / chunkSize))
        message("Integrating signal from raw data files")
        pb <- progress_bar$new(format = paste0("[:bar] :current/:",
                                               "total (:percent) in ",
                                               ":elapsed"),
                               total = length(chunks) + 1L, clear = FALSE)
        res <- lapply(chunks, function(z, ...) {
            pb$tick()
            .h5_xmse_integrate_chrom_peaks(
                .h5_subset_xcms_experiment(
                    object, i = z, keepChromPeaks = TRUE,
                    keepAdjustedRtime = TRUE, ignoreHistory = TRUE),
                pal = pal[z], param = CentWaveParam(), msLevel = msLevel,
                is_filled = FALSE, BPPARAM = BPPARAM)
        })
        pb$tick()
        object@hdf5_mod_count <- max(unlist(res, use.names = FALSE))
        object@chrom_peaks_ms_level <- union(
            object@chrom_peaks_ms_level, msLevel)
        object
    })

#' @rdname hidden_aliases
setMethod(
    "chromPeakSummary",
    signature(object = "XcmsExperimentHdf5", param = "BetaDistributionParam"),
    function(object, param, msLevel = 1L, chunkSize = 2L, BPPARAM = bpparam()) {
        if (length(msLevel) != 1)
            stop("Can only perform peak metrics for one MS level at a time.")
        if (!hasChromPeaks(object, msLevel = msLevel))
            stop("No ChromPeaks definitions for MS level ",
                 msLevel, " present.")
        chunks <- split(seq_along(object),
                        ceiling(seq_along(object) / chunkSize))
        pb <- progress_bar$new(format = paste0("[:bar] :current/:",
                                               "total (:percent) in ",
                                               ":elapsed"),
                               total = length(chunks) + 1L, clear = FALSE)
        pb$tick(0)
        res <- lapply(chunks, function(z) {
            pb$tick()
            x <- .h5_subset_xcms_experiment(
                object, i = z, keepChromPeaks = TRUE, keepAdjustedRtime = TRUE,
                keepFeatures = TRUE, ignoreHistory = TRUE)
            pal <- chromPeaks(x, msLevel = msLevel, bySample = TRUE,
                              columns = c("mzmin", "mzmax", "rtmin", "rtmax"))
            names(pal) <- z
            .h5_xmse_integrate_chrom_peaks(
                x, pal, intFun = .chrom_peak_beta_metrics, msLevel = msLevel,
                BPPARAM = BPPARAM, storeToHdf5 = FALSE)
        })
        res <- do.call(rbind, unlist(res, recursive = FALSE))
        pb$tick()
        res
    })

################################################################################
##
##        RETENTION TIME ALIGNMENT
##
################################################################################

#' @rdname XcmsExperimentHdf5
setMethod(
    "adjustRtimePeakGroups", c("XcmsExperimentHdf5", "PeakGroupsParam"),
    function(object, param = PeakGroupsParam(), msLevel = 1L) {
        if (!hasFeatures(object, msLevel = msLevel))
            stop("No features present. Please run 'groupChromPeaks' first.")
        subs <- param@subset
        if (!length(subs)) subs <- seq_along(object)
        if (!all(subs %in% seq_along(object)))
            stop("Parameter 'subset' is out of bounds.")
        object@sample_id <- object@sample_id[subs] # quick subset hack.
        rts <- .h5_feature_values_ms_level(msLevel, object, method = "maxint",
                                           intensity = "into", value = "rt")
        pres <- apply(rts, 1, function(z) sum(!is.na(z)))
        rts <- rts[pres >= param@minFraction * length(subs), , drop = FALSE]
        colnames(rts) <- basename(fileNames(object))[subs]
        rts[order(rowMedians(rts, na.rm = TRUE)), , drop = FALSE]
    })

#' @rdname hidden_aliases
setMethod("dropAdjustedRtime", "XcmsExperimentHdf5", function(object) {
    if (!hasAdjustedRtime(object))
        return(object)
    ptype <- vapply(object@processHistory, processType, character(1))
    nom <- length(ptype) + 1L
    idx_al <- .match_last(.PROCSTEP.RTIME.CORRECTION, ptype, nomatch = nom)
    idx_co <- .match_last(.PROCSTEP.PEAK.GROUPING, ptype, nomatch = nom)
    if (hasChromPeaks(object)) {
        fidx <- as.factor(fromFile(object))
        object <- updateChromPeaksRtime(
            object, rtraw = split(rtime(object, adjusted = TRUE), fidx),
            rtadj = split(rtime(object, adjusted = FALSE), fidx))
    }
    svs <- unique(c(spectraVariables(object@spectra), "mz", "intensity"))
    object@spectra <- selectSpectraVariables(
        object@spectra, svs[svs != "rtime_adjusted"])
    object@processHistory <- dropProcessHistoriesList(
        object@processHistory, type = .PROCSTEP.RTIME.CORRECTION, num = 1L)
    if (hasFeatures(object) && idx_co > idx_al) {
        warning("Had to remove feature definitions along with the adjusted ",
                "retention times because of the dependency between them.")
        object <- dropFeatureDefinitions(object)
    }
    object
})

################################################################################
##
##        CORRESPONDENCE
##
################################################################################

#' @rdname hidden_aliases
setMethod("hasFeatures", "XcmsExperimentHdf5", function(object,
                                                        msLevel = integer()) {
    .h5_require_rhdf5()
    if (!length(object)) return(FALSE)
    if (length(object@features_ms_level)) {
        if (length(msLevel))
            all(msLevel %in% object@features_ms_level)
        else TRUE
    } else FALSE
})

#' @rdname hidden_aliases
setMethod(
    "dropFeatureDefinitions", "XcmsExperimentHdf5",
    function(object, keepAdjustedRtime = FALSE) {
        if (!hasFeatures(object))
            return(object)
        ptype <- vapply(object@processHistory, processType, character(1))
        nom <- length(ptype) + 1L
        idx_al <- .match_last(.PROCSTEP.RTIME.CORRECTION, ptype, nomatch = nom)
        idx_co <- .match_last(.PROCSTEP.PEAK.GROUPING, ptype, nomatch = nom)
        object@processHistory <- dropProcessHistoriesList(
            object@processHistory, type = .PROCSTEP.PEAK.GROUPING, num = 1L)
        object@features_ms_level <- integer()
        if (.hasFilledPeaks(object))
            object <- dropFilledChromPeaks(object)
        if (!keepAdjustedRtime && hasAdjustedRtime(object) && idx_al > idx_co) {
            object <- dropAdjustedRtime(object)
        }
        object
    })

#' @rdname hidden_aliases
setMethod(
    "groupChromPeaks",
    signature(object = "XcmsExperimentHdf5", param = "Param"),
    function(object, param, msLevel = 1L, add = FALSE) {
        msLevel <- unique(msLevel)
        if (length(msLevel) != 1)
            stop("Can only perform the correspondence analysis on one MS",
                 " level at a time. Please repeat for other MS levels ",
                 "with parameter `add = TRUE`.")
        if (!hasChromPeaks(object, msLevel))
            stop("No chromatographic peak for MS level ", msLevel,
                 " present. Please perform first a peak detection ",
                 "using the 'findChromPeaks' method.", call. = FALSE)
        if (hasFeatures(object) && !add)
            object <- dropFeatureDefinitions(object)
        if (hasFeatures(object, msLevel) && add)
            stop("Adding new features to existing correspondence results of ",
                 "the same MS level is currently not supported.", call. = FALSE)
        cps <- .h5_chrom_peaks(object, msLevel = msLevel,
                               columns = c("mz", "rt"), by_sample = FALSE)
        res <- .xmse_group_cpeaks(cps, param = param)
        if (!nrow(res))
            return(object)
        object@features_ms_level <- as.integer(
            unique(c(object@features_ms_level, msLevel)))
        cpk_idx <- res$peakidx
        res$peakidx <- NULL
        attr(res, "row.names") <- .featureIDs(
            nrow(res), prefix = paste0("FT", msLevel), min_len = 6)
        ## Save features to "/features/ms_<msLevel>/feature_definitions"
        .h5_write_data(object@hdf5_file, list(features = res),
                       name = "feature_definitions", ms_level = msLevel,
                       replace = TRUE, write_rownames = TRUE)
        ## Define the mapping of features to chrom peaks PER SAMPLE
        i <- ncol(cps)
        cps <- unname(cps)
        sample_f <- factor(
            unlist(lapply(cpk_idx, function(idx) cps[idx, i])),
            levels = seq_along(object@sample_id))
        fts_to_cps <- split.data.frame(
            cbind(rep(seq_len(nrow(res)), lengths(cpk_idx)),
                  unlist(cpk_idx, use.names = FALSE)), f = sample_f)
        ## Fix indices relative to chrom peak matrices WITHIN sample
        ## Maybe add an index -1 for features with no assigned peaks? Just
        ## to ensure all features are listed at least once?
        npks <- as.integer(
            table(factor(cps[, i], levels = seq_along(object@sample_id))))
        fts_to_cps <- mapply(function(a, b) {
            if (nrow(a))
                a[, 2L] <- a[, 2L] - b
            a
        }, fts_to_cps, as.list(c(0, cumsum(npks)[-length(npks)])))
        names(fts_to_cps) <- object@sample_id
        mc <- .h5_write_data(
            object@hdf5_file, fts_to_cps, name = "feature_to_chrom_peaks",
            ms_level = rep(msLevel, length(fts_to_cps)), replace = TRUE,
            write_rownames = FALSE, write_colnames = FALSE)
        object@hdf5_mod_count <- mc
        xph <- XProcessHistory(param = param, type. = .PROCSTEP.PEAK.GROUPING,
                               fileIndex = seq_along(object), msLevel = msLevel)
        object <- addProcessHistory(object, xph)
        validObject(object)
        object
    })

#' @rdname hidden_aliases
setMethod(
    "featureArea", "XcmsExperimentHdf5",
    function(object, mzmin = min, mzmax = max, rtmin = min,
             rtmax = max, features = character(), msLevel = 1L,
             minMzWidthPpm = 0.0) {
        if (!hasFeatures(object, msLevel))
            stop("No correspondence results available. Please run ",
                 "'groupChromPeaks' first.", call. = FALSE)
        if (!length(features))
            features <- .h5_feature_definitions_rownames(object, msLevel)[[1L]]
        .h5_features_ms_region(
            object, mzmin = mzmin, mzmax = mzmax, rtmin = rtmin,
            rtmax = rtmax, features, ms_level = msLevel,
            minMzWidthPpm = minMzWidthPpm)
    })

#' @rdname hidden_aliases
setReplaceMethod("featureDefinitions", "XcmsExperimentHdf5",
                 function(object, value) {
                     stop("Not implemented for ", class(object)[1L])
                 })

#' @rdname hidden_aliases
setMethod(
    "featureDefinitions", "XcmsExperimentHdf5",
    function(object, mz = numeric(), rt = numeric(), ppm = 0,
             type = c("any", "within", "apex_within"), msLevel = integer()) {
        if (!hasFeatures(object))
            return(object@featureDefinitions)
        type <- match.arg(type)
        if (length(msLevel))
            msLevel <- intersect(msLevel, object@features_ms_level)
        else msLevel <- object@features_ms_level
        if (!length(msLevel))
            return(object@featureDefinitions)
        fd <- .h5_read_data(object@hdf5_file, rep("features", length(msLevel)),
                            name = "feature_definitions", ms_level = msLevel,
                            read_rownames = TRUE)
        msl <- rep(msLevel, vapply(fd, nrow, 1L))
        fd <- rbindlistWithRownames(fd, fill = TRUE)
        fd$ms_level <- msl
        .subset_feature_definitions(fd, mz = mz, rt = rt,
                                    ppm = ppm, type = type)
    })

#' @rdname featureChromPeaks
setMethod("featureChromPeaks", "XcmsExperimentHdf5",
          function(object, msLevel = integer()) {
              if (!length(msLevel))
                  msLevel <- object@features_ms_level
              if (length(m <- msLevel[!msLevel %in% object@features_ms_level]))
                  stop("No features available for MS level(s) ",
                       paste0(m, collapse = ", "), call. = FALSE)
              map <- do.call(rbind, lapply(msLevel, function(m) {
                  tmp <- .h5_read_data(object@hdf5_file, object@sample_id,
                                       "feature_to_chrom_peak",
                                       rep(m, length(object@sample_id)))
                  fid <- .h5_feature_definitions_rownames(object, m)[[1L]]
                  cid <- .h5_chrom_peaks_rownames(object, m)
                  do.call(rbind, mapply(tmp, cid, FUN = function(a, b) {
                      cbind(fid[a[, 1L]], b[a[, 2L]])
                  }))
              }))
              colnames(map) <- c("feature_id", "chrom_peak_id")
              fids <- unlist(.h5_feature_definitions_rownames(object, msLevel))
              as.data.frame(
                  map[order(match(map[, "feature_id"], fids)), , drop = FALSE])
          })

#' @rdname featureChromPeaks
setMethod("featurePeakidx", "XcmsExperimentHdf5", function(object,
                                                           msLevel = integer()){
    if (!length(msLevel))
        msLevel <- object@features_ms_level
    map <- featureChromPeaks(object, msLevel)
    pids <- unlist(.h5_chrom_peaks_rownames(object, msLevel), use.names = FALSE)
    split(match(map$chrom_peak_id, pids),
          factor(map$feature_id, levels = unique(map$feature_id)))
})

#' @rdname hidden_aliases
setMethod(
    "featureValues", "XcmsExperimentHdf5",
    function(object, method = c("medret", "maxint", "sum"), value = "into",
             intensity = "into", filled = TRUE, missing = NA_real_,
             msLevel = integer()) {
        if (!length(msLevel)) msLevel <- object@features_ms_level
        if (!hasFeatures(object, msLevel = msLevel))
            stop("No feature definitions for MS level(s) ", msLevel," present.")
        method <- match.arg(method)
        if (value == "index")
            stop("'featureValues' for 'XcmsExperimentHdf5' does not support ",
                 "'value = \"index\"'.", call. = FALSE)
        if (is.character(missing) && !(missing %in% c("rowmin_half")))
            stop("if 'missing' is not 'NA' or a numeric it should",
                 " be one of: \"rowmin_half\".")
        msLevel <- intersect(msLevel, object@features_ms_level)
        vals <- rbindFill(
            lapply(msLevel, .h5_feature_values_ms_level, x = object,
                   method = method, value = value, intensity = intensity,
                   filled = filled)
        )
        colnames(vals) <- basename(fileNames(object))
        missing <- missing[1L]
        if (!is.na(missing)) {
            if (is.numeric(missing))
                vals[is.na(vals)] <- missing
            if (missing == "rowmin_half")
                for (i in seq_len(nrow(vals))) {
                    nas <- is.na(vals[i, ])
                    if (any(nas))
                        vals[i, nas] <- min(vals[i, ], na.rm = TRUE) / 2
            }
        }
        vals
    })

#' @rdname XcmsExperimentHdf5
setMethod(
    "filterFeatureDefinitions", "XcmsExperimentHdf5",
    function(object, features = integer()) {
        if (!length(features))
            return(object)
        if (!hasFeatures(object))
            stop("No feature definitions present! Please run ",
                 "'groupChromPeaks' first.")
        fid <- rownames(featureDefinitions(object))
        idx <- .i2index(features, ids = fid, name = "features")
        if (anyDuplicated(idx) || is.unsorted(idx))
            stop("Filtering with duplicated or unsorted indices",
                 " is not supported", call. = FALSE)
        ## filter the feature definitions and all feature mappings.
        mc <- .h5_filter_feature_definitions(
            object, msLevel = object@features_ms_level, feature_id = fid[idx])
        object@hdf5_mod_count <- mc
        object
    })

#' @rdname hidden_aliases
setMethod(
    "manualFeatures", "XcmsExperimentHdf5",
    function(object, peakIdx = list(), msLevel = 1L) {
        if (!length(peakIdx)) return(object)
        if (length(msLevel) > 1L)
            stop("Can only define features for one MS level at a time")
        if (!hasChromPeaks(object))
            stop("No chromatographic peaks present. ",
                 "Please run 'findChromPeaks' first.")
        ## - Get the full chrom peaks matrix with mz and rt.
        pks <- chromPeaks(object, columns = c("rt", "mz"))
        n <- nrow(pks)
        ## - Get the IDs (rownames) for the peakIdx
        pks <- lapply(peakIdx, function(z) {
            if (!is.numeric(z))
                stop("'peakIdx' is expected to be a list of integer indices")
            if (any(z < 1L) || any(z > n))
                stop("Provided peak indices out of bounds")
            pks[sort(z), , drop = FALSE] # peaks are expected to be sorted
        })
        ## - build featureDefinitions and write to HDF5
        fdef <- as.data.frame(do.call(rbind, lapply(pks, function(z) {
            c(mzmed = median(z[, "mz"]),
              mzmin = min(z[, "mz"]),
              mzmax = max(z[, "mz"]),
              rtmed = median(z[, "rt"]),
              rtmin = min(z[, "rt"]),
              rtmax = max(z[, "rt"]),
              npeaks = nrow(z))
        })))
        fdef$npeaks <- as.integer(fdef$npeaks)
        fdef$ms_level <- msLevel
        fdef_have <- featureDefinitions(object, msLevel = msLevel)
        have_features <- hasFeatures(object, msLevel = msLevel)
        if (have_features)
            from_idx <- max(as.integer(sub(paste0("FT", msLevel), "",
                                           rownames(fdef_have)))) + 1L
        else from_idx <- 1L
        fid <- .featureIDs(nrow(fdef), prefix = paste0("FT", msLevel),
                           from = from_idx, min_len = 6)
        rownames(fdef) <- fid
        if (have_features)
            fdef <- MsCoreUtils::rbindFill(fdef_have, fdef)
        fidx <- match(fid, rownames(fdef))
        mc <- .h5_write_data(object@hdf5_file, list(features = fdef),
                             name = "feature_definitions", ms_level = msLevel,
                             replace = TRUE, write_rownames = TRUE)
        object@features_ms_level <- union(object@features_ms_level,
                                          as.integer(msLevel))
        ## - Process feature to chrom peak mapping.
        pks <- do.call(rbind, pks)
        pks <- cbind(pks, feature_idx = rep(fidx, fdef$npeaks[fidx]))
        for (id in object@sample_id) {
            cpk_ids <- rhdf5::h5read(file = object@hdf5_file,
                                     paste0("/", id, "/ms_", msLevel,
                                            "/chrom_peaks_rownames"))
            idx <- match(rownames(pks), cpk_ids)
            keep <- !is.na(idx)
            fmap <- cbind(unname(pks[keep, "feature_idx"]), idx[keep])
            if (have_features) {
                fmap <- rbind(
                    .h5_read_data(object@hdf5_file, id,
                                  "feature_to_chrom_peaks", msLevel)[[1L]],
                    fmap)
            }
            fmap <- list(fmap)
            names(fmap) <- id
            mc <- .h5_write_data(
                object@hdf5_file, fmap, "feature_to_chrom_peaks",
                msLevel, write_colnames = FALSE, write_rownames = FALSE)
        }
        object@hdf5_mod_count <- mc
        object
    })

################################################################################
##
##        OTHER FUNCTIONALIY
##
################################################################################

#' @rdname hidden_aliases
setMethod(
    "chromatogram", "XcmsExperimentHdf5",
    function(object, rt = matrix(nrow = 0, ncol = 2),
             mz = matrix(nrow = 0, ncol = 2), aggregationFun = "sum",
             msLevel = 1L, chunkSize = 2L, isolationWindowTargetMz = NULL,
             return.type = c("XChromatograms", "MChromatograms",
                             "Chromatograms"),
             include = character(),
             chromPeaks = c("apex_within", "any", "none"),
             BPPARAM = bpparam()) {
        if (!is.matrix(rt)) rt <- matrix(rt, ncol = 2L)
        if (!is.matrix(mz)) mz <- matrix(mz, ncol = 2L)
        if (length(include)) {
            warning("Parameter 'include' is deprecated, please use ",
                    "'chromPeaks' instead")
            chromPeaks <- include
        }
        if (nrow(mz) && !nrow(rt))
            rt <- cbind(rep(-Inf, nrow(mz)), rep(Inf, nrow(mz)))
        if (nrow(rt) && !nrow(mz))
            mz <- cbind(rep(-Inf, nrow(rt)), rep(Inf, nrow(rt)))
        return.type <- match.arg(return.type)
        chromPeaks <- match.arg(chromPeaks)
        if (!hasChromPeaks(object, msLevel))
            chromPeaks <- "none"
        if (hasAdjustedRtime(object))
            object <- applyAdjustedRtime(object)
        switch(return.type,
               MChromatograms = .h5_x_chromatograms(
                   object, ms_level = msLevel, chromPeaks = chromPeaks,
                   mz = mz, rt = rt, aggregationFun = aggregationFun,
                   chunkSize = chunkSize, return.type = return.type,
                   isolationWindow = isolationWindowTargetMz,
                   BPPARAM = BPPARAM),
               XChromatograms = .h5_x_chromatograms(
                   object, ms_level = msLevel, chromPeaks = chromPeaks,
                   mz = mz, rt = rt, aggregationFun = aggregationFun,
                   chunkSize = chunkSize, return.type = return.type,
                   isolationWindow = isolationWindowTargetMz,
                   BPPARAM = BPPARAM),
               Chromatograms = .mse_chromatograms_for_ranges(
                   object, rt = rt, mz = mz, aggregationFun = aggregationFun,
                   msLevel = msLevel, isolationWindow = isolationWindowTargetMz)
               )
    })

#' @rdname hidden_aliases
setMethod(
    "chromPeakSpectra", "XcmsExperimentHdf5",
    function(object, method = c("all", "closest_rt", "closest_mz",
                                "largest_tic", "largest_bpi"),
             msLevel = 2L, expandRt = 0, expandMz = 0, ppm = 0,
             skipFilled = FALSE, peaks = character(),
             chromPeakColumns = c("rt", "mz"),
             return.type = c("Spectra", "List"), ...) {
        if (hasAdjustedRtime(object))
            object <- applyAdjustedRtime(object)
        if (!is.character(peaks))
            stop("'peaks' has to be a character vector with the IDs of",
                 " the chromatographic peaks", call. = FALSE)
        method <- match.arg(method)
        return.type <- match.arg(return.type)
        if (msLevel == 1L && method %in% c("closest_mz")) {
            warning("method = \"closest_mz\" is not supported for msLevel = 1.",
                    " Changing to method = \"all\".")
            method <- "all"
        }
        ids <- object@sample_id
        h5f <- object@hdf5_file
        object <- as(object, "MsExperiment") # need only a Spectra container
        ## Need to iterate through files/samples
        res <- lapply(seq_along(object), function(i) {
            .h5_chrom_peak_spectra_sample(
                h5f, ids[i], spectra(object[i]),
                method = method, msLevel = msLevel, expandRt = expandRt,
                expandMz = expandMz, ppm = ppm, skipFilled = skipFilled,
                peaks = unique(peaks), chromPeakColumns = chromPeakColumns)
        })
        res <- concatenateSpectra(res)
        if (return.type == "Spectra") {
            if (length(peaks))
                res[as.matrix(findMatches(peaks, res$chrom_peak_id))[, 2L]]
            else res
        } else {
            if (!length(peaks))
                peaks <- unique(res$chrom_peak_id)
            as(split(res, factor(res$chrom_peak_id))[peaks], "List")
        }
    })

#' @rdname hidden_aliases
setMethod(
    "featureSpectra", "XcmsExperimentHdf5",
    function(object, msLevel = 2L, expandRt = 0, expandMz = 0, ppm = 0,
             skipFilled = FALSE, return.type = c("Spectra", "List"),
             features = character(),
             method = c("all", "closest_rt", "closest_mz",
                        "largest_tic", "largest_bpi"),
             chromPeakColumns = c("rt", "mz"),
             featureColumns = c("rtmed", "mzmed"),
             ...) {
        return.type <- match.arg(return.type)
        method <- match.arg(method)
        chromPeaksMsLevel <- 1L # consider only chrom peaks from that MS level;
        ## manually restricting to MS1 for now.
        if (!hasFeatures(object, msLevel = chromPeaksMsLevel))
            stop("No feature definitions present. Please run ",
                 "'groupChromPeaks' first.")
        if (!is.character(features))
            stop("'features' has to be a character vector with the IDs of",
                 " the LC-MS features", call. = FALSE)
        fd <- featureDefinitions(object, msLevel = chromPeaksMsLevel)
        if (length(features)) {
            fd_idx <- match(features, rownames(fd))
            if (anyNA(fd_idx))
                stop(paste0(features[is.na(fd_idx)], collapse = ", "),
                     " are not valid feature IDs", call. = FALSE)
            fd_idx <- unique(fd_idx)
        }
        else {
            fd_idx <- seq_len(nrow(fd))
            features <- rownames(fd)
        }
        h5f <- object@hdf5_file
        ids <- object@sample_id
        if (hasAdjustedRtime(object))
            object <- applyAdjustedRtime(object)
        object <- as(object, "MsExperiment") # need only a Spectra container
        ## Need to iterate through files/samples
        res <- lapply(seq_along(object), function(i) {
            id <- ids[i]
            ## feature to chrom peak map.
            fmap <- .h5_read_data(h5f, id, "feature_to_chrom_peaks",
                                  chromPeaksMsLevel)[[1L]]
            fmap <- fmap[fmap[, 1L] %in% fd_idx, , drop = FALSE]
            if (nrow(fmap)) {
                s <- .h5_chrom_peak_spectra_sample(
                    h5f, id, spectra(object[i]), method = method,
                    msLevel = msLevel, expandRt = expandRt, expandMz = expandMz,
                    ppm = ppm, skipFilled = skipFilled, peaks=unique(fmap[,2L]),
                    chromPeakColumns = chromPeakColumns)
                ## map chrom peak id back to features.
                cp_ids <- rhdf5::h5read(h5f, paste0("/", id, "/ms_",
                                                    chromPeaksMsLevel,
                                                    "/chrom_peaks_rownames"),
                                        drop = TRUE)
                fmapl <- data.frame(fid = rownames(fd)[fmap[, 1L]],
                                    pid = cp_ids[fmap[, 2L]])
                s$feature_id <- fmapl$fid[match(s$chrom_peak_id, fmapl$pid)]
                s
            } else Spectra()
        })
        res <- concatenateSpectra(res)
        ## add feature columns
        fd <- fd[match(res$feature_id, rownames(fd)),
                 featureColumns, drop = FALSE]
        colnames(fd) <- paste0("feature_", colnames(fd))
        res <- .add_spectra_data(res, fd)
        if (return.type == "List") {
            res <- List(split(res, f = factor(res$feature_id,
                                              levels = unique(features))))
            res[features]
        } else
            res[to(findMatches(features, res$feature_id))]
    })

#' @rdname hidden_aliases
setMethod(
    "featureChromatograms", "XcmsExperimentHdf5",
    function(object, expandRt = 0, expandMz = 0, aggregationFun = "max",
             features = character(),
             return.type = c("XChromatograms","MChromatograms","Chromatograms"),
             chunkSize = 2L, mzmin = min, mzmax = max, rtmin = min,
             rtmax = max, featureArea = TRUE, ...,
             progressbar = TRUE, BPPARAM = bpparam()) {
        return.type <- match.arg(return.type)
        if (hasAdjustedRtime(object))
            object <- applyAdjustedRtime(object)
        switch(return.type,
               Chromatograms = .xmse_chromatograms_for_features(
                   object, expandRt = expandRt, expandMz = expandMz,
                   aggregationFun = aggregationFun, features = features,
                   mzmin = mzmin, mzmax = mzmax, rtmin = rtmin,
                   rtmax = rtmax, featureArea = featureArea),
               MChromatograms = .xmse_mchromatograms_for_features(
                   object, expandRt = expandRt, expandMz = expandMz,
                   aggregationFun = aggregationFun, features = features,
                   chunkSize = chunkSize, mzmin = mzmin, mzmax = mzmax,
                   rtmin = rtmin, rtmax = rtmax, progressbar = progressbar,
                   BPPARAM = BPPARAM),
               XChromatograms = .h5_xchromatograms_for_features(
                   object, expandRt = expandRt, expandMz = expandMz,
                   aggregationFun = aggregationFun, features = features,
                   chunkSize = chunkSize, mzmin = mzmin, mzmax = mzmax,
                   rtmin = rtmin, rtmax = rtmax, progressbar = progressbar,
                   BPPARAM = BPPARAM)
               )
    })

#' Helper method to update the retention times of the chrom peaks matrix
#'
#' @noRd
setMethod("updateChromPeaksRtime", "XcmsExperimentHdf5",
          function(object, rtraw, rtadj) {
              res <- mapply(
                  FUN = .h5_update_rt_chrom_peaks_sample,
                  object@sample_id, rtraw, rtadj,
                  MoreArgs = list(ms_level = object@chrom_peaks_ms_level,
                                  hdf5_file = object@hdf5_file))
              object@hdf5_mod_count <- max(unlist(res, use.names = FALSE))
              object
          })

#' TODO: LLLLLL
#'
#' - Vignette describing the functionality and some notes/properties.
#' - `reconstructChromPeakSpectra()`
#'
#' - SWATH support: need to check if it's already available.
#'   - `findChromPeaksIsolationWindow()`.
#'
#' @noRd
NULL

#' Implementation notes
#'
#' - `chromPeakChromatograms()`: uses the `XcmsExperiment` method. This loads
#'   first the full chrom peaks and chromPeakData. With the current result
#'   object (i.e. the `XChromatograms`) there is not much alternative, since we
#'   need to load that data anyway.
#'
#' - `adjustRtime`: using `XcmsExperiment` implementations. The only difference
#'   is how the retention times of the chromatographic peaks get updated. For
#'   that we have now a `updateChromPeaksRtime()` method.
#' @noRd
NULL
