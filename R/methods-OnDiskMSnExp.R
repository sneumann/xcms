## Methods for MSnbase's OnDiskMSnExp and MSnExp objects.
#' @include DataClasses.R functions-OnDiskMSnExp.R do_findChromPeaks-functions.R


## Main roxygen documentation for the centWace feature detection is in
## DataClasses, before the definition of the CentWaveParam class.

## The centWave peak detection method for OnDiskMSnExp:
#' @title Chromatographic peak detection using the centWave method
#'
#' @description The \code{detectChromPeaks,OnDiskMSnExp,CentWaveParam} method
#'     performs chromatographic peak detection using the \emph{centWave}
#'     algorithm on all samples from an \code{\link{OnDiskMSnExp}}
#'     object. \code{\link{OnDiskMSnExp}} objects encapsule all
#'     experiment specific data and load the spectra data (mz and intensity
#'     values) on the fly from the original files applying also all eventual
#'     data manipulations.
#'
#' @details Parallel processing (one process per sample) is supported and can
#'     be configured either by the \code{BPPARAM} parameter or by globally
#'     defining the parallel processing mode using the
#'     \code{\link{register}} method from the \code{BiocParallel}
#'     package.
#'
#' @param object For \code{findChromPeaks}: an
#'     \code{\link{OnDiskMSnExp}}  object containing the MS- and all
#'     other experiment-relevant data.
#'
#'     For all other methods: a parameter object.
#'
#' @param param An \code{CentWaveParam} object containing all settings for the
#'     centWave algorithm.
#'
#' @param BPPARAM A parameter class specifying if and how parallel processing
#'     should be performed. It defaults to \code{\link{bpparam}}.
#'     See documentation of the \code{BiocParallel} for more details. If
#'     parallel processing is enabled, peak detection is performed in parallel
#'     on several of the input samples.
#'
#' @param return.type Character specifying what type of object the method should
#'     return. Can be either \code{"XCMSnExp"} (default), \code{"list"} or
#'     \code{"xcmsSet"}.
#'
#' @param msLevel \code{integer(1)} defining the MS level on which the peak
#'     detection should be performed. Defaults to \code{msLevel = 1}.
#' 
#' @return For \code{findChromPeaks}: if \code{return.type = "XCMSnExp"} an
#'     \code{\link{XCMSnExp}} object with the results of the peak detection.
#'     If \code{return.type = "list"} a list of length equal to the number of
#'     samples with matrices specifying the identified peaks.
#'     If \code{return.type = "xcmsSet"} an \code{\linkS4class{xcmsSet}} object
#'     with the results of the peak detection.
#'
#' @seealso \code{\link{XCMSnExp}} for the object containing the results of
#'     the peak detection.
#'
#' @rdname findChromPeaks-centWave
setMethod("findChromPeaks",
          signature(object = "OnDiskMSnExp", param = "CentWaveParam"),
          function(object, param, BPPARAM = bpparam(), return.type = "XCMSnExp",
                   msLevel = 1L) {
              return.type <- match.arg(return.type, c("XCMSnExp", "list",
                                                      "xcmsSet"))
              startDate <- date()
              ## Restrict to MS X data.
              if (length(msLevel) > 1)
                  stop("Currently only peak detection in a single MS level is ",
                       "supported")
              ## Restrict to MS level for peak detection, but keep the orignal
              ## object.
              object_mslevel <- filterMsLevel(object, msLevel. = msLevel)
              if (length(object_mslevel) == 0)
                  stop("No MS level ", msLevel, " spectra present to perform ",
                       "peak detection")
              ## Check if the data is centroided
              suppressWarnings(
                  centroided <- isCentroided(object_mslevel[[1]])
              )
              ## issue #181: if there are too few mass peaks the function
              ## returns NA.
              if (is.na(centroided)) {
                  ## check all spectra in the file - takes longer.
                  centroided <- isCentroided(object_mslevel)
                  if (length(which(centroided)) > 0 &
                      length(which(!centroided)) == 0)
                      centroided <- TRUE
              }
              if (!centroided)
                  warning("Your data appears to be not centroided! CentWave",
                          " works best on data in centroid mode.")
              ## (1) split the object per file. Ensure we keep adjusted
              ##     retention times (issue #213).
              args <- list(X = 1:length(fileNames(object_mslevel)),
                           FUN = filterFile, object = object_mslevel)
              if (hasAdjustedRtime(object_mslevel))
                  args$keepAdjustedRtime <- TRUE
              ## (2) use bplapply to do the peak detection.
              resList <- bplapply(do.call("lapply", args),
                                  FUN = findChromPeaks_OnDiskMSnExp,
                                  method = "centWave",
                                  param = param, BPPARAM = BPPARAM)
              ## (3) collect the results.
              res <- .processResultList(resList,
                                        getProcHist = return.type == "xcmsSet",
                                        fnames = fileNames(object_mslevel))

              if (return.type == "list")
                  return(res$peaks)
              object <- .peaks_to_result(res, object, startDate, param, msLevel,
                                         object_mslevel)
              if (return.type == "xcmsSet")
                  as(object, "xcmsSet")
              else object
          })

## The matchedFilter peak detection method for OnDiskMSnExp:
#' @title Peak detection in the chromatographic time domain
#'
#' @description The \code{findChromPeaks,OnDiskMSnExp,MatchedFilterParam}
#'     method performs peak detection using the \emph{matchedFilter} algorithm
#'     on all samples from an \code{\link{OnDiskMSnExp}} object.
#'     \code{\link{OnDiskMSnExp}} objects encapsule all experiment
#'     specific data and load the spectra data (mz and intensity values) on the
#'     fly from the original files applying also all eventual data
#'     manipulations.
#'
#' @details Parallel processing (one process per sample) is supported and can
#'     be configured either by the \code{BPPARAM} parameter or by globally
#'     defining the parallel processing mode using the
#'     \code{\link{register}} method from the \code{BiocParallel}
#'     package.
#' 
#' @param object For \code{findChromPeaks}: an
#'     \code{\link{OnDiskMSnExp}} object containing the MS- and all
#'     other experiment-relevant data.
#'
#'     For all other methods: a parameter object.
#'
#' @param param An \code{MatchedFilterParam} object containing all settings for
#'     the matchedFilter algorithm.
#'
#' @inheritParams findChromPeaks-centWave
#'
#' @return For \code{findChromPeaks}: if \code{return.type = "XCMSnExp"} an
#'     \code{\link{XCMSnExp}} object with the results of the peak detection.
#'     If \code{return.type = "list"} a list of length equal to the number of
#'     samples with matrices specifying the identified peaks.
#'     If \code{return.type = "xcmsSet"} an \code{\linkS4class{xcmsSet}} object
#'     with the results of the peak detection.
#' 
#' @seealso \code{\link{XCMSnExp}} for the object containing the results of
#'     the chromatographic peak detection.
#'
#' @rdname findChromPeaks-matchedFilter
setMethod("findChromPeaks",
          signature(object = "OnDiskMSnExp", param = "MatchedFilterParam"),
          function(object, param, BPPARAM = bpparam(), return.type = "XCMSnExp",
                   msLevel = 1L) {
              return.type <- match.arg(return.type, c("XCMSnExp", "list",
                                                      "xcmsSet"))
              startDate <- date()
              ## Restrict to MS X data.
              if (length(msLevel) > 1)
                  stop("Currently only peak detection in a single MS level is ",
                       "supported")
              object_mslevel <- filterMsLevel(object, msLevel. = msLevel)
              if (length(object_mslevel) == 0)
                  stop("No MS level ", msLevel, " spectra present to perform ",
                       "peak detection")
              ## (1) split the object per file. Ensure we keep adjusted
              ##     retention times (issue #213).
              args <- list(X = 1:length(fileNames(object_mslevel)),
                           FUN = filterFile, object = object_mslevel)
              if (hasAdjustedRtime(object_mslevel))
                  args$keepAdjustedRtime <- TRUE
              ## (2) use bplapply to do the peak detection.
              resList <- bplapply(do.call("lapply", args),
                                  FUN = findChromPeaks_OnDiskMSnExp,
                                  method = "matchedFilter",
                                  param = param,
                                  BPPARAM = BPPARAM)
              ## (3) collect the results.
              res <- .processResultList(resList,
                                        getProcHist = return.type == "xcmsSet",
                                        fnames = fileNames(object_mslevel))
              if (return.type == "list")
                  return(res$peaks)
              object <- .peaks_to_result(res, object, startDate, param, msLevel,
                                         object_mslevel)
              if (return.type == "xcmsSet")
                  as(object, "xcmsSet")
              else object
          })

#' Simple helper function to convert the peak finding results to an XCMSnExp
#' result object.
#'
#' @noRd
.peaks_to_result <- function(res, object, startDate, param, msLevel,
                             object_mslevel) {
    xph <- XProcessHistory(param = param, date. = startDate,
                           type. = .PROCSTEP.PEAK.DETECTION,
                           fileIndex = 1:length(fileNames(object_mslevel)),
                           msLevel = msLevel)
    object <- as(object, "XCMSnExp")
    object@.processHistory <- c(processHistory(object), list(xph))
    ## if (hasAdjustedRtime(object) | hasFeatures(object))
    ##     object@msFeatureData <- new("MsFeatureData")
    pks <- do.call(rbind, res$peaks)
    if (length(pks) > 0)
        chromPeaks(object) <- cbind(pks, is_filled = 0)
    if (validObject(object))
        object
}

## massifquant
## The massifquant peak detection method for OnDiskMSnExp:
#' @title Chromatographic peak detection using the massifquant method
#'
#' @description The \code{findChromPeaks,OnDiskMSnExp,MassifquantParam}
#'     method performs chromatographic peak detection using the
#'     \emph{massifquant} algorithm on all samples from an
#'     \code{\link{OnDiskMSnExp}} object.
#'     \code{\link{OnDiskMSnExp}} objects encapsule all experiment
#'     specific data and load the spectra data (mz and intensity values) on the
#'     fly from the original files applying also all eventual data
#'     manipulations.
#'
#' @details Parallel processing (one process per sample) is supported and can
#'     be configured either by the \code{BPPARAM} parameter or by globally
#'     defining the parallel processing mode using the
#'     \code{\link{register}} method from the \code{BiocParallel}
#'     package.
#'
#' @param object For \code{findChromPeaks}: an
#'     \code{\link{OnDiskMSnExp}} object containing the MS- and all
#'     other experiment-relevant data.
#'
#'     For all other methods: a parameter object.
#'
#' @param param An \code{MassifquantParam} object containing all settings for
#'     the massifquant algorithm.
#'
#' @inheritParams findChromPeaks-centWave
#'
#' @return For \code{findChromPeaks}: if \code{return.type = "XCMSnExp"} an
#'     \code{\link{XCMSnExp}} object with the results of the peak detection.
#'     If \code{return.type = "list"} a list of length equal to the number of
#'     samples with matrices specifying the identified peaks.
#'     If \code{return.type = "xcmsSet"} an \code{\linkS4class{xcmsSet}} object
#'     with the results of the peak detection.
#'
#' @seealso \code{\link{XCMSnExp}} for the object containing the results of
#'     the peak detection.
#'
#' @rdname findChromPeaks-massifquant
setMethod("findChromPeaks",
          signature(object = "OnDiskMSnExp", param = "MassifquantParam"),
          function(object, param, BPPARAM = bpparam(), return.type = "XCMSnExp",
                   msLevel = 1L) {
              return.type <- match.arg(return.type, c("XCMSnExp", "list",
                                                      "xcmsSet"))
              startDate <- date()
              ## Restrict to MS X data.
              if (length(msLevel) > 1)
                  stop("Currently only peak detection in a single MS level is ",
                       "supported")
              object_mslevel <- filterMsLevel(object, msLevel. = msLevel)
              if (length(object_mslevel) == 0)
                  stop("No MS level ", msLevel, " spectra present to perform ",
                       "peak detection")
              ## (1) split the object per file. Ensure we keep adjusted
              ##     retention times (issue #213).
              args <- list(X = 1:length(fileNames(object_mslevel)),
                           FUN = filterFile, object = object_mslevel)
              if (hasAdjustedRtime(object_mslevel))
                  args$keepAdjustedRtime <- TRUE
              ## (2) use bplapply to do the peaks detection.
              resList <- bplapply(do.call("lapply", args),
                                  FUN = findChromPeaks_OnDiskMSnExp,
                                  method = "massifquant", param = param,
                                  BPPARAM = BPPARAM)
              ## (3) collect the results.
              res <- .processResultList(resList,
                                        getProcHist = return.type == "xcmsSet",
                                        fnames = fileNames(object_mslevel))
              if (return.type == "list")
                  return(res$peaks)
              object <- .peaks_to_result(res, object, startDate, param, msLevel,
                                         object_mslevel)
              if (return.type == "xcmsSet")
                  as(object, "xcmsSet")
              else object
          })


## MSW
## The MSW peak detection method for OnDiskMSnExp:
#' @title Single-spectrum non-chromatography MS data peak detection
#'
#' @description The \code{findChromPeaks,OnDiskMSnExp,MSWParam}
#'     method performs peak detection in single-spectrum non-chromatography MS
#'     data using functionality from the \code{MassSpecWavelet} package on all
#'     samples from an \code{\link{OnDiskMSnExp}} object.
#'     \code{\link{OnDiskMSnExp}} objects encapsule all experiment
#'     specific data and load the spectra data (mz and intensity values) on the
#'     fly from the original files applying also all eventual data
#'     manipulations.
#'
#' @details Parallel processing (one process per sample) is supported and can
#'     be configured either by the \code{BPPARAM} parameter or by globally
#'     defining the parallel processing mode using the
#'     \code{\link{register}} method from the \code{BiocParallel}
#'     package.
#'
#' @param object For \code{findChromPeaks}: an
#'     \code{\link{OnDiskMSnExp}} object containing the MS- and all
#'     other experiment-relevant data.
#'
#'     For all other methods: a parameter object.
#'
#' @param param An \code{MSWParam} object containing all settings for
#'     the algorithm.
#'
#' @inheritParams findChromPeaks-centWave
#'
#' @return For \code{findChromPeaks}: if \code{return.type = "XCMSnExp"} an
#'     \code{\link{XCMSnExp}} object with the results of the peak detection.
#'     If \code{return.type = "list"} a list of length equal to the number of
#'     samples with matrices specifying the identified peaks.
#'     If \code{return.type = "xcmsSet"} an \code{\linkS4class{xcmsSet}} object
#'     with the results of the detection.
#'
#' @seealso \code{\link{XCMSnExp}} for the object containing the results of
#'     the peak detection.
#'
#' @rdname findPeaks-MSW
setMethod("findChromPeaks",
          signature(object = "OnDiskMSnExp", param = "MSWParam"),
          function(object, param, BPPARAM = bpparam(), return.type = "XCMSnExp",
                   msLevel = 1L) {
              return.type <- match.arg(return.type, c("XCMSnExp", "list",
                                                      "xcmsSet"))
              startDate <- date()
              ## Restrict to MS X data.
              if (length(msLevel) > 1)
                  stop("Currently only peak detection in a single MS level is ",
                       "supported")
              object_mslevel <- filterMsLevel(object, msLevel. = msLevel)
              if (length(object_mslevel) == 0)
                  stop("No MS level ", msLevel, " spectra present to perform ",
                       "peak detection")

              rts <- split(rtime(object_mslevel), f = fromFile(object_mslevel))
              if (any(lengths(rts) > 1))
                  stop("The MSW method can only be applied to single spectrum,",
                       " non-chromatographic, files (i.e. with a single ",
                       "retention time).")
              ## (1) split the object per file. Ensure we keep adjusted
              ##     retention times (issue #213).
              args <- list(X = 1:length(fileNames(object_mslevel)),
                           FUN = filterFile, object = object_mslevel)
              if (hasAdjustedRtime(object_mslevel))
                  args$keepAdjustedRtime <- TRUE
              ## (2) use bplapply to do the peak detection.
              resList <- bplapply(do.call("lapply", args),
                                  FUN = findPeaks_MSW_OnDiskMSnExp,
                                  method = "MSW", param = param,
                                  BPPARAM = BPPARAM)
              ## (3) collect the results.
              res <- .processResultList(resList,
                                        getProcHist = return.type == "xcmsSet",
                                        fnames = fileNames(object_mslevel))
              if (return.type == "list")
                  return(res$peaks)
              object <- .peaks_to_result(res, object, startDate, param, msLevel,
                                         object_mslevel)
              if (return.type == "xcmsSet")
                  as(object, "xcmsSet")
              else object
          })

## The centWave with predicted isotope peak detection method for OnDiskMSnExp:
#' @title Two-step centWave peak detection considering also isotopes
#'
#' @description The \code{findChromPeaks,OnDiskMSnExp,CentWavePredIsoParam}
#'     method performs a two-step centWave-based chromatographic peak detection
#'     on all samples from an \code{\link{OnDiskMSnExp}} object.
#'     \code{\link{OnDiskMSnExp}} objects encapsule all experiment
#'     specific data and load the spectra data (mz and intensity values) on the
#'     fly from the original files applying also all eventual data
#'     manipulations.
#'
#' @details Parallel processing (one process per sample) is supported and can
#'     be configured either by the \code{BPPARAM} parameter or by globally
#'     defining the parallel processing mode using the
#'     \code{\link{register}} method from the \code{BiocParallel}
#'     package.
#'
#' @param param An \code{CentWavePredIsoParam} object with the settings for the
#'     chromatographic peak detection algorithm.
#' 
#' @inheritParams findChromPeaks-centWave
#' 
#' @return For \code{findChromPeaks}: if \code{return.type = "XCMSnExp"} an
#'     \code{\link{XCMSnExp}} object with the results of the peak detection.
#'     If \code{return.type = "list"} a list of length equal to the number of
#'     samples with matrices specifying the identified peaks.
#'     If \code{return.type = "xcmsSet"} an \code{\linkS4class{xcmsSet}} object
#'     with the results of the peak detection.
#'
#' @seealso \code{\link{XCMSnExp}} for the object containing the results of
#'     the peak detection.
#'
#' @rdname findChromPeaks-centWaveWithPredIsoROIs
setMethod("findChromPeaks",
          signature(object = "OnDiskMSnExp", param = "CentWavePredIsoParam"),
          function(object, param, BPPARAM = bpparam(), return.type = "XCMSnExp",
                   msLevel = 1L) {
              return.type <- match.arg(return.type, c("XCMSnExp", "list",
                                                      "xcmsSet"))
              startDate <- date()
              ## Restrict to MS X data.
              if (length(msLevel) > 1)
                  stop("Currently only peak detection in a single MS level is ",
                       "supported")
              object_mslevel <- filterMsLevel(object, msLevel. = msLevel)
              if (length(object_mslevel) == 0)
                  stop("No MS level ", msLevel, " spectra present to perform ",
                       "peak detection")
              ## Check if the data is centroided
              suppressWarnings(
                  centroided <- isCentroided(object_mslevel[[1]])
              )
              ## issue #181: if there are too few mass peaks the function
              ## returns NA.
              if (is.na(centroided)) {
                  ## check all spectra in the file - takes longer.
                  centroided <- isCentroided(object_mslevel)
                  if (length(which(centroided)) > 0 &
                      length(which(!centroided)) == 0)
                      centroided <- TRUE
              }
              if (!centroided)
                  warning("Your data appears to be not centroided! CentWave",
                          " works best on data in centroid mode.")
              ## (1) split the object per file. Ensure we keep adjusted
              ##     retention times (issue #213).
              args <- list(X = 1:length(fileNames(object_mslevel)),
                           FUN = filterFile, object = object_mslevel)
              if (hasAdjustedRtime(object_mslevel))
                  args$keepAdjustedRtime <- TRUE
              ## (2) use bplapply to do the peak detection.
              resList <- bplapply(do.call("lapply", args),
                                  FUN = findChromPeaks_OnDiskMSnExp,
                                  method = "centWaveWithPredIsoROIs",
                                  param = param, BPPARAM = BPPARAM)
              ## (3) collect the results.
              res <- .processResultList(resList,
                                        getProcHist = return.type == "xcmsSet",
                                        fnames = fileNames(object_mslevel))

              if (return.type == "list")
                  return(res$peaks)
              object <- .peaks_to_result(res, object, startDate, param, msLevel,
                                         object_mslevel)
              if (return.type == "xcmsSet")
                  as(object, "xcmsSet")
              else object
          })

## profMat method for XCMSnExp/OnDiskMSnExp.
#' @description \code{profMat}: creates a \emph{profile matrix}, which
#'     is a n x m matrix, n (rows) representing equally spaced m/z values (bins)
#'     and m (columns) the retention time of the corresponding scans. Each cell
#'     contains the maximum intensity measured for the specific scan and m/z
#'     values. See \code{\link{profMat}} for more details and description of
#'     the various binning methods.
#'
#' @param ... Additional parameters.
#' 
#' @return For \code{profMat}: a \code{list} with a the profile matrix
#'     \code{matrix} (or matrices if \code{fileIndex} was not specified or if
#'     \code{length(fileIndex) > 1}). See \code{\link{profile-matrix}} for
#'     general help and information about the profile matrix.
#'
#' @inheritParams profMat-xcmsSet
#'
#' @rdname XCMSnExp-class
setMethod("profMat", signature(object = "OnDiskMSnExp"), function(object,
                                                                  method = "bin",
                                                                  step = 0.1,
                                                                  baselevel = NULL,
                                                                  basespace = NULL,
                                                                  mzrange. = NULL,
                                                                  fileIndex,
                                                                  ...) {
    ## Subset the object by fileIndex.
    if (!missing(fileIndex)) {
        if (!is.numeric(fileIndex))
            stop("'fileIndex' has to be an integer.")
        if (!all(fileIndex %in% seq_along(fileNames(object))))
            stop("'fileIndex' has to be an integer between 1 and ",
                 length(fileNames(object)), "!")
        object <- filterFile(object, fileIndex)
    }
    ## Split it by file and bplapply over it to generate the profile matrix.
    theF <- factor(seq_along(fileNames(object)))
    theDots <- list(...)
    if (any(names(theDots) == "returnBreaks"))
        returnBreaks <- theDots$returnBreaks
    else
        returnBreaks <- FALSE
    res <- bplapply(splitByFile(object, f = theF), function(z, bmethod, bstep,
                                                            bbaselevel,
                                                            bbasespace,
                                                            bmzrange.,
                                                            breturnBreaks) {
        require(xcms, quietly = TRUE)
        ## Note: this is way faster than spectrapply with
        ## as.data.frame!
        sps <- spectra(z, BPPARAM = SerialParam())
        mzs <- lapply(sps, mz)
        vps <- lengths(mzs, use.names = FALSE)
        return(.createProfileMatrix(mz = unlist(mzs, use.names = FALSE),
                                    int = unlist(lapply(sps, intensity),
                                                 use.names = FALSE),
                                    valsPerSpect = vps,
                                    method = bmethod,
                                    step = bstep,
                                    baselevel = bbaselevel,
                                    basespace = bbasespace,
                                    mzrange. = bmzrange.,
                                    returnBreaks = breturnBreaks)
               )
    }, bmethod = method, bstep = step, bbaselevel = baselevel,
    bbasespace = basespace, bmzrange. = mzrange., breturnBreaks = returnBreaks)
    return(res)
})

#' @rdname adjustRtime-obiwarp
setMethod("adjustRtime",
          signature(object = "OnDiskMSnExp", param = "ObiwarpParam"),
          function(object, param, msLevel = 1L) {
              ## Filter for MS level, perform adjustment and if the object
              ## contains spectra from other MS levels too, adjust all raw
              ## rts based on the difference between adjusted and raw rts.
              object_sub <- filterMsLevel(object, msLevel = msLevel)
              if (length(object_sub) == 0)
                  stop("No spectra of MS level ", msLevel, " present")
              res <- xcms:::.obiwarp(object_sub, param = param)
              ## Adjust the retention time for spectra of all MS levels, if
              ## if there are some other than msLevel (issue #214).
              if (length(unique(msLevel(object))) !=
                  length(unique(msLevel(object_sub)))) {
                  message("Apply retention time correction performed on MS",
                          msLevel, " to spectra from all MS levels")
                  ## I need raw and adjusted rt for the adjusted spectra
                  ## and the raw rt of all.
                  rtime_all <- split(rtime(object), fromFile(object))
                  rtime_sub <- split(rtime(object_sub), fromFile(object_sub))
                  ## For loop is faster than lapply. No sense to do parallel
                  for (i in 1:length(rtime_all)) {
                      n_vals <- length(rtime_sub[[i]])
                      idx_below <- which(rtime_all[[i]] < rtime_sub[[i]][1])
                      if (length(idx_below))
                          vals_below <- rtime_all[[i]][idx_below]
                      idx_above <- which(rtime_all[[i]] >
                                         rtime_sub[[i]][n_vals])
                      if (length(idx_above))
                          vals_above <- rtime_all[[i]][idx_above]
                      ## Adjust the retention time. Note: this should be
                      ## OK even if values are not sorted.
                      adj_fun <- approxfun(x = rtime_sub[[i]], y = res[[i]])
                      rtime_all[[i]] <- adj_fun(rtime_all[[i]])
                      ## Adjust rtime < smallest adjusted rtime.
                      if (length(idx_below)) {
                          rtime_all[[i]][idx_below] <- vals_below +
                              res[[i]][1] - rtime_sub[[i]][1]
                      }
                      ## Adjust rtime > largest adjusted rtime
                      if (length(idx_above)) {
                          rtime_all[[i]][idx_above] <- vals_above +
                              res[[i]][n_vals] - rtime_sub[[i]][n_vals]
                      }
                  }
                  res <- rtime_all
              }
              res <- unlist(res, use.names = FALSE)
              sNames <- unlist(split(featureNames(object), fromFile(object)),
                               use.names = FALSE)
              names(res) <- sNames
              res <- res[featureNames(object)]
              res
          })

#' @rdname extractMsData-method
setMethod("extractMsData", signature(object = "OnDiskMSnExp"),
          function(object, rt, mz, msLevel = 1L) {
              .Deprecated(msg = paste0("Use of 'extractMsData' is deprecated.",
                                       " Please use 'as(x, \"data.frame\")'"))
              .extractMsData(object, rt = rt, mz = mz, msLevel = msLevel)
          })

setMethod("hasAdjustedRtime", signature(object = "OnDiskMSnExp"),
          function(object)
              FALSE
          )
