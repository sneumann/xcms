## Methods for MSnbase's OnDiskMSnExp and MSnExp objects.
#' @include DataClasses.R functions-OnDiskMSnExp.R do_findChromPeaks-functions.R


## Main roxygen documentation for the centWace feature detection is in
## DataClasses, before the definition of the CentWaveParam class.

## The centWave peak detection method for OnDiskMSnExp:
##' @title Chromatographic peak detection using the centWave method
##'
##' @description The \code{detectChromPeaks,OnDiskMSnExp,CentWaveParam} method
##' performs chromatographic peak detection using the \emph{centWave} algorithm
##' on all samples from an \code{\link[MSnbase]{OnDiskMSnExp}} object.
##' \code{\link[MSnbase]{OnDiskMSnExp}} objects encapsule all experiment specific
##' data and load the spectra data (mz and intensity values) on the fly from the
##' original files applying also all eventual data manipulations.
##'
##' @details Parallel processing (one process per sample) is supported and can
##' be configured either by the \code{BPPARAM} parameter or by globally defining
##' the parallel processing mode using the \code{\link[BiocParallel]{register}}
##' method from the \code{BiocParallel} package.
##'
##' @param object For \code{findChromPeaks}: an
##' \code{\link[MSnbase]{OnDiskMSnExp}}  object containing the MS- and all other
##' experiment-relevant data.
##'
##' For all other methods: a parameter object.
##'
##' @param param An \code{CentWaveParam} object containing all settings for the
##' centWave algorithm.
##'
##' @param BPPARAM A parameter class specifying if and how parallel processing
##' should be performed. It defaults to \code{\link[BiocParallel]{bpparam}}.
##' See documentation of the \code{BiocParallel} for more details. If parallel
##' processing is enables, peak detection is performed in parallel on several
##' of the input samples.
##'
##' @param return.type Character specifying what type of object the method should
##' return. Can be either \code{"XCMSnExp"} (default), \code{"list"} or
##' \code{"xcmsSet"}.
##'
##' @return For \code{findChromPeaks}: if \code{return.type = "XCMSnExp"} an
##' \code{\link{XCMSnExp}} object with the results of the peak detection.
##' If \code{return.type = "list"} a list of length equal to the number of
##' samples with matrices specifying the identified peaks.
##' If \code{return.type = "xcmsSet"} an \code{\linkS4class{xcmsSet}} object
##' with the results of the peak detection.
##'
##' @seealso \code{\link{XCMSnExp}} for the object containing the results of
##' the peak detection.
##'
##' @rdname findChromPeaks-centWave
setMethod("findChromPeaks",
          signature(object = "OnDiskMSnExp", param = "CentWaveParam"),
          function(object, param, BPPARAM = bpparam(), return.type = "XCMSnExp") {
              return.type <- match.arg(return.type, c("XCMSnExp", "list",
                                                      "xcmsSet"))
              startDate <- date()
              ## Restrict to MS1 data.
              object <- filterMsLevel(object, msLevel. = 1)
              ## Check if the data is centroided
              if (!isCentroided(object[[1]]))
                  warning("Your data appears to be not centroided! CentWave",
                          " works best on data in centroid mode.")
              ## (1) split the object per file.
              ## (2) use bplapply to do the peak detection.
              resList <- bplapply(lapply(1:length(fileNames(object)),
                                     filterFile, object = object),
                              FUN = findChromPeaks_OnDiskMSnExp,
                              method = "centWave", param = param,
                              BPPARAM = BPPARAM)
              ## (3) collect the results.
              res <- .processResultList(resList,
                                        getProcHist = return.type == "xcmsSet",
                                        fnames = fileNames(object))
              if (return.type == "XCMSnExp") {
                  ## Creating one XProcessHistory for all; eventually change
                  ## that later, but for now seems reasonable to have it in one,
                  ## since we're calling the method once on all.
                  xph <- XProcessHistory(param = param, date. = startDate,
                                         type. = .PROCSTEP.PEAK.DETECTION,
                                         fileIndex = 1:length(fileNames(object)))
                  object <- as(object, "XCMSnExp")
                  object@.processHistory <- list(xph)
                  if (hasAdjustedRtime(object) | hasFeatures(object))
                      object@msFeatureData <- new("MsFeatureData")
                  pks <- do.call(rbind, res$peaks)
                  if (length(pks) > 0)
                      chromPeaks(object) <- cbind(pks, is_filled = 0)
                  if (validObject(object))
                      return(object)
              }
              if (return.type == "list")
                  return(res$peaks)
              if (return.type == "xcmsSet") {
                  xs <- .pSet2xcmsSet(object)
                  peaks(xs) <- do.call(rbind, res$peaks)
                  xs@.processHistory <- res$procHist
                  OK <- .validProcessHistory(object)
                  if (!is.logical(OK))
                      stop(OK)
                  if (!any(colnames(pData(object)) == "class"))
                      message("Note: you might want to set/adjust the",
                              " 'sampclass' of the returned xcmSet object",
                              " before proceeding with the analysis.")
                  return(xs)
              }
          })


## ## The centWave peak detection method for MSnExp:
## ##' @title Chromatographic peak detection using the centWave method
## ##'
## ##' @description The \code{findChromPeaks,MSnExp,CentWaveParam} method performs
## ##' peak detection using the \emph{centWave} algorithm on all samples from
## ##' an \code{\link[MSnbase]{MSnExp}} object. These objects contain mz and
## ##' intensity values of all spectra hence no additional data input from the
## ##' original files is required.
## ##'
## ##' @rdname findChromPeaks-centWave
## setMethod("findChromPeaks",
##           signature(object = "MSnExp", param = "CentWaveParam"),
##           function(object, param, BPPARAM = bpparam(), return.type = "list") {
##               return.type <- match.arg(return.type, c("list", "xcmsSet"))
##               ## Restrict to MS1 data.
##               ## Man, that's too slow! We're doing the MS1 restriction below.
##               ## object <- filterMsLevel(object, msLevel. = 1)
##               ## (1) split the spectra per file - this means we have a second
##               ##     copy of the data, but there is no way around that as
##               ##     filterFile is pretty slow on MSnExp.
##               ms1_idx <- which(unname(msLevel(object)) == 1)
##               if (length(ms1_idx) == 0)
##                   stop("No MS1 spectra available for chromatographic peak",
##                        " detection!")
##               ## Check if the data is centroided
##               if (!isCentroided(object[[ms1_idx[1]]]))
##                   warning("Your data appears to be not centroided! CentWave",
##                           " works best on data in centroid mode.")
##               spect_list <- split(spectra(object)[ms1_idx],
##                                   fromFile(object)[ms1_idx])
##               ## (2) use bplapply to do the peak detection.
##               resList <- bplapply(spect_list, function(z) {
##                   findChromPeaks_Spectrum_list(z,
##                                                method = "centWave",
##                                                param = param)
##               }, BPPARAM = BPPARAM)
##               ## (3) collect the results.
##               res <- .processResultList(resList,
##                                         getProcHist = return.type != "list",
##                                         fnames = fileNames(object))
##               if (return.type == "list")
##                   return(res$peaks)
##               if (return.type == "xcmsSet") {
##                   xs <- .pSet2xcmsSet(object)
##                   peaks(xs) <- do.call(rbind, res$peaks)
##                   xs@.processHistory <- res$procHist
##                   OK <- .validProcessHistory(xs)
##                   if (!is.logical(OK))
##                       stop(OK)
##                   if (!any(colnames(pData(object)) == "class"))
##                       message("Note: you might want to set/adjust the",
##                               " 'sampclass' of the returned xcmSet object",
##                               " before proceeding with the analysis.")
##                   return(xs)
##               }
##           })

## The matchedFilter peak detection method for OnDiskMSnExp:
##' @title Peak detection in the chromatographic time domain
##'
##' @description The \code{findChromPeaks,OnDiskMSnExp,MatchedFilterParam}
##' method performs peak detection using the \emph{matchedFilter} algorithm
##' on all samples from an \code{\link[MSnbase]{OnDiskMSnExp}} object.
##' \code{\link[MSnbase]{OnDiskMSnExp}} objects encapsule all experiment specific
##' data and load the spectra data (mz and intensity values) on the fly from the
##' original files applying also all eventual data manipulations.
##'
##' @details Parallel processing (one process per sample) is supported and can
##' be configured either by the \code{BPPARAM} parameter or by globally defining
##' the parallel processing mode using the \code{\link[BiocParallel]{register}}
##' method from the \code{BiocParallel} package.
##' 
##' @param object For \code{findChromPeaks}: an
##' \code{\link[MSnbase]{OnDiskMSnExp}} object containing the MS- and all other
##' experiment-relevant data.
##'
##' For all other methods: a parameter object.
##'
##' @param param An \code{MatchedFilterParam} object containing all settings for
##' the matchedFilter algorithm.
##'
##' @inheritParams findChromPeaks-centWave
##'
##' @return For \code{findChromPeaks}: if \code{return.type = "XCMSnExp"} an
##' \code{\link{XCMSnExp}} object with the results of the peak detection.
##' If \code{return.type = "list"} a list of length equal to the number of
##' samples with matrices specifying the identified peaks.
##' If \code{return.type = "xcmsSet"} an \code{\linkS4class{xcmsSet}} object
##' with the results of the peak detection.
##'
##' @seealso \code{\link{XCMSnExp}} for the object containing the results of
##' the chromatographic peak detection.
##'
##' @rdname findChromPeaks-matchedFilter
setMethod("findChromPeaks",
          signature(object = "OnDiskMSnExp", param = "MatchedFilterParam"),
          function(object, param, BPPARAM = bpparam(), return.type = "XCMSnExp") {
              return.type <- match.arg(return.type, c("XCMSnExp", "list",
                                                      "xcmsSet"))
              startDate <- date()
              ## Restrict to MS1 data.
              object <- filterMsLevel(object, msLevel. = 1)
              ## (1) split the object per file.
              ## (2) use bplapply to do the peak detection.
              resList <- bplapply(lapply(1:length(fileNames(object)),
                                     filterFile, object = object),
                              FUN = findChromPeaks_OnDiskMSnExp,
                              method = "matchedFilter", param = param,
                              BPPARAM = BPPARAM)
              ## (3) collect the results.
              res <- .processResultList(resList,
                                        getProcHist = return.type == "xcmsSet",
                                        fnames = fileNames(object))
              if (return.type == "XCMSnExp") {
                  ## Creating one XProcessHistory for all; eventually change
                  ## that later, but for now seems reasonable to have it in one,
                  ## since we're calling the method once on all.
                  xph <- XProcessHistory(param = param, date. = startDate,
                                         type. = .PROCSTEP.PEAK.DETECTION,
                                         fileIndex = 1:length(fileNames(object)))
                  object <- as(object, "XCMSnExp")
                  object@.processHistory <- list(xph)
                  if (hasAdjustedRtime(object) | hasFeatures(object))
                      object@msFeatureData <- new("MsFeatureData")
                  pks <- do.call(rbind, res$peaks)
                  if (length(pks) > 0)
                      chromPeaks(object) <- cbind(pks, is_filled = 0)
                  ## ## chromPeaks(object) <- do.call(rbind, res$peaks)
                  ## chromPeaks(object) <- cbind(do.call(rbind, res$peaks),
                  ##                             is_filled = 0)
                  if (validObject(object))
                      return(object)
              }
              if (return.type == "list")
                  return(res$peaks)
              if (return.type == "xcmsSet") {
                  xs <- .pSet2xcmsSet(object)
                  peaks(xs) <- do.call(rbind, res$peaks)
                  xs@.processHistory <- res$procHist
                  OK <- .validProcessHistory(object)
                  if (!is.logical(OK))
                      stop(OK)
                  if (!any(colnames(pData(object)) == "class"))
                      message("Note: you might want to set/adjust the",
                              " 'sampclass' of the returned xcmSet object",
                              " before proceeding with the analysis.")
                  return(xs)
              }
          })

## ##' @title Peak detection in the chromatographic time domain
## ##'
## ##' @description The \code{findChromPeaks,MSnExp,MatchedFilterParam} method
## ##' performs peak detection using the \emph{matchedFilter} method on all
## ##' samples from an \code{\link[MSnbase]{MSnExp}} object. These objects contain
## ##' mz and intensity values of all spectra hence no additional
## ##' data input from the original files is required.
## ##'
## ##' @rdname findChromPeaks-matchedFilter
## setMethod("findChromPeaks",
##           signature(object = "MSnExp", param = "MatchedFilterParam"),
##           function(object, param, BPPARAM = bpparam(), return.type = "list") {
##               return.type <- match.arg(return.type, c("list", "xcmsSet"))
##               ms1_idx <- which(unname(msLevel(object)) == 1)
##               if (length(ms1_idx) == 0)
##                   stop("No MS1 spectra available for chromatographic peak",
##                        " detection!")
##               spect_list <- split(spectra(object)[ms1_idx],
##                                   fromFile(object)[ms1_idx])
##               resList <- bplapply(spect_list, function(z) {
##                   findChromPeaks_Spectrum_list(z,
##                                                method = "matchedFilter",
##                                                param = param)
##               }, BPPARAM = BPPARAM)
##               res <- .processResultList(resList,
##                                         getProcHist = return.type != "list",
##                                         fnames = fileNames(object))
##               if (return.type == "list")
##                   return(res$peaks)
##               if (return.type == "xcmsSet") {
##                   xs <- .pSet2xcmsSet(object)
##                   peaks(xs) <- do.call(rbind, res$peaks)
##                   xs@.processHistory <- res$procHist
##                   OK <- .validProcessHistory(xs)
##                   if (!is.logical(OK))
##                       stop(OK)
##                   if (!any(colnames(pData(object)) == "class"))
##                       message("Note: you might want to set/adjust the",
##                               " 'sampclass' of the returned xcmSet object",
##                               " before proceeding with the analysis.")
##                   return(xs)
##               }
##           })

## massifquant
## The massifquant peak detection method for OnDiskMSnExp:
##' @title Chromatographic peak detection using the massifquant method
##'
##' @description The \code{findChromPeaks,OnDiskMSnExp,MassifquantParam}
##' method performs chromatographic peak detection using the \emph{massifquant}
##' algorithm on all samples from an \code{\link[MSnbase]{OnDiskMSnExp}} object.
##' \code{\link[MSnbase]{OnDiskMSnExp}} objects encapsule all experiment specific
##' data and load the spectra data (mz and intensity values) on the fly from the
##' original files applying also all eventual data manipulations.
##'
##' @details Parallel processing (one process per sample) is supported and can
##' be configured either by the \code{BPPARAM} parameter or by globally defining
##' the parallel processing mode using the \code{\link[BiocParallel]{register}}
##' method from the \code{BiocParallel} package.
##'
##' @param object For \code{findChromPeaks}: an
##' \code{\link[MSnbase]{OnDiskMSnExp}} object containing the MS- and all other
##' experiment-relevant data.
##'
##' For all other methods: a parameter object.
##'
##' @param param An \code{MassifquantParam} object containing all settings for
##' the massifquant algorithm.
##'
##' @inheritParams findChromPeaks-centWave
##'
##' @return For \code{findChromPeaks}: if \code{return.type = "XCMSnExp"} an
##' \code{\link{XCMSnExp}} object with the results of the peak detection.
##' If \code{return.type = "list"} a list of length equal to the number of
##' samples with matrices specifying the identified peaks.
##' If \code{return.type = "xcmsSet"} an \code{\linkS4class{xcmsSet}} object
##' with the results of the peak detection.
##'
##' @seealso \code{\link{XCMSnExp}} for the object containing the results of
##' the peak detection.
##'
##' @rdname findChromPeaks-massifquant
setMethod("findChromPeaks",
          signature(object = "OnDiskMSnExp", param = "MassifquantParam"),
          function(object, param, BPPARAM = bpparam(), return.type = "XCMSnExp") {
              return.type <- match.arg(return.type, c("XCMSnExp", "list",
                                                      "xcmsSet"))
              startDate <- date()
              ## Restrict to MS1 data.
              object <- filterMsLevel(object, msLevel. = 1)
              ## (1) split the object per file.
              ## (2) use bplapply to do the peaks detection.
              resList <- bplapply(lapply(1:length(fileNames(object)),
                                     filterFile, object = object),
                              FUN = findChromPeaks_OnDiskMSnExp,
                              method = "massifquant", param = param,
                              BPPARAM = BPPARAM)
              ## (3) collect the results.
              res <- .processResultList(resList,
                                        getProcHist = return.type == "xcmsSet",
                                        fnames = fileNames(object))
              if (return.type == "XCMSnExp") {
                  ## Creating one XProcessHistory for all; eventually change
                  ## that later, but for now seems reasonable to have it in one,
                  ## since we're calling the method once on all.
                  xph <- XProcessHistory(param = param, date. = startDate,
                                         type. = .PROCSTEP.PEAK.DETECTION,
                                         fileIndex = 1:length(fileNames(object)))
                  object <- as(object, "XCMSnExp")
                  object@.processHistory <- list(xph)
                  if (hasAdjustedRtime(object) | hasFeatures(object))
                      object@msFeatureData <- new("MsFeatureData")
                  pks <- do.call(rbind, res$peaks)
                  if (length(pks) > 0)
                      chromPeaks(object) <- cbind(pks, is_filled = 0)
                  ## ## chromPeaks(object) <- do.call(rbind, res$peaks)
                  ## chromPeaks(object) <- cbind(do.call(rbind, res$peaks),
                  ##                             is_filled = 0)
                  if (validObject(object))
                      return(object)
              }
              if (return.type == "list")
                  return(res$peaks)
              if (return.type == "xcmsSet") {
                  xs <- .pSet2xcmsSet(object)
                  peaks(xs) <- do.call(rbind, res$peaks)
                  xs@.processHistory <- res$procHist
                  OK <- .validProcessHistory(object)
                  if (!is.logical(OK))
                      stop(OK)
                  if (!any(colnames(pData(object)) == "class"))
                      message("Note: you might want to set/adjust the",
                              " 'sampclass' of the returned xcmSet object",
                              " before proceeding with the analysis.")
                  return(xs)
              }
          })


## ##' @title Chromatographic peak detection using the massifquant method
## ##'
## ##' @description The \code{findChromPeaks,MSnExp,MassifquantParam} method
## ##' performs chromatographic peak detection using the \emph{massifquant} method
## ##' on all samples from an \code{\link[MSnbase]{MSnExp}} object. These objects
## ##' contain mz and intensity values of all spectra hence no additional
## ##' data input from the original files is required.
## ##'
## ##' @rdname findChromPeaks-massifquant
## setMethod("findChromPeaks",
##           signature(object = "MSnExp", param = "MassifquantParam"),
##           function(object, param, BPPARAM = bpparam(), return.type = "list") {
##               return.type <- match.arg(return.type, c("list", "xcmsSet"))
##               ms1_idx <- which(unname(msLevel(object)) == 1)
##               if (length(ms1_idx) == 0)
##                   stop("No MS1 spectra available for chromatographic peak",
##                        " detection!")
##               spect_list <- split(spectra(object)[ms1_idx],
##                                   fromFile(object)[ms1_idx])
##               resList <- bplapply(spect_list, function(z) {
##                   findChromPeaks_Spectrum_list(z,
##                                                method = "massifquant",
##                                                param = param)
##               }, BPPARAM = BPPARAM)
##               res <- .processResultList(resList,
##                                         getProcHist = return.type != "list",
##                                         fnames = fileNames(object))
##               if (return.type == "list")
##                   return(res$peaks)
##               if (return.type == "xcmsSet") {
##                   xs <- .pSet2xcmsSet(object)
##                   peaks(xs) <- do.call(rbind, res$peaks)
##                   xs@.processHistory <- res$procHist
##                   OK <- .validProcessHistory(xs)
##                   if (!is.logical(OK))
##                       stop(OK)
##                   if (!any(colnames(pData(object)) == "class"))
##                       message("Note: you might want to set/adjust the",
##                               " 'sampclass' of the returned xcmSet object",
##                               " before proceeding with the analysis.")
##                   return(xs)
##               }
##           })


## MSW
## The MSW peak detection method for OnDiskMSnExp:
##' @title Single-spectrum non-chromatography MS data peak detection
##'
##' @description The \code{findChromPeaks,OnDiskMSnExp,MSWParam}
##' method performs peak detection in single-spectrum non-chromatography MS
##' data using functionality from the \code{MassSpecWavelet} package on all
##' samples from an \code{\link[MSnbase]{OnDiskMSnExp}} object.
##' \code{\link[MSnbase]{OnDiskMSnExp}} objects encapsule all experiment specific
##' data and load the spectra data (mz and intensity values) on the fly from the
##' original files applying also all eventual data manipulations.
##'
##' @details Parallel processing (one process per sample) is supported and can
##' be configured either by the \code{BPPARAM} parameter or by globally defining
##' the parallel processing mode using the \code{\link[BiocParallel]{register}}
##' method from the \code{BiocParallel} package.
##'
##' @param object For \code{findChromPeaks}: an
##' \code{\link[MSnbase]{OnDiskMSnExp}} object containing the MS- and all other
##' experiment-relevant data.
##'
##' For all other methods: a parameter object.
##'
##' @param param An \code{MSWParam} object containing all settings for
##' the algorithm.
##'
##' @inheritParams findChromPeaks-centWave
##'
##' @return For \code{findChromPeaks}: if \code{return.type = "XCMSnExp"} an
##' \code{\link{XCMSnExp}} object with the results of the peak detection.
##' If \code{return.type = "list"} a list of length equal to the number of
##' samples with matrices specifying the identified peaks.
##' If \code{return.type = "xcmsSet"} an \code{\linkS4class{xcmsSet}} object
##' with the results of the detection.
##'
##' @seealso \code{\link{XCMSnExp}} for the object containing the results of
##' the peak detection.
##'
##' @rdname findPeaks-MSW
setMethod("findChromPeaks",
          signature(object = "OnDiskMSnExp", param = "MSWParam"),
          function(object, param, BPPARAM = bpparam(), return.type = "XCMSnExp") {
              return.type <- match.arg(return.type, c("XCMSnExp", "list",
                                                      "xcmsSet"))
              startDate <- date()
              ## Restrict to MS1 data.
              object <- filterMsLevel(object, msLevel. = 1)

              rts <- split(rtime(object), f = fromFile(object))
              if (any(lengths(rts)) > 1)
                  stop("The MSW method can only be applied to single spectrum,",
                       " non-chromatogrphic, files (i.e. with a single ",
                       "retention time).")
              ## (1) split the object per file.
              ## (2) use bplapply to do the peak detection.
              resList <- bplapply(lapply(1:length(fileNames(object)),
                                     filterFile, object = object),
                              FUN = findPeaks_MSW_OnDiskMSnExp,
                              method = "MSW", param = param,
                              BPPARAM = BPPARAM)
              ## (3) collect the results.
              res <- .processResultList(resList,
                                        getProcHist = return.type == "xcmsSet",
                                        fnames = fileNames(object))
              if (return.type == "XCMSnExp") {
                  ## Creating one XProcessHistory for all; eventually change
                  ## that later, but for now seems reasonable to have it in one,
                  ## since we're calling the method once on all.
                  xph <- XProcessHistory(param = param, date. = startDate,
                                         type. = .PROCSTEP.PEAK.DETECTION,
                                         fileIndex = 1:length(fileNames(object)))
                  object <- as(object, "XCMSnExp")
                  object@.processHistory <- list(xph)
                  if (hasAdjustedRtime(object) | hasFeatures(object))
                      object@msFeatureData <- new("MsFeatureData")
                  pks <- do.call(rbind, res$peaks)
                  if (length(pks) > 0)
                      chromPeaks(object) <- cbind(pks, is_filled = 0)
                  ## ## chromPeaks(object) <- do.call(rbind, res$peaks)
                  ## chromPeaks(object) <- cbind(do.call(rbind, res$peaks),
                  ##                             is_filled = 0)
                  if (validObject(object))
                      return(object)
              }
              if (return.type == "list")
                  return(res$peaks)
              if (return.type == "xcmsSet") {
                  xs <- .pSet2xcmsSet(object)
                  peaks(xs) <- do.call(rbind, res$peaks)
                  xs@.processHistory <- res$procHist
                  OK <- .validProcessHistory(object)
                  if (!is.logical(OK))
                      stop(OK)
                  if (!any(colnames(pData(object)) == "class"))
                      message("Note: you might want to set/adjust the",
                              " 'sampclass' of the returned xcmSet object",
                              " before proceeding with the analysis.")
                  return(xs)
              }
          })

## ##' @title Single-spectrum non-chromatography MS data peak detection
## ##'
## ##' @description The \code{findChromPeaks,MSnExp,MSWParam} method
## ##' performs peak detection in single-spectrum non-chromatography MS
## ##' data using functionality from the \code{MassSpecWavelet} package on all
## ##' samples from an \code{\link[MSnbase]{MSnExp}} object. These objects contain
## ##' mz and intensity values of all spectra hence no additional
## ##' data input from the original files is required.
## ##'
## ##' @rdname findPeaks-MSW
## setMethod("findChromPeaks",
##           signature(object = "MSnExp", param = "MSWParam"),
##           function(object, param, BPPARAM = bpparam(), return.type = "list") {
##               return.type <- match.arg(return.type, c("list", "xcmsSet"))
##               ms1_idx <- which(unname(msLevel(object)) == 1)
##               if (length(ms1_idx) == 0)
##                   stop("No MS1 spectra available for chromatographic peak",
##                        " detection!")
##               spect_list <- split(spectra(object)[ms1_idx],
##                                   fromFile(object)[ms1_idx])
##               resList <- bplapply(spect_list, function(z) {
##                   findPeaks_MSW_Spectrum_list(z,
##                                               method = "MSW",
##                                               param = param)
##               }, BPPARAM = BPPARAM)
##               res <- .processResultList(resList,
##                                         getProcHist = return.type != "list",
##                                         fnames = fileNames(object))
##               if (return.type == "list")
##                   return(res$peaks)
##               if (return.type == "xcmsSet") {
##                   xs <- .pSet2xcmsSet(object)
##                   peaks(xs) <- do.call(rbind, res$peaks)
##                   xs@.processHistory <- res$procHist
##                   OK <- .validProcessHistory(xs)
##                   if (!is.logical(OK))
##                       stop(OK)
##                   if (!any(colnames(pData(object)) == "class"))
##                       message("Note: you might want to set/adjust the",
##                               " 'sampclass' of the returned xcmSet object",
##                               " before proceeding with the analysis.")
##                   return(xs)
##               }
##           })

## The centWave with predicted isotope peak detection method for OnDiskMSnExp:
##' @title Two-step centWave peak detection considering also isotopes
##'
##' @description The \code{findChromPeaks,OnDiskMSnExp,CentWavePredIsoParam} method
##' performs a two-step centWave-based chromatographic peak detection on all
##' samples from an \code{\link[MSnbase]{OnDiskMSnExp}} object.
##' \code{\link[MSnbase]{OnDiskMSnExp}} objects encapsule all experiment specific
##' data and load the spectra data (mz and intensity values) on the fly from
##' the original files applying also all eventual data manipulations.
##'
##' @details Parallel processing (one process per sample) is supported and can
##' be configured either by the \code{BPPARAM} parameter or by globally defining
##' the parallel processing mode using the \code{\link[BiocParallel]{register}}
##' method from the \code{BiocParallel} package.
##'
##' @param param An \code{CentWavePredIsoParam} object with the settings for the
##' chromatographic peak detection algorithm.
##' @inheritParams findChromPeaks-centWave
##' 
##' @return For \code{findChromPeaks}: if \code{return.type = "XCMSnExp"} an
##' \code{\link{XCMSnExp}} object with the results of the peak detection.
##' If \code{return.type = "list"} a list of length equal to the number of
##' samples with matrices specifying the identified peaks.
##' If \code{return.type = "xcmsSet"} an \code{\linkS4class{xcmsSet}} object
##' with the results of the peak detection.
##'
##' @seealso \code{\link{XCMSnExp}} for the object containing the results of
##' the peak detection.
##'
##' @rdname findChromPeaks-centWaveWithPredIsoROIs
setMethod("findChromPeaks",
          signature(object = "OnDiskMSnExp", param = "CentWavePredIsoParam"),
          function(object, param, BPPARAM = bpparam(), return.type = "XCMSnExp") {
              return.type <- match.arg(return.type, c("XCMSnExp", "list",
                                                      "xcmsSet"))
              startDate <- date()
              ## Restrict to MS1 data.
              object <- filterMsLevel(object, msLevel. = 1)
              ## Check if the data is centroided
              if (!isCentroided(object[[1]]))
                  warning("Your data appears to be not centroided! CentWave",
                          " works best on data in centroid mode.")
              ## (1) split the object per file.
              ## (2) use bplapply to do the peak detection.
              resList <- bplapply(lapply(1:length(fileNames(object)),
                                     filterFile, object = object),
                              FUN = findChromPeaks_OnDiskMSnExp,
                              method = "centWaveWithPredIsoROIs", param = param,
                              BPPARAM = BPPARAM)
              ## (3) collect the results.
              res <- .processResultList(resList,
                                        getProcHist = return.type == "xcmsSet",
                                        fnames = fileNames(object))
              if (return.type == "XCMSnExp") {
                  ## Creating one XProcessHistory for all; eventually change
                  ## that later, but for now seems reasonable to have it in one,
                  ## since we're calling the method once on all.
                  xph <- XProcessHistory(param = param, date. = startDate,
                                         type. = .PROCSTEP.PEAK.DETECTION,
                                         fileIndex = 1:length(fileNames(object)))
                  object <- as(object, "XCMSnExp")
                  object@.processHistory <- list(xph)
                  if (hasAdjustedRtime(object) | hasFeatures(object))
                      object@msFeatureData <- new("MsFeatureData")
                  pks <- do.call(rbind, res$peaks)
                  if (length(pks) > 0)
                      chromPeaks(object) <- cbind(pks, is_filled = 0)
                  ## ## chromPeaks(object) <- do.call(rbind, res$peaks)
                  ## chromPeaks(object) <- cbind(do.call(rbind, res$peaks),
                  ##                             is_filled = 0)
                  if (validObject(object))
                      return(object)
              }
              if (return.type == "list")
                  return(res$peaks)
              if (return.type == "xcmsSet") {
                  xs <- .pSet2xcmsSet(object)
                  peaks(xs) <- do.call(rbind, res$peaks)
                  xs@.processHistory <- res$procHist
                  OK <- .validProcessHistory(object)
                  if (!is.logical(OK))
                      stop(OK)
                  if (!any(colnames(pData(object)) == "class"))
                      message("Note: you might want to set/adjust the",
                              " 'sampclass' of the returned xcmSet object",
                              " before proceeding with the analysis.")
                  return(xs)
              }
          })


## ## The centWave with predicted isotope peak detection method for MSnExp:
## ##' @title Two-step centWave peak detection considering also isotopes
## ##'
## ##' @description The \code{findChromPeaks,MSnExp,CentWavePredIsoParam} method
## ##' performs a two-step centWave-based peak detection on all samples from
## ##' an \code{\link[MSnbase]{MSnExp}} object. These objects contain mz and
## ##' intensity values of all spectra hence no additional data input from the
## ##' original files is required.
## ##'
## ##' @rdname findChromPeaks-centWaveWithPredIsoROIs
## setMethod("findChromPeaks",
##           signature(object = "MSnExp", param = "CentWavePredIsoParam"),
##           function(object, param, BPPARAM = bpparam(), return.type = "list") {
##               return.type <- match.arg(return.type, c("list", "xcmsSet"))
##               ## Restrict to MS1 data.
##               ## Man, that's too slow! We're doing the MS1 restriction below.
##               ## object <- filterMsLevel(object, msLevel. = 1)
##               ## (1) split the spectra per file - this means we have a second
##               ##     copy of the data, but there is no way around that as
##               ##     filterFile is pretty slow on MSnExp.
##               ms1_idx <- which(unname(msLevel(object)) == 1)
##               if (length(ms1_idx) == 0)
##                   stop("No MS1 spectra available for chromatographic peak",
##                        " detection!")
##               ## Check if the data is centroided
##               if (!isCentroided(object[[ms1_idx[1]]]))
##                   warning("Your data appears to be not centroided! CentWave",
##                           " works best on data in centroid mode.")
##               spect_list <- split(spectra(object)[ms1_idx],
##                                   fromFile(object)[ms1_idx])
##               ## (2) use bplapply to do the peak detection.
##               resList <- bplapply(spect_list, function(z) {
##                   findChromPeaks_Spectrum_list(z,
##                                                method = "centWaveWithPredIsoROIs",
##                                                param = param)
##               }, BPPARAM = BPPARAM)
##               ## (3) collect the results.
##               res <- .processResultList(resList,
##                                         getProcHist = return.type != "list",
##                                         fnames = fileNames(object))
##               if (return.type == "list")
##                   return(res$peaks)
##               if (return.type == "xcmsSet") {
##                   xs <- .pSet2xcmsSet(object)
##                   peaks(xs) <- do.call(rbind, res$peaks)
##                   xs@.processHistory <- res$procHist
##                   OK <- .validProcessHistory(xs)
##                   if (!is.logical(OK))
##                       stop(OK)
##                   if (!any(colnames(pData(object)) == "class"))
##                       message("Note: you might want to set/adjust the",
##                               " 'sampclass' of the returned xcmSet object",
##                               " before proceeding with the analysis.")
##                   return(xs)
##               }
##           })

## profMat method for XCMSnExp/OnDiskMSnExp.
##' @description \code{profMat}: creates a \emph{profile matrix}, which
##' is a n x m matrix, n (rows) representing equally spaced m/z values (bins) and
##' m (columns) the retention time of the corresponding scans. Each cell contains
##' the maximum intensity measured for the specific scan and m/z values. See
##' \code{\link{profMat}} for more details and description of the various binning
##' methods.
##'
##' @param ... Additional parameters.
##' 
##' @return For \code{profMat}: a \code{list} with a the profile matrix
##' \code{matrix} (or matrices if \code{fileIndex} was not specified or if
##' \code{length(fileIndex) > 1}). See \code{\link{profile-matrix}} for general
##' help and information about the profile matrix.
##'
##' @inheritParams profMat-xcmsSet
##'
##' @rdname XCMSnExp-class
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

##' @rdname adjustRtime-obiwarp
setMethod("adjustRtime",
          signature(object = "OnDiskMSnExp", param = "ObiwarpParam"),
          function(object, param) {
              res <- .obiwarp(object, param = param)
              res <- unlist(res, use.names = FALSE)
              sNames <- unlist(split(featureNames(object), fromFile(object)),
                               use.names = FALSE)
              names(res) <- sNames
              res <- res[featureNames(object)]
              return(res)
          })

#' @rdname extractChromatograms-method
setMethod("extractChromatograms",
          signature(object = "OnDiskMSnExp"),
          function(object, rt, mz, aggregationFun = "sum") {
              if (!missing(rt)) {
                  if (is.null(ncol(rt)))
                      rt <- matrix(range(rt), ncol = 2, nrow = 1)
              }
              if (!missing(mz)) {
                  if (is.null(ncol(mz)))
                      mz <- matrix(range(mz), ncol = 2, nrow = 1)
              }
              ## return(.extractChromatogram(x = object, rt = rt, mz = mz,
              ##                             aggregationFun = aggregationFun))
              .extractMultipleChromatograms(object, rt = rt, mz = mz,
                                            aggregationFun = aggregationFun)
          })
