## Methods for MSnbase's OnDiskMSnExp and MSnExp objects.
#' @include functions-OnDiskMSnExp.R


## Main roxygen documentation for the centWace feature detection is in
## DataClasses, before the definition of the CentWaveParam class.

## The centWave feature detection method for OnDiskMSnExp:
##' @title Feature detection using the centWave method
##'
##' @description The \code{detectFeatures,OnDiskMSnExp,CentWaveParam} method
##' performs feature detection using the \emph{centWave} algorithm on all
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
##' @param object For \code{detectFeatures}: Either an
##' \code{\link[MSnbase]{OnDiskMSnExp}} or a \code{\link[MSnbase]{MSnExp}}
##' object containing the MS- and all other experiment-relevant data.
##'
##' For all other methods: a parameter object.
##'
##' @param param An \code{CentWaveParam} object containing all settings for the
##' centWave algorithm.
##'
##' @param BPPARAM A parameter class specifying if and how parallel processing
##' should be performed. It defaults to \code{\link[BiocParallel]{bpparam}}.
##' See documentation of the \code{BiocParallel} for more details. If parallel
##' processing is enables, feature detection is performed in parallel on several
##' of the input samples.
##'
##' @param return.type Character specifying what type of object the method should
##' return. Can be either \code{"list"} or \code{"xcmsSet"}.
##'
##' @return For \code{detectFeatures}: if \code{return.type = "list"} a list of
##' length equal to the number of samples with matrices specifying the
##' identified features/peaks. If \code{return.type = "xcmsSet"} an
##' \code{\linkS4class{xcmsSet}} object with the results of the feature
##' detection.
##'
##' @rdname featureDetection-centWave
setMethod("detectFeatures",
          signature(object = "OnDiskMSnExp", param = "CentWaveParam"),
          function(object, param, BPPARAM = bpparam(), return.type = "list") {
              return.type <- match.arg(return.type, c("list", "xcmsSet"))
              ## Restrict to MS1 data.
              object <- filterMsLevel(object, msLevel. = 1)
              ## Check if the data is centroided
              if (!isCentroided(object[[1]]))
                  warning("Your data appears to be not centroided! CentWave",
                          " works best on data in centroid mode.")
              ## (1) split the object per file.
              ## (2) use bplapply to do the feature detection.
              resList <- bplapply(lapply(1:length(fileNames(object)),
                                     filterFile, object = object),
                              FUN = detectFeatures_OnDiskMSnExp,
                              method = "centWave", param = param)
              ## (3) collect the results.
              res <- .processResultList(resList,
                                        getProcHist = return.type != "list",
                                        fnames = fileNames(object))
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


## The centWave feature detection method for MSnExp:
##' @title Feature detection using the centWave method
##'
##' @description The \code{detectFeatures,MSnExp,CentWaveParam} method performs
##' feature detection using the \emph{centWave} algorithm on all samples from
##' an \code{\link[MSnbase]{MSnExp}} object. These objects contain mz and
##' intensity values of all spectra hence no additional data input from the
##' original files is required.
##'
##' @rdname featureDetection-centWave
setMethod("detectFeatures",
          signature(object = "MSnExp", param = "CentWaveParam"),
          function(object, param, BPPARAM = bpparam(), return.type = "list") {
              return.type <- match.arg(return.type, c("list", "xcmsSet"))
              ## Restrict to MS1 data.
              ## Man, that's too slow! We're doing the MS1 restriction below.
              ## object <- filterMsLevel(object, msLevel. = 1)
              ## (1) split the spectra per file - this means we have a second
              ##     copy of the data, but there is no way around that as
              ##     filterFile is pretty slow on MSnExp.
              ms1_idx <- which(unname(msLevel(object)) == 1)
              if (length(ms1_idx) == 0)
                  stop("No MS1 spectra available for feature detection!")
              ## Check if the data is centroided
              if (!isCentroided(object[[ms1_idx[1]]]))
                  warning("Your data appears to be not centroided! CentWave",
                          " works best on data in centroid mode.")
              spect_list <- split(spectra(object)[ms1_idx],
                                  fromFile(object)[ms1_idx])
              ## (2) use bplapply to do the feature detection.
              resList <- bplapply(spect_list, function(z) {
                  detectFeatures_Spectrum_list(z,
                                               method = "centWave",
                                               param = param)
              }, BPPARAM = BPPARAM)
              ## (3) collect the results.
              res <- .processResultList(resList,
                                        getProcHist = return.type != "list",
                                        fnames = fileNames(object))
              if (return.type == "list")
                  return(res$peaks)
              if (return.type == "xcmsSet") {
                  xs <- .pSet2xcmsSet(object)
                  peaks(xs) <- do.call(rbind, res$peaks)
                  xs@.processHistory <- res$procHist
                  OK <- .validProcessHistory(xs)
                  if (!is.logical(OK))
                      stop(OK)
                  if (!any(colnames(pData(object)) == "class"))
                      message("Note: you might want to set/adjust the",
                              " 'sampclass' of the returned xcmSet object",
                              " before proceeding with the analysis.")
                  return(xs)
              }
          })

## The matchedFilter feature detection method for OnDiskMSnExp:
##' @title Peak detection in the chromatographic time domain
##'
##' @description The \code{detectFeatures,OnDiskMSnExp,MatchedFilterParam}
##' method performs feature detection using the \emph{matchedFilter} algorithm
##' on all samples from an \code{\link[MSnbase]{OnDiskMSnExp}} object.
##' \code{\link[MSnbase]{OnDiskMSnExp}} objects encapsule all experiment specific
##' data and load the spectra data (mz and intensity values) on the fly from the
##' original files applying also all eventual data manipulations.
##'
##' @details Parallel processing (one process per sample) is supported and can
##' be configured either by the \code{BPPARAM} parameter or by globally defining
##' the parallel processing mode using the \code{\link[BiocParallel]{register}}
##' method from the \code{BiocParallel} package.

##' @param object For \code{detectFeatures}: Either an
##' \code{\link[MSnbase]{OnDiskMSnExp}} or a \code{\link[MSnbase]{MSnExp}}
##' object containing the MS- and all other experiment-relevant data.
##'
##' For all other methods: a parameter object.
##'
##' @param param An \code{MatchedFilterParam} object containing all settings for
##' the matchedFilter algorithm.
##'
##' @inheritParams featureDetection-centWave
##'
##' @return For \code{detectFeatures}: if \code{return.type = "list"} a list of
##' length equal to the number of samples with matrices specifying the
##' identified features/peaks. If \code{return.type = "xcmsSet"} an
##' \code{\linkS4class{xcmsSet}} object with the results of the feature
##' detection.
##'
##' @rdname featureDetection-matchedFilter
setMethod("detectFeatures",
          signature(object = "OnDiskMSnExp", param = "MatchedFilterParam"),
          function(object, param, BPPARAM = bpparam(), return.type = "list") {
              return.type <- match.arg(return.type, c("list", "xcmsSet"))
              ## Restrict to MS1 data.
              object <- filterMsLevel(object, msLevel. = 1)
              ## (1) split the object per file.
              ## (2) use bplapply to do the feature detection.
              resList <- bplapply(lapply(1:length(fileNames(object)),
                                     filterFile, object = object),
                              FUN = detectFeatures_OnDiskMSnExp,
                              method = "matchedFilter", param = param)
              ## (3) collect the results.
              res <- .processResultList(resList,
                                        getProcHist = return.type != "list",
                                        fnames = fileNames(object))
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

##' @title Peak detection in the chromatographic time domain
##'
##' @description The \code{detectFeatures,MSnExp,MatchedFilterParam} method
##' performs feature detection using the \emph{matchedFilter} method on all
##' samples from an \code{\link[MSnbase]{MSnExp}} object. These objects contain
##' mz and intensity values of all spectra hence no additional
##' data input from the original files is required.
##'
##' @rdname featureDetection-matchedFilter
setMethod("detectFeatures",
          signature(object = "MSnExp", param = "MatchedFilterParam"),
          function(object, param, BPPARAM = bpparam(), return.type = "list") {
              return.type <- match.arg(return.type, c("list", "xcmsSet"))
              ms1_idx <- which(unname(msLevel(object)) == 1)
              if (length(ms1_idx) == 0)
                  stop("No MS1 spectra available for feature detection!")
              spect_list <- split(spectra(object)[ms1_idx],
                                  fromFile(object)[ms1_idx])
              resList <- bplapply(spect_list, function(z) {
                  detectFeatures_Spectrum_list(z,
                                               method = "matchedFilter",
                                               param = param)
              }, BPPARAM = BPPARAM)
              res <- .processResultList(resList,
                                        getProcHist = return.type != "list",
                                        fnames = fileNames(object))
              if (return.type == "list")
                  return(res$peaks)
              if (return.type == "xcmsSet") {
                  xs <- .pSet2xcmsSet(object)
                  peaks(xs) <- do.call(rbind, res$peaks)
                  xs@.processHistory <- res$procHist
                  OK <- .validProcessHistory(xs)
                  if (!is.logical(OK))
                      stop(OK)
                  if (!any(colnames(pData(object)) == "class"))
                      message("Note: you might want to set/adjust the",
                              " 'sampclass' of the returned xcmSet object",
                              " before proceeding with the analysis.")
                  return(xs)
              }
          })

## massifquant
## The massifquant feature detection method for OnDiskMSnExp:
##' @title Feature detection using the massifquant method
##'
##' @description The \code{detectFeatures,OnDiskMSnExp,MassifquantParam}
##' method performs feature detection using the \emph{massifquant} algorithm
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
##' @param object For \code{detectFeatures}: Either an
##' \code{\link[MSnbase]{OnDiskMSnExp}} or a \code{\link[MSnbase]{MSnExp}}
##' object containing the MS- and all other experiment-relevant data.
##'
##' For all other methods: a parameter object.
##'
##' @param param An \code{MassifquantParam} object containing all settings for
##' the matchedFilter algorithm.
##'
##' @inheritParams featureDetection-centWave
##'
##' @return For \code{detectFeatures}: if \code{return.type = "list"} a list of
##' length equal to the number of samples with matrices specifying the
##' identified features/peaks. If \code{return.type = "xcmsSet"} an
##' \code{\linkS4class{xcmsSet}} object with the results of the feature
##' detection.
##'
##' @rdname featureDetection-massifquant
setMethod("detectFeatures",
          signature(object = "OnDiskMSnExp", param = "MassifquantParam"),
          function(object, param, BPPARAM = bpparam(), return.type = "list") {
              return.type <- match.arg(return.type, c("list", "xcmsSet"))
              ## Restrict to MS1 data.
              object <- filterMsLevel(object, msLevel. = 1)
              ## (1) split the object per file.
              ## (2) use bplapply to do the feature detection.
              resList <- bplapply(lapply(1:length(fileNames(object)),
                                     filterFile, object = object),
                              FUN = detectFeatures_OnDiskMSnExp,
                              method = "massifquant", param = param)
              ## (3) collect the results.
              res <- .processResultList(resList,
                                        getProcHist = return.type != "list",
                                        fnames = fileNames(object))
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


##' @title Feature detection using the massifquant method
##'
##' @description The \code{detectFeatures,MSnExp,MassifquantParam} method
##' performs feature detection using the \emph{massifquant} method on all
##' samples from an \code{\link[MSnbase]{MSnExp}} object. These objects contain
##' mz and intensity values of all spectra hence no additional
##' data input from the original files is required.
##'
##' @rdname featureDetection-massifquant
setMethod("detectFeatures",
          signature(object = "MSnExp", param = "MassifquantParam"),
          function(object, param, BPPARAM = bpparam(), return.type = "list") {
              return.type <- match.arg(return.type, c("list", "xcmsSet"))
              ms1_idx <- which(unname(msLevel(object)) == 1)
              if (length(ms1_idx) == 0)
                  stop("No MS1 spectra available for feature detection!")
              spect_list <- split(spectra(object)[ms1_idx],
                                  fromFile(object)[ms1_idx])
              resList <- bplapply(spect_list, function(z) {
                  detectFeatures_Spectrum_list(z,
                                               method = "massifquant",
                                               param = param)
              }, BPPARAM = BPPARAM)
              res <- .processResultList(resList,
                                        getProcHist = return.type != "list",
                                        fnames = fileNames(object))
              if (return.type == "list")
                  return(res$peaks)
              if (return.type == "xcmsSet") {
                  xs <- .pSet2xcmsSet(object)
                  peaks(xs) <- do.call(rbind, res$peaks)
                  xs@.processHistory <- res$procHist
                  OK <- .validProcessHistory(xs)
                  if (!is.logical(OK))
                      stop(OK)
                  if (!any(colnames(pData(object)) == "class"))
                      message("Note: you might want to set/adjust the",
                              " 'sampclass' of the returned xcmSet object",
                              " before proceeding with the analysis.")
                  return(xs)
              }
          })
