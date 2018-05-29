## Methods for the XCMSnExp object representing untargeted metabolomics
## results
#' @include functions-XCMSnExp.R do_groupChromPeaks-functions.R functions-utils.R
#' do_adjustRtime-functions.R methods-xcmsRaw.R functions-OnDiskMSnExp.R

setMethod("initialize", "XCMSnExp", function(.Object, ...) {
    classVersion(.Object)["XCMSnExp"] <- "0.0.1"
    .Object <- callNextMethod(.Object, ...)
    lockEnvironment(.Object@msFeatureData)
    return(.Object)
})

#' @aliases show,MsFeatureData-method
#' 
#' @rdname XCMSnExp-class
setMethod("show", "XCMSnExp", function(object) {
    callNextMethod()
    ## And not XCMSnExp related stuff.
    cat("- - - xcms preprocessing - - -\n")
    if (hasChromPeaks(object)) {
        cat("Chromatographic peak detection:\n")
        ph <- processHistory(object, type = .PROCSTEP.PEAK.DETECTION)
        cat(" method:", .param2string(ph[[1]]@param), "\n")
        cat(" ", nrow(chromPeaks(object)), " peaks identified in ",
            length(fileNames(object)), " samples.\n", sep = "")
        cat(" On average ",
            format(mean(table(chromPeaks(object)[, "sample"])), digits = 3),
            " chromatographic peaks per sample.\n", sep = "")
    }
    if (hasAdjustedRtime(object)) {
        cat("Alignment/retention time adjustment:\n")
        ph <- processHistory(object, type = .PROCSTEP.RTIME.CORRECTION)
        cat(" method:", .param2string(ph[[1]]@param), "\n")
    }
    if (hasFeatures(object)) {
        cat("Correspondence:\n")
        ph <- processHistory(object, type = .PROCSTEP.PEAK.GROUPING)
        cat(" method:", .param2string(ph[[1]]@param), "\n")
        cat(" ", nrow(featureDefinitions(object)), " features identified.\n",
            sep = "")
        cat(" Median mz range of features: ",
            format(median(featureDefinitions(object)[, "mzmax"] -
                          featureDefinitions(object)[, "mzmin"]), digits = 5),
            "\n", sep = "")
        cat(" Median rt range of features: ",
            format(median(featureDefinitions(object)[, "rtmax"] -
                          featureDefinitions(object)[, "rtmin"]), digits = 5),
            "\n", sep = "")
        if (.hasFilledPeaks(object)) {
            totF <- chromPeaks(object)[, "is_filled"] == 1
            fp <- chromPeaks(object)[totF, , drop = FALSE]
            cat("", sum(totF), "filled peaks (on average",
                mean(table(fp[, "sample"])), "per sample).\n")
        }
    }
})

#' @aliases hasAdjustedRtime hasAdjustedRtime,MsFeatureData-method hasAdjustedRtime,OnDiskMSnExp-method
#'
#' @description
#'
#' \code{hasAdjustedRtime}: whether the object provides adjusted
#' retention times.
#'
#' @rdname XCMSnExp-class
setMethod("hasAdjustedRtime", "XCMSnExp", function(object) {
    hasAdjustedRtime(object@msFeatureData)
})

#' @aliases hasFeatures hasFeatures,MsFeatureData-method
#'
#' @description
#'
#' \code{hasFeatures}: whether the object contains correspondence
#' results (i.e. features).
#'
#' @rdname XCMSnExp-class
setMethod("hasFeatures", "XCMSnExp", function(object) {
    hasFeatures(object@msFeatureData)
})

#' @aliases hasChromPeaks hasChromPeaks,MsFeatureData-method
#'
#' @description
#'
#' \code{hasChromPeaks}: whether the object contains peak
#' detection results.
#'
#' @rdname XCMSnExp-class
setMethod("hasChromPeaks", "XCMSnExp", function(object) {
    hasChromPeaks(object@msFeatureData)
})

#' @aliases adjustedRtime adjustedRtime,MsFeatureData-method
#'
#' @description
#'
#' \code{adjustedRtime},\code{adjustedRtime<-}:
#' extract/set adjusted retention times. \code{adjustedRtime<-} should not
#' be called manually, it is called internally by the
#' \code{\link{adjustRtime}} methods. For \code{XCMSnExp} objects,
#' \code{adjustedRtime<-} does also apply retention time adjustments to
#' eventually present chromatographic peaks. The \code{bySample} parameter
#' allows to specify whether the adjusted retention time should be grouped
#' by sample (file).
#'
#' @return
#'
#' For \code{adjustedRtime}: if \code{bySample = FALSE} a \code{numeric}
#' vector with the adjusted retention for each spectrum of all files/samples
#' within the object. If \code{bySample = TRUE } a \code{list} (length equal
#' to the number of samples) with adjusted retention times grouped by
#' sample. Returns \code{NULL} if no adjusted retention times are present.
#'
#' @rdname XCMSnExp-class
setMethod("adjustedRtime", "XCMSnExp", function(object, bySample = FALSE) {
    res <- adjustedRtime(object@msFeatureData)
    if (length(res) == 0)
        return(res)
    ## Adjusted retention time is a list of retention times.
    if (!bySample) {
        ## Have to re-order the adjusted retention times by spectrum name, such
        ## that rtime are.
        res <- unlist(res, use.names = FALSE)
        sNames <- unlist(split(featureNames(object), fromFile(object)),
                         use.names = FALSE)
        names(res) <- sNames
        res <- res[featureNames(object)]
    }
    return(res)
})
#' @aliases adjustedRtime<- adjustedRtime<-,MsFeatureData-method
#'
#' @rdname XCMSnExp-class
setReplaceMethod("adjustedRtime", "XCMSnExp", function(object, value) {
    if (!is.list(value))
        stop("'value' is supposed to be a list of retention time values!")
    if (hasAdjustedRtime(object))
        object <- dropAdjustedRtime(object)
    ## Check if we have some unsorted retention times (issue #146)
    unsorted <- unlist(lapply(value, is.unsorted), use.names = FALSE)
    if (any(unsorted))
        warning("Adjusted retention times for file(s) ",
                paste(basename(fileNames(object)[unsorted]), collapse = ", "),
                " not sorted increasingly.")
    newFd <- new("MsFeatureData")
    newFd@.xData <- .copy_env(object@msFeatureData)
    adjustedRtime(newFd) <- value
    if (hasChromPeaks(newFd)) {
        ## Change also the retention times reported in the peak matrix.
        if (length(value) != length(rtime(object, bySample = TRUE)))
            stop("The length of 'value' has to match the number of samples!")
        message("Applying retention time adjustment to the identified",
                " chromatographic peaks ... ", appendLF = FALSE)
        fts <- .applyRtAdjToChromPeaks(chromPeaks(newFd),
                                       rtraw = rtime(object, bySample = TRUE),
                                       rtadj = value)
        ## Calling this on the MsFeatureData to avoid all results being removed
        ## again by the chromPeaks<- method.
        chromPeaks(newFd) <- fts
        message("OK")
    }
    lockEnvironment(newFd, bindings = TRUE)
    object@msFeatureData <- newFd
    if (validObject(object))
        return(object)
})

#' @aliases featureDefinitions featureDefinitions,MsFeatureData-method
#'
#' @description
#'
#' \code{featureDefinitions}, \code{featureDefinitions<-}: extract
#' or set the correspondence results, i.e. the mz-rt features (peak groups).
#' Similar to the \code{chromPeaks} it is possible to extract features for
#' specified m/z and/or rt ranges. The function supports also the parameter
#' \code{type} that allows to specify which features to be returned if any
#' of \code{rt} or \code{mz} is specified. For details see help of
#' \code{chromPeaks}.
#' See also \code{\link{featureSummary}} for a function to calculate simple
#' feature summaries.
#' 
#' @return
#'
#' For \code{featureDefinitions}: a \code{DataFrame} with peak grouping
#' information, each row corresponding to one mz-rt feature (grouped peaks
#' within and across samples) and columns \code{"mzmed"} (median mz value),
#' \code{"mzmin"} (minimal mz value), \code{"mzmax"} (maximum mz value),
#' \code{"rtmed"} (median retention time), \code{"rtmin"} (minimal retention
#' time), \code{"rtmax"} (maximal retention time) and \code{"peakidx"}.
#' Column \code{"peakidx"} contains a \code{list} with indices of
#' chromatographic peaks (rows) in the matrix returned by the
#' \code{chromPeaks} method that belong to that feature group. The method
#' returns \code{NULL} if no feature definitions are present.
#'
#' @rdname XCMSnExp-class
setMethod("featureDefinitions", "XCMSnExp", function(object, mz = numeric(),
                                                     rt = numeric(), ppm = 0,
                                                     type = c("any", "within",
                                                              "apex_within")) {
    feat_def <- featureDefinitions(object@msFeatureData)
    type <- match.arg(type)
    ## Select features within rt range.
    if (length(rt) && nrow(feat_def)) {
        rt <- range(rt)
        if (type == "any")
            keep <- which(feat_def$rtmin <= rt[2] & feat_def$rtmax >= rt[1])
        if (type == "within")
            keep <- which(feat_def$rtmin >= rt[1] & feat_def$rtmax <= rt[2])
        if (type == "apex_within")
            keep <- which(feat_def$rtmed >= rt[1] & feat_def$rtmed <= rt[2])
        feat_def <- feat_def[keep, , drop = FALSE]
    }
    ## Select peaks within mz range, considering also ppm
    if (length(mz) && nrow(feat_def)) {
        mz <- range(mz)
        ## Increase mz by ppm.
        if (is.finite(mz[1]))
            mz[1] <- mz[1] - mz[1] * ppm / 1e6
        if (is.finite(mz[2]))
            mz[2] <- mz[2] + mz[2] * ppm / 1e6
        if (type == "any")
            keep <- which(feat_def$mzmin <= mz[2] & feat_def$mzmax >= mz[1])
        if (type == "within")
            keep <- which(feat_def$mzmin >= mz[1] & feat_def$mzmax <= mz[2])
        if (type == "apex_within")
            keep <- which(feat_def$mzmed >= mz[1] & feat_def$mzmed <= mz[2])
        feat_def <- feat_def[keep, , drop = FALSE]
    }
    feat_def
})
#' @aliases featureDefinitions<- featureDefinitions<-,MsFeatureData-method
#'
#' @rdname XCMSnExp-class
setReplaceMethod("featureDefinitions", "XCMSnExp", function(object, value) {
    if (hasFeatures(object))
        object <- dropFeatureDefinitions(object)
    newFd <- new("MsFeatureData")
    newFd@.xData <- .copy_env(object@msFeatureData)
    featureDefinitions(newFd) <- value
    lockEnvironment(newFd, bindings = TRUE)
    object@msFeatureData <- newFd
    if (validObject(object)) {
        ## Lock the environment so that only accessor methods can change values.
        ## lockEnvironment(newFd, bindings = TRUE)
        ## object@msFeatureData <- newFd
        return(object)
    }
})

#' @aliases chromPeaks chromPeaks,MsFeatureData-method
#'
#' @description
#'
#' \code{chromPeaks}, \code{chromPeaks<-}: extract or set
#' the matrix containing the information on identified chromatographic
#' peaks. Parameter \code{bySample} allows to specify whether peaks should
#' be returned ungrouped (default \code{bySample = FALSE}) or grouped by
#' sample (\code{bySample = TRUE}). The \code{chromPeaks<-} method for
#' \code{XCMSnExp} objects removes also all correspondence (peak grouping)
#' and retention time correction (alignment) results. The optional
#' arguments \code{rt}, \code{mz}, \code{ppm} and \code{type} allow to extract
#' only chromatographic peaks overlapping the defined retention time and/or
#' m/z ranges. Argument \code{type} allows to define how \emph{overlapping} is
#' determined: for \code{type == "any"} (the default), all peaks that are even
#' partially overlapping the region are returned, for \code{type == "within"}
#' the full peak has to be within the region and for
#' \code{type == "apex_within"} the peak's apex position (highest signal of the
#' peak) has to be within the region.
#' See description of the return value for details on the returned matrix.
#' Users usually don't have to use the \code{chromPeaks<-} method directly
#' as detected chromatographic peaks are added to the object by the
#' \code{\link{findChromPeaks}} method.
#'
#' @param rt optional \code{numeric(2)} defining the retention time range for
#'     which chromatographic peaks should be returned.
#'
#' @param mz optional \code{numeric(2)} defining the mz range for which
#'     chromatographic peaks should be returned.
#'
#' @param ppm optional \code{numeric(1)} specifying the ppm by which the
#'     \code{mz} range should be extended. For a value of \code{ppm = 10}, all
#'     peaks within \code{mz[1] - ppm / 1e6} and \code{mz[2] + ppm / 1e6} are
#'     returned.
#' 
#' @return
#'
#' For \code{chromPeaks}: if \code{bySample = FALSE} a \code{matrix}
#' with at least the following columns:
#' \code{"mz"} (intensity-weighted mean of mz values of the peak across
#' scans/retention times),
#' \code{"mzmin"} (minimal mz value),
#' \code{"mzmax"} (maximal mz value),
#' \code{"rt"} (retention time for the peak apex),
#' \code{"rtmin"} (minimal retention time),
#' \code{"rtmax"} (maximal retention time),
#' \code{"into"} (integrated, original, intensity of the peak),
#' \code{"maxo"} (maximum intentity of the peak),
#' \code{"sample"} (sample index in which the peak was identified) and
#' \code{"is_filled"} defining whether the chromatographic peak was
#' identified by the peak picking algorithm (\code{0}) or was added by the
#' \code{fillChromPeaks} method (\code{1}).
#' Depending on the employed peak detection algorithm and the
#' \code{verboseColumns} parameter of it additional columns might be
#' returned. For \code{bySample = TRUE} the chronatographic peaks are
#' returned as a \code{list} of matrices, each containing the
#' chromatographic peaks of a specific sample. For samples in which no
#' peaks were detected a matrix with 0 rows is returned.
#'
#' @rdname XCMSnExp-class
setMethod("chromPeaks", "XCMSnExp", function(object, bySample = FALSE,
                                             rt = numeric(), mz = numeric(),
                                             ppm = 0,
                                             type = c("any", "within",
                                                      "apex_within")) {
    pks <- chromPeaks(object@msFeatureData)
    type <- match.arg(type)
    ## Select peaks within rt range.
    if (length(rt)) {
        rt <- range(rt)
        if (type == "any")
            keep <- which(pks[, "rtmin"] <= rt[2] & pks[, "rtmax"] >= rt[1])
        if (type == "within")
            keep <- which(pks[, "rtmin"] >= rt[1] & pks[, "rtmax"] <= rt[2])
        if (type == "apex_within")
            keep <- which(pks[, "rt"] >= rt[1] & pks[, "rt"] <= rt[2])
        pks <- pks[keep, , drop = FALSE]
    }
    ## Select peaks within mz range, considering also ppm
    if (length(mz) && length(pks)) {
        mz <- range(mz)
        ## Increase mz by ppm.
        if (is.finite(mz[1]))
            mz[1] <- mz[1] - mz[1] * ppm / 1e6
        if (is.finite(mz[2]))
            mz[2] <- mz[2] + mz[2] * ppm / 1e6
        if (type == "any")
            keep <- which(pks[, "mzmin"] <= mz[2] & pks[, "mzmax"] >= mz[1])
        if (type == "within")
            keep <- which(pks[, "mzmin"] >= mz[1] & pks[, "mzmax"] <= mz[2])
        if (type == "apex_within")
            keep <- which(pks[, "mz"] >= mz[1] & pks[, "mz"] <= mz[2])
        pks <- pks[keep, , drop = FALSE]
    }
    if (bySample) {
        ## Ensure we return something for each sample in case there is a sample
        ## without detected peaks.
        res <- vector("list", length(fileNames(object)))
        names(res) <- as.character(1:length(res))
        if (length(pks)) {
            tmp <- split.data.frame(pks,
                                    f = pks[, "sample"])
            res[as.numeric(names(tmp))] <- tmp
            if (any(lengths(res) == 0)) {
                emat <- matrix(nrow = 0, ncol = ncol(tmp[[1]]))
                colnames(emat) <- colnames(tmp[[1]])
                for (i in which(lengths(res) == 0))
                    res[[i]] <- emat
            }
        } else {
            for(i in 1:length(res))
                res[[i]] <- pks
        }
        res
    } else
        pks
})
#' @aliases chromPeaks<- chromPeaks<-,MsFeatureData-method
#'
#' @rdname XCMSnExp-class
setReplaceMethod("chromPeaks", "XCMSnExp", function(object, value) {
    newFd <- new("MsFeatureData")
    ## Dropping all alignment results and all retention time corrections.
    suppressMessages(
        object <- dropChromPeaks(object)
    )
    ## Ensure that we remove ALL related process history steps
    newFd@.xData <- .copy_env(object@msFeatureData)
    chromPeaks(newFd) <- value
    lockEnvironment(newFd, bindings = TRUE)
    object@msFeatureData <- newFd
    if (validObject(object)) {
        return(object)
    }
})

#' @description
#'
#' \code{rtime}: extracts the retention time for each
#' scan. The \code{bySample} parameter allows to return the values grouped
#' by sample/file and \code{adjusted} whether adjusted or raw retention
#' times should be returned. By default the method returns adjusted
#' retention times, if they are available (i.e. if retention times were
#' adjusted using the \code{\link{adjustRtime}} method).
#'
#' @param bySample logical(1) specifying whether results should be grouped by
#'     sample.
#'
#' @param adjusted logical(1) whether adjusted or raw (i.e. the original
#'     retention times reported in the files) should be returned.
#' 
#' @return
#'
#' For \code{rtime}: if \code{bySample = FALSE} a numeric vector with
#' the retention times of each scan, if \code{bySample = TRUE} a
#' \code{list} of numeric vectors with the retention times per sample.
#'
#' @rdname XCMSnExp-class
setMethod("rtime", "XCMSnExp", function(object, bySample = FALSE,
                                        adjusted = hasAdjustedRtime(object)) {
    if (adjusted) {
        ## ensure that we DO have adjusted retention times.
        if (hasAdjustedRtime(object)) {
            return(adjustedRtime(object = object, bySample = bySample))
        } else {
            warning("Adjusted retention times requested but none present. ",
                    "returning raw retention times instead.")
        }
    }
    ## Alternative:
    ## theM <- getMethod("rtime", "OnDiskMSnExp")
    ## res <- theM(object)
    res <- callNextMethod(object = object)
    if (bySample) {
        tmp <- split(res, fromFile(object))
        res <- vector("list", length(fileNames(object)))
        names(res) <- as.character(1:length(res))
        res[as.numeric(names(tmp))] <- tmp
    }
    return(res)
})

#' @description
#'
#' \code{mz}: extracts the mz values from each scan of
#' all files within an \code{XCMSnExp} object. These values are extracted
#' from the original data files and eventual processing steps are applied
#' \emph{on the fly}. Using the \code{bySample} parameter it is possible to
#' switch from the default grouping of mz values by spectrum/scan to a
#' grouping by sample/file.
#'
#' @return
#'
#' For \code{mz}: if \code{bySample = FALSE} a \code{list} with the mz
#' values (numeric vectors) of each scan. If \code{bySample = TRUE} a
#' \code{list} with the mz values per sample.
#'
#' @rdname XCMSnExp-class
setMethod("mz", "XCMSnExp", function(object, bySample = FALSE,
                                     BPPARAM = bpparam()) {
    res <- callNextMethod(object = object, BPPARAM = BPPARAM)
    if (bySample) {
        tmp <- lapply(split(res, fromFile(object)), unlist, use.names = FALSE)
        res <- vector("list", length(fileNames(object)))
        names(res) <- as.character(1:length(res))
        res[as.numeric(names(tmp))] <- tmp
    }
    return(res)
})

#' @description
#'
#' \code{intensity}: extracts the intensity values from
#' each scan of all files within an \code{XCMSnExp} object. These values are
#' extracted from the original data files and eventual processing steps are
#' applied \emph{on the fly}. Using the \code{bySample} parameter it is
#' possible to switch from the default grouping of intensity values by
#' spectrum/scan to a grouping by sample/file.
#'
#' @return
#'
#' For \code{intensity}: if \code{bySample = FALSE} a \code{list} with
#' the intensity values (numeric vectors) of each scan. If
#' \code{bySample = TRUE} a \code{list} with the intensity values per
#' sample.
#'
#' @rdname XCMSnExp-class
setMethod("intensity", "XCMSnExp", function(object, bySample = FALSE,
                                            BPPARAM = bpparam()) {
    res <- callNextMethod(object = object, BPPARAM = BPPARAM)
    if (bySample) {
        tmp <- lapply(split(res, fromFile(object)), unlist, use.names = FALSE)
        res <- vector("list", length(fileNames(object)))
        names(res) <- as.character(1:length(res))
        res[as.numeric(names(tmp))] <- tmp
    }
    return(res)
})

#' @description
#'
#' \code{spectra}: extracts the
#' \code{\link{Spectrum}} objects containing all data from
#' \code{object}. The values are extracted from the original data files and
#' eventual processing steps are applied \emph{on the fly}. By setting
#' \code{bySample = TRUE}, the spectra are returned grouped by sample/file.
#' If the \code{XCMSnExp} object contains adjusted retention times, these
#' are returned by default in the \code{Spectrum} objects (can be
#' overwritten by setting \code{adjusted = FALSE}).
#'
#' @param BPPARAM Parameter class for parallel processing. See
#'     \code{\link{bpparam}}.
#' 
#' @return
#'
#' For \code{spectra}: if \code{bySample = FALSE} a \code{list} with
#' \code{\link{Spectrum}} objects. If \code{bySample = TRUE} the
#' result is grouped by sample, i.e. as a \code{list} of \code{lists}, each
#' element in the \emph{outer} \code{list} being the \code{list} of spectra
#' of the specific file.
#'
#' @rdname XCMSnExp-class
setMethod("spectra", "XCMSnExp", function(object, bySample = FALSE,
                                          adjusted = hasAdjustedRtime(object),
                                          BPPARAM = bpparam()) {
    if (adjusted & hasAdjustedRtime(object))
        fData(object)$retentionTime <- rtime(object, adjusted = TRUE)
    else object <- as(object, "OnDiskMSnExp")
    res <- callNextMethod(object = object, BPPARAM = BPPARAM)
    if (bySample) {
        tmp <- split(res, fromFile(object))
        ## That's to ensure that we're always returning something for all files.
        res <- vector("list", length(fileNames(object)))
        names(res) <- as.character(1:length(res))
        res[as.numeric(names(tmp))] <- tmp
    }
    return(res)
})

#' @aliases processHistory
#' 
#' @description
#'
#' \code{processHistory}: returns a \code{list} with
#' \code{\link{ProcessHistory}} objects (or objects inheriting from this
#' base class) representing the individual processing steps that have been
#' performed, eventually along with their settings (\code{Param} parameter
#' class). Optional arguments \code{fileIndex}, \code{type} and
#' \code{msLevel} allow to restrict to process steps of a certain type or
#' performed on a certain file or MS level.
#'
#' @param fileIndex For \code{processHistory}: optional \code{numeric}
#'     specifying the index of the files/samples for which the
#'     \code{\link{ProcessHistory}} objects should be retrieved.
#'
#' @param type For \code{processHistory}: restrict returned
#'     \code{\link{ProcessHistory}} objects to analysis steps of a certain
#'     type. Use the \code{processHistoryTypes} to list all supported values.
#'     For \code{chromPeaks}: \code{character} specifying which peaks to return
#'     if \code{rt} or \code{mz} are defined. For \code{type = "any"} all
#'     chromatographic peaks partially overlapping the range defined by 
#'     \code{mz} and/or \code{rt} are returned, \code{type = "within"} returns
#'     only peaks completely within the region and \code{type = "apex_within"}
#'     peaks for which the peak's apex is within the region.
#' 
#' @return
#'
#' For \code{processHistory}: a \code{list} of
#' \code{\link{ProcessHistory}} objects providing the details of the
#' individual data processing steps that have been performed.
#'
#' @rdname XCMSnExp-class
setMethod("processHistory", "XCMSnExp", function(object, fileIndex, type,
                                                 msLevel) {
    ph <- object@.processHistory
    if (length(ph)) {
        if (!missing(fileIndex)) {
            if (!all(fileIndex %in% 1:length(fileNames(object))))
                stop("'fileIndex' has to be within 1 and the number of samples!")
            gotIt <- unlist(lapply(ph, function(z) {
                any(z@fileIndex %in% fileIndex)
            }))
            ph <- ph[gotIt]
        }
        if (!missing(type) & length(ph)) {
            gotIt <- unlist(lapply(ph, function(z) {
                any(type == processType(z))
            }))
            ph <- ph[gotIt]
        }
        if (!missing(msLevel) & length(ph)) {
            gotIt <- unlist(lapply(ph, function(z) {
                msLevel(z) %in% msLevel
            }))
            ph <- ph[gotIt]
        }
        return(ph)
    }
    list()
})

#' @description
#'
#' \code{addProcessHistory}: adds (appends) a single
#' \code{\link{ProcessHistory}} object to the \code{.processHistory} slot.
#'
#' @return
#'
#' The \code{addProcessHistory} method returns the input object with the
#' provided \code{\link{ProcessHistory}} appended to the process history.
#' 
#' @noRd
setMethod("addProcessHistory", "XCMSnExp", function(object, ph) {
    if (!inherits(ph, "ProcessHistory"))
        stop("Argument 'ph' has to be of type 'ProcessHistory' or a class ",
             "extending it!")
    object@.processHistory[[(length(object@.processHistory) + 1)]] <- ph
    if (validObject(object))
        return(object)
})

#' @aliases dropChromPeaks dropChromPeaks,MsFeatureData-method
#'
#' @description
#'
#' \code{dropChromPeaks}: drops any identified chromatographic
#' peaks and returns the object without that information. Note that for
#' \code{XCMSnExp} objects the method drops by default also results from a
#' correspondence (peak grouping) analysis. Adjusted retention times are
#' removed if the alignment has been performed \emph{after} peak detection.
#' This can be overruled with \code{keepAdjustedRtime = TRUE}.
#'
#' @rdname XCMSnExp-class
setMethod("dropChromPeaks", "XCMSnExp", function(object,
                                                 keepAdjustedRtime = FALSE) {
    if (hasChromPeaks(object)) {
        phTypes <- unlist(lapply(processHistory(object), function(z)
            processType(z)))
        object <- dropProcessHistories(object, type = .PROCSTEP.PEAK.DETECTION)
        ## Make sure we delete all related process history steps
        object <- dropProcessHistories(object, type = .PROCSTEP.PEAK.GROUPING)
        object <- dropProcessHistories(object, type = .PROCSTEP.PEAK.FILLING)
        object <- dropProcessHistories(object, type = .PROCSTEP.CALIBRATION)
        newFd <- new("MsFeatureData")
        newFd@.xData <- .copy_env(object@msFeatureData)
        newFd <- dropChromPeaks(newFd)
        if (hasFeatures(newFd))
            newFd <- dropFeatureDefinitions(newFd)
        ## Drop adjusted retention times if performed AFTER peak detection.
        if (hasAdjustedRtime(newFd) & !keepAdjustedRtime) {
            idx_rt_adj <- which(phTypes == .PROCSTEP.RTIME.CORRECTION)
            idx_pk_det <- which(phTypes == .PROCSTEP.PEAK.DETECTION)
            if (idx_rt_adj > idx_pk_det) {
                object <- dropProcessHistories(object,
                                               type = .PROCSTEP.RTIME.CORRECTION)
                newFd <- dropAdjustedRtime(newFd)
            }
        }
        lockEnvironment(newFd, bindings = TRUE)
        object@msFeatureData <- newFd
    }
    if (validObject(object))
        return(object)
})

#' @aliases dropFeatureDefinitions dropFeatureDefinitions,MsFeatureData-method
#'
#' @description
#' 
#' \code{dropFeatureDefinitions}: drops the results from a
#' correspondence (peak grouping) analysis, i.e. the definition of the mz-rt
#' features and returns the object without that information. Note that for
#' \code{XCMSnExp} objects the method will also by default drop retention
#' time adjustment results, if these were performed after the last peak
#' grouping (i.e. which base on the results from the peak grouping that are
#' going to be removed). All related process history steps are
#' removed too as well as eventually filled in peaks
#' (by \code{\link{fillChromPeaks}}). The parameter \code{keepAdjustedRtime}
#' can be used to avoid removal of adjusted retention times.
#'
#' @param keepAdjustedRtime For \code{dropFeatureDefinitions,XCMSnExp}:
#'     \code{logical(1)} defining whether eventually present retention time
#'     adjustment should not be dropped. By default dropping feature definitions
#'     drops retention time adjustment results too.
#'
#' @param dropLastN For \code{dropFeatureDefinitions,XCMSnExp}:
#'     \code{numeric(1)} defining the number of peak grouping related process
#'     history steps to remove. By default \code{dropLastN = -1}, dropping the
#'     chromatographic peaks removes all process history steps related to peak
#'     grouping. Setting e.g. \code{dropLastN = 1} will only remove the most
#'     recent peak grouping related process history step.
#' 
#' @rdname XCMSnExp-class
setMethod("dropFeatureDefinitions", "XCMSnExp", function(object,
                                                         keepAdjustedRtime = FALSE,
                                                         dropLastN = -1) {
    if (hasFeatures(object)) {
        phTypes <- unlist(lapply(processHistory(object), function(z)
            processType(z)))
        idx_art <- which(phTypes == .PROCSTEP.RTIME.CORRECTION)
        idx_fal <- which(phTypes == .PROCSTEP.PEAK.GROUPING)
        if (length(idx_art) == 0)
            idx_art <- -1L
        if (length(idx_fal) == 0)
            idx_fal <- -1L
        ## 1) drop last related process history step and results
        object <- dropProcessHistories(object,
                                       type = .PROCSTEP.PEAK.GROUPING,
                                       num = 1)
        ## Drop also eventual filterFeatureDefinitions
        object <- dropGenericProcessHistory(object,
                                            fun = "filterFeatureDefinitions")
        newFd <- new("MsFeatureData")
        newFd@.xData <- .copy_env(object@msFeatureData)
        newFd <- dropFeatureDefinitions(newFd)
        if (.hasFilledPeaks(object)) {
            ## Remove filled in peaks
            chromPeaks(newFd) <-
                chromPeaks(newFd)[chromPeaks(newFd)[, "is_filled"] == 0, ,
                                  drop = FALSE]
            object <- dropProcessHistories(object, type = .PROCSTEP.PEAK.FILLING)
        }
        lockEnvironment(newFd, bindings = TRUE)
        object@msFeatureData <- newFd
        ## 2) If retention time correction was performed after the latest peak
        ##    alignment, drop also the retention time correction and all related
        ##    process history steps.
        ##    Otherwise (grouping performed after retention time adjustment) do
        ##    nothing - this keeps eventual alignment related process history
        ##    steps performed before retention time correction.
        if (hasAdjustedRtime(object)) {
            if (max(idx_art) > max(idx_fal) & !keepAdjustedRtime) {
                object <- dropProcessHistories(object,
                                               type = .PROCSTEP.PEAK.GROUPING)
                ## This will ensure that the retention times of the peaks
                ## are restored.
                object <- dropAdjustedRtime(object)
                warning("Removed also correspondence (peak grouping) results as",
                        " these based on the retention time correction results",
                        " that were dropped.")
            }
        }
    }
    if (validObject(object))
        return(object)
})

#' @aliases dropAdjustedRtime dropAdjustedRtime,MsFeatureData-method
#'
#' @description
#'
#' \code{dropAdjustedRtime}: drops any retention time
#' adjustment information and returns the object without adjusted retention
#' time. For \code{XCMSnExp} objects, this also reverts the retention times
#' reported for the chromatographic peaks in the peak matrix to the
#' original, raw, ones (after chromatographic peak detection). Note that
#' for \code{XCMSnExp} objects the method drops also all peak grouping
#' results if these were performed \emph{after} the retention time
#' adjustment. All related process history steps are removed too.
#'
#' @rdname XCMSnExp-class
setMethod("dropAdjustedRtime", "XCMSnExp", function(object) {
    if (hasAdjustedRtime(object)) {
        ## Get the process history types to determine the order of the analysis
        ## steps.
        phTypes <- unlist(lapply(processHistory(object), function(z)
            processType(z)))
        idx_art <- which(phTypes == .PROCSTEP.RTIME.CORRECTION)
        idx_fal <- which(phTypes == .PROCSTEP.PEAK.GROUPING)
        if (length(idx_art) == 0)
            idx_art <- -1L
        if (length(idx_fal) == 0)
            idx_fal <- -1L
        ## Copy the content of the object
        newFd <- new("MsFeatureData")
        newFd@.xData <- .copy_env(object@msFeatureData)        
        ## Revert applied adjustments in peaks:
        if (hasChromPeaks(newFd)) {
            message("Reverting retention times of identified peaks to ",
                    "original values ... ", appendLF = FALSE)
            fts <- .applyRtAdjToChromPeaks(chromPeaks(newFd),
                                           rtraw = adjustedRtime(object,
                                                                 bySample = TRUE),
                                           rtadj = rtime(object,
                                                         bySample = TRUE,
                                                         adjusted = FALSE))
            ## Replacing peaks in MsFeatureData, not in XCMSnExp to avoid
            ## all results being removed.
            chromPeaks(newFd) <- fts
            message("OK")
        }
        ## 1) Drop the retention time adjustment and (the latest) related process
        ##    history
        object <- dropProcessHistories(object,
                                       type = .PROCSTEP.RTIME.CORRECTION,
                                       num = 1)
        newFd <- dropAdjustedRtime(newFd)
        object@msFeatureData <- newFd
        lockEnvironment(newFd, bindings = TRUE)
        ## 2) If grouping has been performed AFTER retention time correction it
        ##    has to be dropped too, including ALL related process histories.
        if (hasFeatures(object)) {
            if (max(idx_fal) > max(idx_art)) {
                object <- dropFeatureDefinitions(object)
                object <- dropProcessHistories(object,
                                               type = .PROCSTEP.PEAK.GROUPING,
                                               num = -1)
            }
        } else {
            ## If there is any peak alignment related process history, but no
            ## peak alignment results, drop them.
            object <- dropProcessHistories(object,
                                           type = .PROCSTEP.PEAK.GROUPING,
                                           num = -1)
        }
    }
    if (validObject(object))
        return(object)
})


#' @description
#'
#' The \code{[} method allows to subset a \code{\link{XCMSnExp}}
#' object by spectra. Be aware that the \code{[} method removes all
#' preprocessing results, except adjusted retention times if
#' \code{keepAdjustedRtime = TRUE} is passed to the method.
#'
#' @param x For \code{[} and \code{[[}: an \code{\link{XCMSnExp}} object.
#'
#' @param i For \code{[}: \code{numeric} or \code{logical} vector specifying to
#'     which spectra the data set should be reduced.
#'     For \code{[[}: a single integer or character.
#'
#' @param j For \code{[} and \code{[[}: not supported.
#'
#' @param drop For \code{[} and \code{[[}: not supported.
#'
#' @rdname XCMSnExp-filter-methods
setMethod("[", "XCMSnExp", function(x, i, j, ..., drop = TRUE) {
    if (!missing(j))
        stop("subsetting by columns ('j') not supported")
    if (missing(i))
        return(x)
    else if (!(is.numeric(i) | is.logical(i)))
        stop("'i' has to be either numeric or logical")
    ## Check if we have keepAdjustedRtime as an additional parameter
    ## in ...
    keepAdjustedRtime <- list(...)$ke
    if (is.null(keepAdjustedRtime))
        keepAdjustedRtime <- FALSE
    if (hasFeatures(x) | hasChromPeaks(x)) {
        suppressMessages(
            x <- dropFeatureDefinitions(x, keepAdjustedRtime =
                                               keepAdjustedRtime))
        suppressMessages(
            x <- dropChromPeaks(x, keepAdjustedRtime =
                                       keepAdjustedRtime))
        warning("Removed preprocessing results")
    }
    if (hasAdjustedRtime(x)) {
        if (keepAdjustedRtime) {
            ## Subset the adjusted rtime
            new_adj <- rtime(x, adjusted = TRUE)[i]
            newFd <- new("MsFeatureData")
            newFd@.xData <- .copy_env(x@msFeatureData)        
            adjustedRtime(newFd) <-
                unname(split(new_adj, f = fromFile(x)[i]))
            lockEnvironment(newFd, bindings = TRUE)
            x@msFeatureData <- newFd
        } else {
            suppressMessages(x <- dropAdjustedRtime(x))
        }
    }
    callNextMethod()
})

#' @description
#'
#' \code{[[} extracts a single \code{\link{Spectrum}}
#' object from an \code{XCMSnExp}. The reported retention time is the
#' adjusted retention time if alignment has been performed on \code{x}.
#' 
#' @rdname XCMSnExp-filter-methods
setMethod("[[", "XCMSnExp",
          function(x, i, j, drop = FALSE) {
              ## If it has adjusted retention times, replace raw ones.
              if (hasAdjustedRtime(x))
                  x@featureData$retentionTime <- rtime(x, adjusted = TRUE)
              x <- as(x, "OnDiskMSnExp")
              callNextMethod()
          })

#' @title XCMSnExp data manipulation methods inherited from MSnbase
#'
#' @description
#'
#' The methods listed on this page are \code{\link{XCMSnExp}}
#' methods inherited from its parent, the
#' \code{\link{OnDiskMSnExp}} class from the \code{MSnbase}
#' package, that alter the raw data or are related to data subsetting. Thus
#' calling any of these methods causes all \code{xcms} pre-processing
#' results to be removed from the \code{\link{XCMSnExp}} object to ensure
#' its data integrity.
#'
#' \code{bin}: allows to \emph{bin} spectra. See
#' \code{\link{bin}} documentation in the \code{MSnbase} package for more
#' details and examples.
#'
#' @param x \code{\link{XCMSnExp}} or \code{\link{OnDiskMSnExp}}
#'     object.
#' 
#' @param object \code{\link{XCMSnExp}} or \code{\link{OnDiskMSnExp}}
#'     object.
#'
#' @param binSize \code{numeric(1)} defining the size of a bin (in Dalton).
#'
#' @param msLevel. For \code{bin}, \code{clean}, \code{filterMsLevel},
#'     \code{removePeaks}: \code{numeric(1)} defining the MS level(s)
#'     to which operations should be applied or to which the object should be
#'     subsetted.
#'
#' @return For all methods: a \code{XCMSnExp} object.
#'
#' @rdname XCMSnExp-inherited-methods
#'
#' @seealso \code{\link{XCMSnExp-filter}} for methods to filter and subset
#'     \code{XCMSnExp} objects.
#'     \code{\link{XCMSnExp}} for base class documentation.
#'     \code{\link{OnDiskMSnExp}} for the documentation of the
#'     parent class.
#'
#' @author Johannes Rainer
setMethod("bin", "XCMSnExp", function(object, binSize = 1L, msLevel.) {
    if (hasAdjustedRtime(object) | hasFeatures(object) |
        hasChromPeaks(object)) {
        ## object@.processHistory <- list()
        ## object@msFeatureData <- new("MsFeatureData")
        object <- dropAdjustedRtime(object)
        object <- dropFeatureDefinitions(object)
        object <- dropChromPeaks(object)
        warning("Removed preprocessing results")
    }
    callNextMethod()
})

#' @description
#'
#' \code{clean}: removes unused \code{0} intensity data
#' points. See \code{\link{clean}} documentation in the \code{MSnbase} package
#' for details and examples.
#'
#' @param all For \code{clean}: \code{logical(1)}, if \code{TRUE} all zeros are
#'     removed.
#'
#' @param verbose \code{logical(1)} whether progress information should be
#'     displayed.
#'
#' @rdname XCMSnExp-inherited-methods
setMethod("clean", "XCMSnExp", function(object, all = FALSE,
                                        verbose = FALSE, msLevel.) {
    if (hasAdjustedRtime(object) | hasFeatures(object) |
        hasChromPeaks(object)) {
        ## object@.processHistory <- list()
        ## object@msFeatureData <- new("MsFeatureData")
        object <- dropAdjustedRtime(object)
        object <- dropFeatureDefinitions(object)
        object <- dropChromPeaks(object)
        warning("Removed preprocessing results")
    }
    callNextMethod()
})

#' @description
#'
#' \code{filterMsLevel}: reduces the \code{\link{XCMSnExp}}
#' object to spectra of the specified MS level(s). See
#' \code{\link{filterMsLevel}} documentation for details and
#' examples. Presently, if \code{msLevel.} is provided, the function
#' removes identified chromatographic peaks and correspondence results
#' while keeping adjusted retention times by default (if present). The
#' latter can be changed setting \code{keepAdjustedRtime = FALSE}.
#' 
#' @rdname XCMSnExp-filter-methods
setMethod("filterMsLevel", "XCMSnExp", function(object, msLevel.,
                                                keepAdjustedRtime =
                                                    hasAdjustedRtime(object)) {
    if (missing(msLevel.)) return(object)
    msLevel. <- as.numeric(msLevel.)
    keep_logical <- msLevel(object) %in% msLevel.
    msg <- paste0("Filter: select MS level(s) ",
                  paste(unique(msLevel.), collapse = " "))
    msg <- paste0(msg, " [", date(), "]")
    if (!any(keep_logical)) {
        res <- new("XCMSnExp")
        res@processingData@processing <- c(res@processingData@processing, msg)
        return(res)
    }

    ## In future we might want to keep also the chromatographic peaks of the
    ## correct MS level.
    if (hasChromPeaks(object))
        warning("Identified chromatographic peaks removed")
    if (hasFeatures(object))
        warning("Feature definitions removed")

    ## Create a new empty MsFeatureData and just add adjusted retention times
    newMfd <- new("MsFeatureData")
    ph <- processHistory(object)
    ## 2) Subset adjusted retention time
    ## if (hasAdjustedRtime(object) & length(keep_fts)) {
    if (hasAdjustedRtime(object)) {
        if (keepAdjustedRtime) {
            ## issue #210: keep adjusted retention times if wanted.
            ## CAVE: we're still removing the chromatographic peaks at the
            ## moment thus we might miss later the peaks and groups on which
            ## alignment has been performed (for peak groups method).
            keep_by_file <- base::split(keep_logical, fromFile(object))
            adj_rt <- base::mapply(FUN = function(y, z) {
                return(y[z])
            }, y = adjustedRtime(object, bySample = TRUE), z = keep_by_file,
            SIMPLIFY = FALSE)
            adjustedRtime(newMfd) <- adj_rt
            ph <- dropProcessHistoriesList(ph,
                                           type = c(.PROCSTEP.PEAK.DETECTION,
                                                    .PROCSTEP.PEAK.GROUPING,
                                                    .PROCSTEP.PEAK.FILLING,
                                                    .PROCSTEP.CALIBRATION))
        } else {
            object <- dropAdjustedRtime(object)
            ph <- dropProcessHistoriesList(ph,
                                           type = .PROCSTEP.RTIME.CORRECTION)
        }
    }
    tmp <- as(object, "OnDiskMSnExp")[base::which(keep_logical)]
    object <- as(tmp, "XCMSnExp")
    ## Put the stuff back
    object@processingData@processing <- c(object@processingData@processing, msg)
    lockEnvironment(newMfd, bindings = TRUE)
    object@msFeatureData <- newMfd
    object@.processHistory <- ph
    if (validObject(object))
        object
})

#' @description \code{filterAcquisitionNum}: filters the
#' \code{\link{XCMSnExp}} object keeping only spectra with the provided
#' acquisition numbers. See \code{\link{filterAcquisitionNum}} for
#' details and examples.
#'
#' @param n For \code{filterAcquisitionNum}: \code{integer} defining the
#'     acquisition numbers of the spectra to which the data set should be
#'     sub-setted.
#'
#' @param file For \code{filterAcquisitionNum}: 
#'     \code{integer} defining the file index within the object to subset the
#'     object by file.
#'
#' @rdname XCMSnExp-inherited-methods
setMethod("filterAcquisitionNum", "XCMSnExp", function(object, n, file) {
    if (hasAdjustedRtime(object) | hasFeatures(object) |
        hasChromPeaks(object)) {
        ## object@.processHistory <- list()
        ## object@msFeatureData <- new("MsFeatureData")
        object <- dropAdjustedRtime(object)
        object <- dropFeatureDefinitions(object)
        object <- dropChromPeaks(object)
        warning("Removed preprocessing results")
    }
    callNextMethod()
})

#' @aliases XCMSnExp-filter
#' 
#' @title XCMSnExp filtering and subsetting
#'
#' @description
#'
#' The methods listed on this page allow to filter and subset
#' \code{\link{XCMSnExp}} objects. Most of them are inherited from the
#' \code{\link{OnDiskMSnExp}} object and have been adapted for
#' \code{\link{XCMSnExp}} to enable subsetting also on the preprocessing
#' results.
#'
#' \code{filterFile}: allows to reduce the
#' \code{\link{XCMSnExp}} to data from only certain files. Identified
#' chromatographic peaks for these files are retained while all eventually
#' present features (peak grouping information) are dropped. By default 
#' also adjusted retention times are removed (if present). This can be
#' overwritten by setting \code{keepAdjustedRtime = TRUE}.
#'
#' @details
#'
#' All subsetting methods try to ensure that the returned data is
#' consistent. Correspondence results for example are removed if the data
#' set is sub-setted by file, since the correspondence results are dependent
#' on the files on which correspondence was performed. Thus, some filter
#' and sub-setting methods drop some of the preprocessing results. An
#' exception are the adjusted retention times: most subsetting methods
#' support the argument \code{keepAdjustedRtime} (even the \code{[} method)
#' that forces the adjusted retention times to be retained even if the
#' default would be to drop them.
#' 
#' @note
#'
#' The \code{filterFile} method removes also process history steps not
#' related to the files to which the object should be sub-setted and updates
#' the \code{fileIndex} attribute accordingly. Also, the method does not
#' allow arbitrary ordering of the files or re-ordering of the files within
#' the object.
#'
#' Note also that most of the filtering methods, and also the subsetting
#' operations \code{[} drop all or selected preprocessing results. To
#' consolidate the alignment results, i.e. ensure that adjusted retention
#' times are always preserved, use the \code{\link{applyAdjustedRtime}}
#' function on the object that contains the alignment results. This replaces
#' the raw retention times with the adjusted ones.
#'
#' @param object A \code{\link{XCMSnExp}} object.
#'
#' @param file For \code{filterFile}: \code{integer} defining the file index
#'     within the object to subset the object by file or \code{character}
#'     specifying the file names to sub set. The indices are expected to be
#'     increasingly ordered, if not they are ordered internally.
#'
#' @param keepAdjustedRtime For \code{filterFile}, \code{filterMsLevel},
#'     \code{[} \code{split}: \code{logical(1)} defining whether the adjusted
#'     retention times should be kept, even if e.g. features are being removed
#'     (and the retention time correction was performed on these features).
#' 
#' @return All methods return an \code{\link{XCMSnExp}} object.
#'
#' @author Johannes Rainer
#'
#' @seealso \code{\link{XCMSnExp}} for base class documentation.
#'
#' @rdname XCMSnExp-filter-methods
#' 
#' @examples
#'
#' ## Load some of the files from the faahKO package.
#' library(faahKO)
#' fs <- c(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
#'         system.file('cdf/KO/ko16.CDF', package = "faahKO"),
#'         system.file('cdf/KO/ko18.CDF', package = "faahKO"))
#' ## Read the files
#' od <- readMSData(fs, mode = "onDisk")
#'
#' ## Perform peak detection on them using the matched filter algorithm. Note
#' ## that we use a large value for binSize to reduce the runtime of the
#' ## example code.
#' mfp <- MatchedFilterParam(binSize = 5)
#' xod <- findChromPeaks(od, param = mfp)
#'
#' ## Subset the dataset to the first and third file.
#' xod_sub <- filterFile(xod, file = c(1, 3))
#'
#' ## The number of chromatographic peaks per file for the full object
#' table(chromPeaks(xod)[, "sample"])
#'
#' ## The number of chromatographic peaks per file for the subset
#' table(chromPeaks(xod_sub)[, "sample"])
#'
#' basename(fileNames(xod))
#' basename(fileNames(xod_sub))
#'
#' ## Filter on mz values; chromatographic peaks and features within the
#' ## mz range are retained (as well as adjusted retention times).
#' xod_sub <- filterMz(xod, mz = c(300, 400))
#' head(chromPeaks(xod_sub))
#' nrow(chromPeaks(xod_sub))
#' nrow(chromPeaks(xod))
#'
#' ## Filter on rt values. All chromatographic peaks and features within the
#' ## retention time range are retained. Filtering is performed by default on
#' ## adjusted retention times, if present.
#' xod_sub <- filterRt(xod, rt = c(2700, 2900))
#'
#' range(rtime(xod_sub))
#' head(chromPeaks(xod_sub))
#' range(chromPeaks(xod_sub)[, "rt"])
#'
#' nrow(chromPeaks(xod))
#' nrow(chromPeaks(xod_sub))
#'
#' ## Extract a single Spectrum
#' xod[[4]]
#'
#' ## Subsetting using [ removes all preprocessing results - using
#' ## keepAdjustedRtime = TRUE would keep adjusted retention times, if present.
#' xod_sub <- xod[fromFile(xod) == 1]
#' xod_sub
#'
#' ## Using split does also remove preprocessing results, but it supports the
#' ## optional parameter keepAdjustedRtime.
#' ## Split the object into a list of XCMSnExp objects, one per file
#' xod_list <- split(xod, f = fromFile(xod))
#' xod_list
setMethod("filterFile", "XCMSnExp", function(object, file,
                                             keepAdjustedRtime = FALSE) {
    if (missing(file)) return(object)
    if (is.character(file)) {
        file <- base::match(file, basename(fileNames(object)))
    }
    ## This will not work if we want to get the files in a different
    ## order (i.e. c(3, 1, 2, 5))
    file <- base::sort(unique(file))
    ## Error checking - seems that's not performed downstream.
    if (!all(file %in% 1:length(fileNames(object))))
        stop("'file' has to be within 1 and the number of files in the object!")
    ## Drop features
    has_features <- hasFeatures(object)
    has_chrom_peaks <- hasChromPeaks(object)
    has_adj_rt <- hasAdjustedRtime(object)
    if (has_features) {
        message("Correspondence results (features) removed.")
        object <- dropFeatureDefinitions(object,
                                         keepAdjustedRtime = keepAdjustedRtime)
    }
    if (has_adj_rt & !keepAdjustedRtime){
        object <- dropAdjustedRtime(object)
        has_adj_rt <- FALSE
    }
    ## Extracting all the XCMSnExp data from the object.
    ph <- processHistory(object)
    newFd <- new("MsFeatureData")
    newFd@.xData <- .copy_env(object@msFeatureData)
    ## Subset original data:
    suppressWarnings(
        object <- callNextMethod(object = object, file = file)
    )
    ## Subset the results per file:
    if (has_adj_rt) {
        adjustedRtime(newFd) <- adjustedRtime(newFd)[file]
    }
    if (has_chrom_peaks) {
        pks <- chromPeaks(newFd)
        pks <- pks[pks[, "sample"] %in% file, , drop = FALSE]
        pks[, "sample"] <- match(pks[, "sample"], file)
        chromPeaks(newFd) <- pks
    }
    ## Remove ProcessHistory not related to any of the files.
    if (length(ph)) {
        kp <- unlist(lapply(ph, function(z) {
            any(fileIndex(z) %in% file)
        }))
        ph <- ph[kp]
    }
    ## Update file index in process histories.
    if (length(ph)) {
        ph <- lapply(ph, function(z) {
            updateFileIndex(z, old = file, new = 1:length(file))
        })
    }
    lockEnvironment(newFd, bindings = TRUE)
    object@msFeatureData <- newFd
    object@.processHistory <- ph
    if (validObject(object))
        object
})

#' @description
#'
#' \code{filterMz}: filters the data set based on the
#' provided mz value range. All chromatographic peaks and features (grouped
#' peaks) falling completely within the provided mz value range are retained
#' (if their minimal mz value is \code{>= mz[1]} and the maximal mz value
#' \code{<= mz[2]}. Adjusted retention times, if present, are not altered by
#' the filtering.
#'
#' @param mz For \code{filterMz}: \code{numeric(2)} defining the lower and upper
#'     mz value for the filtering.
#'
#' @param msLevel. For \code{filterMz}, \code{filterRt}, \code{numeric(1)}
#'     defining the MS level(s) to which operations should be applied or to
#'     which the object should be subsetted. See \code{\link{filterMz}}
#'     for more details
#'
#' @param ... Optional additional arguments.
#'
#' @rdname XCMSnExp-filter-methods
setMethod("filterMz", "XCMSnExp", function(object, mz, msLevel., ...) {
    if (missing(mz))
        return(object)
    if (!is.numeric(mz))
        stop("'mz' has to be a numeric vector of length(2)!")
    mz <- range(mz)
    ## Subset peaks if present.
    object <- callNextMethod()  # just adds to processing queue.

    if (hasChromPeaks(object)) {
        fts <- chromPeaks(object)
        keepIdx <- which(fts[, "mzmin"] >= mz[1] & fts[, "mzmax"] <= mz[2])
        newE <- .filterChromPeaks(object@msFeatureData, idx = keepIdx)
        lockEnvironment(newE, bindings = TRUE)
        object@msFeatureData <- newE
    }
    if (validObject(object))
        return(object)
})

#' @description
#'
#' \code{filterRt}: filters the data set based on the
#' provided retention time range. All chromatographic peaks and features
#' (grouped peaks) the specified retention time window are retained (i.e. if
#' the retention time corresponding to the peak's apex is within the
#' specified rt range). If retention time correction has been performed,
#' the method will by default filter the object by adjusted retention times.
#' The argument \code{adjusted} allows to specify manually whether filtering
#' should be performed by raw or adjusted retention times. Filtering by
#' retention time does not drop any preprocessing results nor does it remove
#' or change alignment results (i.e. adjusted retention times).
#' The method returns an empty object if no spectrum or feature is within
#' the specified retention time range.
#'
#' @param rt For \code{filterRt}: \code{numeric(2)} defining the retention time
#'     window (lower and upper bound) for the filtering.
#'
#' @param adjusted For \code{filterRt}: \code{logical} indicating whether the
#'     object should be filtered by original (\code{adjusted = FALSE}) or
#'     adjusted retention times (\code{adjusted = TRUE}).
#'     For \code{spectra}: whether the retention times in the individual
#'     \code{Spectrum} objects should be the adjusted or raw retention times.
#'
#' @rdname XCMSnExp-filter-methods
setMethod("filterRt", "XCMSnExp", function(object, rt, msLevel.,
                                           adjusted = hasAdjustedRtime(object)) {
    if (missing(rt))
        return(object)
    if (!missing(msLevel.))
        warning("Parameter 'msLevel.' currently ignored.")
    rt <- range(rt)
    ## Get index of spectra within the rt window.
    ## Subset using [
    ## Subset peaks
    ## Subset features
    ## Subset adjusted retention time
    if (!adjusted) {
        have_rt <- rtime(object, adjusted = FALSE, bySample = FALSE)
    } else {
        have_rt <- adjustedRtime(object, bySample = FALSE)
        if (is.null(have_rt))
            stop("No adjusted retention time available!")
    }
    keep_logical <- have_rt >= rt[1] & have_rt <= rt[2]
    msg <- paste0("Filter: select retention time [",
                  paste0(rt, collapse = "-"),
                  "] and MS level(s), ",
                  paste(base::unique(msLevel(object)),
                        collapse = " "))
    msg <- paste0(msg, " [", date(), "]")
    if (!any(keep_logical)) {
        res <- new("XCMSnExp")
        res@processingData@processing <- c(res@processingData@processing, msg)
        return(res)
    }

    ## Extract the stuff we want to keep
    ## mfd <- as(.copy_env(object@msFeatureData), "MsFeatureData")
    newMfd <- new("MsFeatureData")
    ph <- processHistory(object)
    ## 1) Subset peaks within the retention time range and peak groups.
    keep_fts <- numeric()
    if (hasChromPeaks(object)) {
        ftrt <- chromPeaks(object)[, "rt"]
        if (!adjusted & hasAdjustedRtime(object)) {
            ## Have to convert the rt before subsetting.
            fts <- .applyRtAdjToChromPeaks(chromPeaks(object),
                                           rtraw = rtime(object, bySample = TRUE),
                                           rtadj = rtime(object, bySample = TRUE,
                                                         adjusted = FALSE))
            ftrt <- fts[, "rt"]
        }
        keep_fts <- base::which(ftrt >= rt[1] & ftrt <= rt[2])
        if (length(keep_fts))
            newMfd <- .filterChromPeaks(object, idx = keep_fts)
            ## features(newMfd) <- features(object)[keep_fts, , drop = FALSE]
        else
            ph <- dropProcessHistoriesList(ph,
                                           type = c(.PROCSTEP.PEAK.DETECTION,
                                                    .PROCSTEP.PEAK.GROUPING,
                                                    .PROCSTEP.PEAK.FILLING,
                                                    .PROCSTEP.CALIBRATION))
    }
    ## 2) Subset adjusted retention time
    ## if (hasAdjustedRtime(object) & length(keep_fts)) {
    if (hasAdjustedRtime(object)) {
        ## Subset the adjusted retention times (which are stored as a list of
        ## rts by file):
        keep_by_file <- base::split(keep_logical, fromFile(object))
        adj_rt <- base::mapply(FUN = function(y, z) {
            return(y[z])
        }, y = adjustedRtime(object, bySample = TRUE), z = keep_by_file,
        SIMPLIFY = FALSE)
        adjustedRtime(newMfd) <- adj_rt
    }
    ## 3) Subset the OnDiskMSnExp part
    ## Note: this is still slightly faster than dropping the msFeatureData and
    ## calling it on the XCMSnExp!
    tmp <- as(object, "OnDiskMSnExp")[base::which(keep_logical)]
    object <- as(tmp, "XCMSnExp")
    ## Put the stuff back
    object@processingData@processing <- c(object@processingData@processing, msg)
    lockEnvironment(newMfd, bindings = TRUE)
    object@msFeatureData <- newMfd
    object@.processHistory <- ph
    if (validObject(object))
        return(object)
})


#' @description
#'
#' The \code{normalize} method performs basic normalization of
#' spectra intensities. See \code{\link{normalize}} documentation
#' in the \code{MSnbase} package for details and examples.
#'
#' @param method For \code{normalize}: \code{character(1)} specifying the
#'     normalization method. See \code{\link{normalize}} in the \code{MSnbase}
#'     package for details.
#'     For \code{pickPeaks}: \code{character(1)} defining the method. See
#'     \code{\link{pickPeaks}} for options. For \code{smooth}:
#'     \code{character(1)} defining the method. See
#'     \code{\link{smooth}} in the \code{MSnbase} package for options and
#'     details.
#'
#' @rdname XCMSnExp-inherited-methods
setMethod("normalize", "XCMSnExp", function(object, method = c("max", "sum"),
                                            ...) {
    if (hasAdjustedRtime(object) | hasFeatures(object) |
        hasChromPeaks(object)) {
        ## object@.processHistory <- list()
        ## object@msFeatureData <- new("MsFeatureData")
        object <- dropAdjustedRtime(object)
        object <- dropFeatureDefinitions(object)
        object <- dropChromPeaks(object)
        warning("Removed preprocessing results")
    }
    callNextMethod()
})

#' @description
#'
#' The \code{pickPeaks} method performs peak picking. See
#' \code{\link{pickPeaks}} documentation for details and examples.
#'
#' @param halfWindowSize For \code{pickPeaks} and \code{smooth}:
#'     \code{integer(1)} defining the window size for the peak picking. See
#'     \code{\link{pickPeaks}} and \code{\link{smooth}} in the \code{MSnbase}
#'     package for details and options.
#'
#' @param SNR For \code{pickPeaks}: \code{numeric(1)} defining the signal to
#'     noise ratio to be considered. See \code{\link{pickPeaks}}
#'     documentation for details.
#'
#' @param ... Optional additional arguments.
#'
#' @rdname XCMSnExp-inherited-methods
setMethod("pickPeaks", "XCMSnExp", function(object, halfWindowSize = 3L,
                                            method = c("MAD", "SuperSmoother"),
                                            SNR = 0L, ...) {
    if (hasAdjustedRtime(object) | hasFeatures(object) |
        hasChromPeaks(object)) {
        ## object@.processHistory <- list()
        ## object@msFeatureData <- new("MsFeatureData")
        object <- dropAdjustedRtime(object)
        object <- dropFeatureDefinitions(object)
        object <- dropChromPeaks(object)
        warning("Removed preprocessing results")
    }
    callNextMethod()
})

#' @description
#'
#' The \code{removePeaks} method removes mass peaks (intensities)
#' lower than a threshold. Note that these peaks refer to \emph{mass}
#' peaks, which are different to the chromatographic peaks detected and
#' analyzed in a metabolomics experiment! See
#' \code{\link{removePeaks}} documentation for details and
#' examples.
#'
#' @param t For \code{removePeaks}: either a \code{numeric(1)} or \code{"min"}
#'     defining the threshold (method) to be used. See
#'     \code{\link{removePeaks}} for details.
#'
#' @rdname XCMSnExp-inherited-methods
setMethod("removePeaks", "XCMSnExp", function(object, t = "min", verbose = FALSE,
                                              msLevel.) {
    if (hasAdjustedRtime(object) | hasFeatures(object) |
        hasChromPeaks(object)) {
        ## object@.processHistory <- list()
        ## object@msFeatureData <- new("MsFeatureData")
        object <- dropAdjustedRtime(object)
        object <- dropFeatureDefinitions(object)
        object <- dropChromPeaks(object)
        warning("Removed preprocessing results")
    }
    callNextMethod()
})

#' @description
#'
#' The \code{smooth} method smooths spectra. See
#' \code{\link{smooth}} documentation in \code{MSnbase} for details and
#' examples.
#'
#' @rdname XCMSnExp-inherited-methods
setMethod("smooth", "XCMSnExp", function(x, method = c("SavitzkyGolay",
                                                       "MovingAverage"),
                                         halfWindowSize = 2L, verbose = FALSE,
                                         ...) {
    if (hasAdjustedRtime(x) | hasFeatures(x) |
        hasChromPeaks(x)) {
        ## x@.processHistory <- list()
        ## x@msFeatureData <- new("MsFeatureData")
        x <- dropAdjustedRtime(x)
        x <- dropFeatureDefinitions(x)
        x <- dropChromPeaks(x)
        warning("Removed preprocessing results")
    }
    callNextMethod()
})

#' @aliases setAs
#' 
#' @rdname XCMSnExp-class
#'
#' @name XCMSnExp-class
setAs(from = "XCMSnExp", to = "xcmsSet", def = .XCMSnExp2xcmsSet)


#' @title Peak grouping/correspondence based on time dimension peak densities
#'
#' @description
#'
#' \code{groupChromPeaks,XCMSnExp,PeakDensityParam}:
#' performs correspondence (peak grouping within and across samples) within
#' in mz dimension overlapping slices of MS data based on the density
#' distribution of the identified chromatographic peaks in the slice along
#' the time axis.
#'
#' @note
#'
#' Calling \code{groupChromPeaks} on an \code{XCMSnExp} object will cause
#' all eventually present previous correspondence results to be dropped.
#'
#' @param object For \code{groupChromPeaks}: an \code{\link{XCMSnExp}} object
#'     containing the results from a previous peak detection analysis (see
#'     \code{\link{findChromPeaks}}).
#'
#'     For all other methods: a \code{PeakDensityParam} object.
#' 
#' @param param A \code{PeakDensityParam} object containing all settings for
#'     the peak grouping algorithm.
#'
#' @return
#'
#' For \code{groupChromPeaks}: a \code{\link{XCMSnExp}} object with the
#' results of the correspondence analysis. The definition of the resulting
#' mz-rt features can be accessed with the \code{\link{featureDefinitions}}
#' method.
#' 
#' @seealso \code{\link{XCMSnExp}} for the object containing the results of
#'     the correspondence.
#' 
#' @rdname groupChromPeaks-density
setMethod("groupChromPeaks",
          signature(object = "XCMSnExp", param = "PeakDensityParam"),
          function(object, param) {
              if (!hasChromPeaks(object))
                  stop("No chromatographic peak detection results in 'object'! ",
                       "Please perform first a peak detection using the ",
                       "'findChromPeaks' method.")
              ## Get rid of any previous results.
              if (hasFeatures(object))
                  object <- dropFeatureDefinitions(object)
              ## Check if we've got any sample groups:
              if (length(sampleGroups(param)) == 0) {
                  sampleGroups(param) <- rep(1, length(fileNames(object)))
                  message("Empty 'sampleGroups' in 'param', assuming all ",
                          "samples to be in the same group.")
              } else {
                  ## Check that the sampleGroups are OK
                  if (length(sampleGroups(param)) != length(fileNames(object)))
                      stop("The 'sampleGroups' value in the provided 'param' ",
                           "class does not match the number of available files/",
                           "samples!")
              }
              startDate <- date()
              res <- do_groupChromPeaks_density(chromPeaks(object),
                                                sampleGroups = sampleGroups(param),
                                                bw = bw(param),
                                                minFraction = minFraction(param),
                                                minSamples = minSamples(param),
                                                binSize = binSize(param),
                                                maxFeatures = maxFeatures(param))
              xph <- XProcessHistory(param = param, date. = startDate,
                                     type. = .PROCSTEP.PEAK.GROUPING,
                                     fileIndex = 1:length(fileNames(object)))
              object <- addProcessHistory(object, xph)              
              ## Add the results.
              df <- DataFrame(res$featureDefinitions)
              if (nrow(df) == 0)
                  stop("Unable to group any chromatographic peaks. You might ",
                       "have to adapt your settings.")
              df$peakidx <- res$peakIndex
              if (nrow(df) > 0)
                  rownames(df) <- .featureIDs(nrow(df))
              featureDefinitions(object) <- df
              if (validObject(object))
                  return(object)
          })


#' @title Single-spectrum non-chromatography MS data peak grouping
#'
#' @description
#'
#' \code{groupChromPeaks,XCMSnExp,MzClustParam}:
#' performs high resolution peak grouping for single spectrum
#' metabolomics data.
#'
#' @note Calling \code{groupChromPeaks} on an \code{XCMSnExp} object will cause
#'     all eventually present previous correspondence results to be dropped.
#'
#' @param object For \code{groupChromPeaks}: an \code{\link{XCMSnExp}} object
#'     containing the results from a previous chromatographic peak detection
#'     analysis (see \code{\link{findChromPeaks}}).
#'
#'     For all other methods: a \code{MzClustParam} object.
#' 
#' @param param A \code{MzClustParam} object containing all settings for
#'     the peak grouping algorithm.
#'
#' @return
#'
#' For \code{groupChromPeaks}: a \code{\link{XCMSnExp}} object with the
#' results of the peak grouping step (i.e. the features). These can be
#' accessed with the \code{\link{featureDefinitions}} method.
#' 
#' @seealso \code{\link{XCMSnExp}} for the object containing the results of
#'     the peak grouping.
#' 
#' @rdname groupChromPeaks-mzClust
setMethod("groupChromPeaks",
          signature(object = "XCMSnExp", param = "MzClustParam"),
          function(object, param) {
              if (!hasChromPeaks(object))
                  stop("No chromatographic peak detection results in 'object'! ",
                       "Please perform first a peak detection using the ",
                       "'findChromPeak' method.")
              ## I'm expecting a single spectrum per file!
              rtL <- split(rtime(object), f = fromFile(object))
              if (any(lengths(rtL) > 1))
                  stop("'object' contains multiple spectra per sample! This ",
                       "algorithm does only work for single spectra ",
                       "files/samples!")
              ## Get rid of any previous results.
              if (hasFeatures(object))
                  object <- dropFeatureDefinitions(object)
              ## Check if we've got any sample groups:
              if (length(sampleGroups(param)) == 0) {
                  sampleGroups(param) <- rep(1, length(fileNames(object)))
                  message("Empty 'sampleGroups' in 'param', assuming all ",
                          "samples to be in the same group.")
              } else {
                  ## Check that the sampleGroups are OK
                  if (length(sampleGroups(param)) != length(fileNames(object)))
                      stop("The 'sampleGroups' value in the provided 'param' ",
                           "class does not match the number of available files/",
                           "samples!")
              }
              startDate <- date()
              res <- do_groupPeaks_mzClust(chromPeaks(object),
                                           sampleGroups = sampleGroups(param),
                                           ppm = ppm(param),
                                           absMz = absMz(param),
                                           minFraction = minFraction(param),
                                           minSamples = minSamples(param))
              xph <- XProcessHistory(param = param, date. = startDate,
                                     type. = .PROCSTEP.PEAK.GROUPING,
                                     fileIndex = 1:length(fileNames(object)))
              object <- addProcessHistory(object, xph)              
              ## Add the results.
              df <- DataFrame(res$featureDefinitions)
              df$peakidx <- res$peakIndex
              if (nrow(df) == 0)
                  stop("Unable to group any chromatographic peaks. You might ",
                       "have to adapt your settings.")
              if (nrow(df) > 0)
                  rownames(df) <- .featureIDs(nrow(df))
              featureDefinitions(object) <- df
              if (validObject(object))
                  return(object)
          })


#' @title Peak grouping/correspondence based on proximity in the mz-rt space
#'
#' @description
#'
#' \code{groupChromPeaks,XCMSnExp,NearestPeaksParam}:
#' performs peak grouping based on the proximity between chromatographic
#' peaks from different samples in the mz-rt range.
#'
#' @note
#'
#' Calling \code{groupChromPeaks} on an \code{XCMSnExp} object will cause
#' all eventually present previous alignment results to be dropped.
#'
#' @param object For \code{groupChromPeaks}: an \code{\link{XCMSnExp}} object
#'     containing the results from a previous chromatographic peak detection
#'     analysis (see \code{\link{findChromPeaks}}).
#'
#'     For all other methods: a \code{NearestPeaksParam} object.
#' 
#' @param param A \code{NearestPeaksParam} object containing all settings for
#'     the peak grouping algorithm.
#'
#' @return
#'
#' For \code{groupChromPeaks}: a \code{\link{XCMSnExp}} object with the
#' results of the peak grouping/correspondence step (i.e. the mz-rt
#' features). These can be accessed with the
#' \code{\link{featureDefinitions}} method.
#' 
#' @seealso \code{\link{XCMSnExp}} for the object containing the results of
#'     the peak grouping.
#' 
#' @rdname groupChromPeaks-nearest
setMethod("groupChromPeaks",
          signature(object = "XCMSnExp", param = "NearestPeaksParam"),
          function(object, param) {
              if (!hasChromPeaks(object))
                  stop("No chromatographic peak detection results in 'object'! ",
                       "Please perform first a peak detection using the ",
                       "'findChromPeaks' method.")
              ## Get rid of any previous results.
              if (hasFeatures(object))
                  object <- dropFeatureDefinitions(object)
              ## Check if we've got any sample groups:
              if (length(sampleGroups(param)) == 0) {
                  sampleGroups(param) <- rep(1, length(fileNames(object)))
                  message("Empty 'sampleGroups' in 'param', assuming all ",
                          "samples to be in the same group.")
              } else {
                  ## Check that the sampleGroups are OK
                  if (length(sampleGroups(param)) != length(fileNames(object)))
                      stop("The 'sampleGroups' value in the provided 'param' ",
                           "class does not match the number of available files/",
                           "samples!")
              }
              startDate <- date()
              res <- do_groupChromPeaks_nearest(chromPeaks(object),
                                                sampleGroups = sampleGroups(param),
                                                mzVsRtBalance = mzVsRtBalance(param),
                                                absMz = absMz(param),
                                                absRt = absRt(param),
                                                kNN = kNN(param))
              xph <- XProcessHistory(param = param, date. = startDate,
                                     type. = .PROCSTEP.PEAK.GROUPING,
                                     fileIndex = 1:length(fileNames(object)))
              object <- addProcessHistory(object, xph)
              ## Add the results.
              df <- DataFrame(res$featureDefinitions)
              df$peakidx <- res$peakIndex
              if (nrow(df) == 0)
                  stop("Unable to group any chromatographic peaks. You might ",
                       "have to adapt your settings.")
              if (nrow(df) > 0)
                  rownames(df) <- .featureIDs(nrow(df))
              featureDefinitions(object) <- df
              if (validObject(object))
                  return(object)
          })

#' @title Retention time correction based on alignment of house keeping peak
#' groups
#'
#' @description
#'
#' \code{adjustRtime,XCMSnExp,PeakGroupsParam}:
#' performs retention time correction based on the alignment of peak groups
#' (features) found in all/most samples. The correction function identified
#' on these peak groups is applied to the retention time of all spectra in
#' the object, i.e. retention times of all spectra, also MS level > 1 are
#' adjusted.
#'
#' @note
#'
#' This method requires that a correspondence analysis has been performed
#' on the data, i.e. that grouped chromatographic peaks/features are present
#' (see \code{\link{groupChromPeaks}} for details).
#'
#' Calling \code{adjustRtime} on an \code{XCMSnExp} object will cause all
#' peak grouping (correspondence) results and any previous retention time
#' adjustments to be dropped.
#' In some instances, the \code{adjustRtime,XCMSnExp,PeakGroupsParam}
#' re-adjusts adjusted retention times to ensure them being in the same
#' order than the raw (original) retention times.
#'
#' @param object For \code{adjustRtime}: an \code{\link{XCMSnExp}} object
#'     containing the results from a previous chromatographic peak detection
#'     (see \code{\link{findChromPeaks}}) and alignment analysis (see
#'     \code{\link{groupChromPeaks}}).
#'
#'     For all other methods: a \code{PeakGroupsParam} object.
#' 
#' @param param A \code{PeakGroupsParam} object containing all settings for
#'     the retention time correction method..
#'
#' @return
#'
#' For \code{adjustRtime}: a \code{\link{XCMSnExp}} object with the
#' results of the retention time adjustment step. These can be accessed
#' with the \code{\link{adjustedRtime}} method. Retention time correction
#' does also adjust the retention time of the identified chromatographic
#' peaks (accessed \emph{via} \code{\link{chromPeaks}}. Note that retention
#' time correction drops all previous alignment results from the result
#' object.
#' 
#' @seealso \code{\link{XCMSnExp}} for the object containing the results of
#'     the alignment.
#' 
#' @rdname adjustRtime-peakGroups
setMethod("adjustRtime",
          signature(object = "XCMSnExp", param = "PeakGroupsParam"),
          function(object, param) {
              if (hasAdjustedRtime(object)) {
                  message("Removing previous alignment results")
                  object <- dropAdjustedRtime(object)
              }
              if (!hasChromPeaks(object))
                  stop("No chromatographic peak detection results in 'object'! ",
                       "Please perform first a peak detection using the ",
                       "'findChromPeaks' method.")
              if (!hasFeatures(object))
                  stop("No feature definitions found in 'object'! Please ",
                       "perform first a peak grouping using the ",
                       "'groupChromPeak' method.")
              startDate <- date()
              ## If param does contain a peakGroupsMatrix extract that one,
              ## otherwise generate it.
              if (nrow(peakGroupsMatrix(param)))
                  pkGrpMat <- peakGroupsMatrix(param)
              else
                  pkGrpMat <- adjustRtimePeakGroups(object, param = param)
              res <- do_adjustRtime_peakGroups(
                  chromPeaks(object),
                  peakIndex = featureDefinitions(object)$peakidx,
                  rtime = rtime(object, bySample = TRUE),
                  minFraction = minFraction(param),
                  extraPeaks = extraPeaks(param),
                  smooth = smooth(param),
                  span = span(param),
                  family = family(param),
                  peakGroupsMatrix = pkGrpMat
              )
              ## Add the pkGrpMat that's being used to the param object.
              peakGroupsMatrix(param) <- pkGrpMat
              ## Dropping the peak groups but don't remove its process history
              ## step.
              ph <- processHistory(object, type = .PROCSTEP.PEAK.GROUPING)
              object <- dropFeatureDefinitions(object)
              ## Add the results. adjustedRtime<- should also fix the retention
              ## times for the peaks! Want to keep also the latest alignment
              ## information
              adjustedRtime(object) <- res
              if (length(ph)) {
                  object <- addProcessHistory(object, ph[[length(ph)]])
              }
              ## Add the process history step, get the msLevel from the peak
              ## detection step.
              ph <- processHistory(object, type = .PROCSTEP.PEAK.DETECTION)
              xph <- XProcessHistory(param = param, date. = startDate,
                                     type. = .PROCSTEP.RTIME.CORRECTION,
                                     fileIndex = 1:length(fileNames(object)),
                                     msLevel = msLevel(ph[[length(ph)]]))
              object <- addProcessHistory(object, xph)
              if (validObject(object))
                  object
          })


#' @title Align retention times across samples using Obiwarp
#'
#' @description
#'
#' \code{adjustRtime,XCMSnExp,ObiwarpParam}:
#' performs retention time correction/alignment based on the total mz-rt
#' data using the \emph{obiwarp} method.
#'
#' @note
#'
#' Alignment using obiwarp is performed on the retention time of spectra
#' of on MS level. Retention times for spectra of other MS levels are
#' subsequently adjusted based on the adjustment function defined on the
#' retention times of the spectra of MS level \code{msLevel}.
#' 
#' Calling \code{adjustRtime} on an \code{XCMSnExp} object will cause
#' all peak grouping (correspondence) results and any previous retention
#' time adjustment results to be dropped.
#'
#' @param object For \code{adjustRtime}: an \code{\link{XCMSnExp}} object.
#'
#'     For all other methods: a \code{ObiwarpParam} object.
#' 
#' @param param A \code{ObiwarpParam} object containing all settings for
#'     the alignment method.
#'
#' @param msLevel \code{integer} defining the MS level on which the retention
#'     time should be performed.
#' 
#' @return
#'
#' For \code{adjustRtime,XCMSnExp,ObiwarpParam}: a
#' \code{\link{XCMSnExp}} object with the results of the retention time
#' adjustment step. These can be accessed with the
#' \code{\link{adjustedRtime}} method. Retention time correction does also
#' adjust the retention time of the identified chromatographic peaks
#' (accessed \emph{via} \code{\link{chromPeaks}}. Note that retention time
#' correction drops all previous peak grouping results from the result
#' object.
#'
#' For \code{adjustRtime,OnDiskMSnExp,ObiwarpParam}: a \code{numeric} with
#' the adjusted retention times per spectra (in the same order than
#' \code{rtime}).
#' 
#' @seealso \code{\link{XCMSnExp}} for the object containing the results of
#'     the alignment.
#'
#' @references
#' John T. Prince and Edward M. Marcotte. "Chromatographic Alignment of
#' ESI-LC-MS Proteomic Data Sets by Ordered Bijective Interpolated Warping"
#' \emph{Anal. Chem.} 2006, 78 (17), 6140-6152.
#' 
#' @rdname adjustRtime-obiwarp
setMethod("adjustRtime",
          signature(object = "XCMSnExp", param = "ObiwarpParam"),
          function(object, param, msLevel = 1L) {
              ## Drop adjusted retention times if there are some.
              if (hasAdjustedRtime(object))
                  object <- dropAdjustedRtime(object)
              ## We don't require any detected or aligned peaks.
              startDate <- date()
              res <- adjustRtime(as(object, "OnDiskMSnExp"), param = param,
                                 msLevel = msLevel)
              ## res <- .obiwarp(as(object, "OnDiskMSnExp"), param = param)
              ## Dropping the feature groups.
              object <- dropFeatureDefinitions(object)
              ## Add the results. adjustedRtime<- should also fix the retention
              ## times for the peaks! Want to keep also the latest alignment
              ## information
              adjustedRtime(object) <- unname(split(res, fromFile(object)))
              ## Add the process history step.
              xph <- XProcessHistory(param = param, date. = startDate,
                                     type. = .PROCSTEP.RTIME.CORRECTION,
                                     fileIndex = 1:length(fileNames(object)),
                                     msLevel = msLevel)
              object <- addProcessHistory(object, xph)
              if (validObject(object))
                  object
          })

#' @rdname XCMSnExp-class
setMethod("profMat", signature(object = "XCMSnExp"), function(object,
                                                              method = "bin",
                                                              step = 0.1,
                                                              baselevel = NULL,
                                                              basespace = NULL,
                                                              mzrange. = NULL,
                                                              fileIndex,
                                                              ...) {
    ## We want to coerce that as OnDiskMSnExp so we don't slow down in the
    ## filterFile, that would, if rt adjustments are present, revert the whole
    ## thing.
    return(profMat(as(object, "OnDiskMSnExp"), method = method, step = step,
                   baselevel = baselevel, basespace = basespace,
                   mzrange. = mzrange., fileIndex = fileIndex, ...))
})


#' @aliases featureValues
#' 
#' @title Accessing mz-rt feature data values
#' 
#' @description
#'
#' \code{featureValues,XCMSnExp} : extract a \code{matrix} for
#' feature values with rows representing features and columns samples.
#' Parameter \code{value} allows to define which column from the
#' \code{\link{chromPeaks}} matrix should be returned. Multiple
#' chromatographic peaks from the same sample can be assigned to a feature.
#' Parameter \code{method} allows to specify the method to be used in such
#' cases to chose from which of the peaks the value should be returned.
#'
#' @note
#'
#' This method is equivalent to the \code{\link{groupval}} for
#' \code{xcmsSet} objects.
#' 
#' @param object A \code{\link{XCMSnExp}} object providing the feature
#'     definitions.
#' 
#' @param method \code{character} specifying the method to resolve
#'     multi-peak mappings within the same sample, i.e. to define the
#'     \emph{representative} peak for a feature in samples where more than
#'     one peak was assigned to the feature. If \code{"medret"}: select the
#'     peak closest to the median retention time of the feature.
#'     If \code{"maxint"}: select the peak yielding the largest signal. If
#'     \code{"sum"}: sum the values (only if \code{value} is \code{"into"} or
#'     \code{"maxo"}.
#'
#' @param value \code{character} specifying the name of the column in
#'     \code{chromPeaks(object)} that should be returned or \code{"index"} (the
#'     default) to return the index of the peak in the
#'     \code{chromPeaks(object)} matrix corresponding to the
#'     \emph{representative} peak for the feature in the respective sample.
#'
#' @param intensity \code{character} specifying the name of the column in the
#'     \code{chromPeaks(objects)} matrix containing the intensity value of the
#'     peak that should be used for the conflict resolution if
#'     \code{method = "maxint"}.
#'
#' @param filled \code{logical(1)} specifying whether values for filled-in
#'     peaks should be returned or not. If \code{filled = FALSE}, an \code{NA}
#'     is returned in the matrix for the respective peak. See
#'     \code{\link{fillChromPeaks}} for details on peak filling. 
#' 
#' @return
#'
#' For \code{featureValues}: a \code{matrix} with
#' feature values, columns representing samples, rows features. The order
#' of the features matches the order found in the
#' \code{featureDefinitions(object)} \code{DataFrame}. The rownames of the
#' \code{matrix} are the same than those of the \code{featureDefinitions}
#' \code{DataFrame}. \code{NA} is reported for features without
#' corresponding chromatographic peak in the respective sample(s).
#' 
#' @author Johannes Rainer
#' 
#' @seealso
#' \code{\link{XCMSnExp}} for information on the data object.
#' \code{\link{featureDefinitions}} to extract the \code{DataFrame} with the
#' feature definitions.
#' \code{\link{hasFeatures}} to evaluate whether the
#' \code{\link{XCMSnExp}} provides feature definitions.
#' \code{\link{groupval}} for the equivalent method on \code{xcmsSet} objects.
#' 
#' @rdname XCMSnExp-peak-grouping-results
setMethod("featureValues",
          signature(object = "XCMSnExp"),
          function(object, method = c("medret", "maxint", "sum"),
                   value = "index", intensity = "into", filled = TRUE) {
              ## Input argument checkings
              if (!hasFeatures(object))
                  stop("No peak groups present! Use 'groupChromPeaks' first.")
              if (!hasChromPeaks(object))
                  stop("No detected chromatographic peaks present! Use ",
                       "'findChromPeaks' first.")
              method <- match.arg(method)
              if (method == "sum" & !(value %in% c("into", "maxo")))
                  stop("method 'sum' is only allowed if value is set to 'into'",
                       " or 'maxo'")
              fNames <- basename(fileNames(object))
              nSamples <- seq_along(fNames)
              ## Copy all of the objects to avoid costly S4 method calls -
              ## improves speed at the cost of higher memory demand.
              fts <- chromPeaks(object)

              ## issue #157: replace all values for filled-in peaks with NA
              if (!filled)
                  fts[fts[, "is_filled"] == 1, ] <- NA
              grps <- featureDefinitions(object)
              ftIdx <- grps$peakidx
              ## Match columns
              idx_rt <- match("rt", colnames(fts))
              idx_int <- match(intensity, colnames(fts))
              idx_samp <- match("sample", colnames(fts))
              
              vals <- matrix(nrow = length(ftIdx), ncol = length(nSamples))
              
              if (method == "sum") {
                  for (i in seq_along(ftIdx)) {
                      cur_pks <- fts[ftIdx[[i]], c(value, "sample")]
                      int_sum <- split(cur_pks[, value], cur_pks[, "sample"])
                      vals[i, as.numeric(names(int_sum))] <-
                          unlist(lapply(int_sum, base::sum), use.names = FALSE)
                  }
              } else {
                  ## Selecting only a single one.
                  ## Get the indices for the elements.
                  if (method == "medret") {
                      medret <- grps$rtmed
                      for (i in seq_along(ftIdx)) {
                          gidx <- ftIdx[[i]][
                              base::order(base::abs(fts[ftIdx[[i]],
                                                        idx_rt] - medret[i]))]
                          vals[i, ] <- gidx[
                              base::match(nSamples, fts[gidx, idx_samp])]
                      }
                  }
                  if (method == "maxint") {
                      for (i in seq_along(ftIdx)) {
                          gidx <- ftIdx[[i]][
                              base::order(fts[ftIdx[[i]], idx_int],
                                          decreasing = TRUE)]
                          vals[i, ] <- gidx[base::match(nSamples,
                                                        fts[gidx, idx_samp])]
                      }
                  }
                  if (value != "index") {
                      if (!any(colnames(fts) == value))
                          stop("Column '", value, "' not present in the ",
                               "chromatographic peaks matrix!")
                      vals <- fts[vals, value]
                      dim(vals) <- c(length(ftIdx), length(nSamples))
                  }
              }
              
              colnames(vals) <- fNames
              rownames(vals) <- rownames(grps)
              vals
})
## #' @rdname XCMSnExp-peak-grouping-results
## setMethod("groupval",
##           signature(object = "XCMSnExp"),
##           function(object, method = c("medret", "maxint"), value = "index",
##                    intensity = "into") {
##               featureValues(object = object, method = method, value = value,
##                             intensity = intensity)
## })


#' @aliases chromatogram
#' 
#' @title Extracting chromatograms
#'
#' @description
#'
#' \code{chromatogram}: the method allows to extract
#' chromatograms from \code{\link{OnDiskMSnExp}} and
#' \code{\link{XCMSnExp}} objects. See also the
#' \code{\link{chromatogram}} implementation for
#' \code{\link{OnDiskMSnExp}} in the \code{MSnbase} package.
#'
#' @details
#'
#' Arguments \code{rt} and \code{mz} allow to specify the MS
#' data slice from which the chromatogram should be extracted.
#' The parameter \code{aggregationSum} allows to specify the function to be
#' used to aggregate the intensities across the mz range for the same
#' retention time. Setting \code{aggregationFun = "sum"} would e.g. allow
#' to calculate the \emph{total ion chromatogram} (TIC),
#' \code{aggregationFun = "max"} the \emph{base peak chromatogram} (BPC).
#' The length of the extracted \code{\link{Chromatogram}} object,
#' i.e. the number of available data points, corresponds to the number of
#' scans/spectra measured in the specified retention time range. If in a
#' specific scan (for a give retention time) no signal was measured in the
#' specified mz range, a \code{NA_real_} is reported as intensity for the
#' retention time (see Notes for more information). This can be changed
#' using the \code{missing} parameter. 
#'
#' @note
#'
#' \code{\link{Chromatogram}} objects extracted with
#' \code{chromatogram}
#' contain \code{NA_real_} values if, for a given retention time, no 
#' signal was measured in the specified mz range. If no spectrum/scan is
#' present in the defined retention time window a \code{Chromatogram} object
#' of length 0 is returned.
#'
#' For \code{\link{XCMSnExp}} objects, if adjusted retention times are
#' available, the \code{chromatogram} method will by default report
#' and use these (for the subsetting based on the provided parameter
#' \code{rt}). This can be overwritten with the parameter
#' \code{adjustedRtime}.
#' 
#' @param object Either a \code{\link{OnDiskMSnExp}} or
#'     \code{\link{XCMSnExp}} object from which the chromatograms should be
#'     extracted.
#'
#' @param rt \code{numeric(2)} or two-column \code{matrix} defining the lower
#'     and upper boundary for the retention time range(s). If not specified,
#'     the full retention time range of the original data will be used.
#'     It is also possible to submit a \code{numeric(1)} in which case
#'     \code{range} is called on it to transform it to a \code{numeric(2)}.
#'
#' @param mz \code{numeric(2)} or two-column \code{matrix} defining the lower
#'     and upper mz value for the MS data slice(s). If not specified, the
#'     chromatograms will be calculated on the full mz range.
#'     It is also possible to submit a \code{numeric(1)} in which case
#'     \code{range} is called on it to transform it to a \code{numeric(2)}.
#'
#' @param adjustedRtime For \code{chromatogram,XCMSnExp}: whether the
#'     adjusted (\code{adjustedRtime = TRUE}) or raw retention times
#'     (\code{adjustedRtime = FALSE}) should be used for filtering and returned
#'     in the resulting \code{\link{Chromatogram}} object. Adjusted
#'     retention times are used by default if available.
#'
#' @param aggregationFun \code{character} specifying the function to be used to
#'     aggregate intensity values across the mz value range for the same
#'     retention time. Allowed values are \code{"sum"}, \code{"max"},
#'     \code{"mean"} and \code{"min"}.
#'
#' @param missing \code{numeric(1)} allowing to specify the intensity value to
#'     be used if for a given retention time no signal was measured within the
#'     mz range of the corresponding scan. Defaults to \code{NA_real_} (see also
#'     Details and Notes sections below). Use \code{missing = 0} to resemble the
#'     behaviour of the \code{getEIC} from the \code{old} user interface.
#'
#' @param msLevel \code{integer} specifying the MS level from which the
#'     chromatogram should be extracted. Defaults to \code{msLevel = 1L}.
#' 
#' @param BPPARAM Parallelisation backend to be used, which will
#'     depend on the architecture. Default is
#'     \code{BiocParallel::bparam()}.
#'
#' @return
#'
#' \code{chromatogram} returns a \code{\link{Chromatograms}} object with
#' the number of columns corresponding to the number of files in
#' \code{object} and number of rows the number of specified ranges (i.e.
#' number of rows of matrices provided with arguments \code{mz} and/or
#' \code{rt}).
#' 
#' @author Johannes Rainer
#'
#' @seealso \code{\link{XCMSnExp}} for the data object.
#'     \code{\link{Chromatogram}} for the object representing
#'     chromatographic data.
#'
#'     \code{\link{Chromatograms}} for the object allowing to arrange
#'     multiple \code{Chromatogram} objects.
#'
#'     \code{\link{plot}} to plot a \code{Chromatogram} or
#'     \code{Chromatograms} objects.
#'
#'     \code{\link{as}} (\code{as(x, "data.frame")}) in \code{MSnbase}
#'     for a method to extract the MS data as \code{data.frame}.
#'
#' @export
#' 
#' @rdname chromatogram-method
#'
#' @examples
#' ## Read some files from the faahKO package.
#' library(xcms)
#' library(faahKO)
#' faahko_3_files <- c(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
#'                     system.file('cdf/KO/ko16.CDF', package = "faahKO"),
#'                     system.file('cdf/KO/ko18.CDF', package = "faahKO"))
#'
#' od <- readMSData(faahko_3_files, mode = "onDisk")
#'
#' ## Extract the ion chromatogram for one chromatographic peak in the data.
#' chrs <- chromatogram(od, rt = c(2700, 2900), mz = 335)
#'
#' chrs
#' 
#' ## Plot the chromatogram 
#' plot(rtime(chrs[1, 2]), intensity(chrs[1, 2]), type = "l", xlab = "rtime",
#'      ylab = "intensity", col = "000080")
#' for(i in c(1, 3)) {
#'   points(rtime(chrs[1, i]), intensity(chrs[1, i]), type = "l",
#'   col = "00000080")
#' }
#'
#' ## Plot the chromatogram using the dedicated plot method.
#' plot(chrs)
#'
#' ## Extract chromatograms for multiple ranges.
#' mzr <- matrix(c(335, 335, 344, 344), ncol = 2, byrow = TRUE)
#' rtr <- matrix(c(2700, 2900, 2600, 2750), ncol = 2, byrow = TRUE)
#' chrs <- chromatogram(od, mz = mzr, rt = rtr)
#'
#' chrs
#' 
#' ## Plot the extracted chromatograms
#' plot(chrs)
#'
#' ## Get access to all chromatograms for the second mz/rt range
#' chrs[1, ]
#'
#' ## Plot just that one
#' plot(chrs[1, , drop = FALSE])
setMethod("chromatogram",
          signature(object = "XCMSnExp"),
          function(object, rt, mz,
                   aggregationFun = "sum", missing = NA_real_,
                   msLevel = 1L, BPPARAM = bpparam(),
                   adjustedRtime = hasAdjustedRtime(object)) {
              ## Coerce to OnDiskMSnExp.
              if (adjustedRtime)
                  adj_rt <- rtime(object, adjusted = TRUE)
              object <- as(object, "OnDiskMSnExp")
              if (adjustedRtime) {
                  ## Replace the original rtime with adjusted ones...
                  object@featureData$retentionTime <- adj_rt
              }
              MSnbase::chromatogram(object, rt = rt, mz = mz,
                                    aggregationFun = aggregationFun,
                                    missing = missing, msLevel = msLevel,
                                    BPPARAM = BPPARAM)
          })

#' @rdname XCMSnExp-class
#'
#' @description
#'
#' \code{findChromPeaks} performs chromatographic peak detection
#' on the provided \code{XCMSnExp} objects. For more details see the method
#' for \code{\linkS4class{XCMSnExp}}. Note that the \code{findChromPeaks}
#' method for \code{XCMSnExp} objects removes previously identified
#' chromatographic peaks and aligned features. Previous alignment (retention
#' time adjustment) results are kept, i.e. chromatographic peak detection
#' is performed using adjusted retention times if the data was first
#' aligned using e.g. obiwarp (\code{\link{adjustRtime-obiwarp}}). 
#' 
#' @param param A \code{\link{CentWaveParam}}, \code{\link{MatchedFilterParam}},
#'     \code{\link{MassifquantParam}}, \code{\link{MSWParam}} or
#'     \code{\link{CentWavePredIsoParam}} object with the settings for the
#'     chromatographic peak detection algorithm.
#' 
#' @inheritParams findChromPeaks-centWave
setMethod("findChromPeaks",
          signature(object = "XCMSnExp", param = "Param"),
          function(object, param, BPPARAM = bpparam(),
                   return.type = "XCMSnExp", msLevel = 1L) {
              ## Remove previous correspondence results.
              if (hasFeatures(object)) {
                  message("Removed feature definitions.")
                  object <- dropFeatureDefinitions(
                      object,
                      keepAdjustedRtime = hasAdjustedRtime(object))
              }
              ## Remove previous chromatographic peaks.
              if (hasChromPeaks(object)) {
                  message("Removed previously identified chromatographic peaks.")
                  object <- dropChromPeaks(
                      object,
                      keepAdjustedRtime = hasAdjustedRtime(object))
              }
              meth <- selectMethod("findChromPeaks",
                                   signature = c(object = "OnDiskMSnExp",
                                                 param = class(param)))
              object <- do.call(meth, args = list(object = object,
                                                  param = param,
                                                  BPPARAM = BPPARAM,
                                                  return.type = return.type,
                                                  msLevel = msLevel))
              ## object@.processHistory <- list()
              if (validObject(object))
                  object
})

#' @aliases fillChromPeaks
#' 
#' @title Integrate areas of missing peaks
#'
#' @description
#'
#' Integrate signal in the mz-rt area of a feature (chromatographic
#' peak group) for samples in which no chromatographic peak for this
#' feature was identified and add it to the \code{chromPeaks}. Such peaks
#' will have a value of \code{1} in the \code{"is_filled"} column of the
#' \code{\link{chromPeaks}} matrix of the object.
#'
#' @details
#'
#' After correspondence (i.e. grouping of chromatographic peaks across
#' samples) there will always be features (peak groups) that do not include
#' peaks from every sample. The \code{fillChromPeaks} method defines
#' intensity values for such features in the missing samples by integrating
#' the signal in the mz-rt region of the feature. The mz-rt area is defined
#' by the median mz and rt start and end points of the other detected
#' chromatographic peaks for a given feature. Various parameters allow to
#' increase this area, either by a constant value (\code{fixedMz} and
#' \code{fixedRt}) or by a feature-relative amount (\code{expandMz} and
#' \code{expandRt}).
#' 
#' Adjusted retention times will be used if available.
#'
#' Based on the peak finding algorithm that was used to identify the
#' (chromatographic) peaks different internal functions are employed to
#' guarantee that the integrated peak signal matches as much as possible
#' the peak signal integration used during the peak detection. For peaks
#' identified with the \code{\link{matchedFilter}} method, signal
#' integration is performed on the \emph{profile matrix} generated with
#' the same settings used also during peak finding (using the same
#' \code{bin} size for example). For direct injection data and peaks
#' identified with the \code{\link{MSW}} algorithm signal is integrated
#' only along the mz dimension. For all other methods the complete (raw)
#' signal within the area defined by \code{"mzmin"}, \code{"mzmax"},
#' \code{"rtmin"} and \code{"rtmax"} is used. 
#' 
#' @note
#'
#' The reported \code{"mzmin"}, \code{"mzmax"}, \code{"rtmin"} and
#' \code{"rtmax"} for the filled peaks represents the actual MS area from
#' which the signal was integrated.
#' Note that no peak is filled in if no signal was present in a file/sample
#' in the respective mz-rt area. These samples will still show a \code{NA}
#' in the matrix returned by the \code{\link{featureValues}} method. This
#' is in contrast to the \code{\link{fillPeaks.chrom}} method that returned
#' an \code{"into"} and \code{"maxo"} of \code{0} for such peak areas.
#' Growing the mz-rt area using the \code{expandMz} and \code{expandRt}
#' might help to reduce the number of missing peak signals after filling.
#'
#' @param object \code{XCMSnExp} object with identified and grouped
#'     chromatographic peaks.
#'
#' @param param A \code{FillChromPeaksParam} object with all settings.
#' 
#' @param expandMz \code{numeric(1)} defining the value by which the mz width of
#'     peaks should be expanded. Each peak is expanded in mz direction by
#'     \code{expandMz *} their original mz width. A value of \code{0} means no
#'     expansion, a value of \code{1} grows each peak by 1 * the mz width of
#'     the peak resulting in peakswith twice their original size in mz
#'     direction (expansion by half mz width to both sides).
#'
#' @param expandRt \code{numeric(1)}, same as \code{expandRt} but for the
#'     retention time width.
#'
#' @param ppm \code{numeric(1)} optionally specifying a \emph{ppm} by which the
#'     mz width of the peak region should be expanded. For peaks with an mz
#'     width smaller than \code{mean(c(mzmin, mzmax)) * ppm / 1e6}, the
#'     \code{mzmin} will be replaced by
#'     \code{mean(c(mzmin, mzmax)) - (mean(c(mzmin, mzmax)) * ppm / 2 / 1e6)}
#'     and \code{mzmax} by
#'     \code{mean(c(mzmin, mzmax)) + (mean(c(mzmin, mzmax)) * ppm / 2 / 1e6)}.
#'     This is applied before eventually expanding the mz width using the
#'     \code{expandMz} parameter.
#'
#' @param fixedMz \code{numeric(1)} defining a constant factor by which the
#'     m/z width of each feature is to be expanded. The m/z width is expanded
#'     on both sides by \code{fixedMz} (i.e. \code{fixedMz} is subtracted
#'     from the lower m/z and added to the upper m/z). This expansion is
#'     applied \emph{after} \code{expandMz} and \code{ppm}.
#'
#' @param fixedRt \code{numeric(1)} defining a constant factor by which the
#'     retention time width of each factor is to be expanded. The rt width is
#'     expanded on both sides by \code{fixedRt} (i.e. \code{fixedRt} is
#'     subtracted from the lower rt and added to the upper rt). This
#'     expansion is applied \emph{after} \code{expandRt}.
#' 
#' @param BPPARAM Parallel processing settings.
#' 
#' @return
#'
#' A \code{\link{XCMSnExp}} object with previously missing
#' chromatographic peaks for features filled into its \code{chromPeaks}
#' matrix.
#'
#' @rdname fillChromPeaks
#' 
#' @author Johannes Rainer
#' 
#' @seealso \code{\link{groupChromPeaks}} for methods to perform the
#'     correspondence.
#'     \code{\link{dropFilledChromPeaks}} for the method to remove filled in peaks.
#'
#' @examples
#' 
#' ## Perform the peak detection using centWave on some of the files from the
#' ## faahKO package. Files are read using the readMSData from the MSnbase
#' ## package
#' library(faahKO)
#' library(xcms)
#' fls <- dir(system.file("cdf/KO", package = "faahKO"), recursive = TRUE,
#'            full.names = TRUE)
#' raw_data <- readMSData(fls[1:2], mode = "onDisk")
#'
#' ## Create a CentWaveParam object. Note that the noise is set to 10000 to
#' ## speed up the execution of the example - in a real use case the default
#' ## value should be used, or it should be set to a reasonable value.
#' cwp <- CentWaveParam(ppm = 20, noise = 10000, snthresh = 40)
#' 
#' res <- findChromPeaks(raw_data, param = cwp)
#'
#' ## Perform the correspondence. We assign all samples to the same group.
#' res <- groupChromPeaks(res,
#'     param = PeakDensityParam(sampleGroups = rep(1, length(fileNames(res)))))
#'
#' ## For how many features do we lack an integrated peak signal?
#' sum(is.na(featureValues(res)))
#'
#' ## Filling missing peak data using default settings.
#' res <- fillChromPeaks(res)
#'
#' ## Get the peaks that have been filled in:
#' fp <- chromPeaks(res)[chromPeaks(res)[, "is_filled"] == 1, ]
#' head(fp)
#' 
#' ## Did we get a signal for all missing peaks?
#' sum(is.na(featureValues(res)))
#'
#' ## No.
#'
#' ## Get the process history step along with the parameters used to perform
#' ## The peak filling: 
#' ph <- processHistory(res, type = "Missing peak filling")[[1]]
#' ph
#'
#' ## The parameter class:
#' ph@param
#' 
#' ## Drop the filled in peaks:
#' res <- dropFilledChromPeaks(res)
#'
#' ## Perform the peak filling with modified settings: allow expansion of the
#' ## mz range by a specified ppm and expanding the mz range by mz width/2 
#' prm <- FillChromPeaksParam(ppm = 40, expandMz = 0.5)
#' res <- fillChromPeaks(res, param = prm)
#'
#' ## Did we get a signal for all missing peaks?
#' sum(is.na(featureValues(res)))
#'
#' ## Still the same missing peaks.
setMethod("fillChromPeaks",
          signature(object = "XCMSnExp", param = "FillChromPeaksParam"),
          function(object, param, BPPARAM = bpparam()) {
              if (!hasFeatures(object))
                  stop("'object' does not provide feature definitions! Please ",
                       "run 'groupChromPeaks' first.")
              ## Don't do that if we have already filled peaks?
              if (.hasFilledPeaks(object))
                  message("Filled peaks already present, adding still missing",
                          " peaks.")
              
              startDate <- date()
              expandMz <- expandMz(param)
              expandRt <- expandRt(param)
              fixedMz <- fixedMz(param)
              fixedRt <- fixedRt(param)
              ppm <- ppm(param)
              ## Define or extend the peak area from which the signal should be
              ## extracted.
              ## Original code: use the median of the min/max rt and mz per peak.
              fdef <- featureDefinitions(object)
              aggFunLow <- median
              aggFunHigh <- median
              ## Note: we ensure in the downstream function that the rt range is
              ## within the rt range. For the mz range it doesn't matter.
              pkArea <- do.call(
                  rbind,
                  lapply(
                      fdef$peakidx, function(z) {
                          tmp <- chromPeaks(object)[z, c("rtmin", "rtmax",
                                                         "mzmin", "mzmax"),
                                                    drop = FALSE]
                          pa <- c(aggFunLow(tmp[, 1]), aggFunHigh(tmp[, 2]),
                                  aggFunLow(tmp[, 3]), aggFunHigh(tmp[, 4]))
                          ## Check if we have to apply ppm replacement:
                          if (ppm != 0) {
                              mzmean <- mean(pa[3:4])
                              tittle <- mzmean * (ppm / 2) / 1E6
                              if ((pa[4] - pa[3]) < (tittle * 2)) {
                                  pa[3] <- mzmean - tittle
                                  pa[4] <- mzmean + tittle
                              }
                          }
                          ## Expand it.
                          if (expandRt != 0) {
                              diffRt <- (pa[2] - pa[1]) * expandRt / 2
                              pa[1] <- pa[1] - diffRt
                              pa[2] <- pa[2] + diffRt
                          }
                          if (expandMz != 0) {
                              diffMz <- (pa[4] - pa[3]) * expandMz / 2
                              pa[3] <- pa[3] - diffMz
                              pa[4] <- pa[4] + diffMz
                          }
                          if (fixedMz != 0) {
                              pa[3] <- pa[3] - fixedMz
                              pa[4] <- pa[4] + fixedMz
                          }
                          if (fixedRt != 0) {
                              pa[1] <- pa[1] - fixedRt
                              pa[2] <- pa[2] + fixedRt
                          }
                          return(pa)
                      }
                  ))
              colnames(pkArea) <- c("rtmin", "rtmax", "mzmin", "mzmax")
              ## Add mzmed column - needed for MSW peak filling.
              pkArea <- cbind(group_idx = 1:nrow(pkArea), pkArea,
                              mzmed = as.numeric(fdef$mzmed))
              
              pkGrpVal <- featureValues(object)
              ## Check if there is anything to fill...
              if (!any(is.na(rowSums(pkGrpVal)))) {
                  message("No missing peaks present.")
                  return(object)
              }
              ## Split the object by file and define the peaks for which
              objectL <- vector("list", length(fileNames(object)))
              pkAreaL <- objectL
              for (i in 1:length(fileNames(object))) {
                  suppressMessages(
                      objectL[[i]] <- filterFile(object, file = i,
                                                 keepAdjustedRtime = TRUE)
                  )
                  ## Want to extract intensities only for peaks that were not
                  ## found in a sample.
                  pkAreaL[[i]] <- pkArea[is.na(pkGrpVal[, i]), , drop = FALSE]
              }
    
              ## Get to know what algorithm was used for the peak detection.
              ## Special cases are MSWParam (no retention time) and
              ## MatchedFilterParam (integrate from profile matrix).
              ph <- processHistory(object, type = .PROCSTEP.PEAK.DETECTION)
              findPeakMethod <- "unknown"
              mzCenterFun <- "wMean"
              if (length(ph)) {
                  if (is(ph[[1]], "XProcessHistory")) {
                      prm <- ph[[1]]@param
                      findPeakMethod <- .param2string(prm)
                      ## Check if the param class has a mzCenterFun slot
                      if (.hasSlot(prm, "mzCenterFun"))
                          mzCenterFun <- prm@mzCenterFun
                  }
              }
              ## Now rename that to the correct function name in xcms.
              mzCenterFun <- paste("mzCenter",
                                   gsub(mzCenterFun, pattern = "mzCenter.",
                                        replacement = "", fixed = TRUE), sep=".")
              if (findPeakMethod == "MSW") {
                  rts <- rtime(object, bySample = TRUE)
                  ## Ensure that we REALLY have direct injection data.
                  if (any(lengths(rts) > 1))
                      stop("The data is supposed to be direct injection data, ",
                           "but I got files with more than one spectrum/",
                           "retention time!")
                  ## That's not working, because integration uses the rt.
                  res <- bpmapply(FUN = .getMSWPeakData, objectL,
                                  pkAreaL, as.list(1:length(objectL)),
                                  MoreArgs = list(
                                      cn = colnames(chromPeaks(object))),
                                  BPPARAM = BPPARAM, SIMPLIFY = FALSE)
              } else if (findPeakMethod == "matchedFilter") {
                  res <- bpmapply(FUN = .getChromPeakData_matchedFilter,
                                  objectL, pkAreaL, as.list(1:length(objectL)),
                                  MoreArgs = list(
                                      cn = colnames(chromPeaks(object)),
                                      param = prm
                                  ),
                                  BPPARAM = BPPARAM, SIMPLIFY = FALSE)
              } else {
                  res <- bpmapply(FUN = .getChromPeakData, objectL,
                                  pkAreaL, as.list(1:length(objectL)),
                                  MoreArgs = list(
                                      cn = colnames(chromPeaks(object)),
                                      mzCenterFun = mzCenterFun),
                                  BPPARAM = BPPARAM, SIMPLIFY = FALSE)
              }
                            
              res <- do.call(rbind, res)
              if (any(colnames(res) == "is_filled"))
                  res[, "is_filled"] <- 1
              else
                  res <- cbind(res, is_filled = 1)
              ## cbind the group_idx column to track the feature/peak group.
              res <- cbind(res, group_idx = do.call(rbind, pkAreaL)[, "group_idx"])
              ## Remove those without a signal
              res <- res[!is.na(res[, "into"]), , drop = FALSE]
              if (nrow(res) == 0) {
                  warning("Could not integrate any signal for the missing ",
                          "peaks! Consider increasing 'expandMz' and 'expandRt'.")
                  return(object)
              }
              
              ## Get the msFeatureData:
              newFd <- new("MsFeatureData")
              newFd@.xData <- .copy_env(object@msFeatureData)
              incr <- nrow(chromPeaks(object))
              for (i in unique(res[, "group_idx"])) {
                  fdef$peakidx[[i]] <- c(fdef$peakidx[[i]],
                  (which(res[, "group_idx"] == i) + incr))
              }
              
              chromPeaks(newFd) <- rbind(chromPeaks(object), res[, -ncol(res)])
              featureDefinitions(newFd) <- fdef
              lockEnvironment(newFd, bindings = TRUE)
              object@msFeatureData <- newFd
              ## Add a process history step
              ph <- XProcessHistory(param = param,
                                    date. = startDate,
                                    type. = .PROCSTEP.PEAK.FILLING,
                                    fileIndex = 1:length(fileNames(object)))
              object <- addProcessHistory(object, ph) ## this also validates object.
              return(object)
          })

#' @rdname fillChromPeaks
setMethod(
    "fillChromPeaks",
    signature(object = "XCMSnExp", param = "missing"),
              function(object,
                       param,
                       BPPARAM = bpparam()) {
                  fillChromPeaks(object, param = FillChromPeaksParam(),
                                 BPPARAM = BPPARAM)
              })

#' @aliases dropFilledChromPeaks
#'
#' @description
#'
#' \code{dropFilledChromPeaks}: drops any filled-in chromatographic
#' peaks (filled in by the \code{\link{fillChromPeaks}} method) and all
#' related process history steps.
#'
#' @rdname XCMSnExp-class
#' 
#' @seealso \code{\link{fillChromPeaks}} for the method to fill-in eventually
#'     missing chromatographic peaks for a feature in some samples.
setMethod("dropFilledChromPeaks", "XCMSnExp", function(object) {
    if (!.hasFilledPeaks(object))
        return(object)
    keep_pks <- which(chromPeaks(object)[, "is_filled"] == 0)
    newFd <- new("MsFeatureData")
    newFd@.xData <- .copy_env(object@msFeatureData)
    ## Update index in featureDefinitions
    fd <- featureDefinitions(newFd)
    fd <- split(fd, 1:nrow(fd))
    fdL <- lapply(fd, function(z) {
        z$peakidx <- list(z$peakidx[[1]][z$peakidx[[1]] %in% keep_pks])
        return(z)
    })
    featureDefinitions(newFd) <- do.call(rbind, fdL)
    ## Remove peaks
    chromPeaks(newFd) <- chromPeaks(newFd)[keep_pks, , drop = FALSE]
    ## newFd <- .filterChromPeaks(object@msFeatureData, idx = keep_pks)
    object@msFeatureData <- newFd
    object <- dropProcessHistories(object, type = .PROCSTEP.PEAK.FILLING)
    if (validObject(object))
        return(object)
})

#' @aliases extractMsData
#'
#' @title DEPRECATED: Extract a `data.frame` containing MS data
#' 
#' @description
#'
#' **UPDATE**: the `extractMsData` and `plotMsData` functions are deprecated
#' and `as(x, "data.frame")` and `plot(x, type = "XIC")` (`x` being an
#' `OnDiskMSnExp` or `XCMSnExp` object) should be used instead. See examples
#' below. Be aware that filtering the raw object might however drop the
#' adjusted retention times. In such cases it is advisable to use the
#' [applyAdjustedRtime()] function prior to filtering.
#'
#' 
#' Extract a `data.frame` of retention time, mz and intensity
#' values from each file/sample in the provided rt-mz range (or for the full
#' data range if `rt` and `mz` are not defined).
#'
#' @param object A `XCMSnExp` or `OnDiskMSnExp` object.
#'
#' @param rt `numeric(2)` with the retention time range from which the
#'     data should be extracted.
#'
#' @param mz `numeric(2)` with the mz range.
#'
#' @param msLevel `integer` defining the MS level(s) to which the data
#'     should be sub-setted prior to extraction; defaults to
#'     `msLevel = 1L`.
#' 
#' @param adjustedRtime (for `extractMsData,XCMSnExp`): `logical(1)`
#'     specifying if adjusted or raw retention times should be reported.
#'     Defaults to adjusted retention times, if these are present in
#'     `object`.
#'
#' @return
#'
#' A `list` of length equal to the number of samples/files in
#' `object`. Each element being a `data.frame` with columns
#' `"rt"`, `"mz"` and `"i"` with the retention time, mz and
#' intensity tuples of a file. If no data is available for the mz-rt range
#' in a file a `data.frame` with 0 rows is returned for that file.
#'
#' @seealso `XCMSnExp` for the data object.
#'
#' @rdname extractMsData-method
#'
#' @author Johannes Rainer
#'
#' @md
#' 
#' @examples
#' ## Read some files from the test data package.
#' library(faahKO)
#' library(xcms)
#' library(magrittr)
#' fls <- dir(system.file("cdf/KO", package = "faahKO"), recursive = TRUE,
#'            full.names = TRUE)
#' raw_data <- readMSData(fls[1:2], mode = "onDisk")
#'
#' ## Extract the full data as a data.frame
#' ms_all <- as(raw_data, "data.frame")
#' head(ms_all)
#' nrow(ms_all)
#'
#' ## Read the full MS data for a defined mz-rt region.
#' res <- raw_data %>%
#'     filterRt(rt = c(2700, 2900)) %>%
#'     filterMz(mz = c(300, 320)) %>%
#'     as("data.frame")
#'
#' head(res)
#' nrow(res)
setMethod("extractMsData", "XCMSnExp",
          function(object, rt, mz, msLevel = 1L,
                   adjustedRtime = hasAdjustedRtime(object)){
              .Deprecated(msg = paste0("Use of 'extractMsData' is deprecated.",
                                       " Please use 'as(x, \"data.frame\")'"))
              ## Now, this method takes the adjusted rts, casts the object to
              ## an OnDiskMSnExp, eventually replaces the rtime in the
              ## featureData with the adjusted retention times (depending on
              ## adjustedRtime and calls the method for OnDiskMSnExp.
              if (adjustedRtime & hasAdjustedRtime(object)) {
                  fData(object)$retentionTime <- rtime(object, adjusted = TRUE)
              }
              object <- as(object, "OnDiskMSnExp")
              extractMsData(object, rt = rt, mz = mz, msLevel = msLevel)
          })


#' @rdname calibrate-calibrant-mass
#'
#' @param object An [XCMSnExp] object.
#' 
#' @param param The `CalibrantMassParam` object with the calibration settings.
#' 
#' @return
#'
#' The `calibrate` method returns an [XCMSnExp] object with the
#' chromatographic peaks being calibrated. Note that **only** the detected
#' peaks are calibrated, but not the individual mz values in each spectrum.
#' 
#' @md
setMethod("calibrate", "XCMSnExp", function(object, param) {
    if (missing(param)) {
        stop("Argument 'param' is missing")
    } else {
        if (!is(param, "CalibrantMassParam"))
            stop("The calibrate method for 'XCMSnExp' objects requires a ",
                 "'CalibrantMassParam' object to be passed with argument ",
                 "'param'")
    }
    if (isCalibrated(object))
        stop("'object' is already calibrated! Recurrent calibrations are not ",
             "supported")
    mzs <- .mz(param)
    n_samps <- length(fileNames(object))
    if (length(mzs) == 1)
        mzs <- replicate(mzs[[1]], n = n_samps, simplify = FALSE)
    if (length(mzs) != n_samps)
        stop("Number of calibrant mz vectors differs from the number of samples")
    startDate <- date()
    ## Dropping grouping results.
    object <- dropFeatureDefinitions(object)
    pks <- chromPeaks(object)
    adj_models <- vector("list", length = n_samps)
    method <- .method(param)
    for (i in 1:n_samps) {
        ## This could also be done with indices...
        pk_mz <- pks[pks[, "sample"] == i, c("mz", "into")]
        order_mz <- order(pk_mz[, 1])
        close_pks <- .matchpeaks2(pk_mz[order_mz, ], mzs[[i]],
                                  mzabs = .mzabs(param), mzppm = .mzppm(param),
                                  neighbours = .neighbors(param))
        if (nrow(close_pks) == 0) {
            warning("Sample ", i, ": can not calibrate as no peaks are close to",
                    "provided mz values")
            next
        } else {
            if (nrow(close_pks) == 1 & method != "shift") {
                warning("Sample ", i, ": only a single peak found, falling ",
                        "back to method = 'shift'")
                method <- "shift"
            }
        }
        
        prms <- estimate(close_pks, method)
        adj_models[[i]] <- prms
        a <- prms[1]                    # slope
        b <- prms[2]                    # intercept
        mz_ <- pk_mz[order_mz, 1]
        mz_min <- mz_[min(close_pks[, "pos"])]
        mz_max <- mz_[max(close_pks[, "pos"])]
        pk_mz[order_mz, "mz"] <- .calibrate_mz(mz_, method = method,
                                               minMz = mz_min, maxMz = mz_max,
                                               slope = a, intercept = b)
        pks[pks[, "sample"] == i, "mz"] <- pk_mz[, "mz"]
    }

    ## Set the new peak definitions. Careful to not drop additional stuff here.
    newFd <- new("MsFeatureData")
    newFd@.xData <- .copy_env(object@msFeatureData)
    chromPeaks(newFd) <- pks
    lockEnvironment(newFd, bindings = TRUE)
    object@msFeatureData <- newFd
    ## Add param to processHistory
    xph <- XProcessHistory(param = param, date. = startDate,
                           type. = .PROCSTEP.CALIBRATION,
                           fileIndex = 1:n_samps)
    object <- addProcessHistory(object, xph)
    if (validObject(object))
        object
})


#' @description
#'
#' \code{spectrapply} applies the provided function to each
#' \code{Spectrum} in the object and returns its
#' results. If no function is specified the function simply returns the
#' \code{list} of \code{Spectrum} objects.
#'
#' @param FUN For \code{spectrapply}: a function that should be applied to each
#'     spectrum in the object.
#' 
#' @rdname XCMSnExp-class
setMethod("spectrapply", "XCMSnExp", function(object, FUN = NULL,
                                              BPPARAM = bpparam(), ...) {
    ## replace raw with adjusted retention times!
    if (hasAdjustedRtime(object))
        fData(object)$retentionTime <- rtime(object, adjusted = TRUE)
    callNextMethod()
})

#' @description
#'
#' \code{split} splits an \code{XCMSnExp} object into a \code{list}
#' of \code{XCMSnExp} objects based on the provided parameter \code{f}.
#' Note that by default all pre-processing results are removed by the
#' splitting, except adjusted retention times, if the optional argument
#' \code{keepAdjustedRtime = TRUE} is provided.
#'
#' @param f For \code{split} a vector of length equal to the length of x
#'     defining how \code{x} will be splitted. It is converted internally to
#'     a \code{factor}.
#' 
#' @rdname XCMSnExp-filter-methods
setMethod("split", "XCMSnExp", function(x, f,
                                        drop = FALSE, ...) {
    if (drop)
        stop("'drop = TRUE' is not supported")
    if (missing(f))
        stop("required argument 'f' is missing")
    if (length(f) != length(x))
        stop("length of 'f' has to match the length of 'x'")
    keepAdjustedRtime <- list(...)$ke
    if (is.null(keepAdjustedRtime))
        keepAdjustedRtime <- FALSE
    if (!is.factor(f))
        f <- factor(f)
    res <- lapply(levels(f), function(z) {
        suppressWarnings(
            x[f == z, keepAdjustedRtime = keepAdjustedRtime]
        )
    })
    names(res) <- levels(f)
    res
})


#' @description
#'
#' \code{XCMSnExp} objects can be combined with the \code{c} function. This
#' combines identified chromatographic peaks and the objects' pheno data but
#' discards alignment results or feature definitions.
#' 
#' @rdname XCMSnExp-class
c.XCMSnExp <- function(...) {
    .concatenate_XCMSnExp(...)
}

#' @title Generate unique group (feature) names based on mass and retention time
#'
#' @description
#'
#' `groupnames` generates names for the identified features from the
#' correspondence analysis based in their mass and retention time. This
#' generates feature names that are equivalent to the group names of the *old*
#' user interface (aka xcms1).
#'
#' @param object `XCMSnExp` object containing correspondence results.
#'
#' @param mzdec `integer(1)` with the number of decimal places to use for m/z (
#'     defaults to `0`).
#'
#' @param rtdec `integer(1)` with the number of decimal places to use for the
#'     retention time (defaults to `0`).
#'
#' @param template `character` with existing group names whose format should
#'     be emulated.
#'
#' @return
#'
#' `character` with unique names for each feature in `object`. The
#' format is `M(m/z)T(time in seconds)`.
#'
#' @seealso [XCMSnExp].
#'
#' @md
#' 
#' @rdname groupnames-XCMSnExp
setMethod("groupnames", "XCMSnExp", function(object, mzdec = 0, rtdec = 0,
                                             template = NULL) {
    if (!hasFeatures(object))
        stop("No feature data present! Use 'groupChromPeaks' first")
    if (!missing(template)) {
        tempsplit <- strsplit(template[1], "[T_]")
        tempsplit <- strsplit(unlist(tempsplit), "\\.")
        if (length(tempsplit[[1]]) > 1)
            mzdec <- nchar(tempsplit[[1]][2])
        else
            mzdec <- 0
        if (length(tempsplit[[2]]) > 1)
            rtdec <- nchar(tempsplit[[2]][2])
        else
            rtdec <- 0
    }
    mzfmt <- paste0("%.", mzdec, "f")
    rtfmt <- paste0("%.", rtdec, "f")
    gnames <- paste0("M", sprintf(mzfmt, featureDefinitions(object)$mzmed),
                     "T", sprintf(rtfmt, featureDefinitions(object)$rtmed))
    if (any(dup <- duplicated(gnames)))
        for (dupname in unique(gnames[dup])) {
            dupidx <- which(gnames == dupname)
            gnames[dupidx] <- paste(gnames[dupidx], seq(along = dupidx),
                                    sep = "_")
        }
    gnames
})
