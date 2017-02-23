## Methods for the XCMSnExp object representing untargeted metabolomics
## results
#' @include functions-XCMSnExp.R do_groupChromPeaks-functions.R
#' do_adjustRtime-functions.R methods-xcmsRaw.R functions-OnDiskMSnExp.R

setMethod("initialize", "XCMSnExp", function(.Object, ...) {
    classVersion(.Object)["XCMSnExp"] <- "0.0.1"
    .Object <- callNextMethod(.Object, ...)
    lockEnvironment(.Object@msFeatureData)
    return(.Object)
})

##' @rdname XCMSnExp-class
setMethod("show", "XCMSnExp", function(object) {
    callNextMethod()
    ## And not XCMSnExp related stuff.
    cat("- - - xcms preprocessing - - -\n")
    if (hasChromPeaks(object)) {
        cat("Chromatographic peak detection:\n")
        cat(" ", nrow(chromPeaks(object)), " peaks identified in ",
            length(fileNames(object)), " samples.\n", sep = "")
        cat(" On average ",
            format(mean(table(chromPeaks(object)[, "sample"])), digits = 3),
            " chromatographic peaks per sample.\n", sep = "")
    }
    if (hasFeatures(object)) {
        cat("Correspondence:\n")
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
    }
    if (hasAdjustedRtime(object)) {
    }
})

##' @aliases hasAdjustedRtime
##'
##' @description \code{hasAdjustedRtime}: whether the object provides adjusted
##' retention times.
##'
##' @rdname XCMSnExp-class
setMethod("hasAdjustedRtime", "XCMSnExp", function(object) {
    return(hasAdjustedRtime(object@msFeatureData))
})

##' @aliases hasFeatures
##'
##' @description \code{hasFeatures}: whether the object contains correspondence
##' results (i.e. features).
##'
##' @rdname XCMSnExp-class
setMethod("hasFeatures", "XCMSnExp", function(object) {
    return(hasFeatures(object@msFeatureData))
})

##' @aliases hasChromPeaks
##'
##' @description \code{hasChromPeaks}: whether the object contains peak
##' detection results.
##'
##' @rdname XCMSnExp-class
setMethod("hasChromPeaks", "XCMSnExp", function(object) {
    return(hasChromPeaks(object@msFeatureData))
})

##' @aliases adjustedRtime
##'
##' @description \code{adjustedRtime},\code{adjustedRtime<-}:
##' extract/set adjusted retention times. \code{adjustedRtime<-} should not be
##' called manually, it is called internally by the \code{\link{adjustRtime}}
##' methods. For \code{XCMSnExp} objects, \code{adjustedRtime<-} does also apply
##' the retention time adjustment to the chromatographic peaks in the object.
##' The \code{bySample} parameter allows to specify whether the adjusted
##' retention time should be grouped by sample (file).
##'
##' @return For \code{adjustedRtime}: if \code{bySample = FALSE} a \code{numeric}
##' vector with the adjusted retention for each spectrum of all files/samples
##' within the object. If \code{bySample = TRUE } a \code{list} (length equal to
##' the number of samples) with adjusted retention times grouped by sample.
##' Returns \code{NULL} if no adjusted retention times are present.
##'
##' @rdname XCMSnExp-class
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
##' @aliases adjustedRtime<-
##'
##' @rdname XCMSnExp-class
setReplaceMethod("adjustedRtime", "XCMSnExp", function(object, value) {
    if (!is.list(value))
        stop("'value' is supposed to be a list of retention time values!")
    if (hasAdjustedRtime(object))
        object <- dropAdjustedRtime(object)
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

##' @aliases featureDefinitions
##'
##' @description \code{featureDefinitions}, \code{featureDefinitions<-}: extract
##' or set the correspondence results, i.e. the mz-rt features (peak groups).
##'
##' @return For \code{featureDefinitions}: a \code{DataFrame} with peak grouping
##' information, each row corresponding to one mz-rt feature (grouped peaks
##' within and across samples) and columns \code{"mzmed"} (median mz value),
##' \code{"mzmin"} (minimal mz value), \code{"mzmax"} (maximum mz value),
##' \code{"rtmed"} (median retention time), \code{"rtmin"} (minimal retention
##' time), \code{"rtmax"} (maximal retention time) and \code{"peakidx"}.
##' Column \code{"peakidx"} contains a \code{list} with indices of
##' chromatographic peaks (rows) in the matrix returned by the \code{chromPeaks}
##' method that belong to that feature group. The method returns \code{NULL} if
##' no feature definitions are present.
##'
##' @rdname XCMSnExp-class
setMethod("featureDefinitions", "XCMSnExp", function(object) {
    return(featureDefinitions(object@msFeatureData))
})
##' @aliases featureDefinitions<-
##'
##' @rdname XCMSnExp-class
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

##' @aliases chromPeaks
##'
##' @description \code{chromPeaks}, \code{chromPeaks<-}: extract or set
##' the matrix containing the information on identified chromatographic peaks.
##' Parameter \code{bySample} allows to specify whether peaks should be
##' returned ungrouped (default \code{bySample = FALSE}) or grouped by sample (
##' \code{bySample = TRUE}). The \code{chromPeaks<-} method for \code{XCMSnExp}
##' objects removes also all correspondence (peak grouping) and retention time
##' correction (alignment) results.
##' See description of the return value for details on the returned matrix. Users
##' usually don't have to use the \code{chromPeaks<-} method directly as detected
##' chromatographic peaks are added to the object by the
##' \code{\link{findChromPeaks}} method.
##'
##' @return For \code{chromPeaks}: if \code{bySample = FALSE} a \code{matrix} with
##' at least the following columns: \code{"mz"} (mz value for the largest
##' intensity), \code{"mzmin"} (minimal mz value), \code{"mzmax"} (maximal mz
##' value), \code{"rt"} (retention time for the peak apex), \code{"rtmin"}
##' (minimal retention time), \code{"rtmax"} (maximal retention time),
##' \code{"into"} (integrated, original, intensity of the peak) and
##' \code{"sample"} (sample index in which the peak was identified).
##' Depending on the employed peak detection algorithm and the
##' \code{verboseColumns} parameter of it additional columns might be returned.
##' For \code{bySample = TRUE} the chronatographic peaks are returned as a
##' \code{list} of matrices, each containing the chromatographic peak of a
##' specific sample. For sample in which no feastures were detected a matrix
##' with 0 rows is returned.
##'
##' @rdname XCMSnExp-class
setMethod("chromPeaks", "XCMSnExp", function(object, bySample = FALSE) {
    if (bySample) {
        tmp <- split(chromPeaks(object), f = chromPeaks(object)[, "sample"])
        ## Ensure we return something for each sample in case there is a sample
        ## without detected peaks.
        res <- vector("list", length(fileNames(object)))
        names(res) <- as.character(1:length(res))
        tmp <- split.data.frame(chromPeaks(object),
                                f = chromPeaks(object)[, "sample"])
        res[as.numeric(names(tmp))] <- tmp
        if (any(lengths(res) == 0)) {
            emat <- matrix(nrow = 0, ncol = ncol(tmp[[1]]))
            colnamers(emat) <- colnames(tmp[[1]])
            res[lengths(res) == 0] <- emat
        }
        return(res)
    } else {
        return(chromPeaks(object@msFeatureData))
    }
})
##' @aliases chromPeaks<-
##'
##' @rdname XCMSnExp-class
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

##' @description \code{rtime}: extracts the retention time for each
##' scan. The \code{bySample} parameter allows to return the values grouped
##' by sample/file and \code{adjusted} whether adjusted or raw retention times
##' should be returned. By default the method returns adjusted retention times,
##' if they are available (i.e. if retention times were adjusted using the
##' \code{\link{adjustRtime}} method).
##'
##' @param bySample logical(1) specifying whether results should be grouped by
##' sample.
##'
##' @param adjusted logical(1) whether adjusted or raw (i.e. the original
##' retention times reported in the files) should be returned.
##' 
##' @return For \code{rtime}: if \code{bySample = FALSE} a numeric vector with the
##' retention times of each scan, if \code{bySample = TRUE} a \code{list} of
##' numeric vectors with the retention times per sample.
##'
##' @rdname XCMSnExp-class
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

##' @description \code{mz}: extracts the mz values from each scan of
##' all files within an \code{XCMSnExp} object. These values are extracted from
##' the original data files and eventual processing steps are applied
##' \emph{on the fly}. Using the \code{bySample} parameter it is possible to
##' switch from the default grouping of mz values by spectrum/scan to a grouping
##' by sample/file.
##'
##' @return For \code{mz}: if \code{bySample = FALSE} a \code{list} with the mz
##' values (numeric vectors) of each scan. If \code{bySample = TRUE} a
##' \code{list} with the mz values per sample.
##'
##' @rdname XCMSnExp-class
setMethod("mz", "XCMSnExp", function(object, bySample = FALSE) {
    res <- callNextMethod(object = object)
    if (bySample) {
        tmp <- lapply(split(res, fromFile(object)), unlist, use.names = FALSE)
        res <- vector("list", length(fileNames(object)))
        names(res) <- as.character(1:length(res))
        res[as.numeric(names(tmp))] <- tmp
    }
    return(res)
})

##' @description \code{intensity}: extracts the intensity values from
##' each scan of all files within an \code{XCMSnExp} object. These values are
##' extracted from the original data files and eventual processing steps are
##' applied \emph{on the fly}. Using the \code{bySample} parameter it is possible
##' to switch from the default grouping of intensity values by spectrum/scan to
##' a grouping by sample/file.
##'
##' @return For \code{intensity}: if \code{bySample = FALSE} a \code{list} with
##' the intensity values (numeric vectors) of each scan. If
##' \code{bySample = TRUE} a \code{list} with the intensity values per sample.
##'
##' @rdname XCMSnExp-class
setMethod("intensity", "XCMSnExp", function(object, bySample = FALSE) {
    res <- callNextMethod(object = object)
    if (bySample) {
        tmp <- lapply(split(res, fromFile(object)), unlist, use.names = FALSE)
        res <- vector("list", length(fileNames(object)))
        names(res) <- as.character(1:length(res))
        res[as.numeric(names(tmp))] <- tmp
    }
    return(res)
})

##' @description \code{spectra}: extracts the
##' \code{\link[MSnbase]{Spectrum}} objects containing all data from
##' \code{object}. The values are extracted from the original data files and
##' eventual processing steps are applied \emph{on the fly}. By setting
##' \code{bySample = TRUE}, the spectra are returned grouped by sample/file. If
##' the \code{XCMSnExp} object contains adjusted retention times, these are
##' returned by default in the \code{Spectrum} objects (can be overwritten
##' by setting \code{adjusted = FALSE}).
##'
##' @return For \code{spectra}: if \code{bySample = FALSE} a \code{list} with
##' \code{\link[MSnbase]{Spectrum}} objects. If \code{bySample = TRUE} the result
##' is grouped by sample, i.e. as a \code{list} of \code{lists}, each element in
##' the \emph{outer} \code{list} being the \code{list} of spectra of the specific
##' file.
##'
##' @rdname XCMSnExp-class
setMethod("spectra", "XCMSnExp", function(object, bySample = FALSE,
                                          adjusted = hasAdjustedRtime(object)) {
    res <- callNextMethod(object = object)
    ## replace the rtime of these with the adjusted ones - if present.
    if (adjusted & hasAdjustedRtime(object)) {
        rts <- adjustedRtime(object)
        res <- mapply(FUN = function(a, b) {
            a@rt <- b
            return(a)
        }, a = res, b = rts)
    }
    if (bySample) {
        tmp <- split(res, fromFile(object))
        ## That's to ensure that we're always returning something for all files.
        res <- vector("list", length(fileNames(object)))
        names(res) <- as.character(1:length(res))
        res[as.numeric(names(tmp))] <- tmp
    }
    return(res)
})

## processHistory
##' @aliases processHistory
##' @description \code{processHistory}: returns a \code{list} with
##' \code{\link{ProcessHistory}} objects (or objects inheriting from this base
##' class) representing the individual processing steps that have been performed,
##' eventually along with their settings (\code{Param} parameter class). Optional
##' arguments \code{fileIndex} and \code{type} allow to restrict to process steps
##' of a certain type or performed on a certain file.
##'
##' @param fileIndex For \code{processHistory}: optional \code{numeric}
##' specifying the index of the files/samples for which the
##' \code{\link{ProcessHistory}} objects should be retrieved.
##'
##' @param type For \code{processHistory}: restrict returned
##' \code{\link{ProcessHistory}} objects to analysis steps of a certain type.
##' Supported values are \code{"Unknown"}, \code{"Peak detection"},
##' \code{"Peak grouping"} and \code{"Retention time correction"}.
##'
##' @return For \code{processHistory}: a \code{list} of
##' \code{\link{ProcessHistory}} objects providing the details of the individual
##' data processing steps that have been performed.
##'
##' @rdname XCMSnExp-class
setMethod("processHistory", "XCMSnExp", function(object, fileIndex, type) {
    ph <- object@.processHistory
    if (length(ph)) {
        if (!missing(fileIndex)) {
            if (!all(fileIndex %in% 1:length(fileNames(object))))
                stop("'fileIndex' has to be within 1 and the number of samples!")
            gotIt <- unlist(lapply(ph, function(z) {
                return(any(z@fileIndex %in% fileIndex))
            }))
            if (!any(gotIt))
                return(list())
            ph <- ph[gotIt]
        }
        if (!missing(type) & length(ph)) {
            gotIt <- unlist(lapply(ph, function(z) {
                return(any(type == processType(z)))
            }))
            if (!any(gotIt))
                return(list())
            ph <- ph[gotIt]
        }
        return(ph)
    } else {
        return(list())
    }
})

##' @description \code{addProcessHistory}: adds (appends) a single
##' \code{\link{ProcessHistory}} object to the \code{.processHistory} slot.
##'
##' @return The \code{addProcessHistory} method returns the input object with the
##' provided \code{\link{ProcessHistory}} appended to the process history.
##' @noRd
setMethod("addProcessHistory", "XCMSnExp", function(object, ph) {
    if (!inherits(ph, "ProcessHistory"))
        stop("Argument 'ph' has to be of type 'ProcessHistory' or a class ",
             "extending it!")
    object@.processHistory[[(length(object@.processHistory) + 1)]] <- ph
    if (validObject(object))
        return(object)
})

##' @aliases dropChromPeaks
##'
##' @description \code{dropChromPeaks}: drops any identified chromatographic
##' peaks and returns the object without that information. Note that for
##' \code{XCMSnExp} objects the method drops all results from a correspondence
##' (peak grouping) or alignment (retention time adjustment) too.
##' For \code{XCMSnExp} objects the method drops also any related process
##' history steps.
##'
##' @rdname XCMSnExp-class
setMethod("dropChromPeaks", "XCMSnExp", function(object) {
    if (hasChromPeaks(object)) {
        object <- dropProcessHistories(object, type = .PROCSTEP.PEAK.DETECTION)
        ## Make sure we delete all related process history steps
        object <- dropProcessHistories(object, type = .PROCSTEP.RTIME.CORRECTION)
        object <- dropProcessHistories(object, type = .PROCSTEP.PEAK.GROUPING)
        ## idx_fd <- which(unlist(lapply(processHistory(object), processType)) ==
        ##                 .PROCSTEP.PEAK.DETECTION)
        ## if (length(idx_fd) > 0)
        ##     object@.processHistory <- object@.processHistory[-idx_fd]
        newFd <- new("MsFeatureData")
        newFd@.xData <- .copy_env(object@msFeatureData)
        newFd <- dropChromPeaks(newFd)
        ## Dropping other results from the environment (not the object).
        if (hasAdjustedRtime(newFd))
            newFd <- dropAdjustedRtime(newFd)
        if (hasFeatures(newFd))
            newFd <- dropFeatureDefinitions(newFd)
        lockEnvironment(newFd, bindings = TRUE)
        object@msFeatureData <- newFd
    }
    if (validObject(object))
        return(object)
})
##' @aliases dropFeatureDefinitions
##'
##' @description \code{dropFeatureDefinitions}: drops the results from a
##' correspondence (peak grouping) analysis, i.e. the definition of the mz-rt
##' features and returns the object without that information. Note that for
##' \code{XCMSnExp} objects the method will also drop retention time adjustment
##' results, if these were performed after the last peak grouping (i.e. which
##' base on the results from the peak grouping that are going to be removed).
##' For \code{XCMSnExp} objects also all related process history steps are
##' removed.
##'
##' @param keepAdjRtime For \code{dropFeatureDefinitions,XCMSnExp}:
##' \code{logical(1)} defining whether eventually present retention time
##' adjustment should not be dropped. By default dropping feature definitions
##' drops retention time adjustment results too.
##'
##' @param dropLastN For \code{dropFeatureDefinitions,XCMSnExp}:
##' \code{numeric(1)} defining the number of peak grouping related process
##' history steps to remove. By default \code{dropLastN = -1}, dropping the
##' chromatographic peaks removes all process history steps related to peak
##' grouping. Setting e.g. \code{dropLastN = 1} will only remove the most recent
##' peak grouping related process history step.
##' 
##' @rdname XCMSnExp-class
setMethod("dropFeatureDefinitions", "XCMSnExp", function(object,
                                                         keepAdjRtime = FALSE,
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
        newFd <- new("MsFeatureData")
        newFd@.xData <- .copy_env(object@msFeatureData)
        newFd <- dropFeatureDefinitions(newFd)
        lockEnvironment(newFd, bindings = TRUE)
        object@msFeatureData <- newFd
        ## 2) If retention time correction was performed after the latest peak
        ##    alignment, drop also the retention time correction and all related
        ##    process history steps.
        ##    Otherwise (grouping performed after retention time adjustment) do
        ##    nothing - this keeps eventual alignment related process history
        ##    steps performed before retention time correction.
        if (hasAdjustedRtime(object)) {
            if (max(idx_art) > max(idx_fal)) {
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

##' @aliases dropAdjustedRtime
##'
##' @description \code{dropAdjustedRtime}: drops any retention time
##' adjustment information and returns the object without adjusted retention
##' time. For \code{XCMSnExp} object this also reverts the retention times
##' reported for the chromatographic peaks in the peak matrix to the original,
##' raw, ones (after chromatographic peak detection). Note that for
##' \code{XCMSnExp} objects the method drops also all peak grouping results
##' if these were performed \emph{after} the retention time adjustment.
##' For \code{XCMSnExp} objects the method drops also any related process history
##' steps.
##'
##' @rdname XCMSnExp-class
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


##' @title XCMSnExp data manipulation methods inherited from MSnbase
##'
##' @description The methods listed on this page are \code{\link{XCMSnExp}}
##' methods inherited from its parent, the \code{\link[MSnbase]{OnDiskMSnExp}}
##' class from the \code{MSnbase} package, that alter the raw data or are related
##' to data subsetting. Thus calling any of these methods causes all \code{xcms}
##' pre-processing results to be removed from the \code{\link{XCMSnExp}} object
##' to ensure its data integrity.
##'
##' The \code{[} method allows to subset a \code{\link{XCMSnExp}} object by
##' spectra. For more details and examples see the documentation for
##' \code{\link[MSnbase]{OnDiskMSnExp}}.
##'
##' @param x For \code{[}: an \code{\link{XCMSnExp}} object.
##'
##' @param i For \code{[}: \code{numeric} or \code{logical} vector specifying to
##' which spectra the data set should be reduced.
##'
##' @param j For \code{[}: not supported.
##'
##' @param drop For \code{[}: not supported.
##'
##' @return For all methods: a \code{XCMSnExp} object.
##'
##' @rdname XCMSnExp-inherited-methods
##'
##' @seealso \code{\link{XCMSnExp-filter}} for methods to filter and subset
##' \code{XCMSnExp} objects.
##' @seealso \code{\link{XCMSnExp}} for base class documentation.
##' @seealso \code{\link[MSnbase]{OnDiskMSnExp}} for the documentation of the
##' parent class.
##'
##' @author Johannes Rainer
setMethod("[", signature(x = "XCMSnExp", i = "logicalOrNumeric", j = "missing",
                         drop = "missing"),
          function(x, i, j, drop) {
              ## Want to support subsetting of the peaks!
              ## This means that we will also have to adjust the process
              ## history accordingly.
              if (hasAdjustedRtime(x) | hasFeatures(x) |
                  hasChromPeaks(x)) {
                  ## x@.processHistory <- list()
                  ## x@msFeatureData <- new("MsFeatureData")
                  suppressMessages(
                      x <- dropAdjustedRtime(x)
                  )
                  suppressMessages(
                      x <- dropFeatureDefinitions(x)
                  )
                  suppressMessages(
                      x <- dropChromPeaks(x)
                  )
                  warning("Removed preprocessing results")
              }
              callNextMethod()
          })

## setMethod("splitByFile", c("XCMSnExp", "factor"), function(x, f) {
##     if (length(f) != length(fileNames(x)))
##         stop("length of 'f' has to match the length of samples/files in 'object'.")
##     idxs <- lapply(levels(f), function(z) which(f == z))
##     ## Now I can run a filterFile on these.
##     res <- lapply(idxs, function(z) {
##         return(filterFile(x, file = z))
##     })
##     names(res) <- levels(f)
##     return(res)
## })

##' @description \code{bin}: allows to \emph{bin} spectra. See
##' \code{\link[MSnbase]{bin}} documentation for more details and examples.
##'
##' @param object An \code{\link{XCMSnExp}} or \code{OnDiskMSnExp} object.
##'
##' @param binSize \code{numeric(1)} defining the size of a bin (in Dalton).
##'
##' @param msLevel. For \code{bin}, \code{clean}, \code{filterMsLevel},
##' \code{removePeaks}:  \code{numeric(1)} defining the MS level(s)
##' to which operations should be applied or to which the object should be
##' subsetted.
##'
##' @rdname XCMSnExp-inherited-methods
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

##' @description \code{clean}: removes unused \code{0} intensity data
##' points. See \code{\link[MSnbase]{clean}} documentation for details and
##' examples.
##'
##' @param all For \code{clean}: \code{logical(1)}, if \code{TRUE} all zeros are
##' removed.
##'
##' @param verbose \code{logical(1)} whether progress information should be
##' displayed.
##'
##' @rdname XCMSnExp-inherited-methods
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

##' @description \code{filterMsLevel}: reduces the \code{\link{XCMSnExp}}
##' object to spectra of the specified MS level(s). See
##' \code{\link[MSnbase]{filterMsLevel}} documentation for details and examples.
##'
##' @rdname XCMSnExp-inherited-methods
setMethod("filterMsLevel", "XCMSnExp", function(object, msLevel.) {
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

##' @description \code{filterAcquisitionNum}: filters the
##' \code{\link{XCMSnExp}} object keeping only spectra with the provided
##' acquisition numbers. See \code{\link[MSnbase]{filterAcquisitionNum}} for
##' details and examples.
##'
##' @param n For \code{filterAcquisitionNum}: \code{integer} defining the
##' acquisition numbers of the spectra to which the data set should be
##' sub-setted.
##'
##' @param file For \code{filterAcquisitionNum}:
##' \code{integer} defining the file index within the object to subset the
##' object by file.
##'
##' @rdname XCMSnExp-inherited-methods
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

##' @aliases XCMSnExp-filter
##' @title XCMSnExp filtering and subsetting
##'
##' The methods listed on this page allow to filter and subset
##' \code{\link{XCMSnExp}} objects. Most of them are inherited from the
##' \code{\link[MSnbase]{OnDiskMSnExp}} object and have been adapted for
##' \code{\link{XCMSnExp}} to enable subsetting also on the preprocessing
##' results.
##'
##' @description \code{filterFile}: allows to reduce the
##' \code{\link{XCMSnExp}} to data from only certain files. Identified
##' chromatographic peaks for these files are retained while all eventually
##' present features (peak grouping information) are dropped. By default also
##' adjusted retention times are removed. This can be overwritten by setting
##' \code{keepAdjustedRtime = TRUE}, but users should use this option with
##' caution.
##'
##' @note The \code{filterFile} method removes also process history steps not
##' related to the files to which the object should be sub-setted and updates
##' the \code{fileIndex} attribute accordingly. Also, the method does not allow
##' arbitrary ordering of the files or re-ordering of the files within the
##' object.
##'
##' @param object A \code{\link{XCMSnExp}} object.
##'
##' @param file For \code{filterFile}: \code{integer} defining the file index
##' within the object to subset the object by file or \code{character} specifying
##' the file names to sub set. The indices are expected to be increasingly
##' ordered, if not they are ordered internally.
##'
##' @param keepAdjustedRtime For \code{filterFile}: \code{logical(1)} defining
##' whether the adjusted retention times should be kept, even if features are
##' being removed (and the retention time correction being potentially performed
##' on these features).
##' 
##' @return All methods return an \code{\link{XCMSnExp}} object.
##'
##' @author Johannes Rainer
##'
##' @seealso \code{\link{XCMSnExp}} for base class documentation.
##'
##' @rdname XCMSnExp-filter-methods
##' @examples
##'
##' ## Load some of the files from the faahKO package.
##' library(faahKO)
##' fs <- c(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
##'         system.file('cdf/KO/ko16.CDF', package = "faahKO"),
##'         system.file('cdf/KO/ko18.CDF', package = "faahKO"))
##' ## Read the files
##' od <- readMSData2(fs)
##'
##' ## Perform peak detection on them using default matched filter settings.
##' mfp <- MatchedFilterParam()
##' xod <- findChromPeaks(od, param = mfp)
##'
##' ## Subset the dataset to the first and third file.
##' xod_sub <- filterFile(xod, file = c(1, 3))
##'
##' ## The number of chromatographic peaks per file for the full object
##' table(chromPeaks(xod)[, "sample"])
##'
##' ## The number of chromatographic peaks per file for the subset
##' table(chromPeaks(xod_sub)[, "sample"])
##'
##' basename(fileNames(xod))
##' basename(fileNames(xod_sub))
##'
##' ## Filter on mz values; chromatographic peaks and features within the
##' ## mz range are retained (as well as adjusted retention times).
##' xod_sub <- filterMz(xod, mz = c(300, 400))
##' head(chromPeaks(xod_sub))
##' nrow(chromPeaks(xod_sub))
##' nrow(chromPeaks(xod))
##'
##' ## Filter on rt values. All chromatographic peaks and features within the
##' ## retention time range are retained. Filtering is performed by default on
##' ## adjusted retention times, if present.
##' xod_sub <- filterRt(xod, rt = c(2700, 2900))
##'
##' range(rtime(xod_sub))
##' head(chromPeaks(xod_sub))
##' range(chromPeaks(xod_sub)[, "rt"])
##'
##' nrow(chromPeaks(xod))
##' nrow(chromPeaks(xod_sub))
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
    ## Dropping data.
    if (hasAdjustedRtime(object) & !keepAdjustedRtime) {
        message("Adjusted retention times removed.")
        object <- dropAdjustedRtime(object)
    }
    suppressWarnings(
        adjRt <- adjustedRtime(object, bySample = TRUE)
    )
    ## Get the data we want to keep/subset.
    fts <- chromPeaks(object)
    ## Keep also the processHistory
    ph <- processHistory(object)
    if (hasFeatures(object)) {
        message("Correspondence results (features) removed.")
        suppressMessages(
            object <- dropFeatureDefinitions(object)
        )
    }
    ph <- dropProcessHistoriesList(ph, type = .PROCSTEP.PEAK.GROUPING)
    ## ## Process the processing history.
    ## object <- dropProcessHistories(object, type = c(.PROCSTEP.PEAK.GROUPING,
    ##                                                 .PROCSTEP.RTIME.CORRECTION))
    ## ph <- processHistory(object)
    ## The next method will actually clean everything, process history and
    ## msFeatureData
    suppressWarnings(
        object <- callNextMethod(object = object, file = file)
    )
    ## Remove ProcessHistory not related to any of the files.
    if (length(ph)) {
        kp <- unlist(lapply(ph, function(z) {
            any(fileIndex(z) %in% file)
        }))
        ph <- ph[kp]
    }
    if (length(ph)) {
        ph <- lapply(ph, function(z) {
            updateFileIndex(z, old = file, new = 1:length(file))
        })
    }
    ## Process peaks.
    fts <- fts[fts[, "sample"] %in% file, , drop = FALSE]
    fts[, "sample"] <- match(fts[, "sample"], file)
    if (length(adjRt) & keepAdjustedRtime) {
        ## Put all directly into the msFeatureData environment to avoid an
        ## additional correction of the peak retention times by the adjusted rt.
        newFd <- new("MsFeatureData")
        newFd@.xData <- .copy_env(object@msFeatureData)
        chromPeaks(newFd) <- fts
        adjustedRtime(newFd) <- adjRt[file]
        lockEnvironment(newFd, bindings = TRUE)
        object@msFeatureData <- newFd
    } else {
        chromPeaks(object) <- fts
    }
    object@.processHistory <- ph
    return(object)
})

##' @description \code{filterMz}: filters the data set based on the
##' provided mz value range. All chromatographic peaks and features (grouped
##' peaks) falling completely within the provided mz value range are retained
##' (if their minimal mz value is \code{>= mz[1]} and the maximal mz value
##' \code{<= mz[2]}. Adjusted retention times, if present, are not altered by
##' the filtering.
##'
##' @param mz For \code{filterMz}: \code{numeric(2)} defining the lower and upper
##' mz value for the filtering.
##'
##' @param msLevel. For \code{filterMz}, \code{filterRt}, \code{numeric(1)}
##' defining the MS level(s) to which operations should be applied or to which
##' the object should be subsetted.
##'
##' @param ... Optional additional arguments.
##'
##' @rdname XCMSnExp-filter-methods
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

##' @description \code{filterRt}: filters the data set based on the
##' provided retention time range. All chromatographic peaks and features
##' (grouped peaks) the specified retention time window are retained (i.e. if
##' the retention time corresponding to the peak's apex is within the specified
##' rt range). If retention time correction has been performed, the method will
##' by default filter the object by adjusted retention times.
##' The argument \code{adjusted} allows to specify manually whether filtering
##' should be performed by raw or adjusted retention times. Filtering by
##' retention time does not drop any preprocessing results.
##' The method returns an empty object if no spectrum or feature is within the
##' specified retention time range.
##'
##' @param rt For \code{filterRt}: \code{numeric(2)} defining the retention time
##' window (lower and upper bound) for the filtering.
##'
##' @param adjusted For \code{filterRt}: \code{logical} indicating whether the
##' object should be filtered by original (\code{adjusted = FALSE}) or adjusted
##' retention times (\code{adjusted = TRUE}).
##' For \code{spectra}: whether the retention times in the individual
##' \code{Spectrum} objects should be the adjusted or raw retention times.
##'
##' @rdname XCMSnExp-filter-methods
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
                                                    .PROCSTEP.PEAK.GROUPING))
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
    ## suppressWarnings(
    ## Specifically call the [ from the OnDiskMSnExp!
    ## Otherwise we unnecessarily have to drop stuff which has a negative
    ## impact on performance.
    ## theM <- getMethod("[", signature = c(x = "OnDiskMSnExp",
    ##                                      i = "logicalOrNumeric",
    ##                                      j = "missing",
    ##                                      drop = "missing"))
    ## object <- theM(x = object, i = base::which(keep_logical))
    ## )
    ## Fix for issue #124
    ## Now, this casting is not ideal - have to find an easier way to call the
    ## subset method from OnDiskMSnExp...
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


##' The \code{normalize} method performs basic normalization of spectra
##' intensities. See \code{\link[MSnbase]{normalize}} documentation for details
##' and examples.
##'
##' @param method For \code{normalize}: \code{character(1)} specifying the
##' normalization method. See \code{\link[MSnbase]{normalize}} for details.
##' For \code{pickPeaks}: \code{character(1)} defining the method. See
##' \code{\link[MSnbase]{pickPeaks}} for options. For \code{smooth}:
##' \code{character(1)} defining the method. See \code{\link[MSnbase]{smooth}}
##' for options and details.
##'
##' @rdname XCMSnExp-inherited-methods
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

##' The \code{pickPeaks} method performs peak picking. See
##' \code{\link[MSnbase]{pickPeaks}} documentation for details and examples.
##'
##' @param halfWindowSize For \code{pickPeaks} and \code{smooth}:
##' \code{integer(1)} defining the window size for the peak picking. See
##' \code{\link[MSnbase]{pickPeaks}} and \code{\link[MSnbase]{smooth}} for
##' details and options.
##'
##' @param SNR For \code{pickPeaks}: \code{numeric(1)} defining the signal to
##' noise ratio to be considered. See \code{\link[MSnbase]{pickPeaks}}
##' documentation for details.
##'
##' @param ... Optional additional arguments.
##'
##' @rdname XCMSnExp-inherited-methods
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

##' The \code{removePeaks} method removes mass peaks (intensities) lower than a
##' threshold. Note that these peaks refer to \emph{mass} peaks, which are
##' different to the chromatographic peaks detected and analyzed in a
##' metabolomics experiment! See \code{\link[MSnbase]{removePeaks}}
##' documentation for details and examples.
##'
##' @param t For \code{removePeaks}: either a \code{numeric(1)} or \code{"min"}
##' defining the threshold (method) to be used. See
##' \code{\link[MSnbase]{removePeaks}} for details.
##'
##' @rdname XCMSnExp-inherited-methods
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

##' The \code{smooth} method smooths spectra. See \code{\link[MSnbase]{smooth}}
##' documentation for details and examples.
##'
##' @rdname XCMSnExp-inherited-methods
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

## @param from For \code{setAs} and \code{as}: an \code{XCMSnExp} object.
## @param to For \code{setAs} and \code{as}: \code{"xcmsSet"}
##' @title Data container storing xcms preprocessing results
##'
##' @aliases setAs
##' @rdname XCMSnExp-class
##' @name XCMSnExp-class
setAs(from = "XCMSnExp", to = "xcmsSet", def = .XCMSnExp2xcmsSet)


##' @title Peak grouping/correspondence based on time dimension peak densities
##'
##' @description \code{groupChromPeaks,XCMSnExp,PeakDensityParam}:
##' performs correspondence (peak grouping within and across samples) within in
##' mz dimension overlapping slices of MS data based on the density distribution
##' of the identified chromatographic peaks in the slice along the time axis.
##'
##' @note Calling \code{groupChromPeaks} on an \code{XCMSnExp} object will cause
##' all eventually present previous correspondence results to be dropped.
##'
##' @param object For \code{groupChromPeaks}: an \code{\link{XCMSnExp}} object
##' containing the results from a previous peak detection analysis (see
##' \code{\link{findChromPeaks}}).
##'
##' For all other methods: a \code{PeakDensityParam} object.
##' 
##' @param param A \code{PeakDensityParam} object containing all settings for
##' the peak grouping algorithm.
##'
##' @return For \code{groupChromPeaks}: a \code{\link{XCMSnExp}} object with the
##' results of the correspondence analysis. The definition of the resulting mz-rt
##' features can be accessed with the \code{\link{featureDefinitions}} method.
##' 
##' @seealso \code{\link{XCMSnExp}} for the object containing the results of
##' the correspondence.
##' 
##' @rdname groupChromPeaks-density
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
              df$peakidx <- res$peakIndex
              featureDefinitions(object) <- df
              if (validObject(object))
                  return(object)
          })


##' @title Single-spectrum non-chromatography MS data peak grouping
##'
##' @description \code{groupChromPeaks,XCMSnExp,MzClustParam}:
##' performs high resolution peak grouping for single spectrum
##' metabolomics data.
##'
##' @note Calling \code{groupChromPeaks} on an \code{XCMSnExp} object will cause
##' all eventually present previous correspondence results to be dropped.
##'
##' @param object For \code{groupChromPeaks}: an \code{\link{XCMSnExp}} object
##' containing the results from a previous chromatographic peak detection
##' analysis (see \code{\link{findChromPeaks}}).
##'
##' For all other methods: a \code{MzClustParam} object.
##' 
##' @param param A \code{MzClustParam} object containing all settings for
##' the peak grouping algorithm.
##'
##' @return For \code{groupChromPeaks}: a \code{\link{XCMSnExp}} object with the
##' results of the peak grouping step (i.e. the features). These can be accessed
##' with the \code{\link{featureDefinitions}} method.
##' 
##' @seealso \code{\link{XCMSnExp}} for the object containing the results of
##' the peak grouping.
##' 
##' @rdname groupChromPeaks-mzClust
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
              featureDefinitions(object) <- df
              if (validObject(object))
                  return(object)
          })


##' @title Peak grouping/correspondence based on proximity in the mz-rt space
##'
##' @description \code{groupChromPeaks,XCMSnExp,NearestPeaksParam}:
##' performs peak grouping based on the proximity between chromatographic
##' peaks from different samples in the mz-rt range.
##'
##' @note Calling \code{groupChromPeaks} on an \code{XCMSnExp} object will cause
##' all eventually present previous alignment results to be dropped.
##'
##' @param object For \code{groupChromPeaks}: an \code{\link{XCMSnExp}} object
##' containing the results from a previous chromatographic peak detection
##' analysis (see \code{\link{findChromPeaks}}).
##'
##' For all other methods: a \code{NearestPeaksParam} object.
##' 
##' @param param A \code{NearestPeaksParam} object containing all settings for
##' the peak grouping algorithm.
##'
##' @return For \code{groupChromPeaks}: a \code{\link{XCMSnExp}} object with the
##' results of the peak grouping/correspondence step (i.e. the mz-rt features).
##' These can be accessed with the \code{\link{featureDefinitions}} method.
##' 
##' @seealso \code{\link{XCMSnExp}} for the object containing the results of
##' the peak grouping.
##' 
##' @rdname groupChromPeaks-nearest
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
              featureDefinitions(object) <- df
              if (validObject(object))
                  return(object)
          })

##' @title Retention time correction based on alignment of house keeping peak
##' groups
##'
##' @description \code{adjustRtime,XCMSnExp,PeakGroupsParam}:
##' performs retention time correction based on the alignment of peak groups
##' (features) found in all/most samples.
##'
##' @note This method requires that a correspondence has been performed on the
##' data (see \code{\link{groupChromPeaks}}). Calling \code{adjustRtime} on an
##' \code{XCMSnExp} object will cause all peak grouping (correspondence) results
##' to be dropped after rt correction.
##'
##' @param object For \code{adjustRtime}: an \code{\link{XCMSnExp}} object
##' containing the results from a previous chromatographic peak detection (see
##' \code{\link{findChromPeaks}}) and alignment analysis (see
##' \code{\link{groupChromPeaks}}).
##'
##' For all other methods: a \code{PeakGroupsParam} object.
##' 
##' @param param A \code{PeakGroupsParam} object containing all settings for
##' the retention time correction method..
##'
##' @return For \code{adjustRtime}: a \code{\link{XCMSnExp}} object with the
##' results of the retention time adjustment step. These can be accessed with the
##' \code{\link{adjustedRtime}} method. Retention time correction does also adjust
##' the retention time of the identified chromatographic peaks (accessed
##' \emph{via} \code{\link{chromPeaks}}. Note that retention time correction
##' drops all previous alignment results from the result object.
##' 
##' @seealso \code{\link{XCMSnExp}} for the object containing the results of
##' the alignment.
##' 
##' @rdname adjustRtime-peakGroups
setMethod("adjustRtime",
          signature(object = "XCMSnExp", param = "PeakGroupsParam"),
          function(object, param) {
              if (!hasChromPeaks(object))
                  stop("No chromatographic peak detection results in 'object'! ",
                       "Please perform first a peak detection using the ",
                       "'findChromPeaks' method.")
              if (!hasFeatures(object))
                  stop("No feature definitions found in 'object'! Please ",
                       "perform first a peak grouping using the ",
                       "'groupChromPeak' method.")
              startDate <- date()
              res <- do_adjustRtime_peakGroups(
                  chromPeaks(object),
                  peakIndex = featureDefinitions(object)$peakidx,
                  rtime = rtime(object, bySample = TRUE),
                  minFraction = minFraction(param),
                  extraPeaks = extraPeaks(param),
                  smooth = smooth(param),
                  span = span(param),
                  family = family(param)
              )
              ## Dropping the peak groups but don't remove its process history
              ## step.
              ph <- processHistory(object, type = .PROCSTEP.PEAK.GROUPING)
              object <- dropFeatureDefinitions(object)
              ## Add the results. adjustedRtime<- should also fix the retention
              ## times for the peaks! Want to keep also the lates alignment
              ## information
              adjustedRtime(object) <- res
              if (length(ph)) {
                  object <- addProcessHistory(object, ph[[length(ph)]])
              }
              ## Add the process history step.
              xph <- XProcessHistory(param = param, date. = startDate,
                                     type. = .PROCSTEP.RTIME.CORRECTION,
                                     fileIndex = 1:length(fileNames(object)))
              object <- addProcessHistory(object, xph)
              if (validObject(object))
                  return(object)
          })


##' @title Align retention times across samples using Obiwarp
##'
##' @description \code{adjustRtime,XCMSnExp,ObiwarpParam}:
##' performs retention time correction/alignment based on the total mz-rt data
##' using the \emph{obiwarp} method.
##'
##' @note Calling \code{adjustRtime} on an \code{XCMSnExp} object will cause
##' all peak grouping (correspondence) results to be dropped.
##'
##' @param object For \code{adjustRtime}: an \code{\link{XCMSnExp}} object.
##'
##' For all other methods: a \code{ObiwarpParam} object.
##' 
##' @param param A \code{ObiwarpParam} object containing all settings for
##' the alignment method.
##'
##' @return For \code{adjustRtime,XCMSnExp,ObiwarpParam}: a
##' \code{\link{XCMSnExp}} object with the results of the retention time
##' adjustment step. These can be accessed with the \code{\link{adjustedRtime}}
##' method. Retention time correction does also adjust the retention time of the
##' identified chromatographic peaks (accessed \emph{via}
##' \code{\link{chromPeaks}}. Note that retention time correction drops all
##' previous peak grouping results from the result object.
##'
##' For \code{adjustRtime,OnDiskMSnExp,ObiwarpParam}: a \code{numeric} with the
##' adjusted retention times per spectra (in the same order than \code{rtime}).
##' 
##' @seealso \code{\link{XCMSnExp}} for the object containing the results of
##' the alignment.
##'
##' @references
##' John T. Prince and Edward M. Marcotte. "Chromatographic Alignment of
##' ESI-LC-MS Proteomic Data Sets by Ordered Bijective Interpolated Warping"
##' \emph{Anal. Chem.} 2006, 78 (17), 6140-6152.
##' 
##' @rdname adjustRtime-obiwarp
setMethod("adjustRtime",
          signature(object = "XCMSnExp", param = "ObiwarpParam"),
          function(object, param) {
              ## We don't require any detected or aligned peaks.
              startDate <- date()
              res <- .obiwarp(as(object, "OnDiskMSnExp"), param = param)
              ## Dropping the feature groups.
              object <- dropFeatureDefinitions(object)
              ## Add the results. adjustedRtime<- should also fix the retention
              ## times for the peaks! Want to keep also the lates alignment
              ## information
              adjustedRtime(object) <- res
              ## Add the process history step.
              xph <- XProcessHistory(param = param, date. = startDate,
                                     type. = .PROCSTEP.RTIME.CORRECTION,
                                     fileIndex = 1:length(fileNames(object)))
              object <- addProcessHistory(object, xph)
              if (validObject(object))
                  return(object)
          })

## profMat for XCMSnExp
##' @rdname XCMSnExp-class
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


##' @title Accessing mz-rt feature data values
##' 
##' @description \code{groupval,XCMSnExp}: extract a \code{matrix} for feature
##' values with rows representing features and columns samples. Parameter
##' \code{value} allows to define which column from the \code{\link{chromPeaks}}
##' matrix should be returned. Multiple chromatographic peaks from the same
##' sample can be assigned to a feature. Parameter \code{method} allows to
##' specify the method to be used in such cases to chose from which of the peaks
##' the value should be returned.
##'
##' @param object A \code{\link{XCMSnExp}} object providing the feature
##' definitions.
##' 
##' @param method \code{character} specifying the method to resolve
##' multi-peak mappings within the same sample, i.e. to define the
##' \emph{representative} peak for a feature in samples where more than
##' one peak was assigned to the feature. If \code{"medret"}: select the
##' peak closest to the median retention time of the feature.
##' If \code{"maxint"}: select the peak yielding the largest signal.
##'
##' @param value \code{character} specifying the name of the column in
##' \code{chromPeaks(object)} that should be returned or \code{"index"} (the
##' default) to return the index of the peak in the \code{chromPeaks(object)}
##' matrix corresponding to the \emph{representative} peak for the feature
##' in the respective sample.
##'
##' @param intensity \code{character} specifying the name of the column in the
##' \code{chromPeaks(objects)} matrix containing the intensity value of the
##' peak that should be used for the conflict resolution if
##' \code{method = "maxint"}.
##'
##' @return For \code{groupval}: a \code{matrix} with feature values, columns
##' representing samples, rows features. The order of the features
##' matches the order found in the \code{featureDefinitions(object)}
##' \code{DataFrame}. An \code{NA} is reported for features without corresponding
##' chromatographic peak in the respective sample(s).
##' 
##' @author Johannes Rainer
##' 
##' @seealso
##' \code{\link{XCMSnExp}} for information on the data object.
##' \code{\link{featureDefinitions}} to extract the \code{DataFrame} with the
##' feature definitions.
##' \code{\link{hasFeatures}} to evaluate whether the
##' \code{\link{XCMSnExp}} provides feature definitions.
##' 
##' @rdname XCMSnExp-peak-grouping-results
setMethod("groupval",
          signature(object = "XCMSnExp"),
          function(object, method = c("medret", "maxint"), value = "index",
                   intensity = "into") {
              ## Input argument checkings
              if (!hasFeatures(object))
                  stop("No peak groups present! Use 'groupChromPeaks' first.")
              if (!hasChromPeaks(object))
                  stop("No detected chromatographic peaks present! Use ",
                       "'findChromPeaks' first.")
              method <- match.arg(method)
              fNames <- basename(fileNames(object))
              nSamples <- seq_along(fNames)
              ## Copy all of the objects to avoid costly S4 method calls -
              ## improves speed at the cost of higher memory demand.
              fts <- chromPeaks(object)
              grps <- featureDefinitions(object)
              ftIdx <- grps$peakidx
              ## Match columns
              idx_rt <- match("rt", colnames(fts))
              idx_int <- match(intensity, colnames(fts))
              idx_samp <- match("sample", colnames(fts))
              
              vals <- matrix(nrow = length(ftIdx), ncol = length(nSamples))
              
              ## Get the indices for the elements.
              if (method == "medret") {
                  medret <- grps$rtmed
                  for (i in seq_along(ftIdx)) {
                      gidx <- ftIdx[[i]][base::order(base::abs(fts[ftIdx[[i]],
                                                                   idx_rt] -
                                                               medret[i]))]
                      vals[i, ] <- gidx[base::match(nSamples, fts[gidx,
                                                                  idx_samp])]
                  }
              } else {
                  for (i in seq_along(ftIdx)) {
                      gidx <- ftIdx[[i]][base::order(fts[ftIdx[[i]], idx_int],
                                                     decreasing = TRUE)]
                      vals[i, ] <- gidx[base::match(nSamples, fts[gidx, idx_samp])]
                  }
              }
              
              if (value != "index") {
                  if (!any(colnames(fts) == value))
                      stop("Column '", value,
                           "' not present in the chromatographic peaks matrix!")
                  vals <- fts[vals, value]
                  dim(vals) <- c(length(ftIdx), length(nSamples))
              }
              colnames(vals) <- fNames
              ## Let's skip row names for now.
              ## rownames(vals) <- paste(base::round(grps$mzmed, 3),
              ##                         base::round(grps$rtmed), sep = "/")
              return(vals)
})

#' @title Extracting chromatograms
#'
#' @description \code{extractChromatograms}: the method allows to extract
#' chromatograms from \code{\link[MSnbase]{OnDiskMSnExp}} and
#' \code{\link{XCMSnExp}} objects.
#'
#' @details Arguments \code{rt} and \code{mz} allow to specify the MS
#' data slice from which the chromatogram should be extracted. The parameter
#' \code{aggregationSum} allows to specify the function to be used to aggregate
#' the intensities across the mz range for the same retention time. Setting
#' \code{aggregationFun = "sum"} would e.g. allow to calculate the \emph{total
#' ion chromatogram} (TIC), \code{aggregationFun = "max"} the \emph{base peak
#' chromatogram} (BPC).
#'
#' @note
#' \code{Chromatogram} objects extracted with \code{extractChromatogram} contain
#' \code{NA_real_} values if, for a given retention time, no valid measurement
#' was available for the provided mz range.
#'
#' For \code{\link{XCMSnExp}} objects, if adjusted retention times are
#' available, the \code{extractChromatograms} method will by default report and
#' use these (for the subsetting based on the provided parameter \code{rt}). This
#' can be overwritten with the parameter \code{adjustedRtime}.
#' 
#' @param object Either a \code{\link[MSnbase]{OnDiskMSnExp}} or
#' \code{\link{XCMSnExp}} object from which the chromatograms should be extracted.
#'
#' @param rt \code{numeric(2)} defining the lower and upper boundary for the
#' retention time range. If not specified, the full retention time range of the
#' original data will be used. It is also possible to submit a \code{numeric(1)}
#' in which case \code{range} is called on it to transform it to a
#' \code{numeric(2)}.
#'
#' @param mz \code{numeric(2)} defining the lower and upper mz value for the
#' MS data slice. If not specified, the chromatograms will be calculated on the
#' full mz range. It is also possible to submit a \code{numeric(1)} in which case
#' \code{range} is called on it to transform it to a \code{numeric(2)}.
#'
#' @param adjustedRtime For \code{extractChromatograms,XCMSnExp}: whether the
#' adjusted (\code{adjustedRtime = TRUE}) or raw retention times
#' (\code{adjustedRtime = FALSE}) should be used for filtering and returned in
#' the resulting \code{\link{Chromatogram}} object. Adjusted retention times are
#' used by default if available.
#'
#' @param aggregationFun \code{character} specifying the function to be used to
#' aggregate intensity values across the mz value range for the same retention
#' time. Allowed values are \code{"sum"}, \code{"max"}, \code{"mean"} and
#' \code{"min"}.
#' 
#' @author Johannes Rainer
#'
#' @seealso \code{\link{XCMSnExp}} for the data object.
#' \code{\link{Chromatogram}} for the object representing chromatographic data.
#'
#' @export
#' @rdname extractChromatograms-method
#'
#' @examples
#' ## Read some files from the faahKO package.
#' library(xcms)
#' library(faahKO)
#' faahko_3_files <- c(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
#'                     system.file('cdf/KO/ko16.CDF', package = "faahKO"),
#'                     system.file('cdf/KO/ko18.CDF', package = "faahKO"))
#'
#' od <- readMSData2(faahko_3_files)
#'
#' ## Extract the ion chromatogram for one chromatographic peak in the data.
#' chrs <- extractChromatograms(od, rt = c(2700, 2900), mz = 335)
#'
#' ## plot the data
#' plot(rtime(chrs[[2]]), intensity(chrs[[2]]), type = "l", xlab = "rtime",
#'      ylab = "intensity", col = "000080")
#' for(i in c(1, 3)) {
#'   points(rtime(chrs[[i]]), intensity(chrs[[i]]), type = "l", col = "00000080")
#' }
setMethod("extractChromatograms",
          signature(object = "XCMSnExp"),
          function(object, rt, mz, adjustedRtime = hasAdjustedRtime(object),
                   aggregationFun = "sum") {
              return(.extractChromatogram(x = object, rt = rt, mz = mz,
                                          aggregationFun = aggregationFun,
                                          adjusted = adjustedRtime))
          })


