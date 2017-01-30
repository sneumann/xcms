## Methods for the XCMSnExp object representing untargeted metabolomics
## results
#' @include functions-XCMSnExp.R do_groupFeatures-functions.R
#' do_adjustRtime-functions.R methods-xcmsRaw.R

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
    if (hasDetectedFeatures(object)) {
        cat("Feature detection:\n")
        cat(" ", nrow(features(object)), " features identified in ",
            length(fileNames(object)), " samples.\n", sep = "")
        cat(" On average ",
            format(mean(table(features(object)[, "sample"])), digits = 3),
            " features per sample.\n", sep = "")
    }
    if (hasAlignedFeatures(object)) {
        cat("Feature alignment:\n")
        cat(" ", nrow(featureGroups(object)), " featureGroups identified.\n",
            sep = "")
        cat(" Median mz range of feature groups: ",
            format(median(featureGroups(object)[, "mzmax"] -
                          featureGroups(object)[, "mzmin"]), digits = 5),
            "\n", sep = "")
        cat(" Median rt range of feature groups: ",
            format(median(featureGroups(object)[, "rtmax"] -
                          featureGroups(object)[, "rtmin"]), digits = 5),
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

##' @aliases hasAlignedFeatures
##'
##' @description \code{hasAlignedFeatures}: whether the object contains feature
##' alignment results.
##'
##' @rdname XCMSnExp-class
setMethod("hasAlignedFeatures", "XCMSnExp", function(object) {
    return(hasAlignedFeatures(object@msFeatureData))
})

##' @aliases hasDetectedFeatures
##'
##' @description \code{hasDetectedFeatures}: whether the object contains feature
##' detection results.
##'
##' @rdname XCMSnExp-class
setMethod("hasDetectedFeatures", "XCMSnExp", function(object) {
    return(hasDetectedFeatures(object@msFeatureData))
})

##' @aliases adjustedRtime
##'
##' @description \code{adjustedRtime},\code{adjustedRtime<-}:
##' extract/set adjusted retention times. \code{adjustedRtime<-} should not be
##' called manually, it is called internally by the \code{\link{adjustRtime}}
##' methods. For \code{XCMSnExp} objects, \code{adjustedRtime<-} does also apply
##' the retention time adjustment to the features in the object.
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
    if (hasDetectedFeatures(newFd)) {
        ## Change also the retention times reported in the features matrix.
        if (length(value) != length(rtime(object, bySample = TRUE)))
            stop("The length of 'value' has to match the number of samples!")
        message("Applying retention time adjustment to the identified",
                " features ... ", appendLF = FALSE)
        fts <- .applyRtAdjToFeatures(features(newFd),
                                     rtraw = rtime(object, bySample = TRUE),
                                     rtadj = value)
        ## Calling this on the MsFeatureData to avoid all results being removed
        ## again by the features<- method.
        features(newFd) <- fts
        message("OK")
    }
    lockEnvironment(newFd, bindings = TRUE)
    object@msFeatureData <- newFd
    if (validObject(object))
        return(object)
})

##' @aliases featureGroups
##'
##' @description \code{featureGroups}, \code{featureGroups<-}: extract
##' or set the feature alignment results.
##'
##' @return For \code{featureGroups}: a \code{DataFrame} with feature alignment
##' information, each row corresponding to one group of aligned features (across
##' samples) and columns \code{"mzmed"} (median mz value), \code{"mzmin"}
##' (minimal mz value), \code{"mzmax"} (maximum mz value), \code{"rtmed"} (median
##' retention time), \code{"rtmin"} (minimal retention time), \code{"rtmax"}
##' (maximal retention time) and \code{"featureidx"}. Column \code{"featureidx"}
##' contains a \code{list} with indices of features (rows) in the matrix returned
##' by the \code{features} method that belong to that feature group. The method
##' returns \code{NULL} if no aligned feature information is present.
##'
##' @rdname XCMSnExp-class
setMethod("featureGroups", "XCMSnExp", function(object) {
    return(featureGroups(object@msFeatureData))
})
##' @aliases featureGroups<-
##'
##' @rdname XCMSnExp-class
setReplaceMethod("featureGroups", "XCMSnExp", function(object, value) {
    if (hasAlignedFeatures(object))
        object <- dropFeatureGroups(object)
    newFd <- new("MsFeatureData")
    newFd@.xData <- .copy_env(object@msFeatureData)
    featureGroups(newFd) <- value
    lockEnvironment(newFd, bindings = TRUE)
    object@msFeatureData <- newFd
    if (validObject(object)) {
        ## Lock the environment so that only accessor methods can change values.
        ## lockEnvironment(newFd, bindings = TRUE)
        ## object@msFeatureData <- newFd
        return(object)
    }
})

##' @aliases features
##'
##' @description \code{features}, \code{features<-}: extract or set
##' the matrix containing the information on identified features. Parameter
##' \code{bySample} allows to specify whether features should be returned
##' ungrouped (default \code{bySample = FALSE}) or grouped by sample (
##' \code{bySample = TRUE}). The \code{features<-} method for \code{XCMSnExp}
##' objects removes also all feature alignment and retention time correction
##' results.
##' See description of the return value for details on the returned matrix. Users
##' usually don't have to use the \code{features<-} method directly as detected
##' features are added to the object by the \code{\link{detectFeatures}} method.
##'
##' @return For \code{features}: if \code{bySample = FALSE} a \code{matrix} with
##' at least the following columns: \code{"mz"} (mz value for the largest
##' intensity), \code{"mzmin"} (minimal mz value), \code{"mzmax"} (maximal mz
##' value), \code{"rt"} (retention time for the peak apex), \code{"rtmin"}
##' (minimal retention time), \code{"rtmax"} (maximal retention time),
##' \code{"into"} (integrated, original, intensity of the feature) and
##' \code{"sample"} (sample index in which the feature was identified).
##' Depending on the employed feature detection algorithm and the
##' \code{verboseColumns} parameter of it additional columns might be returned.
##' For \code{bySample = TRUE} the features are returned as a \code{list} of
##' matrices, each containing the features of a specific sample. For sample in
##' which no feastures were detected a matrix with 0 rows is returned.
##'
##' @rdname XCMSnExp-class
setMethod("features", "XCMSnExp", function(object, bySample = FALSE) {
    if (bySample) {
        tmp <- split(features(object), f = features(object)[, "sample"])
        ## Ensure we return something for each sample in case there is a sample
        ## without detected features.
        res <- vector("list", length(fileNames(object)))
        names(res) <- as.character(1:length(res))
        tmp <- split.data.frame(features(object), f = features(object)[, "sample"])
        res[as.numeric(names(tmp))] <- tmp
        if (any(lengths(res) == 0)) {
            emat <- matrix(nrow = 0, ncol = ncol(tmp[[1]]))
            colnamers(emat) <- colnames(tmp[[1]])
            res[lengths(res) == 0] <- emat
        }
        return(res)
    } else {
        return(features(object@msFeatureData))
    }
})
##' @aliases features<-
##'
##' @rdname XCMSnExp-class
setReplaceMethod("features", "XCMSnExp", function(object, value) {
    newFd <- new("MsFeatureData")
    ## Dropping all alignment results and all retention time corrections.
    suppressMessages(
        object <- dropFeatures(object)
    )
    ## Ensure that we remove ALL related process history steps
    newFd@.xData <- .copy_env(object@msFeatureData)
    features(newFd) <- value
    lockEnvironment(newFd, bindings = TRUE)
    object@msFeatureData <- newFd
    if (validObject(object)) {
        return(object)
    }
})

##' @description \code{rtime}: extracts the retention time for each
##' scan. The \code{bySample} parameter allows to return the values grouped
##' by sample/file.
##'
##' @param bySample logical(1) specifying whether results should be grouped by
##' sample.
##'
##' @return For \code{rtime}: if \code{bySample = FALSE} a numeric vector with the
##' retention times of each scan, if \code{bySample = TRUE} a \code{list} of
##' numeric vectors with the retention times per sample.
##'
##' @rdname XCMSnExp-class
setMethod("rtime", "XCMSnExp", function(object, bySample = FALSE) {
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
##' \code{object}. These values are extracted from the original data files and
##' eventual processing steps are applied \emph{on the fly}. Setting
##' \code{bySample = TRUE} the spectra are returned grouped by sample/file.
##'
##' @return For \code{spectra}: if \code{bySample = FALSE} a \code{list} with
##' \code{\link[MSnbase]{Spectrum}} objects. If \code{bySample = TRUE} the result
##' is grouped by sample, i.e. as a \code{list} of \code{lists}, each element in
##' the \emph{outer} \code{list} being the \code{list} of spectra of the specific
##' file.
##'
##' @rdname XCMSnExp-class
setMethod("spectra", "XCMSnExp", function(object, bySample = FALSE) {
    res <- callNextMethod(object = object)
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
##' Supported values are \code{"Unknown"}, \code{"Feature detection"},
##' \code{"Feature alignment"} and \code{"Retention time correction"}.
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

##' @aliases dropFeatures
##'
##' @description \code{dropFeatures}: drops any identified features
##' and returns the object without that information. Note that for
##' \code{XCMSnExp} objects the method drops all results from a feature alignment
##' or retention time adjustment too. For \code{XCMSnExp} objects the method
##' drops also any related process history steps.
##'
##' @rdname XCMSnExp-class
setMethod("dropFeatures", "XCMSnExp", function(object) {
    if (hasDetectedFeatures(object)) {
        object <- dropProcessHistories(object, type = .PROCSTEP.FEATURE.DETECTION)
        ## Make sure we delete all related process history steps
        object <- dropProcessHistories(object, type = .PROCSTEP.RTIME.CORRECTION)
        object <- dropProcessHistories(object, type = .PROCSTEP.FEATURE.ALIGNMENT)
        ## idx_fd <- which(unlist(lapply(processHistory(object), processType)) ==
        ##                 .PROCSTEP.FEATURE.DETECTION)
        ## if (length(idx_fd) > 0)
        ##     object@.processHistory <- object@.processHistory[-idx_fd]
        newFd <- new("MsFeatureData")
        newFd@.xData <- .copy_env(object@msFeatureData)
        newFd <- dropFeatures(newFd)
        ## Dropping other results from the environment (not the object).
        if (hasAdjustedRtime(newFd))
            newFd <- dropAdjustedRtime(newFd)
        if (hasAlignedFeatures(newFd))
            newFd <- dropFeatureGroups(newFd)
        lockEnvironment(newFd, bindings = TRUE)
        object@msFeatureData <- newFd
    }
    if (validObject(object))
        return(object)
})
##' @aliases dropFeatureGroups
##'
##' @description \code{dropFeatureGroups}: drops aligned feature
##' results (i.e. feature groups) and returns the object
##' without that information. Note that for \code{XCMSnExp} objects the method
##' will also drop retention time adjustment results, if these were performed
##' after the last feature alignment (i.e. which base on the results from the
##' feature alignment that are going to be removed). For \code{XCMSnExp} objects
##' also all related process history steps are removed.
##'
##' @param keepAdjRtime For \code{dropFeatureGroups,XCMSnExp}: logical(1)
##' defining whether eventually present retention time adjustment should not be
##' dropped. By default dropping feature groups drops retention time adjustment
##' results too.
##'
##' @param dropLastN For \code{dropFeatureGroups,XCMSnExp}: numeric(1) defining
##' the number of feature alignment related process history steps to remove. By
##' default \code{dropLastN = -1}, dropping the features removes all process
##' history steps related to feature alignment. Setting e.g. \code{dropLastN = 1}
##' will only remove the most recent feature alignment related process history
##' step.
##' 
##' @rdname XCMSnExp-class
setMethod("dropFeatureGroups", "XCMSnExp", function(object, keepAdjRtime = FALSE,
                                                    dropLastN = -1) {
    if (hasAlignedFeatures(object)) {
        phTypes <- unlist(lapply(processHistory(object), function(z)
            processType(z)))
        idx_art <- which(phTypes == .PROCSTEP.RTIME.CORRECTION)
        idx_fal <- which(phTypes == .PROCSTEP.FEATURE.ALIGNMENT)
        ## 1) drop last related process history step and results
        object <- dropProcessHistories(object,
                                       type = .PROCSTEP.FEATURE.ALIGNMENT,
                                       num = 1)
        newFd <- new("MsFeatureData")
        newFd@.xData <- .copy_env(object@msFeatureData)
        newFd <- dropFeatureGroups(newFd)
        lockEnvironment(newFd, bindings = TRUE)
        object@msFeatureData <- newFd
        ## 2) If retention time correction was performed after the latest feature
        ##    alignment, drop also the retention time correction and all related
        ##    process history steps.
        ##    Otherwise (grouping performed after retention time adjustment) do
        ##    nothing - this keeps eventual alignment related process history
        ##    steps performed before retention time correction.
        if (hasAdjustedRtime(object)) {
            if (max(idx_art) > max(idx_fal)) {
                object <- dropProcessHistories(object,
                                               type = .PROCSTEP.FEATURE.ALIGNMENT)
                ## This will ensure that the retention times of the features
                ## are restored.
                object <- dropAdjustedRtime(object)
                warning("Removed also feature alignment results as these based",
                        " on the retention time correction results that were",
                        " dropped.")
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
##' reported for the features in the feature matrix to the original, raw, ones
##' (after feature detection). Note that for \code{XCMSnExp} objects the method
##' drops also all feature alignment results if these were performed \emph{after}
##' the retention time adjustment. For \code{XCMSnExp} objects the method drops
##' also any related process history steps.
##'
##' @rdname XCMSnExp-class
setMethod("dropAdjustedRtime", "XCMSnExp", function(object) {
    if (hasAdjustedRtime(object)) {
        ## Get the process history types to determine the order of the analysis
        ## steps.
        phTypes <- unlist(lapply(processHistory(object), function(z)
            processType(z)))
        idx_art <- which(phTypes == .PROCSTEP.RTIME.CORRECTION)
        idx_fal <- which(phTypes == .PROCSTEP.FEATURE.ALIGNMENT)
        ## Copy the content of the object
        newFd <- new("MsFeatureData")
        newFd@.xData <- .copy_env(object@msFeatureData)        
        ## Revert applied adjustments in features:
        if (hasDetectedFeatures(newFd)) {
            message("Reverting retention times of identified features to ",
                    "original values ... ", appendLF = FALSE)
            fts <- .applyRtAdjToFeatures(features(newFd),
                                         rtraw = adjustedRtime(object,
                                                               bySample = TRUE),
                                         rtadj = rtime(object,
                                                       bySample = TRUE))
            ## Replacing features in MsFeatureData, not in XCMSnExp to avoid
            ## all results being removed.
            features(newFd) <- fts
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
        if (hasAlignedFeatures(object)) {
            if (max(idx_fal) > max(idx_art)) {
                object <- dropFeatureGroups(object)
                object <- dropProcessHistories(object,
                                               type = .PROCSTEP.FEATURE.ALIGNMENT,
                                               num = -1)
            }
        } else {
            ## If there is any feature alignment related process history, but no
            ## feature alignment results, drop them.
            object <- dropProcessHistories(object,
                                           type = .PROCSTEP.FEATURE.ALIGNMENT,
                                           num = -1)
        }
    }
    if (validObject(object))
        return(object)
})


############################################################
## Methods inherited from OnDiskMSnExp.
## For some of these methods (altering the raw data or subsetting) we have to
## remove the results.
## Methods to consider:
## o [ subset spectra, return an OnDiskMSnExp.
## o bin remove the results
## o clean remove the results.
## o featureNames returns the names of the spectra.
## o filterAcquisitionNum remove the results.
## o filterFile remove the results.
## o filterMsLevel remove the results.
## o filterMz
## o filterRt
## o fromFile<-
## o normalize remove the results.
## o pickPeaks remove the results.
## o removePeaks remove the results.
## o smooth remove the results.

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
              ## Want to support subsetting of the features!
              ## This means that we will also have to adjust the process
              ## history accordingly.
              if (hasAdjustedRtime(x) | hasAlignedFeatures(x) |
                  hasDetectedFeatures(x)) {
                  ## x@.processHistory <- list()
                  ## x@msFeatureData <- new("MsFeatureData")
                  x <- dropAdjustedRtime(x)
                  x <- dropFeatureGroups(x)
                  x <- dropFeatures(x)
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
    if (hasAdjustedRtime(object) | hasAlignedFeatures(object) |
        hasDetectedFeatures(object)) {
        ## object@.processHistory <- list()
        ## object@msFeatureData <- new("MsFeatureData")
        object <- dropAdjustedRtime(object)
        object <- dropFeatureGroups(object)
        object <- dropFeatures(object)
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
    if (hasAdjustedRtime(object) | hasAlignedFeatures(object) |
        hasDetectedFeatures(object)) {
        ## object@.processHistory <- list()
        ## object@msFeatureData <- new("MsFeatureData")
        object <- dropAdjustedRtime(object)
        object <- dropFeatureGroups(object)
        object <- dropFeatures(object)
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
    if (hasAdjustedRtime(object) | hasAlignedFeatures(object) |
        hasDetectedFeatures(object)) {
        ## object@.processHistory <- list()
        ## object@msFeatureData <- new("MsFeatureData")
        object <- dropAdjustedRtime(object)
        object <- dropFeatureGroups(object)
        object <- dropFeatures(object)
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
    if (hasAdjustedRtime(object) | hasAlignedFeatures(object) |
        hasDetectedFeatures(object)) {
        ## object@.processHistory <- list()
        ## object@msFeatureData <- new("MsFeatureData")
        object <- dropAdjustedRtime(object)
        object <- dropFeatureGroups(object)
        object <- dropFeatures(object)
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
##' \code{\link{XCMSnExp}} to data from only certain files. Identified features
##' for these files are retained while eventually all present feature
##' alignment/grouping information and adjusted retention times are dropped..
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
##' ## Perform feature detection on them using default matched filter settings.
##' mfp <- MatchedFilterParam()
##' xod <- detectFeatures(od, param = mfp)
##'
##' ## Subset the dataset to the first and third file.
##' xod_sub <- filterFile(xod, file = c(1, 3))
##'
##' ## The number of features per file for the full object
##' table(features(xod)[, "sample"])
##'
##' ## The number of features per file for the subset
##' table(features(xod_sub)[, "sample"])
##'
##' basename(fileNames(xod))
##' basename(fileNames(xod_sub))
setMethod("filterFile", "XCMSnExp", function(object, file) {
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
    ## Get the data we want to keep/subset
    fts <- features(object)
    if (hasAdjustedRtime(object)) {
        warning("Adjusted retention times removed.")
        object <- dropAdjustedRtime(object)
    }
    if (hasAlignedFeatures(object)) {
        warning("Feature alignment information removed.")
        object <- dropFeatureGroups(object)
    }
    ## Process the processing history.
    object <- dropProcessHistories(object, type = c(.PROCSTEP.FEATURE.ALIGNMENT,
                                                    .PROCSTEP.RTIME.CORRECTION))
    ph <- processHistory(object)
    ## The next method will actually clean everything, process history and
    ## msFeatureData
    suppressWarnings(
        object <- callNextMethod()
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
    ## Process features.
    fts <- fts[fts[, "sample"] %in% file, , drop = FALSE]
    fts[, "sample"] <- match(fts[, "sample"], file)
    features(object) <- fts
    object@.processHistory <- ph
    return(object)
})

##' @description \code{filterMz}: filters the data set based on the
##' provided mz value range. All features and feature groups (aligned features)
##' falling completely within the provided mz value range are retained (if their
##' minimal mz value is \code{>= mz[1]} and the maximal mz value \code{<= mz[2]}.
##' Adjusted retention times, if present, are not altered by the filtering.
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
    if (!is.numeric(mz) | length(mz) != 2)
        stop("'mz' has to be a numeric vector of length(2)!")
    mz <- range(mz)
    ## Subset features if present.
    object <- callNextMethod()  # just adds to processing queue.

    if (hasDetectedFeatures(object)) {
        fts <- features(object)
        keepIdx <- which(fts[, "mzmin"] >= mz[1] & fts[, "mzmax"] <= mz[2])
        newE <- .filterFeatures(object@msFeatureData, idx = keepIdx)
        lockEnvironment(newE, bindings = TRUE)
        object@msFeatureData <- newE
    }
    if (validObject(object))
        return(object)
})

##' @description \code{filterRt}: filters the data set based on the
##' provided retention time range. All features and feature groups within
##' the specified retention time window are retained. Filtering by retention time
##' does not drop any preprocessing results. The method returns an empty object
##' if no spectrum or feature is within the specified retention time range.
##'
##' @param rt For \code{filterRt}: \code{numeric(2)} defining the retention time
##' window (lower and upper bound) for the filtering.
##'
##' @param adjusted For \code{filterRt}: \code{logical} indicating whether the
##' object should be filtered by original (\code{adjusted = FALSE}) or adjusted
##' retention times (\code{adjusted = TRUE}).
##'
##' @rdname XCMSnExp-filter-methods
setMethod("filterRt", "XCMSnExp", function(object, rt, msLevel.,
                                           adjusted = FALSE) {
    if (missing(rt))
        return(object)
    if (!missing(msLevel.))
        warning("Parameter 'msLevel.' currently ignored.")
    rt <- range(rt)
    ## Get index of spectra within the rt window.
    ## Subset using [
    ## Subset features
    ## Subset feature groups
    ## Subset adjusted retention time
    if (!adjusted) {
        have_rt <- rtime(object)
    } else {
        have_rt <- adjustedRtime(object, bySample = FALSE)
        if (is.null(have_rt))
            stop("No adjusted retention time available!")
    }
    keep_logical <- have_rt >= rt[1] & have_rt <= rt[2]
    msg <- paste0("Filter: select retention time [",
                  paste0(rt, collapse = "-"),
                  "] and MS level(s), ",
                  paste(unique(msLevel(object)),
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
    ## 1) Subset features within the retention time range and feature groups.
    keep_fts <- numeric()
    if (hasDetectedFeatures(object)) {
        keep_fts <- which(features(object)[, "rtmin"] >= rt[1] &
                          features(object)[, "rtmax"] <= rt[2])
        if (length(keep_fts))
            newMfd <- .filterFeatures(object, idx = keep_fts)
            ## features(newMfd) <- features(object)[keep_fts, , drop = FALSE]
        else
            ph <- dropProcessHistoriesList(ph,
                                           type = c(.PROCSTEP.FEATURE.DETECTION,
                                                    .PROCSTEP.FEATURE.ALIGNMENT,
                                                    .PROCSTEP.RTIME.CORRECTION))
    }
    ## 2) Subset adjusted retention time
    if (hasAdjustedRtime(object) & length(keep_fts)) {
        ## Subset the adjusted retention times (which are stored as a list of
        ## rts by file):
        keep_by_file <- split(keep_logical, fromFile(object))
        adj_rt <- mapply(FUN = function(y, z) {
            return(y[z])
        }, y = adjustedRtime(object, bySample = TRUE), z = keep_by_file,
        SIMPLIFY = FALSE)
        adjustedRtime(newMfd) <- adj_rt
    }
    ## 3) Subset the OnDiskMSnExp part
    suppressWarnings(
        object <- object[which(keep_logical)]
    )
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
    if (hasAdjustedRtime(object) | hasAlignedFeatures(object) |
        hasDetectedFeatures(object)) {
        ## object@.processHistory <- list()
        ## object@msFeatureData <- new("MsFeatureData")
        object <- dropAdjustedRtime(object)
        object <- dropFeatureGroups(object)
        object <- dropFeatures(object)
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
    if (hasAdjustedRtime(object) | hasAlignedFeatures(object) |
        hasDetectedFeatures(object)) {
        ## object@.processHistory <- list()
        ## object@msFeatureData <- new("MsFeatureData")
        object <- dropAdjustedRtime(object)
        object <- dropFeatureGroups(object)
        object <- dropFeatures(object)
        warning("Removed preprocessing results")
    }
    callNextMethod()
})

##' The \code{removePeaks} method removes peaks (intensities) lower than a
##' threshold. Note that these peaks are not features! See \code{\link[MSnbase]{removePeaks}} documentation for details and examples.
##'
##' @param t For \code{removePeaks}: either a \code{numeric(1)} or \code{"min"}
##' defining the threshold (method) to be used. See
##' \code{\link[MSnbase]{removePeaks}} for details.
##'
##' @rdname XCMSnExp-inherited-methods
setMethod("removePeaks", "XCMSnExp", function(object, t = "min", verbose = FALSE,
                                              msLevel.) {
    if (hasAdjustedRtime(object) | hasAlignedFeatures(object) |
        hasDetectedFeatures(object)) {
        ## object@.processHistory <- list()
        ## object@msFeatureData <- new("MsFeatureData")
        object <- dropAdjustedRtime(object)
        object <- dropFeatureGroups(object)
        object <- dropFeatures(object)
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
    if (hasAdjustedRtime(x) | hasAlignedFeatures(x) |
        hasDetectedFeatures(x)) {
        ## x@.processHistory <- list()
        ## x@msFeatureData <- new("MsFeatureData")
        x <- dropAdjustedRtime(x)
        x <- dropFeatureGroups(x)
        x <- dropFeatures(x)
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


##' @title Feature alignment based on time dimension feature densities
##'
##' @description \code{groupFeatures,XCMSnExp,FeatureDensityParam}:
##' performs feature alignment (within and across samples) within in mz dimension
##' overlapping slices of MS data based on the density distribution of the
##' identified features in the slice along the time axis.
##'
##' @note Calling \code{groupFeatures} on an \code{XCMSnExp} object will cause
##' all eventually present previous alignment results to be dropped.
##'
##' @param object For \code{groupFeatures}: an \code{\link{XCMSnExp}} object
##' containing the results from a previous feature detection analysis (see
##' \code{\link{detectFeatures}}).
##'
##' For all other methods: a \code{FeatureDensityParam} object.
##' 
##' @param param A \code{FeatureDensityParam} object containing all settings for
##' the feature alignment algorithm.
##'
##' @return For \code{groupFeatures}: a \code{\link{XCMSnExp}} object with the
##' results of the feature alignment step. These can be accessed with the
##' \code{\link{featureGroups}} method.
##' 
##' @seealso \code{\link{XCMSnExp}} for the object containing the results of
##' the feature alignment.
##' 
##' @rdname groupFeatures-density
setMethod("groupFeatures",
          signature(object = "XCMSnExp", param = "FeatureDensityParam"),
          function(object, param) {
              if (!hasDetectedFeatures(object))
                  stop("No feature detection results in 'object'! Please ",
                       "perform first a feature detection using the ",
                       "'detectFeatures' method.")
              ## Get rid of any previous results.
              if (hasAlignedFeatures(object))
                  object <- dropFeatureGroups(object)
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
              res <- do_groupFeatures_density(features(object),
                                              sampleGroups = sampleGroups(param),
                                              bw = bw(param),
                                              minFraction = minFraction(param),
                                              minSamples = minSamples(param),
                                              binSize = binSize(param),
                                              maxFeatures = maxFeatures(param))
              xph <- XProcessHistory(param = param, date. = startDate,
                                     type. = .PROCSTEP.FEATURE.ALIGNMENT,
                                     fileIndex = 1:length(fileNames(object)))
              object <- addProcessHistory(object, xph)              
              ## Add the results.
              df <- DataFrame(res$featureGroups)
              df$featureidx <- res$featureIndex
              featureGroups(object) <- df
              if (validObject(object))
                  return(object)
          })


##' @title Single-spectrum non-chromatography MS data feature detection
##'
##' @description \code{groupFeatures,XCMSnExp,MzClustParam}:
##' performs high resolution feature alignment for single spectrum metabolomics
##' data.
##'
##' @note Calling \code{groupFeatures} on an \code{XCMSnExp} object will cause
##' all eventually present previous alignment results to be dropped.
##'
##' @param object For \code{groupFeatures}: an \code{\link{XCMSnExp}} object
##' containing the results from a previous feature detection analysis (see
##' \code{\link{detectFeatures}}).
##'
##' For all other methods: a \code{MzClustParam} object.
##' 
##' @param param A \code{MzClustParam} object containing all settings for
##' the feature alignment algorithm.
##'
##' @return For \code{groupFeatures}: a \code{\link{XCMSnExp}} object with the
##' results of the feature alignment step. These can be accessed with the
##' \code{\link{featureGroups}} method.
##' 
##' @seealso \code{\link{XCMSnExp}} for the object containing the results of
##' the feature alignment.
##' 
##' @rdname groupFeatures-mzClust
setMethod("groupFeatures",
          signature(object = "XCMSnExp", param = "MzClustParam"),
          function(object, param) {
              if (!hasDetectedFeatures(object))
                  stop("No feature detection results in 'object'! Please ",
                       "perform first a feature detection using the ",
                       "'detectFeatures' method.")
              ## I'm expecting a single spectrum per file!
              rtL <- split(rtime(object), f = fromFile(object))
              if (any(lengths(rtL) > 1))
                  stop("'object' contains multiple spectra per sample! This ",
                       "algorithm does only work for single spectra ",
                       "files/samples!")
              ## Get rid of any previous results.
              if (hasAlignedFeatures(object))
                  object <- dropFeatureGroups(object)
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
              res <- do_groupFeatures_mzClust(features(object),
                                              sampleGroups = sampleGroups(param),
                                              ppm = ppm(param),
                                              absMz = absMz(param),
                                              minFraction = minFraction(param),
                                              minSamples = minSamples(param))
              xph <- XProcessHistory(param = param, date. = startDate,
                                     type. = .PROCSTEP.FEATURE.ALIGNMENT,
                                     fileIndex = 1:length(fileNames(object)))
              object <- addProcessHistory(object, xph)              
              ## Add the results.
              df <- DataFrame(res$featureGroups)
              df$featureidx <- res$featureIndex
              featureGroups(object) <- df
              if (validObject(object))
                  return(object)
          })


##' @title Feature alignment based on proximity in the mz-rt space
##'
##' @description \code{groupFeatures,XCMSnExp,NearestFeaturesParam}:
##' performs feature alignment based on the proximity between features from
##' different samples in the mz-rt range.
##'
##' @note Calling \code{groupFeatures} on an \code{XCMSnExp} object will cause
##' all eventually present previous alignment results to be dropped.
##'
##' @param object For \code{groupFeatures}: an \code{\link{XCMSnExp}} object
##' containing the results from a previous feature detection analysis (see
##' \code{\link{detectFeatures}}).
##'
##' For all other methods: a \code{NearestFeaturesParam} object.
##' 
##' @param param A \code{NearestFeaturesParam} object containing all settings for
##' the feature alignment algorithm.
##'
##' @return For \code{groupFeatures}: a \code{\link{XCMSnExp}} object with the
##' results of the feature alignment step. These can be accessed with the
##' \code{\link{featureGroups}} method.
##' 
##' @seealso \code{\link{XCMSnExp}} for the object containing the results of
##' the feature alignment.
##' 
##' @rdname groupFeatures-nearest
setMethod("groupFeatures",
          signature(object = "XCMSnExp", param = "NearestFeaturesParam"),
          function(object, param) {
              if (!hasDetectedFeatures(object))
                  stop("No feature detection results in 'object'! Please ",
                       "perform first a feature detection using the ",
                       "'detectFeatures' method.")
              ## Get rid of any previous results.
              if (hasAlignedFeatures(object))
                  object <- dropFeatureGroups(object)
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
              res <- do_groupFeatures_nearest(features(object),
                                              sampleGroups = sampleGroups(param),
                                              mzVsRtBalance = mzVsRtBalance(param),
                                              absMz = absMz(param),
                                              absRt = absRt(param),
                                              kNN = kNN(param))
              xph <- XProcessHistory(param = param, date. = startDate,
                                     type. = .PROCSTEP.FEATURE.ALIGNMENT,
                                     fileIndex = 1:length(fileNames(object)))
              object <- addProcessHistory(object, xph)
              ## Add the results.
              df <- DataFrame(res$featureGroups)
              df$featureidx <- res$featureIndex
              featureGroups(object) <- df
              if (validObject(object))
                  return(object)
          })

##' @title Retention time correction based on alignment of house keeping feature
##' groups
##'
##' @description \code{adjustRtime,XCMSnExp,FeatureGroupsParam}:
##' performs retention time correction based on the alignment of feature groups
##' found in all/most samples.
##'
##' @note Calling \code{adjustRtime} on an \code{XCMSnExp} object will cause
##' all feature grouping (alignment) results to be dropped.
##'
##' @param object For \code{adjustRtime}: an \code{\link{XCMSnExp}} object
##' containing the results from a previous feature detection (see
##' \code{\link{detectFeatures}}) and alignment analysis (see
##' \code{\link{groupFeatures}}).
##'
##' For all other methods: a \code{FeatureGroupsParam} object.
##' 
##' @param param A \code{FeatureGroupsParam} object containing all settings for
##' the retention time correction method..
##'
##' @return For \code{adjustRtime}: a \code{\link{XCMSnExp}} object with the
##' results of the retention time adjustment step. These can be accessed with the
##' \code{\link{adjustedRtime}} method. Retention time correction does also adjust
##' the retention time of the identified features (accessed \emph{via}
##' \code{\link{features}}. Note that retention time correction drops
##' all previous alignment results from the result object.
##' 
##' @seealso \code{\link{XCMSnExp}} for the object containing the results of
##' the feature alignment.
##' 
##' @rdname adjustRtime-featureGroups
setMethod("adjustRtime",
          signature(object = "XCMSnExp", param = "FeatureGroupsParam"),
          function(object, param) {
              if (!hasDetectedFeatures(object))
                  stop("No feature detection results in 'object'! Please ",
                       "perform first a feature detection using the ",
                       "'detectFeatures' method.")
              if (!hasAlignedFeatures(object))
                  stop("No feature alignment results in 'object'! Please ",
                       "perform first a feature alignment using the ",
                       "'groupFeatures' method.")
              startDate <- date()
              res <- do_adjustRtime_featureGroups(features(object),
                                                  featureIndex = featureGroups(object)$featureidx,
                                                  rtime = rtime(object, bySample = TRUE),
                                                  minFraction = minFraction(param),
                                                  extraFeatures = extraFeatures(param),
                                                  smooth = smooth(param),
                                                  span = span(param),
                                                  family = family(param)
                                                  )
              ## Dropping the feature groups but don't remove its process history
              ## step.
              ph <- processHistory(object, type = .PROCSTEP.FEATURE.ALIGNMENT)
              object <- dropFeatureGroups(object)
              ## Add the results. adjustedRtime<- should also fix the retention
              ## times for the features! Want to keep also the lates alignment
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
        sps <- spectra(z)
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

## This is for test purposes until we really get to know how to use the function
.adjustRtime_obiwarp <- function(object, profStep = 1, centerSamples = NULL,
                                 response = 1, distFun = "cor_opt",
                                 gapInit = NULL, gapExtend = NULL,
                                 factorDiag = 2, factorGap = 1,
                                 localAlignment = FALSE, initPenalty = 0) {
    ## Load test files...
    faahko_3_files <- c(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
                        system.file('cdf/KO/ko16.CDF', package = "faahKO"),
                        system.file('cdf/KO/ko18.CDF', package = "faahKO"))
    faahko_od <- readMSData2(faahko_3_files)

    ## Feature alignment on those:
    object <- detectFeatures(faahko_od, param = CentWaveParam(noise = 10000,
                                                              snthresh = 40))
    xs <- xcmsSet(faahko_3_files, profparam = list(step = 0),
                  method = "centWave", noise = 10000, snthresh = 40)
    
    profStep = 1
    centerSample = NULL
    response = 1
    distFun = "cor_opt"
    gapInit = NULL
    gapExtend = NULL
    factorDiag = 2
    factorGap = 1
    localAlignment = FALSE
    initPenalty = 0
    ##

    ## Input checking:
    if (!hasDetectedFeatures(object))
        stop("No feature detection results in 'object'! Please ",
             "perform first a feature detection using the ",
             "'detectFeatures' method.")
    ## if (!hasAlignedFeatures(object))
    ##     stop("No feature alignment results in 'object'! Please ",
    ##          "perform first a feature alignment using the ",
    ##          "'groupFeatures' method.")
    startDate <- date()
    
    nSamples <- length(fileNames(object))

    ## centerSample
    if (!is.null(centerSample)) {
        if (length(centerSample) > 1 | !(centerSample %in% 1:nSamples))
            stop("'centerSample' has to be a single integer between 1 and ",
                 nSamples, "!")
    } else {
        tbl <- sort(table(features(object)[, "sample"]), decreasing = TRUE)
        centerSample <- as.numeric(names(tbl)[1])
    }
    ## Note: localAlignment has to be converted to a number.

    ## set default parameter. PUT THAT INTO THE PARAMETER CLASS!
    if (is.null(gapInit)) {
        if (distFunc == "cor") gapInit = 0.3
        if (distFunc == "cor_opt") gapInit = 0.3
        if (distFunc == "cov") gapInit = 0.0
        if (distFunc == "euc") gapInit = 0.9
        if (distFunc == "prd") gapInit = 0.0
    }
    if (is.null(gapExtend)) {
        if (distFunc == "cor") gapExtend = 2.4
        if (distFunc == "cor_opt") gapExtend = 2.4
        if (distFunc == "cov") gapExtend = 11.7
        if (distFunc == "euc") gapExtend = 1.8
        if (distFunc == "prd") gapExtend = 7.8
    }

    ## peakmat: features
    ## rtcor: rtime
    ## samples
    
    ## peakmat <- peaks(object)
    ## samples <- sampnames(object)
    ## classlabel <- as.vector(unclass(sampclass(object)))
    ## N <- length(samples)
    ## corpeaks <- peakmat
    ## plottype <- match.arg(plottype)

    ## Initialize result lists.
    rtimecor <- vector("list", nSamples)
    rtdevsmo <- vector("list", nSamples)

    ## Now for XCMSnExp (or OnDiskMSnExp) objects:
    ## o Have to ensure that the scan range of all objects fits with the center
    ##   sample - let's to that within the lapply though - just like in the
    ##   original function.
    ## o split it by file
    ## o Extract the center sample, remove it from the list
    ## o do an lapply over them:
    
    message("Sample number ", centerSample, " used as center sample.")

    ## Casting the XCMSnExp to an OnDiskMSnExp, so we don't run into any problems
    ## ... essentially, we could run the obiwarp ON A OnDiskMSnExp!!!
    
LLLLL    
###### BEWARE: OLD CODE BELOW.
    N <- nSamples
    center <- centerSample
    scanrange <- NULL
    idx <- which(seq(1,N) != center)
    obj1 <- xcmsRaw(xs@filepaths[center], profmethod="bin", profstep=0)

    ##     ## added t automatically find the correct scan range from the xcmsSet object
    ## if(length(obj1@scantime) != length(xs@rt$raw[[center]])){
    ##     ## This is in case the xcmsSet was read using a scanrange, i.e. if
    ##     ## the data was read in with defining a scan range, then we would have a
    ##     ## mismatch here. This code essentially ensures that the retention time
    ##     ## of the raw object would match the retention time present in the xcmsSet.
    ##     ## This was before the days in which @scanrange was added as a slot to
    ##     ## xcmsSet.
    ##     ##figure out the scan time range
    ##     scantime.start	<-object@rt$raw[[center]][1]
    ##     scantime.end	<-object@rt$raw[[center]][length(object@rt$raw[[center]])]
        
    ##     scanrange.start	<-which.min(abs(obj1@scantime - scantime.start))
    ##     scanrange.end	<-which.min(abs(obj1@scantime - scantime.end))
    ##     scanrange<-c(scanrange.start, scanrange.end)
    ##     obj1 <- xcmsRaw(object@filepaths[center], profmethod="bin", profstep=0, scanrange=scanrange)
    ## } else{
    ##     scanrange<-NULL
    ## }
    
    for (si in 1:length(idx)) {
        s <- idx[si]
        cat(samples[s], " ")

        ##
        ## Might be better to just get the profile matrix from the center object
        ## outside of the for loop and then modifying a internal variable within
        ## the loop - avoids creation of two profile matrices in each iteration.
        xcms:::profStepPad(obj1) <- profStep ## (re-)generate profile matrix, since it might have been modified during previous iteration
		if(is.null(scanrange)){
			obj2 <- xcmsRaw(xs@filepaths[s], profmethod="bin", profstep=0)
		} else{
			obj2 <- xcmsRaw(xs@filepaths[s], profmethod="bin", profstep=0, scanrange=scanrange)
		}
        xcms:::profStepPad(obj2) <- profStep ## generate profile matrix

        mzmin <-  min(obj1@mzrange[1], obj2@mzrange[1])
        mzmax <-  max(obj1@mzrange[2], obj2@mzrange[2])

        mz <- seq(mzmin,mzmax, by=profStep)
        mz <- as.double(mz)
        mzval <- length(mz)

        scantime1 <- obj1@scantime
        scantime2 <- obj2@scantime

        ## median difference between spectras' scan time.
        mstdiff <- median(c(diff(scantime1), diff(scantime2)))

        rtup1 <- c(1:length(scantime1))
        rtup2 <- c(1:length(scantime2))

        ## OK, will have to use that 1:1
        mst1 <- which(diff(scantime1) > 5 * mstdiff)[1]
        if(!is.na(mst1)) {
            rtup1 <- which(rtup1<=mst1)
            cat("Found gaps: cut scantime-vector at ", scantime1[mst1],"seconds", "\n")
        }

        mst2 <- which(diff(scantime2)>5*mstdiff)[1]
        if(!is.na(mst2)) {
            rtup2 <- which(rtup2<=mst2)
            cat("Found gaps: cut scantime-vector at ", scantime2[mst2],"seconds", "\n")
        }

        scantime1 <- scantime1[rtup1]
        scantime2 <- scantime2[rtup2]

        ## Drift of measured scan times - expected to be largest at the end.
        rtmaxdiff <- abs(diff(c(scantime1[length(scantime1)],
                                scantime2[length(scantime2)])))
        ## If the drift is larger than the threshold, cut the matrix up to the
        ## max allowed difference.
        if(rtmaxdiff > (5 * mstdiff)){
            rtmax <- min(scantime1[length(scantime1)],
                         scantime2[length(scantime2)])
            rtup1 <- which(scantime1<=rtmax)
            rtup2 <- which(scantime2<=rtmax)
        }

        scantime1 <- scantime1[rtup1]
        scantime2 <- scantime2[rtup2]
        valscantime1 <- length(scantime1)
        valscantime2 <- length(scantime2)

        ## Restrict the profile matrix to columns 1:valscantime
        if(length(obj1@scantime)>valscantime1) {
            obj1@env$profile <- obj1@env$profile[,-c((valscantime1+1):length(obj1@scantime))]
        }
        if(length(obj2@scantime)>valscantime2) {
            obj2@env$profile <- obj2@env$profile[,-c((valscantime2+1):length(obj2@scantime))]
        }

        ## Now ensure that the nrow of the profile matrix matches.
        ## Add empty rows at the beginning
        if(mzmin < obj1@mzrange[1]) {
            seqlen <- length(seq(mzmin, obj1@mzrange[1], profStep))-1
            x <- matrix(0, seqlen,dim(obj1@env$profile)[2])
            obj1@env$profile <- rbind(x, obj1@env$profile)
        }
        ## Add emtpy rows at the end.
        if(mzmax > obj1@mzrange[2]){
            seqlen <- length(seq(obj1@mzrange[2], mzmax, profStep))-1
            x <- matrix(0, seqlen, dim(obj1@env$profile)[2])
            obj1@env$profile <- rbind(obj1@env$profile, x)
        }
        if(mzmin < obj2@mzrange[1]){
            seqlen <- length(seq(mzmin, obj2@mzrange[1], profStep))-1
            x <- matrix(0, seqlen, dim(obj2@env$profile)[2])
            obj2@env$profile <- rbind(x, obj2@env$profile)
        }
        if(mzmax > obj2@mzrange[2]){
            seqlen <- length(seq(obj2@mzrange[2], mzmax, profStep))-1
            x <- matrix(0, seqlen, dim(obj2@env$profile)[2])
            obj2@env$profile <- rbind(obj2@env$profile, x)
        }

        ## OK, now that the matrices are "aligned" extract the intensities
        intensity1 <- obj1@env$profile
        intensity2 <- obj2@env$profile

        if ((mzval * valscantime1 != length(intensity1)) ||  (mzval * valscantime2 != length(intensity2)))
            stop("Dimensions of profile matrices do not match !\n")

        ## Would it be possible to supply non-binned data too???
        ## valsscantime1, valsscantime2: INTEGER
        ## scantime1, scantime2: REAL
        ## mzval: coerceVector INTSXP?
        ## mz1, mz2: REAL
        ## intensity1, intensity2: REAL
        ## response: INTEGER from 0 to 100
        ## localAlignment: INTEGER
        rtimecor[[s]] <-.Call("R_set_from_xcms",
                              valscantime1,scantime1,mzval,mz,intensity1,
                              valscantime2,scantime2,mzval,mz,intensity2,
                              response, distFunc,
                              gapInit, gapExtend,
                              factorDiag, factorGap,
                              localAlignment, initPenalty)

        ## Hm, silently add the raw retention times if we cut the retention time
        ## vector above - would merit at least a warning I believe.
        if(length(obj2@scantime) > valscantime2) {
            xs@rt$corrected[[s]] <- c(rtimecor[[s]],
                                          obj2@scantime[(max(rtup2)+1):length(obj2@scantime)])
        } else {
            xs@rt$corrected[[s]] <- rtimecor[[s]]
        }

        ## why are we rounding here???
        rtdevsmo[[s]] <- round(rtime[[s]]-object@rt$corrected[[s]],2)

        rm(obj2)
        gc()

        ## updateProgressInfo
        object@progressInfo$retcor.obiwarp <-  si / length(idx)
        xcms:::progressInfoUpdate(object)

    }

    cat("\n")
    rtdevsmo[[center]] <- round(rtime[[center]] - object@rt$corrected[[center]], 2)


    for (i in 1:N) {
        cfun <- stepfun(rtime[[i]][-1] - diff(rtime[[i]])/2, rtime[[i]] - rtdevsmo[[i]])
        rtime[[i]] <- rtime[[i]] - rtdevsmo[[i]]

        sidx <- which(corpeaks[,"sample"] == i)
        corpeaks[sidx, c("rt", "rtmin", "rtmax")] <- cfun(corpeaks[sidx, c("rt", "rtmin", "rtmax")])
    }

    peaks(object) <- corpeaks
    groups(object) <- matrix(nrow = 0, ncol = 0)
    groupidx(object) <- list()
    invisible(object)

}
