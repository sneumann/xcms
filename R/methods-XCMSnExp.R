## Methods for the XCMSnExp object representing untargeted metabolomics
## results
#' @include functions-XCMSnExp.R

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
##' @description The \code{adjustedRtime},\code{adjustedRtime<-} method
##' extract/set adjusted retention times. Retention times are adjusted by
##' retention time correction/adjustment methods. The \code{bySample} parameter
##' allows to specify whether the adjusted retention time should be grouped by
##' sample (file).
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
    newFd <- new("MsFeatureData")
    newFd@.xData <- .copy_env(object@msFeatureData)
    adjustedRtime(newFd) <- value
    lockEnvironment(newFd, bindings = TRUE)
    object@msFeatureData <- newFd
    if (validObject(object)) {
        ## Lock the environment so that only accessor methods can change values.
        ## lockEnvironment(newFd, bindings = TRUE)
        ## object@msFeatureData <- newFd
        return(object)
    }
})

##' @aliases featureGroups
##'
##' @description The \code{featureGroups}, \code{featureGroups<-} methods extract
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
##' @description The \code{features}, \code{features<-} methods extract or set
##' the matrix containing the information on identified features. See description
##' on the return value for details on the matrix columns. Users usually don't
##' have to use the \code{features<-} method directly as detected features are
##' added to the object by the \code{\link{detectFeatures}} method.
##'
##' @return For \code{features}: a \code{matrix} with at least the following
##' columns: \code{"mz"} (mz value for the largest intensity), \code{"mzmin"} (
##' minimal mz value), \code{"mzmax"} (maximal mz value), \code{"rt"} (retention
##' time for the peak apex), \code{"rtmin"} (minimal retention time), \code{"rtmax"}
##' (maximal retention time), \code{"into"} (integrated, original, intensity of
##' the feature) and \code{"sample"} (sample index in which the feature was
##' identified). Depending on the employed feature detection algorithm and the
##' \code{verboseColumns} parameter of it additional columns might be returned.
##'
##' @rdname XCMSnExp-class
setMethod("features", "XCMSnExp", function(object) {
    return(features(object@msFeatureData))
})
##' @aliases features<-
##'
##' @rdname XCMSnExp-class
setReplaceMethod("features", "XCMSnExp", function(object, value) {
    newFd <- new("MsFeatureData")
    newFd@.xData <- .copy_env(object@msFeatureData)
    features(newFd) <- value
    lockEnvironment(newFd, bindings = TRUE)
    object@msFeatureData <- newFd
    if (validObject(object)) {
        ## ## Lock the environment so that only accessor methods can change values.
        ## lockEnvironment(newFd, bindings = TRUE)
        ## object@msFeatureData <- newFd
        return(object)
    }
})

##' @description The \code{rtime} method extracts the retention time for each
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
    if (bySample)
        res <- split(res, fromFile(object))
    return(res)
})

##' @description The \code{mz} method extracts the mz values from each scan of
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
    if (bySample)
        res <- lapply(split(res, fromFile(object)), unlist, use.names = FALSE)
    return(res)
})

##' @description The \code{intensity} method extracts the intensity values from
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
    if (bySample)
        res <- lapply(split(res, fromFile(object)), unlist, use.names = FALSE)
    return(res)
})

## processHistory
##' @aliases processHistory
##' @description The \code{processHistory} method returns a \code{list} with
##' \code{\link{ProcessHistory}} objects (or objects inheriting from this base
##' class) representing the individual processing steps that have been performed,
##' eventually along with their settings (\code{Param} parameter class).
##'
##' @param fileIndex For \code{processHistory}: optional \code{numeric}
##' specifying the index of the files/samples for which the
##' \code{\link{ProcessHistory}} objects should be retrieved.
##'
##' @return For \code{processHistory}: a \code{list} of
##' \code{\link{ProcessHistory}} objects providing the details of the individual
##' data processing steps that have been performed.
##'
##' @rdname XCMSnExp-class
setMethod("processHistory", "XCMSnExp", function(object, fileIndex) {
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
        return(ph)
    } else {
        return(list())
    }
})

##' @description The \code{addProcessHistory} method adds (appends) a single
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
##' @description The \code{dropFeatures} method drops any identified features
##' and returns the object without that information. Note that for
##' \code{XCMSnExp} objects the method drops all results from a feature alignment
##' or retention time adjustment too. For \code{XCMSnExp} objects the method
##' drops also any related process history steps.
##'
##' @rdname XCMSnExp-class
setMethod("dropFeatures", "XCMSnExp", function(object) {
    if (hasDetectedFeatures(object)) {
        object <- dropFeatureGroups(object)
        object <- dropAdjustedRtime(object)
        object <- dropProcessHistories(object, type = .PROCSTEP.FEATURE.DETECTION)
        ## idx_fd <- which(unlist(lapply(processHistory(object), processType)) ==
        ##                 .PROCSTEP.FEATURE.DETECTION)
        ## if (length(idx_fd) > 0)
        ##     object@.processHistory <- object@.processHistory[-idx_fd]
        newFd <- new("MsFeatureData")
        newFd@.xData <- .copy_env(object@msFeatureData)
        newFd <- dropFeatures(newFd)
        lockEnvironment(newFd, bindings = TRUE)
        object@msFeatureData <- newFd
    }
    if (validObject(object))
        return(object)
})
##' @aliases dropFeatureGroups
##'
##' @description The \code{dropFeatureGroups} method drops aligned feature
##' information (i.e. feature groups) and returns the object
##' without that information. Note that for \code{XCMSnExp} objects the method
##' drops also retention time adjustments.
##' For \code{XCMSnExp} objects the method drops also any related process history
##' steps.
##'
##' @rdname XCMSnExp-class
setMethod("dropFeatureGroups", "XCMSnExp", function(object) {
    if (hasAlignedFeatures(object)) {
        ## phTypes <- unlist(lapply(processHistory(object), processType))
        ## idx_fal <- which(phTypes == .PROCSTEP.FEATURE.ALIGNMENT)
        ## idx_art <- which(phTypes == .PROCSTEP.RTIME.CORRECTION)
        ## if (length(idx_fal) > 0)
        ##     object@.processHistory <- object@.processHistory[-idx_fal]
        object <- dropProcessHistories(object, type = .PROCSTEP.FEATURE.ALIGNMENT)
        newFd <- new("MsFeatureData")
        newFd@.xData <- .copy_env(object@msFeatureData)
        newFd <- dropFeatureGroups(newFd)
        lockEnvironment(newFd, bindings = TRUE)
        object@msFeatureData <- newFd
        if (hasAdjustedRtime(object)) {
            ## ALWAYS drop retention time adjustments, since these are performed
            ## after alignment.
            object <- dropAdjustedRtime(object)
        }
    }
    if (validObject(object))
        return(object)
})
##' @aliases dropAdjustedRtime
##'
##' @description The \code{dropAdjustedRtime} method drops any retention time
##' adjustment information and returns the object without adjusted retention
##' time. Note that for \code{XCMSnExp} objects the method drops also all feature
##' alignment results if these were performed after the retention time adjustment.
##' For \code{XCMSnExp} objects the method drops also any related process history
##' steps.
##'
##' @rdname XCMSnExp-class
setMethod("dropAdjustedRtime", "XCMSnExp", function(object) {
    if (hasAdjustedRtime(object)) {
        phTypes <- unlist(lapply(processHistory(object), function(z)
            processType(z)))
        idx_art <- which(phTypes == .PROCSTEP.RTIME.CORRECTION)
        idx_fal <- which(phTypes == .PROCSTEP.FEATURE.ALIGNMENT)
        ## Drop retention time
        object@.processHistory <- object@.processHistory[-idx_art]
        newFd <- new("MsFeatureData")
        newFd@.xData <- .copy_env(object@msFeatureData)
        newFd <- dropAdjustedRtime(newFd)
        object@msFeatureData <- newFd
        lockEnvironment(newFd, bindings = TRUE)
        if (hasAlignedFeatures(object)) {
            if (max(idx_fal) > max(idx_art))
                object <- dropFeatureGroups(object)
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

##' @description The \code{bin} method allows to \emph{bin} spectra. See
##' \code{\link[MSnbase]{bin}} documentation for more details and examples.
##'
##' @param object An \code{\link{XCMSnExp}} object.
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

##' @description The \code{clean} method removes unused \code{0} intensity data
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

##' @description The \code{filterMsLevel} reduces the \code{\link{XCMSnExp}}
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

##' @description The \code{filterAcquisitionNum} method filters the
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
##' @description The \code{filterFile} method allows to reduce the
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

##' @description The \code{filterMz} method filters the data set based on the
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

##' @description The \code{filterRt} method filters the data set based on the
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
