## Methods for the XCMSnExp object representing untargeted metabolomics
## results
#' @include functions-XCMSnExp.R do_groupChromPeaks-functions.R functions-utils.R
#' do_adjustRtime-functions.R methods-xcmsRaw.R functions-OnDiskMSnExp.R

setMethod("initialize", "XCMSnExp", function(.Object, ...) {
    .Object <- callNextMethod(.Object, ...)
    lockEnvironment(.Object@msFeatureData)
    return(.Object)
})

#' @aliases show,MsFeatureData-method
#'
#' @rdname XCMSnExp-class
setMethod("show", "XCMSnExp", function(object) {
    callNextMethod()
    cat("- - - xcms preprocessing - - -\n")
    if (hasChromPeaks(object)) {
        cat("Chromatographic peak detection:\n")
        ph <- processHistory(object, type = .PROCSTEP.PEAK.DETECTION)
        if (length(ph))
            cat(" method:", .param2string(ph[[1]]@param), "\n")
        else cat(" unknown method.\n")
        cat(" ", nrow(chromPeaks(object)), " peaks identified in ",
            length(fileNames(object)), " samples.\n", sep = "")
        cat(" On average ",
            format(mean(table(chromPeaks(object)[, "sample"])), digits = 3),
            " chromatographic peaks per sample.\n", sep = "")
    }
    if (hasAdjustedRtime(object)) {
        cat("Alignment/retention time adjustment:\n")
        ph <- processHistory(object, type = .PROCSTEP.RTIME.CORRECTION)
        if (length(ph))
            cat(" method:", .param2string(ph[[1]]@param), "\n")
        else cat(" unknown method.\n")
    }
    if (hasFeatures(object)) {
        cat("Correspondence:\n")
        ph <- processHistory(object, type = .PROCSTEP.PEAK.GROUPING)
        if (length(ph))
            cat(" method:", .param2string(ph[[1]]@param), "\n")
        else cat(" unknown method.\n")
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
            fp <- chromPeaks(object)[chromPeakData(object)$is_filled, ,
                                     drop = FALSE]
            cat("", sum(chromPeakData(object)$is_filled),
                "filled peaks (on average",
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
setMethod("hasFeatures", "XCMSnExp", function(object, msLevel = integer()) {
    hasFeatures(object@msFeatureData, msLevel = msLevel)
})

#' @aliases hasChromPeaks hasChromPeaks,MsFeatureData-method
#'
#' @description
#'
#' \code{hasChromPeaks}: whether the object contains peak
#' detection results.
#'
#' @rdname XCMSnExp-class
setMethod("hasChromPeaks", "XCMSnExp", function(object, msLevel = integer()) {
    hasChromPeaks(object@msFeatureData, msLevel = msLevel)
})

#' @aliases hasFilledChromPeaks
#'
#' @description
#'
#' \code{hasFilledChromPeaks}: whether the object contains any filled-in
#' chromatographic peaks.
#'
#' @rdname XCMSnExp-class
setMethod("hasFilledChromPeaks", "XCMSnExp", function(object) {
    .hasFilledPeaks(object)
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
        sNames <- unlist(split(featureNames(object),
                               as.factor(fromFile(object))),
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
    validObject(object)
    object
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
#' \code{featureDefinitions} supports also parameters \code{mz}, \code{rt},
#' \code{ppm} and \code{type} to return only features within certain ranges (see
#' description of \code{chromPeaks} for details).
#'
#' @rdname XCMSnExp-class
setMethod("featureDefinitions", "XCMSnExp",
          function(object, mz = numeric(), rt = numeric(), ppm = 0,
                   type = c("any", "within", "apex_within"),
                   msLevel = integer()) {
              feat_def <- featureDefinitions(object@msFeatureData,
                                             msLevel = msLevel)
              type <- match.arg(type)
              ## Select features within rt range.
              if (length(rt) && nrow(feat_def)) {
                  rt <- range(rt)
                  if (type == "any")
                      keep <- which(feat_def$rtmin <= rt[2] &
                                    feat_def$rtmax >= rt[1])
                  if (type == "within")
                      keep <- which(feat_def$rtmin >= rt[1] &
                                    feat_def$rtmax <= rt[2])
                  if (type == "apex_within")
                      keep <- which(feat_def$rtmed >= rt[1] &
                                    feat_def$rtmed <= rt[2])
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
                      keep <- which(feat_def$mzmin <= mz[2] &
                                    feat_def$mzmax >= mz[1])
                  if (type == "within")
                      keep <- which(feat_def$mzmin >= mz[1] &
                                    feat_def$mzmax <= mz[2])
                  if (type == "apex_within")
                      keep <- which(feat_def$mzmed >= mz[1] &
                                    feat_def$mzmed <= mz[2])
                  feat_def <- feat_def[keep, , drop = FALSE]
              }
              feat_def
          })
#' @aliases featureDefinitions<- featureDefinitions<-,MsFeatureData-method
#'
#' @rdname XCMSnExp-class
setReplaceMethod("featureDefinitions", "XCMSnExp", function(object, value) {
    ## if (hasFeatures(object))
    ##     object <- dropFeatureDefinitions(object)
    newFd <- new("MsFeatureData")
    newFd@.xData <- xcms:::.copy_env(object@msFeatureData)
    featureDefinitions(newFd) <- value
    lockEnvironment(newFd, bindings = TRUE)
    object@msFeatureData <- newFd
    validObject(object)
    object
})

#' @aliases chromPeaks chromPeaks,MsFeatureData-method chromPeakData,MsFeatureData-method chromPeakData
#'
#' @description
#'
#' \code{chromPeaks}, \code{chromPeaks<-}: extract or set
#' the matrix containing the information on identified chromatographic
#' peaks. Rownames of the matrix represent unique IDs of the respective peaks
#' within the experiment.
#' Parameter \code{bySample} allows to specify whether peaks should
#' be returned ungrouped (default \code{bySample = FALSE}) or grouped by
#' sample (\code{bySample = TRUE}). The \code{chromPeaks<-} method for
#' \code{XCMSnExp} objects removes also all correspondence (peak grouping)
#' and retention time correction (alignment) results. The optional
#' arguments \code{rt}, \code{mz}, \code{ppm} and \code{type} allow to extract
#' only chromatographic peaks overlapping the defined retention time and/or
#' m/z ranges. Argument \code{type} allows to define how \emph{overlapping} is
#' determined: for \code{type == "any"} (the default), all peaks that are even
#' partially overlapping the region are returned (i.e. for which either
#' \code{"mzmin"} or \code{"mzmax"} of the \code{chromPeaks} or
#' \code{featureDefinitions} matrix are within the provided m/z range), for
#' \code{type == "within"} the full peak has to be within the region (i.e.
#' both \code{"mzmin"} and \code{"mzmax"} have to be within the m/z range) and
#' for \code{type == "apex_within"} the peak's apex position (highest signal
#' of the peak) has to be within the region (i.e. the peak's or features m/z
#' has to be within the m/z range).
#' See description of the return value for details on the returned matrix.
#' Users usually don't have to use the \code{chromPeaks<-} method directly
#' as detected chromatographic peaks are added to the object by the
#' \code{\link{findChromPeaks}} method. Also, \code{chromPeaks<-} will replace
#' any existing \code{chromPeakData}.
#'
#' \code{chromPeakData} and \code{chromPeakData<-} allow to get or set arbitrary
#' chromatographic peak annotations. These are returned or ar returned as a
#' \code{DataFrame}. Note that the number of rows and the rownames of the
#' \code{DataFrame} have to match those of \code{chromPeaks}.
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
#' @param msLevel \code{integer} specifying the MS level(s) for which identified
#'     chromatographic peaks should be returned.
#'
#' @param isFilledColumn \code{logical(1)} whether a column \code{"is_filled"}
#'     is included in the returned \code{"matrix"} providing the information
#'     if a peak was filled in. Alternatively, this information would be
#'     provided by the \code{chromPeakData} data frame.
#'
#' @return
#'
#' For \code{chromPeaks}: if \code{bySample = FALSE} a \code{matrix} (each row
#' being a chromatographic peak, rownames representing unique IDs of the peaks)
#' with at least the following columns:
#' \code{"mz"} (intensity-weighted mean of mz values of the peak across
#' scans/retention times),
#' \code{"mzmin"} (minimal mz value),
#' \code{"mzmax"} (maximal mz value),
#' \code{"rt"} (retention time of the peak apex),
#' \code{"rtmin"} (minimal retention time),
#' \code{"rtmax"} (maximal retention time),
#' \code{"into"} (integrated, original, intensity of the peak),
#' \code{"maxo"} (maximum intentity of the peak),
#' \code{"sample"} (sample index in which the peak was identified) and
#' Depending on the employed peak detection algorithm and the
#' \code{verboseColumns} parameter of it, additional columns might be
#' returned. If parameter \code{isFilledColumn} was set to \code{TRUE} a column
#' named \code{"is_filled"} is also returned.
#' For \code{bySample = TRUE} the chromatographic peaks are
#' returned as a \code{list} of matrices, each containing the
#' chromatographic peaks of a specific sample. For samples in which no
#' peaks were detected a matrix with 0 rows is returned.
#'
#' @rdname XCMSnExp-class
setMethod("chromPeaks", "XCMSnExp", function(object, bySample = FALSE,
                                             rt = numeric(), mz = numeric(),
                                             ppm = 0, msLevel = integer(),
                                             type = c("any", "within",
                                                      "apex_within"),
                                             isFilledColumn = FALSE) {
    type <- match.arg(type)
    pks <- chromPeaks(object@msFeatureData)
    if (isFilledColumn)
        pks <- cbind(pks, is_filled = as.numeric(chromPeakData(object)$is_filled))
    if (length(msLevel))
        pks <- pks[which(chromPeakData(object)$ms_level %in% msLevel), ,
                   drop = FALSE]
    ## Select peaks within rt range.
    if (length(rt)) {
        rt <- range(as.numeric(rt))
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
        mz <- range(as.numeric(mz))
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
        if (nrow(pks)) {
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
#' @aliases chromPeaks<- chromPeaks<-,MsFeatureData-method chromPeakData<-,MsFeatureData-method chromPeakData<-
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
    ## Set rownames if not present
    if (is.null(rownames(value)) && nrow(value))
        rownames(value) <- .featureIDs(nrow(value), prefix = "CP")
    chromPeaks(newFd) <- value
    chromPeakData(newFd) <- DataFrame(ms_level = rep(1L, nrow(value)),
                                      is_filled = rep(FALSE, nrow(value)),
                                      row.names = rownames(value))
    lockEnvironment(newFd, bindings = TRUE)
    object@msFeatureData <- newFd
    validObject(object)
    object
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
        tmp <- split(res, as.factor(fromFile(object)))
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
        tmp <- lapply(split(res, as.factor(fromFile(object))),
                      unlist, use.names = FALSE)
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
        tmp <- lapply(split(res, as.factor(fromFile(object))),
                      unlist, use.names = FALSE)
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
        tmp <- split(res, as.factor(fromFile(object)))
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
#' \code{processHistory}: returns a \code{list} of
#' \code{\link{ProcessHistory}} objects (or objects inheriting from this
#' base class) representing the individual processing steps that have been
#' performed, eventually along with their settings (\code{Param} parameter
#' class). Optional arguments \code{fileIndex}, \code{type} and
#' \code{msLevel} allow to restrict to process steps of a certain type or
#' performed on a certain file or MS level.
#'
#' @param fileIndex For \code{processHistory}: optional \code{integer}
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
        object <- dropProcessHistories(
            object, type = c(.PROCSTEP.PEAK.DETECTION, .PROCSTEP.PEAK.GROUPING,
                             .PROCSTEP.PEAK.FILLING, .PROCSTEP.CALIBRATION,
                             .PROCSTEP.PEAK.REFINEMENT))
        newFd <- new("MsFeatureData")
        if (hasAdjustedRtime(object)) {
            idx_rt_adj <- which(phTypes == .PROCSTEP.RTIME.CORRECTION)
            idx_pk_det <- which(phTypes == .PROCSTEP.PEAK.DETECTION)
            if (length(idx_rt_adj) && length(idx_pk_det) &&
                max(idx_rt_adj) > max(idx_pk_det) && !keepAdjustedRtime) {
                object <- dropProcessHistories(
                    object, type = .PROCSTEP.RTIME.CORRECTION)
            } else
                keepAdjustedRtime <- TRUE
            if (keepAdjustedRtime)
                adjustedRtime(newFd) <- adjustedRtime(object@msFeatureData)
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
        object <- dropProcessHistories(
            object, type = .PROCSTEP.PEAK.GROUPING, num = 1)
        ## Drop also eventual filterFeatureDefinitions
        object <- dropGenericProcessHistory(
            object, fun = "filterFeatureDefinitions")
        newFd <- new("MsFeatureData")
        newFd@.xData <- .copy_env(object@msFeatureData)
        drop_proc_hist <- character()
        dropAdjustedRtime <- FALSE
        if (max(idx_art) > max(idx_fal) & !keepAdjustedRtime) {
            drop_proc_hist <- .PROCSTEP.PEAK.GROUPING
            dropAdjustedRtime <- TRUE
        }
        newFd <- dropFeatureDefinitions(newFd, dropAdjustedRtime)
        if (.hasFilledPeaks(object))
            drop_proc_hist <- c(drop_proc_hist, .PROCSTEP.PEAK.FILLING)
        lockEnvironment(newFd, bindings = TRUE)
        object@msFeatureData <- newFd
        if (!hasChromPeaks(object))
            drop_proc_hist <- c(drop_proc_hist, .PROCSTEP.PEAK.DETECTION,
                                .PROCSTEP.PEAK.REFINEMENT)
        if (length(drop_proc_hist))
            object <- dropProcessHistories(object, drop_proc_hist)
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
        object <- dropProcessHistories(
            object, type = .PROCSTEP.RTIME.CORRECTION, num = 1)
        newFd <- dropAdjustedRtime(
            newFd, rtime(object, bySample = TRUE, adjusted = FALSE))
        ## 2) If grouping has been performed AFTER retention time correction it
        ##    has to be dropped too, including ALL related process histories.
        if (hasFeatures(object)) {
            if (max(idx_fal) > max(idx_art)) {
                newFd <- dropFeatureDefinitions(newFd)
                object <- dropProcessHistories(
                    object, type = .PROCSTEP.PEAK.GROUPING, num = -1)
            }
        }  else {
            ## If there is any peak alignment related process history, but no
            ## peak alignment results, drop them.
            object <- dropProcessHistories(
                object, type = .PROCSTEP.PEAK.GROUPING, num = -1)
        }
        object@msFeatureData <- newFd
        lockEnvironment(newFd, bindings = TRUE)
    }
    if (validObject(object))
        return(object)
})


#' @title XCMSnExp filtering and subsetting
#'
#' @aliases XCMSnExp-filter
#'
#' @description
#'
#' The methods listed on this page allow to filter and subset [XCMSnExp]
#' objects. Most of them are inherited from the [OnDiskMSnExp] object defined
#' in the `MSnbase` package and have been adapted for `XCMSnExp` to enable
#' correct subsetting of preprocessing results.
#'
#' - `[`: subset a `XCMSnExp` object by spectra. Be aware that this removes
#'   **all** preprocessing results, except adjusted retention times if
#'   `keepAdjustedRtime = TRUE` is passed to the method.
#'
#' - `[[`: extracts a single `Spectrum` object (defined in `MSnbase`). The
#'   reported retention time is the adjusted retention time if alignment has
#'   been performed.
#'
#' - `filterChromPeaks`: subset the `chromPeaks` `matrix` in `object`. Parameter
#'   `method` allows to specify how the chromatographic peaks should be
#'   filtered. Currently, only `method = "keep"` is supported which allows to
#'   specify chromatographic peaks to keep with parameter `keep` (i.e. provide
#'   a `logical`, `integer` or `character` defining which chromatographic peaks
#'   to keep). Feature definitions (if present) are updated correspondingly.
#'
#' - `filterFeatureDefinitions`: allows to subset the feature definitions of
#'   an `XCMSnExp` object. Parameter `features` allow to define which features
#'   to keep. It can be a `logical`, `integer` (index of features to keep) or
#'   `character` (feature IDs) vector.
#'
#' - `filterFile`: allows to reduce the `XCMSnExp` to data from only selected
#'   files. Identified chromatographic peaks for these files are retained while
#'   correspondence results (feature definitions) are removed by default. To
#'   force keeping feature definitions use `keepFeatures = TRUE`. Adjusted
#'   retention times (if present) are retained by default if present. Use
#'   `keepAdjustedRtime = FALSE` to drop them.
#'
#' - `filterMsLevel`: reduces the `XCMSnExp` object to spectra of the
#'   specified MS level(s). Chromatographic peaks and identified features are
#'   also subsetted to the respective MS level. See also the `filterMsLevel`
#'   documentation in `MSnbase` for details and examples.
#'
#' - `filterMz`: filters the data set based on the provided m/z value range.
#'   All chromatographic peaks and features (grouped peaks) falling
#'   **completely** within the provided mz value range are retained
#'   (i.e. if their minimal m/z value is `>= mz[1]` and the maximal m/z value
#'   `<= mz[2]`. Adjusted retention times, if present, are kept.
#'
#' - `filterRt`: filters the data set based on the provided retention time
#'   range. All chromatographic peaks and features (grouped peaks)
#'   **completely** within the specified retention time window are retained
#'   (i.e. if the retention time corresponding to the peak's apex is within the
#'   specified rt range). If retention time correction has been performed,
#'   the method will by default filter the object by adjusted retention times.
#'   The argument `adjusted` allows to specify manually whether filtering
#'   should be performed on raw or adjusted retention times. Filtering by
#'   retention time does not drop any preprocessing results nor does it remove
#'   or change alignment results (i.e. adjusted retention times).
#'   The method returns an empty object if no spectrum or feature is within
#'   the specified retention time range.
#'
#' - `split`: splits an `XCMSnExp` object into a `list` of `XCMSnExp` objects
#'   based on the provided parameter `f`. Note that by default all
#'   pre-processing results are removed by the splitting, except adjusted
#'   retention times, if the optional argument `keepAdjustedRtime = TRUE` is
#'   provided.
#'
#' @details
#'
#' All subsetting methods try to ensure that the returned data is
#' consistent. Correspondence results for example are removed by default if the
#' data set is sub-setted by file, since the correspondence results are
#' dependent on the files on which correspondence was performed. This can be
#' changed by setting `keepFeatures = TRUE`.
#' For adjusted retention times, most subsetting methods
#' support the argument `keepAdjustedRtime` (even the `[` method)
#' that forces the adjusted retention times to be retained even if the
#' default would be to drop them.
#'
#' @note
#'
#' The `filterFile` method removes also process history steps not
#' related to the files to which the object should be sub-setted and updates
#' the `fileIndex` attribute accordingly. Also, the method does not
#' allow arbitrary ordering of the files or re-ordering of the files within
#' the object.
#'
#' Note also that most of the filtering methods, and also the subsetting
#' operations `[` drop all or selected preprocessing results. To
#' consolidate the alignment results, i.e. ensure that adjusted retention
#' times are always preserved, use the [applyAdjustedRtime()]
#' function on the object that contains the alignment results. This replaces
#' the raw retention times with the adjusted ones.
#'
#' @param adjusted For `filterRt`: `logical` indicating whether the
#'     object should be filtered by original (`adjusted = FALSE`) or
#'     adjusted retention times (`adjusted = TRUE`).
#'     For `spectra`: whether the retention times in the individual
#'     `Spectrum` objects should be the adjusted or raw retention times.
#'
#' @param drop For `[` and `[[`: not supported.
#'
#' @param f For `split` a vector of length equal to the length of x
#'     defining how `x` should be splitted. It is converted internally to
#'     a `factor`.
#'
#' @param features For `filterFeatureDefinitions`: either a `integer`
#'     specifying the indices of the features (rows) to keep, a `logical`
#'     with a length matching the number of rows of `featureDefinitions`
#'     or a `character` with the feature (row) names.
#'
#' @param file For `filterFile`: `integer` defining the file index
#'     within the object to subset the object by file or `character`
#'     specifying the file names to sub set. The indices are expected to be
#'     increasingly ordered, if not they are ordered internally.
#'
#' @param i For `[`: `numeric` or `logical` vector specifying to
#'     which spectra the data set should be reduced.
#'     For `[[`: a single integer or character.
#'
#' @param j For `[` and `[[`: not supported.
#'
#' @param keep For `filterChromPeaks`: `logical`, `integer` or `character`
#'     defining which chromatographic peaks should be retained.
#'
#' @param keepAdjustedRtime For `filterFile`, `filterMsLevel`,
#'     `[`, `split`: `logical(1)` defining whether the adjusted
#'     retention times should be kept, even if e.g. features are being removed
#'     (and the retention time correction was performed on these features).
#'
#' @param keepFeatures For `filterFile`: `logical(1)` whether
#'     correspondence results (feature definitions) should be kept or dropped.
#'     Defaults to `keepFeatures = FALSE` hence feature definitions are removed
#'     from the returned object by default.
#'
#' @param method For `filterChromPeaks`: `character(1)` allowing to specify the
#'     method by which chromatographic peaks should be filtered. Currently only
#'     `method = "keep"` is supported (i.e. specify with parameter `keep` which
#'     chromatographic peaks should be retained).
#'
#' @param msLevel. For `filterMz`, `filterRt`: `numeric`
#'     defining the MS level(s) to which operations should be applied or to
#'     which the object should be subsetted.
#'
#' @param mz For `filterMz`: `numeric(2)` defining the lower and upper
#'     mz value for the filtering.
#'
#' @param object A [XCMSnExp] object.
#'
#' @param rt For `filterRt`: `numeric(2)` defining the retention time
#'     window (lower and upper bound) for the filtering.
#'
#' @param x For `[` and `[[`: an `XCMSnExp` object.
#'
#' @param ... Optional additional arguments.
#'
#' @return All methods return an [XCMSnExp] object.
#'
#' @author Johannes Rainer
#'
#' @seealso [XCMSnExp] for base class documentation.
#'
#' @seealso [XChromatograms()] for similar filter functions on
#'     `XChromatograms` objects.
#'
#' @rdname XCMSnExp-filter-methods
#'
#' @md
#'
#' @examples
#'
#' ## Loading a test data set with identified chromatographic peaks
#' data(faahko_sub)
#' ## Update the path to the files for the local system
#' dirname(faahko_sub) <- system.file("cdf/KO", package = "faahKO")
#'
#' ## Disable parallel processing for this example
#' register(SerialParam())
#'
#' ## Subset the dataset to the first and third file.
#' xod_sub <- filterFile(faahko_sub, file = c(1, 3))
#'
#' ## The number of chromatographic peaks per file for the full object
#' table(chromPeaks(faahko_sub)[, "sample"])
#'
#' ## The number of chromatographic peaks per file for the subset
#' table(chromPeaks(xod_sub)[, "sample"])
#'
#' basename(fileNames(faahko_sub))
#' basename(fileNames(xod_sub))
#'
#' ## Filter on mz values; chromatographic peaks and features within the
#' ## mz range are retained (as well as adjusted retention times).
#' xod_sub <- filterMz(faahko_sub, mz = c(300, 400))
#' head(chromPeaks(xod_sub))
#' nrow(chromPeaks(xod_sub))
#' nrow(chromPeaks(faahko_sub))
#'
#' ## Filter on rt values. All chromatographic peaks and features within the
#' ## retention time range are retained. Filtering is performed by default on
#' ## adjusted retention times, if present.
#' xod_sub <- filterRt(faahko_sub, rt = c(2700, 2900))
#'
#' range(rtime(xod_sub))
#' head(chromPeaks(xod_sub))
#' range(chromPeaks(xod_sub)[, "rt"])
#'
#' nrow(chromPeaks(faahko_sub))
#' nrow(chromPeaks(xod_sub))
#'
#' ## Extract a single Spectrum
#' faahko_sub[[4]]
#'
#' ## Subsetting using [ removes all preprocessing results - using
#' ## keepAdjustedRtime = TRUE would keep adjusted retention times, if present.
#' xod_sub <- faahko_sub[fromFile(faahko_sub) == 1]
#' xod_sub
#'
#' ## Using split does also remove preprocessing results, but it supports the
#' ## optional parameter keepAdjustedRtime.
#' ## Split the object into a list of XCMSnExp objects, one per file
#' xod_list <- split(faahko_sub, f = fromFile(faahko_sub))
#' xod_list
setMethod("[", "XCMSnExp", function(x, i, j, ..., drop = TRUE) {
    if (!missing(j))
        stop("subsetting by columns ('j') not supported")
    if (missing(i))
        return(x)
    else if (!(is.numeric(i) | is.logical(i)))
        stop("'i' has to be either numeric or logical")
    ## Only result we might eventually keep is adjusted rtimes...
    newFd <- new("MsFeatureData")
    ph <- list()
    ## Check if we have keepAdjustedRtime as an additional parameter
    ## in ...
    keepAdjustedRtime <- list(...)$ke
    if (is.null(keepAdjustedRtime))
        keepAdjustedRtime <- FALSE
    if ((hasFeatures(x) || hasChromPeaks(x)) && length(i) > 30)
        warning("Removed preprocessing results")
    if (hasAdjustedRtime(x) && keepAdjustedRtime) {
        new_adj <- rtime(x, adjusted = TRUE)[i]
        adjustedRtime(newFd) <-
            unname(split(new_adj, f = fromFile(x)[i]))
        p <- processHistory(x, type = .PROCSTEP.RTIME.CORRECTION)
        ph <- p[length(ph)]
    }
    lockEnvironment(newFd, bindings = TRUE)
    x@msFeatureData <- newFd
    x@.processHistory <- ph
    callNextMethod()
})

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
setMethod("bin", "XCMSnExp", function(x, binSize = 1L, msLevel.) {
    if (hasAdjustedRtime(x) | hasFeatures(x) |
        hasChromPeaks(x)) {
        ## object@.processHistory <- list()
        ## object@msFeatureData <- new("MsFeatureData")
        x <- dropAdjustedRtime(x)
        x <- dropFeatureDefinitions(x)
        x <- dropChromPeaks(x)
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

    newMfd <- new("MsFeatureData")
    ph <- processHistory(object)
    if (hasAdjustedRtime(object)) {
        if (keepAdjustedRtime) {
            ## issue #210: keep adjusted retention times if wanted.
            keep_by_file <- base::split(keep_logical, as.factor(fromFile(object)))
            adj_rt <- base::mapply(FUN = function(y, z) {
                return(y[z])
            }, y = adjustedRtime(object, bySample = TRUE), z = keep_by_file,
            SIMPLIFY = FALSE)
            adjustedRtime(newMfd) <- adj_rt
        } else {
            object <- dropAdjustedRtime(object)
            ph <- dropProcessHistoriesList(ph,
                                           type = .PROCSTEP.RTIME.CORRECTION)
        }
    }
    if (hasChromPeaks(object)) {
        newMfd2 <- .filterChromPeaks(
            object@msFeatureData, which(chromPeakData(object)$ms_level %in%
                                                             msLevel.))
        if (hasChromPeaks(newMfd2)) {
            chromPeaks(newMfd) <- chromPeaks(newMfd2)
            chromPeakData(newMfd) <- chromPeakData(newMfd2)
        }
        if (hasFeatures(newMfd2))
            featureDefinitions(newMfd) <- featureDefinitions(newMfd2)
    }
    ## Subset processing history
    keep_ph <- vapply(ph, function(z) {
        if (inherits(z, "XProcessHistory")) {
            is_ok <- any(z@msLevel == msLevel.)
            if (is.na(is_ok) || is_ok) TRUE
            else FALSE
        } else TRUE
    }, logical(1))
    ph <- ph[keep_ph]
    ## Subsetting the object.
    tmp <- as(object, "OnDiskMSnExp")[base::which(keep_logical)]
    object <- as(tmp, "XCMSnExp")
    ## Put the stuff back
    object@processingData@processing <- c(object@processingData@processing, msg)
    lockEnvironment(newMfd, bindings = TRUE)
    object@msFeatureData <- newMfd
    object@.processHistory <- ph
    validObject(object)
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

#' @rdname XCMSnExp-filter-methods
setMethod(
    "filterFile", "XCMSnExp",
    function(object, file, keepAdjustedRtime = hasAdjustedRtime(object),
             keepFeatures = FALSE) {
        object <- .filter_file_XCMSnExp(object, file = file,
                                        keepAdjustedRtime = keepAdjustedRtime,
                                        keepFeatures = keepFeatures)
        validObject(object)
        object
    })

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
        pks <- chromPeaks(object)
        keepIdx <- which(pks[, "mzmin"] >= mz[1] & pks[, "mzmax"] <= mz[2])
        newE <- .filterChromPeaks(object@msFeatureData, idx = keepIdx)
        lockEnvironment(newE, bindings = TRUE)
        object@msFeatureData <- newE
    }
    validObject(object)
    object
})

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
            newMfd <- .filterChromPeaks(object@msFeatureData, idx = keep_fts)
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
        keep_by_file <- base::split(keep_logical, as.factor(fromFile(object)))
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
    object
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
        object@.processHistory <- list()
        object@msFeatureData <- new("MsFeatureData")
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
        object@.processHistory <- list()
        object@msFeatureData <- new("MsFeatureData")
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
        object@.processHistory <- list()
        object@msFeatureData <- new("MsFeatureData")
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
        x@.processHistory <- list()
        x@msFeatureData <- new("MsFeatureData")
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

#' @rdname XCMSnExp-peak-grouping-results
setMethod("quantify", "XCMSnExp", function(object, ...) {
    .XCMSnExp2SummarizedExperiment(object, ...)
})

#' @title Peak grouping/correspondence based on time dimension peak densities
#'
#' @description
#'
#' `groupChromPeaks,XCMSnExp,PeakDensityParam`:
#' performs correspondence (peak grouping within and across samples) within
#' in mz dimension overlapping slices of MS data based on the density
#' distribution of the identified chromatographic peaks in the slice along
#' the time axis.
#'
#' The correspondence analysis can be performed on chromatographic peaks of
#' any MS level (if present and if chromatographic peak detection has been
#' performed for that MS level) defining features combining these peaks. The
#' MS level can be selected with the parameter `msLevel`. By default, calling
#' `groupChromPeaks` will remove any previous correspondence results. This can
#' be disabled with `add = TRUE`, which will add newly defined features to
#' already present feature definitions.
#'
#' @param object For `groupChromPeaks`: an [XCMSnExp] object
#'     containing the results from a previous peak detection analysis (see
#'     [findChromPeaks()]).
#'
#'     For all other methods: a `PeakDensityParam` object.
#'
#' @param param A `PeakDensityParam` object containing all settings for
#'     the peak grouping algorithm.
#'
#' @param msLevel `integer(1)` (default `msLevel = 1L`) defining the MS level
#'     on which the correspondence should be performed. It is required that
#'     chromatographic peaks of the respective MS level are present.
#'
#' @param add `logical(1)` (default `add = FALSE`) allowing to perform an
#'     additional round of correspondence (e.g. on a different MS level) and
#'     add features to the already present feature definitions.
#'
#' @return
#'
#' For `groupChromPeaks`: a [XCMSnExp] object with the
#' results of the correspondence analysis. The definition of the resulting
#' mz-rt features can be accessed with the [featureDefinitions()] method
#'
#' @seealso
#'
#' [XCMSnExp] for the object containing the results of the correspondence.
#'
#' [plotChromPeakDensity()] for plotting chromatographic peak density with the
#' possibility to test different parameter settings.
#'
#' @md
#'
#' @rdname groupChromPeaks-density
setMethod("groupChromPeaks",
          signature(object = "XCMSnExp", param = "PeakDensityParam"),
          function(object, param, msLevel = 1L, add = FALSE) {
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
                           "samples!", call. = FALSE)
              }
              if (hasChromPeaks(object) && !.has_chrom_peak_data(object))
                  object <- updateObject(object)
              if (hasFeatures(object) &&
                  !any(colnames(featureDefinitions(object)) == "ms_level"))
                  object <- updateObject(object)
              startDate <- date()
              res <- do_groupChromPeaks_density(
                  chromPeaks(object, msLevel = msLevel),
                  sampleGroups = sampleGroups(param),
                  bw = bw(param),
                  minFraction = minFraction(param),
                  minSamples = minSamples(param),
                  binSize = binSize(param),
                  maxFeatures = maxFeatures(param))
              xph <- XProcessHistory(param = param, date. = startDate,
                                     type. = .PROCSTEP.PEAK.GROUPING,
                                     fileIndex = 1:length(fileNames(object)),
                                     msLevel = msLevel)
              object <- addProcessHistory(object, xph)
              ## Add the results.
              df <- DataFrame(res)
              if (!nrow(df)) {
                  warning("Unable to group any chromatographic peaks. ",
                          "You might have to adapt your settings.")
                  return(object)
              }
              df$ms_level <- as.integer(msLevel)
              if (!all(chromPeakData(object)$ms_level %in% msLevel))
                  df <- .update_feature_definitions(
                      df, rownames(chromPeaks(object, msLevel = msLevel)),
                      rownames(chromPeaks(object)))
              if (hasFeatures(object)) {
                  startFrom <- max(as.integer(
                      sub("FT", "", rownames(featureDefinitions(object))))) + 1
                  rownames(df) <- .featureIDs(nrow(df), from = startFrom)
                  df <- rbind(featureDefinitions(object), df)
              } else
                  rownames(df) <- .featureIDs(nrow(df))
              featureDefinitions(object) <- df
              validObject(object)
              object
          })


#' @title Single-spectrum non-chromatography MS data peak grouping
#'
#' @description
#'
#' `groupChromPeaks,XCMSnExp,MzClustParam`: performs high resolution peak
#' grouping for single spectrum metabolomics data.
#'
#' @note Calling `groupChromPeaks` on an `XCMSnExp` object will cause
#'     all eventually present previous correspondence results to be dropped.
#'
#' @param object For `groupChromPeaks`: an [XCMSnExp] object containing the
#'     results from a previous chromatographic peak detection analysis (see
#'     [findChromPeaks()]).
#'
#'     For all other methods: a `MzClustParam` object.
#'
#' @param param A `MzClustParam` object containing all settings for
#'     the peak grouping algorithm.
#'
#' @param msLevel `integer(1)` defining the MS level. Currently only MS level
#'     1 is supported.
#'
#' @return
#'
#' For `groupChromPeaks`: a [XCMSnExp] object with the results of the peak
#' grouping step (i.e. the features). These can be accessed with the
#' [featureDefinitions()] method.
#'
#' @seealso [XCMSnExp] for the object containing the results of
#'     the peak grouping.
#'
#' @md
#'
#' @rdname groupChromPeaks-mzClust
setMethod("groupChromPeaks",
          signature(object = "XCMSnExp", param = "MzClustParam"),
          function(object, param, msLevel = 1L) {
              if (!hasChromPeaks(object))
                  stop("No chromatographic peak detection results in 'object'! ",
                       "Please perform first a peak detection using the ",
                       "'findChromPeak' method.")
              if (any(msLevel != 1))
                  stop("Currently peak grouping is only supported for MS level 1")
              ## I'm expecting a single spectrum per file!
              rtL <- split(rtime(object), f = as.factor(fromFile(object)))
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
              if (hasChromPeaks(object) & !.has_chrom_peak_data(object))
                  object <- updateObject(object)
              startDate <- date()
              res <- do_groupPeaks_mzClust(chromPeaks(object, msLevel = msLevel),
                                           sampleGroups = sampleGroups(param),
                                           ppm = ppm(param),
                                           absMz = absMz(param),
                                           minFraction = minFraction(param),
                                           minSamples = minSamples(param))
              xph <- XProcessHistory(param = param, date. = startDate,
                                     type. = .PROCSTEP.PEAK.GROUPING,
                                     fileIndex = 1:length(fileNames(object)),
                                     msLevel = msLevel)
              object <- addProcessHistory(object, xph)
              ## Add the results.
              df <- DataFrame(res$featureDefinitions)
              df$peakidx <- res$peakIndex
              if (nrow(df) == 0)
                  stop("Unable to group any chromatographic peaks. You might ",
                       "have to adapt your settings.")
              if (!all(chromPeakData(object)$ms_level %in% msLevel))
                  df <- .update_feature_definitions(
                      df, rownames(chromPeaks(object, msLevel = msLevel)),
                      rownames(chromPeaks(object)))
              if (nrow(df) > 0)
                  rownames(df) <- .featureIDs(nrow(df))
              featureDefinitions(object) <- df
              validObject(object)
              object
          })


#' @title Peak grouping/correspondence based on proximity in the mz-rt space
#'
#' @description
#'
#' `groupChromPeaks,XCMSnExp,NearestPeaksParam`:
#' performs peak grouping based on the proximity between chromatographic
#' peaks from different samples in the mz-rt range.
#'
#' The correspondence analysis can be performed on chromatographic peaks of
#' any MS level (if present and if chromatographic peak detection has been
#' performed for that MS level) defining features combining these peaks. The
#' MS level can be selected with the parameter `msLevel`. By default, calling
#' `groupChromPeaks` will remove any previous correspondence results. This can
#' be disabled with `add = TRUE`, which will add newly defined features to
#' already present feature definitions.
#'
#' @param object For `groupChromPeaks`: an [XCMSnExp] object containing the
#'     results from a previous chromatographic peak detection
#'     analysis (see [findChromPeaks()]).
#'
#'     For all other methods: a `NearestPeaksParam` object.
#'
#' @param msLevel `integer(1)` defining the MS level on which the correspondence
#'     should be performed. It is required that chromatographic peaks of the
#'     respective MS level are present.
#'
#' @param add `logical(1)` (default `add = FALSE`) allowing to perform an
#'     additional round of correspondence (e.g. on a different MS level) and
#'     add features to the already present feature definitions.
#'
#' @param msLevel `integer(1)` defining the MS level. Currently only MS level
#'     1 is supported.
#'
#' @return
#'
#' For `groupChromPeaks`: a [XCMSnExp] object with the results of the peak
#' grouping/correspondence step (i.e. the mz-rt features). These can be
#' accessed with the [featureDefinitions()] method.
#'
#' @seealso [XCMSnExp] for the object containing the results of
#'     the peak grouping.
#'
#' @md
#'
#' @rdname groupChromPeaks-nearest
setMethod("groupChromPeaks",
          signature(object = "XCMSnExp", param = "NearestPeaksParam"),
          function(object, param, msLevel = 1L, add = FALSE) {
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
              if (hasChromPeaks(object) & !.has_chrom_peak_data(object))
                  object <- updateObject(object)
              if (hasFeatures(object) &&
                  !any(colnames(featureDefinitions(object)) == "ms_level"))
                  object <- updateObject(object)
              startDate <- date()
              res <- do_groupChromPeaks_nearest(
                  chromPeaks(object, msLevel = msLevel),
                  sampleGroups = sampleGroups(param),
                  mzVsRtBalance = mzVsRtBalance(param),
                  absMz = absMz(param),
                  absRt = absRt(param),
                  kNN = kNN(param))
              xph <- XProcessHistory(param = param, date. = startDate,
                                     type. = .PROCSTEP.PEAK.GROUPING,
                                     fileIndex = 1:length(fileNames(object)),
                                     msLevel = msLevel)
              object <- addProcessHistory(object, xph)
              ## Add the results.
              df <- DataFrame(res$featureDefinitions)
              if (!nrow(df)) {
                  warning("Unable to group any chromatographic peaks. ",
                          "You might have to adapt your settings.")
                  return(object)
              }
              df$peakidx <- res$peakIndex
              df$ms_level <- as.integer(msLevel)
              if (!all(chromPeakData(object)$ms_level %in% msLevel))
                  df <- .update_feature_definitions(
                      df, rownames(chromPeaks(object, msLevel = msLevel)),
                      rownames(chromPeaks(object)))
              if (hasFeatures(object)) {
                  startFrom <- max(as.integer(
                      sub("FT", "", rownames(featureDefinitions(object))))) + 1
                  rownames(df) <- .featureIDs(nrow(df), from = startFrom)
                  df <- rbind(featureDefinitions(object), df)
              } else
                  rownames(df) <- .featureIDs(nrow(df))
              featureDefinitions(object) <- df
              validObject(object)
              object
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
#' @param msLevel \code{integer(1)} specifying the MS level. Currently only MS
#'     level 1 is supported.
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
          function(object, param, msLevel = 1L) {
              if (hasAdjustedRtime(object)) {
                  message("Removing previous alignment results")
                  object <- dropAdjustedRtime(object)
              }
              if (any(msLevel != 1))
                  stop("Alignment is currently only supported for MS level 1")
              if (!hasChromPeaks(object))
                  stop("No chromatographic peak detection results in 'object'! ",
                       "Please perform first a peak detection using the ",
                       "'findChromPeaks' method.")
              if (!hasFeatures(object))
                  stop("No feature definitions found in 'object'! Please ",
                       "perform first a peak grouping using the ",
                       "'groupChromPeak' method.")
              if (hasChromPeaks(object) & !.has_chrom_peak_data(object))
                  object <- updateObject(object)
              startDate <- date()
              ## If param does contain a peakGroupsMatrix extract that one,
              ## otherwise generate it.
              if (nrow(peakGroupsMatrix(param)))
                  pkGrpMat <- peakGroupsMatrix(param)
              else
                  pkGrpMat <- adjustRtimePeakGroups(object, param = param)
              res <- do_adjustRtime_peakGroups(
                  chromPeaks(object, msLevel = msLevel),
                  peakIndex = .update_feature_definitions(
                      featureDefinitions(object), rownames(chromPeaks(object)),
                      rownames(chromPeaks(object, msLevel = msLevel)))$peakidx,
                  rtime = rtime(object, bySample = TRUE),
                  minFraction = minFraction(param),
                  extraPeaks = extraPeaks(param),
                  smooth = smooth(param),
                  span = span(param),
                  family = family(param),
                  peakGroupsMatrix = pkGrpMat,
                  subset = subset(param),
                  subsetAdjust = subsetAdjust(param)
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
                                     msLevel = msLevel)
              object <- addProcessHistory(object, xph)
              validObject(object)
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
              if (any(msLevel != 1))
                  stop("Alignment is currently only supported for MS level 1")
              if (hasChromPeaks(object) & !.has_chrom_peak_data(object))
                  object <- updateObject(object)
              startDate <- date()
              res <- adjustRtime(as(object, "OnDiskMSnExp"), param = param,
                                 msLevel = msLevel)
              ## res <- .obiwarp(as(object, "OnDiskMSnExp"), param = param)
              ## Dropping the feature groups.
              object <- dropFeatureDefinitions(object)
              ## Add the results. adjustedRtime<- should also fix the retention
              ## times for the peaks! Want to keep also the latest alignment
              ## information
              adjustedRtime(object) <- unname(
                  split(res, as.factor(fromFile(object))))
              ## Add the process history step.
              xph <- XProcessHistory(param = param, date. = startDate,
                                     type. = .PROCSTEP.RTIME.CORRECTION,
                                     fileIndex = 1:length(fileNames(object)),
                                     msLevel = msLevel)
              object <- addProcessHistory(object, xph)
              validObject(object)
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
    profMat(as(object, "OnDiskMSnExp"), method = method, step = step,
            baselevel = baselevel, basespace = basespace,
            mzrange. = mzrange., fileIndex = fileIndex, ...)
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
#' Parameter `msLevel` allows to choose a specific MS level for which feature
#' values should be returned (given that features have been defined for that MS
#' level).
#'
#' \code{quantify,XCMSnExp}: return the preprocessing results as an
#' \code{\link{SummarizedExperiment}} object containing the feature abundances
#' as assay matrix, the feature definitions (returned by
#' \code{\link{featureDefinitions}}) as \code{rowData} and the phenotype
#' information as \code{colData}. This is an ideal container for further
#' processing of the data. Internally, the \code{\link{featureValues}} method
#' is used to extract the feature abundances, parameters for that method can
#' be passed to \code{quantify} with \code{...}.
#'
#' @note
#'
#' This method is equivalent to the \code{\link{groupval}} for
#' \code{xcmsSet} objects. Note that \code{missing = 0} should be used to
#' get the same behaviour as \code{groupval}, i.e. report missing values as 0
#' after a call to \code{fillPeaks}.
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
#'     \code{chromPeaks(object)} that should be returned. Defaults to
#'     \code{"into"} in which case the integrated peak area is returned. To
#'     get the index of the peak in the \code{chromPeaks(object)} matrix use
#'     \code{"index"}.
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
#' @param missing how missing values should be reported. Allowed values are
#'     \code{NA} (the default), a \code{numeric} or
#'     \code{missing = "rowmin_half"}. The latter replaces any \code{NA} with
#'     half of the row's minimal (non-missing) value.
#'
#' @param msLevel for `featureValues`: `integer` defining the MS level(s) for
#'     which feature values should be returned. By default, values for features
#'     defined for all MS levels are returned.
#'
#' @param ... For \code{quantify}: additional parameters to be passed on to the
#'     \code{\link{featureValues}} method.
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
#' For \code{quantify}: a \code{\link{SummarizedExperiment}} representing
#' the preprocessing results.
#'
#' @author Johannes Rainer
#'
#' @seealso
#' \code{\link{XCMSnExp}} for information on the data object.
#'
#' \code{\link{featureDefinitions}} to extract the \code{DataFrame} with the
#' feature definitions.
#'
#' \code{\link{featureChromatograms}} to extract ion chromatograms for each
#' feature.
#'
#' \code{\link{hasFeatures}} to evaluate whether the
#' \code{\link{XCMSnExp}} provides feature definitions.
#'
#' \code{\link{groupval}} for the equivalent method on \code{xcmsSet} objects.
#'
#' @rdname XCMSnExp-peak-grouping-results
setMethod("featureValues", "XCMSnExp", function(object, method = c("medret",
                                                                   "maxint",
                                                                   "sum"),
                                                value = "into",
                                                intensity = "into",
                                                filled = TRUE, missing = NA,
                                                msLevel = integer()) {
    ## Input argument checkings
    if (!hasFeatures(object, msLevel = msLevel))
        stop("No feature definitions for MS level(s) ", msLevel,
             " present. Call 'groupChromPeaks' first.")
    method <- match.arg(method)
    if (method == "sum" & !(value %in% c("into", "maxo")))
        stop("method 'sum' is only allowed if value is set to 'into'",
             " or 'maxo'")
    if (is.character(missing)) {
        if (!(missing %in% c("rowmin_half")))
            stop("if 'missing' is not 'NA' or a numeric it should",
                 " be one of: \"rowmin_half\".")
    } else {
        if (!is.numeric(missing) & !is.na(missing))
            stop("'missing' should be either 'NA', a numeric or one",
                 " of: \"rowmin_half\".")
    }
    fNames <- basename(fileNames(object))
    pks <- chromPeaks(object)
    ## issue #157: replace all values for filled-in peaks with NA
    if (!filled)
        pks[chromPeakData(object)$is_filled, ] <- NA
    .feature_values(
        pks = pks, fts = featureDefinitions(object, msLevel = msLevel),
        method = method, value = value, intensity = intensity,
        colnames = fNames, missing = missing)
})

#' Internal function to extract feature values based on featureDefinitions
#' `fts` and chromPeaks `pks`.
#'
#' @author Johannes Rainer
#'
#' @noRd
.feature_values <- function(pks, fts, method, value = "into",
                            intensity = "into", colnames,
                            missing = NA) {
    ftIdx <- fts$peakidx
    ## Match columns
    idx_rt <- match("rt", colnames(pks))
    idx_int <- match(intensity, colnames(pks))
    idx_samp <- match("sample", colnames(pks))
    vals <- matrix(nrow = length(ftIdx), ncol = length(colnames))
    nSamples <- seq_along(colnames)
    if (method == "sum") {
        for (i in seq_along(ftIdx)) {
            cur_pks <- pks[ftIdx[[i]], c(value, "sample"), drop=FALSE]
            int_sum <- split(cur_pks[, value],
                             as.factor(as.integer(cur_pks[, "sample"])))
            vals[i, as.numeric(names(int_sum))] <-
                unlist(lapply(int_sum, base::sum), use.names = FALSE)
        }
    } else {
        if (method == "medret") {
            medret <- fts$rtmed
            for (i in seq_along(ftIdx)) {
                gidx <- ftIdx[[i]][
                    base::order(base::abs(pks[ftIdx[[i]],
                                              idx_rt] - medret[i]))]
                vals[i, ] <- gidx[
                    base::match(nSamples, pks[gidx, idx_samp])]
            }
        }
        if (method == "maxint") {
            for (i in seq_along(ftIdx)) {
                gidx <- ftIdx[[i]][
                    base::order(pks[ftIdx[[i]], idx_int],
                                decreasing = TRUE)]
                vals[i, ] <- gidx[base::match(nSamples,
                                              pks[gidx, idx_samp])]
            }
        }
        if (value != "index") {
            if (!any(colnames(pks) == value))
                stop("Column '", value, "' not present in the ",
                     "chromatographic peaks matrix!")
            vals <- pks[vals, value]
            dim(vals) <- c(length(ftIdx), length(nSamples))
        }
    }
    if (value != "index") {
        if (is.numeric(missing)) {
            vals[is.na(vals)] <- missing
        }
        if (!is.na(missing) & missing == "rowmin_half") {
            for (i in seq_len(nrow(vals))) {
                nas <- is.na(vals[i, ])
                if (any(nas))
                    vals[i, nas] <- min(vals[i, ], na.rm = TRUE) / 2
            }
        }
    }
    colnames(vals) <- colnames
    rownames(vals) <- rownames(fts)
    vals
}

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
#' `chromatogram`: extract chromatographic data (such as an extracted ion
#' chromatogram, a base peak chromatogram or total ion chromatogram) from
#' an [OnDiskMSnExp] or [XCMSnExp] objects. See also the help page of the
#' `chromatogram` function in the `MSnbase` package.
#'
#' @details
#'
#' Arguments `rt` and `mz` allow to specify the MS data slice (i.e. the m/z
#' range and retention time window) from which the chromatogram should be
#' extracted. These parameters can be either a `numeric` of length 2 with the
#' lower and upper limit, or a `matrix` with two columns with the lower and
#' upper limits to extract multiple EICs at once.
#' The parameter `aggregationSum` allows to specify the function to be
#' used to aggregate the intensities across the m/z range for the same
#' retention time. Setting `aggregationFun = "sum"` would e.g. allow
#' to calculate the **total ion chromatogram** (TIC),
#' `aggregationFun = "max"` the **base peak chromatogram** (BPC).
#'
#' If for a given retention time no intensity is measured in that spectrum a
#' `NA` intensity value is returned by default. This can be changed with the
#' parameter `missing`, setting `missing = 0` would result in a `0` intensity
#' being returned in these cases.
#'
#' @note
#'
#' For [XCMSnExp] objects, if adjusted retention times are
#' available, the `chromatogram` method will by default report
#' and use these (for the subsetting based on the provided parameter
#' `rt`). This can be changed by setting `adjustedRtime = FALSE`.
#'
#' @param object Either a [OnDiskMSnExp] or [XCMSnExp] object from which the
#'     chromatograms should be extracted.
#'
#' @param rt `numeric(2)` or two-column `matrix` defining the lower
#'     and upper boundary for the retention time range(s). If not specified,
#'     the full retention time range of the original data will be used.
#'
#' @param mz `numeric(2)` or two-column `matrix` defining the lower
#'     and upper mz value for the MS data slice(s). If not specified, the
#'     chromatograms will be calculated on the full mz range.
#'
#' @param adjustedRtime For `chromatogram,XCMSnExp`: whether the
#'     adjusted (`adjustedRtime = TRUE`) or raw retention times
#'     (`adjustedRtime = FALSE`) should be used for filtering and returned
#'     in the resulting [MChromatograms] object. Adjusted
#'     retention times are used by default if available.
#'
#' @param aggregationFun `character(1)` specifying the function to be used to
#'     aggregate intensity values across the mz value range for the same
#'     retention time. Allowed values are `"sum"` (the default), `"max"`,
#'     `"mean"` and `"min"`.
#'
#' @param missing `numeric(1)` allowing to specify the intensity value to
#'     be used if for a given retention time no signal was measured within the
#'     mz range of the corresponding scan. Defaults to `NA_real_` (see also
#'     Details and Notes sections below). Use `missing = 0` to resemble the
#'     behaviour of the `getEIC` from the *old* user interface.
#'
#' @param msLevel `integer(1)` specifying the MS level from which the
#'     chromatogram should be extracted. Defaults to `msLevel = 1L`.
#'
#' @param BPPARAM Parallelisation backend to be used, which will
#'     depend on the architecture. Default is
#'     `BiocParallel::bparam()`.
#'
#' @param filled `logical(1)` whether filled-in peaks should also be
#'     returned. Defaults to `filled = FALSE`, i.e. returns only detected
#'     chromatographic peaks in the result object.
#'
#' @param include `character(1)` defining which chromatographic peaks should be
#'     returned. Supported are `include = "apex_within"` (the default) which
#'     returns chromatographic peaks that have their apex within the `mz` `rt`
#'     range, `include = "any"` to return all chromatographic peaks which
#'     m/z and rt ranges overlap the `mz` and `rt` or `include = "none"` to
#'     not include any chromatographic peaks.
#'
#' @return
#'
#' `chromatogram` returns a [XChromatograms] object with
#' the number of columns corresponding to the number of files in
#' `object` and number of rows the number of specified ranges (i.e.
#' number of rows of matrices provided with arguments `mz` and/or
#' `rt`). All chromatographic peaks with their apex position within the
#' m/z and retention time range are also retained as well as all feature
#' definitions for these peaks.
#'
#' @author Johannes Rainer
#'
#' @seealso [XCMSnExp] for the data object.
#'     [Chromatogram] for the object representing chromatographic data.
#'
#'     [XChromatograms] for the object allowing to arrange
#'     multiple [XChromatogram] objects.
#'
#'     [plot] to plot a [XChromatogram] or [MChromatograms] objects.
#'
#'     `as` (`as(x, "data.frame")`) in `MSnbase` for a method to extract
#'     the MS data as `data.frame`.
#'
#' @export
#'
#' @md
#'
#' @rdname chromatogram-method
#'
#' @examples
#'
#' ## Load a test data set with identified chromatographic peaks
#' data(faahko_sub)
#' ## Update the path to the files for the local system
#' dirname(faahko_sub) <- system.file("cdf/KO", package = "faahKO")
#'
#' ## Disable parallel processing for this example
#' register(SerialParam())
#'
#' ## Extract the ion chromatogram for one chromatographic peak in the data.
#' chrs <- chromatogram(faahko_sub, rt = c(2700, 2900), mz = 335)
#'
#' chrs
#'
#' ## Identified chromatographic peaks
#' chromPeaks(chrs)
#'
#' ## Plot the chromatogram
#' plot(chrs)
#'
#' ## Extract chromatograms for multiple ranges.
#' mzr <- matrix(c(335, 335, 344, 344), ncol = 2, byrow = TRUE)
#' rtr <- matrix(c(2700, 2900, 2600, 2750), ncol = 2, byrow = TRUE)
#' chrs <- chromatogram(faahko_sub, mz = mzr, rt = rtr)
#'
#' chromPeaks(chrs)
#'
#' plot(chrs)
#'
#' ## Get access to all chromatograms for the second mz/rt range
#' chrs[1, ]
#'
#' ## Plot just that one
#' plot(chrs[1, , drop = FALSE])
setMethod(
    "chromatogram", "XCMSnExp",
    function(object, rt, mz, aggregationFun = "sum", missing = NA_real_,
             msLevel = 1L, BPPARAM = bpparam(),
             adjustedRtime = hasAdjustedRtime(object), filled = FALSE,
             include = c("apex_within", "any", "none")) {
        include <- match.arg(include)
        if (adjustedRtime)
            adj_rt <- rtime(object, adjusted = TRUE)
        object_od <- as(object, "OnDiskMSnExp")
        fcs <- c("fileIdx", "spIdx", "seqNum", "acquisitionNum", "msLevel",
                 "polarity", "retentionTime", "precursorScanNum")
        fcs <- intersect(fcs, colnames(.fdata(object)))
        object_od <- selectFeatureData(object_od, fcol = fcs)
        if (adjustedRtime)
            object_od@featureData$retentionTime <- adj_rt
        res <- MSnbase::chromatogram(object_od, rt = rt, mz = mz,
                                     aggregationFun = aggregationFun,
                                     missing = missing, msLevel = msLevel,
                                     BPPARAM = BPPARAM)
        if (!hasChromPeaks(object) | include == "none")
            return(res)
        ## Process peaks
        lvls <- 1:length(fileNames(object))
        if (missing(rt))
            rt <- c(-Inf, Inf)
        if (missing(mz))
            mz <- c(-Inf, Inf)
        if (is.matrix(rt) | is.matrix(mz)) {
            ## Ensure rt and mz are aligned.
            if (!is.matrix(rt))
                rt <- matrix(rt, ncol = 2)
            if (!is.matrix(mz))
                mz <- matrix(mz, ncol = 2)
            if (nrow(rt) == 1)
                rt <- matrix(rep(rt, nrow(mz)), ncol = 2, byrow = TRUE)
            if (nrow(mz) == 1)
                mz <- matrix(rep(mz, nrow(rt)), ncol = 2, byrow = TRUE)
            pk_list <- vector("list", nrow(mz))
            pkd_list <- vector("list", nrow(mz))
            for (i in 1:nrow(mz)) {
                pks <- chromPeaks(object, rt = rt[i, ], mz = mz[i, ],
                                  type = include)
                pkd <- chromPeakData(object)[rownames(pks), , drop = FALSE]
                if (!filled) {
                    pks <- pks[!pkd$is_filled, , drop = FALSE]
                    pkd <- extractROWS(pkd, which(!pkd$is_filled))
                }
                smpls <- factor(pks[, "sample"], levels = lvls)
                pk_list[[i]] <- split.data.frame(pks, smpls)
                pkd_list[[i]] <- split.data.frame(pkd, smpls)
            }
            pks <- do.call(rbind, pk_list)
            pks <- pks[seq_along(pks)]
            pkd <- do.call(rbind, pkd_list)
            pkd <- pkd[seq_along(pkd)]
        } else {
            pks <- chromPeaks(object, rt = rt, mz = mz,
                              type = include)
            pkd <- chromPeakData(object)[rownames(pks), , drop = FALSE]
            if (!filled) {
                pks <- pks[!pkd$is_filled, , drop = FALSE]
                pkd <- extractROWS(pkd, which(!pkd$is_filled))
            }
            smpls <- factor(pks[, "sample"], levels = lvls)
            pks <- split.data.frame(pks, smpls)
            pkd <- split.data.frame(pkd, smpls)
        }
        res <- as(res, "XChromatograms")
        res@.Data <- matrix(
            mapply(unlist(res), pks, pkd, FUN = function(chr, pk, pd) {
                chr@chromPeaks <- pk
                chr@chromPeakData <- pd
                chr
            }), nrow = nrow(res), dimnames = dimnames(res))
        res@.processHistory <- object@.processHistory
        if (hasFeatures(object)) {
            pks_sub <- chromPeaks(res)
            ## Loop through each EIC "row" to ensure all features in
            ## that EIC are retained.
            fts <- lapply(seq_len(nrow(res)), function(r) {
                fdev <- featureDefinitions(object, mz = mz(res)[r, ],
                                           rt = rt)
                if (nrow(fdev)) {
                    fdev$row <- r
                    .subset_features_on_chrom_peaks(
                        fdev, chromPeaks(object), pks_sub)
                } else DataFrame()
            })
            res@featureDefinitions <- do.call(rbind, fts)
        }
        validObject(res)
        res
    })

#' @rdname XCMSnExp-class
#'
#' @aliases faahko_sub
#'
#' @description
#'
#' \code{findChromPeaks} performs chromatographic peak detection
#' on the provided \code{XCMSnExp} objects. For more details see the method
#' for \code{\linkS4class{XCMSnExp}}.
#' Note that by default (with parameter \code{add = FALSE}) previous peak
#' detection results are removed. Use \code{add = TRUE} to perform a second
#' round of peak detection and add the newly identified peaks to the previous
#' peak detection results. Correspondence results (features) are always removed
#' prior to peak detection. Previous alignment (retention
#' time adjustment) results are kept, i.e. chromatographic peak detection
#' is performed using adjusted retention times if the data was first
#' aligned using e.g. obiwarp (\code{\link{adjustRtime-obiwarp}}).
#'
#' @param param A \code{\link{CentWaveParam}}, \code{\link{MatchedFilterParam}},
#'     \code{\link{MassifquantParam}}, \code{\link{MSWParam}} or
#'     \code{\link{CentWavePredIsoParam}} object with the settings for the
#'     chromatographic peak detection algorithm.
#'
#' @param add For \code{findChromPeaks}: if newly identified chromatographic
#'     peaks should be added to the peak matrix with the already identified
#'     chromatographic peaks. By default (\code{add = FALSE}) previous
#'     peak detection results will be removed.
#'
#' @inheritParams findChromPeaks-centWave
setMethod("findChromPeaks",
          signature(object = "XCMSnExp", param = "Param"),
          function(object, param, BPPARAM = bpparam(),
                   return.type = "XCMSnExp", msLevel = 1L, add = FALSE) {
              ## Remove previous correspondence results.
              if (hasFeatures(object)) {
                  message("Removed feature definitions.")
                  object <- dropFeatureDefinitions(
                      object,
                      keepAdjustedRtime = hasAdjustedRtime(object))
              }
              ## Remove previous chromatographic peaks.
              has_peaks <- hasChromPeaks(object)
              if (has_peaks & !add) {
                  message("Removed previously identified chromatographic peaks.")
                  object <- dropChromPeaks(
                      object,
                      keepAdjustedRtime = hasAdjustedRtime(object))
              }
              if (add && has_peaks) {
                  old_cp <- chromPeaks(object)
                  old_cpd <- chromPeakData(object)
              }
              meth <- selectMethod("findChromPeaks",
                                   signature = c(object = "OnDiskMSnExp",
                                                 param = class(param)))
              object <- do.call(meth, args = list(object = object,
                                                  param = param,
                                                  BPPARAM = BPPARAM,
                                                  return.type = return.type,
                                                  msLevel = msLevel))
              if (add && has_peaks) {
                  old_cp <- rbindFill(old_cp, chromPeaks(object))
                  old_cpd <- rbindFill(old_cpd, chromPeakData(object))
                  old_hist <- object@.processHistory
                  chromPeaks(object) <- old_cp
                  rownames(old_cpd) <- rownames(chromPeaks(object))
                  chromPeakData(object) <- old_cpd
                  object@.processHistory <- old_hist
              }
              ## object@.processHistory <- list()
              validObject(object)
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
#' feature was identified and add it to the [chromPeaks()] matrix. Such
#' *filled-in* peaks are indicated with a `TRUE` in column `"is_filled"` in
#' the result object's [chromPeakData()] data frame.
#'
#' Two different gap-filling approaches are implemented:
#'
#' - `param = FillChromPeaksParam()`: the default of the original `xcms` code.
#'   Signal is integrated from the m/z and retention time range as defined in
#'   the [featureDefinitions()] data frame, i.e. from the `"rtmin"`, `"rtmax"`,
#'   `"mzmin"` and `"mzmax"`. See details below for more information and
#'   settings for this method.
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
#' the signal in the mz-rt region of the feature. Two different approaches
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
#' Note that no peak is filled in if no signal was present in a file/sample
#' in the respective mz-rt area. These samples will still show a `NA`
#' in the matrix returned by the [featureValues()] method.
#'
#' @param object `XCMSnExp` object with identified and grouped chromatographic
#'     peaks.
#'
#' @param param `FillChromPeaksParam` or `ChromPeakAreaParam` object
#'     defining which approach should be used (see details section).
#'
#' @param mzmin `function` to be applied to values in the `"mzmin"` column of all
#'     chromatographic peaks of a feature to define the lower m/z value of the
#'     area from which signal for the feature should be integrated. Defaults to
#'     `mzmin = function(z) quantile(z, probs = 0.25)` hence using the 25%
#'     quantile of all values.
#'
#' @param mzmax `function` to be applied to values in the `"mzmax"` column of all
#'     chromatographic peaks of a feature to define the upper m/z value of the
#'     area from which signal for the feature should be integrated. Defaults to
#'     `mzmax = function(z) quantile(z, probs = 0.75)` hence using the 75%
#'     quantile of all values.
#'
#' @param rtmin `function` to be applied to values in the `"rtmin"` column of all
#'     chromatographic peaks of a feature to define the lower rt value of the
#'     area from which signal for the feature should be integrated. Defaults to
#'     `rtmin = function(z) quantile(z, probs = 0.25)` hence using the 25%
#'     quantile of all values.
#'
#' @param rtmax `function` to be applied to values in the `"rtmax"` column of all
#'     chromatographic peaks of a feature to define the upper rt value of the
#'     area from which signal for the feature should be integrated. Defaults to
#'     `rtmax = function(z) quantile(z, probs = 0.75)` hence using the 75%
#'     quantile of all values.
#'
#' @param expandMz for `FillChromPeaksParam`: `numeric(1)` defining the value
#'     by which the mz width of peaks should be expanded. Each peak is expanded
#'     in mz direction by `expandMz *` their original m/z width. A value of
#'     `0` means no expansion, a value of `1` grows each peak by `1 *` the m/z
#'     width of the peak resulting in peaks with twice their original size in
#'     m/z direction (expansion by half m/z width to both sides).
#'
#' @param expandRt for `FillChromPeaksParam`: `numeric(1)`, same as `expandMz`
#'     but for the retention time width.
#'
#' @param ppm for `FillChromPeaksParam`: `numeric(1)` optionally specifying a
#'     *ppm* by which the m/z width of the peak region should be expanded. For
#'     peaks with an m/z width smaller than `mean(c(mzmin, mzmax)) * ppm / 1e6`,
#'     the `mzmin` will be replaced by
#'     `mean(c(mzmin, mzmax)) - (mean(c(mzmin, mzmax)) * ppm / 2 / 1e6)`
#'     `mzmax` by
#'     `mean(c(mzmin, mzmax)) + (mean(c(mzmin, mzmax)) * ppm / 2 / 1e6)`.
#'     This is applied before eventually expanding the m/z width using the
#'     `expandMz` parameter.
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
#'     on one MS level at a time is supported, to fill in peaks for MS level 1
#'     and 2 run first using `msLevel = 1` and then (on the returned
#'     result object) again with `msLevel = 2`.
#'
#' @param BPPARAM Parallel processing settings.
#'
#' @return
#'
#' A `XCMSnExp` object with previously missing chromatographic peaks for
#' features filled into its [chromPeaks()] matrix.
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
setMethod("fillChromPeaks",
          signature(object = "XCMSnExp", param = "FillChromPeaksParam"),
          function(object, param, msLevel = 1L, BPPARAM = bpparam()) {
              if (length(msLevel) != 1)
                  stop("Can only perform peak filling for one MS level at a time")
              if (!hasFeatures(object, msLevel = msLevel))
                  stop("No feature definitions for MS level ", msLevel,
                       " present. Please run 'groupChromPeaks' first.")
              if (.hasFilledPeaks(object))
                  message("Filled peaks already present, adding still missing",
                          " peaks.")
              if (hasChromPeaks(object) & !.has_chrom_peak_data(object))
                  object <- updateObject(object)
              startDate <- date()
              expandMz <- expandMz(param)
              expandRt <- expandRt(param)
              fixedMz <- fixedMz(param)
              fixedRt <- fixedRt(param)
              ppm <- ppm(param)
              message("Defining peak areas for filling-in .",
                      appendLF = FALSE)
              ## Define or extend the peak area from which the signal should be
              ## extracted.
              ## Original code: use the median of the min/max rt and mz per peak.
              fdef <- featureDefinitions(object, msLevel = msLevel)
              aggFunLow <- median
              aggFunHigh <- median
              ## Note: we ensure in the downstream function that the rt range is
              ## within the rt range. For the mz range it doesn't matter.
              tmp_pks <- chromPeaks(object)[, c("rtmin", "rtmax", "mzmin",
                                                "mzmax")]
              pkArea <- do.call(
                  rbind,
                  lapply(
                      fdef$peakidx, function(z) {
                          pa <- c(aggFunLow(tmp_pks[z, 1]),
                                  aggFunHigh(tmp_pks[z, 2]),
                                  aggFunLow(tmp_pks[z, 3]),
                                  aggFunHigh(tmp_pks[z, 4]))
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
                          pa
                      }
                  ))
              rm(tmp_pks)
              message(".", appendLF = FALSE)
              colnames(pkArea) <- c("rtmin", "rtmax", "mzmin", "mzmax")
              ## Add mzmed column - needed for MSW peak filling.
              pkArea <- cbind(group_idx = 1:nrow(pkArea), pkArea,
                              mzmed = as.numeric(fdef$mzmed))
              pkGrpVal <- featureValues(object, value = "index",
                                        msLevel = msLevel)
              message(".", appendLF = FALSE)
              ## Check if there is anything to fill...
              if (!any(is.na(rowSums(pkGrpVal)))) {
                  message("No missing peaks present.")
                  return(object)
              }
              message(".", appendLF = FALSE)
              ## Split the object by file and define the peaks for which
              objectL <- vector("list", length(fileNames(object)))
              pkAreaL <- objectL
              ## We need "only" a list of OnDiskMSnExp, one for each file but
              ## instead of filtering by file we create small objects to keep
              ## memory requirement to a minimum.
              req_fcol <- requiredFvarLabels("OnDiskMSnExp")
              min_fdata <- .fdata(object)[, req_fcol]
              rt_range <- range(pkArea[, c("rtmin", "rtmax")])
              if (hasAdjustedRtime(object))
                  min_fdata$retentionTime <- adjustedRtime(object)
              for (i in 1:length(fileNames(object))) {
                  fd <- min_fdata[min_fdata$fileIdx == i, ]
                  fd$fileIdx <- 1L
                  objectL[[i]] <- new(
                      "OnDiskMSnExp",
                      processingData = new("MSnProcess",
                                           files = fileNames(object)[i]),
                      featureData = new("AnnotatedDataFrame", fd),
                      phenoData = new("NAnnotatedDataFrame",
                                      data.frame(sampleNames = "1")),
                      experimentData = new("MIAPE",
                                           instrumentManufacturer = "a",
                                           instrumentModel = "a",
                                           ionSource = "a",
                                           analyser = "a",
                                           detectorType = "a"))
                  ## Want to extract intensities only for peaks that were not
                  ## found in a sample.
                  pkAreaL[[i]] <- pkArea[is.na(pkGrpVal[, i]), , drop = FALSE]
              }
              rm(pkGrpVal)
              rm(pkArea)
              rm(min_fdata)
              message(" OK\nStart integrating peak areas from original files")
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
              cp_colnames <- colnames(chromPeaks(object))
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
                                      cn = cp_colnames),
                                  BPPARAM = BPPARAM, SIMPLIFY = FALSE)
              } else if (findPeakMethod == "matchedFilter") {
                  res <- bpmapply(FUN = .getChromPeakData_matchedFilter,
                                  objectL, pkAreaL, as.list(1:length(objectL)),
                                  MoreArgs = list(cn = cp_colnames,
                                                  param = prm,
                                                  msLevel = msLevel),
                                  BPPARAM = BPPARAM, SIMPLIFY = FALSE)
              } else {
                  res <- bpmapply(FUN = .getChromPeakData, objectL,
                                  pkAreaL, as.list(1:length(objectL)),
                                  MoreArgs = list(cn = cp_colnames,
                                                  mzCenterFun = mzCenterFun,
                                                  msLevel = msLevel),
                                  BPPARAM = BPPARAM, SIMPLIFY = FALSE)
              }
              rm(objectL)

              res <- do.call(rbind, res)
              ## cbind the group_idx column to track the feature/peak group.
              res <- cbind(
                  res, group_idx = unlist(lapply(pkAreaL,
                                                 function(z) z[, "group_idx"]),
                                          use.names = FALSE))
              ## Remove those without a signal
              res <- res[!is.na(res[, "into"]), , drop = FALSE]
              if (nrow(res) == 0) {
                  warning("Could not integrate any signal for the missing ",
                          "peaks! Consider increasing 'expandMz' and 'expandRt'.")
                  return(object)
              }
              ## Intermediate cleanup of objects.
              rm(pkAreaL)
              gc()

              ## Get the msFeatureData:
              newFd <- new("MsFeatureData")
              newFd@.xData <- .copy_env(object@msFeatureData)
              object@msFeatureData <- new("MsFeatureData")
              incr <- nrow(chromPeaks(newFd))
              for (i in unique(res[, "group_idx"])) {
                  fdef$peakidx[[i]] <- c(fdef$peakidx[[i]],
                  (which(res[, "group_idx"] == i) + incr))
              }
              ## Combine feature data with those from other MS levels
              fdef <- rbind(
                  fdef, featureDefinitions(newFd)[
                           featureDefinitions(newFd)$ms_level != msLevel, ,
                           drop = FALSE])
              if (!any(colnames(fdef) == "ms_level"))
                  fdef$ms_level <- 1L
              else
                  fdef <- fdef[order(fdef$ms_level), ]
              ## Define IDs for the new peaks; include fix for issue #347
              maxId <- max(as.numeric(
                  sub("M", "", sub("^CP", "", rownames(chromPeaks(newFd))))))
              if (maxId < 1)
                  stop("chromPeaks matrix lacks rownames; please update ",
                       "'object' with the 'updateObject' function.")
              toId <- maxId + nrow(res)
              rownames(res) <- sprintf(
                  paste0("CP", "%0", ceiling(log10(toId + 1L)), "d"),
                  (maxId + 1L):toId)
              chromPeaks(newFd) <- rbind(chromPeaks(newFd),
                                         res[, -ncol(res)])
              cpd <- extractROWS(chromPeakData(newFd), rep(1L, nrow(res)))
              cpd[,] <- NA
              cpd$ms_level <- as.integer(msLevel)
              cpd$is_filled <- TRUE
              if (!any(colnames(chromPeakData(newFd)) == "is_filled"))
                  chromPeakData(newFd)$is_filled <- FALSE
              chromPeakData(newFd) <- rbind(chromPeakData(newFd), cpd)
              rownames(chromPeakData(newFd)) <- rownames(chromPeaks(newFd))
              featureDefinitions(newFd) <- fdef
              lockEnvironment(newFd, bindings = TRUE)
              object@msFeatureData <- newFd
              ## Add a process history step
              ph <- XProcessHistory(param = param,
                                    date. = startDate,
                                    type. = .PROCSTEP.PEAK.FILLING,
                                    fileIndex = 1:length(fileNames(object)),
                                    msLevel = msLevel)
              object <- addProcessHistory(object, ph) ## this also validates object.
              object
          })


#' @rdname fillChromPeaks
setMethod("fillChromPeaks",
          signature(object = "XCMSnExp", param = "ChromPeakAreaParam"),
          function(object, param, msLevel = 1L, BPPARAM = bpparam()) {
              if (length(msLevel) != 1)
                  stop("Can only perform peak filling for one MS level at a time")
              if (!hasFeatures(object, msLevel = msLevel))
                  stop("No feature definitions for MS level ", msLevel,
                       " present. Please run 'groupChromPeaks' first.")
              if (.hasFilledPeaks(object))
                  message("Filled peaks already present, adding still missing",
                          " peaks.")
              if (hasChromPeaks(object) & !.has_chrom_peak_data(object))
                  object <- updateObject(object)
              startDate <- date()
              message("Defining peak areas for filling-in .",
                      appendLF = FALSE)
              fts_region <- .features_ms_region(
                  object, mzmin = param@mzmin, mzmax = param@mzmax,
                  rtmin = param@rtmin, rtmax = param@rtmax, msLevel = msLevel)
              fts_region <- cbind(group_idx = seq_len(nrow(fts_region)),
                                  fts_region,
                                  mzmed = featureDefinitions(object)$mzmed)
              message(".", appendLF = FALSE)
              pk_idx <- featureValues(object, value = "index",
                                      msLevel = msLevel)
              message(".", appendLF = FALSE)
              ## Check if there is anything to fill...
              if (!any(is.na(rowSums(pk_idx)))) {
                  message("No missing peaks present.")
                  return(object)
              }
              ## Split the object by file and define the peaks for which
              objectL <- vector("list", length(fileNames(object)))
              pkAreaL <- objectL
              ## We need "only" a list of OnDiskMSnExp, one for each file but
              ## instead of filtering by file we create small objects to keep
              ## memory requirement to a minimum.
              req_fcol <- requiredFvarLabels("OnDiskMSnExp")
              min_fdata <- .fdata(object)[, req_fcol]
              if (hasAdjustedRtime(object))
                  min_fdata$retentionTime <- adjustedRtime(object)
              for (i in 1:length(fileNames(object))) {
                  fd <- min_fdata[min_fdata$fileIdx == i, ]
                  fd$fileIdx <- 1L
                  objectL[[i]] <- new(
                      "OnDiskMSnExp",
                      processingData = new("MSnProcess",
                                           files = fileNames(object)[i]),
                      featureData = new("AnnotatedDataFrame", fd),
                      phenoData = new("NAnnotatedDataFrame",
                                      data.frame(sampleNames = "1")),
                      experimentData = new("MIAPE",
                                           instrumentManufacturer = "a",
                                           instrumentModel = "a",
                                           ionSource = "a",
                                           analyser = "a",
                                           detectorType = "a"))
                  ## Want to extract intensities only for peaks that were not
                  ## found in a sample.
                  pkAreaL[[i]] <- fts_region[is.na(pk_idx[, i]), , drop = FALSE]
              }
              rm(pk_idx)
              rm(fts_region)
              rm(min_fdata)
              message(" OK\nStart integrating peak areas from original files")
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
              cp_colnames <- colnames(chromPeaks(object))
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
                                      cn = cp_colnames),
                                  BPPARAM = BPPARAM, SIMPLIFY = FALSE)
              } else if (findPeakMethod == "matchedFilter") {
                  res <- bpmapply(FUN = .getChromPeakData_matchedFilter,
                                  objectL, pkAreaL, as.list(1:length(objectL)),
                                  MoreArgs = list(cn = cp_colnames,
                                                  param = prm,
                                                  msLevel = msLevel),
                                  BPPARAM = BPPARAM, SIMPLIFY = FALSE)
              } else {
                  res <- bpmapply(FUN = .getChromPeakData, objectL,
                                  pkAreaL, as.list(1:length(objectL)),
                                  MoreArgs = list(cn = cp_colnames,
                                                  mzCenterFun = mzCenterFun,
                                                  msLevel = msLevel),
                                  BPPARAM = BPPARAM, SIMPLIFY = FALSE)
              }
              rm(objectL)

              res <- do.call(rbind, res)
              ## cbind the group_idx column to track the feature/peak group.
              res <- cbind(
                  res, group_idx = unlist(lapply(pkAreaL,
                                                 function(z) z[, "group_idx"]),
                                          use.names = FALSE))
              ## Remove those without a signal
              res <- res[!is.na(res[, "into"]), , drop = FALSE]
              if (nrow(res) == 0) {
                  warning("Could not integrate any signal for the missing ",
                          "peaks!")
                  return(object)
              }
              ## Intermediate cleanup of objects.
              rm(pkAreaL)

              ## Get the msFeatureData:
              newFd <- new("MsFeatureData")
              newFd@.xData <- .copy_env(object@msFeatureData)
              object@msFeatureData <- new("MsFeatureData")
              incr <- nrow(chromPeaks(newFd))
              fdef <- featureDefinitions(newFd, msLevel = msLevel)
              for (i in unique(res[, "group_idx"])) {
                  fdef$peakidx[[i]] <- c(fdef$peakidx[[i]],
                  (which(res[, "group_idx"] == i) + incr))
              }
              ## Combine feature data with those from other MS levels
              fdef <- rbind(fdef,
                            extractROWS(
                                featureDefinitions(newFd),
                                which(featureDefinitions(newFd)$ms_level != msLevel)))
              if (!any(colnames(fdef) == "ms_level"))
                  fdef$ms_level <- 1L
              else
                  fdef <- extractROWS(fdef, order(fdef$ms_level))
              ## Define IDs for the new peaks; include fix for issue #347
              maxId <- max(as.numeric(
                  sub("M", "", sub("^CP", "", rownames(chromPeaks(newFd))))))
              if (maxId < 1)
                  stop("chromPeaks matrix lacks rownames; please update ",
                       "'object' with the 'updateObject' function.")
              toId <- maxId + nrow(res)
              rownames(res) <- sprintf(
                  paste0("CP", "%0", ceiling(log10(toId + 1L)), "d"),
                  (maxId + 1L):toId)
              chromPeaks(newFd) <- rbind(chromPeaks(newFd),
                                         res[, -ncol(res)])
              cpd <- extractROWS(chromPeakData(newFd), rep(1L, nrow(res)))
              cpd[,] <- NA
              cpd$ms_level <- as.integer(msLevel)
              cpd$is_filled <- TRUE
              if (!any(colnames(chromPeakData(newFd)) == "is_filled"))
                  chromPeakData(newFd)$is_filled <- FALSE
              chromPeakData(newFd) <- rbind(chromPeakData(newFd), cpd)
              rownames(chromPeakData(newFd)) <- rownames(chromPeaks(newFd))
              featureDefinitions(newFd) <- fdef
              lockEnvironment(newFd, bindings = TRUE)
              object@msFeatureData <- newFd
              ## Add a process history step
              ph <- XProcessHistory(param = param,
                                    date. = startDate,
                                    type. = .PROCSTEP.PEAK.FILLING,
                                    fileIndex = 1:length(fileNames(object)),
                                    msLevel = msLevel)
              object <- addProcessHistory(object, ph) ## this also validates object.
              object
          })



#' @rdname fillChromPeaks
setMethod(
    "fillChromPeaks",
    signature(object = "XCMSnExp", param = "missing"),
              function(object,
                       param,
                       BPPARAM = bpparam(),
                       msLevel = 1L) {
                  fillChromPeaks(object, param = FillChromPeaksParam(),
                                 BPPARAM = BPPARAM, msLevel = msLevel)
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
    keep_pks <- which(!chromPeakData(object)$is_filled)
    object@msFeatureData <- .filterChromPeaks(object@msFeatureData, keep_pks)
    object <- dropProcessHistories(object, type = .PROCSTEP.PEAK.FILLING)
    validObject(object)
    object
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
#'
#' ## Load a test data set with detected peaks
#' data(faahko_sub)
#' ## Update the path to the files for the local system
#' dirname(faahko_sub) <- system.file("cdf/KO", package = "faahKO")
#'
#' ## Disable parallel processing for this example
#' register(SerialParam())
#'
#' ## Extract the full MS data for a certain retention time range
#' ## as a data.frame
#' tmp <- filterRt(faahko_sub, rt = c(2800, 2900))
#' ms_all <- as(tmp, "data.frame")
#' head(ms_all)
#' nrow(ms_all)
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
    validObject(object)
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

#' @title Export MS data to mzML/mzXML files
#'
#' @description
#'
#' `writeMSData` exports mass spectrometry data in mzML or mzXML format.
#' If adjusted retention times are present, these are used as retention time of
#' the exported spectra.
#'
#' @param object [XCMSnExp] object with the mass spectrometry data.
#'
#' @param file `character` with the file name(s). The length of this parameter
#'     has to match the number of files/samples of `object`.
#'
#' @param outformat `character(1)` defining the format of the output files (
#'     either `"mzml"` or `"mzxml"`).
#'
#' @param copy `logical(1)` if metadata (data processing, software used,
#'     original file names etc) should be copied from the original files.
#'
#' @param software_processing optionally provide specific data processing steps.
#'     See documentation of the `software_processing` parameter of
#'     [mzR::writeMSData()].
#'
#' @param ... Additional parameters to pass down to the [writeMSData()]
#'     function in the `MSnbase` package, such as `outformat` to specify the
#'     output format (`"mzml"` or `"mzxml"`) or `copy` to specify whether
#'     general information from the original MS data files (such as data
#'     processing, software etc) should be copied to the new files.
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @seealso [writeMSData()] function in the `MSnbase` package.
setMethod("writeMSData", signature(object = "XCMSnExp", file = "character"),
          function(object, file, outformat = c("mzml", "mzxml"),
                   copy = FALSE, software_processing = NULL, ...) {
              if (hasAdjustedRtime(object)) {
                  object <- applyAdjustedRtime(object)
                  ## Define the software processing.
                  software_processing <- c(
                      software_processing,
                      list(c("xcms", paste(packageVersion("xcms"),
                                           collapse = "."),
                             "MS:1001582", "MS:1000745")))
              }
              callNextMethod(object = object, file = file,
                             outformat = outformat, copy = copy,
                             software_processing = software_processing, ...)
          })

#' @title Plot chromatographic peak density along the retention time axis
#'
#' @aliases plotChromPeakDensity
#'
#' @description
#'
#' Plot the density of chromatographic peaks along the retention
#' time axis and indicate which peaks would be (or were) grouped into the
#' same feature based using the *peak density* correspondence method.
#' Settings for the *peak density* method can be passed with an
#' [PeakDensityParam] object to parameter `param`. If the `object` contains
#' correspondence results and the correspondence was performed with the
#' *peak groups* method, the results from that correspondence can be
#' visualized setting `simulate = FALSE`.
#'
#' @details
#'
#' The `plotChromPeakDensity` function allows to evaluate
#' different settings for the *peak density* on an mz slice of
#' interest (e.g. containing chromatographic peaks corresponding to a known
#' metabolite).
#' The plot shows the individual peaks that were detected within the
#' specified `mz` slice at their retention time (x-axis) and sample in
#' which they were detected (y-axis). The density function is plotted as a
#' black line. Parameters for the `density` function are taken from the
#' `param` object. Grey rectangles indicate which chromatographic peaks
#' would be grouped into a feature by the `peak density` correspondence
#' method. Parameters for the algorithm are also taken from `param`.
#' See [groupChromPeaks-density()] for more information about the
#' algorithm and its supported settings.
#'
#' @param object A [XCMSnExp] object with identified
#'     chromatographic peaks.
#'
#' @param mz `numeric(2)` defining an mz range for which the peak density
#'     should be plotted.
#'
#' @param rt `numeric(2)` defining an optional rt range for which the
#'     peak density should be plotted. Defaults to the absolute retention time
#'     range of `object`.
#'
#' @param param [PeakDensityParam] from which parameters for the
#'     *peak density* correspondence algorithm can be extracted. If not provided
#'     and if `object` contains feature definitions with the correspondence/
#'     peak grouping being performed by the *peak density* method, the
#'     corresponding parameter class stored in `object` is used.
#'
#' @param simulate `logical(1)` defining whether correspondence should be
#'     simulated within the specified m/z / rt region or (with
#'     `simulate = FALSE`) whether the results from an already performed
#'     correspondence should be shown.
#'
#' @param col Color to be used for the individual samples. Length has to be 1
#'     or equal to the number of samples in `object`.
#'
#' @param xlab `character(1)` with the label for the x-axis.
#'
#' @param ylab `character(1)` with the label for the y-axis.
#'
#' @param xlim `numeric(2)` representing the limits for the x-axis.
#'     Defaults to the range of the `rt` parameter.
#'
#' @param main `character(1)` defining the title of the plot. By default
#'     (for `main = NULL`) the mz-range is used.
#'
#' @param type `character(1)` specifying how peaks are called to be located
#'     within the region defined by `mz` and `rt`. Can be one of `"any"`,
#'     `"within"`, and `"apex_within"` for all peaks that are even partially
#'     overlapping the region, peaks that are completely within the region, and
#'     peaks for which the apex is within the region. This parameter is passed
#'     to the [chromPeaks] function. See related documentation for more
#'     information and examples.
#'
#' @param ... Additional parameters to be passed to the `plot` function. Data
#'     point specific parameters such as `bg` or `pch` have to be of length 1
#'     or equal to the number of samples.
#'
#' @return The function is called for its side effect, i.e. to create a plot.
#'
#' @author Johannes Rainer
#'
#' @seealso [groupChromPeaks-density()] for details on the
#'     *peak density* correspondence method and supported settings.
#'
#' @md
#'
#' @rdname plotChromPeakDensity
#'
#' @examples
#'
#' ## Load a test data set with detected peaks
#' data(faahko_sub)
#' ## Update the path to the files for the local system
#' dirname(faahko_sub) <- system.file("cdf/KO", package = "faahKO")
#'
#' ## Plot the chromatographic peak density for a specific mz range to evaluate
#' ## different peak density correspondence settings.
#' mzr <- c(305.05, 305.15)
#'
#' plotChromPeakDensity(faahko_sub, mz = mzr, pch = 16,
#'     param = PeakDensityParam(sampleGroups = rep(1, length(fileNames(faahko_sub)))))
#'
setMethod("plotChromPeakDensity", "XCMSnExp", .plotChromPeakDensity)

setMethod("updateObject", "XCMSnExp", function(object) {
    newFd <- new("MsFeatureData")
    newFd@.xData <- .copy_env(object@msFeatureData)
    if (hasChromPeaks(newFd)) {
        if (is.null(rownames(chromPeaks(newFd))))
            rownames(chromPeaks(newFd)) <-
                .featureIDs(nrow(chromPeaks(newFd)), "CP")
        if (!.has_chrom_peak_data(newFd)) {
            newFd$chromPeakData <- DataFrame(
                ms_level = rep(1L, nrow(chromPeaks(newFd))),
                row.names = rownames(chromPeaks(newFd)))
            if (any(colnames(chromPeaks(newFd)) == "is_filled")) {
                newFd$chromPeakData$is_filled <- as.logical(
                    chromPeaks(newFd)[, "is_filled"])
                newFd$chromPeaks <-
                    newFd$chromPeaks[, colnames(newFd$chromPeaks) != "is_filled"]
            } else
                newFd$chromPeakData$is_filled <- FALSE
        }
        if (hasFeatures(newFd) &&
            !any(colnames(featureDefinitions(newFd)) == "ms_level"))
            newFd$featureDefinitions$ms_level <- 1L
    }
    lockEnvironment(newFd, bindings = TRUE)
    object@msFeatureData <- newFd
    if (!length(object@.processHistory))
        object@.processHistory <- list()
    validObject(object)
    object
})

#' @rdname XCMSnExp-class
setMethod("chromPeakData", "XCMSnExp", function(object) {
    chromPeakData(object@msFeatureData)
})
#' @rdname XCMSnExp-class
setReplaceMethod("chromPeakData", "XCMSnExp", function(object, value) {
    newFd <- new("MsFeatureData")
    newFd@.xData <- .copy_env(object@msFeatureData)
    chromPeakData(newFd) <- value
    lockEnvironment(newFd, bindings = TRUE)
    object@msFeatureData <- newFd
    validObject(object)
    object
})

#' @description
#'
#' \code{plot} plots the spectrum data (see \code{\link{plot}} for
#' \code{\link{MSnExp}} objects in the \code{MSnbase} package for more details.
#' For \code{type = "XIC"}, identified chromatographic peaks will be indicated
#' as rectangles with border color \code{peakCol}.
#'
#' @param x For \code{plot}: \code{XCMSnExp} object.
#'
#' @param y For \code{plot}: not used.
#'
#' @param peakCol For \code{plot}: the color that should be used to indicate
#'     identified chromatographic peaks (only in combination with
#'     \code{type = "XIC"} and if chromatographic peaks are present).
#'
#' @rdname XCMSnExp-class
setMethod("plot", c("XCMSnExp", "missing"),
          function(x, y, type = c("spectra", "XIC"),
                   peakCol = "#ff000060", ...) {
              type <- match.arg(type)
              if (type == "spectra" || !hasChromPeaks(x))
                  callNextMethod(x = x, type = type, ...)
              else .plot_XIC(x, peakCol = peakCol, ...)
          })

#' @title Remove chromatographic peaks with too large rt width
#'
#' @aliases refineChromPeaks CleanPeaksParam-class show,CleanPeaksParam-method
#'
#' @description
#'
#' Remove chromatographic peaks with a retention time range larger than the
#' provided maximal acceptable width (`maxPeakwidth`).
#'
#' @note
#'
#' `refineChromPeaks` methods will always remove feature definitions, because
#' a call to this method can change or remove identified chromatographic peaks,
#' which may be part of features.
#'
#' @param maxPeakwidth for `CleanPeaksParam`: `numeric(1)` defining the maximal
#'     allowed peak width (in retention time).
#'
#' @param msLevel `integer` defining for which MS level(s) the chromatographic
#'     peaks should be cleaned.
#'
#' @param object [XCMSnExp] object with identified chromatographic peaks.
#'
#' @param param `CleanPeaksParam` object defining the settings for the method.
#'
#' @return `XCMSnExp` object with chromatographic peaks exceeding the specified
#'     maximal retention time width being removed.
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @family chromatographic peak refinement methods
#'
#' @rdname refineChromPeaks-clean
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
#' ## Distribution of chromatographic peak widths
#' quantile(chromPeaks(faahko_sub)[, "rtmax"] - chromPeaks(faahko_sub)[, "rtmin"])
#'
#' ## Remove all chromatographic peaks with a width larger 60 seconds
#' data <- refineChromPeaks(faahko_sub, param = CleanPeaksParam(60))
#'
#' quantile(chromPeaks(data)[, "rtmax"] - chromPeaks(data)[, "rtmin"])
setMethod("refineChromPeaks", c(object = "XCMSnExp", param = "CleanPeaksParam"),
          function(object, param = CleanPeaksParam(),
                   msLevel = 1L) {
              if (!hasChromPeaks(object, msLevel = msLevel)) {
                  warning("No chromatographic peaks present in for MS level ",
                          msLevel, ". Please run 'findChromPeaks' first.")
                  return(object)
              }
              if (hasFeatures(object)) {
                  message("Removing feature definitions.")
                  object <- dropFeatureDefinitions(object)
              }
              validObject(param)
              rtwidths <- chromPeaks(object)[, "rtmax"] -
                  chromPeaks(object)[, "rtmin"]
              sel_ms <- chromPeakData(object)$ms_level %in% msLevel
              sel_rt <- rtwidths < param@maxPeakwidth & sel_ms
              keep <- which(sel_rt | !sel_ms)
              message("Removed ", nrow(chromPeaks(object)) - length(keep),
                      " of ", nrow(chromPeaks(object)),
                      " chromatographic peaks.")
              msf <- new("MsFeatureData")
              msf@.xData <- .copy_env(object@msFeatureData)
              chromPeaks(msf) <- chromPeaks(object)[keep, , drop = FALSE]
              chromPeakData(msf) <- extractROWS(chromPeakData(object), keep)
              object@msFeatureData <- msf
              ph <- processHistory(object, type = .PROCSTEP.PEAK.DETECTION)
              xph <- XProcessHistory(param = param, date. = date(),
                                     type. = .PROCSTEP.PEAK.REFINEMENT,
                                     fileIndex = 1:length(fileNames(object)),
                                     msLevel = msLevel)
              object <- addProcessHistory(object, xph)
              validObject(object)
              object
          })

#' @title Merge neighboring and overlapping chromatographic peaks
#'
#' @aliases MergeNeighboringPeaksParam-class show,MergeNeighboringPeaksParam-method
#'
#' @description
#'
#' Peak detection sometimes fails to identify a chromatographic peak correctly,
#' especially for broad peaks and if the peak shape is irregular (mostly for
#' HILIC data). In such cases several smaller peaks are reported. Also, peak
#' detection can result in partially or completely overlapping peaks. To reduce
#' such peak detection artifacts, this function merges chromatographic peaks
#' which are overlapping or close in rt and m/z dimension considering also the
#' measured signal intensities in the region between them.
#'
#' Chromatographic peaks are first expanded in m/z and retention time dimension
#' (based on parameters `expandMz`, `ppm` and `expandRt`) and subsequently
#' grouped into sets of merge candidates if they are (after expansion)
#' overlapping in both m/z and rt (within the same sample).
#' Candidate peaks are merged if the average intensity of the 3 data
#' points in the middle position between them (i.e. at half the distance between
#' `"rtmax"` of the first and `"rtmin"` of the second peak) is larger than a
#' certain proportion (`minProp`) of the smaller maximal intensity (`"maxo"`)
#' of both peaks. In cases in which this calculated mid point is **not**
#' located between the apexes of the two peaks (e.g. if the peaks are largely
#' overlapping) the average signal intensity at half way between the apexes is
#' used instead. Candidate peaks are not joined if all 3 data points between
#' them have `NA` intensities.
#' The joined peaks get the `"mz"`, `"rt"`, `"sn"` and `"maxo"` values from
#' the peak with the largest signal (`"maxo"`) as well as its row in the
#' metadata data frame of the peak (`chromPeakData`). The `"rtmin"`, `"rtmax"`
#' of the merged peaks are updated and `"into"` is recalculated based on all
#' the signal between `"rtmin"` and `"rtmax"` of the new merged peak. See
#' details for information on the `"mzmin"` and `"mzmax"` values of the merged
#' peak.
#'
#' @note
#'
#' Note that **each** peak gets expanded by `expandMz` and `expandRt`, thus
#' peaks differing by `2 * expandMz` (or `expandRt`) will be identified as
#' *overlapping*. As an example: m/z max of one peak is 12.2, m/z min of
#' another one is 12.4, if `expandMz = 0.1` the m/z max of the first peak
#' will be 12.3 and the m/z min of the second one 12.3, thus both are
#' considered overlapping.
#'
#' `refineChromPeaks` methods will always remove feature definitions, because
#' a call to this method can change or remove identified chromatographic peaks,
#' which may be part of features.
#'
#' Merging of chromatographic peaks is performed along the retention time axis,
#' i.e. candidate peaks are first ordered by their `"rtmin"` value. The signals
#' at half way between the first and the second candidate peak are then compared
#' to the smallest `"maxo"` of both and the two peaks are then merged if the
#' average signal between the peaks is larger `minProp`. For merging any
#' additional peak in a candidate peak list the `"maxo"` of that peak and the
#' newly merged peak are considered.
#'
#' @details
#'
#' For each set of candidate peaks an ion chromatogram is
#' extracted using the range of retention times and m/z values of these peaks.
#' The m/z range for the extracted ion chromatogram is expanded by `expandMz`
#' and `ppm` (on both sides) to reduce the possibility of missing signal
#' intensities between candidate peaks (variance of measured m/z values for
#' lower intensities is larger than for higher intensities and thus data points
#' not being part of identified chromatographic peaks tend to have m/z values
#' outside of the m/z range of the candidate peaks - especially for ToF
#' instruments). This also ensures that all data points from the same ion are
#' considered for the peak integration of merged peaks. The smallest and largest
#' m/z value of all data points used in the peak integration of the merged peak
#' are used as the merged peak's m/z range (i.e. columns `"mzmin"` and `"mzmax"`).
#'
#' @param expandRt `numeric(1)` defining by how many seconds the retention time
#'     window is expanded on both sides to check for overlapping peaks.
#'
#' @param expandMz `numeric(1)` constant value by which the m/z range of each
#'     chromatographic peak is expanded (on both sides!) to check for
#'     overlapping peaks.
#'
#' @param ppm `numeric(1)` defining a m/z relative value (in parts per million)
#'     by which the m/z range of each chromatographic peak is expanded
#'     to check for overlapping peaks.
#'
#' @param minProp `numeric(1)` between `0` and `1` representing the proporion
#'     of intensity to be required for peaks to be joined. See description for
#'     more details. The default (`minProp = 0.75`) means that peaks are only
#'     joined if the signal half way between then is larger 75% of the smallest
#'     of the two peak's `"maxo"` (maximal intensity at peak apex).
#'
#' @param msLevel `integer` defining for which MS level(s) the chromatographic
#'     peaks should be merged.
#'
#' @param object [XCMSnExp] object with identified chromatographic peaks.
#'
#' @param param `MergeNeighboringPeaksParam` object defining the settings for
#'     the method.
#'
#' @param BPPARAM parameter object to set up parallel processing. Uses the
#'     default parallel processing setup returned by `bpparam()`. See
#'     [bpparam()] for details and examples.
#'
#' @return `XCMSnExp` object with chromatographic peaks matching the defined
#'     conditions being merged.
#'
#' @author Johannes Rainer, Mar Garcia-Aloy
#'
#' @md
#'
#' @family chromatographic peak refinement methods
#'
#' @rdname refineChromPeaks-merge
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
setMethod("refineChromPeaks", c(object = "XCMSnExp",
                                param = "MergeNeighboringPeaksParam"),
          function(object, param = MergeNeighboringPeaksParam(),
                   msLevel = 1L, BPPARAM = bpparam()) {
              if (!hasChromPeaks(object, msLevel = msLevel)) {
                  warning("No chromatographic peaks present in for MS level ",
                          msLevel, ". Please run 'findChromPeaks' first.")
                  return(object)
              }
              if (hasFeatures(object)) {
                  message("Removing feature definitions.")
                  object <- dropFeatureDefinitions(object)
              }
              validObject(param)
              peak_count <- nrow(chromPeaks(object))
              res <- bplapply(.split_by_file2(object, msLevel. = msLevel,
                                              to_class = "XCMSnExp",
                                              subsetFeatureData = TRUE,
                                              keep_sample_idx = TRUE),
                              FUN = .merge_neighboring_peaks,
                              expandRt = param@expandRt,
                              expandMz = param@expandMz, ppm = param@ppm,
                              minProp = param@minProp,
                              BPPARAM = BPPARAM)
              pks <- do.call(rbind, lapply(res, "[[", 1))
              pkd <- do.call(rbind, lapply(res, "[[", 2))
              ## Add also peaks for other MS levels!
              other_msl <- !(chromPeakData(object)$ms_level %in% msLevel)
              if (any(other_msl)) {
                  pks <- rbind(pks, chromPeaks(object)[other_msl, , drop = FALSE])
                  pkd <- rbind(pkd, extractROWS(chromPeakData(object), which(other_msl)))
              }
              which_new <- is.na(rownames(pks))
              pkd$merged <- which_new
              max_id <- max(as.numeric(sub("CP", "", rownames(pks))),
                            na.rm = TRUE)
              if (!is.finite(max_id))
                  max_id <- 0
              rownames(pks)[which_new] <- .featureIDs(sum(which_new),
                                                      prefix = "CPM",
                                                      from = max_id + 1)
              rownames(pkd) <- rownames(pks)
              message("Merging reduced ", peak_count, " chromPeaks to ",
                      nrow(pks), ".")
              msf <- new("MsFeatureData")
              msf@.xData <- .copy_env(object@msFeatureData)
              chromPeaks(msf) <- pks
              chromPeakData(msf) <- pkd
              object@msFeatureData <- msf
              ph <- processHistory(object, type = .PROCSTEP.PEAK.DETECTION)
              xph <- XProcessHistory(param = param, date. = date(),
                                     type. = .PROCSTEP.PEAK.REFINEMENT,
                                     fileIndex = 1:length(fileNames(object)),
                                     msLevel = msLevel)
              object <- addProcessHistory(object, xph)
              validObject(object)
              object
          })

#' @title Remove chromatographic peaks based on intensity
#'
#' @aliases FilterIntensityParam-class show,FilterIntensityParam-method
#'
#' @description
#'
#' Remove chromatographic peaks with intensities below the specified threshold.
#' By default, with `nValues = 1`, all peaks with an intensity
#' `>= threshold` are retained. Parameter `value` allows to specify the column of
#' the [chromPeaks()] matrix that should be used for the filtering (defaults to
#' `value = "maxo"` and thus evaluating the maximal intensity for each peak).
#' With `nValues > 1` it is possible to keep only peaks that have `nValues`
#' intensities `>= threshold`. Note that this requires data import from the
#' original MS files and run time of the call can thus be significantly larger.
#' Also, for `nValues > 1` parameter `value` is ignored.
#'
#' @param threshold `numeric(1)` defining the minimal required intensity for
#'     a peak to be retained. Defaults to `threshold = 0`.
#'
#' @param nValues `integer(1)` defining the number of data points (per
#'     chromatographic peak) that have to be `>= threshold`. Defaults to
#'     `nValues = 1`.
#'
#' @param value `character(1)` specifying the column in [chromPeaks()] that
#'     should be used for the comparison. This is ignored for `nValues > 1`.
#'
#' @param msLevel `integer(1)` defining the MS level in which peaks should be
#'     filtered.
#'
#' @param object [XCMSnExp] object with identified chromatographic peaks.
#'
#' @param param `FilterIntensityParam` object defining the settings for
#'     the method.
#'
#' @param BPPARAM parameter object to set up parallel processing. Uses the
#'     default parallel processing setup returned by `bpparam()`. See
#'     [bpparam()] for details and examples.
#'
#' @return `XCMSnExp` object with filtererd chromatographic peaks.
#'
#' @author Johannes Rainer, Mar Garcia-Aloy
#'
#' @md
#'
#' @family chromatographic peak refinement methods
#'
#' @rdname refineChromPeaks-filter-intensity
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
#' ## Remove all peaks with a maximal intensity below 50000
#' res <- refineChromPeaks(faahko_sub, param = FilterIntensityParam(threshold = 50000))
#'
#' nrow(chromPeaks(faahko_sub))
#' nrow(chromPeaks(res))
#'
#' all(chromPeaks(res)[, "maxo"] > 50000)
#'
#' ## Keep only chromatographic peaks that have 3 signals above 20000; we
#' ## perform this on the data of a single file.
#' xdata <- filterFile(faahko_sub)
#'
#' res <- refineChromPeaks(xdata, FilterIntensityParam(threshold = 20000, nValues = 3))
#' nrow(chromPeaks(xdata))
#' nrow(chromPeaks(res))
setMethod("refineChromPeaks", c(object = "XCMSnExp",
                                param = "FilterIntensityParam"),
          function(object, param = FilterIntensityParam(),
                   msLevel = 1L, BPPARAM = bpparam()) {
              if (!hasChromPeaks(object, msLevel = msLevel)) {
                  warning("No chromatographic peaks present in for MS level ",
                          msLevel, ". Please run 'findChromPeaks' first.")
                  return(object)
              }
              if (hasFeatures(object)) {
                  message("Removing feature definitions.")
                  object <- dropFeatureDefinitions(object)
              }
              validObject(param)
              peak_count <- nrow(chromPeaks(object))
              if (param@nValues == 1) {
                  ## Simple subsetting of the chromPeaks matrix.
                  if (!any(colnames(chromPeaks(object)) %in% param@value))
                      stop("Column '", value, "' not found in chromPeaks matrix")
                  keep <- chromPeaks(object)[, param@value] >= param@threshold |
                      !chromPeakData(object)$ms_level %in% msLevel
              } else {
                  res <- bplapply(.split_by_file2(object, to_class = "XCMSnExp",
                                                  msLevel. = 1:10),
                                  FUN = .chrom_peaks_above_threshold,
                                  nValues = param@nValues,
                                  threshold = param@threshold,
                                  msLevel = msLevel, BPPARAM = BPPARAM)
                  keep <- unlist(res, use.names = FALSE)
                  if (length(keep) != nrow(chromPeaks(object)))
                      stop("Length of variable 'keep' does not match number ",
                           "of peaks. Please contact developers.")
              }
              msfd <- .filterChromPeaks(object, idx = which(keep))
              object@msFeatureData <- msfd
              message("Removed ", peak_count - nrow(chromPeaks(object)),
                      " chromatographic peaks.")
              ph <- processHistory(object, type = .PROCSTEP.PEAK.DETECTION)
              xph <- XProcessHistory(param = param, date. = date(),
                                     type. = .PROCSTEP.PEAK.REFINEMENT,
                                     fileIndex = 1:length(fileNames(object)),
                                     msLevel = msLevel)
              object <- addProcessHistory(object, xph)
              validObject(object)
              object
          })

#' @rdname XCMSnExp-filter-methods
setMethod("filterChromPeaks", "XCMSnExp",
          function(object, keep = rep(TRUE, nrow(chromPeaks(object))),
                   method = "keep", ...) {
              method <- match.arg(method)
              object <- switch(
                  method,
                  keep = {
                      idx <- .i2index(keep, ids = rownames(chromPeaks(object)),
                                      name = "keep")
                      object@msFeatureData <- .filterChromPeaks(object, idx)
                      object
                  })
              validObject(object)
              object
          })
