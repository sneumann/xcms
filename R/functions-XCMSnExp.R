#' @include DataClasses.R functions-utils.R

#' @description
#'
#' Takes a XCMSnExp and drops ProcessHistory steps from the
#' `@.processHistory` slot matching the provided type.
#'
#' @param num which should be dropped? If \code{-1} all matching will be dropped,
#'     otherwise just the most recent num.
#'
#' @return The XCMSnExp input object with selected ProcessHistory steps dropped.
#'
#' @noRd
dropProcessHistories <- function(x, type, num = -1) {
    x@.processHistory <- dropProcessHistoriesList(processHistory(x),
                                                  type = type, num = num)
    return(x)
}

#' Drop a generic processing history step based on the provided function name
#'
#' @noRd
dropGenericProcessHistory <- function(x, fun) {
    x@.processHistory <- dropGenericProcessHistoryList(processHistory(x), fun)
    x
}

#' Convert an XCMSnExp to an xcmsSet.
#'
#' @noRd
.XCMSnExp2xcmsSet <- function(from) {
    if (any(msLevel(from) > 1))
        stop("Coercing an XCMSnExp with MS level > 1 is not yet supported!")
    xs <- new("xcmsSet")
    ## @peaks <- chromPeaks
    if (hasChromPeaks(from))
        xs@peaks <- chromPeaks(from)
    ## @groups <- part of featureDefinitions
    ## @groupidx <- featureDefinitions(x)$peakidx
    if (hasFeatures(from)){
        fgs <- featureDefinitions(from)
        xs@groups <- S4Vectors::as.matrix(fgs[,names(fgs)!="peakidx"])
        rownames(xs@groups) <- NULL
        xs@groupidx <- fgs$peakidx
    }
    ## @rt combination from rtime(x) and adjustedRtime(x)
    rts <- list()
    ## Ensure we're getting the raw rt
    rts$raw <- rtime(from, bySample = TRUE, adjusted = FALSE)
    if (hasAdjustedRtime(from))
        rts$corrected <- adjustedRtime(from, bySample = TRUE)
    else
        rts$corrected <- rts$raw
    xs@rt <- rts

    ## @phenoData
    pd <- pData(from)
    if (nrow(pd) != length(fileNames(from))) {
        pd <- data.frame(file_name = basename(fileNames(from)))
        rownames(pd) <- pd$file_name
    }
    xs@phenoData <- pd
    ## @filepaths
    xs@filepaths <- fileNames(from)

    ## @profinfo (list)
    profMethod <- "bin"
    profStep <- 0.1
    profParam <- list()
    ## If we've got any MatchedFilterParam we can take the values from there
    ph <- processHistory(from, type = .PROCSTEP.PEAK.DETECTION)
    if (length(ph)) {
        if (is(ph[[1]], "XProcessHistory")) {
            prm <- processParam(ph[[1]])
            if (is(prm, "MatchedFilterParam")) {
                if (impute(prm) != "none") {
                    profMethod <- paste0("bin", impute(prm))
                    profStep <- binSize(prm)
                    if (profMethod == "binlinbase")
                        warning("Optional parameters for 'binlinbase' can not be",
                                " converted yet, using default parameters.")
                ## distance
                ## baseValue
                }
            }
        }
    }
    profinfo(xs) <- c(list(method = profMethod, step = profStep), profParam)

    ## @mslevel <- msLevel?
    xs@mslevel <- unique(msLevel(from))

    ## @scanrange
    xs@scanrange <- range(scanIndex(from))

    ## .processHistory: just take the processHistory as is.
    xs@.processHistory <- processHistory(from)

    ## Implement later or skip:
    ## @polarity (character), has to be "positive" or "negative"; not clear how
    ## and if this is used at all.

    ## @filled ... not yet.
    if (any(chromPeakData(from)$is_filled)) {
        fld <- which(chromPeakData(from)$is_filled)
        xs@filled <- as.integer(fld)
    }
    ## @dataCorrection (numeric) ? in xcmsSet function, if lockMassFreq.
    ## @progressInfo skip
    ## @progressCallback skip
    if (!any(colnames(pData(from)) == "class"))
        message("Note: you might want to set/adjust the",
                " 'sampclass' of the returned xcmSet object",
                " before proceeding with the analysis.")
    if (validObject(xs))
        return(xs)
}

.XCMSnExp2SummarizedExperiment <- function(x, ...) {
    if (!hasFeatures(x))
        stop("No correspondence analysis results present. Please run ",
             "groupChromPeaks first.")
    SummarizedExperiment(assays = list(raw = featureValues(x, ...)),
                         rowData = featureDefinitions(x),
                         colData = pData(x),
                         metadata = processHistory(x))
}

#' @description
#'
#' Extract a \code{data.frame} of retention time, mz and intensity
#' values from each file/sample in the provided rt-mz range.
#'
#' @note
#'
#' Ideally, \code{x} should be an \code{OnDiskMSnExp} object as subsetting
#' of a \code{XCMSnExp} object is more costly (removing of preprocessing
#' results, restoring data etc). If retention times reported in the
#' featureData are replaced by adjusted retention times, these are set
#' in the Spectrum objects as retention time.
#'
#' @param x An \code{OnDiskMSnExp} object.
#'
#' @param rt \code{numeric(2)} with the retention time range from which the
#'     data should be extracted.
#'
#' @param mz \code{numeric(2)} with the mz range.
#'
#' @param msLevel \code{integer} defining the MS level(s) to which the data
#'     should be restricted prior to data extraction.
#'
#' @return
#'
#' A \code{list} with length equal to the number of files and
#' each element being a \code{data.frame} with the extracted values.
#'
#' @noRd
#'
#' @author Johannes Rainer
.extractMsData <- function(x, rt, mz, msLevel = 1L) {
    if (!missing(rt)) {
        rt <- range(rt, na.rm = TRUE)
        if (length(rt) != 2)
            stop("'rt' has to be a numeric of length 2!")
    }
    if (!missing(mz)) {
        mz <- range(mz, na.rm = TRUE)
        if (length(mz) != 2)
            stop("'mz' has to be a numeric of length 2!")
        fmzr <- mz
    } else fmzr <- c(0, 0)
    ## Subset the object based on MS level rt and mz range.
    subs <- filterMz(filterRt(filterMsLevel(x, msLevel = msLevel), rt = rt),
                     mz = mz)
    if (length(subs) == 0) {
        ## Return a list with empty data.frames
        empty_df <- data.frame(rt = numeric(), mz = numeric(), i = integer())
        return(lapply(1:length(fileNames(x)), FUN = function(z){empty_df}))
    }
    suppressWarnings(
        dfs <- spectrapply(subs, FUN = function(z) {
            if (!z@peaksCount)
                return(data.frame(rt = numeric(), mz = numeric(),
                                  i = integer()))
            data.frame(rt = rep_len(z@rt, length(z@mz)),
                       mz = z@mz, i = z@intensity)
        })
    )
    fns <- fileNames(x)
    fromF <- base::match(fileNames(subs), fns)

    ## Now I want to rbind the spectrum data frames per file
    L <- split(dfs, f = fromFile(subs))
    L <- lapply(L, do.call, what = rbind)
    ## Put them into a vector same length that we have files.
    res <- vector(mode = "list", length = length(fns))
    res[fromF] <- L
    return(res)
}

#' @description
#'
#' Integrates the intensities for chromatograpic peak(s). This is
#' supposed to be called by the fillChromPeaks method.
#'
#' @note This reads the full data first and does the subsetting later in R.
#'
#' @param object An \code{XCMSnExp} object representing a single sample.
#'
#' @param peakArea A \code{matrix} with the peak definition, i.e.
#'     \code{"rtmin"}, \code{"rtmax"}, \code{"mzmin"} and \code{"mzmax"}.
#'
#' @param sample_idx \code{integer(1)} with the index of the sample in the
#'     object.
#'
#' @param mzCenterFun Name of the function to be used to calculate the mz value.
#'     Defaults to \code{weighted.mean}, i.e. the intensity weighted mean mz.
#'
#' @param cn \code{character} with the names of the result matrix.
#'
#' @return
#'
#' A \code{matrix} with at least columns \code{"mz"}, \code{"rt"},
#' \code{"into"} and \code{"maxo"} with the by intensity weighted mean of
#' mz, rt or the maximal intensity in the area, the integrated signal in
#' the area and the maximal signal in the area.
#'
#' @noRd
.getChromPeakData <- function(object, peakArea, sample_idx,
                              mzCenterFun = "weighted.mean",
                              cn = c("mz", "rt", "into", "maxo", "sample"),
                              msLevel = 1L) {
    if (length(fileNames(object)) != 1)
        stop("'object' should be an XCMSnExp for a single file!")
    ncols <- length(cn)
    res <- matrix(ncol = ncols, nrow = nrow(peakArea))
    colnames(res) <- cn
    res[, "sample"] <- sample_idx
    res[, c("mzmin", "mzmax")] <- peakArea[, c("mzmin", "mzmax")]
    ## Load the data
    message("Requesting ", nrow(res), " peaks from ",
            basename(fileNames(object)), " ... ", appendLF = FALSE)
    object <- filterRt(
        object, rt = range(peakArea[, c("rtmin", "rtmax")]) + c(-2, 2))
    object <- filterMsLevel(object, msLevel)
    if (!length(object)) {
        message("FAIL: no MS level ", msLevel, " data available.")
        return(res)
    }
    spctr <- spectra(object, BPPARAM = SerialParam())
    mzs <- lapply(spctr, function(z) z@mz)
    valsPerSpect <- lengths(mzs)
    ints <- unlist(lapply(spctr, function(z) z@intensity), use.names = FALSE)
    rm(spctr)
    mzs <- unlist(mzs, use.names = FALSE)
    mzs_range <- range(mzs)
    rtim <- rtime(object)
    rtim_range <- range(rtim)
    for (i in seq_len(nrow(res))) {
        rtr <- peakArea[i, c("rtmin", "rtmax")]
        mzr <- peakArea[i, c("mzmin", "mzmax")]
        ## If the rt range is completely out; additional fix for #267
        if (rtr[2] < rtim_range[1] | rtr[1] > rtim_range[2] |
            mzr[2] < mzs_range[1] | mzr[1] > mzs_range[2]) {
            res[i, ] <- rep(NA_real_, ncols)
            next
        }
        ## Ensure that the mz and rt region is within the range of the data.
        rtr[1] <- max(rtr[1], rtim_range[1])
        rtr[2] <- min(rtr[2], rtim_range[2])
        mzr[1] <- max(mzr[1], mzs_range[1])
        mzr[2] <- min(mzr[2], mzs_range[2])
        mtx <- .rawMat(mz = mzs, int = ints, scantime = rtim,
                       valsPerSpect = valsPerSpect, rtrange = rtr,
                       mzrange = mzr)
        if (length(mtx)) {
            if (any(!is.na(mtx[, 3]))) {
                ## How to calculate the area: (1)sum of all intensities / (2)by
                ## the number of data points (REAL ones, considering also NAs)
                ## and multiplied with the (3)rt width.
                ## (1) sum(mtx[, 3], na.rm = TRUE)
                ## (2) sum(rtim >= rtr[1] & rtim <= rtr[2]) - 1 ; if we used
                ## nrow(mtx) here, which would correspond to the non-NA
                ## intensities within the rt range we don't get the same results
                ## as e.g. centWave. Using max(1, ... to avoid getting Inf in
                ## case the signal is based on a single data point.
                ## (3) rtr[2] - rtr[1]
                res[i, "into"] <- sum(mtx[, 3], na.rm = TRUE) *
                    ((rtr[2] - rtr[1]) /
                     max(1, (sum(rtim >= rtr[1] & rtim <= rtr[2]) - 1)))
                maxi <- which.max(mtx[, 3])
                res[i, c("rt", "maxo")] <- mtx[maxi[1], c(1, 3)]
                res[i, c("rtmin", "rtmax")] <- rtr
                ## Calculate the intensity weighted mean mz
                meanMz <- do.call(mzCenterFun, list(mtx[, 2], mtx[, 3]))
                if (is.na(meanMz)) meanMz <- mtx[maxi[1], 2]
                res[i, "mz"] <- meanMz
            } else {
                res[i, ] <- rep(NA_real_, ncols)
            }
        } else {
            res[i, ] <- rep(NA_real_, ncols)
        }
    }
    message("got ", sum(!is.na(res[, "into"])), ".")
    res
}

#' @description Same as getChromPeakData, just without retention time.
#'
#' @note
#'
#' The mz and maxo are however estimated differently than for the
#' getChromPeakData: mz is the mz closest to the median mz of the feature and
#' maxo its intensity.
#'
#' @noRd
.getMSWPeakData <- function(object, peakArea, sample_idx,
                              cn = c("mz", "rt", "into", "maxo", "sample")) {
    if (length(fileNames(object)) != 1)
        stop("'object' should be an XCMSnExp for a single file!")
    ncols <- length(cn)
    res <- matrix(ncol = ncols, nrow = nrow(peakArea))
    colnames(res) <- cn
    res[, "sample"] <- sample_idx
    res[, "rt"] <- -1
    res[, "rtmin"] <- -1
    res[, "rtmax"] <- -1
    res[, c("mzmin", "mzmax")] <- peakArea[, c("mzmin", "mzmax")]
    ## Load the data
    message("Requesting ", nrow(res), " peaks from ",
             basename(fileNames(object)), " ... ", appendLF = FALSE)
    spctr <- spectra(object, BPPARAM = SerialParam())
    mzs <- lapply(spctr, mz)
    valsPerSpect <- lengths(mzs)
    ints <- unlist(lapply(spctr, intensity), use.names = FALSE)
    rm(spctr)
    mzs <- unlist(mzs, use.names = FALSE)
    for (i in 1:nrow(res)) {
        mz_area <- which(mzs >= peakArea[i, "mzmin"] &
                         mzs <= peakArea[i, "mzmax"])
        ## Alternative version from original code: but this can also pick up
        ## mzs from outside of the range! See also comments on issue #130
        ## mz_area <- seq(which.min(abs(mzs - peakArea[i, "mzmin"])),
        ##                which.min(abs(mzs - peakArea[i, "mzmax"])))
        mtx <- cbind(time = -1, mz = mzs[mz_area], intensity = ints[mz_area])
        ## mtx <- xcms:::.rawMat(mz = mzs, int = ints, scantime = rtime(object),
        ##                valsPerSpect = valsPerSpect,
        ##                mzrange = peakArea[i, c("mzmin", "mzmax")])
        if (length(mz_area)) {
            if (!all(is.na(mtx[, 3]))) {
                ## How to calculate the area: (1)sum of all intensities
                res[i, "into"] <- sum(mtx[, 3], na.rm = TRUE)
                ## Get the index of the mz value(s) closest to the mzmed of the
                ## feature
                mzDiff <- abs(mtx[, 2] - peakArea[i, "mzmed"])
                mz_idx <- which(mzDiff == min(mzDiff))
                ## Now get the one with the highest intensity.
                maxi <- mz_idx[which.max(mtx[mz_idx, 3])]
                ## Return these.
                res[i, c("mz", "maxo")] <- mtx[maxi, 2:3]
                ## ## mz should be the weighted mean!
                ## res[i, c("mz", "maxo")] <- c(weighted.mean(mtx[, 2], mtx[, 3]),
                ##                              mtx[maxi[1], 3])
            } else {
                res[i, ] <- rep(NA_real_, ncols)
            }
        } else {
            res[i, ] <- rep(NA_real_, ncols)
        }
    }
    message("got ", sum(!is.na(res[, "into"])), ".")
    return(res)
}
## The same version as above, but the maxo is the maximum signal of the peak,
## and the mz the intensity weighted mean mz.
.getMSWPeakData2 <- function(object, peakArea, sample_idx,
                              cn = c("mz", "rt", "into", "maxo", "sample")) {
    if (length(fileNames(object)) != 1)
        stop("'object' should be an XCMSnExp for a single file!")
    ncols <- length(cn)
    res <- matrix(ncol = ncols, nrow = nrow(peakArea))
    colnames(res) <- cn
    res[, "sample"] <- sample_idx
    res[, "rt"] <- -1
    res[, "rtmin"] <- -1
    res[, "rtmax"] <- -1
    res[, c("mzmin", "mzmax")] <- peakArea[, c("mzmin", "mzmax")]
    ## Load the data
    message("Reguesting ", nrow(res), " peaks from ",
             basename(fileNames(object)), " ... ", appendLF = FALSE)
    spctr <- spectra(object, BPPARAM = SerialParam())
    mzs <- lapply(spctr, mz)
    valsPerSpect <- lengths(mzs)
    ints <- unlist(lapply(spctr, intensity), use.names = FALSE)
    rm(spctr)
    mzs <- unlist(mzs, use.names = FALSE)
    for (i in 1:nrow(res)) {
        mtx <- .rawMat(mz = mzs, int = ints, scantime = rtime(object),
                       valsPerSpect = valsPerSpect,
                       mzrange = peakArea[i, c("mzmin", "mzmax")])
        if (length(mtx)) {
            if (!all(is.na(mtx[, 3]))) {
                ## How to calculate the area: (1)sum of all intensities
                res[i, "into"] <- sum(mtx[, 3], na.rm = TRUE)
                res[i, c("mz", "maxo")] <- c(weighted.mean(mtx[, 2], mtx[, 3]),
                                             max(mtx[, 3], na.rm = TRUE))
            } else {
                res[i, ] <- rep(NA_real_, ncols)
            }
        } else {
            res[i, ] <- rep(NA_real_, ncols)
        }
    }
    message("got ", sum(!is.na(res[, "into"])), ".")
    return(res)
}

## Same as .getChromPeakData but for matchedFilter, i.e. using the profile
## matrix instead of the original signal.
.getChromPeakData_matchedFilter <- function(object, peakArea, sample_idx,
                                            mzCenterFun = "weighted.mean",
                                            param = MatchedFilterParam(),
                                            cn = c("mz", "rt", "into", "maxo",
                                                   "sample"), msLevel = 1L) {
    if (length(fileNames(object)) != 1)
        stop("'object' should be an XCMSnExp for a single file!")
    ncols <- length(cn)
    res <- matrix(ncol = ncols, nrow = nrow(peakArea))
    colnames(res) <- cn
    res[, "sample"] <- sample_idx
    res[, c("mzmin", "mzmax")] <-
        peakArea[, c("mzmin", "mzmax")]
    ## Load the data
    message("Requesting ", nrow(res), " peaks from ",
            basename(fileNames(object)), " ... ", appendLF = FALSE)
    object <- filterRt(
        object, rt = range(peakArea[, c("rtmin", "rtmax")]) + c(-2, 2))
    object <- filterMsLevel(object, msLevel)
    if (!length(object)) {
        message("FAIL: no MS level ", msLevel, " data available.")
        return(res)
    }
    spctr <- spectra(object, BPPARAM = SerialParam())
    mzs <- lapply(spctr, mz)
    vps <- lengths(mzs)
    ints <- unlist(lapply(spctr, intensity), use.names = FALSE)
    rm(spctr)
    mzs <- unlist(mzs, use.names = FALSE)
    rtim <- rtime(object)
    rtim_range <- range(rtim)
    ## Now, if we do have "distance" defined it get's tricky:
    basespc <- NULL
    if (length(distance(param)) > 0) {
        mass <- seq(floor(min(mzs) / binSize(param)) * binSize(param),
                    ceiling(max(mzs) / binSize(param)) * binSize(param),
                    by = binSize(param))
        bin_size <- (max(mass) - min(mass)) / (length(mass) - 1)
        basespc <- distance(param) * bin_size
    }
    ## Create the profile matrix:
    pMat <- .createProfileMatrix(mz = mzs, int = ints, valsPerSpect = vps,
                                 method = .impute2method(param),
                                 step = binSize(param),
                                 baselevel = baseValue(param),
                                 basespace = basespc,
                                 returnBreaks = TRUE,
                                 baseValue = NA)  # We want to return NA not 0
                                        # if nothing was found
    brks <- pMat$breaks
    pMat <- pMat$profMat  ## rows are masses, cols are retention times/scans.
    bin_size <- diff(brks[1:2])
    bin_half <- bin_size / 2
    ## Calculate the mean mass per bin using the breaks used for the binning.
    mass <- brks[-length(brks)] + bin_half ## midpoint for the breaks
    mass_range <- range(mass)

    for (i in 1:nrow(res)) {
        rtr <- peakArea[i, c("rtmin", "rtmax")]
        mzr <- peakArea[i, c("mzmin", "mzmax")]
        ## Ensure that the rt region is within the rtrange of the data.
        rtr[1] <- max(rtr[1], rtim_range[1])
        rtr[2] <- min(rtr[2], rtim_range[2])
        mzr[1] <- max(mzr[1], mass_range[1])
        mzr[2] <- min(mzr[2], mass_range[2])
        ## Get the index of rt in rtim that are within the rt range rtr
        ## range_rt <- c(min(which(rtim >= rtr[1])), max(which(rtim <= rtr[2])))
        ## range_rt <- c(which.min(abs(rtim - rtr[1])),
        ##               which.min(abs(rtim - rtr[2])))
        range_rt <- findRange(rtim, rtr, TRUE)
        idx_rt <- range_rt[1]:range_rt[2]
        ## Get the index of the mz in the data that are within the mz range.
        ##range_mz <- c(min(which(brks >= mzr[1])) - 1, max(which(brks <= mzr[2])))
        range_mz <- findRange(mass, c(mzr[1] - bin_half, mzr[2] + bin_half),
                              TRUE)
        idx_mz <- range_mz[1]:range_mz[2]

        if (length(idx_mz) > 0 & length(idx_rt) > 0) {
            intMat <- pMat[idx_mz, idx_rt, drop = FALSE]
            is_na <- is.na(intMat)
            if (all(is_na)) {
                res[i, ] <- rep(NA_real_, ncols)
                next
            }
            intMat_0 <- intMat
            intMat_0[is.na(intMat)] <- 0
            ## Calculate the mean mz value using a intensity weighted mean.
            mz_ints <- rowSums(intMat, na.rm = TRUE)
            ## mz_ints <- Biobase::rowMax(intMat)
            if (length(mz_ints) != length(idx_mz)) {
                ## Take from original code...
                warning("weighted.mean: x and weights have to have same length!")
                mz_ints <- rep(1, length(idx_mz))
            }
            ## mean mz: intensity weighted mean
            mean_mz <- weighted.mean(mass[idx_mz], mz_ints, na.rm = TRUE)
            if (is.nan(mean_mz) || is.na(mean_mz))
                mean_mz <- mean(mzr, na.rm = TRUE)
            res[i, "mz"] <- mean_mz
            ## mean rt: position of the maximum intensity (along rt.)
            rt_ints <- colMax(intMat_0, na.rm = TRUE)
            res[i, c("rt", "rtmin", "rtmax")] <-
                c(rtim[idx_rt][which.max(rt_ints)], rtr)
            ## maxo
            res[i, "maxo"] <- max(rt_ints, na.rm = TRUE)
            ## into
            rt_width <- diff(rtim[range_rt])/diff(range_rt)
            res[i, "into"] <- rt_width * sum(rt_ints, na.rm = TRUE)
        } else
            res[i, ] <- rep(NA_real_, ncols)
    }
    message("got ", sum(!is.na(res[, "into"])), ".")
    return(res)
}

.hasFilledPeaks <- function(object) {
    hasChromPeaks(object) & any(chromPeakData(object)$is_filled, na.rm = TRUE)
}

#' @description
#'
#' Simple helper function to extract the peakidx column from the
#' featureDefinitions DataFrame. The function ensures that the names of the
#' returned list correspond to the rownames of the DataFrame
#'
#' @noRd
.peakIndex <- function(object) {
    if (inherits(object, "DataFrame")) {
        idxs <- object$peakidx
        names(idxs) <- rownames(object)
    } else {
        if (!hasFeatures(object))
            stop("No feature definitions present. Please run groupChromPeaks first.")
        idxs <- featureDefinitions(object)$peakidx
        names(idxs) <- rownames(featureDefinitions(object))
    }
    idxs
}

#' @description
#'
#' \code{adjustRtimePeakGroups} returns the features (peak groups)
#' which would, depending on the provided \code{\link{PeakGroupsParam}}, be
#' selected for alignment/retention time correction.
#'
#' @note
#'
#' \code{adjustRtimePeakGroups} is supposed to be called \emph{before} the
#' sample alignment, but after a correspondence (peak grouping).
#'
#' @return
#'
#' For \code{adjustRtimePeakGroups}: a \code{matrix}, rows being
#' features, columns samples, of retention times. The features are ordered
#' by the median retention time across columns.
#'
#' @rdname adjustRtime-peakGroups
adjustRtimePeakGroups <- function(object, param = PeakGroupsParam(),
                                  msLevel = 1L) {
    if (!is(object, "XCMSnExp"))
        stop("'object' has to be an 'XCMSnExp' object.")
    if (!hasFeatures(object))
        stop("No features present. Please run 'groupChromPeaks' first.")
    if (hasAdjustedRtime(object))
        warning("Alignment/retention time correction was already performed, ",
                "returning a matrix with adjusted retention times.")
    subs <- subset(param)
    if (!length(subs))
        subs <- seq_along(fileNames(object))
    nSamples <- length(subs)
    pkGrp <- .getPeakGroupsRtMatrix(
        peaks = chromPeaks(object, msLevel = msLevel),
        peakIndex = .peakIndex(
            .update_feature_definitions(
                featureDefinitions(object), rownames(chromPeaks(object)),
                rownames(chromPeaks(object, msLevel = msLevel)))),
        sampleIndex = subs,
        missingSample = nSamples - (nSamples * minFraction(param)),
        extraPeaks = extraPeaks(param)
    )
    colnames(pkGrp) <- basename(fileNames(object))[subs]
    pkGrp
}

#' @title Visualization of alignment results
#'
#' @description
#'
#' Plot the difference between the adjusted and the raw retention
#' time (y-axis) for each file along the (adjusted or raw) retention time
#' (x-axis). If alignment was performed using the
#' \code{\link{adjustRtime-peakGroups}} method, also the features (peak
#' groups) used for the alignment are shown.
#'
#' @param object A \code{\link{XCMSnExp}} object with the alignment results.
#'
#' @param col colors to be used for the lines corresponding to the individual
#'     samples.
#'
#' @param lwd line width to be used for the lines of the individual samples.
#'
#' @param lty line type to be used for the lines of the individual samples.
#'
#' @param type plot type to be used. See help on the \code{par} function for
#'     supported values.
#'
#' @param adjustedRtime logical(1) whether adjusted or raw retention times
#'     should be shown on the x-axis.
#'
#' @param xlab the label for the x-axis.
#'
#' @param ylab the label for the y-axis.
#'
#' @param peakGroupsCol color to be used for the peak groups (only used if
#'     alignment was performed using the \code{\link{adjustRtime-peakGroups}}
#'     method.
#'
#' @param peakGroupsPch point character (\code{pch}) to be used for the peak
#'     groups (only used if alignment was performed using the
#'     \code{\link{adjustRtime-peakGroups}} method.
#'
#' @param peakGroupsLty line type (\code{lty}) to be used to connect points for
#'     each peak groups (only used if alignment was performed using the
#'     \code{\link{adjustRtime-peakGroups}} method.
#'
#' @param ylim optional \code{numeric(2)} with the upper and lower limits on
#'     the y-axis.
#'
#' @param ... Additional arguments to be passed down to the \code{plot}
#'     function.
#'
#' @seealso \code{\link{adjustRtime}} for all retention time correction/
#'     alignment methods.
#'
#' @author Johannes Rainer
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
#' ## Performing the peak grouping using the "peak density" method.
#' p <- PeakDensityParam(sampleGroups = c(1, 1, 1))
#' res <- groupChromPeaks(faahko_sub, param = p)
#'
#' ## Perform the retention time adjustment using peak groups found in both
#' ## files.
#' fgp <- PeakGroupsParam(minFraction = 1)
#' res <- adjustRtime(res, param = fgp)
#'
#' ## Visualize the impact of the alignment.
#' plotAdjustedRtime(res, adjusted = FALSE)
#' grid()
plotAdjustedRtime <- function(object, col = "#00000080", lty = 1, lwd = 1,
                              type = "l", adjustedRtime = TRUE,
                              xlab = ifelse(adjustedRtime,
                                            yes = expression(rt[adj]),
                                            no = expression(rt[raw])),
                              ylab = expression(rt[adj]-rt[raw]),
                              peakGroupsCol = "#00000060",
                              peakGroupsPch = 16,
                              peakGroupsLty = 3, ylim, ...) {
    if (!is(object, "XCMSnExp"))
        stop("'object' has to be an 'XCMSnExp' object.")
    if (!hasAdjustedRtime(object))
        warning("No alignment/retention time correction results present.")
    diffRt <- rtime(object, adjusted = TRUE) - rtime(object, adjusted = FALSE)
    diffRt <- split(diffRt, fromFile(object))
    ## Define the rt that is shown on x-axis
    xRt <- rtime(object, adjusted = adjustedRtime, bySample = TRUE)
    ## Check colors.
    if (length(col) == 1)
        col <- rep(col, length(diffRt))
    if (length(lty) == 1)
        lty <- rep(lty, length(diffRt))
    if (length(lwd) == 1)
        lwd <- rep(lwd, length(diffRt))
    if (length(col) != length(diffRt)) {
        warning("length of 'col' does not match the number of samples! Will ",
                "use 'col[1]' for all samples.")
        col <- rep(col[1], length(diffRt))
    }
    if (length(lty) != length(lty)) {
        warning("length of 'lty' does not match the number of samples! Will ",
                "use 'lty[1]' for all samples.")
        lty <- rep(lty[1], length(diffRt))
    }
    if (length(lwd) != length(lwd)) {
        warning("length of 'lwd' does not match the number of samples! Will ",
                "use 'lwd[1]' for all samples.")
        lwd <- rep(lwd[1], length(diffRt))
    }
    ## Initialize plot.
    if (missing(ylim))
        ylim <- range(diffRt, na.rm = TRUE)
    plot(3, 3, pch = NA, xlim = range(xRt, na.rm = TRUE),
         ylim = ylim, xlab = xlab, ylab = ylab, ...)
    ## Plot all.
    for (i in 1:length(diffRt))
        points(x = xRt[[i]], y = diffRt[[i]], col = col[i], lty = lty[i],
               type = type, lwd = lwd[i])
    ## If alignment was performed using the peak groups method highlight also
    ## those in the plot.
    ph <- processHistory(object, type = .PROCSTEP.RTIME.CORRECTION)
    if (length(ph)) {
        ph <- ph[[length(ph)]]
        if (is(ph, "XProcessHistory")) {
            ## Check if we've got a PeakGroupsParam parameter class
            prm <- processParam(ph)
            if (is(prm, "PeakGroupsParam")) {
                rm(diffRt)
                rm(xRt)
                rawRt <- rtime(object, adjusted = FALSE, bySample = TRUE)
                adjRt <- rtime(object, adjusted = TRUE, bySample = TRUE)
                pkGroup <- peakGroupsMatrix(prm)
                subs <- subset(prm)
                if (!length(subs))
                    subs <- 1:ncol(pkGroup)
                rawRt <- rawRt[subs]
                adjRt <- adjRt[subs]
                ## Have to "adjust" these:
                pkGroupAdj <- pkGroup
                for (i in 1:ncol(pkGroup)) {
                    pkGroupAdj[, i] <- .applyRtAdjustment(pkGroup[, i],
                                                                 rawRt[[i]],
                                                                 adjRt[[i]])
                }
                diffRt <- pkGroupAdj - pkGroup
                if (adjustedRtime)
                    xRt <- pkGroupAdj
                else
                    xRt <- pkGroup
                ## Loop through the rows and plot points - ordered by diffRt!
                for (i in 1:nrow(xRt)) {
                    idx <- order(diffRt[i, ])
                    points(x = xRt[i, ][idx], diffRt[i, ][idx],
                           col = peakGroupsCol, type = "b",
                           pch = peakGroupsPch, lty = peakGroupsLty)
                }
            }
        }
    }
}

.plotChromPeakDensity <- function(object, mz, rt, param, simulate = TRUE,
                                 col = "#00000080", xlab = "retention time",
                                 ylab = "sample", xlim = range(rt),
                                 main = NULL, type = c("any", "within",
                                                       "apex_within"), ...) {
    .Deprecated(
        msg = paste0("Use of 'plotChromPeakDensity' on 'XCMSnExp' is",
                     "discouraged. Please extract chromatographic ",
                     "data first and call 'plotChromPeakDensity' ",
                     "directly on the 'XChromatograms' object. See ",
                     "?XChromatograms, section 'Correspondence ",
                     "analysis' for more details."))
    type <- match.arg(type)
    if (missing(object))
        stop("Required parameter 'object' is missing")
    if (!is(object, "XCMSnExp"))
        stop("'object' must be an XCMSnExp object")
    if (!hasChromPeaks(object))
        stop("No chromatographic peaks present in 'object'")
    if (missing(mz))
        mz <- c(-Inf, Inf)
    if (missing(rt))
        rt <- range(rtime(object))
    if (!simulate & !hasFeatures(object)) {
        warning("Falling back to 'simulate = TRUE' because no correspondence ",
                "results are present")
        simulate <- TRUE
    }
    if (missing(param)) {
        ## Use the internal parameter - if available.
        if (hasFeatures(object)) {
            ph <- processHistory(object, type = .PROCSTEP.PEAK.GROUPING)
            param <- processParam(ph[[length(ph)]])
        } else {
            param = PeakDensityParam(
                sampleGroups = rep(1, length(fileNames(object))))
        }
    }
    if (!is(param, "PeakDensityParam"))
        stop("'param' has to be a 'PeakDensityParam'")
    mz <- range(mz)
    rt <- range(rt)
    ## Get all the data we require.
    nsamples <- length(fileNames(object))
    if (length(col) != nsamples)
        col <- rep_len(col[1], nsamples)
    pks <- chromPeaks(object, mz = mz, rt = rt, type = type, msLevel = 1L)
    if (nrow(pks)) {
        ## Extract parameters from the param object
        bw = bw(param)
        ## Ensure the density is calculated exactly as it would in the real
        ## case.
        full_rt_range <- range(chromPeaks(object)[, "rt"])
        dens_from <- full_rt_range[1] - 3 * bw
        dens_to <- full_rt_range[2] + 3 * bw
        ## That's Jan Stanstrup's fix (issue #161).
        densN <- max(512, 2 * 2^(ceiling(log2(diff(full_rt_range) / (bw / 2)))))
        sample_groups <- sampleGroups(param)
        if (length(sample_groups) == 0)
            sample_groups <- rep(1, nsamples)
        if (length(sample_groups) != nsamples)
            stop("If provided, the 'sampleGroups' parameter in the 'param' ",
                 "class has to have the same length than there are samples ",
                 "in 'object'")
        sample_groups_table <- table(sample_groups)
        dens <- density(pks[, "rt"], bw = bw, from = dens_from, to = dens_to,
                        n = densN)
        yl <- c(0, max(dens$y))
        ypos <- seq(from = yl[1], to = yl[2], length.out = nsamples)
        if (is.null(main))
            main <- paste0(format(mz, digits = 7), collapse = " - ")
        ## Plot the peaks as points.
        ## Fix assignment of point types for each sample.
        dots <- list(...)
        if (any(names(dots) == "bg")) {
            bg <- dots$bg
            if (length(bg) != nsamples)
                bg <- rep_len(bg[1], nsamples)
            dots$bg <- bg[pks[, "sample"]]
        }
        if (any(names(dots) == "pch")) {
            pch <- dots$pch
            if (length(pch) != nsamples)
                pch <- rep_len(pch[1], nsamples)
            dots$pch <- pch[pks[, "sample"]]
        }
        do.call("plot", args = c(list(x = pks[, "rt"],
                                      y = ypos[pks[, "sample"]], xlim = xlim,
                                      col = col[pks[, "sample"]], xlab = xlab,
                                      yaxt = "n", ylab = ylab,
                                      main = main, ylim = yl), dots))
        axis(side = 2, at = ypos, labels = 1:nsamples)
        points(x = dens$x, y = dens$y, type = "l")
        ## Estimate what would be combined to a feature
        ## Code is taken from do_groupChromPeaks_density
        dens_max <- max(dens$y)
        dens_y <- dens$y
        if (simulate) {
            snum <- 0
            while(dens_y[max_y <- which.max(dens_y)] > dens_max / 20 &&
                  snum < maxFeatures(param)) {
                      feat_range <- descendMin(dens_y, max_y)
                      dens_y[feat_range[1]:feat_range[2]] <- 0
                      feat_idx <- which(pks[, "rt"] >= dens$x[feat_range[1]] &
                                        pks[, "rt"] <= dens$x[feat_range[2]])
                      tt <- table(sample_groups[pks[feat_idx, "sample"]])
                      if (!any(tt / sample_groups_table[names(tt)] >=
                               minFraction(param) & tt >= minSamples(param)))
                          next
                      rect(xleft = min(pks[feat_idx, "rt"]), ybottom = yl[1],
                           xright = max(pks[feat_idx, "rt"]), ytop = yl[2],
                           border = "#00000040", col = "#00000020")
                  }
        } else {
            ## Plot all features in the region.
            fts <- featureDefinitions(object, mz = mz, rt = rt, type = type)
            if (nrow(fts)) {
                rect(xleft = fts$rtmin, xright = fts$rtmax,
                     ybottom = rep(yl[1], nrow(fts)),
                     ytop = rep(yl[2], nrow(fts)),
                     border = "#00000040", col = "#00000020")
                abline(v = fts$rtmed, col = "#00000040", lty = 2)
            }
        }
    } else {
        plot(3, 3, pch = NA, xlim = rt, xlab = xlab,
             main = paste0(format(mz, digits = 7), collapse = " - "))
    }
}

#' @title Add definition of chromatographic peaks to an extracted chromatogram
#'     plot
#'
#' @description
#'
#' The \code{highlightChromPeaks} function adds chromatographic
#' peak definitions to an existing plot, such as one created by the
#' \code{plot} method on a \code{\link{Chromatogram}} or
#' \code{\link{MChromatograms}} object.
#'
#' @param x For \code{highlightChromPeaks}: \code{XCMSnExp} object with the
#'     detected peaks.
#'
#' @param rt For \code{highlightChromPeaks}: \code{numeric(2)} with the
#'     retention time range from which peaks should be extracted and plotted.
#'
#' @param mz \code{numeric(2)} with the mz range from which the peaks should
#'     be extracted and plotted.
#'
#' @param peakIds \code{character} defining the IDs (i.e. rownames of the peak
#'     in the \code{chromPeaks} table) of the chromatographic peaks to be
#'     highlighted in a plot.
#'
#' @param border colors to be used to color the border of the rectangles/peaks.
#'     Has to be equal to the number of samples in \code{x}.
#'
#' @param lwd \code{numeric(1)} defining the width of the line/border.
#'
#' @param col For \code{highlightChromPeaks}: color to be used to fill the
#'     rectangle (if \code{type = "rect"}) or the peak
#'     (for \code{type = "polygon"}).
#'
#' @param type the plotting type. See \code{\link{plot}} in base grapics for
#'     more details.
#'     For \code{highlightChromPeaks}: \code{character(1)} defining how the peak
#'     should be highlighted: \code{type = "rect"} draws a rectangle
#'     representing the peak definition, \code{type = "point"} indicates a
#'     chromatographic peak with a single point at the position of the peak's
#'     \code{"rt"} and \code{"maxo"} and \code{type = "polygon"} will highlight
#'     the peak shape. For \code{type = "polygon"} the color of the border and
#'     area can be defined with parameters \code{"border"} and \code{"col"},
#'     respectively.
#'
#' @param whichPeaks \code{character(1)} specifying how peaks are called to be
#'     located within the region defined by \code{mz} and \code{rt}. Can be
#'     one of \code{"any"}, \code{"within"}, and \code{"apex_within"} for all
#'     peaks that are even partially overlapping the region, peaks that are
#'     completely within the region, and peaks for which the apex is within
#'     the region. This parameter is passed to the \code{type} argument of the
#'     \code{\link{chromPeaks}} function. See related documentation for more
#'     information and examples.
#'
#' @param ... additional parameters to the \code{\link{matplot}} or \code{plot}
#'     function.
#'
#' @author Johannes Rainer
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
#' ## Extract the ion chromatogram for one chromatographic peak in the data.
#' chrs <- chromatogram(faahko_sub, rt = c(2700, 2900), mz = 335)
#'
#' plot(chrs)
#'
#' ## Extract chromatographic peaks for the mz/rt range (if any).
#' chromPeaks(faahko_sub, rt = c(2700, 2900), mz = 335)
#'
#' ## Highlight the chromatographic peaks in the area
#' ## Show the peak definition with a rectangle
#' highlightChromPeaks(faahko_sub, rt = c(2700, 2900), mz = 335)
#'
#' ## Color the actual peak
#' highlightChromPeaks(faahko_sub, rt = c(2700, 2900), mz = 335,
#'     col = c("#ff000020", "#00ff0020"), type = "polygon")
highlightChromPeaks <- function(x, rt, mz, peakIds = character(),
                                border = rep("00000040", length(fileNames(x))),
                                lwd = 1, col = NA,
                                type = c("rect", "point", "polygon"),
                                whichPeaks = c("any", "within", "apex_within"),
                                ...) {
    type <- match.arg(type)
    msLevel <- 1L
    whichPeaks <- match.arg(whichPeaks)
    n_samples <- length(fileNames(x))
    if (!hasChromPeaks(x))
        stop("'x' does not contain detected peaks")
    if (length(peakIds)) {
        if (!missing(rt) | !missing(mz))
            warning("Ignoring 'rt' and 'mz' because peakIds were provided")
        if (!all(peakIds %in% rownames(chromPeaks(x, msLevel = msLevel))))
            stop("'peakIds' do not match rownames of 'chromPeaks(x)'")
        rt <- range(chromPeaks(x)[peakIds, c("rtmin", "rtmax")])
        mz <- range(chromPeaks(x)[peakIds, c("mzmin", "mzmax")])
    }
    if (missing(rt))
        rt <- c(-Inf, Inf)
    if (missing(mz))
        mz <- c(-Inf, Inf)
    if (!is(x, "XCMSnExp"))
        stop("'x' has to be a XCMSnExp object")
    if (length(peakIds))
        pks <- chromPeaks(x)[peakIds, , drop = FALSE]
    else pks <- chromPeaks(x, rt = rt, mz = mz, ppm = 0, type = whichPeaks,
                           msLevel = msLevel)
    if (length(col) != n_samples)
        col <- rep(col[1], n_samples)
    if (length(border) != n_samples)
        border <- rep(border[1], n_samples)
    if (length(pks)) {
        switch(type,
               rect = rect(xleft = pks[, "rtmin"], xright = pks[, "rtmax"],
                           ybottom = rep(0, nrow(pks)), ytop = pks[, "maxo"],
                           border = border[pks[, "sample"]], lwd = lwd,
                           col = col[pks[, "sample"]]),
               point = {
                   ## Fix assignment of point types for each sample.
                   dots <- list(...)
                   if (any(names(dots) == "bg")) {
                       bg <- dots$bg
                       if (length(bg) != n_samples)
                           bg <- rep_len(bg[1], n_samples)
                       dots$bg <- bg[pks[, "sample"]]
                   }
                   if (any(names(dots) == "pch")) {
                       pch <- dots$pch
                       if (length(pch) != n_samples)
                           pch <- rep_len(pch[1], n_samples)
                       dots$pch <- pch[pks[, "sample"]]
                   }
                   if (any(is.na(col)))
                       col <- border
                       ## Draw a point at the position defined by the "rt" column
                       do.call("points",
                               args = c(list(x = pks[, "rt"],
                                             y = pks[, "maxo"],
                                             col = col[pks[, "sample"]]), dots))
               },
               polygon = {
                   if (nrow(pks)) {
                       chrs <- chromatogram(
                           x, rt = range(pks[, c("rtmin", "rtmax")]), mz = mz)
                       pks <- pks[order(pks[, "maxo"], decreasing = TRUE), ,
                                  drop = FALSE]
                       for (j in seq_len(nrow(pks))) {
                           i <- pks[j, "sample"]
                           chr <- filterRt(chrs[1, i],
                                           rt = pks[j, c("rtmin", "rtmax")])
                           xs <- rtime(chr)
                           xs <- c(xs, xs[length(xs)], xs[1])
                           ys <- c(intensity(chr), 0, 0)
                           nona <- !is.na(ys)
                           polygon(xs[nona], ys[nona], border = border[i],
                                   col = col[i])
                       }
                   }
               })
    }
}

#' @title General visualizations of peak detection results
#'
#' @description
#'
#' \code{plotChromPeaks} plots the identified chromatographic
#' peaks from one file into the plane spanned by the retention time and mz
#' dimension (x-axis representing the retention time and y-axis mz).
#' Each chromatographic peak is plotted as a rectangle representing its
#' width in rt and mz dimension.
#'
#' This plot is supposed to provide some initial overview of the
#' chromatographic peak detection results.
#'
#' @details
#'
#' The width and line type of the rectangles indicating the detected
#' chromatographic peaks for the \code{plotChromPeaks} function can be
#' specified using the \code{par} function, i.e. with \code{par(lwd = 3)}
#' and \code{par(lty = 2)}, respectively.
#'
#' @param x \code{\link{XCMSnExp}} object.
#'
#' @param file For \code{plotChromPeaks}: \code{numeric(1)} specifying the
#'     index of the file within \code{x} for which the plot should be created.
#'     Defaults to \code{1}.
#'
#' @param xlim \code{numeric(2)} specifying the x-axis limits (retention time
#'     dimension). Defaults to \code{NULL} in which case the full retention
#'     time range of the file is used.
#'
#' @param ylim For \code{plotChromPeaks}: \code{numeric(2)} specifying the
#'     y-axis limits (mz dimension). Defaults to \code{NULL} in which case the
#'     full mz range of the file is used.
#'
#' @param add For \code{plotChromPeaks}: \code{logical(1)} whether the plot
#'     should be added or created as a new plot.
#'
#' @param border For \code{plotChromPeaks}: the color for the rectangles'
#'     border.
#'
#' @param col For \code{plotChromPeaks}: the color to be used to fill the
#'     rectangles.
#'
#' @param xlab \code{character(1)} defining the x-axis label.
#'
#' @param ylab For \code{plotChromPeaks}: \code{character(1)} defining the
#'     y-axis label.
#'
#' @param main \code{character(1)} defining the plot title. By default (i.e.
#'     \code{main = NULL} the name of the file will be used as title.
#'
#' @param ... Additional arguments passed to the \code{plot} (for
#'     \code{plotChromPeaks}) and \code{image} (for
#'     \code{plotChromPeakImage}) functions. Ignored if \code{add = TRUE}.
#'
#' @author Johannes Rainer
#'
#' @seealso \code{\link{highlightChromPeaks}} for the function to highlight
#'     detected chromatographic peaks in extracted ion chromatogram plots.
#'
#' @examples
#'
#' ## Load a test data set with detected peaks
#' data(faahko_sub)
#' ## Update the path to the files for the local system
#' dirname(faahko_sub) <- system.file("cdf/KO", package = "faahKO")
#'
#' ## plotChromPeakImage: plot an image for the identified peaks per file
#' plotChromPeakImage(faahko_sub)
#'
#' ## Show all detected chromatographic peaks from the first file
#' plotChromPeaks(faahko_sub)
#'
#' ## Plot all detected peaks from the second file and restrict the plot to a
#' ## mz-rt slice
#' plotChromPeaks(faahko_sub, file = 2, xlim = c(3500, 3600), ylim = c(400, 600))
plotChromPeaks <- function(x, file = 1, xlim = NULL, ylim = NULL,
                               add = FALSE, border = "#00000060", col = NA,
                               xlab = "retention time", ylab = "mz",
                               main = NULL, ...) {
    if (!is(x, "XCMSnExp"))
        stop("'x' is supposed to be an 'XCMSnExp' object, but I got a ",
             class(x))
    suppressMessages(
        x_file <- filterFile(x, file = file[1], keepAdjustedRtime = TRUE)
    )
    if (is.null(xlim))
        xlim <- range(rtime(x_file))
    if (is.null(ylim))
        ylim <- range(mz(x_file))
    if (is.null(main))
        main <- basename(fileNames(x_file))
    ## Get the peaks from the file, restricting to the current limits (might
    ## speed up things).
    pks <- chromPeaks(x_file, mz = ylim, rt = xlim, msLevel = 1L)
    ## Initialize plot
    if (!add)
        plot(3, 3, pch = NA, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab,
             main = main, ...)
    if (nrow(pks))
        rect(xleft = pks[, "rtmin"], xright = pks[, "rtmax"],
             ybottom = pks[, "mzmin"], ytop = pks[, "mzmax"], col = col,
             border = border)
}

#' @description
#'
#' \code{plotChromPeakImage} plots the number of detected peaks for
#' each sample along the retention time axis as an \emph{image} plot, i.e.
#' with the number of peaks detected in each bin along the retention time
#' represented with the color of the respective cell.
#'
#' @param binSize For \code{plotChromPeakImage}: \code{numeric(1)} defining the
#'     size of the bins along the x-axis (retention time). Defaults to
#'     \code{binSize = 30}, peaks within each 30 seconds will thus counted and
#'     plotted.
#'
#' @param log For \code{plotChromPeakImage}: \code{logical(1)} whether the peak
#'     counts should be log2 transformed before plotting.
#'
#' @param yaxt For \code{plotChromPeakImage}: \code{character(1)} defining
#'     whether y-axis labels should be added. To disable the y-axis use
#'     \code{yaxt = "n"}. For any other value of \code{yaxt} the axis will be
#'     drawn. See \code{par} help page for more details.
#'
#' @rdname plotChromPeaks
plotChromPeakImage <- function(x, binSize = 30, xlim = NULL, log = FALSE,
                               xlab = "retention time", yaxt = par("yaxt"),
                               main = "Chromatographic peak counts", ...) {
    if (!is(x, "XCMSnExp"))
        stop("'x' is supposed to be an 'XCMSnExp' object, but I got a ",
             class(x))
    if (is.null(xlim))
        xlim <- c(floor(min(rtime(x))), ceiling(max(rtime(x))))
    brks <- seq(xlim[1], xlim[2], by = binSize)
    if (brks[length(brks)] < xlim[2])
        brks <- c(brks, brks[length(brks)] + binSize)
    pks <- chromPeaks(x, rt = xlim, msLevel = 1L)
    if (nrow(pks)) {
        rts <- split(pks[, "rt"], as.factor(as.integer(pks[, "sample"])))
        cnts <- lapply(rts, function(z) {
            hst <- hist(z, breaks = brks, plot = FALSE)
            hst$counts
        })
        ## Add 0 vectors for samples in which no peaks were found.
        n_samples <- length(fileNames(x))
        sample_idxs <- 1:n_samples
        sample_idxs <- sample_idxs[!(as.character(sample_idxs) %in% names(rts))]
        if (length(sample_idxs)) {
            all_cnts <- vector("list", n_samples)
            all_cnts[as.numeric(names(cnts))] <- cnts
            zeros <- rep(0, (length(brks) - 1))
            all_cnts[sample_idxs] <- list(zeros)
            cnts <- all_cnts
        }
        cnts <- t(do.call(rbind, cnts))
        if (log)
            cnts <- log2(cnts)
        image(z = cnts, x = brks - (brks[2] - brks[1]) / 2, xaxs = "r",
              xlab = xlab, yaxt = "n", ...)
        if (yaxt != "n")
            axis(side = 2, at = seq(0, 1, length.out = n_samples),
                 labels = basename(fileNames(x)), las = 2)
    }
}

#' @rdname calibrate-calibrant-mass
#'
#' @description
#'
#' The `isCalibrated` function returns `TRUE` if chromatographic
#' peaks of the [XCMSnExp] object `x` were calibrated and `FALSE` otherwise.
#'
#' @md
isCalibrated <- function(object) {
    if (length(processHistory(object, type = .PROCSTEP.CALIBRATION)))
        TRUE
    else
        FALSE
}

#' @title Replace raw with adjusted retention times
#'
#' @description
#'
#' Replaces the raw retention times with the adjusted retention
#' time or returns the object unchanged if none are present.
#'
#' @details
#'
#' Adjusted retention times are stored *in parallel* to the adjusted
#' retention times in the `XCMSnExp`. The `applyAdjustedRtime` replaces the
#' raw retention times (stored in the *feature data* (`fData` `data.frame`))
#' with the adjusted retention times.
#'
#' @note
#'
#' Replacing the raw retention times with adjusted retention times
#' disables the possibility to restore raw retention times using the
#' [dropAdjustedRtime()] method. This function does **not** remove the
#' retention time processing step with the settings of the alignment from
#' the [processHistory()] of the `object` to ensure that the processing
#' history is preserved.
#'
#' @param object An [XCMSnExp] object.
#'
#' @md
#'
#' @return
#'
#' A `XCMSnExp` with the raw retention times being replaced with the
#' adjusted retention time.
#'
#' @author Johannes Rainer
#'
#' @seealso [adjustRtime()] for the function to perform the alignment (retention
#'     time correction).
#'
#'     [adjustedRtime()] for the method to extract adjusted retention times from
#'     an [XCMSnExp] object.
#'
#'     [dropAdjustedRtime] for the method to delete alignment results and to
#'     restore the raw retention times.
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
#' xod <- adjustRtime(faahko_sub, param = ObiwarpParam())
#'
#' hasAdjustedRtime(xod)
#'
#' ## Replace raw retention times with adjusted retention times.
#' xod <- applyAdjustedRtime(xod)
#'
#' ## No adjusted retention times present
#' hasAdjustedRtime(xod)
#'
#' ## Raw retention times have been replaced with adjusted retention times
#' plot(split(rtime(faahko_sub), fromFile(faahko_sub))[[1]] -
#'     split(rtime(xod), fromFile(xod))[[1]], type = "l")
#'
#' ## And the process history still contains the settings for the alignment
#' processHistory(xod)
applyAdjustedRtime <- function(object) {
    if (!is(object, "XCMSnExp"))
        stop("'object' has to be an 'XCMSnExp' object")
    if (!hasAdjustedRtime(object))
        return(object)
    ## Replace retention times.
    fData(object)$retentionTime <- rtime(object, adjusted = TRUE)
    ## Copy the data
    newFd <- new("MsFeatureData")
    newFd@.xData <- .copy_env(object@msFeatureData)
    rm(list = "adjustedRtime", envir = newFd)
    object@msFeatureData <- newFd
    object
}

## Find mz ranges with multiple peaks per sample.
## Use the density distribution for that? with a bandwidth = 0.001, check
## density method for that...


.concatenate_XCMSnExp <- function(...) {
    x <- list(...)
    if (length(x) == 0)
        return(NULL)
    if (length(x) == 1)
        return(x[[1]])
    ## Check that all are XCMSnExp objects.
    if (!all(unlist(lapply(x, function(z) is(z, "XCMSnExp")))))
        stop("All passed objects should be 'XCMSnExp' objects")
    new_x <- as(.concatenate_OnDiskMSnExp(...), "XCMSnExp")
    ## If any of the XCMSnExp has alignment results or detected features drop
    ## them!
    x <- lapply(x, function(z) {
        if (hasAdjustedRtime(z)) {
            z <- dropAdjustedRtime(z)
            warning("Adjusted retention times found, had to drop them.")
        }
        if (hasFeatures(z)) {
            z <- dropFeatureDefinitions(z)
            warning("Feature definitions found, had to drop them.")
        }
        z
    })
    ## Combine peaks
    fls <- lapply(x, fileNames)
    startidx <- cumsum(lengths(fls))
    pks <- lapply(x, chromPeaks)
    procH <- lapply(x, processHistory)
    for (i in 2:length(fls)) {
        pks[[i]][, "sample"] <- pks[[i]][, "sample"] + startidx[i - 1]
        procH[[i]] <- lapply(procH[[i]], function(z) {
            z@fileIndex <- as.integer(z@fileIndex + startidx[i - 1])
            z
            })
    }
    pks <- do.call(rbind, pks)
    rownames(pks) <- NULL
    new_x@.processHistory <- unlist(procH)
    chromPeaks(new_x) <- pks
    rm(pks)
    cpd <- do.call(rbind, lapply(x, chromPeakData))
    rownames(cpd) <- rownames(chromPeaks(new_x))
    chromPeakData(new_x) <- cpd
    validObject(new_x)
    new_x
}

#' @rdname XCMSnExp-filter-methods
filterFeatureDefinitions <- function(x, features) {
    if (!is(x, "XCMSnExp"))
        stop("'x' is expected to be an 'XCMSnExp' object")
    if (missing(features))
        return(x)
    if (!hasFeatures(x))
        stop("No feature definitions present! Run 'groupChromPeaks' first.")
    fts <- featureDefinitions(x)
    ## features input parameter checking:
    if (is.logical(features))
        features <- which(features)
    if (is.character(features))
        features <- match(features, rownames(fts))
    if (is.numeric(features)) {
        if (any(is.na(features)))
            stop("no 'NA' in 'features' allowed!")
        if (!all(features %in% 1:nrow(fts)))
            stop("specified 'features' out of range")
    } else {
        stop("'features' has to be either a 'integer' with the indices of the ",
             "features, a 'logical' or a 'character' matching rownames of the ",
             "'featureDefinitions' data frame.")
    }
    ## Actual sub-setting...
    newFd <- new("MsFeatureData")
    newFd@.xData <- .copy_env(x@msFeatureData)
    featureDefinitions(newFd) <- fts[features, , drop = FALSE]
    lockEnvironment(newFd, bindings = TRUE)
    x@msFeatureData <- newFd
    ## Add a generic filtering process history.
    x <- addProcessHistory(x, GenericProcessHistory(
                                  fun = "filterFeatureDefinitions",
                                  args = list(features = features),
                                  fileIndex. = 1:length(fileNames(x))))
    if (validObject(x))
        x
}

#' @title Simple feature summaries
#'
#' @description
#'
#' Simple function to calculate feature summaries. These include counts and
#' percentages of samples in which a chromatographic peak is present for each
#' feature and counts and percentages of samples in which more than one
#' chromatographic peak was annotated to the feature. Also relative standard
#' deviations (RSD) are calculated for the integrated peak areas per feature
#' across samples. For `perSampleCounts = TRUE` also the individual
#' chromatographic peak counts per sample are returned.
#'
#' @param x `XCMSnExp` object with correspondence results.
#'
#' @param group `numeric`, `logical`, `character` or `factor` with the same
#'     length than `x` has samples to aggregate counts by the groups defined
#'     in `group`.
#'
#' @param perSampleCounts `logical(1)` whether feature wise individual peak
#'     counts per sample should be returned too.
#'
#' @param method `character` passed to the [featureValues()] function. See
#'     respective help page for more information.
#'
#' @param skipFilled `logical(1)` whether filled-in peaks should be excluded
#'     (default) or included in the summary calculation.
#'
#' @return
#'
#' `matrix` with one row per feature and columns:
#'
#' - `"count"`: the total number of samples in which a peak was found.
#' - `"perc"`: the percentage of samples in which a peak was found.
#' - `"multi_count"`: the total number of samples in which more than one peak
#'   was assigned to the feature.
#' - `"multi_perc"`: the percentage of those samples in which a peak was found,
#'   that have also multiple peaks annotated to the feature. Example: for a
#'   feature, at least one peak was detected in 50 samples. In 5 of them 2 peaks
#'   were assigned to the feature. `"multi_perc"` is in this case 10%.
#' - `"rsd"`: relative standard deviation (coefficient of variation) of the
#'   integrated peak area of the feature's peaks.
#' - The same 4 columns are repeated for each unique element (level) in `group`
#'   if `group` was provided.
#'
#' If `perSampleCounts = TRUE` also one column for each sample is returned
#' with the peak counts per sample.
#'
#' @author Johannes Rainer
featureSummary <- function(x, group, perSampleCounts = FALSE,
                           method = "maxint", skipFilled = TRUE) {
    if (!is(x, "XCMSnExp"))
        stop("'x' is expected to be an 'XCMSnExp' object")
    if (!hasFeatures(x))
        stop("No feature definitions found in 'x'. Please perform first a ",
             "correspondence analysis with the 'groupChromPeaks' function.")
    if (!missing(group)) {
        if (length(group) != length(fileNames(x)))
            stop("length of 'group' does not match the number of ",
                 "samples in 'x'")
    }
    if (skipFilled && .hasFilledPeaks(x))
        x <- dropFilledChromPeaks(x)
    ## First determine the number of peaks per sample
    smpls <- seq_along(fileNames(x))
    pks_per_sample <- lapply(featureDefinitions(x)$peakidx, function(z)
        as.numeric(table(factor(chromPeaks(x)[z, "sample"],
                                levels = smpls))))
    pks_per_sample <- do.call(rbind, pks_per_sample)
    rownames(pks_per_sample) <- rownames(featureDefinitions(x))
    sum_fun <- function(z) {
        cnt <- sum(z > 0)
        cnt_multi <- sum(z > 1)
        sm <- c(count = cnt,
                 perc = cnt * 100 / length(z),
                 multi_count = cnt_multi,
                 multi_perc = cnt_multi * 100 / cnt
                 )
        sm[is.na(sm)] <- 0
        sm
    }
    fts_sum <- t(apply(pks_per_sample, MARGIN = 1, sum_fun))
    ## Calculate RSD
    rsd <- function(z) {sd(z, na.rm = TRUE) / abs(mean(z, na.rm = TRUE))}
    fvals <- featureValues(x, value = "into", method = method)
    fts_sum <- cbind(fts_sum, rsd = apply(fvals, MARGIN = 1, rsd))
    if (!missing(group)) {
        per_group <- lapply(unique(group), function(grp) {
            idx <- which(group == grp)
            res <- cbind(t(apply(pks_per_sample[, idx, drop = FALSE],
                                 MARGIN = 1, sum_fun)),
                         rsd = apply(fvals[, idx, drop = FALSE], MARGIN = 1,
                                     rsd))
            colnames(res) <- paste0(grp, "_", colnames(res))
            res
        })
        fts_sum <- cbind(fts_sum, do.call(cbind, per_group))
    }
    if (perSampleCounts) {
        colnames(pks_per_sample) <- basename(fileNames(x))
        fts_sum <- cbind(fts_sum, pks_per_sample)
    }
    fts_sum
}

#' @title Identify overlapping features
#'
#' @description
#'
#' `overlappingFeatures` identifies features that are overlapping or close in
#' the m/z - rt space.
#'
#' @param x `XCMSnExp` with the features.
#'
#' @param expandMz `numeric(1)` with the value to expand each feature (on each
#'     side) in m/z dimension before identifying overlapping features.
#'     The resulting `"mzmin"` for the feature is thus `mzmin - expandMz` and
#'     the `"mzmax"` `mzmax + expandMz`.
#'
#' @param expandRt `numeric(1)` with the value to expand each feature (on each
#'     side) in retention time dimension before identifying overlapping
#'     features. The resulting `"rtmin"` for the
#'     feature is thus `rtmin - expandRt` and the `"rtmax"` `rtmax + expandRt`.
#'
#' @param ppm `numeric(1)` to grow the m/z width of the feature by a relative
#'     value: `mzmin - mzmin * ppm / 2e6`, `mzmax + mzmax * ppm / 2e6`. Each
#'     feature is thus expanded in m/z dimension by ppm/2 on each side before
#'     identifying overlapping features.
#'
#' @return `list` with indices of features (in [featureDefinitions()]) that
#'     are overlapping.
#'
#' @md
#'
#' @author Johannes Rainer
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
#' ## Correspondence analysis
#' xdata <- groupChromPeaks(faahko_sub, param = PeakDensityParam(sampleGroups = c(1, 1, 1)))
#'
#' ## Identify overlapping features
#' overlappingFeatures(xdata)
#'
#' ## Identify features that are separated on retention time by less than
#' ## 2 minutes
#' overlappingFeatures(xdata, expandRt = 60)
overlappingFeatures <- function(x, expandMz = 0, expandRt = 0, ppm = 0) {
    if (!is(x, "XCMSnExp"))
        stop("'x' is expected to be an 'XCMSnExp' object")
    if (!hasFeatures(x))
        stop("No feature definitions found in 'x'. Please perform first a ",
             "correspondence analysis with the 'groupChromPeaks' function.")
    xl <- featureDefinitions(x)$rtmin
    xr <- featureDefinitions(x)$rtmax
    yb <- featureDefinitions(x)$mzmin
    yt <- featureDefinitions(x)$mzmax
    ## Expand them?
    if (expandMz != 0) {
        yb <- yb - expandMz
        yt <- yt + expandMz
    }
    if (ppm != 0) {
        yb <- yb - yb * ppm / 2e6
        yt <- yt + yt * ppm / 2e6
    }
    if (expandRt != 0) {
        rng <- range(rtime(x))
        xl <- xl - expandRt
        xl[xl < rng[1]] <- rng[1]
        xr <- xr + expandRt
        xr[xr > rng[2]] <- rng[2]
    }
    .rect_overlap(xleft = xl, xright = xr, ybottom = yb, ytop = yt)
}


#' @title Export data for use in MetaboAnalyst
#'
#' @description
#'
#' Export the feature table for further analysis in the MetaboAnalyst
#' software (or the `MetaboAnalystR` R package).
#'
#' @param x [XCMSnExp] object with identified chromatographic peaks grouped
#'     across samples.
#'
#' @param file `character(1)` defining the file name. If not specified, the
#'     `matrix` with the content is returned.
#'
#' @param label either `character(1)` specifying the phenodata column in `x`
#'     defining the sample grouping or a vector with the same length than
#'     samples in `x` defining the group assignment of the samples.
#'
#' @param value `character(1)` specifying the value to be returned for each
#'     feature. See [featureValues()] for more details.
#'
#' @param digits `integer(1)` defining the number of significant digits to be
#'     used for numeric. The default `NULL` uses `getOption("digits")`. See
#'     [format()] for more information.
#'
#' @param groupnames `logical(1)` whether row names of the resulting matrix
#'     should be the feature IDs (`groupnames = FALSE`; default) or IDs that
#'     are composed of the m/z and retention time of the features (in the
#'     format `M<m/z>T<rt>` (`groupnames = TRUE`). See help of the [groupnames]
#'     function for details.
#'
#' @param ... additional parameters to be passed to the [featureValues()]
#'     function.
#'
#' @return If `file` is not specified, the function returns the `matrix` in
#'     the format supported by MetaboAnalyst.
#'
#' @export
#'
#' @author Johannes Rainer
#'
#' @md
exportMetaboAnalyst <- function(x, file = NULL, label,
                                value = "into", digits = NULL,
                                groupnames = FALSE, ...) {
    if (!is(x, "XCMSnExp"))
        stop("'x' is supposed to be an XCMSnExp object")
    fv <- featureValues(x, value = value, ...)
    nas <- is.na(fv)
    fv <- format(fv, trim = TRUE, digits = digits)
    fv[nas] <- NA
    if (missing(label))
        stop("Please provide the group assignment of the samples with the ",
             "'label' parameter")
    if (groupnames)
        rownames(fv) <- groupnames(x)
    if (length(label) == 1) {
        if (any(colnames(pData(x)) == label))
            label <- as.character(pData(x)[, label])
        else
            stop("No pheno data column \"", label, "\" found. If length ",
                 "'label' is 1 it is expected to specify the column in ",
                 "the object's phenodata data.frame containing the",
                 " group assignment.")
    }
    if (length(label) != ncol(fv))
        stop("Length of 'label' has to match the number of samples (",
             ncol(fv), ").")
    fv <- rbind(Sample = colnames(fv),
                Label = as.character(label),
                fv)
    if (!is.null(file))
        write.table(fv, file = file, sep = ",",
                    qmethod = "double", col.names = FALSE, row.names = TRUE)
    else
        fv
}

#' @description
#'
#' Identifies for all peaks of an `XCMSnExp` object MS2 spectra and returns
#' them. Identification is performed separately (but parallel) for each file.
#'
#' @author Johannes Rainer
#'
#' @noRd
ms2_mspectrum_for_all_peaks <- function(x, expandRt = 0, expandMz = 0,
                                  ppm = 0, method = c("all",
                                                      "closest_rt",
                                                      "closest_mz",
                                                      "signal"),
                                  skipFilled = FALSE, subset = NULL,
                                  BPPARAM = bpparam()) {
    ## DEPRECATE THIS IN BIOC3.14
    method <- match.arg(method)
    pks <- chromPeaks(x)
    if (ppm != 0)
        ppm <- pks[, "mz"] * ppm / 1e6
    if (expandMz != 0 || length(ppm) > 1) {
        pks[, "mzmin"] <- pks[, "mzmin"] - expandMz - ppm
        pks[, "mzmax"] <- pks[, "mzmax"] + expandMz + ppm
    }
    if (expandRt != 0) {
        pks[, "rtmin"] <- pks[, "rtmin"] - expandRt
        pks[, "rtmax"] <- pks[, "rtmax"] + expandRt
    }
    if (length(subset)) {
        if (!(min(subset) >= 1 && max(subset) <= nrow(pks)))
            stop("If 'subset' is defined it has to be >= 1 and <= ",
                 nrow(pks), ".")
        not_subset <- rep(TRUE, nrow(pks))
        not_subset[subset] <- FALSE
        pks[not_subset, "mz"] <- NA
        rm(not_subset)
    }
    if (skipFilled && any(chromPeakData(x)$is_filled))
        pks[chromPeakData(x)$is_filled, ] <- NA
    ## Split data per file
    file_factor <- factor(pks[, "sample"])
    peak_ids <- rownames(pks)
    pks <- split.data.frame(pks, f = file_factor)
    x <- .split_by_file2(
        x, msLevel. = 2L, subsetFeatureData = FALSE)[as.integer(levels(file_factor))]
    res <- bpmapply(ms2_mspectrum_for_peaks_from_file, x, pks,
                    MoreArgs = list(method = method), SIMPLIFY = FALSE,
                    USE.NAMES = FALSE, BPPARAM = BPPARAM)
    res <- unsplit(res, file_factor)
    names(res) <- peak_ids
    res
}

#' @description
#'
#' Identify for all chromatographic peaks of a single file MS2 spectra with
#' an precursor m/z within the peak's m/z and retention time within the peak's
#' retention time width and return them.
#'
#' @note
#'
#' If needed, the m/z and rt width of the peaks should be increased previously.
#' No MS2 spectra are identified for peaks with an `"mz"` of `NA`, thus, to
#' skip identification for some (e.g. filled-in) peaks their value in the
#' `"mz"` column of `pks` should be set to `NA` before passing the parameter
#' to the function.
#'
#' @param x `OnDiskMSnExp` with (only) MS2 spectra of a single file.
#'
#' @param pks `matrix` with chromatographic peaks of a single file.
#'
#' @param method `character` defining the method to optionally select a single
#'     MS2 spectrum for the peak.
#'
#' @return `list` with length equal to the number of rows of `pks` with
#'     `Spectrum2` objects
#'
#' @author Johannes Rainer
#'
#' @noRd
ms2_mspectrum_for_peaks_from_file <- function(x, pks, method = c("all",
                                                                 "closest_rt",
                                                                 "closest_mz",
                                                                 "signal")) {
    ## DEPRECATE THIS IN BIOC3.14
    res <- vector(mode = "list", nrow(pks))
    if (nrow(pks) == 0 || !any(msLevel(x) == 2))
        return(res)
    method <- match.arg(method)
    fromFile <- as.integer(pks[1, "sample"])
    sps <- spectra(x)
    pmz <- precursorMz(x)
    rtm <- rtime(x)
    for (i in 1:nrow(pks)) {
        if (is.na(pks[i, "mz"]))
            next
        idx <- which(pmz >= pks[i, "mzmin"] & pmz <= pks[i, "mzmax"] &
                     rtm >= pks[i, "rtmin"] & rtm <= pks[i, "rtmax"])
        if (length(idx)) {
            if (length(idx) > 1 & method != "all") {
                if (method == "closest_rt")
                    idx <- idx[order(abs(rtm[idx] - pks[i, "rt"]))][1]
                if (method == "closest_mz")
                    idx <- idx[order(abs(pmz[idx] - pks[i, "mz"]))][1]
                if (method == "signal") {
                    sps_sub <- sps[idx]
                    ints <- vapply(sps_sub, function(z) sum(intensity(z)),
                                   numeric(1))
                    idx <- idx[order(abs(ints - pks[i, "maxo"]))][1]
                }
                if (method == "largest_tic") {
                    sps_sub <- sps[idx]
                    ints <- vapply(sps_sub, function(z) sum(intensity(z)),
                                   numeric(1))
                    idx <- idx[order(ints, decreasing = TRUE)][1L]
                }
                if (method == "largest_bpi") {
                    sps_sub <- sps[idx]
                    ints <- vapply(sps_sub, function(z) max(intensity(z)),
                                   numeric(1))
                    idx <- idx[order(ints, decreasing = TRUE)][1L]
                }
            }
            res[[i]] <- lapply(sps[idx], function(z) {
                z@fromFile = fromFile
                z
            })
        }
    }
    names(res) <- rownames(pks)
    res
}

#' given an XCMSnExp this function identifies for each MS1 chromatographic
#' peak all MS2 spectra with a precursor m/z within the peak region and returns
#' a `Spectra` object with all of these spectra.
#'
#' @return a `list`, same length than there are `chromPeaks` with a `Spectra`
#'     in each, if found.
#'
#' @noRd
.spectra_for_peaks <- function(x, method = c("all", "closest_rt",
                                                  "closest_mz", "signal",
                                                  "largest_tic", "largest_bpi"),
                               msLevel = 2L, expandRt = 0, expandMz = 0,
                               ppm = 0, skipFilled = FALSE,
                               peaks = character()) {
    if (is(x, "XCMSnExp") && hasAdjustedRtime(x))
        fData(x)$retentionTime <- rtime(x)
    ## from_msl <- 1L
    method <- match.arg(method)
    if (msLevel == 1L && method %in% c("closest_mz", "signal")) {
        warning("method = \"closest_mz\" and method = \"signa;\" are not",
                " supported for msLevel = 1. Changing to method = \"all\".")
        method <- "all"
    }
    pks <- as.data.frame(chromPeaks(x))[, c("mz", "mzmin", "mzmax", "rt",
                                            "rtmin", "rtmax", "maxo", "sample")]
    if (ppm != 0)
        expandMz <- expandMz + pks$mz * ppm / 1e6
    if (expandMz[1L] != 0) {
        pks$mzmin <- pks$mzmin - expandMz
        pks$mzmax <- pks$mzmax + expandMz
    }
    if (expandRt != 0) {
        pks$rtmin <- pks$rtmin - expandRt
        pks$rtmax <- pks$rtmax + expandRt
    }
    if (length(peaks)) {
        peaks <- .i2index(peaks, rownames(pks), "peaks")
        keep <- rep(FALSE, nrow(pks))
        keep[peaks] <- TRUE
    } else {
        keep <- rep(TRUE, nrow(pks))
        if (skipFilled && any(chromPeakData(x)$is_filled))
            keep <- !chromPeakData(x)$is_filled
        ## if (any(chromPeakData(x)$ms_level))
        ##     keep <- keep & chromPeakData(x)$ms_level == from_msl
    }
    ## maybe subset by peak ID.
    fns <- fileNames(x)
    res <- vector("list", length(keep))
    be <- new("MsBackendMzR")
    sps <- new("Spectra")
    for (i in which(keep)) {
        sel <- .fdata(x)$msLevel == msLevel &
               .fdata(x)$retentionTime >= pks$rtmin[i] &
               .fdata(x)$retentionTime <= pks$rtmax[i] &
               .fdata(x)$fileIdx == pks$sample[i]
        if (msLevel > 1L)
            sel <- sel &
                   .fdata(x)$precursorMZ >= pks$mzmin[i] &
                   .fdata(x)$precursorMZ <= pks$mzmax[i]
        fd <- .fdata(x)[which(sel), ]
        if (nrow(fd)) {
            fd$peak_index <- i
            fd$peak_id <- rownames(pks)[i]
            sp <- .fData2MsBackendMzR(fd, fns, be)
            sp <- switch(
                method,
                all = sp,
                closest_rt = {
                    sp[which.min(abs(pks$rt[i] - rtime(sp)))[1L]]
                },
                closest_mz = {
                    sp[which.min(abs(pks$mz[i] - precursorMz(sp)))[1L]]
                },
                signal = {
                    ints <- vapply(intensity(sp), sum, numeric(1))
                    sp[which.min(abs(ints - pks$maxo[i]))[1L]]
                },
                largest_tic = {
                    ints <- vapply(intensity(sp), sum, numeric(1))
                    sp[which.max(ints)[1L]]
                },
                largest_bpi = {
                    ints <- vapply(intensity(sp), max, numeric(1))
                    sp[which.max(ints)[1L]]
                })
            slot(sps, "backend", check = FALSE) <- sp
            res[[i]] <- sps
        }
    }
    names(res) <- rownames(pks)
    if (length(peaks))
        res[peaks]
    else res
}

#' @title Extract spectra associated with chromatographic peaks
#'
#' @description
#'
#' Extract (MS1 or MS2) spectra from an [XCMSnExp] object for each identified
#' chromatographic peak. The function returns by default spectra for
#' chromatographic peaks of **all** MS levels, but parameter `peaks` allows to
#' restrict the result to selected chromatographic peaks.
#' For `msLevel = 1L` (only supported for `return.type = "Spectra"` or
#' `return.type = "List"`) MS1 spectra within the retention time boundaries
#' (in the file in which the peak was detected) are returned. For
#' `msLevel = 2L` MS2 spectra are returned for a chromatographic
#' peak if their precursor m/z is within the retention time and m/z range of
#' the chromatographic peak. Parameter `method` allows to define whether all
#' or a single spectrum should be returned:
#'
#' - `method = "all"`: (default): return all spectra for each peak.
#' - `method = "closest_rt"`: return the spectrum with the retention time
#'   closest to the peak's retention time (at apex).
#' - `method = "closest_mz"`: return the spectrum with the precursor m/z
#'   closest to the peaks's m/z (at apex); only supported for `msLevel = 2L`.
#' - `method = "signal"`: return the spectrum with the sum of intensities most
#'   similar to the peak's apex signal (`"maxo"`); only supported for
#'   `msLevel = 2L`.
#' - `method = "largest_tic"`: return the spectrum with the largest total
#'   signal (sum of peaks intensities).
#' - `method = "largest_bpi"`: return the spectrum with the largest peak
#'   intensity (maximal peak intensity).
#'
#' Parameter `return.type` allows to specify the *type* of the result object.
#' Please use `return.type = "Spectra"` or `return.type = "List"`,
#' `return.type = "list"` or the default `return.type = "MSpectra"` will be
#' deprecated (also, they do not support extracting MS1 spectra).
#'
#' See also the *LC-MS/MS data analysis* vignette for more details and examples.
#'
#' @param x [XCMSnExp] object with identified chromatographic peaks.
#'
#' @param msLevel `integer(1)` defining whether MS1 or MS2 spectra should be
#'     returned. `msLevel = 1` is currently only supported for `return.type`
#'     being `"Spectra"` or `"List"`.
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
#'     parameter overrides `skipFilled` and is only supported for `return.type`
#'     being either `"Spectra"` or `"List"`.
#'
#' @param skipFilled `logical(1)` whether spectra for filled-in peaks should
#'     be reported or not.
#'
#' @param return.type `character(1)` defining the result type. Defaults to
#'     `return.type = "MSpectra"` but `return.type = "Spectra"` or
#'     `return.type = "List"` are preferred. See below for more information.
#'
#' @return
#'
#' parameter `return.type` allow to specify the type of the returned object:
#'
#' - `return.type = "MSpectra"`: a [MSpectra] object with elements being
#'   [Spectrum-class] objects. The result objects contains all spectra
#'   for all peaks. Metadata column `"peak_id"` provides the ID of the
#'   respective peak (i.e. its rowname in [chromPeaks()]).
#' - `return.type = "Spectra"`: a `Spectra` object (defined in the `Spectra`
#'   package). The result contains all spectra for all peaks. Metadata column
#'   `"peak_id"` provides the ID of the respective peak (i.e. its rowname in
#'   [chromPeaks()] and `"peak_index"` its index in the object's `chromPeaks`
#'   matrix.
#' - `return.type = "list"`: `list` of `list`s that are either of length
#'   0 or contain [Spectrum2-class] object(s) within the m/z-rt range. The
#'   length of the list matches the number of peaks.
#' - `return.type = "List"`: `List` of length equal to the number of
#'   chromatographic peaks is returned with elements being either `NULL` (no
#'   spectrum found) or a `Spectra` object.
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @examples
#'
#' ## Read a file with DDA LC-MS/MS data
#' fl <- system.file("TripleTOF-SWATH/PestMix1_DDA.mzML", package = "msdata")
#' dda <- readMSData(fl, mode = "onDisk")
#'
#' ## Subset the object to reduce runtime of the example
#' dda <- filterRt(dda, c(200, 400))
#'
#' ## Perform MS1 peak detection
#' dda <- findChromPeaks(dda, CentWaveParam(peakwidth = c(5, 15), prefilter = c(5, 1000)))
#'
#' ## Load the required Spectra package and return all MS2 spectro for each
#' ## chromatographic peaks as a Spectra object
#' ms2_sps <- chromPeakSpectra(dda, return.type = "Spectra")
#' ms2_sps
#'
#' ## columns peak_id or peak_index assign spectra to the chromatographic peaks
#' ms2_sps$peak_id
#' ms2_sps$peak_index
#' chromPeaks(dda)
#'
#' ## Alternatively, return the result as a List of Spectra objects. This list
#' ## is parallel to chromPeaks hence the mapping between chromatographic peaks
#' ## and MS2 spectra is easier.
#' ms2_sps <- chromPeakSpectra(dda, return.type = "List")
#' ms2_sps[[1L]]
#' length(ms2_sps)
#'
#' ## In addition to MS2 spectra we could also return the MS1 spectrum for each
#' ## chromatographic peak which is closest to the peak's apex position.
#' ms1_sps <- chromPeakSpectra(dda, msLevel = 1L, method = "closest_rt",
#'     return.type = "Spectra")
#' ms1_sps
#'
#' ## Parameter peaks would allow to extract spectra for specific peaks only
#' chromPeakSpectra(dda, msLevel = 1L, method = "closest_rt", peaks = c(3, 5))
chromPeakSpectra <- function(x, msLevel = 2L, expandRt = 0, expandMz = 0,
                             ppm = 0, method = c("all", "closest_rt",
                                                 "closest_mz", "signal",
                                                 "largest_tic", "largest_bpi"),
                             skipFilled = FALSE,
                             return.type = c("MSpectra", "Spectra",
                                             "list", "List"),
                             peaks = character()) {
    method <- match.arg(method)
    return.type <- match.arg(return.type)
    if (!is(x, "XCMSnExp"))
        stop("'x' is supposed to be an 'XCMSnExp' object.")
    if (!hasChromPeaks(x))
        stop("No chromatographic peaks present. Please run 'findChromPeaks' ",
             "first")
    if (return.type %in% c("Spectra", "List")) {
        .require_spectra()
        if (length(x@spectraProcessingQueue))
            warning("Lazy evaluation queue is not empty. Will ignore any",
                    " processing steps as 'return.type = \"Spectra\"' and",
                    " 'return.type = \"List\"' currently don't support a ",
                    "non-empty processing queue.")
        res <- .spectra_for_peaks(x, msLevel = msLevel, method = method,
                                  expandRt = expandRt, expandMz = expandMz,
                                  ppm = ppm, skipFilled = skipFilled,
                                  peaks = peaks)
        if (return.type == "Spectra") {
            res <- do.call(c, unname(res[lengths(res) > 0]))
            if (is(res, "Spectra"))
                res@processing <- character()
        } else res <- List(res)
    } else {
        ## DEPRECATE THIS IN BIOC 3.14
        if (length(peaks))
            warning("Ignoring parameter 'peaks' which is only supported for ",
                    "'return.type = \"Spectra\"' and 'return.type = \"List\"'.")
        if (msLevel != 2 || (msLevel == 2 & !any(msLevel(x) == 2))) {
            res <- vector(mode = "list", length = nrow(chromPeaks(x)))
            names(res) <- rownames(chromPeaks(x))
            if (msLevel != 2)
                warning("msLevel = 1 is currently not supported for",
                        "return.type = \"MSpectra\" or return.type = \"list\".")
            if (msLevel == 2 & !any(msLevel(x) == 2))
                warning("No MS2 spectra available in 'x'.")
        } else {
            res <- ms2_mspectrum_for_all_peaks(x, expandRt = expandRt,
                                               expandMz = expandMz,
                                               ppm = ppm, method = method,
                                               skipFilled = skipFilled)
        }
        if (return.type == "MSpectra") {
            pids <- rep(names(res), lengths(res))
            res <- res[lengths(res) > 0]
            if (length(res))
                res <- unlist(res)
            res <- MSpectra(res, elementMetadata = DataFrame(peak_id = pids))
        }
    }
    res
}

#' For information and details see featureSpectra
#'
#' @noRd
ms2_mspectrum_for_features <- function(x, expandRt = 0, expandMz = 0, ppm = 0,
                                       skipFilled = FALSE, ...) {
    idxs <- featureDefinitions(x)$peakidx
    sp_pks <- ms2_mspectrum_for_all_peaks(x, expandRt = expandRt,
                                          expandMz = expandMz, ppm = ppm,
                                          skipFilled = skipFilled,
                                          subset = unique(unlist(idxs)), ...)
    res <- lapply(idxs, function(z) unlist(sp_pks[z]))
    names(res) <- rownames(featureDefinitions(x))
    res
}

#' @title Extract spectra associated with features
#'
#' @description
#'
#' This function returns spectra associated with the identified features in the
#' input object. By default, spectra are returned for all features (from all
#' MS levels), but parameter `features` allows to specify selected features for
#' which the result should be returned.
#' Parameter `msLevel` allows to define whether MS level 1 or 2
#' spectra should be returned. For `msLevel = 1L` all MS1 spectra within the
#' retention time range of each chromatographic peak (in that respective data
#' file) associated with a feature are returned. Note that for samples in which
#' no peak was identified (or even filled-in) no spectra are returned.
#' For `msLevel = 2L` all MS2
#' spectra with a retention time within the retention time range and their
#' precursor m/z within the m/z range of any chromatographic peak of a feature
#' are returned. See also [chromPeakSpectra()] (used internally to extract
#' spectra for each chromatographic peak of a feature) for additional
#' information.
#'
#' In contrast to the [chromPeakSpectra()] function, selecting a `method`
#' different than `"all"` will not return a single spectrum per feature, but
#' one spectrum per **chromatographic peak** assigned to the feature.
#'
#' Note also that `msLevel = 1L` is only supported for `return.type = "List"`
#' or `return.type = "Spectra"`.
#'
#' @param x [XCMSnExp] object with feature defitions available.
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
#' parameter `return.type` allow to specify the type of the returned object:
#'
#' - `return.type = "MSpectra"`: a [MSpectra] object with elements being
#'   [Spectrum-class] objects. The result objects contains all spectra
#'   for all features. Metadata column `"feature_id"` provides the ID of the
#'   respective feature (i.e. its rowname in [featureDefinitions()]).
#' - `return.type = "Spectra"`: a `Spectra` object (defined in the `Spectra`
#'   package). The result contains all spectra for all features. Metadata column
#'   `"feature_id"` provides the ID of the respective feature (i.e. its rowname
#'   in [featureDefinitions()].
#' - `return.type = "list"`: `list` of `list`s that are either of length
#'   0 or contain [Spectrum2-class] object(s) within the m/z-rt range. The
#'   length of the list matches the number of features.
#' - `return.type = "List"`: `List` of length equal to the number of
#'   features with MS level `msLevel` is returned with elements being either
#'   `NULL` (no spectrum found) or a `Spectra` object.
#'
#' @author Johannes Rainer
#'
#' @md
featureSpectra <- function(x, msLevel = 2L, expandRt = 0, expandMz = 0,
                           ppm = 0, skipFilled = FALSE,
                           return.type = c("MSpectra", "Spectra",
                                           "list", "List"),
                           features = character(), ...) {
    if (!is(x, "XCMSnExp"))
        stop("'x' is supposed to be an 'XCMSnExp' object.")
    return.type <- match.arg(return.type)
    if (!hasFeatures(x))
        stop("No feature definitions present. Please run 'groupChromPeaks' ",
             "first.")
    if (return.type %in% c("Spectra", "List")) {
        .require_spectra()
        if (length(x@spectraProcessingQueue))
            warning("Lazy evaluation queue is not empty. Will ignore any",
                    " processing steps as 'return.type = \"Spectra\"' and",
                    " 'return.type = \"List\"' currently don't support a ",
                    "non-empty processing queue.")
        res <- .spectra_for_features(x, msLevel = msLevel, expandRt = expandRt,
                                     expandMz = expandMz, ppm = ppm,
                                     skipFilled = skipFilled,
                                     features = features, ...)
        if (return.type == "Spectra") {
            res <- do.call(c, unname(res[lengths(res) > 0]))
            if (!length(res)) {
                warning("No MS level ", msLevel, " spectra found")
                if (!is(res, "Spectra"))
                    res <- Spectra::Spectra()
            }
            res@processing <- character()
        } else res <- List(res)
    } else {
        ## DEPRECATE IN BIOC3.14
        if (length(features))
            warning("Ignoring parameter 'features': this is only supported",
                    " for 'return.type = \"Spectra\"' and ",
                    "'return.type = \"List\"'.")
        if (msLevel != 2 || (msLevel == 2 & !any(msLevel(x) == 2))) {
            res <- vector(mode = "list", length = nrow(featureDefinitions(x)))
            names(res) <- rownames(featureDefinitions(x))
            if (msLevel != 2)
                warning("msLevel = 1 is currently only supported for ",
                        "'return.type = \"Spectra\"' and ",
                        "'return.type = \"List\"'.")
            if (msLevel == 2 & !any(msLevel(x) == 2))
                warning("No MS2 spectra available in 'x'.")
        } else {
            res <- ms2_mspectrum_for_features(
                x, expandRt = expandRt, expandMz = expandMz,
                ppm = ppm, skipFilled = skipFilled, ...)
        }
        if (return.type == "MSpectra") {
            fids <- rep(names(res), lengths(res))
            res <- res[lengths(res) > 0]
            if (length(res)) {
                res <- unlist(res)
                pids <- vapply(strsplit(names(res), ".", TRUE),
                               `[`, character(1), 2)
            } else {
                pids <- character()
            }
            res <- MSpectra(res, elementMetadata = DataFrame(feature_id = fids,
                                                             peak_id = pids))
        }
    }
    res
}

#' get spectra per feature:
#' 1) call chromPeakSpectra to get spectra for all peaks (as List)
#' 2) iterate over features to extract all
#' @noRd
.spectra_for_features <- function(x, msLevel = 2L, expandRt = 0,
                                  expandMz = 0, ppm = 0, skipFilled = TRUE,
                                  features = character(), ...) {
    fids <- rownames(featureDefinitions(x))
    idx <- featureDefinitions(x)$peakidx
    if (length(features)) {
        findex <- .i2index(features, fids, "features")
        fids <- fids[findex]
        idx <- idx[findex]
    }
    ## Get spectra for all peaks of these features
    pkidx <- sort(unique(unlist(idx, use.names = FALSE)))
    peak_sp <- vector("list", nrow(chromPeaks(x)))
    peak_sp[pkidx] <- xcms:::.spectra_for_peaks(
        x, msLevel = msLevel, expandRt = expandRt, expandMz = expandMz,
        ppm = ppm, skipFilled = skipFilled, peaks = pkidx, ...)
    res <- lapply(seq_along(fids), function(i) {
        z <- peak_sp[idx[[i]]]
        if (any(lengths(z))) {
            z <- Spectra::concatenateSpectra(z)
            z@backend@spectraData <- cbind(z@backend@spectraData,
                                           DataFrame(feature_id = fids[i]))
            z@processing <- character()
            z
        }
    })
    names(res) <- fids
    res
}

#' @title Extract ion chromatograms for each feature
#'
#' @description
#'
#' Extract ion chromatograms for features in an [XCMSnExp-class] object. The
#' function returns for each feature its extracted ion chromatogram and all
#' associated peaks with it. The chromatogram is extracted from the m/z - rt
#' region including all chromatographic peaks of that features (i.e. based on
#' the ranges of `"mzmin"`, `"mzmax"`, `"rtmin"`, `"rtmax"` of all
#' chromatographic peaks of the feature).
#'
#' By default only chromatographic peaks associated with a feature are included
#' for an extracted ion chromatogram (parameter `include = "feature_only"`). By
#' setting `include = "apex_within"` all chromatographic peaks (and eventually
#' the feature which they are part of - if feature definitions are present)
#' that have their apex position within the m/z - rt range from which the
#' chromatogram is extracted are returned too.
#' With `include = "any"` or `include = "all"` all chromatographic peaks (and
#' eventually the feature in which they are present) overlapping the m/z and rt
#' range will be returned.
#'
#' @note
#'
#' When extracting EICs from only the top `n` samples it can happen that one
#' or more of the features specified with `features` are dropped because they
#' have no detected peak in the *top n* samples. The chance for this to happen
#' is smaller if `x` contains also filled-in peaks (with `fillChromPeaks`).
#'
#' @param x `XCMSnExp` object with grouped chromatographic peaks.
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
#' @param aggregationFun `character(1)` specifying the name that should be
#'     used to aggregate intensity values across the m/z value range for
#'     the same retention time. The default `"sum"` returns a base peak
#'     chromatogram.
#'
#' @param features `integer`, `character` or `logical` defining a subset of
#'     features for which chromatograms should be returned. Can be the index
#'     of the features in `featureDefinitions`, feature IDs (row names of
#'     `featureDefinitions`) or a logical vector.
#'
#' @param include `character(1)` defining which chromatographic peaks (and
#'     related feature definitions) should be included in the returned
#'     [XChromatograms()]. Defaults to `"feature_only"`; See description above
#'     for options and details.
#'
#' @param filled `logical(1)` whether filled-in peaks should be included in
#'     the result object. The default is `filled = FALSE`, i.e. only detected
#'     peaks are reported.
#'
#' @param n `integer(1)` to optionally specify the number of *top n* samples
#'     from which the EIC should be extracted.
#'
#' @param value `character(1)` specifying the column to be used to sort the
#'     samples. Can be either `"maxo"` (the default) or `"into"` to use the
#'     maximal peak intensity or the integrated peak area, respectively.
#'
#' @param ... optional arguments to be passed along to the [chromatogram()]
#'     function.
#'
#' @return [XChromatograms()] object.
#'
#' @md
#'
#' @seealso [filterColumnsKeepTop()] to filter the extracted EICs keeping only
#'     the *top n* columns (samples) with the highest intensity.
#'
#' @author Johannes Rainer
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
#' ## Subset the object to a smaller retention time range
#' xdata <- filterRt(faahko_sub, c(2500, 3500))
#'
#' xdata <- groupChromPeaks(xdata,
#'     param = PeakDensityParam(minFraction = 0.8, sampleGroups = rep(1, 3)))
#'
#' ## Get the feature definitions
#' featureDefinitions(xdata)
#'
#' ## Extract ion chromatograms for the first 3 features. Parameter
#' ## `features` can be either the feature IDs or feature indices.
#' chrs <- featureChromatograms(xdata, features = 1:3)
#'
#' ## Plot the XIC for the first feature using different colors for each file
#' plot(chrs[1, ], col = c("red", "green", "blue"))
featureChromatograms <- function(x, expandRt = 0, aggregationFun = "max",
                                 features,
                                 include = c("feature_only", "apex_within",
                                             "any", "all"),
                                 filled = FALSE,
                                 n = length(fileNames(x)),
                                 value = c("maxo", "into"),
                                 expandMz = 0, ...) {
    include <- match.arg(include)
    value <- match.arg(value)
    if (!hasFeatures(x))
        stop("No feature definitions present. Please run first 'groupChromPeaks'")
    if (!missing(features)) {
        features <- .i2index(features, ids = rownames(featureDefinitions(x)),
                             "features")
    } else features <- seq_len(nrow(featureDefinitions(x)))
    if (!length(features))
        return(MChromatograms(ncol = length(fileNames(x)), nrow = 0))
    ## If we want to get chromatograms only in a reduced number of samples
    n <- ceiling(n[1])
    if (n != length(fileNames(x))) {
        fids <- rownames(featureDefinitions(x))[features]
        if (n > length(fileNames(x)) || n < 1)
            stop("'n' has to be a positive integer between 1 and ",
                 length(fileNames(x)))
        pks <- chromPeaks(x)
        if (!filled)
            pks[chromPeakData(x)$is_filled, ] <- NA
        vals <- apply(.feature_values(
            pks, extractROWS(featureDefinitions(x), features),
            method = "maxint", value = value, intensity = "maxo",
            colnames = basename(fileNames(x))),
            MARGIN = 2, sum, na.rm = TRUE)
        sample_idx <- order(vals, decreasing = TRUE)[seq_len(n)]
        x <- .filter_file_XCMSnExp(x, sample_idx, keepFeatures = TRUE)
        features <- match(fids, rownames(featureDefinitions(x)))
        features <- features[!is.na(features)]
        if (length(features) < length(fids))
            warning(length(fids) - length(features), " of ", length(fids),
                    " features not present in the selected samples")
        if (!length(features))
            return(MChromatograms(ncol = length(fileNames(x)), nrow = 0))
    }
    pks <- chromPeaks(x)
    if (length(unique(rownames(pks))) != nrow(pks)) {
        rownames(pks) <- .featureIDs(nrow(pks), "CP")
        rownames(chromPeaks(x)) <- rownames(pks)
    }
    pk_idx <- featureDefinitions(x)$peakidx[features]
    rt_cols <- c("rtmin", "rtmax")
    mz_cols <- c("mzmin", "mzmax")
    mat <- do.call(rbind, lapply(pk_idx, function(z) {
        pks_current <- pks[z, , drop = FALSE]
        c(range(pks_current[, rt_cols]),
          range(pks_current[, mz_cols]))
    }))
    include_peaks <- "apex_within"
    if (include %in% c("any", "all"))
        include_peaks <- "any"
    mat[, 1] <- mat[, 1] - expandRt
    mat[, 2] <- mat[, 2] + expandRt
    mat[, 3] <- mat[, 3] - expandMz
    mat[, 4] <- mat[, 4] + expandMz
    colnames(mat) <- c("rtmin", "rtmax", "mzmin", "mzmax")
    chrs <- chromatogram(x, rt = mat[, 1:2], mz = mat[, 3:4],
                         aggregationFun = aggregationFun, filled = filled,
                         include = include_peaks, ...)
    if (include == "feature_only") {
        ## Loop over rows/features:
        ## subset to peaks of a feature.
        fts_all <- featureDefinitions(chrs)
        pks_all <- chromPeaks(chrs)
        chrs@featureDefinitions <- fts_all[integer(), ]
        nr <- nrow(chrs)
        nc <- ncol(chrs)
        ft_defs <- vector("list", nr)
        ft_ids <- rownames(featureDefinitions(x))[features]
        ## Keep only a single feature per row
        ## Keep only peaks for the features of interest.
        for (i in seq_len(nr)) {
            is_feature <- which(fts_all$row == i & rownames(fts_all) == ft_ids[i])
            if (!length(is_feature))
                next
            ft_def <- extractROWS(fts_all, is_feature)
            ft_defs[[i]] <- ft_def
            pk_ids <- rownames(pks_all)[ft_def$peakidx[[1]]]
            for (j in seq_len(nc)) {
                cur_pks <- chrs@.Data[i, j][[1]]@chromPeaks
                if (nrow(cur_pks)) {
                    keep <- which(rownames(cur_pks) %in% pk_ids)
                    chrs@.Data[i, j][[1]]@chromPeaks <-
                        cur_pks[keep, , drop = FALSE]
                    chrs@.Data[i, j][[1]]@chromPeakData <-
                        extractROWS(chrs@.Data[i, j][[1]]@chromPeakData, keep)
                }
            }
        }
        pks_sub <- chromPeaks(chrs)
        ## Update the index/mapping between features and peaks (in a loop to
        ## support duplicated features).
        fts <- lapply(seq_len(nr), function(r) {
            .subset_features_on_chrom_peaks(
                ft_defs[[r]], pks_all, pks_sub)
        })
        chrs@featureDefinitions <- do.call(rbind, fts)
    }
    if (validObject(chrs))
        chrs
}

#' @description
#'
#' \code{hasFilledChromPeaks}: whether filled-in peaks are present or not.
#'
#' @rdname XCMSnExp-class
hasFilledChromPeaks <- function(object) {
    .hasFilledPeaks(object)
}

#' Process the results from a peak detection in SWATH pockets.
#'
#' @param x `list` of `XCMSnExp` objects.
#'
#' @param msf `MsFeatureData` of the original object
#'
#' @param fileNames `character` with the file names of the original object. This
#'     is required to ensure that column `"sample"` in the chrom peaks matrix
#'     contains the correct indices.
#'
#' @return `MsFeatureData` with the `chromPeaks` and `chromPeakData` updated.
#'
#' @author Johannes Rainer
#'
#' @noRd
.swath_collect_chrom_peaks <- function(x, msf, fileNames) {
    pks <- do.call(rbind, lapply(x, function(z) {
        suppressWarnings(cpks <- chromPeaks(z))
        if (!is.null(cpks) && nrow(cpks))
            cpks[, "sample"] <- match(fileNames(z)[cpks[, "sample"]], fileNames)
        cpks
    }))
    cpd <- do.call(rbind, lapply(x, function(z) {
        if (hasChromPeaks(z) && nrow(chromPeakData(z))) {
            ret <- chromPeakData(z)
            target_mz <- isolationWindowTargetMz(z)[1]
            ret$isolationWindow <- .fdata(z)$isolationWindow[1]
            ret$isolationWindowTargetMZ <- target_mz
            ret$isolationWindowLowerMz <-
                target_mz - .fdata(z)$isolationWindowLowerOffset[1]
            ret$isolationWindowUpperMz <-
                target_mz + .fdata(z)$isolationWindowUpperOffset[1]
            ret
        } else DataFrame()
    }))
    if (!nrow(cpd))
        return(msf)
    if (hasChromPeaks(msf)) {
        idx_start <- max(nrow(chromPeaks(msf)),
                         as.numeric(sub("CP", "", rownames(chromPeaks(msf)))))
        rownames(pks) <- rownames(cpd) <- .featureIDs(nrow(pks),
                                                      from = idx_start + 1,
                                                      prefix = "CP")
        chromPeaks(msf) <- .rbind_fill(chromPeaks(msf), pks)
        chromPeakData(msf) <- .rbind_fill(chromPeakData(msf), cpd)
    } else {
        rownames(pks) <- rownames(cpd) <- .featureIDs(nrow(pks), prefix = "CP")
        chromPeaks(msf) <- pks
        chromPeakData(msf) <- cpd
    }
    msf
}

#' @title Data independent acquisition (DIA): peak detection in isolation windows
#'
#' @description
#'
#' The `findChromPeaksIsolationWindow` function allows to perform a
#' chromatographic peak detection in MS level > 1 spectra of certain isolation
#' windows (e.g. SWATH pockets). The function performs a peak detection,
#' separately for all spectra belonging to the same isolation window and adds
#' them to the [chromPeaks()] matrix of the result object, information about
#' the isolation window they were detected in is added to [chromPeakData()].
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
#' @param object `OnDiskMSnExp` or `XCMSnExp` object with the DIA data.
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
#' @param ... currently not used.
#'
#' @return
#'
#' An `XCMSnExp` object with the chromatographic peaks identified in spectra of
#' each isolation window from each file added to the `chromPeaks` matrix.
#' Isolation window definition for each identified peak are stored as additional
#' columns in [chromPeakData()].
#'
#' @author Johannes Rainer, Michael Witting
#'
#' @seealso [reconstructChromPeakSpectra()] for the function to reconstruct
#'     MS2 spectra for each MS1 chromatographic peak.
#'
#' @md
findChromPeaksIsolationWindow <-
    function(object, param, msLevel = 2L,
             isolationWindow = isolationWindowTargetMz(object), ...) {
        startDate <- date()
        if (!is.factor(isolationWindow))
            isolationWindow <- factor(isolationWindow)
        if (length(isolationWindow) != length(object))
            stop("length of 'isolationWindow' has to match length of 'object'")
        if (all(is.na(isolationWindow)))
            stop("all isolation windows in 'isolationWindow' are NA")
        if (!inherits(object, "OnDiskMSnExp"))
            stop("'object' should be an 'OnDiskMSnExp' or 'XCMSnExp' object")
        fData(object)$isolationWindow <- isolationWindow
        obj_sub <- selectFeatureData(as(object, "OnDiskMSnExp"),
                                     fcol = c(MSnbase:::.MSnExpReqFvarLabels,
                                              "centroided",
                                              "isolationWindow",
                                              "isolationWindowTargetMZ",
                                              "isolationWindowLowerOffset",
                                              "isolationWindowUpperOffset"))
        if (inherits(object, "XCMSnExp"))
            fData(obj_sub)$retentionTime <- rtime(object)
        res <- lapply(split(obj_sub, f = isolationWindow),
                      FUN = findChromPeaks, param = param, msLevel = msLevel)
        if (!inherits(object, "XCMSnExp"))
            object <- as(object, "XCMSnExp")
        msf <- new("MsFeatureData")
        msf@.xData <- .copy_env(object@msFeatureData)
        msf <- .swath_collect_chrom_peaks(res, msf, fileNames(object))
        lockEnvironment(msf, bindings = TRUE)
        object@msFeatureData <- msf
        xph <- XProcessHistory(param = param, date. = startDate,
                               type. = .PROCSTEP.PEAK.DETECTION,
                           fileIndex = 1:length(fileNames(object)),
                           msLevel = msLevel)
        object@.processHistory <- c(processHistory(object), list(xph))
        validObject(object)
        object
    }

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
#' The resulting `MSpectra` object provides also the peak IDs of the MS2
#' chromatographic peaks for each spectrum as well as their correlation value.
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
#'     Can be either `return.type = "MSpectra"` (the default) to return a
#'     `MSnbase::MSpectra` object or `return.type = "Spectra"` for the newer
#'     `Spectra::Spectra` object.
#'
#' @return
#'
#' Depending on `return.type`:
#'
#' - [MSpectra()] with the reconstructed MS2 spectra for all MS1 peaks
#'   in `object`. Contains empty [Spectrum2-class] objects for MS1 peaks for
#'   which reconstruction was not possible (either no MS2 signal was recorded
#'   or the correlation of the MS2 chromatographic peaks with the MS1
#'   chromatographic peak was below threshold `minCor`. `MSpectra` metadata
#'   columns `"ms2_peak_id"` and `"ms2_peak_cor"` (of type [CharacterList()]
#'   and [NumericList()] with length equal to the number of peaks per
#'   reconstructed MS2 spectrum) providing the IDs and the correlation of the
#'   MS2 chromatographic peaks from which the MS2 spectrum was reconstructed.
#'   As retention time the median retention times of all MS2 chromatographic
#'   peaks used for the spectrum reconstruction is reported. The MS1
#'   chromatographic peak intensity is reported as the reconstructed
#'   spectrum's `precursorIntensity` value (see parameter `intensity` above).
#' - `Spectra` object (defined in the `Spectra` package). The same content and
#'   information than above.
#'
#' @author Johannes Rainer, Michael Witting
#'
#' @md
#'
#' @seealso [findChromPeaksIsolationWindow()] for the function to perform MS2
#'     peak detection in DIA isolation windows and for examples.
reconstructChromPeakSpectra <- function(object, expandRt = 0, diffRt = 2,
                                        minCor = 0.8, intensity = "maxo",
                                        peakId = rownames(
                                            chromPeaks(object, msLevel = 1L)),
                                        BPPARAM = bpparam(),
                                        return.type = c("MSpectra", "Spectra")) {
    if (!inherits(object, "XCMSnExp") || !hasChromPeaks(object))
        stop("'object' should be an 'XCMSnExp' object with identified ",
             "chromatographic peaks")
    if (!is.character(peakId))
        stop("'peakId' has to be of type character")
    return.type <- match.arg(return.type)
    n_peak_id <- length(peakId)
    peakId <- intersect(peakId, rownames(chromPeaks(object, msLevel = 1L)))
    if (!length(peakId))
        stop("None of the provided 'peakId' matches IDs of MS1 ",
             "chromatographic peaks")
    if (length(peakId) < n_peak_id)
        warning("Only ", length(peakId), " of the provided",
                " identifiers match IDs of MS1 chromatographic peaks")
    object <- selectFeatureData(object,
                                fcol = c(MSnbase:::.MSnExpReqFvarLabels,
                                         "centroided",
                                         "polarity",
                                         "isolationWindow",
                                         "isolationWindowTargetMZ",
                                         "isolationWindowLowerOffset",
                                         "isolationWindowUpperOffset"))
    sps <- bplapply(
        .split_by_file2(
            object, subsetFeatureData = FALSE, to_class = "XCMSnExp"),
        FUN = function(x, files, expandRt, diffRt, minCor, col, pkId,
                       return.type) {
            .reconstruct_ms2_for_peaks_file(
                x, expandRt = expandRt, diffRt = diffRt,
                minCor = minCor, fromFile = match(fileNames(x), files),
                column = col, peakId = pkId, return.type = return.type)
        },
        files = fileNames(object), expandRt = expandRt, diffRt = diffRt,
        minCor = minCor, col = intensity, pkId = peakId, BPPARAM = BPPARAM,
        return.type = return.type)
    do.call(c, sps)
}

#' This function *overwrites* the `MSnbase` .plot_XIC function by adding also
#' a rectangle with the identified chromatographic peak.
#'
#' @param x `XCMSnExp` object with identifie chromatographic peaks.
#'
#' @param peakCol color for the border of the rectangle.
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @noRd
.plot_XIC <- function(x, peakCol = "#00000060", ...) {
    peakCol <- peakCol[1]
    x <- filterMsLevel(x, 1L)
    if (!length(x))
        stop("No MS1 data available")
    pks <- chromPeaks(x)
    fls <- basename(fileNames(x))
    x <- .split_by_file2(x)
    x <- lapply(x, as, "data.frame")
    ## Check if we are greedy and plot a too large area
    if (any(unlist(lapply(x, nrow)) > 20000))
        warning("The MS area to be plotted seems rather large. It is suggested",
                " to restrict the data first using 'filterRt' and 'filterMz'. ",
                "See also ?chromatogram and ?Chromatogram for more efficient ",
                "functions to plot a total ion chromatogram or base peak ",
                "chromatogram.",
                immediate = TRUE, call = FALSE)
    ## Define the layout.
    dots <- list(...)
    if (any(names(dots) == "layout")) {
        if (!is.null(dots$layout))
            layout(layout)
        dots$layout <- NULL
    } else
        layout(MSnbase:::.vertical_sub_layout(length(x)))
    for (i in seq_along(x)) {
        do.call(MSnbase:::.plotXIC,
                c(list(x = x[[i]], main = fls[i], layout = NULL), dots))
        pks_current <- pks[pks[, "sample"] == i, , drop = FALSE]
        if (length(pks_current) && nrow(pks_current)) {
            do.call(rect, c(list(xleft = pks_current[, "rtmin"],
                                 ybottom = pks_current[, "mzmin"],
                                 xright = pks_current[, "rtmax"],
                                 ytop = pks_current[, "mzmax"],
                                 border = peakCol), dots))
        }
    }
}

#' @title Group chromatographic peaks based on m/z or retention time
#'
#' @description
#'
#' Group chromatographic peaks if they are overlapping on m/z (independently
#' of their retention time) or *vice versa*.
#'
#' @param x `matrix` with columns `"mzmin"` and `"mzmax"` or `"rtmin"` and
#'     `"rtmax"`.
#'
#' @param min_col
#'
#' @param max_col `character(1)` with the name of the column with the upper
#'     range (e.g. `"mzmax"` or `"rtmax"`).
#'
#' @param min_col `character(1)` with the name of the column with the lower
#'     range (e.g. `"mzmin"` or `"rtmin"`).
#'
#' @param expand `numeric(1)` defining a constant value by which each e.g. m/z
#'     range is supposed to be expanded. Note that the mz range will be expanded
#'     by `expandMz` in both dimensions (i.e. `"mzmin"` - `expandMz` and
#'     `"mzmax"` + `expandMz`.
#'
#' @param ppm `numeric(1)` defining an m/z relative value by which the m/z range
#'     should be expanded.
#'
#' @note
#'
#' `x` is supposed to be a `chromPeaks` matrix of a single file, otherwise we're
#' grouping chromatographic peaks across samples.
#'
#' Note also that **each** peak gets expanded by `expandMz`, thus
#' peaks differing by `2 * expand` will be overlapping. As an example: m/z max
#' of one peak is 12.2, m/z min of another one is 12.4, if `expand = 0.1` is
#' used the m/z max of the first peak will be 12.3 and the m/z min of the second
#' one 12.3, thus both are considered *overlapping*.
#'
#' @author Johannes Rainer
#'
#' @return `list` with rownames (chromatographic peak IDs) of peak groups.
#'
#' @noRd
#'
#' @examples
#'
#' mat <- cbind(rtmin = c(10, 13, 16, 18), rtmax = c(12, 15, 17, 20),
#'     mzmin = c(2, 3, 4, 7), mzmax = c(2.5, 3.5, 4.2, 7.6))
#' rownames(mat) <- c("a", "b", "c", "d")
#' .group_overlapping_peaks(mat)
#'
#' .group_overlapping_peaks(mat, expand = 1)
#'
#' .group_overlapping_peaks(mat, expand = 0.25)
.group_overlapping_peaks <- function(x, min_col = "mzmin", max_col = "mzmax",
                                     expand = 0, ppm = 0) {
    x[, min_col] <- x[, min_col] - expand - x[, min_col] * ppm / 1e6
    x[, max_col] <- x[, max_col] + expand + x[, max_col] * ppm / 1e6
    reduced_ranges <- .reduce(x[, min_col], x[, max_col])
    res <- vector("list", nrow(reduced_ranges))
    tolerance <- sqrt(.Machine$double.eps)
    for (i in seq_along(res)) {
        res[[i]] <- rownames(x)[
            x[, min_col] >= reduced_ranges[i, 1] - tolerance &
            x[, max_col] <= reduced_ranges[i, 2] + tolerance
        ]
    }
    res
}

#' @description
#'
#' Identify chromatographic peaks overlapping in m/z dimension and being close
#' on retention time to combine them if they fulfill the additional criteria:
#' intensity at `"rtmax"` for the first chromatographic peak is
#' `> prop * "maxo"` (`"maxo"` being the maximal intensity of the first peak)
#' **and** intensity at `"rtmin"` for the second chromatographic peak is
#' `> prop * "maxo"` of the second peak.
#'
#' @details
#'
#' The function first identifies chromatographic peaks within the same sample
#' that are overlapping on their m/z range. The m/z range can be expanded with
#' parameter `expandMz` or `ppm`. Note that both the upper and lower m/z is
#' expanded by these resulting in m/z ranges that are of size *original m/z
#' range* `+ 2 * expandMz`.
#'
#' All peaks are first ordered by theyr `"mzmin"` and subsequently expanded by
#' `expandMz` and `ppm`. Peaks are grouped if their expanded m/z ranges
#' (`"mzmin" - expandMz - ppm("mzmin")` to `"mzmax + expandMz + ppm("mzmax")`)
#' and rt ranges (`"rtmin" - expandRt` to `"rtmax" + expandRt`) are overlapping.
#'
#' For overlapping peak-candidates a chromatogram is extracted, with the
#' m/z range being the range of the individual chromatographic peak's m/z range
#' expanded by `expandMz` and `ppm` (on both sides). This is to avoid data
#' points in between peaks being `NA`.
#'
#' @param x `XCMSnExp` object with chromatographic peaks of a **single** file or
#'     an `OnDiskMSnExp` object, in which case parameters `pks` and `pkd` have
#'     to be provided.
#'
#' @param pks `chromPeaks` matrix.
#'
#' @param pkd `chromPeakData` data frame.
#'
#' @param sample_index `integer(1)` representing the index of the sample in the
#'     original object. To be used in column `"sample"` of the new peaks.
#'
#' @param expandRt `numeric(1)` defining by how many seconds the retention time
#'     window is expanded on both sides to check for overlapping peaks.
#'
#' @param expandMz `numeric(1)` constant value by which the m/z range of each
#'     chromatographic peak should be expanded (on both sides!) to check for
#'     overlapping peaks.
#'
#' @param ppm `numeric(1)` defining a m/z relative value (in parts per million)
#'     by which the m/z range of each chromatographic peak should be expanded
#'     to check for overlapping peaks.
#'
#' @param minProp `numeric(1)` between `0` and `1` representing the proporion
#'     of intensity to be required for peaks to be joined. See description for
#'     more details. The default (`minProp = 0.75`) means that peaks are only
#'     joined if the signal half way between then is larger 75% of the smallest
#'     of the two peak's `"maxo"` (maximal intensity at peak apex).
#'
#' @return `list` with element `"chromPeaks"`, that contains the peaks `matrix`
#'     containing newly merged peaks and original peaks if they could not be
#'     merged and `"chromPeakData"` that represents the `DataFrame` with the
#'     corresponding metadata information. The merged peaks will have a row
#'     name of `NA`.
#'
#' @author Johannes Rainer, Mar Garcia-Aloy
#'
#' @md
#'
#' @noRd
#'
#' @examples
#'
#' xd <- readMSData(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
#'     mode = "onDisk")
#' xd <- findChromPeaks(xd, param = CentWaveParam())
#'
#' xchr <- chromatogram(xd, mz = c(-0.5, 0.5) + 305.1)
#' plot(xchr)
#'
#' res <- xcms:::.merge_neighboring_peaks(xd, expandRt = 4)
#'
#' res_sub <- res[res[, "mz"] >= 305.05 & res[, "mz"] <= 305.15, ]
#' rect(res_sub[, "rtmin"], 0, res_sub[, "rtmax"], res_sub[, "maxo"],
#'     border = "red")
#'
#' xchr <- chromatogram(xd, mz = c(-0.5, 0.5) + 496.2)
#' plot(xchr)
#'
#' res <- xcms:::.merge_neighboring_peaks(xd, expandRt = 4)
#'
#' res_sub <- res[res[, "mz"] >= 496.15 & res[, "mz"] <= 496.25, ]
#' rect(res_sub[, "rtmin"], 0, res_sub[, "rtmax"], res_sub[, "maxo"],
#'     border = "red")
.merge_neighboring_peaks <- function(x, pks = chromPeaks(x),
                                     pkd = chromPeakData(x),
                                     expandRt = 2, expandMz = 0, ppm = 10,
                                     minProp = 0.75) {
    pks <- force(pks)
    pkd <- force(pkd)
    if (is.null(rownames(pks)))
        stop("Chromatographic peak IDs are required.")
    ms_level <- unique(pkd$ms_level)
    if (length(ms_level) != 1)
        stop("Got chromatographic peaks from different MS levels.", call. = FALSE)
    mz_groups <- .group_overlapping_peaks(pks, expand = expandMz, ppm = ppm)
    mz_groups <- mz_groups[lengths(mz_groups) > 1]
    message("Evaluating ", length(mz_groups), " peaks in file ",
            basename(fileNames(x)), " for merging ... ", appendLF = FALSE)
    if (!length(mz_groups)) {
        message("OK")
        return(list(chromPeaks = pks, chromPeakData = pkd))
    }
    ## Defining merge candidates
    pk_groups <- list()
    chr_def_mat <- list()
    current_group <- 1
    for (i in seq_along(mz_groups)) {
        rt_groups <- .group_overlapping_peaks(
            pks[mz_groups[[i]], , drop = FALSE],
            expand = expandRt, min_col = "rtmin",
            max_col = "rtmax"
        )
        rt_groups <- rt_groups[lengths(rt_groups) > 1]
        for (j in seq_along(rt_groups)) {
            rt_group <- rt_groups[[j]]
            pk_groups[[current_group]] <- rt_group
            pks_sub <- pks[rt_group, ]
            mzr_sub <- range(pks_sub[, c("mzmin", "mzmax")])
            ## Expand the mz range in a similar fashion than used for checking
            ## if chrom peaks are overlapping. This fixes an issue with very
            ## low intensities in between two peaks, that tend to have shifted
            ## m/z value (because their intensities are so low).
            mzr_sub <- mzr_sub + c(-1, 1) * mzr_sub * ppm * 1e-6 + expandMz
            chr_def_mat[[current_group]] <-
                c(mzr_sub,
                  range(pks_sub[, c("rtmin", "rtmax")]))
            current_group <- current_group + 1
        }
    }
    if (!length(chr_def_mat)) {
        message("OK")
        return(list(chromPeaks = pks, chromPeakData = pkd))
    }
    chr_def_mat <- do.call(rbind, chr_def_mat)
    chrs <- chromatogram(x, mz = chr_def_mat[, c(1, 2)],
                         rt = chr_def_mat[, c(3, 4)],
                         msLevel = ms_level)
    ## Now proceed to process them.
    res_list <- pkd_list <- vector("list", length(pk_groups))
    for (i in seq_along(pk_groups)) {
        pk_group <- pk_groups[[i]]
        res <- .chrom_merge_neighboring_peaks(
            chrs[i, 1], pks = pks[pk_group, , drop = FALSE],
            extractROWS(pkd, pk_group), diffRt = 2 * expandRt,
            minProp = minProp)
        res_list[[i]] <- res$chromPeaks
        pkd_list[[i]] <- res$chromPeakData
    }
    pks_new <- do.call(rbind, res_list)
    pks_new[, "sample"] <- pks[1, "sample"]
    pkd_new <- do.call(rbind, pkd_list)
    ## drop peaks that were candidates, but that were not returned (i.e.
    ## were either merged or dropped)
    keep <- which(!(rownames(pks) %in% setdiff(unlist(pk_groups,
                                                      use.names = FALSE),
                                               rownames(pks_new))))
    pks <- pks[keep, , drop = FALSE]
    pkd <- extractROWS(pkd, keep)
    ## add merged peaks
    idx_new <- is.na(rownames(pks_new))
    message("OK")
    list(chromPeaks = rbind(pks, pks_new[idx_new, , drop = FALSE]),
         chromPeakData = rbind(pkd, extractROWS(pkd_new, idx_new)))
}

.filter_file_XCMSnExp <- function(object, file,
                                  keepAdjustedRtime = hasAdjustedRtime(object),
                                  keepFeatures = FALSE) {
    if (missing(file)) return(object)
    if (is.character(file))
        file <- base::match(file, basename(fileNames(object)))
    ## This will not work if we want to get the files in a different
    ## order (i.e. c(3, 1, 2, 5))
    file <- base::sort(unique(file))
    ## Error checking - seems that's not performed downstream.
    if (!all(file %in% seq_along(fileNames(object))))
        stop("'file' has to be within 1 and the number of files in object.")
    has_features <- hasFeatures(object)
    has_chrom_peaks <- hasChromPeaks(object)
    has_adj_rt <- hasAdjustedRtime(object)
    ph <- processHistory(object)
    if (has_features && !keepFeatures) {
        has_features <- FALSE
        ph <- dropProcessHistoriesList(ph, .PROCSTEP.PEAK.GROUPING, num = 1)
    }
    ## Extracting all the XCMSnExp data from the object.
    newFd <- new("MsFeatureData")
    ## Subset original data:
    nobject <- as(filterFile(as(object, "OnDiskMSnExp"), file = file),
                 "XCMSnExp")
    if (has_adj_rt)
        adjustedRtime(newFd) <- adjustedRtime(object, bySample = TRUE)[file]
    if (has_chrom_peaks) {
        pks <- chromPeaks(object)
        idx <- pks[, "sample"] %in% file
        pks <- pks[idx, , drop = FALSE]
        pks[, "sample"] <- match(pks[, "sample"], file)
        if (has_features) {
            featureDefinitions(newFd) <- .update_feature_definitions(
                featureDefinitions(object),
                original_names = rownames(chromPeaks(object)),
                subset_names = rownames(pks))
        }
        chromPeaks(newFd) <- pks
        chromPeakData(newFd) <- extractROWS(chromPeakData(object), idx)
    }
    if (hasAdjustedRtime(newFd) && !keepAdjustedRtime)
        newFd <- dropAdjustedRtime(newFd, rtime(nobject, bySample = TRUE,
                                                adjusted = FALSE))
    ## Remove ProcessHistory not related to any of the files.
    ## if (length(ph)) {
    ##     kp <- unlist(lapply(ph, function(z) {
    ##         any(fileIndex(z) %in% file)
    ##     }))
    ##     ph <- ph[kp]
    ## }
    ## ## Update file index in process histories.
    ## if (length(ph)) {
    ##     ph <- lapply(ph, function(z) {
    ##         updateFileIndex(z, old = file, new = 1:length(file))
    ##     })
    ## }
    lockEnvironment(newFd, bindings = TRUE)
    nobject@msFeatureData <- newFd
    nobject@.processHistory <- ph
    nobject
}

#' Define the MS region (m/z - rt range) for each feature based on the rtmin,
#' rtmax, mzmin, mzmax of the corresponding detected peaks.
#'
#' @param x `XCMSnExp` object
#'
#' @param mzmin, mzmax, rtmin, rtmax `function` to be applied to the values
#'     (rtmin, ...) of the chrom peaks. Defaults to `median` but would also
#'     work with `mean` etc.
#'
#' @return `matrix` with columns `"mzmin"`, `"mzmax"`, `"rtmin"`, `"rtmax"`
#'     defining the range of
#'
#' @author Johannes Rainer
#'
#' @noRd
.features_ms_region <- function(x, mzmin = median, mzmax = median,
                                rtmin = median, rtmax = median,
                                msLevel = unique(msLevel(x)),
                                features = character()) {
    pk_idx <- featureValues(x, value = "index", method = "maxint",
                            msLevel = msLevel)
    if (length(features)) {
        if (!all(features %in% rownames(pk_idx)))
            stop(sum(!features %in% rownames(pk_idx)), " IDs defined with ",
                 "'features' are not available in 'object' for MS level ",
                 msLevel)
        pk_idx <- pk_idx[features, , drop = FALSE]
    }
    n_ft <- nrow(pk_idx)
    rt_min <- rt_max <- mz_min <- mz_max <- numeric(n_ft)
    for (i in seq_len(n_ft)) {
        idx <- pk_idx[i, ]
        tmp_pks <- chromPeaks(x)[idx[!is.na(idx)], , drop = FALSE]
        rt_min[i] <- rtmin(tmp_pks[, "rtmin"])
        rt_max[i] <- rtmax(tmp_pks[, "rtmax"])
        mz_min[i] <- mzmin(tmp_pks[, "mzmin"])
        mz_max[i] <- mzmax(tmp_pks[, "mzmax"])
    }
    res <- cbind(mzmin = mz_min, mzmax = mz_max, rtmin = rt_min, rtmax = rt_max)
    rownames(res) <- rownames(pk_idx)
    res
}

#' @rdname XCMSnExp-class
#'
#' @description
#'
#' \code{featureArea} extracts the m/z - retention time region for each feature.
#' This area is defined by the m/z - retention time regions of all
#' chromatographic peaks associated with a feature. Parameters \code{mzmin},
#' \code{mzmax}, \code{rtmin} and \code{rtmax} allow to define functions how
#' the corresponding value is calculated from the individual values (such as
#' the \code{"rtmin"}) of all chromatographic peaks of that feature. By default
#' the median \code{"rtmin"}, \code{"rtmax"}, \code{"mzmin"} and \code{"mzmax"}
#' is reported. Parameter \code{features} allows to provide feature IDs for
#' which the area should be extracted. By default it is extracted for all
#' features.
#'
#' @param features for \code{featureArea}: IDs of features for which the area
#'     should be extracted.
#'
#' @param mzmin for \code{featureArea}: \code{function} to be applied to values
#'     in the \code{"mzmin"} column of all chromatographic peaks of a feature
#'     to define the lower m/z value of the feature area.
#'     Defaults to \code{median}.
#'
#' @param mzmax for \code{featureArea}: \code{function} same as \code{mzmin}
#'     but for the \code{"mzmax"} column.
#'
#' @param rtmin for \code{featureArea}: \code{function} same as \code{mzmin}
#'     but for the \code{"rtmin"} column.
#'
#' @param rtmax for \code{featureArea}: \code{function} same as \code{mzmin}
#'     but for the \code{"rtmax"} column.
featureArea <- function(object, mzmin = median, mzmax = median, rtmin = median,
                        rtmax = median, msLevel = unique(msLevel(object)),
                        features = character()) {
    if (!hasFeatures(object))
        stop("No correspondence results available. Please run ",
             "'groupChromPeaks' first.")
    .features_ms_region(object, mzmin = mzmin, mzmax = mzmax, rtmin = rtmin,
                        rtmax = rtmax, msLevel = force(msLevel),
                        features = features)
}

#' @param x `XCMSnExp` object of a single file.
#'
#' @param nValues `integer(1)` defining the number of values that have to be above
#'     threshold.
#'
#' @param threshold `numeric(1)` with the threshold value.
#'
#' @param msLevel `integer(1)` with the MS level.
#'
#' @return `integer` with the index of the peaks that should be retained.
#'
#' @author Johannes Rainer
#'
#' @noRd
.chrom_peaks_above_threshold <- function(x, nValues = 1L, threshold = 0,
                                         msLevel = 1L) {
    if (length(msLevel) > 1)
        stop("Currently only filtering of a single MS level at a ",
             "time is supported")
    x <- applyAdjustedRtime(x)
    chrs <- chromatogram(
        filterMsLevel(as(x, "OnDiskMSnExp"), msLevel = msLevel),
        rt = chromPeaks(x, msLevel = msLevel)[, c("rtmin", "rtmax")],
        mz = chromPeaks(x, msLevel = msLevel)[, c("mzmin", "mzmax")],
        msLevel = msLevel)
    keep <- vapply(chrs@.Data, function(z) {
        sum(z@intensity >= threshold, na.rm = TRUE) >= nValues
    }, logical(1))
    keep_all <- rep(TRUE, nrow(chromPeaks(x)))
    keep_all[chromPeakData(x)$ms_level == msLevel] <- keep
    keep_all
}

#' @title Manual peak integration and feature definition
#'
#' @description
#'
#' The `manualChromPeaks` function allows to manually define chromatographic
#' peaks which are added to the object's `chromPeaks` matrix. In contrast to
#' [findChromPeaks()], no *peak detection* is performed (e.g. using an
#' algorithm such as *centWave*) but the peak is added as defined by the user.
#' Note that a peak will not be added if no signal (intensity) was found in a
#' sample within the provided boundaries.
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
#' @param object `XCMSnExp` or `OnDiskMSnExp` object.
#'
#' @param chromPeaks `matrix` defining the boundaries of the chromatographic
#'     peaks, one row per chromatographic peak, columns `"mzmin"`, `"mzmax"`,
#'     `"rtmin"` and `"rtmax"` defining the m/z and retention time region of
#'     each peak.
#'
#' @param peakIdx for `nabbyakFeatyres`: `list` of `integer` vectors with the
#'     indices of chromatographic peaks in the object's `chromPeaks` matrix
#'     that should be grouped into features.
#'
#' @param samples optional `integer` to select samples in which the peak
#'     integration should be performed. By default performed in all samples.
#'
#' @param BPPARAM parallel processing settings (see [bpparam()] for details).
#'
#' @param msLevel `integer(1)` defining the MS level in which peak integration
#'     should be performed.
#'
#' @return `XCMSnExp` with the manually added chromatographic peaks or features.
#'
#' @author Johannes Rainer
#'
#' @md
manualChromPeaks <- function(object, chromPeaks = matrix(),
                             samples = seq_along(fileNames(object)),
                             BPPARAM = bpparam(), msLevel = 1L) {
    if (length(msLevel) > 1L)
        stop("Length 'msLevel' is > 1: can only add peaks for one MS level",
             " at a time.")
    if (!inherits(object, "OnDiskMSnExp"))
        stop("'object' has to be either an OnDiskMSnExp or XCMSnExp object")
    if (is(object, "XCMSnExp") && hasFeatures(object))
        object <- dropFeatureDefinitions(object)
    else object <- as(object, "XCMSnExp")
    if (!all(c("mzmin", "mzmax", "rtmin", "rtmax") %in% colnames(chromPeaks)))
        stop("'chromPeaks' lacks one or more of the required columns: 'mzmin',",
             " 'mzmax', 'rtmin' and 'rtmax'")
    if (is.data.frame(chromPeaks)) chromPeaks <- as.matrix(chromPeaks)
    if (!all(samples %in% seq_along(fileNames(object))))
        stop("'samples' out of bounds")
    if (hasChromPeaks(object))
        cn <- colnames(chromPeaks(object))
    else cn <- c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "into",
                 "intb", "maxo", "sn", "sample")
    ## Do integration
    res <- bpmapply(xcms:::.split_by_file2(object)[samples], samples,
                    FUN = function(obj, idx, peakArea, msLevel, cn) {
                        xcms:::.getChromPeakData(obj, peakArea = peakArea,
                                                 sample_idx = idx,
                                                 msLevel = msLevel,
                                                 cn = cn)
                    }, MoreArgs = list(peakArea = chromPeaks,
                                       msLevel = msLevel, cn = cn),
                    BPPARAM = BPPARAM, SIMPLIFY = FALSE, USE.NAMES = FALSE)
    res <- do.call(rbind, res)
    res <- res[!is.na(res[, "sample"]), , drop = FALSE]
    if (nrow(res) == 0)
        return(object)
    newFd <- new("MsFeatureData")
    if (hasChromPeaks(object)) {
        newFd@.xData <- .copy_env(object@msFeatureData)
        object@msFeatureData <- new("MsFeatureData")
        incr <- nrow(chromPeaks(newFd))
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
        chromPeaks(newFd) <- rbind(chromPeaks(newFd), res)
        cpd <- extractROWS(chromPeakData(newFd), rep(1L, nrow(res)))
        cpd[,] <- NA
        cpd$ms_level <- as.integer(msLevel)
        cpd$is_filled <- FALSE
        if (!any(colnames(chromPeakData(newFd)) == "is_filled"))
            chromPeakData(newFd)$is_filled <- FALSE
        chromPeakData(newFd) <- rbind(chromPeakData(newFd), cpd)
        rownames(chromPeakData(newFd)) <- rownames(chromPeaks(newFd))
        lockEnvironment(newFd, bindings = TRUE)
        object@msFeatureData <- newFd
    } else {
        nr <- nrow(res)
        rownames(res) <- sprintf(
            paste0("CP", "%0", ceiling(log10(nr + 1L)), "d"),
            seq_len(nr))
        chromPeaks(newFd) <- res
        cpd <- DataFrame(ms_level = rep(msLevel, nr), is_filled = rep(FALSE, nr))
        rownames(cpd) <- rownames(res)
        chromPeakData(newFd) <- cpd
        lockEnvironment(newFd, bindings = TRUE)
        object@msFeatureData <- newFd
    }
    object
}

#' @rdname manualChromPeaks
manualFeatures <- function(object, peakIdx = list(), msLevel = 1L) {
    if (!length(peakIdx))
        return(object)
    if (length(msLevel) > 1L)
        stop("Length 'msLevel' is > 1: can only add peaks for one MS level",
             " at a time.")
    if (!inherits(object, "XCMSnExp"))
        stop("'object' has to be an XCMSnExp object")
    if (!hasChromPeaks(object))
        stop("No features present. Please run 'findChromPeaks' first.")
    peakIdx <- lapply(peakIdx, as.integer)
    newFd <- new("MsFeatureData")
    newFd@.xData <- xcms:::.copy_env(object@msFeatureData)
    if (hasFeatures(newFd)) {
        fnew <- as.data.frame(featureDefinitions(newFd))[1L, ]
        fnew[] <- NA
        rownames(fnew) <- NULL
    } else {
        cn <- c("mzmed", "mzmin", "mzmax", "rtmed", "rtmin", "rtmax", "npeaks")
        fnew <- as.data.frame(
            matrix(ncol = 7, nrow = 1, dimnames = list(character(), cn)))
    }
    res <- lapply(peakIdx, function(z) {
        cp <- chromPeaks(newFd)[z, , drop = FALSE]
        if (any(is.na(cp[, "sample"])))
            stop("Some of the provided indices are out of bounds. 'peakIdx' ",
                 "needs to be a list of valid indices in the 'chromPeaks' ",
                 "matrix.", call. = FALSE)
        newf <- fnew
        newf$mzmed <- median(cp[, "mz"])
        newf$mzmin <- min(cp[, "mz"])
        newf$mzmax <- max(cp[, "mz"])
        newf$rtmed <- median(cp[, "rt"])
        newf$rtmin <- min(cp[, "rt"])
        newf$rtmax <- max(cp[, "rt"])
        newf$npeaks <- length(z)
        newf
    })
    res <- DataFrame(do.call(rbind, res))
    res$peakidx <- peakIdx
    res$ms_level <- msLevel
    ## Define feature IDs
    if (hasFeatures(newFd))
        max_id <- max(
            as.integer(sub("FT", "", rownames(featureDefinitions(newFd)))))
    else max_id <- 0
    rownames(res) <- sprintf(
        paste0("FT", "%0", ceiling(log10(max_id + nrow(res) + 1L)), "d"),
        (max_id + 1L):(max_id + nrow(res)))
    suppressWarnings(
        featureDefinitions(newFd) <- rbind(featureDefinitions(newFd), res)
    )
    lockEnvironment(newFd, bindings = TRUE)
    object@msFeatureData <- newFd
    object
}
