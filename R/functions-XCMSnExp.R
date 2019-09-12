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
        xs@groups <- S4Vectors::as.matrix(fgs[, -ncol(fgs)])
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
    xs@phenoData <- pData(from)
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


## #' @description Integrates the intensities for chromatograpic peak(s). This is
## #'     supposed to be called by the fillChromPeaks method.
## #'
## #' @note This reads the full data first and does the subsetting later in R.
## #'
## #' @param object An \code{XCMSnExp} object representing a single sample.
## #'
## #' @param peakArea A \code{matrix} with the peak definition, i.e. \code{"rtmin"},
## #'     \code{"rtmax"}, \code{"mzmin"} and \code{"mzmax"}.
## #'
## #' @noRd
## .getPeakInt2 <- function(object, peakArea) {
##     if (length(fileNames(object)) != 1)
##         stop("'object' should be an XCMSnExp for a single file!")
##     res <- numeric(nrow(peakArea))
##     spctr <- spectra(object, BPPARAM = SerialParam())
##     mzs <- lapply(spctr, mz)
##     valsPerSpect <- lengths(mzs)
##     ints <- unlist(lapply(spctr, intensity), use.names = FALSE)
##     rm(spctr)
##     mzs <- unlist(mzs, use.names = FALSE)
##     rtim <- rtime(object)
##     for (i in 1:length(res)) {
##         rtr <- peakArea[i, c("rtmin", "rtmax")]
##         mtx <- .rawMat(mz = mzs, int = ints, scantime = rtim,
##                        valsPerSpect = valsPerSpect, rtrange = rtr,
##                        mzrange = peakArea[i, c("mzmin", "mzmax")])
##         if (length(mtx)) {
##             if (!all(is.na(mtx[, 3]))) {
##                 ## How to calculate the area: (1)sum of all intensities / (2)by
##                 ## the number of data points (REAL ones, considering also NAs)
##                 ## and multiplied with the (3)rt width.
##                 ## (1) sum(mtx[, 3], na.rm = TRUE)
##                 ## (2) sum(rtim >= rtr[1] & rtim <= rtr[2]) - 1 ; if we used
##                 ## nrow(mtx) here, which would correspond to the non-NA
##                 ## intensities within the rt range we don't get the same results
##                 ## as e.g. centWave.
##                 ## (3) rtr[2] - rtr[1]
##                 res[i] <- sum(mtx[, 3], na.rm = TRUE) *
##                     ((rtr[2] - rtr[1]) /
##                      (sum(rtim >= rtr[1] & rtim <= rtr[2]) - 1))
##             } else {
##                 res[i] <- NA_real_
##             }
##         } else {
##             res[i] <- NA_real_
##         }
##     }
##     return(unname(res))
## }

## #' @description Integrates the intensities for chromatograpic peak(s). This is
## #'     supposed to be called by the fillChromPeaks method.
## #'
## #' @note This reads the full data first and does the subsetting later in R. This
## #'     function uses the C getEIC function.
## #'
## #' @param object An \code{XCMSnExp} object representing a single sample.
## #'
## #' @param peakArea A \code{matrix} with the peak definition, i.e. \code{"rtmin"},
## #'     \code{"rtmax"}, \code{"mzmin"} and \code{"mzmax"}.
## #'
## #' @noRd
## .getPeakInt3 <- function(object, peakArea) {
##     if (length(fileNames(object)) != 1)
##         stop("'object' should be an XCMSnExp for a single file!")
##     if (nrow(peakArea) == 0) {
##         return(numeric())
##     }
##     res <- matrix(ncol = 4, nrow = nrow(peakArea))
##     res <- numeric(nrow(peakArea))
##     spctr <- spectra(object, BPPARAM = SerialParam())
##     mzs <- lapply(spctr, mz)
##     valsPerSpect <- lengths(mzs)
##     scanindex <- valueCount2ScanIndex(valsPerSpect) ## Index vector for C calls
##     ints <- unlist(lapply(spctr, intensity), use.names = FALSE)
##     rm(spctr)
##     mzs <- unlist(mzs, use.names = FALSE)
##     rtim <- rtime(object)
##     for (i in 1:length(res)) {
##         rtr <- peakArea[i, c("rtmin", "rtmax")]
##         sr <- c(min(which(rtim >= rtr[1])), max(which(rtim <= rtr[2])))
##         eic <- .Call("getEIC", mzs, ints, scanindex,
##                      as.double(peakArea[i, c("mzmin", "mzmax")]),
##                      as.integer(sr), as.integer(length(scanindex)),
##                      PACKAGE = "xcms")
##         if (length(eic$intensity)) {
##             ## How to calculate the area: (1)sum of all intensities / (2)by
##             ## the number of data points (REAL ones, considering also NAs)
##             ## and multiplied with the (3)rt width.
##             if (!all(is.na(eic$intensity)) && !all(eic$intensity == 0)) {
##                 res[i] <- sum(eic$intensity, na.rm = TRUE) *
##                     ((rtr[2] - rtr[1]) / (length(eic$intensity) - 1))
##             } else {
##                 res[i] <- NA_real_
##             }
##         } else {
##             res[i] <- NA_real_
##         }
##     }
##     return(unname(res))
## }


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
                              cn = c("mz", "rt", "into", "maxo", "sample")) {
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
    spctr <- spectra(object, BPPARAM = SerialParam())
    mzs <- lapply(spctr, mz)
    valsPerSpect <- lengths(mzs)
    ints <- unlist(lapply(spctr, intensity), use.names = FALSE)
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
    return(res)
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
        if (length(mtx)) {
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
                                                   "sample")) {
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
            .update_feature_definitions(featureDefinitions(object),
                                        rownames(chromPeaks(object)),
                                        rownames(chromPeaks(object,
                                                            msLevel = msLevel)))),
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
#' ## Below we perform first a peak detection (using the matchedFilter
#' ## method) on some of the test files from the faahKO package followed by
#' ## a peak grouping and retention time adjustment using the "peak groups"
#' ## method
#' library(faahKO)
#' library(xcms)
#' fls <- dir(system.file("cdf/KO", package = "faahKO"), recursive = TRUE,
#'            full.names = TRUE)
#'
#' ## Reading 2 of the KO samples
#' raw_data <- readMSData(fls[1:2], mode = "onDisk")
#'
#' ## Perform the peak detection using the matchedFilter method.
#' mfp <- MatchedFilterParam(snthresh = 20, binSize = 1)
#' res <- findChromPeaks(raw_data, param = mfp)
#'
#' ## Performing the peak grouping using the "peak density" method.
#' p <- PeakDensityParam(sampleGroups = c(1, 1))
#' res <- groupChromPeaks(res, param = p)
#'
#' ## Perform the retention time adjustment using peak groups found in both
#' ## files.
#' fgp <- PeakGroupsParam(minFraction = 1)
#' res <- adjustRtime(res, param = fgp)
#'
#' ## Visualize the impact of the alignment. We show both versions of the plot,
#' ## with the raw retention times on the x-axis (top) and with the adjusted
#' ## retention times (bottom).
#' par(mfrow = c(2, 1))
#' plotAdjustedRtime(res, adjusted = FALSE)
#' grid()
#' plotAdjustedRtime(res)
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
#' \code{\link{Chromatograms}} object.
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
#' ## Read some files from the faahKO package.
#' library(xcms)
#' library(faahKO)
#' faahko_3_files <- c(system.file('cdf/KO/ko16.CDF', package = "faahKO"),
#'                     system.file('cdf/KO/ko18.CDF', package = "faahKO"))
#'
#' od <- readMSData(faahko_3_files, mode = "onDisk")
#'
#' ## Peak detection using the 'matchedFilter' method. Note that we are using a
#' ## larger binSize to reduce the runtime of the example.
#' xod <- findChromPeaks(od, param = MatchedFilterParam(binSize = 0.3, snthresh = 20))
#'
#' ## Extract the ion chromatogram for one chromatographic peak in the data.
#' chrs <- chromatogram(xod, rt = c(2700, 2900), mz = 335)
#'
#' plot(chrs)
#'
#' ## Extract chromatographic peaks for the mz/rt range (if any).
#' chromPeaks(xod, rt = c(2700, 2900), mz = 335)
#'
#' ## Highlight the chromatographic peaks in the area
#' ## Show the peak definition with a rectangle
#' highlightChromPeaks(xod, rt = c(2700, 2900), mz = 335)
#'
#' ## Color the actual peak
#' highlightChromPeaks(xod, rt = c(2700, 2900), mz = 335,
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
#' ## Perform peak detection on two files from the faahKO package.
#' library(xcms)
#' library(faahKO)
#' faahko_file <- c(system.file('cdf/KO/ko16.CDF', package = "faahKO"),
#'                  system.file('cdf/KO/ko18.CDF', package = "faahKO"))
#'
#' od <- readMSData(faahko_file, mode = "onDisk")
#'
#' ## Peak detection using the 'matchedFilter' method. Note that we are using a
#' ## larger binSize to reduce the runtime of the example.
#' xod <- findChromPeaks(od, param = MatchedFilterParam(binSize = 0.3, snthresh = 20))
#'
#' ## plotChromPeakImage: plot an image for the identified peaks per file
#' plotChromPeakImage(xod)
#'
#' ## Show all detected chromatographic peaks from the first file
#' plotChromPeaks(xod)
#'
#' ## Plot all detected peaks from the second file and restrict the plot to a
#' ## mz-rt slice
#' plotChromPeaks(xod, file = 2, xlim = c(3500, 3600), ylim = c(400, 600))
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
        rts <- split(pks[, "rt"], pks[, "sample"])
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
#' ## Load test data
#' files <- c(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
#'     system.file('cdf/KO/ko16.CDF', package = "faahKO"),
#'     system.file('cdf/KO/ko18.CDF', package = "faahKO"))
#'
#' od <- readMSData(files, mode = "onDisk")
#'
#' ## Apply obiwarp retention time adjustment. We have to convert the
#' ## OnDiskMSnExp first to an XCMSnExp
#' xod <- as(od, "XCMSnExp")
#' xod <- adjustRtime(xod, param = ObiwarpParam())
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
#' plot(split(rtime(od), fromFile(od))[[1]] -
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
    newFd <- dropAdjustedRtime(newFd)
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

#' @description
#'
#' \code{filterFeatureDefinitions} allows to subset the feature definitions of
#' an \code{XCMSnExp} object. Which feature definitions should be kept can be
#' specified with the \code{features} argument that can be a \code{logical},
#' \code{integer} or \code{character} vector. The function returns the
#' \code{XCMSnExp} with the reduced \code{featureDefinitions} data frame.
#'
#' @param features For \code{filterFeatureDefinitions}: either a \code{integer}
#'     specifying the indices of the features (rows) to keep, a \code{logical}
#'     with a length matching the number of rows of \code{featureDefinitions}
#'     or a \code{character} with the feature (row) names.
#'
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
#' ## Load 2 test files.
#' data <- readMSData(c(system.file("cdf/KO/ko15.CDF", package = "faahKO"),
#'                      system.file("cdf/KO/ko16.CDF", package = "faahKO")),
#'                    mode = "onDisk")
#'
#' ## Perform peak detection; parameters set to reduce processing speed
#' data <- findChromPeaks(data, CentWaveParam(noise = 10000, snthresh = 40))
#'
#' ## Correspondence analysis
#' data <- groupChromPeaks(data, param = PeakDensityParam(sampleGroups = c(1, 1)))
#'
#' ## Identify overlapping features
#' overlappingFeatures(data)
#'
#' ## Identify features that are separated on retention time by less than
#' ## 2 minutes
#' overlappingFeatures(data, expandRt = 60)
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
#' software (or the `MetaboAnalystR` R package.
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
#'     format *M<m/z>T<rt>* (`groupnames = TRUE`). See help of the [groupnames]
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
ms2_spectra_for_all_peaks <- function(x, expandRt = 0, expandMz = 0,
                                  ppm = 0, method = c("all",
                                                      "closest_rt",
                                                      "closest_mz",
                                                      "signal"),
                                  skipFilled = FALSE, subset = NULL,
                                  BPPARAM = bpparam()) {
    method <- match.arg(method)
    if (hasAdjustedRtime(x))
        x <- applyAdjustedRtime(x)
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
        pks[chromPeakData(x)$is_filled] <- NA
    x <- filterMsLevel(as(x, "OnDiskMSnExp"), 2L)
    ## Split data per file
    file_factor <- factor(pks[, "sample"])
    peak_ids <- rownames(pks)
    pks <- split.data.frame(pks, f = file_factor)
    x <- lapply(as.integer(levels(file_factor)), filterFile, object = x)
    res <- bpmapply(ms2_spectra_for_peaks_from_file, x, pks,
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
ms2_spectra_for_peaks_from_file <- function(x, pks, method = c("all",
                                                               "closest_rt",
                                                               "closest_mz",
                                                               "signal")) {
    if (nrow(pks) == 0)
        return(list())
    method <- match.arg(method)
    fromFile <- as.integer(pks[1, "sample"])
    sps <- spectra(x)
    pmz <- precursorMz(x)
    rtm <- rtime(x)
    res <- vector(mode = "list", nrow(pks))
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

#' @title Extract (MS2) spectra associated with chromatographic peaks
#'
#' @description
#'
#' Extract (MS2) spectra from an [XCMSnExp] object that represent ions within
#' the rt and m/z range of each chromatographic peak (in the same file
#' /sample in which the peak was detected). All MS2 spectra are returned for
#' chromatographic peak `i` for which the precursor m/z is
#' `>= chromPeaks(x)[i, "mzmin"]` and `<= chromPeaks(x)[i, "mzmax"]` and the
#' retention time is `>= chromPeaks(x)[i, "rtmin"]` and
#' `<= chromPeaks(x)[i, "rtmax"]`.
#'
#' @param x [XCMSnExp] object with identified chromatographic peaks.
#'
#' @param msLevel `integer(1)` defining whether MS1 or MS2 spectra should be
#'     returned. Currently only `msLevel = 2` is supported.
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
#' @param method `character(1)` specifying which MS2 spectra should be included.
#'     Defaults to `"all"` in which all MS2 spectra within the rt and m/z range
#'     of a chromatographic peak are returned. `"closest_rt"` returns the one
#'     MS2 spectrum with the retention time closest to the chromatographic
#'     peak's apex rt. `"closest_mz"` returns the MS2 spectrum with the
#'     precursor m/z closest to the chromatographic peak's m/z. `"signal"`
#'     returns the MS2 spectrum which total signal is closest to the
#'     chromatographic peak's maximal signal (`"maxo"`).
#'
#' @param skipFilled `logical(1)` whether no spectra for filled-in peaks should
#'     be reported.
#'
#' @param return.type `character(1)` defining whether the result should be a
#'     [Spectra] object or a simple `list`. See below for more information.
#'
#' @return
#'
#' Which object is returned depends on the value of `return.type`:
#'
#' - For `return.type = "Spectra"`: a [Spectra] object with elements being
#'   [Spectrum-class] objects. The result objects contains all spectra
#'   for all peaks. Metadata column `"peak_id"` provides the ID of the
#'   respective peak (i.e. its rowname in [chromPeaks()]).
#' - If `return.type = "list"`: `list` of `list`s that are either of length
#'   0 or contain [Spectrum2-class] object(s) within the m/z-rt range. The
#'   length of the list matches the number of peaks.
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
#' ## Perform MS1 peak detection
#' dda <- findChromPeaks(dda, CentWaveParam(peakwidth = c(5, 15)))
#' ms2_sps <- chromPeakSpectra(dda)
#' ms2_sps
#'
#' ## Metadata column `peak_id` contains the ID of the chromatographic peak
#' ## of the MS2 spectrum
chromPeakSpectra <- function(x, msLevel = 2L, expandRt = 0, expandMz = 0,
                             ppm = 0, method = c("all", "closest_rt",
                                                 "closest_mz", "signal"),
                             skipFilled = FALSE,
                             return.type = c("Spectra", "list")) {
    method <- match.arg(method)
    return.type <- match.arg(return.type)
    if (!is(x, "XCMSnExp"))
        stop("'x' is supposed to be an 'XCMSnExp' object.")
    if (!hasChromPeaks(x))
        stop("No chromatographic peaks present. Please run 'findChromPeaks' ",
             "first")
    if (msLevel != 2 || (msLevel == 2 & !any(msLevel(x) == 2))) {
        res <- vector(mode = "list", length = nrow(chromPeaks(x)))
        names(res) <- rownames(chromPeaks(x))
        if (msLevel != 2)
            warning("msLevel = 1 is currently not supported.")
        if (msLevel == 2 & !any(msLevel(x) == 2))
            warning("No MS2 spectra available in 'x'.")
    } else {
        res <- ms2_spectra_for_all_peaks(x, expandRt = expandRt,
                                         expandMz = expandMz,
                                         ppm = ppm, method = method,
                                         skipFilled = skipFilled)
    }
    if (return.type == "Spectra") {
        pids <- rep(names(res), lengths(res))
        res <- res[lengths(res) > 0]
        if (length(res))
            res <- unlist(res)
        res <- Spectra(res, elementMetadata = DataFrame(peak_id = pids))
    }
    res
}

#' For information and details see featureSpectra
#'
#' @noRd
ms2_spectra_for_features <- function(x, expandRt = 0, expandMz = 0, ppm = 0,
                                      skipFilled = FALSE, ...) {
    idxs <- featureDefinitions(x)$peakidx
    sp_pks <- ms2_spectra_for_all_peaks(x, expandRt = expandRt,
                                        expandMz = expandMz, ppm = ppm,
                                        skipFilled = skipFilled,
                                        subset = unique(unlist(idxs)), ...)
    res <- lapply(idxs, function(z) unlist(sp_pks[z]))
    names(res) <- rownames(featureDefinitions(x))
    res
}

#' @title Extract (MS2) spectra associated with features
#'
#' @description
#'
#' Return (MS2) spectra for all chromatographic peaks associated with the
#' features in `x` (defined by [featureDefinitions()]).
#'
#' @details
#'
#' The function identifies all MS2 spectra with their precursor m/z within the
#' m/z range of a chromatographic peak (i.e. `>=` mzmin and `<=` mzmax) of a
#' feature and their retention time within the rt range of the same peak
#' (`>=` rtmin and `<=` rtmax).
#'
#' The optional parameter `method` allows to ensure that for each
#' chromatographic peak in one sample only one MS2 spectrum is returned.
#' See [chromPeakSpectra()] for more details.
#'
#' @param x [XCMSnExp] object with feature defitions available.
#'
#' @inheritParams chromPeakSpectra
#'
#' @param ... additional arguments to be passed along to [chromPeakSpectra()],
#'     such as `method`.
#'
#' @return
#'
#' Which object is returned depends on the value of `return.type`:
#'
#' - For `return.type = "Spectra"` (the default): a [Spectra] object with data
#'   only for features for which a [Spectrum2-class] was found. The
#'   ID of the feature and of the chromatographic peak to which the
#'   spectrum is associated are provided with the `"feature_id"` and
#'   `"peak_id"` metadata columns.
#' - For `return.type = "list"`: a `list, same length than there are
#'   features in `x`, each element being a `list` of [Spectrum2-class]
#'   objects with the MS2 spectra that potentially represent ions measured
#'   by each chromatographic peak associated with the feature, or `NULL`
#'   if none are present.
#'
#' @author Johannes Rainer
#'
#' @md
featureSpectra <- function(x, msLevel = 2, expandRt = 0, expandMz = 0,
                           ppm = 0, skipFilled = FALSE,
                           return.type = c("Spectra", "list"), ...) {
    if (!is(x, "XCMSnExp"))
        stop("'x' is supposed to be an 'XCMSnExp' object.")
    return.type <- match.arg(return.type)
    if (!hasFeatures(x))
        stop("No feature definitions present. Please run 'groupChromPeaks' ",
             "first.")
    if (msLevel != 2 || (msLevel == 2 & !any(msLevel(x) == 2))) {
        res <- vector(mode = "list", length = nrow(featureDefinitions(x)))
        names(res) <- rownames(featureDefinitions(x))
        if (msLevel != 2)
            warning("msLevel = 1 is currently not supported.")
        if (msLevel == 2 & !any(msLevel(x) == 2))
            warning("No MS2 spectra available in 'x'.")
    } else {
        res <- ms2_spectra_for_features(x, expandRt = expandRt,
                                        expandMz = expandMz,
                                        ppm = ppm, skipFilled = skipFilled,
                                        ...)
    }
    if (return.type == "Spectra") {
        fids <- rep(names(res), lengths(res))
        res <- res[lengths(res) > 0]
        if (length(res)) {
            res <- unlist(res)
            pids <- vapply(strsplit(names(res), ".", TRUE),
                           `[`, character(1), 2)
        } else {
            pids <- character()
        }
        res <- Spectra(res, elementMetadata = DataFrame(feature_id = fids,
                                                        peak_id = pids))
    }
    res
}

#' @title Extract ion chromatograms for each feature
#'
#' @description
#'
#' Extract ion chromatograms for features in an [XCMSnExp-class] object. The
#' function returns for each feature its extracted ion chromatogram and all
#' associated peaks with it.
#'
#' By default only chromatographic peaks associated with a feature are included
#' for an extracted ion chromatogram. Setting `include = "all"` (instead of
#' the default `include = "feature_only"`) will return all chromatographic peaks
#' identified in the m/z - rt data slice of a feature (and eventually also other
#' features within that region).
#'
#' @param x `XCMSnExp` object with grouped chromatographic peaks.
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
#' @param include `character(1)` defining which chromatographic peaks and
#'     feature definitions should be included in the returned
#'     [XChromatograms()]. See description above for details.
#'
#' @param filled `logical(1)` whether filled-in peaks should be included in
#'     the result object. The default is `filled = FALSE`, i.e. only detected
#'     peaks are reported.
#'
#' @param ... optional arguments to be passed along to the [chromatogram()]
#'     function.
#'
#' @return [XChromatograms()] object.
#'
#' @md
#'
#' @author Johannes Rainer
#'
#' @examples
#'
#' library(xcms)
#' library(faahKO)
#' faahko_3_files <- c(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
#'                     system.file('cdf/KO/ko16.CDF', package = "faahKO"),
#'                     system.file('cdf/KO/ko18.CDF', package = "faahKO"))
#'
#' ## Do a simple and fast preprocessing of the test data
#' od <- readMSData(faahko_3_files, mode = "onDisk")
#' od <- findChromPeaks(od, param = CentWaveParam(peakwidth = c(30, 80),
#'     noise = 1000))
#' od <- adjustRtime(od, param = ObiwarpParam(binSize = 0.6))
#' od <- groupChromPeaks(od,
#'     param = PeakDensityParam(minFraction = 0.8, sampleGroups = rep(1, 3)))
#'
#' ## Extract ion chromatograms for each feature
#' chrs <- featureChromatograms(od)
#'
#' ## Plot the XIC for the first feature using different colors for each file
#' par(mfrow = c(1, 2))
#' plot(chrs[1, ], col = c("red", "green", "blue"))
featureChromatograms <- function(x, expandRt = 0, aggregationFun = "max",
                                 features, include = c("feature_only", "all"),
                                 filled = FALSE, ...) {
    include <- match.arg(include)
    if (!hasFeatures(x))
        stop("No feature definitions present. Please run first 'groupChromPeaks'")
    if (!missing(features)) {
        if (is.logical(features)) {
            if (length(features) != nrow(featureDefinitions(x)))
                stop("If 'features' is a logical vector, its length has to ",
                     "be equal to the number of features.")
            features <- which(features)
        }
        if (is.character(features)) {
            features <- match(features, rownames(featureDefinitions(x)))
            if (any(is.na(features)))
                stop("Some of the feature IDs specified with 'features' ",
                     "can not be found.")
        }
        features <- as.integer(features)
        if (any(features < 1) || any(features > nrow(featureDefinitions(x))))
            stop("'features' are out of range")
    } else features <- seq_len(nrow(featureDefinitions(x)))
    if (!length(features))
        return(Chromatograms(ncol = length(fileNames(x)), nrow = 0))
    pks <- chromPeaks(x)
    if (length(unique(rownames(pks))) != nrow(pks)) {
        rownames(pks) <- .featureIDs(nrow(pks), "CP")
        rownames(chromPeaks(x)) <- rownames(pks)
    }
    pk_idx <- featureDefinitions(x)$peakidx[features]
    mat <- do.call(rbind, lapply(pk_idx, function(z) {
        pks_current <- pks[z, , drop = FALSE]
        c(range(pks_current[, c("rtmin", "rtmax")]),
          range(pks_current[, c("mzmin", "mzmax")]))
    }))
    mat[, 1] <- mat[, 1] - expandRt
    mat[, 2] <- mat[, 2] + expandRt
    colnames(mat) <- c("rtmin", "rtmax", "mzmin", "mzmax")
    chrs <- chromatogram(x, rt = mat[, 1:2], mz = mat[, 3:4],
                         aggregationFun = aggregationFun, filled = filled, ...)
    if (include == "feature_only") {
        pk_ids <- lapply(pk_idx, function(z) rownames(pks)[z])
        fts_all <- featureDefinitions(chrs)
        pks_all <- chromPeaks(chrs)
        chrs@featureDefinitions <- fts_all[integer(), ]
        for (i in 1:nrow(chrs)) {
            for (j in 1:ncol(chrs)) {
                cur_pks <- chrs@.Data[i, j][[1]]@chromPeaks
                if (nrow(cur_pks)) {
                    keep <- rownames(cur_pks) %in% pk_ids[[i]]
                    chrs@.Data[i, j][[1]]@chromPeaks <- cur_pks[keep, ,
                                                                drop = FALSE]
                    chrs@.Data[i, j][[1]]@chromPeakData <-
                        chrs@.Data[i, j][[1]]@chromPeakData[keep, , drop = FALSE]
                }
            }
        }
        chrs@featureDefinitions <- .subset_features_on_chrom_peaks(
            fts_all, pks_all, chromPeaks(chrs))
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
            ret$isolationWindow <- fData(z)$isolationWindow[1]
            ret$isolationWindowTargetMZ <- target_mz
            ret$isolationWindowLowerMz <-
                target_mz - fData(z)$isolationWindowLowerOffset[1]
            ret$isolationWindowUpperMz <-
                target_mz + fData(z)$isolationWindowUpperOffset[1]
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
#' Reconstructs MS2 spectra for each MS1 chromatographic peak (if possible) for
#' data independent acquisition (DIA) data (such as SWATH).
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
#' The resulting `Spectra` object provides also the peak IDs of the MS2
#' chromatographic peaks for each spectrum as well as their correlation value.
#'
#' @param object `XCMSnExp` with identified chromatographic peaks.
#'
#' @param expandRt `numeric(1)` allowing to expand the retention time range
#'     for extracted ion chromatograms by a constant value (for the peak
#'     shape correlation).
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
#'     spectra's peaks.
#'
#' @param peakId optional `character` vector with peak IDs (i.e. rownames of
#'     `chromPeaks`) of MS1 peaks for which MS2 spectra should be reconstructed.
#'     By default they are reconstructed for all MS1 chromatographic peaks.
#'
#' @param BPPARAM parallel processing setup. See [bpparam()] for more
#'     information.
#'
#' @return [Spectra()] with the reconstructed MS2 spectra for all MS1 peaks
#'     in `object`. Contains empty [Spectrum2-class] objects for MS1 peaks for
#'     which reconstruction was not possible (either no MS2 signal was recorded
#'     or the correlation of the MS2 chromatographic peaks with the MS1
#'     chromatographic peak was below threshold `minCor`. `Spectra` metadata
#'     columns `"ms2_peak_id"` and `"ms2_peak_cor"` (of type [CharacterList()]
#'     and [NumericList()] with length equal to the number of peaks per
#'     reconstructed MS2 spectrum) providing the IDs and the correlation of the
#'     MS2 chromatographic peaks from which the MS2 spectrum was reconstructed.
#'
#' @author Johannes Rainer, Micheal Witting
#'
#' @md
#'
#' @seealso [findChromPeaksIsolationWindow()] for the function to perform MS2
#'     peak detection in DIA isolation windows and for examples.
reconstructChromPeakSpectra <- function(object, expandRt = 1, diffRt = 2,
                                        minCor = 0.8, intensity = "maxo",
                                        peakId = rownames(
                                            chromPeaks(object, msLevel = 1L)),
                                        BPPARAM = bpparam()) {
    if (!inherits(object, "XCMSnExp") || !hasChromPeaks(object))
        stop("'object' should be an 'XCMSnExp' object with identified ",
             "chromatographic peaks")
    if (!is.character(peakId))
        stop("'peakId' has to be of type character")
    n_peak_id <- length(peakId)
    peakId <- intersect(peakId, rownames(chromPeaks(object, msLevel = 1L)))
    if (!length(peakId))
        stop("None of the provided 'peakId' matches IDs of MS1 ",
             "chromatographic peaks")
    if (length(peakId) < n_peak_id)
        warning("Only ", length(peakId), " of the provided",
                " identifiers match IDs of MS1 chromatographic peaks")
    ## Drop featureDefinitions if present
    if (hasFeatures(object))
        suppressMessages(object <- dropFeatureDefinitions(
                             object, keepAdjustedRtime = TRUE))
    object <- selectFeatureData(object,
                                fcol = c(MSnbase:::.MSnExpReqFvarLabels,
                                         "centroided",
                                         "polarity",
                                         "isolationWindow",
                                         "isolationWindowTargetMZ",
                                         "isolationWindowLowerOffset",
                                         "isolationWindowUpperOffset"))
    sps <- bplapply(
        lapply(seq_len(length(fileNames(object))), filterFile, object = object,
               keepAdjustedRtime = TRUE),
        FUN = function(x, files, expandRt, diffRt, minCor, col, pkId) {
            .reconstruct_ms2_for_peaks_file(
                x, expandRt = expandRt, diffRt = diffRt,
                minCor = minCor, fromFile = match(fileNames(x), files),
                column = col, peakId = pkId)
        },
        files = fileNames(object), expandRt = expandRt, diffRt = diffRt,
        minCor = minCor, col = intensity, pkId = peakId, BPPARAM = BPPARAM)
    do.call(c, sps)
}
