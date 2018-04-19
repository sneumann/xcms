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
    if (any(chromPeaks(from)[, "is_filled"] == 1)) {
        fld <- which(chromPeaks(from)[, "is_filled"] == 1)
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
    res[, c("mzmin", "mzmax")] <-
        peakArea[, c("mzmin", "mzmax")]
    ## Load the data
    message("Requesting ", nrow(res), " missing peaks from ",
             basename(fileNames(object)), " ... ", appendLF = FALSE)
    spctr <- spectra(object, BPPARAM = SerialParam())
    mzs <- lapply(spctr, mz)
    valsPerSpect <- lengths(mzs)
    ints <- unlist(lapply(spctr, intensity), use.names = FALSE)
    rm(spctr)
    mzs <- unlist(mzs, use.names = FALSE)
    rtim <- rtime(object)
    rtim_range <- range(rtim)
    for (i in 1:nrow(res)) {
        rtr <- peakArea[i, c("rtmin", "rtmax")]
        ## Ensure that the rt region is within the rtrange of the data.
        rtr[1] <- max(rtr[1], rtim_range[1])
        rtr[2] <- min(rtr[2], rtim_range[2])
        mtx <- .rawMat(mz = mzs, int = ints, scantime = rtim,
                       valsPerSpect = valsPerSpect, rtrange = rtr,
                       mzrange = peakArea[i, c("mzmin", "mzmax")])
        if (length(mtx)) {
            if (!all(is.na(mtx[, 3]))) {
                ## How to calculate the area: (1)sum of all intensities / (2)by
                ## the number of data points (REAL ones, considering also NAs)
                ## and multiplied with the (3)rt width.
                ## (1) sum(mtx[, 3], na.rm = TRUE)
                ## (2) sum(rtim >= rtr[1] & rtim <= rtr[2]) - 1 ; if we used
                ## nrow(mtx) here, which would correspond to the non-NA
                ## intensities within the rt range we don't get the same results
                ## as e.g. centWave.
                ## (3) rtr[2] - rtr[1]
                res[i, "into"] <- sum(mtx[, 3], na.rm = TRUE) *
                    ((rtr[2] - rtr[1]) /
                     (sum(rtim >= rtr[1] & rtim <= rtr[2]) - 1))
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
    message("Requesting ", nrow(res), " missing peaks from ",
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
    message("Reguesting ", nrow(res), " missing peaks from ",
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
    message("Requesting ", nrow(res), " missing peaks from ",
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
    if (hasChromPeaks(object))
        if (any(colnames(chromPeaks(object)) == "is_filled"))
            return(any(chromPeaks(object)[, "is_filled"] == 1))
    FALSE
}

#' @description
#'
#' Simple helper function to extract the peakidx column from the
#' featureDefinitions DataFrame. The function ensures that the names of the
#' returned list correspond to the rownames of the DataFrame
#' 
#' @noRd
.peakIndex <- function(object) {
    if (!hasFeatures(object))
        stop("No feature definitions present. Please run groupChromPeaks first.")
    idxs <- featureDefinitions(object)$peakidx
    names(idxs) <- rownames(featureDefinitions(object))
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
adjustRtimePeakGroups <- function(object, param = PeakGroupsParam()) {
    if (!is(object, "XCMSnExp"))
        stop("'object' has to be an 'XCMSnExp' object.")
    if (!hasFeatures(object))
        stop("No features present. Please run 'groupChromPeaks' first.")
    if (hasAdjustedRtime(object))
        warning("Alignment/retention time correction was already performed, ",
                "returning a matrix with adjusted retention times.")
    nSamples <- length(fileNames(object))
    pkGrp <- .getPeakGroupsRtMatrix(
        peaks = chromPeaks(object),
        peakIndex = .peakIndex(object),
        nSamples = nSamples,
        missingSample = nSamples - (nSamples * minFraction(param)),
        extraPeaks = extraPeaks(param)
    )
    colnames(pkGrp) <- basename(fileNames(object))
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
plotAdjustedRtime <- function(object, col = "#00000080", lty = 1, type = "l",
                              adjustedRtime = TRUE,
                              xlab = ifelse(adjustedRtime,
                                            yes = expression(rt[adj]),
                                            no = expression(rt[raw])),
                              ylab = expression(rt[adj]-rt[raw]),
                              peakGroupsCol = "#00000060",
                              peakGroupsPch = 16,
                              peakGroupsLty = 3, ...) {
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
    ## Initialize plot.
    plot(3, 3, pch = NA, xlim = range(xRt, na.rm = TRUE),
         ylim = range(diffRt, na.rm = TRUE), xlab = xlab, ylab = ylab, ...)
    ## Plot all.
    for (i in 1:length(diffRt))
        points(x = xRt[[i]], y = diffRt[[i]], col = col[i], lty = lty[i],
               type = type)
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

#' @title Plot chromatographic peak density along the retention time axis
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
#' @examples
#'
#' ## Below we perform first a peak detection (using the centWave
#' ## method) on some of the test files from the faahKO package.
#' library(faahKO)
#' library(xcms)
#' fls <- dir(system.file("cdf/KO", package = "faahKO"), recursive = TRUE,
#'            full.names = TRUE)
#' 
#' ## Reading 2 of the KO samples
#' raw_data <- readMSData(fls[1:2], mode = "onDisk")
#'
#' ## Perform the peak detection using the centWave method (settings are tuned
#' ## to speed up example execution)
#' res <- findChromPeaks(raw_data, param = CentWaveParam(noise = 3000, snthresh = 40))
#'
#' ## Align the samples using obiwarp
#' res <- adjustRtime(res, param = ObiwarpParam())
#'
#' ## Plot the chromatographic peak density for a specific mz range to evaluate
#' ## different peak density correspondence settings.
#' mzr <- c(305.05, 305.15)
#'
#' plotChromPeakDensity(res, mz = mzr, pch = 16,
#'     param = PeakDensityParam(sampleGroups = rep(1, length(fileNames(res)))))
#'
#' ## Use a larger bandwidth
#' plotChromPeakDensity(res, mz = mzr, param = PeakDensityParam(bw = 60,
#'     sampleGroups = rep(1, length(fileNames(res)))), pch = 16)
#' ## Neighboring peaks are now fused into one.
#'
#' ## Require the chromatographic peak to be present in all samples of a group
#' plotChromPeakDensity(res, mz = mzr, pch = 16,
#'     param = PeakDensityParam(minFraction = 1,
#'     sampleGroups = rep(1, length(fileNames(res)))))
plotChromPeakDensity <- function(object, mz, rt, param, simulate = TRUE,
                                 col = "#00000080", xlab = "retention time",
                                 ylab = "sample", xlim = range(rt),
                                 main = NULL, ...) {
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
    pks <- chromPeaks(object, mz = mz, rt = rt)
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
            fts <- featureDefinitions(object, mz = mz, rt = rt)
            rect(xleft = fts$rtmin, xright = fts$rtmax,
                 ybottom = rep(yl[1], nrow(fts)), ytop = rep(yl[2], nrow(fts)),
                 border = "#00000040", col = "#00000020")
            abline(v = fts$rtmed, col = "#00000040", lty = 2)
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
#' @param border colors to be used to color the border of the rectangles. Has to
#'     be equal to the number of samples in \code{x}.
#' 
#' @param lwd \code{numeric(1)} defining the width of the line/border.
#'
#' @param col For \code{highlightChromPeaks}: color to be used to fill the
#'     rectangle.
#'
#' @param type the plotting type. See \code{\link{plot}} in base grapics for
#'     more details.
#'     For \code{highlightChromPeaks}: \code{character(1)} defining how the peak
#'     should be highlighted: \code{type = "rect"} draws a rectangle
#'     representing the peak definition, \code{type = "point"} indicates a
#'     chromatographic peak with a single point at the position of the peak's
#'     \code{"rt"} and \code{"maxo"}.
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
#' highlightChromPeaks(xod, rt = c(2700, 2900), mz = 335)
highlightChromPeaks <- function(x, rt, mz,
                                border = rep("00000040", length(fileNames(x))),
                                lwd = 1, col = NA, type = c("rect", "point"),
                                ...) {
    type <- match.arg(type)
    n_samples <- length(fileNames(x))
    if (missing(rt))
        rt <- c(-Inf, Inf)
    if (missing(mz))
        mz <- c(-Inf, Inf)
    if (!is(x, "XCMSnExp"))
        stop("'x' has to be a XCMSnExp object")
    if (!hasChromPeaks(x))
        stop("'x' does not contain any detected peaks")
    pks <- chromPeaks(x, rt = rt, mz = mz, ppm = 0)
    if (length(col) != n_samples)
        col <- rep(col[1], n_samples)
    if (length(border) != n_samples)
        border <- rep(border[1], n_samples)
    if (length(pks)) {
        if (type == "rect")
            rect(xleft = pks[, "rtmin"], xright = pks[, "rtmax"],
                 ybottom = rep(0, nrow(pks)), ytop = pks[, "maxo"],
                 border = border[pks[, "sample"]], lwd = lwd,
                 col = col[pks[, "sample"]])
        if (type == "point") {
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
            do.call("points", args = c(list(x = pks[, "rt"], y = pks[, "maxo"],
                                            col = col[pks[, "sample"]]), dots))
        }
    }
}


#' @title General visualizations of peak detection results
#' 
#' @description
#'
#' \code{plotChromPeakImage} plots the identified chromatographic
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
    pks <- chromPeaks(x_file, mz = ylim, rt = xlim)
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
    pks <- chromPeaks(x, rt = xlim)
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
    new_x@.processHistory <- unlist(procH)
    chromPeaks(new_x) <- pks
    if (validObject(new_x))
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

