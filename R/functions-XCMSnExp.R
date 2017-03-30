#' @include DataClasses.R functions-utils.R

##' Takes a XCMSnExp and drops ProcessHistory steps from the @.processHistory
##' slot matching the provided type.
##'
##' @param num which should be dropped? If \code{-1} all matching will be dropped,
##' otherwise just the most recent num.
##' @return The XCMSnExp input object with selected ProcessHistory steps dropped.
##' @noRd
dropProcessHistories <- function(x, type, num = -1) {
    x@.processHistory <- dropProcessHistoriesList(processHistory(x),
                                                  type = type, num = num)
    return(x)
}

dropProcessHistoriesList <- function(x, type, num = -1) {
    if (!missing(type)) {
        toRem <- unlist(lapply(x, function(z) {
            return(processType(z) %in% type)
        }))
        if (any(toRem)) {
            if (num < 0) {
                x <- x[!toRem]
            } else {
                idx <- which(toRem)
                idx <- tail(idx, n = num)
                if (length(idx))
                x <- x[-idx]
            }
        }
    }
    return(x)
}

##' Convert an XCMSnExp to an xcmsSet.
##' @noRd
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

##' @title Extract spectra subsets
##'
##' Extract subsets of spectra matching the given mz and retention time ranges.
##' @noRd
.spectraSubsets <- function(x, rtrange, mzrange) {
    ## o Should allow to provide single ranges, but also matrices of rtrange
    ##   and/or mzranges.
    ## o Perform the data fetching by file.
    ## o Return the (mz subsetted) spectra by file and by ranges.
    ## o Since data should be processed on a by-file basis representation of the
    ##   result as a list (files) of list(ranges) of list(spectra) seems to be
    ##   best.
    ## SEE runit.XCMSnExp.R,
}

## This is somewhat similar to the getEIC, just that it extracts for each
## mz/rt range pair a data.frame with rt, mz, intensity per sample.
## This version works on a single rt/mz range pair at a time.
## CHECK:
## 1) mz range outside.
## 2) rt range outside.
.extractMsData <- function(x, rtrange, mzrange) {
    ## Subset the OnDiskMSnExp
    fns <- fileNames(x)
    subs <- filterMz(filterRt(x, rt = rtrange), mz = mzrange)
    if (length(subs) == 0) {
        return(NULL)
    }
    fromF <- base::match(fileNames(subs), fns)
    ## Now extract mz-intensity pairs from each spectrum.
    ## system.time(
    ## suppressWarnings(
    ##     dfs <- spectrapply(tmp, as.data.frame)
    ## )
    ## ) ## 0.73sec
    ## system.time(
    suppressWarnings(
        dfs <- spectrapply(subs, FUN = function(z) {
            if (peaksCount(z))
                return(base::data.frame(rt = rep_len(rtime(z), length(z@mz)),
                                        as.data.frame(z)))
            else
                return(base::data.frame(rt = numeric(), mz = numeric(),
                                        i = integer()))
        })
    )
    ## Now I want to rbind the spectrum data frames per file
    ff <- fromFile(subs)
    L <- split(dfs, f = ff)
    L <- lapply(L, do.call, what = rbind)
    ## Put them into a vector same length that we have files.
    res <- vector(mode = "list", length = length(fns))
    res[fromF] <- L
    return(res)
}

## Same as above, but we're applying a function - or none.
.sliceApply <- function(x, FUN = NULL, rtrange, mzrange) {
    fns <- fileNames(x)
    if (is(x, "XCMSnExp")) {
        ## Now, the filterRt might get heavy for XCMSnExp objects if we're
        ## filtering also the chromatographic peaks and features!
        msfd <- new("MsFeatureData")
        if (hasAdjustedRtime(x)) {
            ## just copy over the retention time.
            msfd$adjustedRtime <- x@msFeatureData$adjustedRtime
        }
        lockEnvironment(msfd, bindings = TRUE)
        x@msFeatureData <- msfd
        ## with an XCMSnExp without chrom. peaks and features filterRt
        ## should be faster.
    }
    subs <- filterMz(filterRt(x, rt = rtrange), mz = mzrange)
    if (base::length(subs) == 0) {
        return(list())
    }
    fromF <- base::match(fileNames(subs), fns)
    suppressWarnings(
        res <- spectrapply(subs, FUN = FUN)
    )
    ## WARN: fromFile reported within the Spectra might not be correct!!!
    return(res)
}

##' @description Extract a chromatogram from an \code{OnDiskMSnExp} or
##' \code{XCMSnExp} subsetting to the provided retention time range
##' (\code{rt}) and using the function \code{aggregationFun} to aggregate
##' intensity values for the same retention time across the mz range
##' (\code{mz}).
##'
##' @param x An \code{OnDiskMSnExp} or \code{XCMSnExp} object.
##'
##' @param rt \code{numeric(2)} providing the lower and upper retention time. It
##' is also possible to submit a \code{numeric(1)} in which case \code{range} is
##' called on it to transform it to a \code{numeric(2)}.
##' 
##' @param mz \code{numeric(2)} providing the lower and upper mz value for
##' the mz range. It is also possible to submit a \code{numeric(1)} in which case
##' \code{range} is called on it to transform it to a \code{numeric(2)}.
##' 
##' @param aggregationFun The function to be used to aggregate intensity values
##' across the mz range for the same retention time.
##'
##' @param ... Additional arguments to be passed to the object's \code{rtime}
##' call.
##' 
##' @return A \code{list} with the \code{Chromatogram} objects. If no data was
##' present for the specified \code{rtrange} and \code{mzrange} the function
##' returns a \code{list} of length \code{0}.
##'
##' @author Johannes Rainer
##' @noRd
.extractChromatogram <- function(x, rt, mz, aggregationFun = "sum", ...) {
    if (!any(.SUPPORTED_AGG_FUN_CHROM == aggregationFun))
        stop("'aggregationFun' should be one of ",
             paste0("'", .SUPPORTED_AGG_FUN_CHROM, "'", collapse = ", "))
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
    ## Subset the object based on rt and mz range.
    subs <- filterMz(filterRt(x, rt = rt), mz = mz)
    if (length(subs) == 0) {
        return(list())
    }
    ## Now, call spectrapply on the object to return the data we need from each
    ## Spectrum: the aggregated intensity values per spectrum and the mz value
    ## range.
    ## Note: we're returning NA in case we don't have a valid measurement for a
    ## retention time within the specified mz
    suppressWarnings(
        res <- spectrapply(subs, FUN = function(z) {
            if (!z@peaksCount)
                return(c(NA_real_, NA_real_, NA_real_))
            ## return(list())
            return(c(range(z@mz), do.call(aggregationFun, list(z@intensity))))
        })
    )
    ## Do I want to drop the names?
    nas <- unlist(lapply(res, function(z) is.na(z[3])), use.names = FALSE)
    if (all(nas))
        return(list())
    ## not_empty <- base::which(base::lengths(res) > 0)
    ## if (length(not_empty)) {
    ## res <- split(res[not_empty], f = fromFile(subs)[not_empty])
    ## rtm <- split(rtime(subs, ...)[not_empty], f = fromFile(subs)[not_empty])
    res <- split(res, f = fromFile(subs))
    rtm <- split(rtime(subs, ...), f = fromFile(subs))
    ## We want to have one Chromatogram per file.
    ## Let's use a simple for loop here - no need for an mapply (yet).
    resL <- vector("list", length(res))
    for (i in 1:length(res)) {
        allVals <- unlist(res[[i]], use.names = FALSE)
        idx <- seq(3, length(allVals), by = 3)
        mzr <- range(allVals[-idx], na.rm = TRUE, finite = TRUE)
        ## Or should we drop the names completely?
        ints <- allVals[idx]
        names(ints) <- names(rtm[[i]])
        resL[[i]] <- Chromatogram(rtime = rtm[[i]],
                                  intensity = ints, mz = mzr,
                                  filterMz = fmzr,
                                  fromFile = as.integer(names(res)[i]),
                                  aggregationFun = aggregationFun)
    }
    return(resL)
    ## } else {
    ##     return(list())
    ## }
}

#' Integrates the intensities for chromatograpic peak(s). This is supposed to be
#' called by the fillChromPeaks method.
#'
#' @note Use one of .getPeakInt2 or .getPeakInt3 instead!
#' 
#' @param object An \code{XCMSnExp} object representing a single sample.
#' @param peakArea A \code{matrix} with the peak definition, i.e. \code{"rtmin"},
#' \code{"rtmax"}, \code{"mzmin"} and \code{"mzmax"}.
#' @noRd
.getPeakInt <- function(object, peakArea) {
    if (length(fileNames(object)) != 1)
        stop("'object' should be an XCMSnExp for a single file!")
    res <- numeric(nrow(peakArea))
    for (i in 1:length(res)) {
        rtr <- peakArea[i, c("rtmin", "rtmax")]
        chr <- extractChromatograms(object,
                                    rt = rtr,
                                    mz = peakArea[i, c("mzmin", "mzmax")])[[1]]
        if (length(chr))
            res[i] <- sum(intensity(chr), na.rm = TRUE) *
                ((rtr[2] - rtr[1]) / (length(chr) - 1))
        else
            res[i] <- NA_real_
    }
    return(unname(res))
}

#' Integrates the intensities for chromatograpic peak(s). This is supposed to be
#' called by the fillChromPeaks method.
#'
#' @note This reads the full data first and does the subsetting later in R.
#' 
#' @param object An \code{XCMSnExp} object representing a single sample.
#' @param peakArea A \code{matrix} with the peak definition, i.e. \code{"rtmin"},
#' \code{"rtmax"}, \code{"mzmin"} and \code{"mzmax"}.
#' @noRd
.getPeakInt2 <- function(object, peakArea) {
    if (length(fileNames(object)) != 1)
        stop("'object' should be an XCMSnExp for a single file!")
    res <- numeric(nrow(peakArea))
    spctr <- spectra(object, BPPARAM = SerialParam())
    mzs <- lapply(spctr, mz)
    valsPerSpect <- lengths(mzs)
    ints <- unlist(lapply(spctr, intensity), use.names = FALSE)
    rm(spctr)
    mzs <- unlist(mzs, use.names = FALSE)
    rtim <- rtime(object)
    for (i in 1:length(res)) {
        rtr <- peakArea[i, c("rtmin", "rtmax")]
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
                res[i] <- sum(mtx[, 3], na.rm = TRUE) *
                    ((rtr[2] - rtr[1]) /
                     (sum(rtim >= rtr[1] & rtim <= rtr[2]) - 1))
            } else {
                res[i] <- NA_real_
            }
        } else {
            res[i] <- NA_real_
        }
    }
    return(unname(res))
}

#' Integrates the intensities for chromatograpic peak(s). This is supposed to be
#' called by the fillChromPeaks method.
#'
#' @note This reads the full data first and does the subsetting later in R. This
#' function uses the C getEIC function.
#' 
#' @param object An \code{XCMSnExp} object representing a single sample.
#' @param peakArea A \code{matrix} with the peak definition, i.e. \code{"rtmin"},
#' \code{"rtmax"}, \code{"mzmin"} and \code{"mzmax"}.
#'
#' @noRd
.getPeakInt3 <- function(object, peakArea) {
    if (length(fileNames(object)) != 1)
        stop("'object' should be an XCMSnExp for a single file!")
    if (nrow(peakArea) == 0) {
        return(numeric())
    }
    res <- matrix(ncol = 4, nrow = nrow(peakArea))
    res <- numeric(nrow(peakArea))
    spctr <- spectra(object, BPPARAM = SerialParam())
    mzs <- lapply(spctr, mz)
    valsPerSpect <- lengths(mzs)
    scanindex <- valueCount2ScanIndex(valsPerSpect) ## Index vector for C calls
    ints <- unlist(lapply(spctr, intensity), use.names = FALSE)
    rm(spctr)
    mzs <- unlist(mzs, use.names = FALSE)
    rtim <- rtime(object)
    for (i in 1:length(res)) {
        rtr <- peakArea[i, c("rtmin", "rtmax")]
        sr <- c(min(which(rtim >= rtr[1])), max(which(rtim <= rtr[2])))
        eic <- .Call("getEIC", mzs, ints, scanindex,
                     as.double(peakArea[i, c("mzmin", "mzmax")]),
                     as.integer(sr), as.integer(length(scanindex)),
                     PACKAGE = "xcms")
        if (length(eic$intensity)) {
            ## How to calculate the area: (1)sum of all intensities / (2)by
            ## the number of data points (REAL ones, considering also NAs)
            ## and multiplied with the (3)rt width.
            if (!all(is.na(eic$intensity)) && !all(eic$intensity == 0)) {
                res[i] <- sum(eic$intensity, na.rm = TRUE) *
                    ((rtr[2] - rtr[1]) / (length(eic$intensity) - 1))
            } else {
                res[i] <- NA_real_
            }
        } else {
            res[i] <- NA_real_
        }
    }
    return(unname(res))
}


#' Integrates the intensities for chromatograpic peak(s). This is supposed to be
#' called by the fillChromPeaks method.
#'
#' @note This reads the full data first and does the subsetting later in R.
#' 
#' @param object An \code{XCMSnExp} object representing a single sample.
#' 
#' @param peakArea A \code{matrix} with the peak definition, i.e. \code{"rtmin"},
#' \code{"rtmax"}, \code{"mzmin"} and \code{"mzmax"}.
#'
#' @param sample_idx \code{integer(1)} with the index of the sample in the
#' object.
#'
#' @param mzCenterFun Name of the function to be used to calculate the mz value.
#' Defaults to \code{weighted.mean}, i.e. the intensity weighted mean mz.
#' 
#' @param cn \code{character} with the names of the result matrix.
#' 
#' @return A \code{matrix} with at least columns \code{"mz"}, \code{"rt"},
#' \code{"into"} and \code{"maxo"} with the by intensity weighted mean of mz,
#' rt or the maximal intensity in the area, the integrated signal in the area
#' and the maximal signal in the area.
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

#' Same as getChromPeakData, just without retention time.
#' @note The mz and maxo are however estimated differently than for the
#' getChromPeakData: mz is the mz closest to the median mz of the feature and
#' maxo its intensity.
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

#' @description Simple helper function to extract the peakidx column from the
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

#' @description \code{adjustRtimePeakGroups} returns the features (peak groups)
#'     which would, depending on the provided \code{\link{PeakGroupsParam}}, be
#'     selected for alignment/retention time correction.
#'
#' @note \code{adjustRtimePeakGroups} is supposed to be called \emph{before} the
#'     sample alignment, but after a correspondence (peak grouping).
#'
#' @return For \code{adjustRtimePeakGroups}: a \code{matrix}, rows being
#'     features, columns samples, of retention times. The features are ordered
#'     by the median retention time across columns.
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
#' @description Plot the difference between the adjusted and the raw retention
#'     time (y-axis) for each file along the (adjusted) retention time (x-axis).
#'     It the alignment was performed using the
#'     \code{\link{adjustRtime-peakGroups}} method, also the features (peak
#'     groups) used for the alignment are shown.
#'
#' @param object A \code{\link{XCMSnExp}} object with the alignment results.
#'
#' @param ... Additional arguments to be passed down to the plotting function.
#' 
#' @seealso \code{\link{adjustRtime}} for all retention time correction/
#'     alignment methods.
#' 
#' @noRd
plotAdjustedRtime <- function(object, col = "#00000080", lty = 1, type = "l",
                              adjRt = TRUE,
                              xlab = ifelse(adjRt, yes = expression(rt[adj]),
                                           no = expression(rt[raw])),
                              ylab = expression(rt[adj]-rt[raw]), ...) {
    if (!is(object, "XCMSnExp"))
        stop("'object' has to be an 'XCMSnExp' object.")
    if (!hasAdjustedRtime(object))
        warning("No alignment/retention time correction results present.")
    ## No legend or anything - just plotting.
    diffRt <- rtime(object, adjusted = TRUE) - rtime(object, adjusted = FALSE)
    diffRt <- split(diffRt, fromFile(object))
    ## Define the rt that is shown on x-axis
    rawRt <- rtime(object, adjusted = FALSE, bySample = TRUE)
    adjRt <- rtime(object, adjusted = TRUE, bySample = TRUE)
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
    plot(3, 3, pch = NA, xlim = range(rawRt, na.rm = TRUE),
         ylim = range(diffRt, na.rm = TRUE), xlab = xlab, ylab = ylab, ...)
    ## Plot all.
    for (i in 1:length(diffRt))
        points(x = rawRt[[i]], y = diffRt[[i]], col = col[i], lty = lty[i],
               type = type)
    ## If alignment was performed using the peak groups method highlight also
    ## those in the plot. Problem is these peak groups don't necessarily have
    ## to be present!
    ## Thus:
    ## 1) Check if grouping was performed for alignment.
    ## 2) If TRUE, check that grouping was performed using the peak groups
    ##    method.
    ## 3) If TRUE: drop the adjusted retention times, perform the grouping
    ##    using the provided setting and return the peak groups matrix.
    procStepTypes <- lapply(processHistory(object), processType)
    idx_algn <- which(procStepTypes == xcms:::.PROCSTEP.RTIME.CORRECTION)
    if (!length(idx_algn))
        idx_algn <- 0
    idx_grp <- which(procStepTypes == xcms:::.PROCSTEP.PEAK.GROUPING)
    if (!length(idx_grp))
        idx_grp <- 0
    if (idx_grp < idx_algn) {
        ph <- processHistory(object)[[idx_grp]]
        if (is(ph, "XProcessHistory")) {
            ## Check what param class we've got
            prm <- processParam(ph)
            if (is(prm, "PeakGroupsParam")) {
                ## Now, replay the grouping.
                message("Replaying correspondence analysis prior to alignment.")
                message(" - ", appendLF = FALSE)
                tmp <- dropAdjustedRtime(object)
                message(" - ", appendLF = FALSE)
                tmp <- groupChromPeaks(tmp, param = prm)
                pkGroup <- adjustRtimePeakGroups(tmp)
                ## OK, I've got the raw retention times now, need to get the
                ## adjusted too.
                pkGroupAdj <- pkGroup
                for (i in 1:ncol(pkGroupAdj)) {
                    pkGroupAdj[, i] <- xcms:::.applyRtAdjustment(pkGroup[, i],
                                                          rawRt[[i]],
                                                          adjRt[[i]])
                }
                ## matplot?, but before that, re-order each row.
            }
        }
    }
}

