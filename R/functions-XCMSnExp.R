#' @include DataClasses.R

##' Takes a XCMSnExp and drops ProcessHistory steps from the @.processHistory
##' slot matching the provided type.
##'
##' @param num which should be dropped? If \code{-1} all matching will be dropped,
##' otherwise just the most recent num.
##' @return The XCMSnExp input object with selected ProcessHistory steps dropped.
##' @noRd
dropProcessHistories <- function(x, type, num = -1) {
    ## ## Drop processing history steps by type.
    ## if (!missing(type)) {
    ##     toRem <- unlist(lapply(processHistory(x), function(z) {
    ##         return(processType(z) %in% type)
    ##     }))
    ##     if (any(toRem))
    ##         x@.processHistory <- processHistory(x)[!toRem]
    ## }
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
    spctr <- spectra(object)
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
    spctr <- spectra(object)
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
#' @param cn \code{character} with the names of the result matrix.
#' 
#' @return A \code{matrix} with at least columns \code{"mz"}, \code{"rt"},
#' \code{"into"} and \code{"maxo"} with the mz and rt or the maximal intensity
#' in the area, the integrated signal in the area and the maximal signal in the
#' area.
#' @noRd
.getChromPeakData <- function(object, peakArea, sample_idx,
                              cn = c("mz", "rt", "into", "maxo", "sample")) {
    if (length(fileNames(object)) != 1)
        stop("'object' should be an XCMSnExp for a single file!")
    ncols <- length(cn)
    res <- matrix(ncol = ncols, nrow = nrow(peakArea))
    colnames(res) <- cn
    res[, "sample"] <- sample_idx
    res[, c("mzmin", "mzmax", "rtmin", "rtmax")] <-
        peakArea[, c("mzmin", "mzmax", "rtmin", "rtmax")]
    ## Load the data
    message("Reguesting ", nrow(res), " missing peaks from ",
             basename(fileNames(object)), " ... ", appendLF = FALSE)
    spctr <- spectra(object)
    mzs <- lapply(spctr, mz)
    valsPerSpect <- lengths(mzs)
    ints <- unlist(lapply(spctr, intensity), use.names = FALSE)
    rm(spctr)
    mzs <- unlist(mzs, use.names = FALSE)
    rtim <- rtime(object)
    for (i in 1:nrow(res)) {
        rtr <- peakArea[i, c("rtmin", "rtmax")]
        mtx <- xcms:::.rawMat(mz = mzs, int = ints, scantime = rtim,
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
                res[i, c("rt", "mz", "maxo")] <- mtx[maxi[1], ]
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

