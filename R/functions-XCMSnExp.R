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
    ## @peaks <- features
    if (hasDetectedFeatures(from))
        xs@peaks <- features(from)
    ## @groups <- part of featureGroups
    ## @groupidx <- featureGroups(x)$featureidx
    if (hasAlignedFeatures(from)){
        fgs <- featureGroups(from)
        xs@groups <- S4Vectors::as.matrix(fgs[, -ncol(fgs)])
        rownames(xs@groups) <- NULL
        xs@groupidx <- fgs$featureidx
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
    ph <- processHistory(from, type = .PROCSTEP.FEATURE.DETECTION)
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
        ## filtering also the features and featureGroups!
        msfd <- new("MsFeatureData")
        if (hasAdjustedRtime(x)) {
            ## just copy over the retention time.
            msfd$adjustedRtime <- x@msFeatureData$adjustedRtime
        }
        lockEnvironment(msfd, bindings = TRUE)
        x@msFeatureData <- msfd
        ## with an XCMSnExp without features and featureGroups filterRt should be
        ## faster.
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
    res <- spectrapply(subs, FUN = function(z) {
        if (!z@peaksCount)
            return(list())
        return(c(range(z@mz), do.call(aggregationFun, list(z@intensity))))
    })
    ## Do I want to drop the names?
    not_empty <- base::which(base::lengths(res) > 0)
    if (length(not_empty)) {
        res <- split(res[not_empty], f = fromFile(subs)[not_empty])
        rtm <- split(rtime(subs)[not_empty], f = fromFile(subs)[not_empty])
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
    } else {
        return(list())
    }
}
