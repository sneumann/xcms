#' @include DataClasses.R

##' Takes a XCMSnExp and drops ProcessHistory steps from the @.processHistory
##' slot matching the provided type.
##'
##' @return The XCMSnExp input object with selected ProcessHistory steps dropped.
##' @noRd
dropProcessHistories <- function(x, type) {
    ## ## Drop processing history steps by type.
    ## if (!missing(type)) {
    ##     toRem <- unlist(lapply(processHistory(x), function(z) {
    ##         return(processType(z) %in% type)
    ##     }))
    ##     if (any(toRem))
    ##         x@.processHistory <- processHistory(x)[!toRem]
    ## }
    x@.processHistory <- dropProcessHistoriesList(processHistory(x), type = type)
    return(x)
}

dropProcessHistoriesList <- function(x, type) {
    if (!missing(type)) {
        toRem <- unlist(lapply(x, function(z) {
            return(processType(z) %in% type)
        }))
        if (any(toRem))
            x <- x[!toRem]
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
        xs@groups <- as.matrix(fgs[, -ncol(fgs)])
        xs@groupidx <- fgs$featureidx
    }
    ## @rt combination from rtime(x) and adjustedRtime(x)
    rts <- list()
    rts$raw <- rtime(from, bySample = TRUE)
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

