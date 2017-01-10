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
.XCMSnExp2xcmsSet <- function(x) {
    xs <- new("xcmsSet")
    ## @peaks <- features
    if (hasDetectedFeatures(x))
        xs@peaks <- features(x)
    ## @groups <- part of featureGroups
    ## @groupidx <- featureGroups(x)$featureidx
    if (hasAlignedFeatures(x)){
        fgs <- featureGroups(x)
        xs@groups <- as.matrix(fgs[, -ncol(fgs)])
        xs@groupidx <- fgs$featureidx
    }
    ## @rt combination from rtime(x) and adjustedRtime(x)
    rts <- list()
    rts$raw <- rtime(x, bySample = TRUE)
    if (hasAdjustedRtime(x))
        rts$corrected <- adjustedRtime(x, bySample = TRUE)
    else
        rts$corrected <- rts$raw

    xs@rt <- rts

    ## @filled ... not yet.
    ## @phenoData <- phenoData?
    ## @filepaths <- fileNames(x) ?
    ## @profinfo (list)
    profMethod <- "bin"
    profStep <- 0.1
    profParam <- list()
    ## If we've got any MatchedFilterParam we can take the values from there
    profinfo(xs) <- c(list(method = profMethod, step = profStep), profParam)
    ## @dataCorrection (numeric) ? in xcmsSet function, if lockMassFreq.
    ## @polarity (character) ?
    ## @progressInfo skip
    ## @progressCallback skip
    ## @mslevel <- msLevel?
    ## @scanrange <- ?
    ## .processHistory <- ? or coerced from XCMSnExp.
    if (validObject(xs))
        return(xs)
}
