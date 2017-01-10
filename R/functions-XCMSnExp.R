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
