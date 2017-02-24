## Unsorted utility functions.
#' @include DataClasses.R

############################################################
## valueCount2ScanIndex
##
## @description Simple helper function that converts the number of values
## per scan/spectrum to an integer vector that can be passed to the base
## xcms functions/downstream C functions.
##
## @title Create index vector for internal C calls
## @param valCount Numeric vector representing the number of values per
## spectrum.
## @return An integer vector with the index (0-based) in the mz or intensity
## vectors indicating the start of a spectrum.
## @author Johannes Rainer
valueCount2ScanIndex <- function(valCount){
    ## Convert into 0 based.
    valCount <- cumsum(valCount)
    return(as.integer(c(0, valCount[-length(valCount)])))
}

############################################################
## useOriginalCode
##
## Simple function allowing the user to enable using the orignal
## code instead of the new implementations.
## This sets options.
##
##' @title Enable usage of old xcms code
##'
##' @description This function allows to enable the usage of old, partially
##' deprecated code from xcms by setting a corresponding global option. See
##' details for functions affected.
##'
##' @note Usage of old code is strongly dicouraged. This function is thought
##' to be used mainly in the transition phase from xcms to xcms version 3.
##' @details The functions/methods that will be affected by this are:
##' \itemize{
##' \item \code{\link{do_findChromPeaks_matchedFilter}}
##' }
##' @param x logical(1) to specify whether or not original
##' old code should be used in corresponding functions. If not provided the
##' function simply returns the value of the global option.
##' @return logical(1) indicating whether old code is being
##' used.
##' @author Johannes Rainer
useOriginalCode <- function(x) {
    if (missing(x)) {
        res <- options()$BioC$xcms$useOriginalCode
        if (is.null(res))
            return(FALSE)
        return(res)
    }
    if (!is.logical(x))
        stop("'x' has to be logical.")
    b_opts <- getOption("BioC")
    x_opts <- b_opts$xcms
    x_opts$useOriginalCode <- x[1]
    b_opts$xcms <-  x_opts
    options(BioC = b_opts)
    return(options()$BioC$xcms$useOriginalCode)
}

## .getOriginalFunction <- function(x) {
##     if (!any(names(.ORIGINAL_FUNCTIONS)))
## }
## .ORIGINAL_FUNCTIONS <- c(
##     matchedFilter = ".matchedFilter_orig"
## )

##' @title Copy the content from an environment to another one
##' This function copies the content of an environment into another one.
##' @param env environment from which to copy.
##' @param inheritLocks logical(1) whether the locking status should be copied
##' too
##' @return an env.
##' @noRd
.copy_env <- function(env, inheritLocks = FALSE) {
    new_e <- new.env(parent = emptyenv())
    eNames <- ls(env, all.names = TRUE)
    if (length(eNames) > 0) {
        for (eN in eNames) {
            new_e[[eN]] <- env[[eN]]
        }
    }
    if (inheritLocks) {
        if (environmentIsLocked(env))
            lockEnvironment(new_e)
    }
    return(new_e)
}
