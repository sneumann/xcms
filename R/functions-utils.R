## Unsorted utility functions.
#' @include DataClasses.R

############################################################
## valueCount2ScanIndex
##
##' @description Simple helper function that converts the number of values
##' per scan/spectrum to an integer vector that can be passed to the base
##' xcms functions/downstream C functions.
##'
##' @title Create index vector for internal C calls
##' @param valCount Numeric vector representing the number of values per
##' spectrum.
##' @return An integer vector with the index (0-based) in the mz or intensity
##' vectors indicating the start of a spectrum.
##' @author Johannes Rainer
valueCount2ScanIndex <- function(valCount){
    ## Convert into 0 based.
    valCount <- cumsum(valCount)
    return(as.integer(c(0, valCount[-length(valCount)])))
}
