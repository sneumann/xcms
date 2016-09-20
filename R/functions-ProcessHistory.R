############################################################
## Functions for ProcessHistory objects

############################################################
## Constructor
ProcessHistory <- function(type., date., info., error., fileIndex.) {
    if (missing(type.))
        type. <- .PROCSTEP.UNKNOWN
    if (missing(info.))
        info. <- character()
    if (missing(error.))
        error. <- NULL
    if (missing(date.))
        date. <- date()
    if (missing(fileIndex.))
        fileIndex. <- integer()
    return(new("ProcessHistory", type = type., info = info.,
               date = date., error = error.,
               fileIndex = as.integer(fileIndex.)))
}

############################################################
## updateFileIndex
## Update the file index mapping index in 'old' to index in 'new' dropping
## all indices that are not in 'new'
updateFileIndex <- function(x, old = integer(), new = integer()) {
    if (length(old) == 0 & length(new) == 0)
        return(x)
    if (length(old) != length(new))
        stop("Lengths of 'old' and 'new' don't match!")
    fidx <- x@fileIndex
    fidx <- fidx[fidx %in% old]
    if (length(fidx) == 0)
        stop("None of the file indices matches 'old'!")
    for (i in 1:length(fidx)) {
        fidx[i] <- new[old == fidx[i]]
    }
    x@fileIndex <- fidx
    return(x)
}
