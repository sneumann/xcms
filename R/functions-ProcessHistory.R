############################################################
## Functions for ProcessHistory objects
#' @include DataClasses.R

#' @description \code{processHistoryTypes} returns the available \emph{types} of
#'     process histories. These can be passed with argument \code{type} to the
#'     \code{processHistory} method to extract specific process step(s).
#'
#' @rdname XCMSnExp-class
processHistoryTypes <- function() {
    .PROCSTEPS
}

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
## To remove indices, use NA in new
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
    x@fileIndex <- as.integer(fidx[!is.na(fidx)])
    return(x)
}

############################################################
## XProcessHistory
XProcessHistory <- function(param = NULL, msLevel = NA_integer_, ...) {
    obj <- ProcessHistory(...)
    obj <- as(obj, "XProcessHistory")
    obj@param <- param
    obj@msLevel <- as.integer(msLevel)
    classVersion(obj)["XProcessHistory"] <- "0.0.2"
    OK <- validObject(obj)
    if (is.character(OK))
        stop(OK)
    return(obj)
}

#' Create a simple generic process history object for a function call.
#'
#' @param fun `character` specifying the function name.
#'
#' @param args `list` with arguments to the function.
#'
#' @param date the date/time when the function was called.
#'
#' @param fileIndex `integer` with the indices of the files on which the
#'     function was called. Usually `1:length(fileNames(xs))`.
#'
#' @md
#' 
#' @noRd
GenericProcessHistory <- function(fun, args = list(), msLevel = NA_integer_,
                                  date. = date(), fileIndex. = NA_integer_) {
    gp <- new("GenericParam", fun = fun, args = args)
    xcms:::XProcessHistory(param = gp, msLevel = msLevel,
                           date. = date., fileIndex. = fileIndex.)
}

#' Remove a generic process history step based on the name of the function.
#'
#' @param x `list` of process history object (such as returned by
#'     `processHistory(xs)`).
#'
#' @param fun `character` with the (exact) name of the function.
#'
#' @return `list` of process history objects.
#'
#' @md
#'
#' @noRd
dropGenericProcessHistoryList <- function(x, fun) {
    got_it <- unlist(lapply(x, function(z) {
        if (is(z, "XProcessHistory")) {
            prm <- processParam(z)
            if (is(prm, "GenericParam"))
                prm@fun == fun
            else FALSE
        } else FALSE
    }))
    if (any(got_it))
        x <- x[!got_it]
    x
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

