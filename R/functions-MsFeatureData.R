## Functions for MsFeatureData classes.
#' @include DataClasses.R

#' Validates a 'chromPeaks' matrix or data.frame and ensures that it contains all
#' required columns and that all columns are of numeric data type.
#' 
#' @return \code{TRUE} or a \code{character} with the error message.
#' 
#' @noRd
.validChromPeaksMatrix <- function(x) {
    msg <- character()
    if (length(x)) {
        if (!(is.matrix(x) | is.data.frame(x)))
            return(paste0("'chromPeaks' has to be a matrix or a data.frame!"))
        hasReqCols <- .REQ_PEAKS_COLS %in% colnames(x)
        if (any(!hasReqCols))
            return(paste0("Required columns ",
                          paste0("'", .REQ_PEAKS_COLS[!hasReqCols],
                                 "'", collapse = ", "), " not",
                          " present in 'chromPeaks' matrix!"))
        ## Check data.types - all have to be numeric.
        typeOK <- apply(x[, .REQ_PEAKS_COLS, drop = FALSE], MARGIN = 2,
                        is.numeric)
        if (any(!typeOK))
            return(paste0("Values in column(s) ",
                          paste0("'", names(typeOK)[!typeOK], "'", collapse = ", ")),
                   " of the 'chromPeaks' matrix are not numeric!")
    }
    return(TRUE)
}


#' @description Performs a validation check of all elements within the object:
#' 1) Allowed are: chromPeaks (matrix), featureDefinitions (DataFrame) and
#'    adjustedRtime (list).
#' 
#' @author Johannes Rainer
#' 
#' @return \code{TRUE} if object is valid, or a message with the error message.
#' 
#' @noRd
validateMsFeatureData <- function(x) {
    msg <- character()
    ks <- ls(x)
    if (length(ks)) {
        validKeys <- ks %in% c("chromPeaks", "featureDefinitions",
                               "adjustedRtime")
        if (!all(validKeys)) {
            msg <- c(msg, paste0("Only elements named 'chromPeaks', ",
                                 "'featureDefinitions' and 'adjustedRtime' ",
                                 "are allowed but I got ",
                                 paste0("'", ks[!validKeys],"'")))
        }
        haveFts <- any(ks == "chromPeaks")
        if (haveFts) {
            OK <- .validChromPeaksMatrix(x$chromPeaks)
            if (is.character(OK))
                msg <- c(msg, OK)
        }
        haveFGs <- any(ks == "featureDefinitions")
        if (haveFGs) {
            if (is(x$featureDefinitions, "DataFrame")) {
                ## Check required columns.
                hasReqCols <- .REQ_PEAKG_COLS %in% colnames(x$featureDefinitions)
                if (any(!hasReqCols)) {
                    msg <- c(msg,
                             paste0("Required columns ",
                                    paste0("'", .REQ_PEAKG_COLS[!hasReqCols],
                                           "'", collapse = ", "), " not",
                                    " present in 'featureDefinitions'!"))
                } else {
                    ## Check content!
                    if (!is(x$featureDefinitions$peakidx, "list"))
                        msg <- c(msg,
                                 paste0("Column 'peakidx' in '",
                                        "featureDefinitions' is not a list!"))
                    for (col in .REQ_PEAKG_COLS[-length(.REQ_PEAKG_COLS)]) {
                        if (!is.numeric(x$featureDefinitions[, col]))
                            msg <- c(msg, paste0("Column '", col, "' has",
                                                 " to be numeric!"))
                    }
                    if (haveFts) {
                        ## Check that indices are within 1:nrow(x$chromPeaks)
                        if (!all(unlist(x$featureDefinitions$peakidx) %in%
                                 1:nrow(x$chromPeaks)))
                            msg <- c(msg,
                                     paste0("Some of the indices in column",
                                            " 'peakidx' of element ",
                                            "'featureDefinitions' do not match ",
                                            "rows of the 'chromPeaks' matrix!"))
                    }
                }

            } else {
                msg <- c(msg, paste0("The 'featureDefinitions' element has to",
                                     " be of type 'DataFrame' and not '",
                                     class(x$featureDefinitions), "'!"))
            }
            if (!haveFts) {
                msg <- c(msg, paste0("Can not have element 'featureDefinitions'",
                                     " without element 'chromPeaks'!"))
            }
        }
        haveRts <- any(ks == "adjustedRtime")
        if (haveRts) {
            ## Not true, since obiwarp works without peaks.
            ## if (!haveFts)
            ##     msg <- c(msg, paste0("Can not have element 'adjustedRtime'",
            ##                          " without element 'chromPeaks'!"))
            ## adjustedRtime has to be a list of numerics.
            if (!is.list(x$adjustedRtime)) {
                msg <- c(msg, paste0("The 'alignedRtime' element has to ",
                                     "be of type 'list' and not '",
                                     class(x$adjustedRtime), "'!"))
            } else {
                areNum <- unlist(lapply(x$adjustedRtime, function(z) {
                    return(is.numeric(z))
                }))
                if (!all(areNum))
                    msg <- c(msg, paste0("The 'alignedRtime' element has",
                                         " to be a list of numeric ",
                                         "vectors!"))
            }
        }
    }
    ## if (length(msg) == 0)
    ##     return(TRUE)
    ## else return(msg)
    return(msg)
}

#' @description Filter chromPeaks and sync them with with the present
#'     filterGroups, i.e. update their peakidx column or remove them.
#'
#' @param x A \code{MsFeatureData} or an \code{XCMSnExp} object.
#' 
#' @param idx \code{numeric} with the indices of the chromatographic peaks to
#'     keep.
#'
#' @return A \code{MsFeatureData}.
#' 
#' @author Johannes Rainer
#' 
#' @noRd
.filterChromPeaks <- function(x, idx) {
    if (missing(idx))
        return(x)
    if (!hasChromPeaks(x))
        return(x)
    fts <- chromPeaks(x)
    idx <- sort(idx)
    if (!all(idx %in% 1:nrow(fts)))
        stop("All indices in 'idx' have to be within 1 and nrow of the peak",
             " matrix.")
    new_e <- new("MsFeatureData")
    chromPeaks(new_e) <- fts[idx, , drop = FALSE]
    if (hasFeatures(x)) {
        af <- featureDefinitions(x)
        af <- split(af, 1:nrow(af))
        afL <- lapply(af, function(z) {
            if(all(z$peakidx[[1]] %in% idx)) {
                z$peakidx <- list(match(z$peakidx[[1]], idx))
                return(z)
            } else {
                return(NULL)
            }
        })
        af <- do.call(rbind, afL)
        if (length(af) > 0)
            featureDefinitions(new_e) <- af
    }
    if (hasAdjustedRtime(x)) {
        if (is(x, "XCMSnExp"))
            adjustedRtime(new_e) <- adjustedRtime(x, bySample = TRUE)
        else
            adjustedRtime(new_e) <- adjustedRtime(x)

    }
    return(new_e)
}
