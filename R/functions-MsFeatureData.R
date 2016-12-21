## Functions for MsFeatureData classes.
#' @include DataClasses.R

##' @description Performs a validation check of all elements within the object:
##' 1) Allowed are: features (matrix), featureGroups (DataFrame) and
##'    adjustedRtime (list).
##' @author Johannes Rainer
##' @return \code{TRUE} if object is valid, or a message with the error message.
##' @noRd
validateMsFeatureData <- function(x) {
    msg <- validMsg(NULL, NULL)
    ks <- ls(x)
    if (length(ks)) {
        validKeys <- ks %in% c("features", "featureGroups", "adjustedRtime")
        if (!all(validKeys)) {
            msg <- validMsg(msg, paste0("Only elements named 'features', ",
                                        "'featureGroups' and 'adjustedRtime' ",
                                        "are allowed but I got ",
                                        paste0("'", ks[!validKeys],"'")))
        }
        haveFts <- any(ks == "features")
        if (haveFts) {
            if (is.matrix(x$features)) {
                hasReqCols <- .XCMS_REQ_FEATS_COLS %in% colnames(x$features)
                if (any(!hasReqCols))
                    msg <- validMsg(msg,
                                    paste0("Required columns ",
                                           paste0("'", .XCMS_REQ_FEATS_COLS[!hasReqCols],
                                                  "'", collapse = ", "), " not",
                                           " present in 'features' matrix!"))
            } else {
                msg <- validMsg(msg, paste0("The 'features' element has to be ",
                                            "of type 'matrix' and not '",
                                            class(x$features), "'!"))
            }
        }
        haveFGs <- any(ks == "featureGroups")
        if (haveFGs) {
            if (is(x$featureGroups, "DataFrame")) {
                ## Check required columns.
                hasReqCols <- .XCMS_REQ_FEATG_COLS %in% colnames(x$featureGroups)
                if (any(!hasReqCols)) {
                    msg <- validMsg(msg,
                                    paste0("Required columns ",
                                           paste0("'", .XCMS_REQ_FEATG_COLS[!hasReqCols],
                                                  "'", collapse = ", "), " not",
                                           " present in 'featureGroups'!"))
                } else {
                    ## Check content!
                    if (!is(x$featureGroups$featureidx, "list"))
                        msg <- validMsg(msg,
                                        paste0("Column 'featureidx' in '",
                                               "featureGroups' is not a list!"))
                    for (col in .XCMS_REQ_FEATG_COLS[-length(.XCMS_REQ_FEATG_COLS)]) {
                        if (!is.numeric(x$featureGroups[, col]))
                            msg <- validMsg(msg, paste0("Column '", col, "' has",
                                                        " to be numeric!"))
                    }
                    if (haveFts) {
                        ## Check that indices are within 1:nrow(x$features)
                        if (!all(unlist(x$featureGroups$featureidx) %in%
                                 1:nrow(x$features)))
                            msg <- validMsg(msg,
                                            paste0("Some of the indices in column",
                                                   " 'featureidx' of element ",
                                                   "'featureGroups' do not match ",
                                                   "rows of the 'features' matrix!"))
                    }
                }

            } else {
                msg <- validMsg(msg, paste0("The 'featureGroups' element has to",
                                            " be of type 'DataFrame' and not '",
                                            class(x$featureGroups), "'!"))
            }
            if (!haveFts) {
                msg <- validMsg(msg, paste0("Can not have element 'featureGroups'",
                                            " without element 'features'!"))
            }
        }
        haveRts <- any(ks == "adjustedRtime")
        if (haveRts) {
            if (!haveFts)
                msg <- validMsg(msg, paste0("Can not have element 'adjustedRtime'",
                                            " without element 'features'!"))
            ## adjustedRtime has to be a list of numerics.
            if (!is.list(x$adjustedRtime)) {
                msg <- validMsg(msg, paste0("The 'alignedRtime' element has to ",
                                            "be of type 'list' and not '",
                                            class(x$adjustedRtime), "'!"))
            } else {
                areNum <- unlist(lapply(x$adjustedRtime, function(z) {
                    return(is.numeric(z))
                }))
                if (!all(areNum))
                    msg <- validMsg(msg, paste0("The 'alignedRtime' element has",
                                                " to be a list of numeric ",
                                                "vectors!"))
            }
        }
    }
    if (is.null(msg))
        return(TRUE)
    else return(msg)
}
