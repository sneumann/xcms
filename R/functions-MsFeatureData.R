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

##' @description Filter features and sync them with with the present
##' filterGroups, i.e. update their featureidx column or remove them.
##'
##' @author Johannes Rainer
##' @noRd
.filterFeatures <- function(x, idx) {
    if (missing(idx))
        return(x)
    if (!hasDetectedFeatures(x))
        return(x)
    fts <- features(x)
    idx <- sort(idx)
    if (!all(idx %in% 1:nrow(fts)))
        stop("All indices in 'idx' have to be within 1 and nrow of the feature",
             " matrix.")
    new_e <- new("MsFeatureData")
    features(new_e) <- fts[idx, , drop = FALSE]
    if (hasAlignedFeatures(x)) {
        af <- featureGroups(x)
        af <- split(af, 1:nrow(af))
        afL <- lapply(af, function(z) {
            if(all(z$featureidx[[1]] %in% idx)) {
                z$featureidx <- list(match(z$featureidx[[1]], idx))
                return(z)
            } else {
                return(NULL)
            }
        })
        af <- do.call(rbind, afL)
        featureGroups(new_e) <- af
    }
    if (hasAdjustedRtime(x))
        adjustedRtime(new_e) <- adjustedRtime(x)
    return(new_e)
}
