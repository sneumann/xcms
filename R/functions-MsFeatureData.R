## Functions for MsFeatureData classes.
#' @include DataClasses.R

##' Validates a 'feature' matrix or data.frame and ensures that it contains all
##' required columns and that all columns are of numeric data type.
##' @return \code{TRUE} or a \code{character} with the error message.
##' @noRd
.validFeatureMatrix <- function(x) {
    msg <- validMsg(NULL, NULL)
    if (length(x)) {
        if (!(is.matrix(x) | is.data.frame(x)))
            return(paste0("'features' has to be a matrix or a data.frame!"))
        hasReqCols <- .XCMS_REQ_FEATS_COLS %in% colnames(x)
        if (any(!hasReqCols))
            return(paste0("Required columns ",
                          paste0("'", .XCMS_REQ_FEATS_COLS[!hasReqCols],
                                 "'", collapse = ", "), " not",
                          " present in 'features' matrix!"))
        ## Check data.types - all have to be numeric.
        typeOK <- apply(x[, .XCMS_REQ_FEATS_COLS, drop = FALSE], MARGIN = 2,
                        is.numeric)
        if (any(!typeOK))
            return(paste0("Values in column(s) ",
                          paste0("'", names(typeOK)[!typeOK], "'", collapse = ", ")),
                   " of the 'features' matrix are not numeric!")
    }
    return(TRUE)
}


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
            OK <- .validFeatureMatrix(x$features)
            if (is.character(OK))
                msg <- validMsg(msg, OK)
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
##' @param x A \code{MsFeatureData} or an \code{XCMSnExp} object.
##' @param idx \code{numeric} with the indices of the features to keep.
##'
##' @return A \code{MsFeatureData}.
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
    if (hasAdjustedRtime(x)) {
        if (is(x, "XCMSnExp"))
            adjustedRtime(new_e) <- adjustedRtime(x, bySample = TRUE)
        else
            adjustedRtime(new_e) <- adjustedRtime(x)

    }
    return(new_e)
}
