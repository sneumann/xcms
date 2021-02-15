## Methods for the MsFeatureData class.
#' @include functions-MsFeatureData.R do_adjustRtime-functions.R

setValidity("MsFeatureData", function(object) {
    validateMsFeatureData(object)
})

#' @noRd
#'
#' @rdname XCMSnExp-class
setMethod("show", "MsFeatureData", function(object) {
    cat("Object of class: ", class(object), "\n")
    ks <- ls(object)
    if (length(ks)) {
        cat(" Available feature data:\n")
        for (k in ks)
            cat(" -", k, "\n")
    }
})

## (features) chromPeaks: getter and setter for the chromatographic peaks matrix.
## (featureDefinitions) featureDefinitions: getter and setter for the features DataFrame.
## adjustedRtime: getter and setter for the adjustedRtime list.

#' @noRd
#'
#' @rdname XCMSnExp-class
setMethod("hasAdjustedRtime", "MsFeatureData", function(object) {
    !is.null(object$adjustedRtime)
})

#' @noRd
#'
#' @rdname XCMSnExp-class
setMethod("hasFeatures", "MsFeatureData", function(object,
                                                   msLevel = integer()) {
    if (length(msLevel) && !is.null(object$featureDefinitions) &&
        any(colnames(object$featureDefinitions) == "ms_level"))
        any(msLevel %in% object$featureDefinitions$ms_level)
    else !is.null(object$featureDefinitions)
})

#' @noRd
#'
#' @rdname XCMSnExp-class
setMethod("hasChromPeaks", "MsFeatureData", function(object,
                                                     msLevel = integer()) {
    if (length(msLevel) && !is.null(object$chromPeaks) &&
        any(colnames(object$chromPeakData) == "ms_level"))
        any(msLevel %in% object$chromPeakData$ms_level)
    else !is.null(object$chromPeaks)
})

#' @noRd
#'
#' @rdname XCMSnExp-class
setMethod("adjustedRtime", "MsFeatureData", function(object) {
    if (hasAdjustedRtime(object))
        return(object$adjustedRtime)
    warning("No adjusted retention times available.")
    return(NULL)
})
#' @noRd
#'
#' @rdname XCMSnExp-class
setReplaceMethod("adjustedRtime", "MsFeatureData", function(object, value) {
    object$adjustedRtime <- value
    object
})
#' @noRd
#'
#' @rdname XCMSnExp-class
setMethod("dropAdjustedRtime", "MsFeatureData", function(object, rtraw) {
    if (hasAdjustedRtime(object)) {
        if (hasChromPeaks(object)) {
            message("Reverting retention times of identified peaks to ",
                    "original values ... ", appendLF = FALSE)
            fts <- .applyRtAdjToChromPeaks(
                chromPeaks(object), rtraw = adjustedRtime(object),
                rtadj = rtraw)
            chromPeaks(object) <- fts
            message("OK")
        }
        rm(list = "adjustedRtime", envir = object)
    }
    return(object)
})

#' @noRd
#'
#' @rdname XCMSnExp-class
setMethod("featureDefinitions", "MsFeatureData", function(object,
                                                          msLevel = integer()) {
    if (length(object$featureDefinitions)) {
        if (any(colnames(object$featureDefinitions) == "ms_level") &&
            length(msLevel))
            object$featureDefinitions[object$featureDefinitions$ms_level %in%
                                      msLevel, ]
        else object$featureDefinitions
    } else {
        warning("No aligned feature information available.", call. = FALSE)
        DataFrame()
    }
})
#' @noRd
#'
#' @rdname XCMSnExp-class
setReplaceMethod("featureDefinitions", "MsFeatureData", function(object, value) {
    object$featureDefinitions <- value
    object
})
#' @noRd
#'
#' @rdname XCMSnExp-class
setMethod("dropFeatureDefinitions", "MsFeatureData",
          function(object, dropAdjustedRtime = FALSE) {
              if (hasFeatures(object)) {
                  if (.hasFilledPeaks(object)) {
                      ## Remove filled in peaks
                      chromPeaks(object) <- chromPeaks(
                          object)[!chromPeakData(object)$is_filled, , drop = FALSE]
                      chromPeakData(object) <- extractROWS(
                          chromPeakData(object), which(!chromPeakData(object)$is_filled))
                  }
                  if (dropAdjustedRtime) {
                      ## This will ensure that the retention times of the peaks
                      ## are restored.
                      object <- dropAdjustedRtime(object)
                  }
                  rm(list = "featureDefinitions", envir = object)
              }
              object
})

#' @noRd
#'
#' @rdname XCMSnExp-class
setMethod("chromPeaks", "MsFeatureData", function(object) {
    object$chromPeaks
})
#' @noRd
#'
#' @rdname XCMSnExp-class
setReplaceMethod("chromPeaks", "MsFeatureData", function(object, value) {
    object$chromPeaks <- value
    object
})
#' @noRd
#'
#' @rdname XCMSnExp-class
setMethod("dropChromPeaks", "MsFeatureData", function(object) {
    if (hasChromPeaks(object))
        rm(list = "chromPeaks", envir = object)
    if (.has_chrom_peak_data(object))
        rm(list = "chromPeakData", envir = object)
    object
})

setMethod("chromPeakData", "MsFeatureData", function(object) {
    .chrom_peak_data(object)
})
setReplaceMethod("chromPeakData", "MsFeatureData", function(object, value) {
    object$chromPeakData <- value
    object
})
