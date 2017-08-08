## Methods for the MsFeatureData class.
#' @include functions-MsFeatureData.R do_adjustRtime-functions.R

setMethod("initialize", "MsFeatureData", function(.Object, ...) {
    classVersion(.Object)["MsFeatureData"] <- "0.0.1"
    callNextMethod(.Object, ...)
})

setValidity("MsFeatureData", function(object) {
    return(validateMsFeatureData(object))
})

##' @rdname XCMSnExp-class
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


#' @rdname XCMSnExp-class
setMethod("hasAdjustedRtime", "MsFeatureData", function(object) {
    !is.null(object$adjustedRtime)
})

#' @rdname XCMSnExp-class
setMethod("hasFeatures", "MsFeatureData", function(object) {
    !is.null(object$featureDefinitions)
})

#' @rdname XCMSnExp-class
setMethod("hasChromPeaks", "MsFeatureData", function(object) {
    !is.null(object$chromPeaks)
})

#' @rdname XCMSnExp-class
setMethod("adjustedRtime", "MsFeatureData", function(object) {
    if (hasAdjustedRtime(object))
        return(object$adjustedRtime)
    warning("No adjusted retention times available.")
    return(NULL)
})
#' @rdname XCMSnExp-class
setReplaceMethod("adjustedRtime", "MsFeatureData", function(object, value) {
    object$adjustedRtime <- value
    if (validObject(object))
        return(object)
})
#' @rdname XCMSnExp-class
setMethod("dropAdjustedRtime", "MsFeatureData", function(object) {
    if (hasAdjustedRtime(object)) {
        rm(list = "adjustedRtime", envir = object)
    }
    return(object)
})

#' @rdname XCMSnExp-class
setMethod("featureDefinitions", "MsFeatureData", function(object) {
    if (hasFeatures(object))
        return(object$featureDefinitions)
    warning("No aligned feature information available.")
    return(NULL)
})
#' @rdname XCMSnExp-class
setReplaceMethod("featureDefinitions", "MsFeatureData", function(object, value) {
    object$featureDefinitions <- value
    if (validObject(object))
        return(object)
})
#' @rdname XCMSnExp-class
setMethod("dropFeatureDefinitions", "MsFeatureData", function(object) {
    if (hasFeatures(object))
        rm(list = "featureDefinitions", envir = object)
    return(object)
})

#' @rdname XCMSnExp-class
setMethod("chromPeaks", "MsFeatureData", function(object) {
    if (hasChromPeaks(object))
        return(object$chromPeaks)
    warning("No chromatographic peaks available.")
    return(NULL)
})
#' @rdname XCMSnExp-class
setReplaceMethod("chromPeaks", "MsFeatureData", function(object, value) {
    object$chromPeaks <- value
    if (validObject(object))
        return(object)
})
#' @rdname XCMSnExp-class
setMethod("dropChromPeaks", "MsFeatureData", function(object) {
    if (hasChromPeaks(object))
        rm(list = "chromPeaks", envir = object)
    return(object)
})
