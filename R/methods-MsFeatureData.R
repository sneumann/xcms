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

## features: getter and setter for the features matrix.
## featureGroups: getter and setter for the featureGroups DataFrame.
## adjustedRtime: getter and setter for the adjustedRtime list.


##' @rdname XCMSnExp-class
setMethod("hasAdjustedRtime", "MsFeatureData", function(object) {
    return(any(ls(object) == "adjustedRtime"))
})

##' @rdname XCMSnExp-class
setMethod("hasAlignedFeatures", "MsFeatureData", function(object) {
    return(any(ls(object) == "featureGroups"))
})

##' @rdname XCMSnExp-class
setMethod("hasDetectedFeatures", "MsFeatureData", function(object) {
    return(any(ls(object) == "features"))
})

##' @rdname XCMSnExp-class
setMethod("adjustedRtime", "MsFeatureData", function(object) {
    if (hasAdjustedRtime(object))
        return(object$adjustedRtime)
    warning("No adjusted retention times available.")
    return(NULL)
})
##' @rdname XCMSnExp-class
setReplaceMethod("adjustedRtime", "MsFeatureData", function(object, value) {
    object$adjustedRtime <- value
    if (validObject(object))
        return(object)
})
##' @rdname XCMSnExp-class
setMethod("dropAdjustedRtime", "MsFeatureData", function(object) {
    if (hasAdjustedRtime(object)) {
        rm(list = "adjustedRtime", envir = object)
    }
    return(object)
})

##' @rdname XCMSnExp-class
setMethod("featureGroups", "MsFeatureData", function(object) {
    if (hasAlignedFeatures(object))
        return(object$featureGroups)
    warning("No aligned feature information available.")
    return(NULL)
})
##' @rdname XCMSnExp-class
setReplaceMethod("featureGroups", "MsFeatureData", function(object, value) {
    object$featureGroups <- value
    if (validObject(object))
        return(object)
})
##' @rdname XCMSnExp-class
setMethod("dropFeatureGroups", "MsFeatureData", function(object) {
    if (hasAlignedFeatures(object))
        rm(list = "featureGroups", envir = object)
    return(object)
})

##' @rdname XCMSnExp-class
setMethod("features", "MsFeatureData", function(object) {
    if (hasDetectedFeatures(object))
        return(object$features)
    warning("No detected features available.")
    return(NULL)
})
##' @rdname XCMSnExp-class
setReplaceMethod("features", "MsFeatureData", function(object, value) {
    object$features <- value
    if (validObject(object))
        return(object)
})
##' @rdname XCMSnExp-class
setMethod("dropFeatures", "MsFeatureData", function(object) {
    if (hasDetectedFeatures(object))
        rm(list = "features", envir = object)
    return(object)
})
