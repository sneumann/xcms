## Methods for the XCMSnExp object representing untargeted metabolomics
## results

setMethod("initialize", "XCMSnExp", function(.Object, ...) {
    classVersion(.Object)["XCMSnExp"] <- "0.0.1"
    callNextMethod(.Object, ...)
})

##' @rdname XCMSnExp-class
setMethod("show", "XCMSnExp", function(object) {
    cat("Object of class: ", class(object), "\n")
    callNextMethod()
})

##' @aliases hasAdjustedRtime
##'
##' @description \code{hasAdjustedRtime}: whether the object provides adjusted
##' retention times.
##'
##' @rdname XCMSnExp-class
setMethod("hasAdjustedRtime", "XCMSnExp", function(object) {
    return(hasAdjustedRtime(object@msFeatureData))
})

##' @aliases hasAlignedFeatures
##'
##' @description \code{hasAlignedFeatures}: whether the object contains feature
##' alignment results.
##'
##' @rdname XCMSnExp-class
setMethod("hasAlignedFeatures", "XCMSnExp", function(object) {
    return(hasAlignedFeatures(object@msFeatureData))
})

##' @aliases hasDetectedFeatures
##'
##' @description \code{hasDetectedFeatures}: whether the object contains feature
##' detection results.
##'
##' @rdname XCMSnExp-class
setMethod("hasDetectedFeatures", "XCMSnExp", function(object) {
    return(hasDetectedFeatures(object@msFeatureData))
})

##' @aliases adjustedRtime
##'
##' @description The \code{adjustedRtime},\code{adjustedRtime<-} method
##' extract/set adjusted retention times. Retention times are adjusted by
##' retention time correction/adjustment methods.
##'
##' @return For \code{adjustedRtime}: a \code{list} (length equal to the number
##' of samples) with numeric vectors with the adjusted retention time for each
##' scan or \code{NULL} if no such data is present.
##'
##' @rdname XCMSnExp-class
setMethod("adjustedRtime", "XCMSnExp", function(object) {
    return(adjustedRtime(object@msFeatureData))
})
##' @aliases adjustedRtime<-
##'
##' @rdname XCMSnExp-class
setReplaceMethod("adjustedRtime", "XCMSnExp", function(object, value) {
    newFd <- new("MsFeatureData")
    newFd@.xData <- .copy_env(object@msFeatureData)
    adjustedRtime(newFd) <- value
    object@msFeatureData <- newFd
    if (validObject(object)) {
        ## Lock the environment so that only accessor methods can change values.
        lockEnvironment(newFd, bindings = TRUE)
        return(object)
    }
})

##' @aliases featureGroups
##'
##' @description The \code{featureGroups}, \code{featureGroups<-} methods extract
##' or set the feature alignment results.
##'
##' @return For \code{featureGroups}: a \code{DataFrame} with feature alignment
##' information, each row corresponding to one group of aligned features (across
##' samples) and columns \code{"mzmed"} (median mz value), \code{"mzmin"}
##' (minimal mz value), \code{"mzmax"} (maximum mz value), \code{"rtmed"} (median
##' retention time), \code{"rtmin"} (minimal retention time), \code{"rtmax"}
##' (maximal retention time) and \code{"featureidx"}. Column \code{"featureidx"}
##' contains a \code{list} with indices of features (rows) in the matrix returned
##' by the \code{features} method that belong to that feature group. The method
##' returns \code{NULL} if no aligned feature information is present.
##'
##' @rdname XCMSnExp-class
setMethod("featureGroups", "XCMSnExp", function(object) {
    return(featureGroups(object@msFeatureData))
})
##' @aliases featureGroups<-
##'
##' @rdname XCMSnExp-class
setReplaceMethod("featureGroups", "XCMSnExp", function(object, value) {
    newFd <- new("MsFeatureData")
    newFd@.xData <- .copy_env(object@msFeatureData)
    featureGroups(newFd) <- value
    object@msFeatureData <- newFd
    if (validObject(object)) {
        ## Lock the environment so that only accessor methods can change values.
        lockEnvironment(newFd, bindings = TRUE)
        return(object)
    }
})

##' @aliases features
##'
##' @description The \code{features}, \code{features<-} methods extract or set
##' the matrix containing the information on identified features. See description
##' on the return value for details on the matrix columns. Users usually don't
##' have to use the \code{features<-} method directly as detected features are
##' added to the object by the \code{\link{detectFeatures}} method.
##'
##' @return For \code{features}: a \code{matrix} with at least the following
##' columns: \code{"mz"} (mz value for the largest intensity), \code{"mzmin"} (
##' minimal mz value), \code{"mzmax"} (maximal mz value), \code{"rt"} (retention
##' time for the peak apex), \code{"rtmin"} (minimal retention time), \code{"rtmax"}
##' (maximal retention time), \code{"into"} (integrated, original, intensity of
##' the feature) and \code{"sample"} (sample index in which the feature was
##' identified). Depending on the employed feature detection algorithm and the
##' \code{verboseColumns} parameter of it additional columns might be returned.
##'
##' @rdname XCMSnExp-class
setMethod("features", "XCMSnExp", function(object) {
    return(features(object@msFeatureData))
})
##' @aliases features<-
##'
##' @rdname XCMSnExp-class
setReplaceMethod("features", "XCMSnExp", function(object, value) {
    newFd <- new("MsFeatureData")
    newFd@.xData <- .copy_env(object@msFeatureData)
    features(newFd) <- value
    object@msFeatureData <- newFd
    if (validObject(object)) {
        ## Lock the environment so that only accessor methods can change values.
        lockEnvironment(newFd, bindings = TRUE)
        return(object)
    }
})
