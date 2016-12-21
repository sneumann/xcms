## Methods for the MsFeatureData class.
#' @include functions-MsFeatureData.R

setMethod("initialize", "MsFeatureData", function(.Object, ...) {
    classVersion(.Object)["MsFeatureData"] <- "0.0.1"
    callNextMethod(.Object, ...)
})

setMethod("show", "MsFeatureData", function(object) {
    cat("Object of class: ", class(object), "\n")
    ks <- ls(object)
    if (length(ks)) {
        cat(" Available feature data:\n")
        for (k in ks)
            cat(" -", k, "\n")
    }
})

setValidity("MsFeatureData", function(object) {
    return(validateMsFeatureData(object))
})

## features: getter and setter for the features matrix.
## featureGroups: getter and setter for the featureGroups DataFrame.
## adjustedRtime: getter and setter for the adjustedRtime list.
