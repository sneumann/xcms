## Methods for the Param class and sub-classes
#' @include functions-Params.R

############################################################
## Param
setMethod("as.list", signature(x = "Param"), function(x, ...) {
    return(.param2list(x))
})
## The 'setAs' method.
setAs("Param" ,"list", function(from){
    return(.param2list(from))
})
setMethod("initialize", "Param", function(.Object, ...) {
    classVersion(.Object)["Param"] <- "0.0.1"
    callNextMethod(.Object, ...)
})


############################################################
## CentWaveParam
setMethod("initialize", "CentWaveParam", function(.Object, ...) {
    classVersion(.Object)["CentWaveParam"] <- "0.0.1"
    callNextMethod(.Object, ...)
})
