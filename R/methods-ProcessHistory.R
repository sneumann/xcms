## Methods for ProcessHistory and XProcessHistory.

setMethod("initialize", "ProcessHistory", function(.Object, ...) {
    classVersion(.Object)["ProcessHistory"] <- "0.0.2"
    callNextMethod(.Object, ...)
})
setMethod("initialize", "XProcessHistory", function(.Object, ...) {
    classVersion(.Object)["XProcessHistory"] <- "0.0.1"
    callNextMethod(.Object, ...)
})


##' @aliases param
##'
##' @description Get or set the parameter class from an \code{XProcessHistory}
##' object.
##'
##' @param object For \code{param}: the \code{XProcessHistory} object from which
##' the parameter class should be returned.
##'
##' @return A parameter object extending the \code{Param} class.
##'
##' @author Johannes Rainer
##' @noRd
setMethod("param", "XProcessHistory", function(object) {
    return(object@param)
})
##' @aliases param<-
##'
##' @param value An object extending the \code{Param} class.
##' @noRd
setReplaceMethod("param", "XProcessHistory", function(object, value) {
    object@param <- value
    if (validObject(object))
        return(object)
})
