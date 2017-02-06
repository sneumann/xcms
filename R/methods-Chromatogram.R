#' @include DataClasses.R functions-Chromatogram.R

##' @rdname Chromatogram-class
setMethod("show", "Chromatogram", function(object) {
    cat("Object of class: ", class(object), "\n", sep = "")
    if (length(object@aggregationFun))
        cat(names(.SUPPORTED_AGG_FUN_CHROM)[.SUPPORTED_AGG_FUN_CHROM ==
                                            object@aggregationFun], "\n")
    cat("length of object: ", length(object@rtime), "\n", sep = "")
    cat("mzrange: [", object@mzrange[1], ", ", object@mzrange[2], "]\n", sep = "")
})

## Methods:

## rtime
##' @description \code{rtime} returns the retention times for the rentention time
##' - intensity pairs stored in the chromatogram.
##'
##' @param object A \code{Chromatogram} object.
##' 
##' @rdname Chromatogram-class
setMethod("rtime", "Chromatogram", function(object) {
    return(object@rtime)
})

## intensity
##' @description \code{intensity} returns the intensity for the rentention time
##' - intensity pairs stored in the chromatogram.
##' 
##' @rdname Chromatogram-class
setMethod("intensity", "Chromatogram", function(object) {
    return(object@intensity)
})

## mzrange
##' @description \code{mzrange} get or set the mz range of the
##' chromatogram.
##' 
##' @rdname Chromatogram-class
setMethod("mzrange", "Chromatogram", function(object) {
    return(object@mzrange)
})
## ##' @rdname Chromatogram-class
## setReplaceMethod("mzrange", "CentWaveParam", function(object, value) {
##     object@mzrange <- value
##     if (validObject(object))
##         return(object)
## })

## aggregationFun
##' @aliases aggregationFun
##' @description \code{aggregationFun,aggregationFun<-} get or set the
##' aggregation function.
##' 
##' @rdname Chromatogram-class
setMethod("aggregationFun", "Chromatogram", function(object) {
    return(object@aggregationFun)
})
## ##' @rdname Chromatogram-class
## setReplaceMethod("aggregationFun", "CentWaveParam", function(object, value) {
##     object@aggregationFun <- value
##     if (validObject(object))
##         return(object)
## })

## fromFile
##' @description \code{fromFile} returns the value from the \code{fromFile} slot.
##' 
##' @rdname Chromatogram-class
setMethod("fromFile", "Chromatogram", function(object) {
    return(object@fromFile)
})

## length
##' @description \code{length} returns the length (number of retention time -
##' intensity pairs) of the chromatogram.
##' 
##' @param x For \code{as.data.frame} and \code{length}: a \code{Chromatogram}
##' object.
##'
##' @rdname Chromatogram-class
setMethod("length", "Chromatogram", function(x) {
    return(length(x@rtime))
})

## as.data.frame
##' @description \code{as.data.frame} returns the \code{rtime} and
##' \code{intensity} values from the object as \code{data.frame}.
##' 
##' @rdname Chromatogram-class
setMethod("as.data.frame", "Chromatogram", function(x) {
    return(data.frame(rtime = x@rtime, intensity = x@intensity))
})
