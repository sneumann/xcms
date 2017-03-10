#' @include DataClasses.R functions-Chromatogram.R

setMethod("initialize", "Chromatogram", function(.Object, ...) {
    classVersion(.Object)["Chromatogram"] <- "0.0.1"
    callNextMethod(.Object, ...)
})


##' @rdname Chromatogram-class
setMethod("show", "Chromatogram", function(object) {
    cat("Object of class: ", class(object), "\n", sep = "")
    if (length(object@aggregationFun))
        cat(names(.SUPPORTED_AGG_FUN_CHROM)[.SUPPORTED_AGG_FUN_CHROM ==
                                            object@aggregationFun], "\n")
    cat("length of object: ", length(object@rtime), "\n", sep = "")
    cat("from file: ", object@fromFile, "\n", sep = "")
    cat("mz range: [", object@mz[1], ", ", object@mz[2], "]\n", sep = "")
    rtr <- range(object@rtime)
    cat("rt range: [", rtr[1], ", ", rtr[2], "]\n", sep = "")
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

## mz
##' @description \code{mz} get the mz (range) of the chromatogram. The
##' function returns a \code{numeric(2)} with the lower and upper mz value.
##'
##' @param filter For \code{mz}: whether the mz range used to filter the
##' original object should be returned (\code{filter = TRUE}), or the mz range
##' calculated on the real data (\code{filter = FALSE}).
##' 
##' @rdname Chromatogram-class
setMethod("mz", "Chromatogram", function(object, filter = FALSE) {
    if (filter)
        return(object@filterMz)
    return(object@mz)
})
## ##' @rdname Chromatogram-class
## setReplaceMethod("mz", "CentWaveParam", function(object, value) {
##     object@mzrange <- value
##     if (validObject(object))
##         return(object)
## })

##' @description \code{precursorMz} get the mz of the precursor ion. The
##' function returns a \code{numeric(2)} with the lower and upper mz value.
##' 
##' @rdname Chromatogram-class
setMethod("precursorMz", "Chromatogram", function(object) {
    return(object@precursorMz)
})

##' @aliases productMz
##' @description \code{productMz} get the mz of the product chromatogram/ion. The
##' function returns a \code{numeric(2)} with the lower and upper mz value.
##' 
##' @rdname Chromatogram-class
setMethod("productMz", "Chromatogram", function(object) {
    return(object@productMz)
})

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
