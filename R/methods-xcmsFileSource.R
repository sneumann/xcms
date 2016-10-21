##' @description Return the file name from an \code{xcmsFileSource} object.
##' @param object The \code{xcmsFileSource} object from which the file name
##' should be returned.
##' @param ... Not used.
##' @return A character representing the file name.
##' @noRd
setMethod("fileName", "xcmsFileSource", function(object, ...) {
    return(as.character(object))
})

