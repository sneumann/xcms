##' @description Return the file name from an \code{xcmsFileSource} object.
##' @param object The \code{xcmsFileSource} object from which the file name
##' should be returned.
##' @param ... Not used.
##' @return A character representing the file name.
##' @noRd
setMethod("fileName", "xcmsFileSource", function(object, ...) {
    return(as.character(object))
})

##' Reads the MS1 and/or MSn data from an netCDF or mzML file.
##' @param object The \code{xcmsFileSource} object representing the file name.
##' @param includeMSn Whether MSn data should be imported too.
##' @return Returns a \code{list} with rt, tic, scanindex, mz and intensity.
##' @noRd
setMethod("loadRaw", "xcmsFileSource", function(object, includeMSn = FALSE) {
    return(readRawData(fileName(object), includeMSn = includeMSn))
})
