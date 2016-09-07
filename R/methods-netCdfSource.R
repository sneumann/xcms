## Methods for netCdfSource class.
#' @include AllGenerics.R DataClasses.R

setMethod("initialize", "netCdfSource", function(.Object, path) {
    .Object@cdf <- mzR:::netCDFOpen(path)
    .Object <- setDataPart(.Object, path)
    validObject(.Object)
    gc()
    .Object
})

setMethod("loadRaw", "netCdfSource", function(object, includeMSn) {
    if (!missing(includeMSn) && includeMSn) {
        warning("Reading of MSn spectra for NetCDF not supported")
    }
    on.exit(mzR:::netCDFClose(object@cdf))
    on.exit(gc(), add = TRUE)
    res <- mzR:::netCDFRawData(object@cdf)
    return(res)
})

