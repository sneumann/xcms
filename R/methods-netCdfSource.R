## Methods for netCdfSource class.
#' @include AllGenerics.R DataClasses.R

setMethod("initialize", "netCdfSource", function(.Object, path) {
    .Object@cdf <- mzR:::netCDFOpen(path)
    .Object <- setDataPart(.Object, path)
    validObject(.Object)
    .Object
})

setMethod("loadRaw", "netCdfSource", function(object, includeMSn) {
    if (!missing(includeMSn) && includeMSn) {
        warning("Reading of MSn spectra for NetCDF not supported")
    }

    on.exit(mzR:::netCDFClose(object@cdf))
    mzR:::netCDFRawData(object@cdf)
})

