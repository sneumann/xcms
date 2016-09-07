## Methods for rampSource
#' @include AllGenerics.R DataClasses.R

setMethod("initialize", "rampSource", function(.Object, path) {
    .Object@rampid <- mzR:::rampOpen(path)
    .Object <- setDataPart(.Object, path)
    validObject(.Object)
    gc()
    .Object
})

setMethod("loadRaw", "rampSource", function(object, includeMSn) {
    on.exit(mzR:::rampClose(object@rampid))
    on.exit(gc(), add = TRUE)
    rawdata <- mzR:::rampRawData(object@rampid)
    if (!missing(includeMSn) && includeMSn) {
        rawdata$MSn <- mzR:::rampRawDataMSn(object@rampid)
    }
    rawdata
})

