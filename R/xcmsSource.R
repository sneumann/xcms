setMethod("xcmsSource", "character", function(object) {
    if (mzR:::netCDFIsFile(object)) {
        new("netCdfSource", object)
    } else if (mzR:::rampIsFile(object)) {
        new("rampSource", object)
    } else {
        stop("Could not determine file type")
    }
})
