setMethod("xcmsSource", "character", function(object) {
    if (! file.exists(object)) {
        stop("xcmsSource: file not found: ", object)
    } else if (mzR:::netCDFIsFile(object)) {
        new("netCdfSource", object)
    } else if (mzR:::rampIsFile(object)) {
        new("rampSource", object)
    } else {
        stop("xcmsSource: Could not determine file type for: ", object)
    }
})
