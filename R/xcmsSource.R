setMethod("xcmsSource", "character", function(object) {
    if (useOriginalCode()) {
        if (! file.exists(object)) {
            stop("xcmsSource: file not found: ", object)
        } else if (mzR:::netCDFIsFile(object)) {
            new("netCdfSource", object)
        } else if (mzR:::rampIsFile(object)) {
            new("rampSource", object)
        } else {
            stop("xcmsSource: Could not determine file type for: ", object)
        }
    } else {
        if (! file.exists(object)) {
            stop("xcmsSource: file not found: ", object)
        } else if (isCdfFile(object)) {
            new("netCdfSource", object)
        } else if (isMzMLFile(object)) {
            new("rampSource", object)
        } else {
            stop("xcmsSource: Could not determine file type for: ", object)
        }
    }
})
