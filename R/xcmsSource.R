## Methods related to xcmsSource
#' @include AllGenerics.R DataClasses.R

setMethod("xcmsSource", "character", function(object) {
    if (!file.exists(object)) {
        stop("xcmsSource: file not found: ", object)
    } else if (isCdfFile(object)) {
        new("netCdfSource", object)
    } else if (isMzMLFile(object)) {
        new("pwizSource", object)
    } else if (isRampFile(object)) {
        new("rampSource", object)
    } else {
        stop("xcmsSource: Could not determine file type for: ", object)
    }
})
