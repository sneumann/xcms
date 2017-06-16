## Methods related to xcmsSource
#' @include AllGenerics.R DataClasses.R

setMethod("xcmsSource", "character", function(object, filetype = character()) {
    if (!file.exists(object))
        stop("xcmsSource: file not found: ", object)
    ## For now we're checking the backend/file by the file name or manually
    ## defined file type (issue #188).
    if (length(filetype) == 1)
        backend <- .mzRBackendFromFiletype(filetype)
    else
        backend <- .mzRBackendFromFilename(object)
    ## issue #188
    if (is.na(backend))
        backend <- .mzRBackendFromFilecontent(object)
    ## Throw error if still NA
    if (is.na(backend))
        stop("xcmsSource: Could not determine file type for: ", object)
    if (backend == "pwiz") {
        return(new("pwizSource", object))
    } else if (backend == "netCDF") {
        return(new("netCdfSource", object))
    } else if (backend == "Ramp") {
        return(new("rampSource", object))
    }
})
