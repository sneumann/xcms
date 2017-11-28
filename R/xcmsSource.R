## Methods related to xcmsSource
#' @include AllGenerics.R DataClasses.R

setMethod("xcmsSource", "character", function(object, filetype = character()) {
    if (!file.exists(object))
        stop("xcmsSource: file not found: ", object)
    return(new("xcmsFileSource", object))
})
