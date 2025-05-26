## Methods for ProcessHistory and XProcessHistory.
#' @include functions-ProcessHistory.R

#' @rdname hidden_aliases
setMethod("show", "ProcessHistory", function(object) {
    cat("Object of class \"", class(object), "\"\n", sep = "")
    cat(" type:", object@type, "\n")
    cat(" date:", object@date, "\n")
    cat(" info:", object@info, "\n")
    cat(" fileIndex:", paste0(object@fileIndex, collapse = ","), "\n")
})
#' @rdname hidden_aliases
setMethod("show", "XProcessHistory", function(object) {
    callNextMethod()
    pcLabel <- "-none-"
    if (length(object@param))
        pcLabel <- class(object@param)
    cat(" Parameter class:", pcLabel, "\n")
    if (!is.na(msLevel(object)))
        cat(" MS level(s)", paste(msLevel(object), sep = " "), "\n")
})

#' @aliases processParam
#'
#' @description `processParam()`, `processParam<-`: get or set the
#'     parameter class from an `XProcessHistory` object.
#'
#' @param object A `ProcessHistory` or `XProcessHistory` object.
#'
#' @return For `processParam`: a parameter object extending the
#'     `Param` class.
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @rdname ProcessHistory-class
setMethod("processParam", "XProcessHistory", function(object) {
    return(object@param)
})
#' @aliases processParam<-
#'
#' @param value An object extending the `Param` class.
#'
#' @md
#'
#' @noRd
setReplaceMethod("processParam", "XProcessHistory", function(object, value) {
    object@param <- value
    if (validObject(object))
        return(object)
})
#' @description `msLevel()`: returns the MS level on which a certain analysis
#'     has been performed, or `NA` if not defined.
#'
#' @md
#'
#' @rdname ProcessHistory-class
setMethod("msLevel", "XProcessHistory", function(object) {
    if (.hasSlot(object, "msLevel"))
        return(object@msLevel)
    NA_integer_
})


## Methods:
#' @aliases processType
#'
#' @description The `processType()` method returns a character specifying the
#'     processing step *type*.
#'
#' @return The `processType()` method returns a character string with the
#'     processing step type.
#'
#' @md
#'
#' @rdname ProcessHistory-class
setMethod("processType", "ProcessHistory", function(object) {
    return(object@type)
})
#' @noRd
setReplaceMethod("processType", "ProcessHistory", function(object, value) {
    object@type <- value
    if (validObject(object))
        return(object)
})

#' @aliases processDate
#'
#' @description The `processDate()` extracts the start date of the processing
#'     step.
#'
#' @return The `processDate()` method returns a character string with the
#'     time stamp of the processing step start.
#'
#' @md
#'
#' @rdname ProcessHistory-class
setMethod("processDate", "ProcessHistory", function(object) {
    return(object@date)
})
#' @noRd
setReplaceMethod("processDate", "ProcessHistory", function(object, value) {
    object@date <- value
    if (validObject(object))
        return(object)
})

#' @aliases processInfo
#'
#' @description The `processInfo()` extracts optional additional information
#'     on the processing step.
#'
#' @return The `processInfo()` method returns a character string with
#'     optional additional informations.
#'
#' @md
#'
#' @rdname ProcessHistory-class
setMethod("processInfo", "ProcessHistory", function(object) {
    return(object@info)
})
#' @noRd
setReplaceMethod("processInfo", "ProcessHistory", function(object, value) {
    object@info <- value
    if (validObject(object))
        return(object)
})

#' @aliases fileIndex
#'
#' @description The `fileIndex()` extracts the indices of the files on which
#'     the processing step was applied.
#'
#' @return The `fileIndex()` method returns a integer vector with the index
#'     of the files/samples on which the processing step was applied.
#'
#' @md
#'
#' @rdname ProcessHistory-class
setMethod("fileIndex", "ProcessHistory", function(object) {
    return(object@fileIndex)
})
#' @noRd
setReplaceMethod("fileIndex", "ProcessHistory", function(object, value) {
    object@fileIndex <- as.integer(value)
    if (validObject(object))
        return(object)
})
