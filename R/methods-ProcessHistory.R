## Methods for ProcessHistory and XProcessHistory.
#' @include functions-ProcessHistory.R

#' @rdname ProcessHistory-class
setMethod("show", "ProcessHistory", function(object) {
    cat("Object of class \"", class(object), "\"\n", sep = "")
    cat(" type:", object@type, "\n")
    cat(" date:", object@date, "\n")
    cat(" info:", object@info, "\n")
    cat(" fileIndex:", paste0(object@fileIndex, collapse = ","), "\n")
})
#' @rdname ProcessHistory-class
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
#' @description \code{processParam}, \code{processParam<-}: get or set the
#'     parameter class from an \code{XProcessHistory} object.
#'
#' @param object A \code{ProcessHistory} or \code{XProcessHistory} object.
#'
#' @return For \code{processParam}: a parameter object extending the
#'     \code{Param} class.
#'
#' @author Johannes Rainer
#'
#' @rdname ProcessHistory-class
setMethod("processParam", "XProcessHistory", function(object) {
    return(object@param)
})
#' @aliases processParam<-
#'
#' @param value An object extending the \code{Param} class.
#'
#' @noRd
setReplaceMethod("processParam", "XProcessHistory", function(object, value) {
    object@param <- value
    if (validObject(object))
        return(object)
})
#' @description \code{msLevel}: returns the MS level on which a certain analysis
#'     has been performed, or \code{NA} if not defined.
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
#' @description The \code{processType} method returns a character specifying the
#'     processing step \emph{type}.
#'
#' @return The \code{processType} method returns a character string with the
#'     processing step type.
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
#' @description The \code{processDate} extracts the start date of the processing
#'     step.
#'
#' @return The \code{processDate} method returns a character string with the
#'     time stamp of the processing step start.
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
#' @description The \code{processInfo} extracts optional additional information
#'     on the processing step.
#'
#' @return The \code{processInfo} method returns a character string with
#'     optional additional informations.
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
#' @description The \code{fileIndex} extracts the indices of the files on which
#'     the processing step was applied.
#'
#' @return The \code{fileIndex} method returns a integer vector with the index
#'     of the files/samples on which the processing step was applied.
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
