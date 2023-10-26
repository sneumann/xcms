#' @title Store `XcmsExperiment` object as .RData file
#'
#' @name RDataParam
#'
#' @export
#'
#' @description
#' The `RDataParam` class and method allow users to save an `XcmsExperiment`
#' object as an .RData file with a chosen filename. This new `param` class and
#' method are part of the possible dispatch of the generic function
#' `storeResults`. 
#' 
#' Other available `param` classes and linked methods include:
#'
#' - `PlainTextParam`: ...
#'
#' - `MzTabMParam`: ...
#' 
#' @param param A parameter defining the format in which the object should be
#' saved. For the `RDataParam`, it has one slot for the `fileName` in which the
#' object is going to be saved. The default will be `tempfile()`.
#'
#' @param object An object of class `XcmsExperiment` that will be saved as an
#' .RData file.
#' 
#' @return The saved object as an .RData file. 
#'
#' @author Philippine Louail
#'
#' @examples
#'
#'## Get an `XcmsExperiment` object
#' x <- ...
#' 
#' ## Define param 
#' param <- RDataParam(fileName = "example_xcms_object")
#' 
#' ## Save as RData
#' storeResults(object = x, param = param)
#' 
NULL

#' @noRd
setClass("RDataParam",
         slots = c(fileName = "character"), 
         contains = "Param",
         prototype = prototype(
             fileName = character()),
         validity = function(object) {
             msg <- NULL
             if (length(object@fileName) != 1)  
                 msg <- c("'fileName' has to be a character string of length 1")
             msg
         })

#' @rdname RDataParam
#'
#' @export
RDataParam <- function(fileName = tempfile()) {
    new("RDataParam", fileName = fileName)
}
 
#' @rdname RDataParam
setMethod("storeResults",
          signature(object = "XcmsExperiment", 
                    param = "RDataParam"), 
          function(object, param){
              save(object, file = param@fileName)
              }
          )
