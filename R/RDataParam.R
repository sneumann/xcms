#' @title Store `XcmsExperiment` object as .RData file
#'
#' @name RDataParam
#'
#' @export
#' 
#' @family xcms result export formats.
#' 
#' @description
#' The `RDataParam` class and method allow users to save an `XcmsExperiment`
#' object as an .RData file with a chosen filename. The object gets exported
#' using [`save()`] function.  This `param` class and method are part of the
#' possible dispatch of the generic function `storeResults`. 
#' 
#' @param fileName A parameter of `RDataParam`, that define its `fileName` slot
#' and therefore the name of the saved .RData file.The default will be
#' `tempfile()`.
#' 
#' @param param A parameter of `storeResults` defining the format in which the
#' object should be saved.
#'
#' @param object An object of class `XcmsExperiment` that will be saved as an
#' .RData file by the function `storeResults`.
#' 
#' @return for `RDataParam`: a `RDataParam` class. `storeResults` does not
#' return anything but saves the object to a *RData* file.
#'
#' @author Philippine Louail
#'
#' @examples
#'
#' ## Load a test data set with detected peaks
#' faahko_sub <- loadXcmsData("faahko_sub2")
#' 
#' ## Define param 
#' param <- RDataParam(fileName = "example_xcms_object")
#' 
#' ## Save as RData
#' storeResults(object = faahko_sub, param = param)
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

