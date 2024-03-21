#' @title Store/Load an `XcmsExperiment` object as/from .RData file
#'
#' @name RDataParam
#'
#' @export
#'
#' @family xcms result export and import formats.
#'
#' @description
#' The `RDataParam` class and method allow users to save or load an
#' `XcmsExperiment` object as/from an .RData file with a defined filename. The
#' object gets exported using [`save()`] function and imported using the
#' [`load()`] function. This `param` class and method are part of the
#' possible dispatch of the generic functions `storeResults()` and
#' `loadResults()`.
#'
#' @param fileName for `RDataParam` `character(1)`, defining the file name. The
#' default will be `tempfile()`.
#'
#' @inheritParams storeResults
#'
#' @return for `RDataParam`: a `RDataParam` class. `storeResults` does not
#' return anything but saves the object to a RData file. `loadResults` returns
#' an object of class `XcmsExperiment`
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
#' ## Load this saved dataset
#' xcmse <- loadResults(param = param)
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

#' @rdname RDataParam
setMethod("loadResults",
          signature(param = "RDataParam"),
          function(param){
              load(file = param@fileName)
          }
)

