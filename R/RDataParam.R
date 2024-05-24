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
#' @param spectraFilePath for `loadResults` `character(1)`, defining the
#' absolute path where the spectra files should be imported from when loading
#' the object. The default will be set using the common file path of all the
#' spectra files when exporting. This is only supported if the backend of the
#' object loaded is `MsBackendMzr()`
#'
#' @inheritParams storeResults
#'
#' @inheritParams loadResults
#'
#' @importFrom Spectra dataStorageBasePath<-
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
#' xcmse <- loadResults(object = XcmsExperiment(), param = param)
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
          signature(object = "XcmsExperiment",
                    param = "RDataParam"),
          function(object, param, spectraFilePath = character()){
              env <- new.env()
              load(file = param@fileName, envir = env)
              res <- get(ls(env)[1], envir = env)
              if (!length(spectraFilePath) == 0 &&
                  inherits(spectra(res)@backend, "MsBackendMzR")) {
                  dataStorageBasePath(spectra(res)) <- spectraFilePath
              }

              res
          }
)

