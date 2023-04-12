#' @title LC-MS preprocessing result test data sets
#'
#' @description
#'
#' Data sets with `xcms` preprocessing results are provided within the `xcms`
#' package and can be loaded with the `loadXcmsData` function. The available
#' Test data sets are:
#'
#' - `xdata`: an [XCMSnExp()] object with the results from a `xcms`-based
#'   pre-processing of an LC-MS untargeted metabolomics data set. The raw data
#'   files are provided in the `faahKO` R package.
#'
#' - `xmse`: an [XcmsExperiment()] object with the results from an `xcms`-based
#'   pre-processing of an LC-MS untargeted metabolomics data set (same original
#'   data set and pre-processing settings as for the `xdata` data set).
#'   The pre-processing of this data set is described in detail in the *xcms*
#'   vignette of the `xcms` package.
#'
#' Data sets can also be loaded using `data`, which would however require to
#' update objects to point to the location of the raw data files. The
#' `loadXcmsData` loads the data and ensures that all paths are updated
#' accordingly.
#'
#' @param x For `loadXcmsData`: `character(1)` with the name of the data file
#'     (object) to load.
#'
#' @name loadXcmsData
#'
#' @aliases xdata
#'
#' @aliases xmse
#'
#' @examples
#'
#' library(xcms)
#' xdata <- loadXcmsData()
loadXcmsData <- function(x = c("xmse", "xdata")) {
    x <- match.arg(x)
    e <- new.env()
    data(list = x, envir = e)
    obj <- get(x, e)
    switch(x,
           "xdata" = {
               dirname(obj) <- c(
                   rep(system.file("cdf", "KO", package = "faahKO"), 4),
                   rep(system.file("cdf", "WT", package = "faahKO"), 4))
               obj
          },
           "xmse" = {
               fls <- dir(system.file("cdf", package = "faahKO"),
                          recursive = TRUE, full.names = TRUE)
               if (!length(fls))
                   stop("Package \"faahKO\" not available. Please install ",
                        "using 'BiocManager::install(\"faahKO\")'")
               idx <- match(basename(dataStorage(spectra(obj))), basename(fls))
               if (anyNA(idx))
                   stop("Some of the original data files not found")
               obj@spectra$dataStorage <- fls[idx]
               obj
           })
}
