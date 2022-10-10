#' @title Next Generation `xcms` Result Object
#'
#' @aliases XcmsExperiment-class
#'
#' @description
#'
#' The `XcmsExperiment` is a data container for `xcms` preprocessing results
#' (i.e. results from chromatographic peak detection, alignment and
#' correspondence analysis).
#'
#' It provides the same functionality than the [XCMSnExp] object, but uses the
#' more advanced and modern MS infrastructure provided by the `MsExperiment`
#' and `Spectra` Bioconductor packages. With this comes a higher flexibility on
#' how and where to store the data.
#'
#' @name XcmsExperiment
#'
#' @author Johannes Rainer
#'
#' @exportClass XcmsExperiment
#'
#' @md
NULL

setClass("XcmsExperiment",
         contains = "MsExperiment",
         slots = c(chromPeaks = "matrix",
                   chromPeakData = "data.frame"),
         prototype = prototype(chromPeaks = matrix(),
                               chromPeakData = data.frame()))

#' validation:
#' - nrow of chromPeaks and chromPeakData have to match.
#' - chromPeaks required columns; all being numeric matrix.
#' - chromPeakData required columns
#' @noRd