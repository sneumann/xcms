#' @rdname findChromPeaks-Chromatogram-CentWaveParam
#'
#' @aliases findChromPeaks-Chromatogram-CentWaveParam
#'
#' @examples
#'
#' ## Perform peak detection on an Chromatograms object
#' od3 <- readMSData(c(system.file("cdf/KO/ko15.CDF", package = "faahKO"),
#'     system.file("cdf/KO/ko16.CDF", package = "faahKO"),
#'     system.file("cdf/KO/ko18.CDF", package = "faahKO")),
#'     mode = "onDisk")
#'
#' ## Extract chromatograms for a m/z - retention time slice
#' chrs <- chromatogram(od3, mz = 344, rt = c(2500, 3500))
#'
#' ## Perform peak detection using CentWave
#' xchrs <- findChromPeaks(chrs, param = CentWaveParam())
#' xchrs
#'
#' ## Extract the identified chromatographic peaks
#' chromPeaks(xchrs)
#'
#' ## plot the result
#' plot(xchrs)
setMethod("findChromPeaks", signature(object = "Chromatograms",
                                      param = "CentWaveParam"),
          function(object, param, BPPARAM = bpparam(), ...) {
              .findChromPeaks_XChromatograms(object = object, param = param,
                                             BPPARAM = BPPARAM, ...)
          })

#' @rdname findChromPeaks-Chromatogram-CentWaveParam
setMethod("findChromPeaks", signature(object = "Chromatograms",
                                      param = "MatchedFilterParam"),
          function(object, param, BPPARAM = BPPARAM, ...) {
              .findChromPeaks_XChromatograms(object = object,
                                             param = param,
                                             BPPARAM = BPPARAM, ...)
          })

.findChromPeaks_XChromatograms <- function(object, param, BPPARAM, ...) {
    startDate <- date()
    if (missing(BPPARAM))
        BPPARAM <- bpparam()
    object <- as(object, "XChromatograms")
    object@.Data <- matrix(bplapply(c(object@.Data), FUN = findChromPeaks,
                                    param = param, BPPARAM = BPPARAM),
                           ncol = ncol(object),
                           dimnames = dimnames(object@.Data))
    ph_len <- length(object@.processHistory)
    if (ph_len && processType(object@.processHistory[[ph_len]]) ==
        .PROCSTEP.PEAK.DETECTION)
        object@.processHistory <-
            object@.processHistory[seq_len(ph_len - 1)]
    object <- addProcessHistory(
        object, XProcessHistory(param = param, date. = startDate,
                                type. = .PROCSTEP.PEAK.DETECTION,
                                fileIndex = seq_len(ncol(object))))
    if (validObject(object)) object
}

#' @rdname align-Chromatogram
setMethod("align", signature = c(x = "Chromatograms", y = "Chromatogram"),
          function(x, y, method = c("matchRtime", "approx"), ...) {
              x@.Data <- matrix(lapply(x@.Data, .align_chromatogram, y = y,
                                       method = method, ...),
                                nrow = nrow(x), dimnames = dimnames(x))
              validObject(x)
              x
          })
