#' @include DataClasses.R do_findChromPeaks-functions.R

#' @title centWave-based peak detection in purely chromatographic data
#'
#' @description
#'
#' `findChromPeaks` on a [Chromatogram] or [Chromatograms] object with a
#' [CentWaveParam] parameter object performs centWave-based peak detection
#' on purely chromatographic data. See [centWave] for details on the method
#' and [CentWaveParam] for details on the parameter class.
#' Note that not all settings from the `CentWaveParam` will be used.
#' See [peaksWithCentWave()] for the arguments used for peak detection
#' on purely chromatographic data.
#'
#' @param object a [Chromatogram] or [Chromatograms] object.
#'
#' @param param a [CentWaveParam] object specifying the settings for the
#'     peak detection. See [peaksWithCentWave()] for the description of
#'     arguments used for peak detection.
#'
#' @param ... currently ignored.
#' 
#' @return
#'
#' If called on a `Chromatogram` object, the method returns a `matrix` with
#' the identified peaks. See [peaksWithCentWave()] for details on the matrix
#' content.
#' 
#' @seealso [peaksWithCentWave()] for the downstream function and [centWave]
#'     for details on the method.
#'
#' @author Johannes Rainer
#'
#' @rdname findChromPeaks-Chromatogram-CentWaveParam
#' 
#' @md
#'
#' @examples
#' 
#' od <- readMSData(system.file("cdf/KO/ko15.CDF", package = "faahKO"),
#'     mode = "onDisk")
#'
#' ## Extract chromatographic data for a small m/z range
#' chr <- chromatogram(od, mz = c(272.1, 272.3))[1, 1]
#'
#' ## Identify peaks with default settings
#' pks <- findChromPeaks(chr, CentWaveParam())
#' pks
#'
#' ## Plot the identified peaks
#' plot(chr, type = "h")
#' rect(xleft = pks[, "rtmin"], xright = pks[, "rtmax"],
#'     ybottom = rep(0, nrow(pks)), ytop = pks[, "maxo"], col = "#ff000020")
#'
#' ## Modify the settings
#' cwp <- CentWaveParam(snthresh = 5, peakwidth = c(10, 60))
#' pks <- findChromPeaks(chr, cwp)
#' pks
#'
#' plot(chr, type = "h")
#' rect(xleft = pks[, "rtmin"], xright = pks[, "rtmax"],
#'     ybottom = rep(0, nrow(pks)), ytop = pks[, "maxo"], col = "#00ff0020")
setMethod("findChromPeaks", signature(object = "Chromatogram",
                                      param = "CentWaveParam"),
          function(object, param, ...) {
              do.call("peaksWithCentWave",
                      args = c(list(int = intensity(object),
                                    rt = rtime(object)),
                               as(param, "list")))
          })

#' @title matchedFilter-based peak detection in purely chromatographic data
#'
#' @description
#'
#' `findChromPeaks` on a [Chromatogram] or [Chromatograms] object with a
#' [MatchedFilterParam] parameter object performs matchedFilter-based peak
#' detection on purely chromatographic data. See [matchedFilter] for details
#' on the method and [MatchedFilterParam] for details on the parameter class.
#' Note that not all settings from the `MatchedFilterParam` will be used.
#' See [peaksWithMatchedFilter()] for the arguments used for peak detection
#' on purely chromatographic data.
#'
#' @param object a [Chromatogram] or [Chromatograms] object.
#'
#' @param param a [MatchedFilterParam] object specifying the settings for the
#'     peak detection. See [peaksWithMatchedFilter()] for the description of
#'     arguments used for peak detection.
#'
#' @param ... currently ignored.
#' 
#' @return
#'
#' If called on a `Chromatogram` object, the method returns a `matrix` with
#' the identified peaks. See [peaksWithMatchedFilter()] for details on the
#' matrix content.
#' 
#' @seealso [peaksWithMatchedFilter()] for the downstream function and
#'     [matchedFilter] for details on the method.
#'
#' @author Johannes Rainer
#'
#' @rdname findChromPeaks-Chromatogram-MatchedFilter
#' 
#' @md
#'
#' @examples
#' 
#' od <- readMSData(system.file("cdf/KO/ko15.CDF", package = "faahKO"),
#'     mode = "onDisk")
#'
#' ## Extract chromatographic data for a small m/z range
#' chr <- chromatogram(od, mz = c(272.1, 272.3))[1, 1]
#'
#' ## Identify peaks with default settings
#' pks <- findChromPeaks(chr, MatchedFilterParam())
#' pks
#'
#' ## Plot the identified peaks
#' plot(chr, type = "h")
#' rect(xleft = pks[, "rtmin"], xright = pks[, "rtmax"],
#'     ybottom = rep(0, nrow(pks)), ytop = pks[, "maxo"], col = "#ff000020")
#'
#' ## Modify the settings
#' mfp <- MatchedFilterParam(fwhm = 60)
#' pks <- findChromPeaks(chr, mfp)
#' pks
#'
#' plot(chr, type = "h")
#' rect(xleft = pks[, "rtmin"], xright = pks[, "rtmax"],
#'     ybottom = rep(0, nrow(pks)), ytop = pks[, "maxo"], col = "#00ff0020")
setMethod("findChromPeaks", signature(object = "Chromatogram",
                                      param = "MatchedFilterParam"),
          function(object, param, ...) {
              do.call("peaksWithMatchedFilter",
                      args = c(list(int = intensity(object),
                                    rt = rtime(object)),
                               as(param, "list")))
          })
