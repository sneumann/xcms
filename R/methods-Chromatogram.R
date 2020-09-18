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
#' After chromatographic peak detection, identified peaks can also be *refined*
#' with the [refineChromPeaks()] method, which can help to reduce peak
#' detection artifacts.
#'
#' @param object a [Chromatogram] or [Chromatograms] object.
#'
#' @param param a [CentWaveParam] object specifying the settings for the
#'     peak detection. See [peaksWithCentWave()] for the description of
#'     arguments used for peak detection.
#'
#' @param BPPARAM a parameter class specifying if and how parallel processing
#'     should be performed (only for `XChromatograms` objects). It defaults to
#'     `bpparam()`. See [bpparam()] for more information.
#'
#' @param ... currently ignored.
#'
#' @return
#'
#' If called on a `Chromatogram` object, the method returns an [XChromatogram]
#' object with the identified peaks. See [peaksWithCentWave()] for details on
#' the peak matrix content.
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
#' ## Disabling parallel processing in this example
#' register(SerialParam())
#'
#' ## Extract chromatographic data for a small m/z range
#' chr <- chromatogram(od, mz = c(272.1, 272.3))[1, 1]
#'
#' ## Identify peaks with default settings
#' xchr <- findChromPeaks(chr, CentWaveParam())
#' xchr
#'
#' ## Plot data and identified peaks.
#' plot(xchr)
#'
#' ## Modify the settings
#' cwp <- CentWaveParam(snthresh = 5, peakwidth = c(10, 60))
#' xchr <- findChromPeaks(chr, cwp)
#' xchr
#'
#' plot(xchr)
setMethod("findChromPeaks", signature(object = "Chromatogram",
                                      param = "CentWaveParam"),
          function(object, param, ...) {
              res <- do.call("peaksWithCentWave",
                             args = c(list(int = intensity(object),
                                           rt = rtime(object)),
                                      as(param, "list")))
              object <- as(object, "XChromatogram")
              chromPeaks(object) <- res
              object
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
#' ## Disabling parallel processing in this example
#' register(SerialParam())
#'
#' ## Extract chromatographic data for a small m/z range
#' chr <- chromatogram(od, mz = c(272.1, 272.3))[1, 1]
#'
#' ## Identify peaks with default settings
#' xchr <- findChromPeaks(chr, MatchedFilterParam())
#'
#' ## Plot the identified peaks
#' plot(xchr)
#'
#' ## Modify the settings
#' mfp <- MatchedFilterParam(fwhm = 60)
#' xchr <- findChromPeaks(chr, mfp)
#'
#' plot(xchr)
setMethod("findChromPeaks", signature(object = "Chromatogram",
                                      param = "MatchedFilterParam"),
          function(object, param, ...) {
              res <- do.call("peaksWithMatchedFilter",
                             args = c(list(int = intensity(object),
                                           rt = rtime(object)),
                                      as(param, "list")))
              object <- as(object, "XChromatogram")
              chromPeaks(object) <- res
              object
          })

#' @title Aligning chromatographic data
#'
#' @aliases align
#'
#' @rdname align-Chromatogram
#'
#' @description
#'
#' Align chromatogram `x` against chromatogram `y`. The resulting chromatogram
#' has the same length (number of data points) than `y` and the same retention
#' times thus allowing to perform any pair-wise comparisons between the
#' chromatograms. If `x` is a [Chromatograms()] object, each `Chromatogram` in
#' it is aligned against `y`.
#'
#' Parameter `method` allows to specify which alignment method
#' should be used. Currently there are the following options:
#'
#' - `method = "matchRtime"` (the default): match data points in the first
#'   chromatogram (`x`) to those of the second (`y`) based on the difference
#'   between their retention times: each data point in `x` is assigned to the
#'   data point in `y` with the smallest difference in their retention times
#'   if their difference is smaller than the minimum average difference between
#'   retention times in `x` or `y`.
#'
#' - `method = "approx"`: uses the base R `approx` function to approximate
#'   intensities in `x` to the retention times in `y` (using linear
#'   interpolation). This should only be used for chromatograms that were
#'   measured in the same measurement run (e.g. MS1 and corresponding MS2
#'   chromatograms from SWATH experiments).
#'
#' @param x [Chromatogram()] or [Chromatograms()] to be aligned against `y`.
#'
#' @param y [Chromatogram()] to which `x` should be aligned to.
#'
#' @param method `character(1)`
#'
#' @param ... additional parameters to be passed along to the alignment
#'     functions.
#'
#' @return `Chromatogram` (or `Chromatograms`) representing `x` aligned
#'     against `y`.
#'
#' @author Johannes Rainer, Michael Witting
#'
#' @md
#'
#' @examples
#'
#' chr1 <- Chromatogram(rtime = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
#'     intensity = c(3, 5, 14, 30, 24, 6, 2, 1, 1, 0))
#' chr2 <- Chromatogram(rtime = c(2.5, 3.42, 4.5, 5.43, 6.5),
#'     intensity = c(5, 12, 15, 11, 5))
#'
#' plot(chr1, col = "black")
#' points(rtime(chr2), intensity(chr2), col = "blue", type = "l")
#'
#' ## Align chr2 to chr1 without interpolation
#' res <- align(chr2, chr1)
#' rtime(res)
#' intensity(res)
#' points(rtime(res), intensity(res), col = "#00ff0080", type = "l")
#'
#' ## Align chr2 to chr1 with interpolation
#' res <- align(chr2, chr1, method = "approx")
#' points(rtime(res), intensity(res), col = "#ff000080", type = "l")
#'
#' legend("topright", col = c("black", "blue", "#00ff0080", "#ff000080"), lty = 1,
#'     legend = c("chr1", "chr2", "chr2 matchRtime", "chr2 approx"))
#'
setMethod("align", signature = c(x = "Chromatogram", y = "Chromatogram"),
          function(x, y, method = c("matchRtime", "approx"), ...) {
              .align_chromatogram(x = x, y = y, method = method, ...)
          })

#' @title Correlate chromatograms
#'
#' @aliases correlate
#'
#' @rdname correlate-Chromatogram
#'
#' @description
#'
#' Correlate intensities of two chromatograms with each other. If the two
#' `Chromatogram` objects have different retention times they are first
#' *aligned* to match data points in the first to data points in the second
#' chromatogram. See [align()] for more details.
#'
#' @param x [Chromatogram()] object.
#'
#' @param y [Chromatogram()] object.
#'
#' @param use `character(1)` passed to the `cor` function. See [cor()] for
#'     details.
#'
#' @param method `character(1)` passed to the `cor` function. See [cor()] for
#'     details.
#'
#' @param align `character(1)` defining the alignment method to be used. See
#'     [align()] for details. The value of this parameter is passed to the
#'     `method` parameter of `align`.
#'
#' @param ... optional parameters passed along to the `align` method.
#'
#' @return `numeric(1)` with the correlation coefficient.
#'
#' @author Michael Witting, Johannes Rainer
#'
#' @md
setMethod("correlate", signature = c(x = "Chromatogram", y = "Chromatogram"),
          function(x, y, use = "pairwise.complete.obs",
                   method = c("pearson", "kendall", "spearman"),
                   align = c("matchRtime", "approx"), ...) {
              .correlate_chromatogram(x, y, use = use, method = method,
                                      align = align, ...)
          })
