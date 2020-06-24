#' @include DataClasses.R do_findChromPeaks-functions.R

#' @title centWave-based peak detection in purely chromatographic data
#'
#' @description
#'
#' `findChromPeaks` on a [Chromatogram] or [MChromatograms] object with a
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
#' @param object a [Chromatogram] or [MChromatograms] object.
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
#' `findChromPeaks` on a [Chromatogram] or [MChromatograms] object with a
#' [MatchedFilterParam] parameter object performs matchedFilter-based peak
#' detection on purely chromatographic data. See [matchedFilter] for details
#' on the method and [MatchedFilterParam] for details on the parameter class.
#' Note that not all settings from the `MatchedFilterParam` will be used.
#' See [peaksWithMatchedFilter()] for the arguments used for peak detection
#' on purely chromatographic data.
#'
#' @param object a [Chromatogram] or [MChromatograms] object.
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
#' chromatograms. If `x` is a [MChromatograms()] object, each `Chromatogram` in
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
#' @param x [Chromatogram()] or [MChromatograms()] to be aligned against `y`.
#'
#' @param y [Chromatogram()] to which `x` should be aligned to.
#'
#' @param method `character(1)`
#'
#' @param ... additional parameters to be passed along to the alignment
#'     functions.
#'
#' @return `Chromatogram` (or `MChromatograms`) representing `x` aligned
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
#' If `correlate` is called on a single [MChromatograms()] object a pairwise
#' correlation of each chromatogram with each other is performed and a `matrix`
#' with the correlation coefficients is returned.
#'
#' Note that the correlation of two chromatograms depends also on their order,
#' e.g. `correlate(chr1, chr2)` might not be identical to
#' `correlate(chr2, chr1)`. The lower and upper triangular part of the
#' correlation matrix might thus be different.
#'
#' For correlating elements of a `MChromatograms` with each other it might be
#' sufficient to calculate just the upper triangular matrix. This can be done
#' by setting `full = FALSE`.
#'
#' @param x [Chromatogram()] or [MChromatograms()] object.
#'
#' @param y [Chromatogram()] or [MChromatograms()] object.
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
#' @param full `logical(1)` for `correlate` on a single `MChromatograms` object:
#'     whether the *full* correlation matrix should be calculated (default) or
#'     just the upper triangular matrix (and diagonal).
#'
#' @param ... optional parameters passed along to the `align` method.
#'
#' @return `numeric(1)` or `matrix` (if called on `MChromatograms` objects)
#'     with the correlation coefficient. If a `matrix` is returned, the rows
#'     represent the chromatograms in `x` and the columns the chromatograms in
#'     `y`.
#'
#' @author Michael Witting, Johannes Rainer
#'
#' @md
#'
#' @examples
#'
#' chr1 <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
#'     intensity = c(5, 29, 50, NA, 100, 12, 3, 4, 1, 3))
#' chr2 <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
#'     intensity = c(80, 50, 20, 10, 9, 4, 3, 4, 1, 3))
#' chr3 <- Chromatogram(rtime = 3:9 + rnorm(7, sd = 0.3),
#'     intensity = c(53, 80, 130, 15, 5, 3, 2))
#'
#' chrs <- MChromatograms(list(chr1, chr2, chr3))
#'
#' correlate(chr1, chr2)
#' correlate(chr2, chr1)
#'
#' correlate(chrs, chrs)
setMethod("correlate", signature = c(x = "Chromatogram", y = "Chromatogram"),
          function(x, y, use = "pairwise.complete.obs",
                   method = c("pearson", "kendall", "spearman"),
                   align = c("matchRtime", "approx"), ...) {
              .correlate_chromatogram(x, y, use = use, method = method,
                                      align = align, ...)
          })

#' @title Remove intensities from chromatographic data
#'
#' @aliases removeIntensity
#'
#' @rdname removeIntensity-Chromatogram
#'
#' @description
#'
#' `removeIntensities` allows to remove intensities from chromatographic data
#' matching certain conditions (depending on parameter `which`). The
#' intensities are actually not *removed* but replaced with `NA_real_`. To
#' actually **remove** the intensities (and the associated retention times)
#' use [clean()] afterwards.
#'
#' Parameter `which` allows to specify which intensities should be replaced by
#' `NA_real_`. By default (`which = "below_threshod"` intensities below
#' `threshold` are removed. If `x` is a `XChromatogram` or `XChromatograms`
#' object (and hence provides also chromatographic peak definitions within the
#' object) `which = "outside_chromPeak"` can be selected which removes any
#' intensity which is outside the boundaries of identified chromatographic
#' peak(s) in the chromatographic data.
#'
#' @param object an object representing chromatographic data. Can be a
#'     [Chromatogram()], [MChromatograms()], [XChromatogram()] or
#'     [XChromatograms()] object.
#'
#' @param which `character(1)` defining the condition to remove intensities.
#'     See description for details and options.
#'
#' @param threshold `numeric(1)` defining the threshold below which intensities
#'     are removed (if `which = "below_threshold"`).
#'
#' @return the input object with matching intensities being replaced by `NA`.
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @examples
#'
#' chr <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
#'     intensity = c(5, 29, 50, NA, 100, 12, 3, 4, 1, 3))
#'
#' ## Remove all intensities below 20
#' res <- removeIntensity(chr, threshold = 20)
#' intensity(res)
setMethod("removeIntensity", "Chromatogram",
          function(object, which = "below_threshold", threshold = 0) {
              which <- match.arg(which)
              if (which == "below_threshold")
                  object@intensity[which(object@intensity < threshold)] <- NA_real_
              object
          })

#' @title Normalize chromatographic data
#'
#' @aliases normalize,Chromatogram-method
#'
#' @name normalize-Chromatogram
#'
#' @description
#'
#' `normalize` *normalizes* the intensities of a chromatogram by dividing them
#' either by the maximum intensity (`method = "max"`) or total intensity
#' (`method = "sum"`) of the chromatogram.
#'
#' @param object [Chromatogram()] or [MChromatograms()] object.
#'
#' @param method `character(1)` defining whether each chromatogram should be
#'     normalized to its maximum signal or total signal.
#'
#' @return the normalized input object.
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @rdname normalize-Chromatogram
setMethod("normalize", "Chromatogram",
          function(object, method = c("max", "sum")) {
              method <- match.arg(method)
              .normalize_chromatogram(object, method)
})
