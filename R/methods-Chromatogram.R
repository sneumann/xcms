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
#' ## Loading a test data set with identified chromatographic peaks
#' data(faahko_sub)
#' ## Update the path to the files for the local system
#' dirname(faahko_sub) <- system.file("cdf/KO", package = "faahKO")
#' faahko_sub <- filterRt(faahko_sub, c(2500, 3700))
#'
#' ##
#' od <- as(filterFile(faahko_sub, 1L), "OnDiskMSnExp")
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
#' ## Loading a test data set with identified chromatographic peaks
#' data(faahko_sub)
#' ## Update the path to the files for the local system
#' dirname(faahko_sub) <- system.file("cdf/KO", package = "faahKO")
#' faahko_sub <- filterRt(faahko_sub, c(2500, 3700))
#'
#' ##
#' od <- as(filterFile(faahko_sub, 1L), "OnDiskMSnExp")
#'
#' ## Extract chromatographic data for a small m/z range
#' chr <- chromatogram(od, mz = c(272.1, 272.3))[1, 1]
#'
#' ## Identify peaks with default settings
#' xchr <- findChromPeaks(chr, MatchedFilterParam())
#'
#' ## Plot the identified peaks
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

#' @title Correlate chromatograms
#'
#' @aliases correlate
#'
#' @rdname correlate-Chromatogram
#'
#' @description
#'
#' **For `xcms` >= 3.15.3 please use [compareChromatograms()] instead of
#' `correlate`**
#'
#' Correlate intensities of two chromatograms with each other. If the two
#' `Chromatogram` objects have different retention times they are first
#' *aligned* to match data points in the first to data points in the second
#' chromatogram. See help on `alignRt` in [MSnbase::Chromatogram()] for more
#' details.
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
#'     help on `alignRt` in [MSnbase::Chromatogram()] for details. The value of
#'     this parameter is passed to the `method` parameter of `alignRt`.
#'
#' @param ... optional parameters passed along to the `alignRt` method such as
#'     `tolerance` that, if set to `0` requires the retention times to be
#'     identical.
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
#' ## Using `compareChromatograms` instead of `correlate`.
#' compareChromatograms(chr1, chr2)
#' compareChromatograms(chr2, chr1)
#'
#' compareChromatograms(chrs, chrs)
setMethod("correlate", signature = c(x = "Chromatogram", y = "Chromatogram"),
          function(x, y, use = "pairwise.complete.obs",
                   method = c("pearson", "kendall", "spearman"),
                   align = c("closest", "approx"), ...) {
              .Deprecated(new = "compareChromatograms")
              lst <- list(...)
              compareChromatograms(
                  x, y, ALIGNFUN = alignRt, FUN = cor,
                  ALIGNFUNARGS = c(list(method = align), lst),
                  FUNARGS = c(list(method = method, use = use), lst))
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
#' Note that [filterIntensity()] might be a better approach to subset/filter
#' chromatographic data.
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
                  object@intensity[which(object@intensity < threshold)] <-
                      NA_real_
              object
          })
