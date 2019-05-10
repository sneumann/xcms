#' @description
#'
#' *Align* chromatogram `x` to chromatogram `y`.
#' The retention times of the first `Chromatogram` (`x`) will be replaced
#' by the retention times of the second (`y`) estimated by the `approx`
#' function.
#'
#' This type of alignment should be used if the chromatograms were measured in
#' the same run (e.g. in SWATH experiments).
#'
#' @param x [Chromatogram()] object that will be aligned against `y`
#'
#' @param y [Chromatogram()] object.
#'
#' @param na.value optional parameter allowing to specify the value with which
#'     `NA`s in the aligned chromatogram `x` will be replaced.
#'
#' @return aligned `x` `Chromatogram` with the same number of data points and
#'     same retention times than `y`.
#'
#' @author Michael Witting
#'
#' @noRd
.align_chromatogram_approx <- function(x, y, na.value = NA, ...) {
    x_aligned <- approx(rtime(x), intensity(x), rtime(y))
    if (!is.na(na.value)) {
        x_aligned$y[is.na(x_aligned$y)] <- na.value
    }
    ## correct rtime and int in chromatogram object
    x@rtime <- x_aligned$x
    x@intensity <- x_aligned$y
    x
}

#' @description
#'
#' *Align* chromatogram `x` against chromatogram `y` matching each data point
#' in `x` to the data point in `y` with the smallest difference between their
#' retention times, but only if the difference in retention times is `<=` the
#' minimal average difference between retention times in `x` or `y` (otherwise
#' the data point will be discarded).
#'
#' @param x [Chromatogram()] object that will be aligned against `y`
#'
#' @param y [Chromatogram()] object.
#'
#' @param na.value optional parameter allowing to specify the value with which
#'     `NA`s in the aligned chromatogram `x` will be replaced.
#'
#' @param ... optional parameters to be passed along to the `.match_closest`
#'     function.
#'
#' @return aligned `x` `Chromatogram` with the same number of data points and
#'     same retention times than `y`.
#'
#' @author Johannes Rainer
#'
#' @noRd
.align_chromatogram_match_rtime <- function(x, y, na.value = NA, ...) {
    idx <- .match_closest(rtime(x), rtime(y), ...)
    not_na <- !is.na(idx)
    x@rtime <- rtime(y)
    new_int <- rep(na.value, length(y))
    if (any(not_na))
        new_int[idx[not_na]] <- x@intensity[not_na]
    x@intensity <- new_int
    x
}

#' @title Align two chromatograms
#'
#' @description
#'
#' Align chromatogram `x` against chromatogram `y`. The resulting chromatogram
#' has the same length (number of data points) than `y` and the same retention
#' times thus allowing to perform any pair-wise comparisons between the
#' chromatograms. Parameter `method` allows to specify which alignment method
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
#' @param x `Chromatogram` that should be aligned against `y`.
#'
#' @param y `Chromatogram` to which `x` should be aligned to.
#'
#' @param method `character(1)`
#'
#' @param na.value optional value with which `NA` intensities in the resulting
#'     `Chromatogram` should be replaced with.
#'
#' @param ... additional parameters to be passed along to the alignment
#'     functions.
#'
#' @return `Chromatogram` representing `x` aligned against `y`.
#'
#' @author Johannes Rainer, Michael Witting
#'
#' @noRd
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
#' res <- .align_chromatogram(chr2, chr1)
#' rtime(res)
#' intensity(res)
#' points(rtime(res), intensity(res), col = "#00ff0080", type = "l")
#'
#' ## Align chr2 to chr1 with interpolation
#' res <- .align_chromatogram(chr2, chr1, method = "approx")
#' points(rtime(res), intensity(res), col = "#ff000080", type = "l")
#'
#' legend("topright", col = c("black", "blue", "#00ff0080", "#ff000080"), lty = 1,
#'     legend = c("chr1", "chr2", "chr2 matchRtime", "chr2 approx"))
#'
.align_chromatogram <- function(x, y, method = c("matchRtime", "approx"),
                                na.value = NA, ...) {
    method <- match.arg(method)
    switch(method,
           matchRtime = .align_chromatogram_match_rtime(
               x, y, na.value = na.value, ...),
           approx = .align_chromatogram_approx(x, y, na.value = na.value, ...))
}



#' @title Correlate chromatograms
#'
#' @description
#'
#' Correlate intensities of two chromatograms with each other. If the two
#' `Chromatogram` objects have different retention times they are first
#' *aligned* to match data points in the first to data points in the second
#' chromatogram.
#'
#' @details
#'
#' Chromatograms will almost never have the same retention times which makes
#' correlation of their intensities a non-trivial task. The default is to
#'
#' Alternatively, with `interpolate = TRUE` the
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
#' @param interpolate `logical(1)` defining whether interpolation should take
#'     place during *alignment* of the chromatograms. See details for more
#'     information.
#'
#' @param ... optional parameters passed along to the `.align_chromatogram`
#'     function.
#'
#' @return `numeric(1)` with the correlation coefficient.
#'
#' @author Michael Witting, Johannes Rainer
#'
#' @noRd
.correlate_chromatogram <- function(x, y, use = "pairwise.complete.obs",
                                    method = c("pearson", "kendall",
                                               "spearman"),
                                    align = c("matchRtime", "approx"),
                                    ...) {
    align <- match.arg(align)
    if(length(x) != length(y) || !all(rtime(x) %in% rtime(y)))
        x <- .align_chromatogram(x, y, interpolate = interpolate, ...)
    cor(intensity(x), intensity(y), use = use, method = method)
}
