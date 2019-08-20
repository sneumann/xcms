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
#' @noRd
.correlate_chromatogram <- function(x, y, use = "pairwise.complete.obs",
                                    method = c("pearson", "kendall",
                                               "spearman"),
                                    align = c("matchRtime", "approx"),
                                    ...) {
    align <- match.arg(align)
    if(length(x) != length(y) || !all(rtime(x) %in% rtime(y)))
        x <- .align_chromatogram(x, y, method = align, ...)
    cor(intensity(x), intensity(y), use = use, method = method)
}
