#' @description
#'
#' *Align* two chromatograms with similar, but not identical retention times.
#' The retention times of the second `Chromatogram` (`y`) will be replaced
#' by the retention times of the first (`x`) estimated by the `approx` function.
#'
#' @param x [Chromatogram()] object.
#'
#' @param y [Chromatogram()] object.
#'
#' @param na.value optional parameter allowing to specify the value with which
#'     `NA`s in the aligned chromatogram `y` will be replaced.
#'
#' @return `y` with the same retention times (and length) than `x`.
#'
#' @author Michael Witting
#'
#' @noRd
.align_chromatogram <- function(x, y, na.value = NA) {
    ## linear interpolation of y on x
    y_aligned <- approx(rtime(y), intensity(y), rtime(x))

    ## deal with NA
    if (!is.na(na.value)) {
        y_aligned$y[is.na(y_aligned$y)] <- na.value
    }

    ## correct rtime and int in chromatogram object
    y@rtime <- y_aligned$x
    y@intensity <- y_aligned$y

    y
}

#' @description
#'
#' Correlate intensities of two chromatograms with each other. If the two
#' `Chromatogram` objects have different retention times `.align_chromatogram`
#' is first called on them.
#'
#' @param x [Chromatogram()] object.
#'
#' @param y [Chromatogram()] object.
#'
#' @param na.value see `.align_chromagram`.
#'
#' @param use `character(1)` passed to the `cor` function. See [cor()] for
#'     details.
#'
#' @param method `character(1)` passed to the `cor` function. See [cor()] for
#'     details.
#'
#' @return `numeric(1)` with the correlation coefficient.
#'
#' @author Michael Witting
#'
#' @noRd
.correlate_chromatogram <- function(x, y,
                                    na.value = NA,
                                    use = "pairwise.complete.obs",
                                    method = c("pearson", "kendall",
                                               "spearman")) {

    ## check if x and y are on the same time scale
    if(!all(rtime(x) == rtime(y))) {
        y <- .align_chromatogram(x, y, na.value = na.value)
    }

    ## correlate values
    cor(intensity(x), intensity(y), use = use, method = method)
}

## wrapper function for multiple chromatograms
# are the MS1 chroms, y are the MS2
.correlate_chromatograms <-function(x, y) {

  # check if same amount of chroms are in x and y
  if(!ncol(x) == ncol(y)) {
    stop("Number of samples different between MS1 and MS2 chroms")
  }

  # make matrix to store correlation values
  corMatrix <- matrix(0, nrow = nrow(y), ncol = ncol(y))

  # iterate and make correlation
  # TODO
  # use biocParallel backend???

  # return
  return(corMatrix)

}
