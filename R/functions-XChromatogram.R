.validXChromatogram <- function(object) {
    txt <- character()
    if (nrow(object@chromPeaks)) {
        if (!all(.CHROMPEAKS_REQ_NAMES %in% colnames(object@chromPeaks)))
            txt <- c(txt, paste0("chromPeaks matrix does not have all required",
                                 " columns: ", paste(.CHROMPEAKS_REQ_NAMES,
                                                     collapse = ",")))
        if (!is.numeric(object@chromPeaks))
            txt <- c(txt, "chromPeaks should be a numeric matrix")
        if (!is.null(colnames(object@chromPeaks)) &&
            any(object@chromPeaks[, "rtmax"] < object@chromPeaks[, "rtmin"]))
            txt <- c(txt, "rtmax has to be larger than rtmin")
    }
    if (length(txt)) txt
    else TRUE
}

#' @title Containers for chromatographic and peak detection data
#'
#' @aliases XChromatogram-class
#'
#' @description
#'
#' The `XChromatogram` object allows to store chromatographic data (e.g.
#' an extracted ion chromatogram) along with identified chromatographic peaks
#' within that data. The object inherits all functions from the [Chromatogram()]
#' object in the `MSnbase` package.
#'
#' All functions are described (grouped into topic-related sections) after the
#' **Arguments** section.
#'
#' @section Creation of objects:
#'
#' Objects can be created with the contructor function `XChromatogram`.
#'
#' @param x An `XChromatogram` object.
#'
#' @param pks `matrix` with required columns `"rt"`, `"rtmin"`, `"rtmax"`,
#'     `"into"`, `"maxo"` and `"sn"`.
#'
#' @param object An `XChromatogram` object.
#'
#' @return
#'
#' For `chromPeaks`: a `matrix` with columns:
#' - `"rt"`: the retention time of the peak apex.
#' - `"rtmin"`: the lower peak boundary.
#' - `"rtmax"`: the upper peak boundary.
#' - `"into"`: the integrated peak signal (area of the peak).
#' - `"maxo"`: the maximum intensity of the peak.
#' Note that the matrix could contain also more columnds depending on the
#' peak detection method.
#'
#' All other methods return an `XChromatogram` object.
#'
#' @md
#'
#' @author Johannes Rainer
#'
#' @rdname XChromatogram
#'
#' @examples
#'
#' ## Create a XChromatogram object
#' chr <- Chromatogram(rtime = 1:10,
#'     intensity = c(4, 8, 14, 19, 18, 12, 9, 8, 5, 2))
#' pks <- matrix(nrow = 1, ncol = 6)
#' colnames(pks) <- c("rt", "rtmin", "rtmax", "into", "maxo", "sn")
#' pks[, "rtmin"] <- 2
#' pks[, "rtmax"] <- 9
#' pks[, "rt"] <- 4
#' pks[, "maxo"] <- 19
#' pks[, "into"] <- 93
#'
#' xchr <- XChromatogram(chr, pks)
#' xchr
XChromatogram <- function(x, pks) {
    if (missing(x))
        x <- Chromatogram()
    if (!is(x, "Chromatogram"))
        stop("'x' has to be a 'Chromatogram' object")
    if (missing(pks))
        pks <- matrix(ncol = length(.CHROMPEAKS_REQ_NAMES), nrow = 0,
                      dimnames = list(character(), .CHROMPEAKS_REQ_NAMES))
    else if (!is.matrix(pks))
        stop("'x' has to be a 'matrix'")
    x <- as(x, "XChromatogram")
    x@chromPeaks <- pks
    if (validObject(x)) x
}
