#' @rdname xchromatograms-class
XChromatograms <- function(data, phenoData, featureData, chromPeaks, ...) {
    chr <- Chromatograms(data, phenoData, featureData, ...)
    chr <- as(chr, "XChromatograms")
    if (!missing(chromPeaks)) {
        if (length(chromPeaks) != length(data))
            stop("Length of 'data' and 'chromPeaks' has to match")
        chr@chromPeaks <- matrix(chromPeaks, ncol = ncol(chr))
    }
    if (validObject(chr))
        chr
}
