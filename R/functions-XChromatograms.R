.validXChromatograms <- function(object) {
    txt <- character()
    if (length(object@.processHistory))
        if (!all(vapply(object@.processHistory,
                        function(z) inherits(z, "ProcessHistory"), logical(1))))
            txt <- c(txt, paste0("Only 'ProcessHistory' objects are allowed ",
                                 "in slot .processHistory"))
    if (!all(vapply(object, function(z)
        inherits(z, "XChromatogram"), logical(1))))
        txt <- c(txt, paste0("'object' should only contain 'XChromatogram' ",
                             "objects"))
    if (nrow(object@featureDefinitions)) {
        ## Check that we have a column "row" and that in that column we have
        ## indices from 1:nrow(object).
    }
    if (length(txt)) txt
    else TRUE
}

#' @rdname XChromatogram
#'
#' @param data For `XChromatograms`: `list` of `Chromatogram` or
#'     `XChromatogram` objects.
#'
#' @param phenoData For `XChromatograms`: either a `data.frame`,
#'     `AnnotatedDataFrame` or `NAnnotatedDataFrame` describing the
#'     phenotypical information of the samples.
#'
#' @param featureData For `XChromatograms`: either a `data.frame` or
#'     `AnnotatedDataFrame` with additional information for each row of
#'     chromatograms.
#'
#' @md
#'
#' @examples
#'
#' ## ---- Creation of XChromatograms ----
#' ##
#' ## Create a XChromatograms from Chromatogram objects
#' dta <- list(Chromatogram(rtime = 1:7, c(3, 4, 6, 12, 8, 3, 2)),
#'     Chromatogram(1:10, c(4, 6, 3, 4, 7, 13, 43, 34, 23, 9)))
#'
#' ## Create an XChromatograms without peak data
#' xchrs <- XChromatograms(dta)
#'
#' ## Create an XChromatograms with peaks data
#' pks <- list(matrix(c(4, 2, 5, 30, 12, NA), nrow = 1,
#'     dimnames = list(NULL, c("rt", "rtmin", "rtmax", "into", "maxo", "sn"))),
#'     NULL)
#' xchrs <- XChromatograms(dta, chromPeaks = pks)
#'
#' ## Create an XChromatograms from XChromatogram objects
#' dta <- lapply(dta, as, "XChromatogram")
#' chromPeaks(dta[[1]]) <- pks[[1]]
#'
#' xchrs <- XChromatograms(dta, nrow = 1)
#'
#' hasChromPeaks(xchrs)
#'
#' ## Load test files and extract chromatograms for a data slice
#' od <- readMSData(c(system.file("cdf/KO/ko15.CDF", package = "faahKO"),
#'     system.file("cdf/KO/ko16.CDF", package = "faahKO"),
#'     system.file("cdf/KO/ko18.CDF", package = "faahKO")),
#'     mode = "onDisk")
#'
#' ## Extract chromatograms for a m/z - retention time slice
#' chrs <- chromatogram(od, mz = 344, rt = c(2500, 3500))
#' chrs
#'
#' ## Perform peak detection using CentWave
#' xchrs <- findChromPeaks(chrs, param = CentWaveParam())
#' xchrs
#'
#' ## Do we have chromatographic peaks?
#' hasChromPeaks(xchrs)
#'
#' ## Process history
#' processHistory(xchrs)
#'
#' ## The chromatographic peaks, columns "row" and "column" provide information
#' ## in which sample the peak was identified.
#' chromPeaks(xchrs)
#'
#' ## Spectifically extract chromatographic peaks for one sample/chromatogram
#' chromPeaks(xchrs[1, 2])
#'
#' ## Plot the results
#' plot(xchrs)
#'
#' ## Plot the results using a different color for each sample
#' sample_colors <- c("#ff000040", "#00ff0040", "#0000ff40")
#' cols <- sample_colors[chromPeaks(xchrs)[, "column"]]
#' plot(xchrs, col = sample_colors, peakBg = cols)
#'
#' ## Indicate the peaks with a rectangle
#' plot(xchrs, col = sample_colors, peakCol = cols, peakType = "rectangle",
#'     peakBg = NA)
XChromatograms <- function(data, phenoData, featureData, chromPeaks, ...) {
    if (missing(data))
        return(new("XChromatograms"))
    if (!missing(chromPeaks)) {
        if (!is.list(chromPeaks) || length(chromPeaks) != length(data))
            stop("If provided, 'chromPeaks' has to be a list same length than",
                 " 'data'.")
        data <- mapply(data, chromPeaks, FUN = function(z, pks) {
            if (is(z, "Chromatogram"))
                z <- as(z, "XChromatogram")
            if (is.matrix(pks) && length(pks))
                chromPeaks(z) <- pks
            z
        })
    }
    object <- Chromatograms(data = data, phenoData = phenoData,
                            featureData = featureData, ...)
    object <- as(object, "XChromatograms")
    if (validObject(object)) object
}
