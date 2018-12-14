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
        if (!all(object@featureDefinitions$row %in% seq_len(nrow(object))))
            txt <- c(txt, paste0("Elements in column 'row' are outside of the",
                                 " number of rows of 'object'"))
        if (!all(unlist(object@featureDefinitions$peakidx) %in%
                 seq_len(nrow(chromPeaks(object)))))
            txt <- c(txt, paste0("peakidx in feature data does not match ",
                                 "the number of present chromatographic peaks"))
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
#'
#' ## Group chromatographic peaks across samples
#' prm <- PeakDensityParam(sampleGroup = rep(1, 3))
#' res <- groupChromPeaks(xchrs, param = prm)
#'
#' hasFeatures(res)
#' featureDefinitions(res)
#'
#' ## Delete the identified feature definitions
#' res <- dropFeatureDefinitions(res)
#' hasFeatures(res)
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

#' Subset the featureDefinitions `DataFrame` fts based on `pks` and `pks_sub`
#' being the `chromPeaks` before and after filtering.
#'
#' @author Johannes Rainer
#'
#' @noRd
.subset_features_on_chrom_peaks <- function(fts, pks, pks_sub) {
    if (nrow(fts)) {
        if (!is.null(rownames(pks)) && !is.null(rownames(pks_sub))) {
            ids_orig <- rownames(pks)
            ids_sub <- rownames(pks_sub)
        } else {
            cns <- intersect(colnames(pks), colnames(pks_sub))
            ids_orig <- apply(pks[, cns, drop = FALSE], 1, paste,
                              collapse = "-")
            if (length(ids_orig) != length(unique(ids_orig)))
                stop("Can not uniquely identify chromatographic peaks.")
            ids_sub <- apply(pks_sub[, cns, drop = FALSE], 1, paste,
                             collapse = "-")
        }
        fts$peakidx <- lapply(fts$peakidx, function(z) {
            newidx <- match(ids_orig[z], ids_sub)
            newidx[!is.na(newidx)]
        })
        fts <- fts[lengths(fts$peakidx) > 0, , drop = FALSE]
    }
    fts
}

#' Subset the chromPeaks matrix from an `XChromatograms` object. The
#' `chromPeaks` matrix is generated dynamically from the `chromPeaks` matrices
#' of each internal `XChromatogram` object, so there is not really a need to
#' subset the `chromPeaks` from an `XChromatograms` - only that we need this
#' to update the `"peakidx"` column of the `featureDefinitions`.
#'
#' Note: the chromPeaks matrix is extracted ordered by row.
#'
#' @author Johannes Rainer
#'
#' @noRd
.subset_chrom_peaks_xchromatograms <- function(x, i, j) {
    if (missing(i) & missing(j))
        return(x)
    if (missing(i)) i <- seq_len(nrow(x))
    if (missing(j)) j <- seq_len(ncol(x))
    x <- x[x[, "row"] %in% i & x[, "column"] %in% j, , drop = FALSE]
    if (nrow(x)) {
        x[, "row"] <- match(x[, "row"], i)
        x[, "column"] <- match(x[, "column"], j)
        x[order(x[, "row"], x[, "column"]), , drop = FALSE]
    } else x
}
