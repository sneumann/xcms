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
    else lapply(object, validObject)
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
#' ## Loading a test data set with identified chromatographic peaks
#' data(faahko_sub)
#' ## Update the path to the files for the local system
#' dirname(faahko_sub) <- system.file("cdf/KO", package = "faahKO")
#'
#' ## Subset the dataset to the first and third file.
#' xod_sub <- filterFile(faahko_sub, file = c(1, 3))
#'
#' od <- as(xod_sub, "OnDiskMSnExp")
#'
#' ## Extract chromatograms for a m/z - retention time slice
#' chrs <- chromatogram(od, mz = 344, rt = c(2500, 3500))
#' chrs
#'
#' ## --------------------------------------------------- ##
#' ##       Chromatographic peak detection                ##
#' ## --------------------------------------------------- ##
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
#' ## --------------------------------------------------- ##
#' ##       Correspondence analysis                       ##
#' ## --------------------------------------------------- ##
#' ## Group chromatographic peaks across samples
#' prm <- PeakDensityParam(sampleGroup = rep(1, 2))
#' res <- groupChromPeaks(xchrs, param = prm)
#'
#' hasFeatures(res)
#' featureDefinitions(res)
#'
#' ## Plot the correspondence results. Use simulate = FALSE to show the
#' ## actual results. Grouped chromatographic peaks are indicated with
#' ## grey shaded rectangles.
#' plotChromPeakDensity(res, simulate = FALSE)
#'
#' ## Simulate a correspondence analysis based on different settings. Larger
#' ## bw will increase the smoothing of the density estimate hence grouping
#' ## chromatographic peaks that are more apart on the retention time axis.
#' prm <- PeakDensityParam(sampleGroup = rep(1, 3), bw = 60)
#' plotChromPeakDensity(res, param = prm)
#'
#' ## Delete the identified feature definitions
#' res <- dropFeatureDefinitions(res)
#' hasFeatures(res)
XChromatograms <- function(data, phenoData, featureData, chromPeaks,
                           chromPeakData, ...) {
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
    if (!missing(chromPeakData)) {
        if (missing(chromPeaks))
            stop("If 'chromPeakData' is provided, also 'chromPeaks' is required")
        if (!is.list(chromPeakData) || length(chromPeakData) != length(data))
            stop("If provided, 'chromPeakData' has to be a list same length ",
                 "than 'data'.")
        data <- mapply(data, chromPeakData, FUN = function(z, pkd) {
            if (length(pkd))
                chromPeakData(z) <- pkd
        })
    }
    object <- MChromatograms(data = data, phenoData = phenoData,
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
            cns <- cns[!(cns %in% c("row", "column"))]
            ids_orig <- apply(pks[, cns, drop = FALSE], 1, paste,
                              collapse = "-")
            ## if (length(ids_orig) != length(unique(ids_orig)))
            ##     stop("Can not uniquely identify chromatographic peaks.")
            ids_sub <- apply(pks_sub[, cns, drop = FALSE], 1, paste,
                             collapse = "-")
        }
        for (i in seq_len(nrow(fts))) {
            fts$peakidx[[i]] <- unname(
                which(ids_sub %in% ids_orig[fts$peakidx[[i]]] &
                      pks_sub[, "row"] == fts$row[i]))
        }
        fts <- extractROWS(fts, which(lengths(fts$peakidx) > 0))
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

#' Convenience function to plot a peak density given provided chromPeaks and
#' a PeakDensityParam object.
#'
#' @param pks `matrix` with chromatographic peaks.
#'
#' @param param `PeakDensityParam`
#'
#' @param xlim optional definition of the x-axis limits.
#'
#' @param main optional title.
#'
#' @param xlab, ylab x- and y-axis labels.
#'
#' @param peakCol foreground color definition for peaks. Either 1 or length
#'     equal to `ncol(pks`).
#'
#' @param peakBg background color definition for peaks.
#'
#' @param peakPch point character.
#'
#' @author Johannes Rainer
#'
#' @noRd
.plot_chrom_peak_density <- function(pks, fts, param, xlim = range(pks[, "rt"]),
                                     main = NA, xlab = "retention time",
                                     ylab = "sample", peakCol = "#00000060",
                                     peakBg = "#00000020", peakPch = 1,
                                     simulate = TRUE, col = "black",
                                     ylim = range(pks[, "column"]), ...) {
    pks_count <- nrow(pks)
    if (pks_count) {
        smpl_col <- which(colnames(pks) == "sample")
        if (!length(smpl_col))
            smpl_col <- which(colnames(pks) == "column")
        if (length(peakCol) == 1)
            peakCol <- rep(peakCol, pks_count)
        if (length(peakBg) == 1)
            peakBg <- rep(peakBg, pks_count)
        if (length(peakPch) == 1)
            peakPch <- rep(peakPch, pks_count)
        if (length(peakCol) != pks_count) {
            warning("Length of 'peakCol' does not match the number of peaks. ",
                    "Using peakCol[1] for all.")
            peakCol <- rep(peakCol[1], pks_count)
        }
        if (length(peakBg) != pks_count) {
            warning("Length of 'peakBg' does not match the number of peaks. ",
                    "Using peakBg[1] for all.")
            peakBg <- rep(peakBg[1], pks_count)
        }
        if (length(peakPch) != pks_count) {
            warning("Length of 'peakPch' does not match the number of peaks. ",
                    "Using peakPch[1] for all.")
            peakPch <- rep(peakPch[1], pks_count)
        }
        bw <- bw(param)
        full_rt_range <- range(pks[, "rt"])
        dens_from <- full_rt_range[1] - 3 * bw
        dens_to <- full_rt_range[2] + 3 * bw
        densN <- max(512, 2 * 2^(ceiling(log2(diff(full_rt_range) / (bw / 2)))))
        sample_groups <- sampleGroups(param)
        dens <- density(pks[, "rt"], bw = bw, from = dens_from, to = dens_to,
                        n = densN)
        yl <- c(0, max(dens$y))
        min_max_smple <- ylim
        ypos <- seq(from = yl[1], to = yl[2],
                    length.out = diff(min_max_smple) + 1)
        plot(pks[, "rt"], ypos[pks[, smpl_col]], xlim = xlim, ylim = yl,
             main = main, yaxt = "n", ylab = ylab, xlab = xlab,
             col = peakCol, bg = peakBg, pch = peakPch, ...)
        points(x = dens$x, y = dens$y, type = "l", ...)
        axis(side = 2, at = ypos, labels = seq(from = min_max_smple[1],
                                               to = min_max_smple[2]))
        sample_groups <- sampleGroups(param)
        sample_groups_table <- table(sample_groups)
        dens_max <- max(dens$y)
        if (simulate) {
            snum <- 0
            while(dens$y[max_y <- which.max(dens$y)] > dens_max / 20 &&
                  snum < maxFeatures(param)) {
                      feat_range <- descendMin(dens$y, max_y)
                      dens$y[feat_range[1]:feat_range[2]] <- 0
                      feat_idx <- which(pks[, "rt"] >= dens$x[feat_range[1]] &
                                        pks[, "rt"] <= dens$x[feat_range[2]])
                      tt <- table(sample_groups[pks[feat_idx, smpl_col]])
                      if (!any(tt / sample_groups_table[names(tt)] >=
                               minFraction(param) & tt >= minSamples(param)))
                          next
                      rect(xleft = min(pks[feat_idx, "rt"]), ybottom = yl[1],
                           xright = max(pks[feat_idx, "rt"]), ytop = yl[2],
                           border = "#00000040", col = "#00000020")
                  }
        } else {
            if (!missing(fts) && nrow(fts)) {
                rect(xleft = fts$rtmin, xright = fts$rtmax,
                     ybottom = rep(yl[1], nrow(fts)),
                     ytop = rep(yl[2], nrow(fts)),
                     border = "#00000040", col = "#00000020")
                abline(v = fts$rtmed, col = "#00000040", lty = 2)
            } else warning("No feature definitions present. Either use ",
                           "'groupChromPeaks' first or set 'simulate = TRUE'")
        }
    } else {
        plot(3, 3, pch = NA, xlim = xlim, main = main, xlab = xlab, ylab = ylab)
    }
}
