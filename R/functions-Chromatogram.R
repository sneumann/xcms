.SUPPORTED_AGG_FUN_CHROM <- c("sum", "max", "min", "mean")
names(.SUPPORTED_AGG_FUN_CHROM) <-
    c("Total ion chromatogram (TIC).", "Base peak chromatogram (BPC).",
      "Intensity representing the minimum intensity across the mz range.",
      "Intensity representing the mean intensity across the mz range.")

#' @title Plot Chromatogram objects
#' 
#' @description \code{plotChromatogram} creates a chromatogram plot for a
#'     single \code{\link[MSnbase]{Chromatogram}} object or a \code{list} of
#'     \code{\link[MSnbase]{Chromatogram}} objects (one line for each
#'     \code{\link[MSnbase]{Chromatogram}}/sample).
#'
#' @details The \code{plotChromatogram} function allows to efficiently plot
#'     the chromatograms of several samples into a single plot.
#' 
#' @param x For \code{plotChromatogram}: \code{list} of
#'     \code{\link[MSnbase]{Chromatogram}} objects. Such as extracted from an
#'     \code{\link{XCMSnExp}} object by the \code{\link{chromatogram}}
#'     method.
#'     For \code{highlightChromPeaks}: \code{XCMSnExp} object with the detected
#'     peaks.
#'
#' @param rt For \code{plotChromatogram}: \code{numeric(2)}, optional parameter
#'     to subset each \code{Chromatogram} by retention time prior to plotting.
#'     Alternatively, the plot could be subsetted by passing a \code{xlim}
#'     parameter.
#'     For \code{highlightChromPeaks}: \code{numeric(2)} with the
#'     retention time range from which peaks should be extracted and plotted.
#' 
#' @param col For \code{plotChromatogram}: color definition for each
#'     line/sample. Has to have the same length as samples/elements in \code{x},
#'     otherwise \code{col[1]} is recycled to generate a vector of
#'     \code{length(x)}.
#'     For \code{highlightChromPeaks}: color to be used to fill the
#'     rectangle.
#'
#' @param lty the line type. See \code{\link[graphics]{plot}} for more details.
#'
#' @param type the plotting type. See \code{\link[graphics]{plot}} for more
#'     details.
#'     For \code{highlightChromPeaks}: \code{character(1)} defining how the peak
#'     should be highlighted: \code{type = "rect"} draws a rectangle
#'     representing the peak definition, \code{type = "point"} indicates a
#'     chromatographic peak with a single point at the position of the peak's
#'     \code{"rt"} and \code{"maxo"}.
#' 
#' @param xlab \code{character(1)} with the label for the x-axis.
#'
#' @param ylab \code{character(1)} with the label for the y-axis.
#'
#' @param main The title for the plot. For \code{plotChromatogram}: if
#'     \code{main = NULL} the mz range of the \code{Chromatogram} object(s) will
#'     be used as the title.
#'
#' @param ... additional parameters to the \code{\link{matplot}} or \code{plot}
#'     function.
#'
#' @seealso \code{\link{chromatogram}} for how to extract a list of
#'     \code{\link[MSnbase]{Chromatogram}} objects from an
#'     \code{\link{XCMSnExp}} objects.
#' 
#' @author Johannes Rainer
#'
#' @examples
#'
#' ## Perform a fast peak detection.
#' library(xcms)
#' library(faahKO)
#' faahko_3_files <- c(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
#'                     system.file('cdf/KO/ko16.CDF', package = "faahKO"),
#'                     system.file('cdf/KO/ko18.CDF', package = "faahKO"))
#'
#' od <- readMSData2(faahko_3_files)
#'
#' od <- findChromPeaks(od, param = CentWaveParam(snthresh = 20, noise = 10000))
#'
#' rtr <- c(2600, 2750)
#' mzr <- c(344, 344)
#' chrs <- chromatogram(od, rt = rtr, mz = mzr)
#'
#' ## Plot a single chromatogram
#' plotChromatogram(chrs[[1]])
#'
#' ## Plot all chromatograms at once, using different colors for each.
#' plotChromatogram(chrs, col = c("#FF000080", "#00FF0080", "#0000FF80"), lwd = 2)
#'
#' ## Highlight identified chromatographic peaks.
#' highlightChromPeaks(od, rt = rtr, mz = mzr,
#'     col = c("#FF000005", "#00FF0005", "#0000FF05"),
#'     border = c("#FF000040", "#00FF0040", "#0000FF40"))
#'
plotChromatogram <- function(x, rt, col = "#00000060",
                             lty = 1, type = "l", xlab = "retention time",
                             ylab = "intensity", main = NULL, ...) {
    if (!is.list(x) & !is(x, "Chromatogram"))
        stop("'x' should be a Chromatogram object or a list of Chromatogram",
             " objects.")
    if (is(x, "Chromatogram"))
        x <- list(x)
    isOK <- lapply(x, function(z) {
        if (is(z, "Chromatogram")) {
            return(TRUE)
        } else {
            if (is.na(z))
                return(TRUE)
        }
        FALSE
    })
    if (any(!unlist(isOK)))
        stop("if 'x' is a list it should only contain Chromatogram objects")
    ## Subset the Chromatogram objects if rt provided.
    if (!missing(rt)) {
        rt <- range(rt)
        x <- lapply(x, function(z) {
            if (is(z, "Chromatogram"))
                filterRt(z, rt = rt)
        })
    }
    if (length(col) != length(x)) {
        col <- rep(col[1], length(x))
    }
    ## If main is NULL use the mz range.
    if (is.null(main)) {
        mzr <- range(lapply(x, mz), na.rm = TRUE, finite = TRUE)
        main <- paste0(format(mzr, digits = 7), collapse = " - ")
    }
    ## Number of measurements we've got per chromatogram. This can be different
    ## between samples, from none (if not a single measurement in the rt/mz)
    ## to the number of data points that were actually measured.
    lens <- unique(lengths(x))
    max_len <- max(lens)
    max_len_vec <- rep_len(NA, max_len)
    ## Generate the matrix of rt values, columns are samples, rows retention
    ## time values. Fill each column with NAs up to the maximum number of values
    ## we've got in a sample/file.
    rts <- do.call(cbind, lapply(x, function(z) {
        cur_len <- length(z)
        if (cur_len == 0)
            max_len_vec
        else {
            ## max_len_vec[,] <- NA  ## don't need that. get's copied.
            max_len_vec[seq_len(cur_len)] <- rtime(z)
            max_len_vec
        }
    }))
    ## Same for the intensities.
    ints <- do.call(cbind, lapply(x, function(z) {
        cur_len <- length(z)
        if (length(z) == 0)
            max_len_vec
        else {
            ## max_len_vec[,] <- NA  ## don't need that. get's copied.
            max_len_vec[seq_len(cur_len)] <- intensity(z)
            max_len_vec
        }
    }))
    ## Define the x and y limits
    x_lim <- c(0, 1)
    y_lim <- c(0, 1)
    if (all(is.na(rts)))
        if (!missing(rt))
            x_lim <- range(rt)
    else
        x_lim <- range(rts, na.rm = TRUE, finite = TRUE)
    if (!all(is.na(ints)))
        y_lim <- range(ints, na.rm = TRUE, finite = TRUE)
    ## Identify columns that have only NAs in either intensity or rt - these
    ## will not be plotted.
    keepCol <- which(apply(ints, MARGIN = 2, function(z) any(!is.na(z))) |
                     apply(rts, MARGIN = 2, function(z) any(!is.na(z))))
    ## Finally plot the data.
    if (length(keepCol)) {
        matplot(x = rts[, keepCol, drop = FALSE],
                y = ints[, keepCol, drop = FALSE], type = type, lty = lty,
                col = col[keepCol], xlab = xlab, ylab = ylab, main = main,
                ...)
    } else
        plot(x = 3, y = 3, pch = NA, xlab = xlab, ylab = ylab, main = main,
             xlim = x_lim, ylim = y_lim)
}

