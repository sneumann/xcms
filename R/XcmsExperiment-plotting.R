## Plotting functions for XcmsExperiment.

## plot: don't think that would be possible... has lower priority, maybe just
## provide a XIC plot.
## plotXIC

#' @title Visualization of Alignment Results
#'
#' @description
#'
#' The `plotAdjustedRtime` function plots the difference between the adjusted
#' and *raw* retention times on the y-axis against the raw retention times on
#' the x-axis. Each line represents the results for one sample (file).
#' If alignment was performed using the *peak groups* method (see
#' [adjustRtime()] for more infromation) also the peak groups used in the
#' alignment are visualized.
#'
#' @param object A [XcmsExperiment()] or [XCMSnExp()] object with the alignment
#'     results.
#'
#' @param col color(s) for the individual lines. Has to be of length 1 or equal
#'     to the number of samples.
#'
#' @param lwd line width for the lines of the individual samples.
#'
#' @param lty line type for the lines of the individual samples.
#'
#' @param type plot *type* (see [par()] for options; defaults to `type = "l"`).
#'
#' @param adjustedRtime `logical(1)` whether adjusted or raw retention times
#'     should be shown on the x-axis.
#'
#' @param xlab the label for the x-axis.
#'
#' @param ylab the label for the y-axis.
#'
#' @param peakGroupsCol color to be used for the peak groups (only if
#'     alignment was performed using the *peak groups* method.
#'
#' @param peakGroupsPch point character (`pch`) to be used for the peak
#'     groups (only if alignment was performed using the *peak groups*
#'     method.
#'
#' @param peakGroupsLty line type (`lty`) to be used to connect points for
#'     each peak groups (only if alignment was performed using the
#'     *peak groups* method.
#'
#' @param ylim optional `numeric(2)` with the upper and lower limits on
#'     the y-axis.b
#'
#' @param ... Additional arguments to be passed down to the `plot`
#'     function.
#'
#' @author Johannes Rainer
#'
#' @examples
#'
#' ## Load a test data set with detected peaks
#' faahko_sub <- loadXcmsData("faahko_sub2")
#'
#' ## Disable parallel processing for this example
#' register(SerialParam())
#'
#' ## Performing the peak grouping using the "peak density" method.
#' p <- PeakDensityParam(sampleGroups = c(1, 1, 1))
#' res <- groupChromPeaks(faahko_sub, param = p)
#'
#' ## Perform the retention time adjustment using peak groups found in both
#' ## files.
#' fgp <- PeakGroupsParam(minFraction = 1)
#' res <- adjustRtime(res, param = fgp)
#'
#' ## Visualize the impact of the alignment.
#' plotAdjustedRtime(res, adjusted = FALSE)
#' grid()
plotAdjustedRtime <- function(object, col = "#00000080", lty = 1, lwd = 1,
                              type = "l", adjustedRtime = TRUE,
                              xlab = ifelse(adjustedRtime,
                                            yes = expression(rt[adj]),
                                            no = expression(rt[raw])),
                              ylab = expression(rt[adj]-rt[raw]),
                              peakGroupsCol = "#00000060",
                              peakGroupsPch = 16, peakGroupsLty = 3,
                              ylim, ...) {
    if (!(inherits(object, "XCMSnExp") |
          inherits(object, "XcmsExperiment")))
        stop("'object' needs to be either a 'XcmsExperiment' or ",
             "'XCMSnExp' object.")
    if (!hasAdjustedRtime(object))
        warning("No alignment/retention time correction results present.")
    rt <- rtime(object, adjusted = FALSE)
    rtadj <- rtime(object, adjusted = TRUE)
    .plot_adjusted_rtime(rt, rtadj, fromFile(object), col = col, lty = lty,
                         lwd = lwd, type = type, xlab = xlab, ylab = ylab,
                         adjustedRtime = adjustedRtime, ...)
    ph <- processHistory(object, type = .PROCSTEP.RTIME.CORRECTION)
    if (length(ph)) {
        ph <- ph[[length(ph)]]
        if (inherits(ph, "XProcessHistory")) {
            prm <- processParam(ph)
            if (is(prm, "PeakGroupsParam")) {
                rt <- split(rt, fromFile(object))
                rtadj <- split(rtadj, fromFile(object))
                peak_group <- peakGroupsMatrix(prm)
                idx_sub <- subset(prm)
                if (length(idx_sub)) {
                    rt <- rt[idx_sub]
                    rtadj <- rtadj[idx_sub]
                }
                .plot_peak_groups(rt, rtadj, peak_group,
                                  col = peakGroupsCol, pch = peakGroupsPch,
                                  lty = peakGroupsLty,
                                  adjustedRtime = adjustedRtime, ...)
            }
        }
    }
}

.plot_peak_groups <- function(rt, rtadj, peak_group, col = "#00000060",
                              pch = 16, lty = 3, adjustedRtime = TRUE, ...) {
    dev.hold()
    on.exit(dev.flush())
    peak_group_adj <- peak_group
    for (i in seq_len(ncol(peak_group)))
        peak_group_adj[, i] <- .applyRtAdjustment(peak_group[, i],
                                                  rt[[i]], rtadj[[i]])
    diff_rt <- peak_group_adj - peak_group
    if (adjustedRtime)
        xrt <- peak_group_adj
    else xrt <- peak_group
    for (i in 1:nrow(xrt)) {
        idx <- order(diff_rt[i, ])
        plot.xy(xy.coords(x = xrt[i, ][idx], y = diff_rt[i, ][idx]),
                col = col, type = "b", pch = pch, lty = lty, ...)
    }
}

.plot_adjusted_rtime <- function(rt, rtadj, from_file,
                                 col = "#00000080", lty = 1, lwd = 1,
                                 type = "l", adjustedRtime = TRUE,
                                 xlab = ifelse(adjustedRtime,
                                               yes = expression(rt[adj]),
                                               no = expression(rt[raw])),
                                 ylab = expression(rt[adj]-rt[raw]),
                                 ylim, main = "", ...) {
    diffRt <- rtadj - rt
    diffRt <- split(diffRt, from_file)
    l <- length(diffRt)
    if (adjustedRtime)
        rt <- rtadj
    xlim <- range(rt, na.rm = TRUE)
    rt <- split(rt, from_file)
    if (length(col) != l) {
        if (length(col) > 1)
            warning("length of 'col' does not match number of samples. Will ",
                    "use 'col[1]' for all.")
        col <- rep(col[1L], l)
    }
    if (length(lty) != l) {
        if (length(lty) > 1)
            warning("length of 'lty' does not match number of samples. Will ",
                    "use 'lty[1]' for all.")
        lty <- rep(lty[1L], l)
    }
    if (length(lwd) != l) {
        if (length(lwd) > 1)
            warning("length of 'lwd' does not match number of samples. Will ",
                    "use 'lwd[1]' for all.")
        lwd <- rep(lwd[1L], l)
    }
    if (missing(ylim))
        ylim <- range(diffRt, na.rm = TRUE)
    dev.hold()
    on.exit(dev.flush())
    plot.new()
    plot.window(xlim = xlim, ylim = ylim, ...)
    axis(side = 1, ...)
    axis(side = 2, ...)
    box()
    title(main = main, xlab = xlab, ylab = ylab, ...)
    for (i in seq_along(diffRt))
        plot.xy(xy.coords(x = rt[[i]], y = diffRt[[i]]), col = col[i],
                lty = lty[i], type = type, lwd = lwd[i])
}

#' @title General visualizations of peak detection results
#'
#' @description
#'
#' `plotChromPeaks` plots the identified chromatographic
#' peaks from one file into the plane spanned by the retention time (x-axis)
#' and m/z (y-axis) dimension. Each chromatographic peak is plotted as a
#' rectangle representing its width in RT and m/z dimension.
#'
#' `plotChromPeakImage` plots the number of detected peaks for
#' each sample along the retention time axis as an *image* plot, i.e.
#' with the number of peaks detected in each bin along the retention time
#' represented with the color of the respective cell.
#'
#' @details
#'
#' The width and line type of the rectangles indicating the detected
#' chromatographic peaks for the `plotChromPeaks` function can be
#' specified using the `par` function, i.e. with `par(lwd = 3)`
#' and `par(lty = 2)`, respectively.
#'
#' @param add For `plotChromPeaks`: `logical(1)` whether the plot
#'     should be added to an existing plot or if a new plot should be created.
#'
#' @param binSize For `plotChromPeakImage`: `numeric(1)` defining the
#'     size of the bins along the x-axis (retention time). Defaults to
#'     `binSize = 30`, peaks within each 30 seconds will thus counted and
#'     plotted.
#'
#' @param border For `plotChromPeaks`: the color for the rectangles'
#'     border.
#'
#' @param col For `plotChromPeaks`: the color to be used to fill the
#'     rectangles.
#'
#' @param file For `plotChromPeaks`: `integer(1)` specifying the
#'     index of the file within `x` for which the plot should be created.
#'     Defaults to `file = 1`.
#'
#' @param log For `plotChromPeakImage`: `logical(1)` whether the peak
#'     counts should be log2 transformed before plotting.
#'
#' @param main `character(1)` defining the plot title. By default (i.e.
#'     `main = NULL`) the name of the file will be used as title.
#'
#' @param msLevel `integer(1)` defining the MS level from which the peaks
#'     should be visualized.
#'
#' @param x A [XcmsExperiment()] or [XCMSnExp()] object.
#'
#' @param xlim `numeric(2)` specifying the x-axis limits (retention time
#'     dimension). Defaults to `xlim = NULL` in which case the full retention
#'     time range of the file is used.
#'
#' @param yaxt For `plotChromPeakImage`: `character(1)` defining
#'     whether y-axis labels should be added. To disable the y-axis use
#'     `yaxt = "n"`. For any other value of `yaxt` the axis will be
#'     drawn. See [par()] help page for more details.
#'
#' @param ylim For `plotChromPeaks`: `numeric(2)` specifying the
#'     y-axis limits (m/z dimension). Defaults to `ylim = NULL` in which
#'     case the full m/z range of the file is used.
#'
#' @param xlab `character(1)` defining the x-axis label.
#'
#' @param ylab For `plotChromPeaks`: `character(1)` defining the
#'     y-axis label.
#'
#' @param ... Additional arguments passed to the `plot` (for
#'     `plotChromPeaks`) and `image` (for
#'     `plotChromPeakImage`) functions. Ignored for `add = TRUE`.
#'
#' @author Johannes Rainer
#'
#' @examples
#'
#' ## Load a test data set with detected peaks
#' faahko_sub <- loadXcmsData("faahko_sub2")
#'
#' ## plotChromPeakImage: plot an image for the identified peaks per file
#' plotChromPeakImage(faahko_sub)
#'
#' ## Show all detected chromatographic peaks from the first file
#' plotChromPeaks(faahko_sub)
#'
#' ## Plot all detected peaks from the second file and restrict the plot to a
#' ## mz-rt slice
#' plotChromPeaks(faahko_sub, file = 2, xlim = c(3500, 3600), ylim = c(400, 600))
plotChromPeaks <- function(x, file = 1, xlim = NULL, ylim = NULL,
                           add = FALSE, border = "#00000060", col = NA,
                           xlab = "retention time", ylab = "mz",
                           main = NULL, msLevel = 1L, ...) {
    if (!(is(x, "XCMSnExp") | inherits(x, "XcmsExperiment")))
        stop("'x' is supposed to be an 'XcmsExperiment' or 'XCMSnExp' object.")
    suppressMessages(
        x_file <- filterFile(x, file = file[1], keepAdjustedRtime = TRUE)
    )
    if (is.null(xlim))
        xlim <- range(rtime(x_file))
    if (is.null(ylim)) {
        if (is(x, "XcmsExperiment"))
            ylim <- range(unlist(mz(spectra(x_file)), use.names = FALSE))
        else ylim <- range(mz(x_file))
    }
    if (is.null(main))
        main <- basename(fileNames(x_file))
    pks <- chromPeaks(x_file, mz = ylim, rt = xlim, msLevel = msLevel)
    ## Initialize plot
    if (!add)
        plot(3, 3, pch = NA, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab,
             main = main, ...)
    if (nrow(pks))
        rect(xleft = pks[, "rtmin"], xright = pks[, "rtmax"],
             ybottom = pks[, "mzmin"], ytop = pks[, "mzmax"], col = col,
             border = border)
    invisible(TRUE)
}

#' @rdname plotChromPeaks
plotChromPeakImage <- function(x, binSize = 30, xlim = NULL, log = FALSE,
                               xlab = "retention time", yaxt = par("yaxt"),
                               main = "Chromatographic peak counts",
                               msLevel = 1L, ...) {
    if (!(is(x, "XCMSnExp") | inherits(x, "XcmsExperiment")))
        stop("'x' is supposed to be an 'XcmsExperiment' or 'XCMSnExp' object.")
    if (is.null(xlim))
        xlim <- c(floor(min(rtime(x))), ceiling(max(rtime(x))))
    brks <- seq(xlim[1], xlim[2], by = binSize)
    if (brks[length(brks)] < xlim[2])
        brks <- c(brks, brks[length(brks)] + binSize)
    pks <- chromPeaks(x, rt = xlim, msLevel = msLevel)
    if (nrow(pks)) {
        rts <- split(pks[, "rt"], as.factor(as.integer(pks[, "sample"])))
        cnts <- lapply(rts, function(z) {
            hst <- hist(z, breaks = brks, plot = FALSE)
            hst$counts
        })
        ## Add 0 vectors for samples in which no peaks were found.
        n_samples <- length(fileNames(x))
        sample_idxs <- 1:n_samples
        sample_idxs <- sample_idxs[!(as.character(sample_idxs) %in% names(rts))]
        if (length(sample_idxs)) {
            all_cnts <- vector("list", n_samples)
            all_cnts[as.numeric(names(cnts))] <- cnts
            zeros <- rep(0, (length(brks) - 1))
            all_cnts[sample_idxs] <- list(zeros)
            cnts <- all_cnts
        }
        cnts <- t(do.call(rbind, cnts))
        if (log)
            cnts <- log2(cnts)
        image(z = cnts, x = brks - (brks[2] - brks[1]) / 2, xaxs = "r",
              xlab = xlab, yaxt = "n", ...)
        if (yaxt != "n")
            axis(side = 2, at = seq(0, 1, length.out = n_samples),
                 labels = basename(fileNames(x)), las = 2)
    }
    invisible(TRUE)
}

#' @rdname XcmsExperiment
setMethod(
    "plot", c("MsExperiment", "missing"),
    function(x, y, msLevel = 1L, peakCol = "#ff000060", ...) {
        if (length(msLevel) > 1)
            warning("'plot' does support only a single MS level. ",
                    "Will use msLevel[1].")
        msLevel <- msLevel[1L]
        .xmse_plot_xic(x, msLevel = msLevel, peakCol = peakCol, ...)
    })

#' same as plot( type = "XIC") on a `XCMSnExp`. This supports both
#' `MsExperiment` and `XcmsExperiment` objects.
#'
#' @noRd
.xmse_plot_xic <- function(x, msLevel = 1L, peakCol = "#00000060", ...) {
    dots <- list(...)
    dots$type <- NULL
    if (any(names(dots) == "layout")) {
        if (!is.null(dots$layout))
            layout(layout)
        dots$layout <- NULL
    } else layout(MSnbase:::.vertical_sub_layout(length(x)))
    dev.hold()
    on.exit(dev.flush())
    peakCol <- peakCol[1L]
    if (inherits(x, "XcmsExperiment")) {
        if (hasAdjustedRtime(x))
            x <- applyAdjustedRtime(x)
        pkl <- chromPeaks(x, msLevel = msLevel)
        pkl <- split.data.frame(
            pkl, factor(pkl[, "sample"], levels = seq_along(x)))
        x <- as(x, "MsExperiment")
    } else pkl <- vector("list", length(x))
    fns <- basename(fileNames(x))
    for (i in seq_along(x)) {
        z <- x[i]
        flt <- filterMsLevel(spectra(z), msLevel = msLevel)
        lst <- as(flt, "list")
        lns <- lengths(lst) / 2
        lst <- do.call(rbind, lst)
        if (!length(lst))
            next
        df <- data.frame(rt = rep(rtime(flt), lns), lst)
        colnames(df)[colnames(df) == "intensity"] <- "i"
        do.call(MSnbase:::.plotXIC,
                c(list(x = df, main = fns[i], layout = NULL), dots))
        mzr <- range(df$mz, na.rm = TRUE)
        rtr <- range(df$rt, na.rm = TRUE)
        pks <- pkl[[i]]
        if (length(pks)) {
            pks <- pks[(pks[, "rtmax"] > rtr[1L] | pks[, "rtmin"] < rtr[2L]) &
                       (pks[, "mzmax"] > mzr[1L] | pks[, "mzmin"] < mzr[2L]),
                     , drop = FALSE]
            if (nrow(pks))
                do.call(rect, c(list(xleft = pks[, "rtmin"],
                                     ybottom = pks[, "mzmin"],
                                     xright = pks[, "rtmax"],
                                     ytop = pks[, "mzmax"],
                                     border = peakCol), dots))
        }
    }
}
