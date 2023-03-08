## Plotting functions for XcmsExperiment.

## plot: don't think that would be possible... has lower priority, maybe just
## provide a XIC plot.

## plotChromPeakDensity: modify functions-XCMSnExp/.plotChromPeakDensity
## Wait - that's deprecated! Use that on an XChromatograms instead!

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
#' data(faahko_sub)
#' ## Update the path to the files for the local system
#' dirname(faahko_sub) <- system.file("cdf/KO", package = "faahKO")
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



## plotChromPeaks: modify functions-XCMSnExp/plotChromPeaks

## plotChromPeakImage: modify functions-XCMSnExp/plotChromPeakImage
