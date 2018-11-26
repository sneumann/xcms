#' @rdname XChromatogram
#'
#' @md
setMethod("show", "XChromatogram", function(object) {
    callNextMethod()
    cat("Identified chromatographic peaks (", nrow(object@chromPeaks),"):\n",
        sep = "")
    cat(" rt\trtmin\trtmax\tinto\tmaxo\tsn\n")
    for (i in seq_len(nrow(object@chromPeaks)))
        cat(" ", paste(object@chromPeaks[i, ], collapse = "\t"), "\n", sep = "")
})

#' @rdname XChromatogram
#'
#'
#' @section Accessing data:
#'
#' - `chromPeaks`, `chromPeaks<-`: extract or set the matrix with the
#'   chromatographic peak definitions. Parameter `rt` allows to specify a
#'   retention time range for which peaks should be returned along with
#'   parameter `type` that defines how *overlapping* is defined (parameter
#'   description for details).
#'
#' - `hasChromPeaks`: infer whether a `XChromatogram` (or `XChromatograms`)
#'   has chromatographic peaks. For `XChromatogram`: returns a `logical(1)`,
#'   for `XChromatograms`: returns a `matrix`, same dimensions than `object`
#'   with either `true` or `false` if chromatographic peaks are available in
#'   the chromatogram at the respective position.
#'
#' @param rt For `chromPeaks`: `numeric(2)` defining the retention time range
#'     for which chromatographic peaks should be returned.
#'     For `filterRt`: `numeric(2)` defining the retention time range to
#'     reduce `object` to.
#'
#' @param ppm For `chromPeaks`: `numeric(1)` defining a ppm to expand the
#'     provided m/z range
#'
#' @param type For `chromPeaks`: `character(1)` defining which peaks to return
#'     if `rt` is provided: `"any"` (default) return all peaks that are even
#'     partially overlapping with `rt`, `"within"` return peaks that are
#'     completely within `rt` and `"apex_within"` return peaks which apex
#'     is within `rt`.
#'     For `plot`: what type of plot should be used for the
#'     chromatogram (such as `"l"` for lines, `"p"` for points etc), see help
#'     of [plot()] in the `graphics` package for more details.
#'
#' @param value For `chromPeaks<-`: a numeric `matrix` with required columns
#'     `"rt"`, `"rtmin"`, `"rtmax"`, `"into"` and `"maxo"`.
#'
#' @md
#'
#' @examples
#'
#' ## Extract the chromatographic peaks
#' chromPeaks(xchr)
setMethod("chromPeaks", "XChromatogram", function(object, rt = numeric(),
                                                  mz = numeric(), ppm = 0,
                                                  type = c("any", "within",
                                                           "apex_within")) {
    type <- match.arg(type)
    pks <- object@chromPeaks
    if (length(rt) && nrow(pks)) {
        rt <- range(rt)
        pks <- switch(type,
                      any = pks[which(pks[, "rtmin"] <= rt[2] &
                                      pks[, "rtmax"] >= rt[1]), , drop = FALSE],
                      within = pks[which(pks[, "rtmin"] >= rt[1] &
                                         pks[, "rtmax"] <= rt[2]), ,
                                   drop = FALSE],
                      apex_within = pks[which(pks[, "rt"] >= rt[1] &
                                              pks[, "rt"] <= rt[2]), ,
                                        drop = FALSE]
                      )
    }
    if (length(mz) && nrow(pks) & all(c("mz", "mzmin", "mzmax")
                                      %in% colnames(pks))) {
        mz <- .ppm_range(mz, ppm = ppm)
        pks <- switch(type,
                      any = pks[which(pks[, "mzmin"] <= mz[2] &
                                      pks[, "mzmax"] >= mz[1]), , drop = FALSE],
                      within = pks[which(pks[, "mzmin"] >= mz[1] &
                                         pks[, "mzmax"] <= mz[2]), ,
                                   drop = FALSE],
                      apex_within = pks[which(pks[, "mz"] >= mz[1] &
                                              pks[, "mz"] <= mz[2]), ,
                                        drop = FALSE]
                      )
    }
    pks
})


#' @rdname XChromatogram
setReplaceMethod("chromPeaks", "XChromatogram", function(object, value) {
    if (!is.matrix(value))
        stop("'value' should be a numeric matrix")
    object@chromPeaks <- value
    if (validObject(object)) object
})

#' @rdname XChromatogram
#'
#' @section Plotting and visualizing:
#'
#' - `plot` draws the chromatogram and highlights in addition any
#'   chromatographic peaks present in the `XChromatogram` (unless
#'   `peakType = "none"` was specified). To draw peaks in different colors
#'   a vector of color definitions with length equal to `nrow(chromPeaks(x))`
#'   has to be submitted  with `peakCol` and/or `peakBg` defining one color
#'   for each peak (in the order as peaks are in `chromPeaks(x))`. For base
#'   peak chromatograms or total ion chromatograms it might be better to set
#'   `peakType = "none"` to avoid generating busy plots.
#'
#' @note
#'
#' Highlighting the peak area(s) in an `XChromatogram` object (`plot` with
#' `peakType = "polygon"`) draws a polygon representing the shown chromatogram
#' from the peak's minimal retention time to the maximal retention time. If the
#' `XChromatogram` was extracted from an [XCMSnExp()] object with the
#' [chromatogram()] function this might therefore not represent the actual
#' identified peak area if the m/z range that was used to extract the
#' chromatogram was larger than the peak's m/z.
#'
#' @param col For `plot`: the color to be used to draw the chromatogram.
#'
#' @param lty For `plot`: the line type.
#'
#' @param xlab For `plot`: the x axis label.
#'
#' @param ylab For `plot`: the y axis label.
#'
#' @param main For `plot`: an optional title for the plot.
#'
#' @param peakType For `plot`: `character(1)` defining how (and if) identified
#'     chromatographic peak within the chromatogram should be plotted. Options
#'     are `"polygon"` (default): draw the peak borders with the `peakCol` color
#'     and fill the peak area with the `peakBg` color, `"point"`: indicate the
#'     peak's apex with a point, `"rectangle"`: draw a rectangle around the
#'     identified peak and `"none"`: don't draw peaks.
#'
#' @param peakCol For `plot`: the foreground color for the peaks.
#'     For `peakType = "polygon"` and `peakType = "rectangle"` this is the
#'     color for the border. Use `NA` to not use a foreground color. This
#'     should either be a single color or a vector of colors with the same
#'     length than `chromPeaks(x)` has rows.
#'
#' @param peakBg For `plot`: the background color for the peaks. For
#'     `peakType = "polygon"` and `peakType = "rectangle"` the peak are or
#'     rectangle will be filled with this color. Use `NA` to skip. This should
#'     be either a single color or a vector of colors with the same length than
#'     `chromPeaks(x)` has rows.
#'
#' @param peakPch For `plot`: the point character to be used for
#'     `peakType = "point"`. See [plot()] in the `graphics` package for more
#'     details.
#'
#' @md
#'
#' @examples
#'
#' ## Plotting of a single XChromatogram object
#' ## o Don't highlight chromatographic peaks
#' plot(xchr, peakType = "none")
#'
#' ## o Indicate peaks with a polygon
#' plot(xchr)
#'
#' ## Add a second peak to the data.
#' pks <- rbind(chromPeaks(xchr), c(7, 7, 10, NA, 15, NA))
#' chromPeaks(xchr) <- pks
#'
#' ## Plot the peaks in different colors
#' plot(xchr, peakCol = c("#ff000080", "#0000ff80"),
#'     peakBg = c("#ff000020", "#0000ff20"))
#'
#' ## Indicate the peaks as rectangles
#' plot(xchr, peakCol = c("#ff000060", "#0000ff60"), peakBg = NA,
#'     peakType = "rectangle")
setMethod("plot", "XChromatogram", function(x, col = "#00000060", lty = 1,
                                            type = "l",
                                            xlab = "retention time",
                                            ylab = "intensity",
                                            main = NULL,
                                            peakType = c("polygon",
                                                         "point",
                                                         "rectangle",
                                                         "none"),
                                            peakCol = "#00000060",
                                            peakBg = "#00000020",
                                            peakPch = 1, ...) {
    peakType <- match.arg(peakType)
    callNextMethod(x = x, col = col, lty = lty, type = type, xlab = xlab,
                   ylab = ylab, main = main, ...)
    pks <- chromPeaks(x)
    nr <- nrow(pks)
    if (nr && peakType != "none") {
        if (length(peakCol) != nr)
            peakCol <- rep(peakCol[1], nr)
        if (length(peakBg) != nr)
            peakBg <- rep(peakBg[1], nr)
        if (length(peakPch) != nr)
            peakPch <- rep(peakPch[1], nr)
        suppressWarnings(.add_chromatogram_peaks(x, pks, col = peakCol,
                                                 bg = peakBg, type = peakType,
                                                 pch = peakPch, ...))
    }
})

.add_chromatogram_peaks <- function(x, pks, col, bg, type, pch, ...) {
    switch(type,
           point = {
               points(pks[, "rt"], pks[, "maxo"], pch = pch, col = col,
                      bg = bg, ...)
           },
           rectangle = {
               rect(xleft = pks[, "rtmin"], xright = pks[, "rtmax"],
                    ybottom = rep(0, nrow(pks)), ytop = pks[, "maxo"],
                    col = bg, border = col, ...)
           },
           polygon = {
               ordr <- order(pks[, "maxo"], decreasing = TRUE)
               pks <- pks[ordr, , drop = FALSE]
               col <- col[ordr]
               bg <- bg[ordr]
               xs_all <- numeric()
               ys_all <- numeric()
               for (i in seq_len(nrow(pks))) {
                   chr <- filterRt(x, rt = pks[i, c("rtmin", "rtmax")])
                   xs <- rtime(chr)
                   if (!length(xs)) {
                       next
                       col <- col[-i]
                       bg <- bg[-i]
                   }
                   xs <- c(xs[1], xs, xs[length(xs)])
                   ys <- c(0, intensity(chr), 0)
                   nona <- !is.na(ys)
                   if (length(xs_all)) {
                       xs_all <- c(xs_all, NA)
                       ys_all <- c(ys_all, NA)
                   }
                   xs_all <- c(xs_all, xs[nona])
                   ys_all <- c(ys_all, ys[nona])
                   ## polygon(xs[nona], ys[nona], border = col[i], col = bg[i],
                   ##         ...)
               }
               polygon(xs_all, ys_all, border = col, col = bg, ...)
           })
}

#' @rdname XChromatogram
#'
#' @section Filtering and subsetting:
#'
#' - `filterMz` filters the chromatographic peaks within an `XChromatogram`, if
#'   a column `"mz"` is present in the `chromPeaks` matrix. This would be the
#'   case if the `XChromatogram` was extracted from an [XCMSnExp()] object with
#'   the [chromatogram()] function. All chromatographic peaks with their m/z
#'   within the m/z range defined by `mz` will be retained.
#'
#' - `filterRt` filters the chromatogram by the provided retention time range.
#'   All eventually present chromatographic peaks with their apex within the
#'   retention time range specified with `rt` will be retained.
#'
#' @md
setMethod("filterMz", "XChromatogram", function(object, mz, ...) {
    if (missing(mz) || length(mz) == 0)
        return(object)
    pks <- chromPeaks(object)
    if (nrow(pks) && any(colnames(pks) == "mz")) {
        mz <- range(mz)
        pks <- pks[which(pks[, "mz"] >= mz[1] & pks[, "mz"] <= mz[2]), ,
                   drop = FALSE]
        chromPeaks(object) <- pks
    }
    if (validObject(object)) object
})

#' @rdname XChromatogram
#'
#' @md
#'
#' @examples
#'
#' ## Filter the XChromatogram by retention time
#' xchr_sub <- filterRt(xchr, rt = c(4, 6))
#' xchr_sub
#' plot(xchr_sub)
setMethod("filterRt", "XChromatogram", function(object, rt, ...) {
    if (missing(rt) || length(rt) == 0) return(object)
    pks <- chromPeaks(object)
    object <- callNextMethod()
    if (nrow(pks)) {
        rt <- range(rt)
        pks <- pks[which(pks[, "rt"] >= rt[1] & pks[, "rt"] <= rt[2]), ,
                   drop = FALSE]
        chromPeaks(object) <- pks
    }
    if (validObject(object)) object
})

#' @rdname XChromatogram
#'
#' @md
setMethod("hasChromPeaks", "XChromatogram", function(object) {
    as.logical(nrow(object@chromPeaks))
})
