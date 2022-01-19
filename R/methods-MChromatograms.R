#' @rdname findChromPeaks-Chromatogram-CentWaveParam
#'
#' @aliases findChromPeaks-Chromatogram-CentWaveParam
#'
#' @examples
#'
#' ## Perform peak detection on an MChromatograms object
#' od3 <- readMSData(c(system.file("cdf/KO/ko15.CDF", package = "faahKO"),
#'     system.file("cdf/KO/ko16.CDF", package = "faahKO"),
#'     system.file("cdf/KO/ko18.CDF", package = "faahKO")),
#'     mode = "onDisk")
#'
#' ## Disable parallel processing for this example
#' register(SerialParam())
#'
#' ## Extract chromatograms for a m/z - retention time slice
#' chrs <- chromatogram(od3, mz = 344, rt = c(2500, 3500))
#'
#' ## Perform peak detection using CentWave
#' xchrs <- findChromPeaks(chrs, param = CentWaveParam())
#' xchrs
#'
#' ## Extract the identified chromatographic peaks
#' chromPeaks(xchrs)
#'
#' ## plot the result
#' plot(xchrs)
setMethod("findChromPeaks", signature(object = "MChromatograms",
                                      param = "CentWaveParam"),
          function(object, param, BPPARAM = bpparam(), ...) {
              .findChromPeaks_XChromatograms(object = object, param = param,
                                             BPPARAM = BPPARAM, ...)
          })

#' @rdname findChromPeaks-Chromatogram-CentWaveParam
setMethod("findChromPeaks", signature(object = "MChromatograms",
                                      param = "MatchedFilterParam"),
          function(object, param, BPPARAM = BPPARAM, ...) {
              .findChromPeaks_XChromatograms(object = object,
                                             param = param,
                                             BPPARAM = BPPARAM, ...)
          })

.findChromPeaks_XChromatograms <- function(object, param, BPPARAM, ...) {
    startDate <- date()
    if (missing(BPPARAM))
        BPPARAM <- bpparam()
    object <- as(object, "XChromatograms")
    object@.Data <- matrix(bplapply(c(object@.Data), FUN = findChromPeaks,
                                    param = param, BPPARAM = BPPARAM),
                           ncol = ncol(object),
                           dimnames = dimnames(object@.Data))
    ph_len <- length(object@.processHistory)
    if (ph_len && processType(object@.processHistory[[ph_len]]) ==
        .PROCSTEP.PEAK.DETECTION)
        object@.processHistory <-
            object@.processHistory[seq_len(ph_len - 1)]
    object <- addProcessHistory(
        object, XProcessHistory(param = param, date. = startDate,
                                type. = .PROCSTEP.PEAK.DETECTION,
                                fileIndex = seq_len(ncol(object))))
    if (validObject(object)) object
}

#' @rdname correlate-Chromatogram
setMethod("correlate",
          signature = c(x = "MChromatograms", y = "missing"),
          function(x, y = NULL, use = "pairwise.complete.obs",
                   method = c("pearson", "kendall", "spearman"),
                   align = c("closest", "approx"),
                   ...) {
              .Deprecated(new = "compareChromatograms")
              dots <- list(...)
              compareChromatograms(
                  x, x, ALIGNFUN = alignRt, FUN = cor,
                  ALIGNFUNARGS = c(list(method = align), dots),
                  FUNARGS = c(list(method = method, use = use), dots))
          })

#' @rdname correlate-Chromatogram
setMethod("correlate",
          signature = c(x = "MChromatograms", y = "MChromatograms"),
          function(x, y = NULL, use = "pairwise.complete.obs",
                   method = c("pearson", "kendall", "spearman"),
                   align = c("closest", "approx"), ...) {
              .Deprecated(new = "compareChromatograms")
              dots <- list(...)
              compareChromatograms(
                  x, y, ALIGNFUN = alignRt, FUN = cor,
                  ALIGNFUNARGS = c(list(method = align), dots),
                  FUNARGS = c(list(method = method, use = use), dots))
          })

#' @rdname removeIntensity-Chromatogram
setMethod("removeIntensity", "MChromatograms",
          function(object, which = "below_threshold", threshold = 0) {
              object@.Data <- matrix(lapply(c(object@.Data),
                                            FUN = removeIntensity,
                                            which = which,
                                            threshold = threshold),
                                     ncol = ncol(object),
                                     dimnames = dimnames(object@.Data))
              object
          })

#' @title Filtering sets of chromatographic data
#'
#' @aliases filterColumnsIntensityAbove filterColumnsKeepTop
#'
#' @rdname filter-MChromatograms
#'
#' @description
#'
#' These functions allow to filter (subset) [MChromatograms()] or
#' [XChromatograms()] objects, i.e. sets of chromatographic data, without
#' changing the data (intensity and retention times) within the individual
#' chromatograms ([Chromatogram()] objects).
#'
#' - `filterColumnsIntensityAbove`: subsets a `MChromatograms` objects keeping
#'   only columns (samples) for which `value` is larger than the provided
#'   `threshold` in `which` rows (i.e. if `which = "any"` a
#'   column is kept if **any** of the chromatograms in that column have a
#'   `value` larger than `threshold` or with `which = "all"` **all**
#'   chromatograms in that column fulfill this criteria). Parameter `value`
#'   allows to define on which value the comparison should be performed, with
#'   `value = "bpi"` the maximum intensity of each chromatogram is compared to
#'   `threshold`, with `value = "tic" the total sum of intensities of each
#'   chromatogram is compared to `threshold`. For `XChromatograms` object,
#'   `value = "maxo"` and `value = "into"` are supported which compares the
#'   largest intensity of all identified chromatographic peaks in the
#'   chromatogram with `threshold`, or the integrated peak area, respectively.
#'
#' - `filterColumnsKeepTop`: subsets a `MChromatograms` object keeping the top
#'   `n` columns sorted by the value specified with `sortBy`. In detail, for
#'   each column the value defined by `sortBy` is extracted from each
#'   chromatogram and aggregated using the `aggregationFun`. Thus, by default,
#'   for each chromatogram the maximum intensity is determined
#'   (`sortBy = "bpi"`) and these values are summed up for chromatograms in the
#'   same column (`aggregationFun = sum`). The columns are then sorted by these
#'   values and the top `n` columns are retained in the returned
#'   `MChromatograms`. Similar to the `filterColumnsIntensityAbove` function,
#'   this function allows to use for `XChromatograms` objects to sort the
#'   columns by column `sortBy = "maxo"` or `sortBy = "into"` of the
#'   `chromPeaks` matrix.
#'
#' @param aggregationFun for `filterColumnsKeepTop`: function to be used to
#'     aggregate (combine) the values from all chromatograms in each column.
#'     Defaults to `aggregationFun = sum` in which case the sum of the values
#'     is used to rank the columns. Alternatively the `mean`, `median` or
#'     similar function can be used.
#'
#' @param n for `filterColumnsKeepTop`: `integer(1)` specifying the number of
#'     columns that should be returned. `n` will be rounded to the closest
#'     (larger) integer value.
#'
#' @param object [MChromatograms()] or [XChromatograms()] object.
#'
#' @param sortBy for `filterColumnsKeepTop`: the value by which columns should
#'     be ordered to determine the top n columns. Can be either `sortBy = "bpi"`
#'     (the default), in which case the maximum intensity of each column's
#'     chromatograms is used, or `sortBy = "tic"` to use the total intensity
#'     sum of all chromatograms. For [XChromatograms()] objects also
#'     `value = "maxo"` and `value = "into"` is supported to use the maximum
#'     intensity or the integrated area of identified chromatographic peaks
#'     in each chromatogram.
#'
#' @param threshold for `filterColumnsIntensityAbove`: `numeric(1)` with the
#'     threshold value to compare against.
#'
#' @param value `character(1)` defining which value should be used in the
#'     comparison or sorting. Can be `value = "bpi"` (default) to use the
#'     maximum intensity per chromatogram or `value = "tic"` to use the sum
#'     of intensities per chromatogram. For [XChromatograms()] objects also
#'     `value = "maxo"` and `value = "into"` is supported to use the maximum
#'     intensity or the integrated area of identified chromatographic peaks
#'     in each chromatogram.
#'
#' @param which for `filterColumnsIntensityAbove`: `character(1)` defining
#'     whether **any** (`which = "any"`, default) or **all** (`which = "all"`)
#'     chromatograms in a column have to fulfill the criteria for the column
#'     to be kept.
#'
#' @return a filtered `MChromatograms` (or `XChromatograms`) object with the
#'     same number of rows (EICs) but eventually a lower number of columns
#'     (samples).
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @examples
#'
#' chr1 <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
#'     intensity = c(5, 29, 50, NA, 100, 12, 3, 4, 1, 3))
#' chr2 <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
#'     intensity = c(80, 50, 20, 10, 9, 4, 3, 4, 1, 3))
#' chr3 <- Chromatogram(rtime = 3:9 + rnorm(7, sd = 0.3),
#'     intensity = c(53, 80, 130, 15, 5, 3, 2))
#'
#' chrs <- MChromatograms(list(chr1, chr2, chr1, chr3, chr2, chr3),
#'     ncol = 3, byrow = FALSE)
#' chrs
#'
#' #### filterColumnsIntensityAbove
#' ##
#' ## Keep all columns with for which the maximum intensity of any of its
#' ## chromatograms is larger 90
#' filterColumnsIntensityAbove(chrs, threshold = 90)
#'
#' ## Require that ALL chromatograms in a column have a value larger 90
#' filterColumnsIntensityAbove(chrs, threshold = 90, which = "all")
#'
#' ## If none of the columns fulfills the criteria no columns are returned
#' filterColumnsIntensityAbove(chrs, threshold = 900)
#'
#' ## Filtering XChromatograms allow in addition to filter on the columns
#' ## "maxo" or "into" of the identified chromatographic peaks within each
#' ## chromatogram.
#'
#' #### filterColumnsKeepTop
#' ##
#' ## Keep the 2 columns with the highest sum of maximal intensities in their
#' ## chromatograms
#' filterColumnsKeepTop(chrs, n = 1)
#'
#' ## Keep the 50 percent of columns with the highest total sum of signal. Note
#' ## that n will be rounded to the next larger integer value
#' filterColumnsKeepTop(chrs, n = 0.5 * ncol(chrs), sortBy = "tic")
setMethod("filterColumnsIntensityAbove", "MChromatograms",
          function(object, threshold = 0, value = c("bpi", "tic"),
                   which = c("any", "all")) {
              value <- match.arg(value)
              which <- match.arg(which)
              if (length(threshold) > 1 || !is.numeric(threshold))
                  stop("'threshold' should be a 'numeric' of length 1")
              nc <- ncol(object)
              if (value == "bpi")
                  FUN <- max
              else FUN <- sum
              which_fun <- getFunction(which)
              keep <- rep(FALSE, nc)
              for (i in seq_len(nc)) {
                  vals <- vapply(object[, i], function(z) {
                      FUN(z@intensity, na.rm = TRUE)
                  }, FUN.VALUE = NA_real_, USE.NAMES = FALSE)
                  keep[i] <- which_fun(vals > threshold)
              }
              object[, keep]
          })

#' @rdname filter-MChromatograms
setMethod("filterColumnsKeepTop", "MChromatograms",
          function(object, n = 1L, sortBy = c("bpi", "tic"),
                   aggregationFun = sum) {
              sortBy <- match.arg(sortBy)
              if (length(n) > 1 || !is.numeric(n))
                  stop("'n' should be an 'integer' of length 1")
              n <- ceiling(n)
              nc <- ncol(object)
              nr <- nrow(object)
              if (n > nc)
                  stop("'n' should be smaller or equal than the number of ",
                       "columns (", nc, ")")
              if (sortBy == "bpi")
                  FUN <- max
              else FUN <- sum
              colval <- numeric(nc)
              for (i in seq_len(nc)) {
                  if (nr == 1)
                      vals <- FUN(object[1, i]@intensity, na.rm = TRUE)
                  else
                      vals <- vapply(object[, i], function(z) {
                          FUN(z@intensity, na.rm = TRUE)
                      }, FUN.VALUE = NA_real_, USE.NAMES = FALSE)
                  colval[i] <- aggregationFun(vals, na.rm = TRUE)
              }
              idx <- order(colval, decreasing = TRUE)[seq_len(n)]
              object[, sort(idx)]
          })

#' @title Plot multiple chromatograms into the same plot
#'
#' @aliases plotChromatogramsOverlay plotChromatogramsOverlay,MChromatograms-method
#'
#' @description
#'
#' `plotOverlay` draws chromatographic peak data from multiple (different)
#' extracted ion chromatograms (EICs) into the same plot. This allows to
#' directly compare the peak shape of these EICs in the same sample. In
#' contrast to the `plot` function for [MChromatograms()] object, which draws
#' the data from the same EIC across multiple samples in the same plot, this
#' function draws the different EICs from the same sample into the same plot.
#'
#' If `plotChromatogramsOverlay` is called on a `XChromatograms` object any
#' present chromatographic peaks will also be highlighted/drawn depending on the
#' parameters `peakType`, `peakCol`, `peakBg` and `peakPch` (see also help on
#' the `plot` function for `XChromatogram()` object for details).
#'
#' @param col definition of the color in which the chromatograms should be
#'     drawn. Can be of length 1 or equal to `nrow(object)` to plot each
#'     overlayed chromatogram in a different color.
#'
#' @param main optional title of the plot. If not defined, the range of m/z
#'     values is used.
#'
#' @param object [MChromatograms()] or [XChromatograms()] object.
#'
#' @param peakBg if `object` is a `XChromatograms` object: definition of
#'     background color(s) for each chromatographic peak. Has to be either of
#'     length 1 or equal to the number of peaks in `object`. If not specified,
#'     the peak will be drawn in the color defined by `col`.
#'
#' @param peakCol if `object` is a `XChromatograms` object: definition of
#'     color(s) for each chromatographic peak. Has to be either of length 1 or
#'     equal to the number of peaks in `object`. If not specified, the peak will
#'     be drawn in the color defined by `col`.
#'
#' @param peakPch if `object` is a `XChromatograms` object: *point character* to
#'     be used to label the apex position of the chromatographic peak if
#'     `peakType = "point"`.
#'
#' @param peakType if `object` is a `XChromatograms` object: how chromatographic
#'     peaks should be drawn: `peakType = "polygon"` (the default): label the
#'     full chromatographic peak area, `peakType = "rectangle"`: indicate the
#'     chromatographic peak by a rectangle and `peakType  = "point"`: label the
#'     chromatographic peaks' apex position with a point.
#'
#' @param stacked `numeric(1)` defining the part (proportion) of the y-axis to
#'     use to *stack* EICs depending on their m/z values. If `stacked = 0` (the
#'     default) no stacking is performed. With `stacked = 1` half of the y-axis
#'     is used for stacking and half for the intensity y-axis (i.e. the ratio
#'     between stacking and intensity y-axis is 1:1). Note that if `stacking`
#'     is different from 0 no y-axis and label are drawn.
#'
#' @param transform `function` to transform the intensity values before
#'     plotting. Defaults to `transform = identity` which plots the data as it
#'     is. With `transform = log10` intensity values would be log10 transformed
#'     before plotting.
#'
#' @param type `character(1)` defing the type of the plot. By default
#'     (`type = "l"`) each chromatogram is drawn as a line.
#'
#' @param xlab `character(1)` defining the x-axis label.
#'
#' @param xlim optional `numeric(2)` defining the x-axis limits.
#'
#' @param ylab `character(1)` defining the y-axis label.
#'
#' @param ylim optional `numeric(2)` defining the y-axis limits.
#'
#' @param ... optional arguments to be passed to the plotting functions (see
#'     help on the base R `plot` function.
#'
#' @return silently returns a `list` (length equal to `ncol(object)` of
#'     `numeric` (length equal to `nrow(object)`) with the y position of
#'     each EIC.
#'
#' @md
#'
#' @author Johannes Rainer
#'
#' @name plotChromatogramsOverlay
#'
#' @examples
#'
#' ## Load preprocessed data and extract EICs for some features.
#' library(xcms)
#' data(xdata)
#' ## Update the path to the files for the local system
#' dirname(xdata) <- c(rep(system.file("cdf", "KO", package = "faahKO"), 4),
#'                     rep(system.file("cdf", "WT", package = "faahKO"), 4))
#' ## Subset to the first 3 files.
#' xdata <- filterFile(xdata, 1:3, keepFeatures = TRUE)
#'
#' ## Define features for which to extract EICs
#' fts <- c("FT097", "FT163", "FT165")
#' chrs <- featureChromatograms(xdata, features = fts)
#'
#' plotChromatogramsOverlay(chrs)
#'
#' ## plot the overlay of EICs in the first sample
#' plotChromatogramsOverlay(chrs[, 1])
#'
#' ## Define a different color for each feature (row in chrs). By default, also
#' ## all chromatographic peaks of a feature is labeled in the same color.
#' plotChromatogramsOverlay(chrs[, 1],
#'     col = c("#ff000040", "#00ff0040", "#0000ff40"))
#'
#' ## Alternatively, we can define a color for each individual chromatographic
#' ## peak and provide this with the `peakBg` and `peakCol` parameters.
#' chromPeaks(chrs[, 1])
#'
#' ## Use a color for each of the two identified peaks in that sample
#' plotChromatogramsOverlay(chrs[, 1],
#'     col = c("#ff000040", "#00ff0040", "#0000ff40"),
#'     peakBg = c("#ffff0020", "#00ffff20"))
#'
#' ## Plotting the data in all samples.
#' plotChromatogramsOverlay(chrs,
#'     col = c("#ff000040", "#00ff0040", "#0000ff40"))
#'
#' ## Creating a "stacked" EIC plot: the EICs are placed along the y-axis
#' ## relative to their m/z value. With `stacked = 1` the y-axis is split in
#' ## half, the lower half being used for the stacking of the EICs, the upper
#' ## half being used for the *original* intensity axis.
#' res <- plotChromatogramsOverlay(chrs[, 1], stacked = 1,
#'     col = c("#ff000040", "#00ff0040", "#0000ff40"))
#' ## add horizontal lines for the m/z values of each EIC
#' abline(h = res[[1]], col = "grey", lty = 2)
#'
#' ## Note that this type of visualization is different than the conventional
#' ## plot function for chromatographic data, which will draw the EICs for
#' ## multiple samples into the same plot
#' plot(chrs)
#'
#' ## Converting the object to a MChromatograms without detected peaks
#' chrs <- as(chrs, "MChromatograms")
#'
#' plotChromatogramsOverlay(chrs,
#'     col = c("#ff000040", "#00ff0040", "#0000ff40"))
setMethod("plotChromatogramsOverlay", "MChromatograms",
          function(object, col = "#00000060", type = "l", main = NULL,
                   xlab = "rtime", ylab = "intensity", xlim = numeric(),
                   ylim = numeric(), stacked = 0, transform = identity, ...) {
              nsam <- ncol(object)
              transform <- match.fun(transform)
              if (nsam > 1)
                  par(mfrow = n2mfrow(nsam, 1))
              res <- vector("list", nsam)
              for (i in seq_len(nsam)) {
                  res[[i]] <- .plot_xchromatograms_overlay(
                      object[, i, drop = FALSE], main = main,
                      xlab = xlab, ylab = ylab, xlim = xlim,
                      ylim = ylim, col = col, type = type,
                      stacked = stacked, transform = transform, ...)
              }
              invisible(res)
          })

#' @rdname plotChromatogramsOverlay
setMethod("plotChromatogramsOverlay", "XChromatograms",
          function(object, col = "#00000060", type = "l", main = NULL,
                   xlab = "rtime", ylab = "intensity", xlim = numeric(),
                   ylim = numeric(), peakType = c("polygon", "point",
                                                  "rectangle", "none"),
                   peakBg = NULL, peakCol = NULL, peakPch = 1,
                   stacked = 0, transform = identity, ...) {
              transform <- match.fun(transform)
              nsam <- ncol(object)
              peakType <- match.arg(peakType)
              if (nsam > 1)
                  par(mfrow = n2mfrow(nsam, 1))
              res <- vector("list", nsam)
              for (i in seq_len(nsam)) {
                  res[[i]] <- .plot_xchromatograms_overlay(
                      object[, i, drop = FALSE], main = main, xlab = xlab,
                      ylab = ylab, xlim = xlim, ylim = ylim, col = col,
                      type = type, peakType = peakType, peakCol = peakCol,
                      peakBg = peakBg, peakPch = peakPch,
                      stacked = stacked, transform = transform, ...)
              }
              invisible(res)
          })

.plot_single_xchromatograms <- function(x, type = "l", col = "#00000060",
                                        peakType = c("polygon", "point",
                                                     "rectangle", "none"),
                                        peakCol = NULL, peakBg = NULL,
                                        peakPch = 1, yoffset = 0,
                                        fill = NA, transform = identity, ...) {
    peakType <- match.arg(peakType)
    .plot_single_chromatograms(x, type = type, col = col,
                               yoffset = yoffset, fill = fill,
                               transform = transform, ...)
    pks <- chromPeaks(x)
    if (nrow(pks) && peakType != "none") {
        .add_chromatogram_peaks(as(x, "Chromatogram"), pks, col = peakCol,
                                bg = peakBg, type = peakType,
                                pch = peakPch, yoffset = yoffset,
                                transform = transform, ...)
    }
}

.plot_single_chromatograms <- function(x, type = "l", col = "#00000060",
                                       fill = NA, yoffset = 0,
                                       transform = identity, ...) {
    ints <- transform(intensity(x)) + yoffset
    ints[is.infinite(ints)] <- 0
    if (!is.na(fill)) {
        nnas <- !is.na(ints)
        rts <- rtime(x)[nnas]
        if (any(nnas))
            polygon(c(rts[1], rts, rts[length(rts)]),
                    c(yoffset, ints[nnas], yoffset), border = NA, col = fill)
    }
    plot.xy(xy.coords(rtime(x), ints),
            type = type, col = col, ...)
}

.plot_xchromatograms_overlay <- function(x, xlab = "rtime", ylab = "intensity",
                                         type = "l", col = "#00000060",
                                         xlim = numeric(), ylim = numeric(),
                                         main = NULL, axes = TRUE,
                                         frame.plot = axes,
                                         peakType = c("polygon",
                                                      "point",
                                                      "rectangle",
                                                      "none"),
                                         peakCol = NULL,
                                         peakBg = NULL,
                                         peakPch = 1,
                                         fill = NA,
                                         yoffset = 0,
                                         stacked = 0,
                                         transform = identity, ...) {
    if (ncol(x) > 1)
        stop(".plot_chromatograms_overlay supports only single column",
             " XChromatograms")
    peakType <- match.arg(peakType)
    stacked <- stacked[1L]
    nchr <- nrow(x)
    if (length(col) != nchr)
        col <- rep(col[1], nchr)
    if (length(fill) != nchr)
        fill <- rep(fill[1], nchr)
    if(!length(xlim))
        xlim <- suppressWarnings(range(lapply(x, rtime), na.rm = TRUE))
    if(!length(ylim)) {
        ylim <- transform(range(c(lapply(x, intensity), 0), na.rm = TRUE))
        ylim[is.infinite(ylim)] <- 0
    }
    if (any(is.infinite(xlim)))
        xlim <- c(0, 0)
    if (any(is.infinite(ylim)))
        ylim <- c(0, 0)
    ## yposition and stacking
    if (stacked != 0) {
        mzs <- mz(x)
        if (any(is.na(mzs)))
            mzs <- seq_len(nchr)
        else mzs <- rowMeans(mzs)
        ymz_range <- stacked * ylim
        mzr <- range(mzs)
        ypos <- ymz_range[2] + (mzs - mzr[2]) * (ymz_range[2] - ymz_range[1]) /
            (mzr[2] - mzr[1])
        ylim[2] <- ylim[2] + ymz_range[2]
        ylab <- ""
    } else
        ypos <- rep(yoffset, nchr)
    dev.hold()
    on.exit(dev.flush())
    plot.new()
    plot.window(xlim = xlim, ylim = ylim)
    if (axes) {
        axis(side = 1, ...)
        if (stacked == 0)
            axis(side = 2, ...)
    }
    if (frame.plot)
        box(...)
    if (!length(main)) {
        mzs <- mz(x)
        mzs <- c(suppressWarnings(min(mzs[, 1], na.rm = TRUE)),
                 suppressWarnings(max(mzs[, 2], na.rm = TRUE)))
        main <- paste0("m/z: ", format(mzs[1], digits = 6), " - ",
                       format(mzs[2], digits = 6))
    }
    title(main = main, xlab = xlab, ylab = ylab, ...)
    ## Define colors for peaks - but only if it IS a XChromatograms.
    if (is(x, "XChromatograms") && peakType != "none") {
        pk_row <- chromPeaks(x)[, "row"]
        np <- length(pk_row)
        f <- factor(pk_row, levels = seq_len(nchr))
        if (!length(peakCol))
            peakCol <- col[f]
        if (!length(peakBg))
            peakBg <- col[f]
        if (length(peakCol) != np)
            peakCol <- rep(peakCol[1L], np)
        if (length(peakBg) != np)
            peakBg <- rep(peakBg[1L], np)
        for (i in order(ypos, decreasing = TRUE)) {
            .plot_single_xchromatograms(x[i, 1], col = col[i], type = type,
                                        peakType = peakType,
                                        peakCol = peakCol[pk_row == i],
                                        peakBg = peakBg[pk_row == i],
                                        peakPch = peakPch, fill = fill[i],
                                        yoffset = ypos[i],
                                        transform = transform, ...)
        }
    } else {
        for (i in order(ypos, decreasing = TRUE)) {
            .plot_single_chromatograms(x[i, 1], col = col[i], type = type,
                                       fill = fill[i], yoffset = ypos[i],
                                       transform = transform, ...)
        }
    }
    ypos
}
