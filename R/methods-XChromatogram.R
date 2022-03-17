#' @rdname XChromatogram
#'
#' @md
setMethod("show", "XChromatogram", function(object) {
    callNextMethod()
    cat("Identified chromatographic peaks (", nrow(object@chromPeaks),"):\n",
        sep = "")
    cat(" rt\trtmin\trtmax\tinto\tmaxo\tsn ")
    nc <- ncol(object@chromPeaks)
    if (nc > 6)
        cat("(", nc - 6, " more column(s))", sep = "")
    cat("\n")
    for (i in seq_len(nrow(object@chromPeaks)))
        cat(" ", paste(object@chromPeaks[i, .CHROMPEAKS_REQ_NAMES],
                       collapse = "\t"), "\n", sep = "")
})

#' @rdname XChromatogram
#'
#' @aliases filterChromPeaks
#'
#' @section Accessing data:
#'
#' See also help of [Chromatogram] in the `MSnbase` package for general
#' information and data access. The methods listed here are specific for
#' `XChromatogram` and `XChromatograms` objects.
#'
#' - `chromPeaks`, `chromPeaks<-`: extract or set the matrix with the
#'   chromatographic peak definitions. Parameter `rt` allows to specify a
#'   retention time range for which peaks should be returned along with
#'   parameter `type` that defines how *overlapping* is defined (parameter
#'   description for details). For `XChromatogram` objects the function returns
#'   a `matrix` with columns `"rt"` (retention time of the peak apex),
#'   `"rtmin"` (the lower peak boundary), `"rtmax"` (the upper peak boundary),
#'   `"into"` (the ingegrated peak signal/area of the peak), `"maxo"` (the
#'   maximum instensity of the peak and `"sn"` (the signal to noise ratio).
#'   Note that, depending on the peak detection algorithm, the matrix may
#'   contain additional columns.
#'   For `XChromatograms` objects the `matrix` contains also columns `"row"`
#'   and `"column"` specifying in which chromatogram of `object` the peak was
#'   identified. Chromatographic peaks are ordered by row.
#'
#' - `chromPeakData`, `chromPeakData<-`: extract or set the [DataFrame()] with
#'   optional chromatographic peak annotations.
#'
#' - `hasChromPeaks`: infer whether a `XChromatogram` (or `XChromatograms`)
#'   has chromatographic peaks. For `XChromatogram`: returns a `logical(1)`,
#'   for `XChromatograms`: returns a `matrix`, same dimensions than `object`
#'   with either `TRUE` or `FALSE` if chromatographic peaks are available in
#'   the chromatogram at the respective position.
#'
#' - `hasFilledChromPeaks`: whether a `XChromatogram` (or a `XChromatogram` in
#'   a `XChromatograms`) has filled-in chromatographic peaks.
#'   For `XChromatogram`: returns a `logical(1)`,
#'   for `XChromatograms`: returns a `matrix`, same dimensions than `object`
#'   with either `TRUE` or `FALSE` if chromatographic peaks are available in
#'   the chromatogram at the respective position.
#'
#' - `dropFilledChromPeaks`: removes filled-in chromatographic peaks. See
#'   [dropFilledChromPeaks()] help for [XCMSnExp()] objects for more
#'   information.
#'
#' - `hasFeatures`: for `XChromatograms` objects only: if correspondence
#'   analysis has been performed and m/z-rt feature definitions are present.
#'   Returns a `logical(1)`.
#'
#' - `dropFeatureDefinitions`: for `XChrmomatograms` objects only: delete any
#'   correspondence analysis results (and related process history).
#'
#' - `featureDefinitions`: for `XChromatograms` objects only. Extract the
#'   results from the correspondence analysis (performed with
#'   `groupChromPeaks`). Returns a `DataFrame` with the properties of the
#'   defined m/z-rt features: their m/z and retention time range. Columns
#'   `peakidx` and `row` contain the index of the chromatographic peaks in the
#'   `chromPeaks` matrix associated with the feature and the row in the
#'   `XChromatograms` object in which the feature was defined. Similar to the
#'   `chromPeaks` method it is possible to filter the returned feature matrix
#'   with the `mz`, `rt` and `ppm` parameters.
#'
#' - `featureValues`: for `XChromatograms` objects only. Extract the abundance
#'   estimates for the individuals features. Note that by default (with
#'   parameter `value = "index"` a `matrix` of indices of the peaks in the
#'   `chromPeaks` matrix associated to the feature is returned. To extract the
#'   integrated peak area use `value = "into"`. The function returns a `matrix`
#'   with one row per feature (in `featureDefinitions`) and each column being
#'   a sample (i.e. column of `object`). For features without a peak associated
#'   in a certain sample `NA` is returned. This can be changed with the
#'   `missing` argument of the function.
#'
#' - `filterChromPeaks`: *filters* chromatographic peaks in `object` depending
#'   on parameter `method` and method-specific parameters passed as additional
#'   arguments with `...`. Available methods are:
#'   - `method = "keepTop"`: keep top `n` (default `n = 1L`) peaks in each
#'     chromatogram ordered by column `order` (defaults to `order = "maxo"`).
#'     Parameter `decreasing` (default `decreasing = TRUE`) can be used to
#'     order peaks in descending (`decreasing = TRUE`) or ascending (
#'     `decreasing = FALSE`) order to keep the top `n` peaks with largest or
#'     smallest values, respectively.
#'
#' - `processHistory`: returns a `list` of [ProcessHistory] objects representing
#'   the individual performed processing steps. Optional parameters `type` and
#'   `fileIndex` allow to further specify which processing steps to return.
#'
#' @section Manipulating data:
#'
#' - `transformIntensity`: transforms the intensity values of the chromatograms
#'   with provided function `FUN`. See [transformIntensity()] in the `MSnbase`
#'   package for details. For `XChromatogram` and `XChromatograms` in addition
#'   to the intensity values also columns `"into"` and `"maxo"` in the object's
#'   `chromPeaks` matrix are transformed by the same function.
#'
#' @param i For `[`: `integer` with the row indices to subset the
#'     `XChromatograms` object.
#'
#' @param j For `[`: `integer` with the column indices to subset the
#'     `XChromatograms` object.
#'
#' @param drop For `[`: `logical(1)` whether the dimensionality should be
#'     dropped (if possible). Defaults to `drop = TRUE`, thus, if length of `i`
#'     and `j` is 1 a `XChromatogram` is returned. Note that `drop` is ignored
#'     if length of `i` or `j` is larger than 1, thus a `XChromatograms` is
#'     returned.
#'
#' @param FUN For `transformIntensity`: a function to transform the intensity
#'     values of `object`.
#'
#' @param method For `featureValues`: `character(1)` specifying the method to
#'     resolve multi-peak mappings within the sample sample, i.e. to select
#'     the *representative* peak for a feature for which more than one peak
#'     was assigned in one sample. Options are `"medret"` (default): select the
#'     peak closest to the median retention time of the feature, `"maxint"`:
#'     select the peak with the largest signal and `"sum"`: sum the values
#'     of all peaks (only if `value` is `"into"` or `"maxo"`).
#'     For `filterChromPeaks`: `character(1)` defining the method that should
#'     be used to filter chromatographic peaks. See help on `filterChromPeaks`
#'     below for details.
#'
#' @param missing For `featureValues`: how missing values should be reported.
#'     Allowed values are `NA` (default), a `numeric(1)` to replace `NA`s with
#'     that value or `missing = "rowmin_half"` to replace `NA`s with half
#'     of the row's minimal (non-missing) value.
#'
#' @param rt For `chromPeaks` and `featureDefinitions`: `numeric(2)` defining
#'     the retention time range for which chromatographic peaks or features
#'     should be returned.
#'     For `filterRt`: `numeric(2)` defining the retention time range to
#'     reduce `object` to.
#'
#' @param ppm For `chromPeaks` and `featureDefinitions`: `numeric(1)` defining
#'     a ppm to expand the provided m/z range.
#'
#' @param type For `chromPeaks` and `featureDefinitions`: `character(1)`
#'     defining which peaks or features to return if `rt` or `mz` is provided:
#'     `"any"` (default) return all peaks that are even
#'     partially overlapping with `rt`, `"within"` return peaks that are
#'     completely within `rt` and `"apex_within"` return peaks which apex
#'     is within `rt`.
#'
#'     For `plot`: what type of plot should be used for the
#'     chromatogram (such as `"l"` for lines, `"p"` for points etc), see help
#'     of [plot()] in the `graphics` package for more details.
#'     For `processHistory`: restrict returned processing steps to specific
#'     types. Use [processHistoryTypes()] to list all supported values.
#'
#' @param value For `chromPeaks<-`: a numeric `matrix` with required columns
#'     `"rt"`, `"rtmin"`, `"rtmax"`, `"into"` and `"maxo"`.
#'
#'     For `featureValues`: `character(1)` specifying the name of the column in
#'     `chromPeaks(object)` that should be returned or `"index"` (default) to
#'     return the index of the peak associated with the feature in each sample.
#'     To return the integrated peak area instead of the index use
#'     `value = "into"`.
#'
#' @md
#'
#' @seealso
#'
#' [findChromPeaks-centWave][findChromPeaks-Chromatogram-CentWaveParam] for peak
#' detection on [MChromatograms()] objects.
#'
#' @examples
#'
#' ## Extract the chromatographic peaks
#' chromPeaks(xchr)
setMethod("chromPeaks", "XChromatogram", function(object, rt = numeric(),
                                                  mz = numeric(), ppm = 0,
                                                  type = c("any", "within",
                                                           "apex_within"),
                                                  msLevel) {
    type <- match.arg(type)
    pks <- object@chromPeaks
    if (!missing(msLevel))
        pks <- pks[chromPeakData(object)$ms_level %in% msLevel, , drop = FALSE]
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
    object@chromPeakData <- DataFrame(ms_level = rep(1L, nrow(value)),
                                      is_filled = rep(FALSE, nrow(value)),
                                      row.names = rownames(value))
    validObject(object)
    object
})

#' @rdname XChromatogram
#'
#' @section Plotting and visualizing:
#'
#' - `plot` draws the chromatogram and highlights in addition any
#'   chromatographic peaks present in the `XChromatogram` or `XChromatograms`
#'   (unless `peakType = "none"` was specified). To draw peaks in different
#'   colors a vector of color definitions with length equal to
#'   `nrow(chromPeaks(x))` has to be submitted  with `peakCol` and/or `peakBg`
#'   defining one color for each peak (in the order as peaks are in
#'   `chromPeaks(x))`. For base peak chromatograms or total ion chromatograms
#'   it might be better to set `peakType = "none"` to avoid generating busy
#'   plots.
#'
#' - `plotChromPeakDensity`: visualize *peak density*-based correspondence
#'   analysis results. See section *Correspondence analysis* for more details.
#'
#' @note
#'
#' Highlighting the peak area(s) in an `XChromatogram` or `XChromatograms`
#' object (`plot` with `peakType = "polygon"`) draws a polygon representing
#' the displayed chromatogram from the peak's minimal retention time to the
#' maximal retention time. If the `XChromatograms` was extracted from an
#' [XCMSnExp()] object with the [chromatogram()] function this might not
#' represent the actual identified peak area if the m/z range that was
#' used to extract the chromatogram was larger than the peak's m/z.
#'
#' @param x For `plot`: an `XChromatogram` or `XChromatograms` object.
#'
#' @param col For `plot`: the color to be used to draw the chromatogram.
#'
#' @param lty For `plot` and `plotChromPeakDensity`: the line type.
#'
#' @param xlab For `plot` and `plotChromPeakDensity`: the x axis label.
#'
#' @param ylab For `plot`: the y axis label.
#'
#' @param main For `plot` and `plotChromPeakDensity`: an optional title for
#'     the plot.
#'
#' @param peakType For `plot` and `plotChromPeakDensity`:
#'     `character(1)` defining how (and if) identified chromatographic peak
#'     within the chromatogram should be plotted. Options
#'     are `"polygon"` (default): draw the peak borders with the `peakCol` color
#'     and fill the peak area with the `peakBg` color, `"point"`: indicate the
#'     peak's apex with a point, `"rectangle"`: draw a rectangle around the
#'     identified peak and `"none"`: don't draw peaks.
#'
#' @param peakCol For `plot` and `plotChromPeakDensity`: the foreground color
#'     for the peaks. For `peakType = "polygon"` and `peakType = "rectangle"`
#'     this is the color for the border. Use `NA` to not use a foreground
#'     color. This should either be a single color or a vector of colors with
#'     the same length than `chromPeaks(x)` has rows.
#'
#' @param peakBg For `plot` and `plotChromPeakDensity`: the background color
#'     for the peaks. For `peakType = "polygon"` and `peakType = "rectangle"`
#'     the peak are or rectangle will be filled with this color. Use `NA` to
#'     skip. This should be either a single color or a vector of colors with
#'     the same length than `chromPeaks(x)` has rows.
#'
#' @param peakPch For `plot` and `plotChromPeakDensity`: the point character
#'     to be used for `peakType = "point"`. See [plot()] in the `graphics`
#'     package for more details.
#'
#' @param param For `groupChromPeaks` and `plotChromPeakDensity`: a
#'     [PeakDensityParam()] object with the settings for the *peak density*
#'     correspondence analysis algorithm.
#'
#' @param simulate For `plotChromPeakDensity`: `logical(1)` whether a
#'     correspondence analysis should be *simulated* based on the available
#'     data and the provided [PeakDensityParam()] `param` argument. See
#'     section *Correspondence analysis* for details.
#'
#' @param ... For `filterChromPeaks`: additional parameters defining how to
#'     filter chromatographic peaks. See function description below for details.
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

#' @rdname XChromatogram
#'
#' @section Filtering and subsetting:
#'
#' - `[` allows to subset a `XChromatograms` object by row (`i`) and column
#'    (`j`), with `i` and `j` being of type `integer`. The `featureDefinitions`
#'   will also be subsetted accordingly and the `peakidx` column updated.
#'
#' - `filterMz` filters the chromatographic peaks within an `XChromatogram` or
#'   `XChromatograms`, if a column `"mz"` is present in the `chromPeaks` matrix.
#'   This would be the case if the `XChromatogram` was extracted from an
#'   [XCMSnExp()] object with the [chromatogram()] function. All
#'   chromatographic peaks with their m/z within the m/z range defined by `mz`
#'   will be retained. Also feature definitions (if present) will be subset
#'   accordingly. The function returns a filtered `XChromatogram` or
#'   `XChromatograms` object.
#'
#' - `filterRt` filters chromatogram(s) by the provided retention time range.
#'   All eventually present chromatographic peaks with their apex within the
#'   retention time range specified with `rt` will be retained. Also feature
#'   definitions, if present, will be filtered accordingly. The function
#'   returns a filtered `XChromatogram` or `XChromatograms` object.
#'
#' @md
setMethod("filterMz", "XChromatogram", function(object, mz, ...) {
    if (missing(mz) || length(mz) == 0)
        return(object)
    pks <- chromPeaks(object)
    if (nrow(pks) && any(colnames(pks) == "mz")) {
        mz <- range(mz)
        keep <- which(pks[, "mz"] >= mz[1] & pks[, "mz"] <= mz[2])
        object@chromPeaks <- pks[keep, , drop = FALSE]
        object@chromPeakData <- extractROWS(object@chromPeakData, keep)
        validObject(object)
    }
    object
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
        keep <- which(pks[, "rt"] >= rt[1] & pks[, "rt"] <= rt[2])
        object@chromPeaks <- pks[keep, , drop = FALSE]
        object@chromPeakData <- extractROWS(object@chromPeakData, keep)
        validObject(object)
    }
    object
})

#' @rdname XChromatogram
#'
#' @md
setMethod("hasChromPeaks", "XChromatogram", function(object) {
    as.logical(nrow(object@chromPeaks))
})

#' @rdname XChromatogram
#'
#' @md
setMethod("dropFilledChromPeaks", "XChromatogram", function(object) {
    if (!.hasFilledPeaks(object))
        return(object)
    not_fld <- which(!object@chromPeakData$is_filled)
    object@chromPeaks <- object@chromPeaks[not_fld, , drop = FALSE]
    object@chromPeakData <- extractROWS(object@chromPeakData, not_fld)
    validObject(object)
    object
})

#' @rdname XChromatogram
setMethod("chromPeakData", "XChromatogram", function(object) {
    object@chromPeakData
})
#' @rdname XChromatogram
setReplaceMethod("chromPeakData", "XChromatogram", function(object, value) {
    object@chromPeakData <- value
    validObject(object)
    object
})

#' @rdname XChromatogram
setMethod("refineChromPeaks", c(object = "XChromatogram",
                                param = "MergeNeighboringPeaksParam"),
          function(object, param = MergeNeighboringPeaksParam()) {
              object <- .xchrom_merge_neighboring_peaks(
                  object, minProp = param@minProp, diffRt = 2 * param@expandRt)
              validObject(object)
              object
          })

#' @rdname removeIntensity-Chromatogram
setMethod("removeIntensity", "XChromatogram",
          function(object, which = c("below_threshold", "outside_chromPeak"),
                   threshold = 0) {
              which <- match.arg(which)
              if (which == "outside_chromPeak") {
                  cps <- chromPeaks(object)
                  if (nrow(cps)) {
                      keep <- rep(FALSE, length(object@rtime))
                      for (i in seq_len(nrow(cps)))
                          keep[which(object@rtime >= cps[i, "rtmin"] &
                                     object@rtime <= cps[i, "rtmax"])] <- TRUE
                      object@intensity[!keep] <- NA_real_
                  } else
                      warning("No chromatographic peaks present. ",
                              "Returning data as is")
                  return(object)
              } else callNextMethod(object, which = which, threshold = threshold)
          })

#' @rdname XChromatogram
setMethod("filterChromPeaks", "XChromatogram",
          function(object, method = c("keepTop"), ...) {
              method <- match.arg(method)
              switch(method,
                     keepTop = .filter_chrom_peaks_keep_top(object, ...))
          })

#' @rdname XChromatogram
setMethod("transformIntensity", "XChromatogram", function(object,
                                                          FUN = identity) {
    object <- callNextMethod()
    object@chromPeaks[, "into"] <- FUN(object@chromPeaks[, "into"])
    object@chromPeaks[, "maxo"] <- FUN(object@chromPeaks[, "maxo"])
    validObject(object)
    object
})
