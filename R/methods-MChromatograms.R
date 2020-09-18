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

#' @rdname align-Chromatogram
setMethod("align", signature = c(x = "MChromatograms", y = "Chromatogram"),
          function(x, y, method = c("matchRtime", "approx"), ...) {
              x@.Data <- matrix(lapply(x@.Data, .align_chromatogram, y = y,
                                       method = method, ...),
                                nrow = nrow(x), dimnames = dimnames(x))
              validObject(x)
              x
          })

#' @rdname correlate-Chromatogram
setMethod("correlate", signature = c(x = "MChromatograms", y = "missing"),
          function(x, y = NULL, use = "pairwise.complete.obs",
                   method = c("pearson", "kendall", "spearman"),
                   align = c("matchRtime", "approx"), full = TRUE, ...) {
              .correlate_chromatograms(x, x, use = use, method = method,
                                       align = align, full = full, ...)
          })

#' @rdname correlate-Chromatogram
setMethod("correlate", signature = c(x = "MChromatograms", y = "MChromatograms"),
          function(x, y = NULL, use = "pairwise.complete.obs",
                   method = c("pearson", "kendall", "spearman"),
                   align = c("matchRtime", "approx"), ...) {
              .correlate_chromatograms(x, y, use = use, method = method,
                                       align = align, ...)
          })

#' @title Correlate chromatograms within an `MChromatograms` with each other
#'
#' @description
#'
#' The `correlate_chromatograms_self` function correlates each chromatogram
#' in an (single column) [MChromatograms()] object with each other.
#'
#' @note
#'
#' This is an optimized version that does not calculate the lower triangular
#' part of the correlation matrix - but might thus also not be correct,
#' @param x `MChromatograms` object or a `list` of [Chromatogram] objects.
#'
#' @param ... Additional arguments to be passed down to the [correlate()]
#'     function.
#'
#' @return `matrix` with the correlation coefficients.
#'
#' @importMethodsFrom xcms correlate
#'
#' @importFrom MSnbase Chromatogram
#'
#' @author Johannes Rainer
#'
#' @noRd
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
#' correlate_chromatograms_self(list(chr1, chr2, chr3))
#'
#' chrs <- MChromatograms(list(chr1, chr2, chr3))
#' correlate_chromatograms_self(chrs)
.correlate_chromatograms_self <- function(x, ...) {
    if (inherits(x, "MChromatograms")) {
        if (ncol(x) > 1)
            stop("Currently only single column MChromatograms are supported.")
        x <- unlist(x)
    }
    nx <- length(x)
    m <- matrix(NA_real_, nrow = nx, ncol = nx)
    cb <- which(lower.tri(m, diag = TRUE), arr.ind = TRUE)
    for (i in seq_len(nrow(cb))) {
        cur <- cb[i, 2L]
        if (i == 1L || cb[i - 1L, 2L] != cur)
            py <- px <- x[[cur]]
        else
            py <- x[[cb[i, 1L]]]
        m[cb[i, 1L], cur] <- m[cur, cb[i, 1L]] <-
            correlate(px, py, ...)
    }
    m
}

.correlate_chromatograms <- function(x, y, full = TRUE, ...) {
    if (inherits(x, "MChromatograms")) {
        if (ncol(x) > 1)
            stop("Currently only single column MChromatograms are supported.")
        x <- unlist(x)
    }
    if (inherits(y, "MChromatograms")) {
        if (ncol(y) > 1)
            stop("Currently only single column MChromatograms are supported.")
        y <- unlist(y)
    }
    nx <- length(x)
    ny <- length(y)

    m <- matrix(NA_real_, nrow = nx, ncol = ny)

    for (i in seq_len(nx)) {
        for (j in seq_len(ny)) {
            if (!full && i > j)
                next
            m[i, j] <- correlate(x[[i]], y[[j]], ...)
        }
    }
    m
}

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

#' @rdname normalize-Chromatogram
setMethod("normalize", "MChromatograms",
          function(object, method = c("max", "sum")) {
              method <- match.arg(method)
              object@.Data <- matrix(lapply(c(object@.Data),
                                            FUN = .normalize_chromatogram,
                                            method = method),
                                     ncol = ncol(object),
                                     dimnames = dimnames(object@.Data))
              object
          })
