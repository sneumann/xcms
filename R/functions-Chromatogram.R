#' @include DataClasses.R
.SUPPORTED_AGG_FUN_CHROM <- c("sum", "max", "min", "mean")
names(.SUPPORTED_AGG_FUN_CHROM) <-
    c("Total ion chromatogram (TIC).", "Base peak chromatogram (BPC).",
      "Intensity representing the minimum intensity across the mz range.",
      "Intensity representing the mean intensity across the mz range.")

#' @title Validation function for Chromatogram objects
#'
#' @description This function can be used instead of the \code{validObject} to
#'     check if the chromatogram is valid, without having to call the validity
#'     method on all super classes.
#'
#' @param object A \code{Chromatogram} object.
#' 
#' @return \code{TRUE} if the \code{object} is valid and the error messages
#'     otherwise (i.e. a \code{character}).
#' 
#' @author Johannes Rainer
#' 
#' @noRd
validChromatogram <- function(object) {
    msg <- character()
    if (length(object@rtime) != length(object@intensity))
        msg <- c(msg, "Length of 'rt' and 'intensity' have to match!")
    if (is.unsorted(object@mz))
        msg <- c(msg, "'mz' has to be increasingly ordered!")
    if (is.unsorted(object@rtime))
        msg <- c(msg, paste0("'rtime' has to be increasingly ordered!"))
    if (length(object@mz) > 0 & length(object@mz) != 2)
        msg <- c(msg, paste0("'mz' is supposed to contain the ",
                             "minimum and maximum mz values for the ",
                             "chromatogram."))
    if (length(object@filterMz) > 0 & length(object@filterMz) != 2)
        msg <- c(msg, paste0("'filterMz' is supposed to contain the ",
                             "minimum and maximum mz values of the filter",
                             " used to create the chromatogram."))
    if (length(object@precursorMz) > 0 & length(object@precursorMz) != 2)
        msg <- c(msg, paste0("'precursorMz' is supposed to be a numeric of",
                             " length 2."))
    if (length(object@productMz) > 0 & length(object@productMz) != 2)
        msg <- c(msg, paste0("'productMz' is supposed to be a numeric of",
                             " length 2."))
    if (length(object@fromFile) > 1 | any(object@fromFile < 0))
        msg <- c(msg, paste0("'fromFile' is supposed to be a single ",
                             "positive integer!"))
    if (length(object@aggregationFun) > 1)
        msg <- c(msg, "Length of 'aggregationFun' has to be 1!")
    if (length(object@aggregationFun)) {
        if (!object@aggregationFun %in% .SUPPORTED_AGG_FUN_CHROM)
            msg <- c(msg, paste0("Invalid value for 'aggregationFun'! only ",
                                 paste0("'", .SUPPORTED_AGG_FUN_CHROM,"'",
                                        collapse = ","), " are allowed!"))
    }
    if (length(msg) == 0) TRUE
    else msg
}

#' @description \code{Chromatogram}: create an instance of the
#'     \code{Chromatogram} class.
#'
#' @param rtime \code{numeric} with the retention times (length has to be equal
#'     to the length of \code{intensity}).
#'
#' @param intensity \code{numeric} with the intensity values (length has to be
#'     equal to the length of \code{rtime}).
#'
#' @param mz \code{numeric(2)} representing the mz value range (min, max)
#'     on which the chromatogram was created. This is supposed to contain the
#'     \emph{real} range of mz values in contrast to the \code{filterMz} below.
#'     If not applicable use \code{mzrange = c(0, 0)}.
#'
#' @param filterMz \code{numeric(2)} representing the mz value range (min,
#'     max) that was used to filter the original object on mz dimension. If not
#'     applicable use \code{filterMz = c(0, 0)}.
#'
#' @param precursorMz \code{numeric(2)} for SRM/MRM transitions.
#'     Represents the mz of the precursor ion. See details for more information.
#' 
#' @param productMz \code{numeric(2)} for SRM/MRM transitions.
#'     Represents the mz of the product. See details for more information.
#' 
#' @param fromFile \code{integer(1)} the index of the file within the
#'     \code{\link{OnDiskMSnExp}} or \code{\link{XCMSnExp}} from which the
#'     chromatogram was extracted.
#'
#' @param aggregationFun \code{character} string specifying the function that
#'     was used to aggregate intensity values for the same retention time across
#'     the mz range. Supported are \code{"sum"} (total ion chromatogram),
#'     \code{"max"} (base peak chromatogram), \code{"min"} and \code{"mean"}.
#' 
#' @slot .__classVersion__,rtime,intensity,mz,filterMz,precursorMz,productMz,fromFile,aggregationFun See corresponding parameter above.
#' 
#' @rdname Chromatogram-class
Chromatogram <- function(rtime = numeric(), intensity = numeric(),
                         mz = c(0, 0), filterMz = c(0, 0),
                         precursorMz = c(NA_real_, NA_real_),
                         productMz = c(NA_real_, NA_real_),
                         fromFile = integer(),
                         aggregationFun = character()) {
    ## Check if we have to re-order the data (issue #145).
    if (is.unsorted(rtime)) {
        idx <- order(rtime)
        rtime <- rtime[idx]
        intensity <- intensity[idx]
    }
    return(new("Chromatogram", rtime = rtime, intensity = intensity,
               mz = range(mz), filterMz = range(filterMz),
               precursorMz = range(precursorMz), productMz = range(productMz),
               fromFile = as.integer(fromFile), aggregationFun = aggregationFun))
}

#' @title Plot Chromatogram objects
#' 
#' @description \code{plotChromatogram} creates a chromatogram plot for a
#'     single \code{Chromatogram} object or a \code{list} of
#'     \code{\link{Chromatogram}} objects (one line for each
#'     \code{\link{Chromatogram}}/sample).
#'
#' @details The \code{plotChromatogram} function allows to efficiently plot
#'     the chromatograms of several samples into a single plot.
#' 
#' @param x For \code{plotChromatogram}: \code{list} of
#'     \code{\link{Chromatogram}} objects. Such as extracted from an
#'     \code{\link{XCMSnExp}} object by the \code{\link{extractChromatograms}}
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
#' @seealso \code{\link{extractChromatograms}} for how to extract a list of
#'     \code{\link{Chromatogram}} objects from an \code{\link{XCMSnExp}}
#'     objects.
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
#' chrs <- extractChromatograms(od, rt = rtr, mz = mzr)
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
    rts <- do.call(cbind, lapply(x, rtime))
    ints <- do.call(cbind, lapply(x, intensity))
    ## ## Fix problem if we've got only NAs:
    ## dotList <- list(...)
    ## if (is.null(ylim)) {
    ##     if (all(is.na(ints)))
    ##         ylim <- c(0, 1)
    ##     else
    ##         ylim <- range(ints, na.rm = TRUE, finite = TRUE)
    ## }
    ## Skip columns that have only NAs?
    keepCol <- which(apply(ints, MARGIN = 2, function(z) any(!is.na(z))))
    if (length(keepCol)) {
        matplot(x = rts[, keepCol, drop = FALSE],
                y = ints[, keepCol, drop = FALSE], type = type, lty = lty,
                col = col[keepCol], xlab = xlab, ylab = ylab, main = main,
                ...)
    } else
        plot(x = 3, y = 3, pch = NA, xlab = xlab, ylab = ylab, main = main,
             xlim = range(rts, na.rm = TRUE), ylim = c(0, 1))
}

#' @description The \code{highlightChromPeaks} function adds chromatographic
#'     peak definitions to an existing plot, such as one created by the
#'     \code{plotChromatograms} function.
#'
#' @param mz \code{numeric(2)} with the mz range from which the peaks should
#'     be extracted and plotted.
#'
#' @param border colors to be used to color the border of the rectangles. Has to
#'     be equal to the number of samples in \code{x}.
#'
#' @param lwd \code{numeric(1)} defining the width of the line/border.
#'
#' @rdname plotChromatogram
highlightChromPeaks <- function(x, rt, mz,
                                border = rep("00000040", length(fileNames(x))),
                                lwd = 1, col = NA, type = c("rect", "point"),
                                ...) {
    type <- match.arg(type)
    if (!is(x, "XCMSnExp"))
        stop("'x' has to be a XCMSnExp object")
    if (!hasChromPeaks(x))
        stop("'x' does not contain any detected peaks")
    pks <- chromPeaks(x, rt = rt, mz = mz, ppm = 0)
    if (length(col) != length(fileNames(x)))
        col <- rep(col[1], length(fileNames(x)))
    if (length(border) != length(fileNames(x)))
        border <- rep(border[1], length(fileNames(x)))
    if (length(pks)) {
        if (type == "rect")
            rect(xleft = pks[, "rtmin"], xright = pks[, "rtmax"],
                 ybottom = rep(0, nrow(pks)), ytop = pks[, "maxo"],
                 border = border[pks[, "sample"]], lwd = lwd,
                 col = col[pks[, "sample"]])
        if (type == "point") {
            if (any(is.na(col)))
                col <- border
            ## Draw a star at the position defined by the "rt" column
            points(x = pks[, "rt"], y = pks[, "maxo"],
                   col = col[pks[, "sample"]], ...)
        }
    }
}

