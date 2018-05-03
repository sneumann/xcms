## Unsorted utility functions.
#' @include DataClasses.R

############################################################
## valueCount2ScanIndex
##
#' @title Create index vector for internal C calls
#' 
#' @description Simple helper function that converts the number of values
#'     per scan/spectrum to an integer vector that can be passed to the base
#'     xcms functions/downstream C functions.
#'
#' @param valCount Numeric vector representing the number of values per
#'     spectrum.
#' 
#' @return An integer vector with the index (0-based) in the mz or intensity
#'     vectors indicating the start of a spectrum.
#' 
#' @author Johannes Rainer
#'
#' @noRd
valueCount2ScanIndex <- function(valCount){
    ## Convert into 0 based.
    valCount <- cumsum(valCount)
    return(as.integer(c(0, valCount[-length(valCount)])))
}

############################################################
## useOriginalCode
##
## Simple function allowing the user to enable using the orignal
## code instead of the new implementations.
## This sets options.
##
#' @title Enable usage of old xcms code
#'
#' @description This function allows to enable the usage of old, partially
#'     deprecated code from xcms by setting a corresponding global option. See
#'     details for functions affected.
#'
#' @note Usage of old code is strongly dicouraged. This function is thought
#'     to be used mainly in the transition phase from xcms to xcms version 3.
#' 
#' @details The functions/methods that will be affected by this are:
#'     \itemize{
#'     \item \code{\link{do_findChromPeaks_matchedFilter}}
#'     }
#' 
#' @param x logical(1) to specify whether or not original
#'     old code should be used in corresponding functions. If not provided the
#'     function simply returns the value of the global option.
#' 
#' @return logical(1) indicating whether old code is being used.
#' 
#' @author Johannes Rainer
useOriginalCode <- function(x) {
    if (missing(x)) {
        res <- options()$BioC$xcms$useOriginalCode
        if (is.null(res))
            return(FALSE)
        return(res)
    }
    if (!is.logical(x))
        stop("'x' has to be logical.")
    b_opts <- getOption("BioC")
    x_opts <- b_opts$xcms
    x_opts$useOriginalCode <- x[1]
    b_opts$xcms <-  x_opts
    options(BioC = b_opts)
    return(options()$BioC$xcms$useOriginalCode)
}

## .getOriginalFunction <- function(x) {
##     if (!any(names(.ORIGINAL_FUNCTIONS)))
## }
## .ORIGINAL_FUNCTIONS <- c(
##     matchedFilter = ".matchedFilter_orig"
## )

#' @title Copy the content from an environment to another one
#' 
#' @description This function copies the content of an environment into another
#'     one.
#' 
#' @param env environment from which to copy.
#' 
#' @param inheritLocks logical(1) whether the locking status should be copied
#'     too.
#' 
#' @return an env.
#' 
#' @noRd
.copy_env <- function(env, inheritLocks = FALSE) {
    new_e <- new.env(parent = emptyenv())
    eNames <- ls(env, all.names = TRUE)
    if (length(eNames) > 0) {
        for (eN in eNames) {
            new_e[[eN]] <- env[[eN]]
        }
    }
    if (inheritLocks) {
        if (environmentIsLocked(env))
            lockEnvironment(new_e)
    }
    return(new_e)
}

## #' Simulates the \code{findRange} function.
## #' @noRd
## findRangeR <- function(x, values) {
##     start <- min(which(x >= values[1]))
##     end <- max(which(x <= values[2]))
##     return(c(start, end))
## }

############################################################
## .createProfileMatrix
#' @title Create the profile matrix
#'
#' @description This function creates a \emph{profile} matrix, i.e. a rt times
#'     m/z matrix of aggregated intensity values with values aggregated within
#'     bins along the m/z dimension.
#'
#' @details This is somewhat the successor function for the deprecated
#'     \code{profBin} methods (\code{profBinM}, \code{profBinLinM},
#'     \code{profBinLinBaseM} and \code{profIntLin}).
#'
#' @param mz Numeric representing the m/z values across all scans/spectra.
#'
#' @param int Numeric representing the intensity values across all
#'     scans/spectra.
#'
#' @param valsPerSpect Numeric representing the number of measurements for each
#'     scan/spectrum.
#'
#' @param method A character string specifying the profile matrix generation
#'     method. Allowed are \code{"bin"}, \code{"binlin"},
#'     \code{"binlinbase"} and \code{"intlin"}.
#'
#' @param step Numeric specifying the size of the m/z bins.
#'
#' @param baselevel Numeric specifying the base value.
#'
#' @param basespace Numeric.
#'
#' @param mzrange. numeric(2) optionally specifying the mz value range
#'     for binning. This is to adopt the old profStepPad<- method used for
#'     obiwarp retention time correction that did the binning from
#'     whole-number limits.
#'
#' @param returnBreaks logical(1): hack to return the breaks of the bins.
#'     Setting this to TRUE causes the function to return a \code{list} with
#'     elements \code{"$profMat"} and \code{"breaks"}.
#'
#' @param baseValue numeric(1) defining the value to be returned if no signal
#'     was found in the corresponding bin. Defaults to 0 for backward
#'     compatibility.
#'
#' @noRd
.createProfileMatrix <- function(mz, int, valsPerSpect,
                                 method, step = 0.1, baselevel = NULL,
                                 basespace = NULL,
                                 mzrange. = NULL,
                                 returnBreaks = FALSE,
                                 baseValue = 0) {
    profMeths <- c("bin", "binlin", "binlinbase", "intlin")
    names(profMeths) <- c("none", "lin", "linbase", "intlin")
    method <- match.arg(method, profMeths)
    impute <- names(profMeths)[profMeths == method]
    brks <- NULL
    
    if (length(mzrange.) != 2) {
        mrange <- range(mz, na.rm = TRUE)
        mzrange. <- c(floor(mrange[1] / step) * step,
                      ceiling(mrange[2] / step) * step)
    }
    mass <- seq(mzrange.[1], mzrange.[2], by = step)
    mlength <- length(mass)
    ## Calculate the "real" bin size; old xcms code oddity that that's different
    ## from step.
    bin_size <- (mass[mlength] - mass[1]) / (mlength - 1)
    ## Define the breaks.
    toIdx <- cumsum(valsPerSpect)
    fromIdx <- c(1L, toIdx[-length(toIdx)] + 1L)
    shiftBy <- TRUE
    binFromX <- min(mass)
    binToX <- max(mass)
    brks <- breaks_on_nBins(fromX = binFromX, toX = binToX,
                            nBins = mlength, shiftByHalfBinSize = TRUE)
    ## for profIntLinM we have to use the old code.
    if (impute == "intlin") {
        profFun <- "profIntLinM"
        profp <- list()
        scanindex <- valueCount2ScanIndex(valsPerSpect)
        buf <- do.call(profFun, args = list(mz, int,
                                            scanindex, mlength,
                                            mass[1], mass[mlength],
                                            TRUE))
    } else {
        ## Binning the data.
        binRes <- binYonX(mz, int,
                          breaks = brks,
                          fromIdx = fromIdx,
                          toIdx = toIdx,
                          baseValue = ifelse(impute == "none", yes = baseValue,
                                             no = NA),
                          sortedX = TRUE,
                          returnIndex = FALSE,
                          returnX = FALSE
                          )
        if (length(toIdx) == 1)
            binRes <- list(binRes)
        ## Missing value imputation.
        if (impute == "linbase") {
            ## need arguments distance and baseValue.
            if (length(basespace) > 0) {
                if (!is.numeric(basespace))
                    stop("'basespace' has to be numeric!")
                distance <- floor(basespace[1] / bin_size)
            } else {
                distance <- floor(0.075 / bin_size)
            }
            if (length(baselevel) > 0) {
                if (!is.numeric(baselevel))
                    stop("'baselevel' has to be numeric!")
                baseValue <- baselevel
            } else {
                baseValue <- min(int, na.rm = TRUE) / 2
            }
        } else {
            distance <- 0
            baseValue <- 0
        }
        if (method == "none") {
            ## binVals <- lapply(binRes, function(z) z$y)
            binVals <- binRes
        } else {
            binVals <- lapply(binRes, function(z) {
                imputeLinInterpol(z$y, method = impute, distance = distance,
                                  noInterpolAtEnds = TRUE,
                                  baseValue = baseValue)
            })
        }
        buf <- base::do.call(cbind, binVals)
    }
    if (returnBreaks)
        buf <- list(profMat = buf, breaks = brks)
    buf
}

#' @description This function creates arbitrary IDs for features.
#' 
#' @param x integer(1) with the number of IDs that should be generated.
#'
#' @noRd
.featureIDs <- function(x) {
    sprintf(paste0("FT%0", ceiling(log10(x + 1L)), "d"), 1:x)
}

## #' @description Expands stretches of TRUE values in \code{x} by one on both
## #'     sides.
## #'
## #' @note The return value for a \code{NA} is always \code{FALSE}.
## #' 
## #' @param x \code{logical} vector.
## #'
## #' @author Johannes Rainer
## #' 
## #' @noRd
## .grow_trues <- function(x) {
##     previous <- NA
##     x_new <- rep_len(FALSE, length(x))
##     for (i in 1:length(x)) {
##         if (is.na(x[i])) {
##             previous <- NA
##             next
##         }
##         ## If current element is TRUE
##         if (x[i]) {
##             x_new[i] <- TRUE
##             ## if last element was FALSE, set last element to TRUE
##             if (!is.na(previous) && !previous)
##                 x_new[i - 1] <- TRUE
##         } else {
##             ## if previous element was TRUE, set current to TRUE.
##             if (!is.na(previous) && previous)
##                 x_new[i] <- TRUE
##         }
##         previous <- x[i]
##     }
##     x_new
## }

#' @title Weighted mean around maximum
#'
#' @description Calculate a weighted mean of the values around the value with
#'     the largest weight. \code{x} could e.g. be mz values and \code{w} the
#'     corresponding intensity values.
#'
#' @param x \code{numeric} vector from which the weighted mean should be
#'     calculated.
#'
#' @param w \code{numeric} of same length than \code{x} with the weights.
#'
#' @param i \code{integer(1)} defining the number of data points left and right
#'     of the index with the largest weight that should be considered for the
#'     weighted mean calculation.
#'
#' @return The weighted mean value.
#'
#' @author Johannes Rainer
#'
#' @noRd
#'
#' @examples
#'
#' mz <- c(124.0796, 124.0812, 124.0828, 124.0843, 124.0859, 124.0875,
#'     124.0890, 124.0906, 124.0922, 124.0938, 124.0953, 124.0969)
#' ints <- c(10193.8, 28438.0, 56987.6, 85107.6, 102531.6, 104262.6,
#'     89528.8, 61741.2, 33485.8, 14146.6, 5192.2, 1630.2)
#'
#' plot(mz, ints)
#'
#' ## What would be found by the max:
#' abline(v = mz[which.max(ints)], col = "grey")
#' ## What does the weighted mean around apex return:
#' abline(v = weightedMeanAroundApex(mz, ints, i = 2), col = "blue")
weightedMeanAroundApex <- function(x, w = rep(1, length(x)), i = 1) {
    max_idx <- which.max(w)
    seq_idx <- max(1, max_idx - i):min(length(x), max_idx + i)
    weighted.mean(x[seq_idx], w[seq_idx])
}



#' @title DEPRECATED: Create a plot that combines a XIC and a mz/rt 2D plot for one sample
#'
#' @description
#'
#' **UPDATE**: please use `plot(x, type = "XIC")` from the `MSnbase` package
#' instead. See examples below.
#' 
#' The `plotMsData` creates a plot that combines an (base peak )
#' extracted ion chromatogram on top (rt against intensity) and a plot of
#' rt against m/z values at the bottom.
#' 
#' @param x `data.frame` such as returned by the [extractMsData()] function.
#'     Only a single `data.frame` is supported.
#'
#' @param main `character(1)` specifying the title.
#'
#' @param cex `numeric(1)` defining the size of points. Passed directly to the
#'     `plot` function.
#' 
#' @param mfrow `numeric(2)` defining the plot layout. This will be passed
#'     directly to `par(mfrow = mfrow)`. See `par` for more information. Setting
#'     `mfrow = NULL` avoids calling `par(mfrow = mfrow)` hence allowing to
#'     pre-define the plot layout.
#'
#' @param grid.color a color definition for the grid line (or `NA` to skip
#'     creating them).
#'
#' @param colramp a *color ramp palette* to be used to color the data points
#'     based on their intensity. See argument `col.regions` in
#'     [lattice::level.colors] documentation.
#' 
#' @author Johannes Rainer
#' 
#' @md
#'
#' @examples
#'
#' ## Read two files from the faahKO package
#' library(faahKO)
#' library(magrittr)
#' cdfs <- dir(system.file("cdf", package = "faahKO"), full.names = TRUE,
#'     recursive = TRUE)[1:2]
#' raw_data <- readMSData(cdfs, mode = "onDisk")
#'
#' ## Subset the object to a rt and mz range and plot the data.
#' raw_data %>%
#'     filterRt(rt = c(2700, 2900)) %>%
#'     filterMz(mz = c(334.9, 335.1)) %>%
#'     plot(type = "XIC")
plotMsData <- function(x, main = "", cex = 1, mfrow = c(2, 1),
                       grid.color = "lightgrey",
                       colramp = colorRampPalette(
                           rev(brewer.pal(9, "YlGnBu")))) {
    .Deprecated(msg = paste0("'plotMsData' is deprecated. Please use ",
                             "'plot(x, type = \"XIC\") instead."))
    if (length(mfrow) == 2)
        par(mfrow = mfrow)
    par(mar = c(0, 4, 2, 1))
    x_split <- split(x$i, f = x$rt)
    ints <- unlist(lapply(x_split, function(z) max(z)))
    brks <- do.breaks(range(x$i), nint = 256)
    cols <- level.colors(ints, at = brks, col.regions = colramp)
    plot(as.numeric(names(ints)), ints, main = main, xlab = "", xaxt = "n",
         ylab = "", las = 2, pch = 21, bg = cols, col = "grey", cex = cex)
    mtext(side = 4, line = 0, "intensity", cex = par("cex.lab"))
    grid(col = grid.color)
    par(mar = c(3.5, 4, 0, 1))
    cols <- level.colors(x$i, at = brks, col.regions = colramp)
    plot(x$rt, x$mz, main = "", pch = 21, bg = cols, col = "grey",
         xlab = "", ylab = "", yaxt = "n", cex = cex)
    axis(side = 2, las = 2)
    grid(col = grid.color)
    mtext(side = 1, line = 2.5, "retention time", cex = par("cex.lab"))
    mtext(side = 4, line = 0, "mz", cex = par("cex.lab"))
}

#' @title Calculate relative log abundances
#' 
#' `rla` calculates the relative log abundances (RLA, see reference) on a
#' `numeric` vector.
#'
#' @details The RLA is defines as the (log) abundance of an analyte relative
#'     to the median across all abundances of the same group.
#' 
#' @param x `numeric` (for `rla`) or `matrix` (for `rowRla`) with the
#'     abundances (in natural scale) on which the RLA should be calculated.
#'
#' @param group `factor`, `numeric` or `character` with the same length
#'     than `x` that groups values in `x`. If omitted all values are considered
#'     to be from the same group.
#'
#' @param log.transform `logical(1)` whether `x` should be log2 transformed.
#'     Set to `log.transform = FALSE` if `x` is already in log scale.
#' 
#' @return `numeric` of the same length than `x` (for `rla`) or `matrix` with
#'     the same dimensions than `x` (for `rowRla`).
#'
#' @rdname rla
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @references
#'
#' De Livera AM, Dias DA, De Souza D, Rupasinghe T, Pyke J, Tull D, Roessner U,
#' McConville M, Speed TP. Normalizing and integrating metabolomics data.
#' *Anal Chem* 2012 Dec 18;84(24):10768-76.
#' 
#' @examples
#'
#' x <- c(3, 4, 5, 1, 2, 3, 7, 8, 9)
#'
#' grp <- c(1, 1, 1, 2, 2, 2, 3, 3, 3)
#'
#' rla(x, grp)
rla <- function(x, group, log.transform = TRUE) {
    if (missing(group))
        group <- rep_len(1, length(x))
    if (length(x) != length(group))
        stop("length of 'x' has to match length of 'group'")
    if (!is.factor(group))
	group <- factor(group, levels = unique(group))
    ## Calculate group medians.
    if (log.transform)
        log_x <- log2(x)
    grp_meds <- unlist(lapply(split(log_x, group), median, na.rm = TRUE))
    log_x - grp_meds[group]
}

#' `rowRla` calculates row-wise RLAs.
#'
#' @rdname rla
#'
#' @md
rowRla <- function(x, group, log.transform = TRUE) {
    t(apply(x, MARGIN = 1, rla, group = group, log.transform = log.transform))
}
