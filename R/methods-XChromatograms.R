setAs("Chromatograms", "XChromatograms", function(from) {
    res <- new("XChromatograms")
    res@.Data <- matrix(lapply(from, function(z) {
        if (is(z, "Chromatogram"))
            as(z, "XChromatogram")
        else z
    }), nrow = nrow(from), ncol = ncol(from), dimnames = dimnames(from))
    res@phenoData <- from@phenoData
    res@featureData <- from@featureData
    if (validObject(resetClass)) res
})

#' @rdname XChromatogram
setMethod("show", "XChromatograms", function(object) {
    nr <- nrow(object)
    nc <- ncol(object)
    cat(class(object), " with ",
        nr, ifelse(nr == 1, " row and ", " rows and "),
        nc, ifelse(nc == 1, " column\n", " columns\n"),
        sep = "")
    sumFun <- function(z) {
        paste0("peaks: ", nrow(z[[1]]@chromPeaks))
    }
    if (nr > 0 && nc > 0) {
        if (nr <= 4) {
            out <- apply(object, MARGIN = c(1, 2), sumFun)
            rownames(out) <- paste0("[", 1:nrow(out), ",]")
        }
        else {
            out <- rbind(
                apply(object[c(1, 2), , drop = FALSE], MARGIN = c(1, 2), sumFun),
                rep(" ... ", ncol(object)),
                apply(object[nrow(object) - c(1, 0), , drop = FALSE],
                      MARGIN = c(1, 2), sumFun)
            )
            rownames(out) <- c("[1,]", "[2,]", "...",
                               paste0("[", c(nrow(object) - c(1, 0)), ",]"))
        }
        rn <- rownames(out)
        out <- rbind(rep("<XChromatogram>", ncol(out)), out)
        rownames(out) <- c("", rn)
        print(out, quote = FALSE, right = TRUE)
    }
    cat("phenoData with", length(varLabels(object@phenoData)), "variables\n")
    cat("featureData with", length(fvarLabels(object)), "variables\n")
})

#' @rdname XChromatogram
setMethod("hasChromPeaks", "XChromatograms", function(object) {
    matrix(vapply(object, hasChromPeaks, logical(1)), ncol = ncol(object),
           dimnames = dimnames(object))
})

#' @rdname XChromatogram
#'
#' @examples
#'
#' ## Extract the chromatographic peaks
#' chromPeaks(xchrs)
setMethod("chromPeaks", "XChromatograms", function(object, rt = numeric(),
                                                   mz = numeric(), ppm = 0,
                                                   type = c("any", "within",
                                                            "apex_within")){
    type <- match.arg(type)
    res <- lapply(object, chromPeaks, rt = rt, mz = mz, ppm = ppm, type = type)
    nrs <- vapply(res, nrow, integer(1))
    row_idx <- rep(seq_len(nrow(object)), ncol(object))
    col_idx <- rep(seq_len(ncol(object)), each = nrow(object))
    res <- do.call(rbind, res)
    cbind(res, row = rep(row_idx, nrs), column = rep(col_idx, nrs))
})

#' @rdname XChromatogram
setMethod("filterMz", "XChromatograms", function(object, mz, ...) {
    if (missing(mz) || length(object) == 0)
        return(object)
    object@.Data <- matrix(lapply(object, filterMz, mz = mz, ...),
                           nrow = nrow(object), dimnames = dimnames(object))
    if (validObject(object)) object
})

#' @rdname XChromatogram
setMethod("filterRt", "XChromatograms", function(object, rt, ...) {
    if (missing(rt) || length(object) == 0)
        return(object)
    object@.Data <- matrix(lapply(object, filterRt, rt = rt, ...),
                           nrow = nrow(object), dimnames = dimnames(object))
    if (validObject(object)) object
})

setMethod("addProcessHistory", "XChromatograms", function(object, ph) {
    if (!inherits(ph, "ProcessHistory"))
        stop("Argument 'ph' has to be of type 'ProcessHistory' or a class ",
             "extending it!")
    object@.processHistory[[(length(object@.processHistory) + 1)]] <- ph
    object
})

#' @rdname XChromatogram
setMethod("plot", "XChromatograms", function(x, col = "#00000060", lty = 1,
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
    nr <- nrow(x)
    if (nr > 1)
        par(mfrow = c(round(sqrt(nr)), ceiling(sqrt(nr))))
    pks_all <- chromPeaks(x)
    pks_nr <- nrow(pks_all)
    if (length(peakCol) != pks_nr)
        peakCol <- rep(peakCol[1], pks_nr)
    if (length(peakBg) != pks_nr)
        peakBg <- rep(peakBg[1], pks_nr)
    if (length(peakPch) != pks_nr)
        peakPch <- rep(peakPch[1], pks_nr)
    for (i in seq_len(nr)) {
        x_sub <- x[i, ]
        plot(as(x_sub, "Chromatograms"), col = col, lty = lty, type = type,
             xlab = xlab, ylab = ylab, main = main, ...)
        idx <- which(pks_all[, "row"] == i)
        if (length(idx) && peakType != "none") {
            pks <- chromPeaks(x_sub)
            .add_chromatogram_peaks(x_sub, pks, col = peakCol[idx],
                                    bg = peakBg[idx], type = peakType,
                                    pch = peakPch[idx], ...)
        }
    }
})
