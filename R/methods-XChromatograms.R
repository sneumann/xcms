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
    cat("- - - xcms preprocessing - - -\n")
    if (any(hasChromPeaks(object))) {
        cat("Chromatographic peak detection:\n")
        ph <- processHistory(object, type = .PROCSTEP.PEAK.DETECTION)
        if (length(ph))
            cat(" method:", .param2string(ph[[1]]@param), "\n")
        else cat(" unknown method.\n")
    }
    if (hasFeatures(object)) {
        cat("Correspondence:\n")
        ph <- processHistory(object, type = .PROCSTEP.PEAK.GROUPING)
        if (length(ph))
            cat(" method:", .param2string(ph[[1]]@param), "\n")
        else cat(" unknown method.\n")
        cat(" ", nrow(object@featureDefinitions), " feature(s) identified.\n",
            sep = "")
        ## if (.hasFilledPeaks(object)) {
        ##     totF <- chromPeaks(object)[, "is_filled"] == 1
        ##     fp <- chromPeaks(object)[totF, , drop = FALSE]
        ##     cat("", sum(totF), "filled peaks (on average",
        ##         mean(table(fp[, "sample"])), "per sample).\n")
        ## }
    }
})

#' @rdname XChromatogram
setMethod("hasChromPeaks", "XChromatograms", function(object) {
    matrix(vapply(object, hasChromPeaks, logical(1)), ncol = ncol(object),
           dimnames = dimnames(object))
})

#' @rdname XChromatogram
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
    res <- cbind(res, row = rep(row_idx, nrs), column = rep(col_idx, nrs))
    res[order(res[, "row"]), , drop = FALSE]
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

#' @rdname XChromatogram
#'
#' @param fileIndex For `processHistory`: optional `integer` specifying the
#'     index of the files/samples for which the [ProcessHistory] objects should
#'     be returned.
#'
#' @md
setMethod("processHistory", "XChromatograms", function(object, fileIndex,
                                                       type) {
    ph <- object@.processHistory
    if (length(ph)) {
        if (!missing(fileIndex)) {
            if (!all(fileIndex %in% seq_len(ncol(object))))
                stop("'fileIndex' has to be within 1 and the number of samples!")
            gotIt <- vapply(ph, function(z) any(z@fileIndex %in% fileIndex),
                            logical(1))
            ph <- ph[gotIt]
        }
        if (!missing(type) && length(ph)) {
            gotIt <- vapply(ph, function(z) any(type == processType(z)),
                            logical(1))
            ph <- ph[gotIt]
        }
        return(ph)
    }
    list()
})

#' @rdname XChromatogram
#'
#' @md
setMethod("hasFeatures", "XChromatograms", function(object, ...) {
    nrow(object@featureDefinitions) > 0
})

#' @rdname XChromatogram
#'
#' @md
setMethod("dropFeatureDefinitions", "XChromatograms", function(object, ...) {
    if (!hasFeatures(object))
        return(object)
    object <- dropProcessHistories(object, type = .PROCSTEP.PEAK.GROUPING, 1)
    object@featureDefinitions <- DataFrame()
    if (validObject(object))
        object
})

#' @rdname XChromatogram
#'
#' @section Correspondence analysis:
#'
#' Identified chromatographic peaks in an `XChromatograms` object can be grouped
#' into *features* with the `groupChromPeaks` function. Currently, such a
#' correspondence analysis can be performed with the *peak density* method
#' (see [groupChromPeaks] for more details) specifying the algorithm settings
#' with a [PeakDensityParam()] object. A correspondence analysis is performed
#' separately for each row in the `XChromatograms` object grouping
#' chromatographic peaks across samples (columns).
#'
#' The analysis results are stored in the returned `XChromatograms` object
#' and can be accessed with the `featureDefinitions` method which returns a
#' `DataFrame` with one row for each feature. Column `"row"` specifies in
#' which row of the `XChromatograms` object the feature was identified.
#'
#' @param param For `groupChromPeaks`: a [PeakDensityParam()] object with the
#'     settings for the *peak density* correspondence analysis algorithm.
#'
#' @md
setMethod("groupChromPeaks",
          signature(object = "XChromatograms", param = "PeakDensityParam"),
          function(object, param) {
              if (!any(hasChromPeaks(object)))
                  stop("No chromatographic peak detection results in 'object'! ",
                       "Please perform first a peak detection using the ",
                       "'findChromPeaks' method.")
              if (hasFeatures(object))
                  object <- dropFeatureDefinitions(object)
              ## Check if we've got any sample groups:
              if (length(sampleGroups(param)) == 0) {
                  sampleGroups(param) <- rep(1, ncol(object))
                  message("Empty 'sampleGroups' in 'param', assuming all ",
                          "samples to be in the same group.")
              } else {
                  ## Check that the sampleGroups are OK
                  if (length(sampleGroups(param)) != ncol(object))
                      stop("The 'sampleGroups' value in the provided 'param' ",
                           "class does not match the number of available files/",
                           "samples!")
              }
              startDate <- date()
              cpks <- chromPeaks(object)
              cpks <- cbind(cpks, index = seq_len(nrow(cpks)))
              if (!any(colnames(cpks) == "sample"))
                  colnames(cpks)[colnames(cpks) == "column"] <- "sample"
              if (!any(colnames(cpks) == "mz"))
                  cpks <- cbind(cpks, mz = NA)
              nr <- nrow(object)
              res <- vector("list", nr)
              bw <- bw(param)
              sgrps <- sampleGroups(param)
              sgrps_tbl <- table(sgrps)
              minfr <- minFraction(param)
              minsam <- minSamples(param)
              maxf <- maxFeatures(param)
              for (i in seq_len(nr)) {
                  cur_pks <- cpks[cpks[, "row"] == i, , drop = FALSE]
                  if (nrow(cur_pks) == 0)
                      next
                  rtr <- range(lapply(object[i, ], rtime), na.rm = TRUE)
                  densFrom <- rtr[1] - 3 * bw
                  densTo <- rtr[2] + 3 * bw
                  densN <- max(512, 2 * 2^(ceiling(log2(diff(rtr) / (bw / 2)))))
                  tmp <- .group_peaks_density(cur_pks, bw = bw,
                                              densFrom = densFrom,
                                              densTo = densTo,
                                              densN = densN,
                                              sampleGroups = sgrps,
                                              sampleGroupTable = sgrps_tbl,
                                              minFraction = minfr,
                                              minSamples = minsam,
                                              maxFeatures = maxf)
                  tmp$row <- rep(i, nrow(tmp))
                  res[[i]] <- tmp
              }
              res <- DataFrame(do.call(rbind, res))
              xph <- XProcessHistory(param = param, date. = startDate,
                                     type. = .PROCSTEP.PEAK.GROUPING,
                                     fileIndex = 1:ncol(object))
              object <- addProcessHistory(object, xph)
              ## Add the results.
              if (nrow(res) == 0)
                  warning("Unable to group any chromatographic peaks. ",
                          "You might have to adapt your settings.")
              if (nrow(res) > 0)
                  rownames(res) <- .featureIDs(nrow(res))
              object@featureDefinitions <- res
              if (validObject(object))
                  return(object)
          })

#' @rdname XChromatogram
setMethod("featureDefinitions", "XChromatograms",
          function(object, mz = numeric(), rt = numeric(), ppm = 0,
                   type = c("any", "within", "apex_within")) {
              if (!hasFeatures(object))
                  return(DataFrame())
              feat_def <- object@featureDefinitions
              type <- match.arg(type)
              ## Select features within rt range.
              if (length(rt) && nrow(feat_def)) {
                  rt <- range(rt)
                  keep <- switch(
                      type,
                      any = which(feat_def$rtmin <= rt[2] &
                                  feat_def$rtmax >= rt[1]),
                      within = which(feat_def$rtmin >= rt[1] &
                                     feat_def$rtmax <= rt[2]),
                      apex_within = which(feat_def$rtmed >= rt[1] &
                                          feat_def$rtmed <= rt[2])
                  )
                  feat_def <- feat_def[keep, , drop = FALSE]
              }
              if (length(mz) && nrow(feat_def)) {
                  mz <- range(mz)
                  ## Increase mz by ppm.
                  if (is.finite(mz[1]))
                      mz[1] <- mz[1] - mz[1] * ppm / 1e6
                  if (is.finite(mz[2]))
                      mz[2] <- mz[2] + mz[2] * ppm / 1e6
                  keep <- switch(
                      type,
                      any = which(feat_def$mzmin <= mz[2] &
                                  feat_def$mzmax >= mz[1]),
                      within = which(feat_def$mzmin >= mz[1] &
                                     feat_def$mzmax <= mz[2]),
                      apex_within = which(feat_def$mzmed >= mz[1] &
                                          feat_def$mzmed <= mz[2])
                  )
                  feat_def <- feat_def[keep, , drop = FALSE]
              }
              feat_def
          })

#' @rdname XChromatogram
setMethod("[", "XChromatograms", function(x, i, j, drop = FALSE) {
    if (missing(i) && missing(j))
        return(x)
    if (missing(i))
        i <- seq_len(nrow(x))
    if (missing(j))
        j <- seq_len(ncol(x))
    if (is.logical(i))
        i <- which(i)
    if (is.logical(j))
        j <- which(j)
    if (length(i) == 1 && length(j) == 1)
        return(x@.Data[i, j, drop = TRUE][[1]])
    cpeaks_orig <- chromPeaks(x)
    fts <- featureDefinitions(x)
    x <- callNextMethod()
    if (nrow(fts)) {
        rownames(cpeaks_orig) <- as.character(seq_len(nrow(cpeaks_orig)))
        cpks <- .subset_chrom_peaks_xchromatograms(cpeaks_orig, i = i, j = j)
        idxs <- as.integer(rownames(cpks))
        fts <- fts[fts$row %in% i, , drop = FALSE]
        fts$row <- match(fts$row, i)
        fts$peakidx <- lapply(fts$peakidx, function(z) {
            newidx <- match(z, idxs)
            newidx[!is.na(newidx)]
        })
        fts <- fts[lengths(fts$peakidx) > 0, , drop = FALSE]
        x@featureDefinitions <- fts[order(fts$row), , drop = FALSE]
    }
    if (validObject(x))
        x
})

#' Subset the chromPeaks matrix from an `XChromatograms` object. The
#' `chromPeaks` matrix is generated dynamically from the `chromPeaks` matrices
#' of each internal `XChromatogram` object, so there is not really a need to
#' subset the `chromPeaks` from an `XChromatograms` - only that we need this
#' to update the `"peakidx"` column of the `featureDefinitions`.
#'
#' Note: the chromPeaks matrix is extracted ordered by row.
#'
#' @author Johannes Rainer
#'
#' @noRd
.subset_chrom_peaks_xchromatograms <- function(x, i, j) {
    if (missing(i) & missing(j))
        return(x)
    if (missing(i)) i <- seq_len(nrow(x))
    if (missing(j)) j <- seq_len(ncol(x))
    x <- x[x[, "row"] %in% i & x[, "column"] %in% j, , drop = FALSE]
    if (nrow(x)) {
        ## x <- x[order(match(x[, "row"], i), match(x[, "column"], j)), ,
        ##        drop = FALSE]
        x[, "row"] <- match(x[, "row"], i)
        x[, "column"] <- match(x[, "column"], j)
        x[order(x[, "row"], x[, "column"]), , drop = FALSE]
    } else x
}