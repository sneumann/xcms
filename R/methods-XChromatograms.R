#' @include methods-MChromatograms.R

setAs("MChromatograms", "XChromatograms", function(from) {
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
setMethod("hasFilledChromPeaks", "XChromatograms", function(object) {
    matrix(vapply(object, .hasFilledPeaks, logical(1)), ncol = ncol(object),
           dimnames = dimnames(object))
})

#' @rdname XChromatogram
setMethod("chromPeaks", "XChromatograms", function(object, rt = numeric(),
                                                   mz = numeric(), ppm = 0,
                                                   type = c("any", "within",
                                                            "apex_within"),
                                                   msLevel) {
    type <- match.arg(type)
    res <- lapply(object, chromPeaks, rt = rt, mz = mz, ppm = ppm, type = type,
                  msLevel = msLevel)
    nrs <- vapply(res, nrow, integer(1))
    row_idx <- rep(seq_len(nrow(object)), ncol(object))
    col_idx <- rep(seq_len(ncol(object)), each = nrow(object))
    res <- do.call(rbind, res)
    res <- cbind(res, row = rep(row_idx, nrs), column = rep(col_idx, nrs))
    res[order(res[, "row"]), , drop = FALSE]
})

#' @rdname XChromatogram
setMethod("chromPeakData", "XChromatograms", function(object) {
    res <- lapply(object, chromPeakData)
    nrs <- vapply(res, nrow, integer(1))
    row_idx <- rep(seq_len(nrow(object)), ncol(object))
    col_idx <- rep(seq_len(ncol(object)), each = nrow(object))
    res <- do.call(rbind, res)
    res$row <- rep(row_idx, nrs)
    res$column <- rep(col_idx, nrs)
    extractROWS(res, order(res[, "row"]))
})

#' @rdname XChromatogram
setMethod("filterMz", "XChromatograms", function(object, mz, ...) {
    if (missing(mz) || length(object) == 0)
        return(object)
    pks_orig <- chromPeaks(object)
    object@.Data <- matrix(lapply(object, filterMz, mz = mz, ...),
                           nrow = nrow(object), dimnames = dimnames(object))
    pks_sub <- chromPeaks(object)
    if (hasFeatures(object)) {
        fts <- .subset_features_on_chrom_peaks(
            object@featureDefinitions, pks_orig, pks_sub)
        fts$row <- vapply(fts$peakidx, function(z) {
            as.integer(pks_sub[z, "row"][1])
        }, integer(1))
        object@featureDefinitions <- fts
    }
    validObject(object)
    object
})

#' @rdname XChromatogram
setMethod("filterRt", "XChromatograms", function(object, rt, ...) {
    if (missing(rt) || length(object) == 0)
        return(object)
    pks_orig <- chromPeaks(object)
    object@.Data <- matrix(lapply(object, filterRt, rt = rt, ...),
                           nrow = nrow(object), dimnames = dimnames(object))
    pks_sub <- chromPeaks(object)
    if (hasFeatures(object)) {
        fts <- .subset_features_on_chrom_peaks(
            object@featureDefinitions, pks_orig, pks_sub)
        fts$row <- vapply(fts$peakidx, function(z) {
            as.integer(pks_sub[z, "row"][1])
        }, integer(1))
        object@featureDefinitions <- fts
    }
    validObject(object)
    object
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
        x_sub <- x[i, , drop = FALSE]
        plot(as(x_sub, ifelse(is(x_sub, "XChromatograms"),
                              "MChromatograms", "Chromatogram")),
             col = col, lty = lty, type = type,
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
    validObject(object)
    object
})

#' @rdname XChromatogram
#'
#' @section Chromatographic peak detection:
#'
#' See [findChromPeaks-Chromatogram-CentWaveParam] for information.
#'
#' After chromatographic peak detection it is also possible to *refine*
#' identified chromatographic peaks with the `refineChromPeaks` method (e.g. to
#' reduce peak detection artifacts). Currently, only peak refinement using the
#' *merge neighboring peaks* method is available (see
#' [MergeNeighboringPeaksParam()] for a detailed description of the approach.
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
#' The `plotChromPeakDensity` method can be used to visualize *peak density*
#' correspondence results, or to *simulate* a peak density correspondence
#' analysis on chromatographic data. The resulting plot consists of two panels,
#' the upper panel showing the chromatographic data as well as the identified
#' chromatographic peaks, the lower panel the distribution of peaks (the peak
#' density) along the retention time axis. This plot shows each peak as a point
#' with it's peak's retention time on the x-axis, and the sample in which it
#' was found on the y-axis. The distribution of peaks along the retention time
#' axis is visualized with a density estimate. Grouped chromatographic peaks
#' are indicated with grey shaded rectangles. Parameter `simulate` allows to
#' define whether the correspondence analysis should be simulated (
#' `simulate=TRUE`, based on the available data and the provided
#' [PeakDensityParam()] parameter class) or not (`simulate=FALSE`). For the
#' latter it is assumed that a correspondence analysis has been performed with
#' the *peak density* method on the `object`.
#' See examples below.
#'
#' Abundance estimates for each feature can be extracted with the
#' `featureValues` function using parameter `value = "into"` to extract the
#' integrated peak area for each feature. The result is a `matrix`, columns
#' being samples and rows features.
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
              validObject(object)
              object
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
                  feat_def <- extractROWS(feat_def, keep)
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
                  feat_def <- extractROWS(feat_def, keep)
              }
              feat_def
          })

#' @rdname XChromatogram
setMethod("[", "XChromatograms", function(x, i, j, drop = TRUE) {
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
    if (length(i) > 1 || length(j) > 1)
        drop <- FALSE
    if (length(i) == 1 && length(j) == 1 && drop)
        return(x@.Data[i, j, drop = TRUE][[1]])
    cpeaks_orig <- chromPeaks(x)
    fts_orig <- featureDefinitions(x)
    ## The following code replicates the [,MChromatograms
    ph <- x@.processHistory
    pd <- x@phenoData
    fd <- x@featureData
    xclass <- class(x)
    x <- as(x@.Data[i = i, j = j, drop = FALSE], xclass)
    pd <- pd[j, ]
    pData(pd) <- droplevels(pData(pd))
    x@phenoData <- pd
    fd <- fd[i, ]
    pData(fd) <- droplevels(pData(fd))
    x@featureData <- fd
    if (nrow(fts_orig)) {
        cpeaks_sub <- chromPeaks(x)
        ## re-order and duplicate fts based on i.
        fts <- vector("list", length(i))
        for (el in seq_along(i)) {
            fts_row <- fts_orig[fts_orig$row == i[el], , drop = FALSE]
            if (nrow(fts_row)) {
                fts_row$row <- el
                fts_row <- .subset_features_on_chrom_peaks(
                    fts_row, cpeaks_orig, cpeaks_sub)
                fts[[el]] <- fts_row
            } else fts[[el]] <- DataFrame()
        }
        x@featureDefinitions <- do.call(rbind, fts)
    }
    x@.processHistory <- .process_history_subset_samples(ph, j = j)
    validObject(x)
    x
})

#' @rdname XChromatogram
setMethod("featureValues", "XChromatograms",
          function(object, method = c("medret", "maxint", "sum"),
                   value = "into", intensity = "into", missing = NA, ...) {
              if (!any(hasChromPeaks(object)))
                  stop("No chromatographic peaks present! Please use ",
                       "'findChromPeaks' first.")
              if (!hasFeatures(object))
                  stop("No features (grouped peaks) present! Please use ",
                       "'groupChromPeaks' first.")
              method = match.arg(method)
              if (method == "sum" & !(value %in% c("into", "maxo")))
                  stop("method 'sum' is only allowed if value is set to 'into'",
                       " or 'maxo'")
              if (is.character(missing)) {
                  if (!(missing %in% c("rowmin_half")))
                      stop("if 'missing' is not 'NA' or a numeric it should",
                           " be one of: \"rowmin_half\".")
              } else {
                  if (!is.numeric(missing) & !is.na(missing))
                      stop("'missing' should be either 'NA', a numeric or one",
                           " of: \"rowmin_half\".")
              }
              cnames <- colnames(object)
              pks <- chromPeaks(object)
              if (any(colnames(pks) == "sample"))
                  pks[, "sample"] <- pks[, "column"]
              else
                  pks <- cbind(pks, sample = pks[, "column"])
              .feature_values(pks = pks, fts = featureDefinitions(object),
                              method = method, value = value,
                              intensity = intensity, colnames = cnames,
                              missing = missing)
})

#' @rdname XChromatogram
setMethod("plotChromPeakDensity", "XChromatograms",
          function(object, param, col = "#00000060", xlab = "retention time",
                   main = NULL, peakType = c("polygon", "point", "rectangle",
                                             "none"), peakCol = "#00000060",
                   peakBg = "#00000020", peakPch = 1, simulate = TRUE, ...) {
              peakType <- match.arg(peakType)
              if (!any(hasChromPeaks(object)))
                  stop("No chromatographic peaks present. Please run ",
                       "'findChromPeaks' first.", call. = FALSE)
              if (nrow(object) > 1)
                  stop("Currently only plotting of a single chromatogram in ",
                       "multiple samples is supported. Please subset 'object' ",
                       "to one row.", call. = FALSE)
              if (missing(param)) {
                  param <- NULL
                  if (hasFeatures(object)) {
                      ph <- processHistory(object,
                                           type = .PROCSTEP.PEAK.GROUPING)
                      if (length(ph)) {
                          ph <- ph[[length(ph)]]
                          if (is(ph, "XProcessHistory") &&
                              is(ph@param, "PeakDensityParam"))
                              param <- ph@param
                      }
                  }
              }
              if (!length(param))
                  stop("Object 'param' is missing", call. = FALSE)
              fts <- NULL
              if (!simulate && hasFeatures(object))
                  fts <- featureDefinitions(object)
              mr <- par("mar")
              mr_1 <- mr[1]
              mr_3 <- mr[3]
              mr[1] <- 0
              xl <- range(lapply(object, function(z) range(rtime(z))))
              par(mfrow = c(2, 1), mar = mr)
              plot(object, col = col, xlab = NA, xaxt = "n", main = main,
                   peakType = peakType, peakCol = peakCol, peakBg = peakBg,
                   peakPch = peakPch, xlim = xl, ...)
              mr[1] <- mr_1
              mr[3] <- 0
              par(mar = mr)
              .plot_chrom_peak_density(chromPeaks(object), fts = fts, col = col,
                                       param = param, xlab = xlab, xlim = xl,
                                       peakCol = peakCol, peakBg = peakBg,
                                       peakPch = peakPch, simulate = simulate,
                                       ylim = c(1, ncol(object)),
                                       ...)
              mr[1] <- mr_1
              mr[3] <- mr_3
              par(mar = mr)
})

#' @rdname XChromatogram
setMethod("dropFilledChromPeaks", "XChromatograms", function(object) {
    pks_orig <- chromPeaks(object)
    object@.Data <- matrix(lapply(object, dropFilledChromPeaks),
                           nrow = nrow(object), dimnames = dimnames(object))
    pks_sub <- chromPeaks(object)
    if (hasFeatures(object)) {
        fts <- .subset_features_on_chrom_peaks(
            object@featureDefinitions, pks_orig, pks_sub)
        ## fts$row <- vapply(fts$peakidx, function(z) {
        ##     as.integer(pks_sub[z, "row"][1])
        ## }, integer(1))
        object@featureDefinitions <- fts
    }
    object <- dropProcessHistories(object, type = .PROCSTEP.PEAK.FILLING)
    validObject(object)
    object
})

#' @rdname XChromatogram
setMethod("refineChromPeaks", c(object = "XChromatograms",
                                param = "MergeNeighboringPeaksParam"),
          function(object, param = MergeNeighboringPeaksParam()) {
              object@.Data <- matrix(
                  lapply(object, .xchrom_merge_neighboring_peaks,
                         diffRt = 2 * param@expandRt,
                         minProp = param@minProp),
                  ncol = ncol(object), nrow = nrow(object),
                  dimnames = dimnames(object))
              xph <- XProcessHistory(param = param, date. = date(),
                                     type. = .PROCSTEP.PEAK.REFINEMENT,
                                     fileIndex = 1:ncol(object))
              object <- addProcessHistory(object, xph)
              validObject(object)
              object
          })

#' @rdname filter-MChromatograms
setMethod("filterColumnsIntensityAbove", "XChromatograms",
          function(object, threshold = 0,
                   value = c("bpi", "tic", "maxo", "into"),
                   which = c("any", "all")) {
              value <- match.arg(value)
              which <- match.arg(which)
              if (length(threshold) > 1 || !is.numeric(threshold))
                  stop("'threshold' should be a 'numeric' of length 1")
              if (value %in% c("maxo", "into")) {
                  nc <- ncol(object)
                  rws <- seq_len(nrow(object))
                  cps <- chromPeaks(object)
                  keep <- rep(FALSE, nc)
                  for (i in seq_len(nc)) {
                      vals <- cps[cps[, "column"] == i &
                                  cps[, value] > threshold, "row"]
                      if (length(vals)) {
                          if (which == "any")
                              keep[i] <- TRUE
                          else keep[i] <- all(rws %in% vals)
                      }
                  }
                  object[, keep]
              } else
                  callNextMethod(object, threshold = threshold, value = value,
                                 which = which)
          })

#' @rdname filter-MChromatograms
setMethod("filterColumnsKeepTop", "XChromatograms",
          function(object, n = 1L, sortBy = c("bpi", "tic", "maxo", "into"),
                   aggregationFun = sum) {
              sortBy <- match.arg(sortBy)
              if (length(n) > 1 || !is.numeric(n))
                  stop("'n' should be an 'integer' of length 1")
              if (sortBy %in% c("maxo", "into")) {
                  n <- ceiling(n)
                  nc <- ncol(object)
                  if (n > nc)
                      stop("'n' should be smaller or equal than the number of ",
                           "columns (", nc, ")")
                  colval <- numeric(nc)
                  cps <- chromPeaks(object)
                  for (i in seq_len(nc)) {
                      vals <- cps[cps[, "column"] == i, c("row", sortBy),
                                  drop = FALSE]
                      if (nrow(vals)) {
                          vals <- sapply(split(vals[, sortBy], vals[, "row"]),
                                         max, na.rm = TRUE)
                          colval[i] <- aggregationFun(vals, na.rm = TRUE)
                      }
                  }
                  idx <- order(colval, decreasing = TRUE)[seq_len(n)]
                  object[, sort(idx)]
              } else
                  callNextMethod(object, n = n, sortBy = sortBy,
                                 aggregationFun = aggregationFun)
          })

#' @rdname XChromatogram
setMethod("filterChromPeaks", "XChromatograms",
          function(object, method = c("keepTop"), ...) {
              method <- match.arg(method)
              pks_orig <- chromPeaks(object)
              object@.Data <- matrix(lapply(object, filterChromPeaks,
                                            method = method, ...),
                                     nrow = nrow(object),
                                     dimnames = dimnames(object))
              pks_sub <- chromPeaks(object)
              if (hasFeatures(object)) {
                  fts <- .subset_features_on_chrom_peaks(
                      object@featureDefinitions, pks_orig, chromPeaks(object))
                  object@featureDefinitions <- fts
              }
              validObject(object)
              object
          })

#' @rdname XChromatogram
setMethod("transformIntensity", "XChromatograms", function(object,
                                                           FUN = identity) {
    object@.Data <- matrix(lapply(object, FUN = transformIntensity, FUN),
                           nrow = nrow(object),
                           dimnames = dimnames(object))
    object
})
