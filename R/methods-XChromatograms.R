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
        cat(" ", nrow(object@featureDefinitions), " features identified.\n",
            sep = "")
        if (.hasFilledPeaks(object)) {
            totF <- chromPeaks(object)[, "is_filled"] == 1
            fp <- chromPeaks(object)[totF, , drop = FALSE]
            cat("", sum(totF), "filled peaks (on average",
                mean(table(fp[, "sample"])), "per sample).\n")
        }
    }
})

#' @rdname XChromatogram
setMethod("hasChromPeaks", "XChromatograms", function(object) {
    matrix(vapply(object, hasChromPeaks, logical(1)), ncol = ncol(object),
           dimnames = dimnames(object))
})

#' @rdname XChromatogram
#'
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
#' with a [PeakDensityParam()] object. The correspondence analysis results are
#' stored in the returned `XChromatograms` object and can be accessed with the
#' [featureDefinitions()] method.
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
