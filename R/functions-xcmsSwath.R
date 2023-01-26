## functions for SWATH/DIA analysis.

#' Get all MS2 peaks from the isolation window containing the MS1 peaks'
#' m/z.
#'
#' @noRd
.which_mz_in_range <- function(mz, lowerMz, upperMz) {
    if (length(mz) > 1)
        return(lapply(mz, .which_mz_in_range, lowerMz, upperMz))
    which(mz >= lowerMz & mz <= upperMz)
}

#' @description
#'
#' Which rows (peaks) in `pks` have overlapping retention time ranges with
#' the peak `x`?
#'
#' @param x `numeric` representing one row from the `chromPeaks` matrix.
#'
#' @param pks `matrix` representing a `chromPeaks` matrix.
#'
#' @return `integer` with the index of rows (peaks) in `pks` overlapping the
#'     retention time range of `x`.
#'
#' @author Johannes Rainer
#'
#' @noRd
.which_chrom_peak_overlap_rt <- function(x, pks) {
    if (is.matrix(x))
        x <- x[1, ]
    which(pks[, "rtmin"] <= x["rtmax"] & pks[, "rtmax"] >= x["rtmin"])
}

#' @description
#'
#' Which rows (peaks) in `pks` have a retention time (of the apex) that is
#' close to the apex of the specified peak `x`. Peaks are considerd close
#' if the difference between their apex retention time and the retention time
#' of the input peak is smaller than `diffRt`.
#'
#' @param x `numeric` representing one row from the `chromPeaks` matrix.
#'
#' @param pks `matrix` representing a `chromPeaks` matrix.
#'
#' @return `integer` with the index of the rows (peaks) in `pks` that are close
#'     to the specified peak `x`.
#'
#' @author Johannes Rainer
#'
#' @noRd
.which_chrom_peak_diff_rt <- function(x, pks, diffRt = 2) {
    if (is.matrix(x))
        x <- x[1, ]
    which(abs(pks[, "rt"] - x["rt"]) <= diffRt)
}

#' *Reconstruct* MS2 spectra for DIA data:
#' For each MS1 chromatographic peak:
#'
#' - find all MS2 chrom peaks from the same isolation window
#' - reduce to MS2 chrom peaks with an rt similar to the one from the MS1
#'   chrom peak
#' - remove EICs with 2 or less data points
#' - create an MS2 spectrum from all MS2 chrom peaks with peak shape
#'   correlation > `minCor`.
#'
#' @param object `XCMSnExp` with data from a **single** file.
#'
#' @note
#'
#' this function first extracts EICs for all chromatographic peaks, thus it
#' will not be efficient for predicting MS2 spectra for selected MS1 peaks.
#'
#' Be aware that this function does only support returning a `Spectra`!
#'
#' @noRd
.reconstruct_dia_ms2 <-
    function(object, expandRt = 2, diffRt = 5, minCor = 0.8, fromFile = 1L,
             column = "maxo",
             peakId = rownames(chromPeaks(object, msLevel = 1L))) {
        if (hasAdjustedRtime(object))
            fData(object)$retentionTime <- rtime(object)
        message("Reconstructing MS2 spectra for ", length(peakId),
                " chrom peaks ...", appendLF = FALSE)
        pks <- chromPeaks(object)[, c("mz", "mzmin", "mzmax", "rt", "rtmin",
                                      "rtmax", column)]
        pks[, "rtmin"] <- pks[, "rtmin"] - expandRt
        pks[, "rtmax"] <- pks[, "rtmax"] + expandRt
        ord <- order(pks[, "mz"])       # m/z need to be ordered in a Spectra
        pks <- pks[ord, ]
        ilmz <- chromPeakData(object)$isolationWindowLowerMz[ord]
        iumz <- chromPeakData(object)$isolationWindowUpperMz[ord]
        ## Get EICs for all chrom peaks (all MS levels)
        object <- filterRt(object, rt = range(pks[, c("rtmin", "rtmax")]))
        chrs <- .chromatograms_for_peaks(
            lapply(spectra(object),
                   function(z) cbind(mz = z@mz, intensity = z@intensity)),
            rt = rtime(object), msl = msLevel(object), file_idx = fromFile,
            tmz = isolationWindowTargetMz(object), pks = pks,
            pks_msl = chromPeakData(object)$ms_level[ord],
            pks_tmz = chromPeakData(object)$isolationWindowTargetMZ[ord])
        idx <- match(peakId, rownames(pks)) # MS1 peaks to loop over
        res <- data.frame(
            peak_id = peakId, precursorMz = pks[idx, "mz"],
            rtime = pks[idx, "rt"], msLevel = 2L,
            polarity = polarity(object)[1L],
            precursorIntensity = pks[idx, column],
            fromFile = fromFile)
        ms2_peak_id <- lapply(idx, function(z) character())
        mzs <- ints <- ms2_peak_cor <-
            lapply(seq_len(nrow(res)), function(z) numeric())
        for (i in seq_along(idx)) {
            ii <- idx[i]
            imz <- .which_mz_in_range(pks[ii, "mz"], ilmz, iumz)
            irt <- .which_chrom_peak_diff_rt(pks[ii, c("rt", "rtmax")],
                                             pks, diffRt = diffRt)
            ix <- intersect(imz, irt)
            if (!length(ix))
                next
            ## Filter empty or sparse chromatograms
            ix <- ix[vapply(
                chrs@.Data[ix], function(z) sum(!is.na(z@intensity)),
                integer(1)) > 2]
            if (!length(ix))
                next
            ## Correlate
            cors <- vapply(
                chrs@.Data[ix], compareChromatograms, numeric(1),
                y = chrs[[ii]], ALIGNFUNARGS = list(method = "approx"))
            keep <- which(cors >= minCor)
            if (length(keep)) {
                ix <- ix[keep]
                res$rtime[i] <- median(pks[ix, "rt"])
                ms2_peak_id[[i]] <- rownames(pks)[ix]
                ms2_peak_cor[[i]] <- unname(cors[keep])
                mzs[[i]] <- unname(pks[ix, "mz"])
                ints[[i]] <- unname(pks[ix, column])
            }
        }
        res$mz <- mzs
        res$intensity <- ints
        res$ms2_peak_id <- ms2_peak_id
        res$ms2_peak_cor <- ms2_peak_cor
        message(" OK")
        .require_spectra()
        Spectra::Spectra(res)
    }

#' Add a function that takes two chromatograms and returns a `logical(1)`
#' whether the two should be merged? Or better maybe some data.frame? That
#' function could then be called by the *framework* function for SWATH data
#' (i.e. identify potentially matching chromatograms and calculating stuff
#' on them.
#'
#' @noRd
NULL
