## functions for SWATH/DIA analysis.
## Chromatogram alignment and correlation functions have been moved to
## functions-Chromatogram.R

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

#' @description
#'
#' Reconstructs MS2 spectra for each MS1 chromatographic peak (if possible) for
#' DIA data (such as SWATH).
#'
#' @param object `XCMSnExp` with identified chromatographic peaks.
#'
#' @return [Spectra()] with the reconstructed MS2 spectra for all MS1 peaks
#'     in `object`. Contains empty [Spectrum2-class] objects for MS1 peaks for
#'     which reconstruction was not possible.
#'
#' @author Johannes Rainer, Micheal Witting
#'
#' @noRd
reconstructChromPeakSpectra <- function(object, expandRt = 2, diffRt = 4,
                                        minCor = 0.8, BPPARAM = bpparam()) {
    ## Input check: need MS1 peaks and MS2 peaks.

    ## Drop featureDefinitions if present
    if (hasFeatures(object))
        object <- dropFeatureDefinitions(object, keepAdjustedRtime = TRUE)
    object <- selectFeatureData(object,
                                fcol = c(MSnbase:::.MSnExpReqFvarLabels,
                                         "centroided",
                                         "polarity",
                                         "isolationWindow",
                                         "isolationWindowTargetMZ",
                                         "isolationWindowLowerOffset",
                                         "isolationWindowUpperOffset"))
    sps <- bplapply(
        lapply(seq_len(length(fileNames(object))), filterFile, object = object,
               keepAdjustedRtime = TRUE),
        FUN = function(x, files, expandRt, diffRt, minCor) {
            .reconstruct_ms2_for_peaks_file(
                x, expandRt = expandRt, diffRt = diffRt,
                minCor = minCor, fromFile = match(fileNames(x), files))
        },
        BPPARAM = BPPARAM, files = fileNames(object), expandRt = expandRt,
        diffRt = diffRt, minCor = minCor)
    do.call(c, sps)
}

#' @description
#'
#' Reconstruct MS2 spectra for all MS1 chromatographic peaks of one file.
#' See `.reconstruct_ms2_for_chrom_peak` for details and parameters.
#'
#' @return `Spectra` with the reconstructed MS2 spectra (emtpy `Spectra` if
#'     no spectrum was reconstructed for any peak or if no MS1 chromatographic
#'     peaks are present).
#'
#' @author Johannes Rainer, Micheal Witting
#'
#' @noRd
.reconstruct_ms2_for_peaks_file <- function(object, expandRt = 2, diffRt = 5,
                                            minCor = 0.8, fromFile = 1L) {
    if (hasAdjustedRtime(object))
        fData(object)$retentionTime <- rtime(object)
    idx <- which(chromPeakData(object)$ms_level == 1L)
    if (!length(idx)) {
        return(Spectra(elementMetadata = DataFrame(peak_id = character(),
                                                   ms2_peak_id = SimpleList())))
    }
    res <- vector("list", length(idx))
    message("Reconstructing MS2 spectra for ", length(idx), " chrom peaks ...",
            appendLF = FALSE)
    for (i in idx) {
        res[[i]] <- .reconstruct_ms2_for_chrom_peak(
            chromPeaks(object)[i, ], object = object, expandRt = expandRt,
            diffRt = diffRt, minCor = minCor, fromFile = fromFile)
    }
    message(" OK")
    res <- do.call(c, res)
    mcols(res)$peak_id <- rownames(chromPeaks(object))[idx]
    ## res[!isEmpty(res)]
    res
}

#' @description
#'
#' Reconstruct the MS2 spectrum for a chromatographic peak from a set of MS2
#' chromatographic peaks within the isolation window matching the MS1 peak's
#' m/z.
#'
#' @param x `numeric` representing one line of a `chromPeaks` `matrix`.
#'
#' @param object `XCMSnExp` with M2 chromatographic peaks identified in
#'     isolation windows.
#'
#' @param fromFile optional `integer(1)` to be used for slot `fromFile` in the
#'     `Spectrum2` object.
#'
#' @param expandRt `numeric(1)` allowing to expand the retention time range
#'     for extracted ion chromatograms by a constant value.
#'
#' @param diffRt `numeric(1)` defining the maximal allowed difference between
#'     the retention time of the chromatographic peak (apex) and the retention
#'     times of MS2 chromatographic peaks (apex) to consider them as
#'     representing candidate fragments of the original ion.
#'
#' @param minCor `numeric(1)` defining the minimal required correlation
#'     coefficient for MS2 chromatographic peaks to be considered for MS2
#'     spectrum reconstruction.
#'
#' @return [Spectra] object with the reconstructed MS2 spectrum. The spectrum
#'     is empty if no MS2 chromatographic peak with a good enough correlation
#'     was found for the MS1 chromatographic peak.
.reconstruct_ms2_for_chrom_peak <- function(x, object, fromFile = 1L,
                                            expandRt = 2, diffRt = 5,
                                            minCor = 0.8) {
    if (is.matrix(x))
        x <- x[1, ]
    ## find MS2 chrom peaks from the isolation window matching the peak's m/z
    idx_mz <- .which_mz_in_range(x["mz"],
                                 chromPeakData(object)$isolationWindowLowerMz,
                                 chromPeakData(object)$isolationWindowUpperMz)
    ## reduce to peaks with overlapping retention times.
    idx_rt <- .which_chrom_peak_diff_rt(x, chromPeaks(object), diffRt = diffRt)
    idx <- intersect(idx_mz, idx_rt)
    ## Only look into the peaks/spectra of the correct isolation window
    od_object <- filterIsolationWindow(as(object, "OnDiskMSnExp"), mz = x["mz"])
    if (!length(idx) | !length(od_object))
        return(Spectra(
            new("Spectrum2", fromFile = fromFile),
            elementMetadata = DataFrame(ms2_peak_id = SimpleList(NA_character_))))
    ## extract ion chromatograms for all chromatographic peaks
    pks <- chromPeaks(object)[idx, , drop = FALSE]
    rts_1 <- c(x["rtmin"] - expandRt, x["rtmax"] + expandRt)
    mzs_1 <- x[c("mzmin", "mzmax")]
    rts_2 <- cbind(pks[, "rtmin"] - expandRt, pks[, "rtmax"] + expandRt)
    mzs_2 <- pks[, c("mzmin", "mzmax")]
    chr_1 <- chromatogram(as(object, "OnDiskMSnExp"), mz = mzs_1, rt = rts_1)
    chr_2 <- chromatogram(od_object, mz = mzs_2, rt = rts_2, msLevel = 2L)
    cors <- vapply(chr_2@.Data, xcms:::.correlate_chromatogram, y = chr_1[1, 1],
                   numeric(1), align = "approx")
    pks <- pks[which(cors >= minCor), , drop = FALSE]
    if (nrow(pks)) {
        warning("Need to plug in the reconstruction function.")
        sp <- new("Spectrum2", fromFile = fromFile)
        df <- DataFrame(matrix(ncol = 0, nrow = 1))
        df$ms2_peak_id <- SimpleList(rownames(pks))
        Spectra(sp, elementMetadata = df)
    } else {
        Spectra(
            new("Spectrum2", fromFile = fromFile),
            elementMetadata = DataFrame(ms2_peak_id = SimpleList(NA_character_)))
    }
}
