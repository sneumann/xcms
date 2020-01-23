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
                                            minCor = 0.8, fromFile = 1L,
                                            column = "maxo",
                                            peakId = rownames(
                                                chromPeaks(object, msLevel = 1L))) {
    if (hasAdjustedRtime(object))
        fData(object)$retentionTime <- rtime(object)
    idx <- match(peakId, rownames(chromPeaks(object)))
    idx <- idx[!is.na(idx)]
    if (!length(idx)) {
        return(Spectra(elementMetadata = DataFrame(
                           peak_id = character(),
                           ms2_peak_id = CharacterList(),
                           ms2_peak_cor = NumericList())))
    }
    res <- vector("list", length(idx))
    message("Reconstructing MS2 spectra for ", length(idx), " chrom peaks ...",
            appendLF = FALSE)
    for (i in seq_along(idx)) {
        res[[i]] <- .reconstruct_ms2_for_chrom_peak(
            chromPeaks(object)[idx[i], ], object = object, expandRt = expandRt,
            diffRt = diffRt, minCor = minCor, fromFile = fromFile,
            column = column)
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
#' @param column `character(1)` specifying the column from the chrom peak matrix
#'     from which the intensity values for the reconstructed spectrum should be
#'     taken (recommended, either `"maxo"` or `"into"`).
#'
#' @return [Spectra] object with the reconstructed MS2 spectrum. The spectrum
#'     is empty if no MS2 chromatographic peak with a good enough correlation
#'     was found for the MS1 chromatographic peak. The `Spectra` contains
#'     metadata columns `"ms2_peak_id"` and `"ms2_peak_cor"` of type
#'     [CharacterList()] and [NumericList()] (length equal to number of peaks
#'     per spectrum) providing the IDs and the correlation of the MS2
#'     chromatographic peaks from which the MS2 spectrum was reconstructed.
#'
#' @author Johannes Rainer, Micheal Witting
#'
#' @noRd
.reconstruct_ms2_for_chrom_peak <- function(x, object, fromFile = 1L,
                                            expandRt = 2, diffRt = 5,
                                            minCor = 0.8, column = "maxo") {
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
            elementMetadata = DataFrame(
                ms2_peak_id = CharacterList(character()),
                ms2_peak_cor = NumericList(numeric()))))
    ## extract ion chromatograms for all chromatographic peaks
    pks <- chromPeaks(object)[idx, , drop = FALSE]
    pks <- pks[order(pks[, "mz"]), , drop = FALSE]
    rts_1 <- c(x["rtmin"] - expandRt, x["rtmax"] + expandRt)
    mzs_1 <- x[c("mzmin", "mzmax")]
    rts_2 <- cbind(pks[, "rtmin"] - expandRt, pks[, "rtmax"] + expandRt)
    mzs_2 <- pks[, c("mzmin", "mzmax")]
    chr_1 <- chromatogram(as(object, "OnDiskMSnExp"), mz = mzs_1, rt = rts_1)
    chr_2 <- chromatogram(od_object, mz = mzs_2, rt = rts_2, msLevel = 2L)
    ## Filter for chromatograms that are empty or have too few data points
    chr_2 <- chr_2[vapply(chr_2, function(z) sum(!is.na(z@intensity)),
                          integer(1)) > 3, ]
    if (!length(chr_2))
        return(Spectra(
            new("Spectrum2", fromFile = fromFile),
            elementMetadata = DataFrame(
                ms2_peak_id = CharacterList(character()),
                ms2_peak_cor = NumericList(numeric()))))
    cors <- vapply(chr_2@.Data, .correlate_chromatogram, y = chr_1[1, 1],
                   numeric(1), align = "approx")
    pks <- pks[which(cors >= minCor), , drop = FALSE]
    if (nrow(pks)) {
        sp <- new("Spectrum2", fromFile = fromFile, centroided = TRUE,
                  mz = pks[, "mz"], intensity = pks[, column],
                  precursorMz = x["mz"])
        df <- DataFrame(matrix(ncol = 0, nrow = 1))
        df$ms2_peak_id <- CharacterList(rownames(pks), compress = FALSE)
        df$ms2_peak_cor <- NumericList(cors[cors >= minCor], compress = FALSE)
        Spectra(sp, elementMetadata = df)
    } else {
        Spectra(
            new("Spectrum2", fromFile = fromFile),
            elementMetadata = DataFrame(
                ms2_peak_id = CharacterList(character()),
                ms2_peak_cor = NumericList(numeric())))
    }
}
