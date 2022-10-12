#' Adding chrom peaks to an XcmsExperiment. If chrom peaks are present they
#' will be added. Rownames will be fixed/added too.
#'
#' @noRd
.mse_add_chrom_peaks <- function(object, pks, pkd) {
    ## Keeping also old chromPeaks.
    nr <- nrow(object@chromPeaks)
    if (nr) {
        ## Keep the original chrom peak names.
        idx <- max(as.integer(sub("CP", "", rownames(object@chromPeaks))))
        rownames(pks) <- rownames(pkd) <- .featureIDs(nrow(pks), prefix = "CP",
                                                      from = idx + 1L)
        object@chromPeaks <- suppressWarnings(
            rbindFill(object@chromPeaks, pks))
        object@chromPeakData <- suppressWarnings(
            rbindFill(object@chromPeakData, pkd))
    } else {
        rownames(pks) <- rownames(pkd) <- .featureIDs(nrow(pks), prefix = "CP")
        object@chromPeaks <- pks
        object@chromPeakData <- pkd
    }
    object
}

.mse_valid_chrom_peaks <- function(x) {
    if (!all(.REQ_PEAKS_COLS %in% colnames(x)))
        "One or more required columns missing in 'chromPeaks'"
    else NULL
}

.mse_valid_chrom_peak_data <- function(x) {
    if (!all(c("ms_level", "is_filled") %in% colnames(x)))
        "One or more required columns missing in 'chromPeakData'"
    else NULL
}

.mse_same_rownames <- function(a, b) {
    if (!all(rownames(a) == rownames(b)))
        "Row names of 'chromPeaks' and 'chromPeakData' don't match"
    else NULL
}
