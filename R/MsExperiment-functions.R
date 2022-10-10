.param_to_fun <- function(x) {
    p2f <- c(CentWaveParam = "do_findChromPeaks_centWave",
             MatchedFilterParam = "do_findChromPeaks_matchedFilter",
             MassifquantParam = "do_findChromPeaks_massifquant",
             MSWParam = "do_findPeaks_MSW",
             CentWavePredIsoParam = "do_findChromPeaks_centWaveWithPredIsoROIs")
    fun <- p2f[class(x)[1L]]
    if (is.na(fun))
        stop("No peak detection function for parameter class ", class(x)[1L])
    unname(fun)
}

#' Apply a function FUN to the data of each sample in an MsExperiment.
#'
#' @author Johannes Rainer
#'
#' @noRd
.mse_sample_apply <- function(x, FUN = identity, ..., BPPARAM = bpparam()) {
    bplapply(split(x, seq_along(x)), function(x, ...) {
        FUN(x, ...)
    }, ..., BPPARAM = BPPARAM)
}

#' A faster way to apply a function to data from one file/sample than
#' `.mse_sample_apply` that does **only** split the `Spectra` and no other data
#' the `MsExperiment`.
#'
#' This function splits the spectra by x@sampleDataLinks [, 1] which is ALWAYS
#' ordered from 1 to number of samples. The result will thus be a `list` in the
#' same order (sample 1:length(x)).
#'
#' @author Johannes Rainer
#'
#' @noRd
.mse_sample_spectra_apply <- function(x, FUN = identity, ...,
                                      BPPARAM = SerialParam()) {
    if (!length(spectra(x)))
        stop("No spectra available.")
    if (!any(names(x@sampleDataLinks) == "spectra"))
        stop("No link between samples and spectra found. Please use ",
             "'linkSampleData' to define which spectra belong to which ",
             "samples.")
    if (nrow(x@sampleDataLinks[["spectra"]]) != length(spectra(x)))
        stop("All spectra in 'x' need to be associated with a sample.")
    Spectra::spectrapply(
                 spectra(x),
                 f = as.factor(x@sampleDataLinks[["spectra"]][, 1L]),
                 FUN = FUN, ..., BPPARAM = BPPARAM)
}

#' @title Perform peak detection on a `Spectra` object
#'
#' @description
#'
#' Performs peak detection on a `Spectra` object.
#'
#' @param x `Spectra` with spectra from a **single sample/file**.
#'
#' @param msLevel `integer(1)` defining on which MS level to do peak detection.
#'
#' @param param parameter object with the settings for the peak detection.
#'
#' @param ... ignored.
#'
#' @author Johannes Rainer
#' @noRd
.mse_sample_find_chrom_peaks <- function(x, msLevel = 1L, param, ...) {
    x <- filterMsLevel(x, msLevel)
    pkd <- Spectra::peaksData(x, c("mz", "intensity"), BPPARAM = SerialParam())
    vals_per_spect <- vapply(pkd, nrow, integer(1), USE.NAMES = FALSE)
    ## Open questions:
    ## - What to do with empty spectra? Remove them? MatchFilter does not like
    ##   them. Maybe add a matrix with m/z = 0, rt = 0.
    if (any(vals_per_spect == 0))
        warning("Found empty spectra. Please run 'filterEmptySpectra' first.",
                call. = FALSE)
    pkd <- do.call(rbind, pkd)
    if (!length(pkd))
        return(matrix(numeric(), ncol =  9, nrow = 0,
                      dimnames = list(character(), c("mz", "mzmin", "mzmax",
                                                     "rt", "rtmin", "rtmax",
                                                     "into", "maxo", "sn"))))
    if (inherits(param, "CentWaveParam")) {
        centroided <- all(centroided(x))
        if (is.na(centroided)) {
            centroided <- isCentroided(x[ceiling(length(x) / 3)])
            if (is.na(centroided) || !centroided)
                warning("Your data appears to be not centroided! CentWave",
                        " works best on data in centroid mode.")
        }
    }
    rts <- rtime(x)
    if (is.unsorted(rts))
        stop("Spectra are not ordered by retention time", .call = FALSE)
    do.call(.param_to_fun(param),
            args = c(list(mz = pkd[, 1L], int = pkd[, 2L], scantime = rts,
                          valsPerSpect = vals_per_spect), as(param, "list")))
}

#' Perform peak detection on an MsExperiment object and returns the `matrix`
#' with identified chromatographic peaks.
#'
#' @author Johannes Rainer
#'
#' @noRd
.mse_find_chrom_peaks <- function(x, msLevel = 1L, param, ...,
                                  BPPARAM = bpparam()) {
    res <- .mse_sample_spectra_apply(x, FUN = .mse_sample_find_chrom_peaks,
                                     msLevel = msLevel, param = param,
                                     BPPARAM = BPPARAM)
    sidx <- vapply(res, nrow, integer(1), USE.NAMES = FALSE)
    cbind(do.call(rbind, res), sample = rep(seq_along(x), sidx))
}
