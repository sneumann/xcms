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
    .mse_check_spectra_sample_mapping(x)
    Spectra::spectrapply(
                 spectra(x),
                 f = as.factor(x@sampleDataLinks[["spectra"]][, 1L]),
                 FUN = FUN, ..., BPPARAM = BPPARAM)
}

#' @title Perform peak detection on a `Spectra` object
#'
#' @description
#'
#' Performs peak detection on a `Spectra` object which is supposed to contain
#' spectra from a single file/sample.
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
.mse_find_chrom_peaks_sample <- function(x, msLevel = 1L, param, ...) {
    x <- filterMsLevel(x, msLevel)
    pkd <- Spectra::peaksData(x, columns = c("mz", "intensity"),
                              BPPARAM = SerialParam())
    vals_per_spect <- vapply(pkd, nrow, integer(1), USE.NAMES = FALSE)
    ## Open questions:
    ## - What to do with empty spectra? Remove them? MatchFilter does not like
    ##   them. Maybe add a matrix with m/z = 0, rt = 0.
    if (any(vals_per_spect == 0))
        warning("Found empty spectra. Please run 'filterEmptySpectra' first.",
                call. = FALSE)
    pkd <- do.call(rbind, pkd)
    if (!length(pkd))
        return(NULL)                    # not returning matrix because of rbind
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
    res <- .mse_sample_spectra_apply(x, FUN = .mse_find_chrom_peaks_sample,
                                     msLevel = msLevel, param = param,
                                     BPPARAM = BPPARAM)
    sidx <- vapply(res, function(z) if (is.matrix(z)) nrow(z) else length(z),
                   integer(1), USE.NAMES = FALSE)
    res <- cbind(do.call(rbind, res), sample = rep(seq_along(x), sidx))
    if (!nrow(res))
        .empty_chrom_peaks()
    else res
}

#' Helper function to process spectra in an MsExperiment in chunks, rather than
#' directly in parallel.
#' This function assigns each spectrum the index of the sample it belongs
#' (as a new spectra variable). This can be used by FUN to split the spectra
#' by sample and process them separately (and in parallel).
#'
#' @author Johannes Rainer
#'
#' @noRd
.mse_spectrapply_chunks <- function(x, FUN, ..., chunkSize = 1L,
                                    BPPARAM = bpparam) {
    if (!length(spectra(x)))
        stop("No spectra available.")
    if (!any(names(x@sampleDataLinks) == "spectra"))
        stop("No link between samples and spectra found. Please use ",
             "'linkSampleData' to define which spectra belong to which ",
             "samples.")
    idx <- seq_along(x)
    chunks <- split(idx, ceiling(idx / chunkSize))
    sps <- spectra(x)[x@sampleDataLinks[["spectra"]][, 2L]]
    sps$.SAMPLE_IDX <- x@sampleDataLinks[["spectra"]][, 1L] # or as.factor?
    lapply(chunks, function(z, ...) {
        FUN(sps[sps$.SAMPLE_IDX %in% z], ...)
    }, ..., BPPARAM = BPPARAM)
    ## res <- vector("list", length(chunks))
    ## for (i in seq_along(chunks)) {
    ##     res[[i]] <- FUN(sps[sps$.SAMPLE_IDX %in% chunks[[i]]], ...,
    ##                     BPPARAM = BPPARAM)
    ## }
    ## res
}

#' @title Perform peak detection on chunks of data
#'
#' @description
#'
#' Perform the peak detection in chunks of samples along `x`. Data retrieval is
#' performed on `chunkSize` samples simultaneous and peak detection is then
#' performed in parallel (separate per file). The value of `chunkSize`
#' determines thus how much memory will be needed in each iteration. Also,
#' parallel processing will only be performed on at maximum `chunkSize` cores.
#' Thus, this parameter should be chosed wisely to ensure that a) memory usage
#' is within the available memory on the system and that b) the parallel
#' processing setup defined by `BPPARAM` is efficiently used. Defining a
#' parallel processing setup with 10 cores but using `chunkSize = 2` would for
#' example be equivalent with a parallel processing setup with 2 cores since
#' only 2 can be used at each iteration.
#'
#' This function is ideal for `Spectra` backends that don't allow parallel
#' retrieval of peaks data (such as backends using SQL database backends).
#'
#' @param x `MsExperiment`.
#'
#' @param msLevel `integer(1)` with the MS level
#'
#' @param param defining the peak finding algorithm
#'
#' @param chunkSize `integer(1)` with the number of samples for which the peaks
#'     data should be loaded in memory.
#'
#' @param BPPARAM parallel processing setting.
#'
#' @author Johannes Rainer
#'
#' @noRd
.mse_find_chrom_peaks_chunks <- function(x, msLevel = 1L, param, ...,
                                         chunkSize = 1L,
                                         BPPARAM = bpparam()) {
    res <- unlist(
        .mse_spectrapply_chunks(x, FUN = .mse_find_chrom_peaks_chunk,
                                msLevel = msLevel, param = param,
                                BPPARAM = BPPARAM),
        recursive = FALSE, use.names = FALSE)
    sidx <- vapply(res, function(z) if (is.matrix(z)) nrow(z) else length(z),
                   integer(1), USE.NAMES = FALSE)
    res <- cbind(do.call(rbind, res), sample = rep(seq_along(x), sidx))
    if (!length(res))
        .empty_chrom_peaks()
    else res
}

## .mse_find_chrom_peaks_chunks <- function(x, msLevel = 1L, param, ...,
##                                          chunkSize = 1L,
##                                          BPPARAM = bpparam()) {
##     if (!length(spectra(x)))
##         stop("No spectra available.")
##     if (!any(names(x@sampleDataLinks) == "spectra"))
##         stop("No link between samples and spectra found. Please use ",
##              "'linkSampleData' to define which spectra belong to which ",
##              "samples.")
##     idx <- seq_along(x)
##     chunks <- split(idx, ceiling(idx / chunkSize))
##     sps <- spectra(x)[x@sampleDataLinks[["spectra"]][, 2L]]
##     sps$.SAMPLE_IDX <- x@sampleDataLinks[["spectra"]][, 1L] # or as.factor?
##     res <- unlist(lapply(chunks, function(z) {
##         .mse_find_chrom_peaks_chunk(
##             sps[sps$.SAMPLE_IDX %in% z],
##             msLevel = msLevel, param = param,
##             BPPARAM = BPPARAM)
##     }), recursive = FALSE, use.names = FALSE)
##     sidx <- vapply(res, function(z) if (is.matrix(z)) nrow(z) else length(z),
##                    integer(1), USE.NAMES = FALSE)
##     res <- cbind(do.call(rbind, res), sample = rep(seq_along(x), sidx))
##     if (!length(res))
##         .empty_chrom_peaks()
##     else res
## }

#' Perform the peak detection on a *chunk* of samples. Data realization (i.e.
#' loading of the peak data) is performed without parallel processing while
#' peak detection is performed in parallel.
#'
#' @param x `Spectra` representing a chunk of samples. Needs a spectra variable
#'     `.SAMPLE_IDX` that defined the sample to which the spectra belong to.
#'
#' @author Johannes Rainer
#'
#' @noRd
.mse_find_chrom_peaks_chunk <- function(x, msLevel = 1L, param, ...,
                                        BPPARAM = bpparam()) {
    sidx <- unique(x$.SAMPLE_IDX)
    x <- filterMsLevel(x, msLevel = msLevel)
    if (length(x))
        f <- factor(x$.SAMPLE_IDX, levels = sidx)
    else f <- factor(integer(), levels = sidx)
    bpmapply(
        split(Spectra::peaksData(x, columns = c("mz", "intensity"),
                                 BPPARAM = SerialParam()), f),
        split(rtime(x), f),
        FUN = function(p, rt, prm, msl) {
            vals_per_spect <- vapply(p, nrow, integer(1), USE.NAMES = FALSE)
            p <- do.call(rbind, p)
            if (!length(p))
                return(NULL)            # not returning matrix because of rbind
            if (inherits(prm, "CentWaveParam")) {
                centroided <- all(centroided(x))
                if (is.na(centroided)) {
                    centroided <- isCentroided(x[ceiling(length(x) / 3)])
                    if (is.na(centroided) || !centroided)
                        warning("Your data appears to be not centroided! ",
                                "CentWave works best on data in centroid mode.")
                }
            }
            if (is.unsorted(rt))
                stop("Spectra are not ordered by retention time", .call = FALSE)
            do.call(
                .param_to_fun(prm),
                args = c(list(mz = p[, 1L], int = p[, 2L], scantime = rt,
                              valsPerSpect = vals_per_spect), as(prm, "list")))

        }, MoreArgs = list(prm = param, msl = msLevel), SIMPLIFY = FALSE,
        USE.NAMES = FALSE, BPPARAM = BPPARAM)
}

#' generic method to apply a filtering to the spectra data. The function will
#' apply the filtering and (most importantly) keep/update the link between
#' spectra and samples.
#'
#' @importMethodsFrom Spectra selectSpectraVariables
#'
#' @param x `MsExperiment`.
#'
#' @param FUN filter function.
#'
#' @param ... parameters for `FUN`.
#'
#' @author Johannes Rainer
#'
#' @noRd
.mse_filter_spectra <- function(x, FUN, ...) {
    ls <- length(spectra(x))
    have_links <- length(x@sampleDataLinks[["spectra"]]) > 0
    if (have_links)
        x@spectra$._SAMPLE_IDX <- seq_len(ls)
    x@spectra <- FUN(x@spectra, ...)
    if (have_links) {
        if (ls != length(spectra(x))) {
            sdl <- x@sampleDataLinks[["spectra"]]
            idx <- match(sdl[, 2L], x@spectra$._SAMPLE_IDX)
            keep <- !is.na(idx)
            sdl <- sdl[keep, , drop = FALSE]
            sdl[, 2L] <- idx[keep]
            x@sampleDataLinks[["spectra"]] <- sdl
        }
        svs <- spectraVariables(spectra(x))
        x@spectra <- selectSpectraVariables(
            x@spectra, svs[svs != "._SAMPLE_IDX"])
    }
    x
}

#' Ensure that each spectrum is assigned to a sample and that we only have 1:1
#' mappings. That is important for most code involving splitting of samples
#' etc.
#'
#' @noRd
.mse_check_spectra_sample_mapping <- function(x) {
    if (!length(x@sampleDataLinks[["spectra"]]))
        stop("No links between samples and spectra are present.", call. = FALSE)
    if (nrow(x@sampleDataLinks[["spectra"]]) != length(x@spectra))
        stop("This functionality requires that all spectra in 'object' are ",
             "assigned to a sample.", call. = FALSE)
    if (anyDuplicated(x@sampleDataLinks[["spectra"]][, 2L]))
        stop("This functionality requires that each spectrum is only ",
             "assigned to a single sample.", call. = FALSE)
    NULL
}
