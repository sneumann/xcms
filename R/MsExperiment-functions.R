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
                              f = factor(), BPPARAM = SerialParam())
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
                                    progressbar = TRUE, BPPARAM = bpparam()) {
    if (!length(spectra(x)))
        stop("No spectra available.")
    if (!any(names(x@sampleDataLinks) == "spectra"))
        stop("No link between samples and spectra found. Please use ",
             "'linkSampleData' to define which spectra belong to which ",
             "samples.")
    idx <- seq_along(x)
    chunks <- split(idx, ceiling(idx / chunkSize))
    if (progressbar) {
        pb <- progress_bar$new(format = paste0("[:bar] :current/:",
                                               "total (:percent) in ",
                                               ":elapsed"),
                               total = length(chunks),
                               clear = FALSE, show_after = 0)
        pb$tick(0)
    }
    sps <- spectra(x)[x@sampleDataLinks[["spectra"]][, 2L]]
    sps$.SAMPLE_IDX <- x@sampleDataLinks[["spectra"]][, 1L] # or as.factor?
    lapply(chunks, function(z, ..., pb) {
        suppressMessages(
            res <- FUN(sps[sps$.SAMPLE_IDX %in% z], ...)
        )
        if (progressbar) pb$tick()
        res
    }, ..., pb = pb, BPPARAM = BPPARAM)
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
                                chunkSize = chunkSize, BPPARAM = BPPARAM),
        recursive = FALSE, use.names = FALSE)
    sidx <- vapply(res, function(z) if (is.matrix(z)) nrow(z) else length(z),
                   integer(1), USE.NAMES = FALSE)
    res <- cbind(do.call(rbind, res), sample = rep(seq_along(x), sidx))
    if (!length(res))
        .empty_chrom_peaks()
    else res
}

#' Perform the peak detection on a *chunk* of samples. Data realization (i.e.
#' loading of the peak data) is performed without parallel processing while
#' peak detection is performed in parallel. Note: we could even perform the
#' data realization in parallel using the supported BPPARAM (i.e. using the
#' `backendBpparam` function). It would however not help much since
#' parallelization if performed on `dataStorage`.
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
    lx <- length(x)
    if (lx)
        f <- factor(x$.SAMPLE_IDX, levels = sidx)
    else f <- factor(integer(), levels = sidx)
    ## Check for random number of spectra if they are centroided. NOT all.
    if (inherits(param, "CentWaveParam")) {
        cntr <- all(centroided(x[sort(sample(seq_along(x), min(c(100, lx))))]))
        if (is.na(cntr)) {
            cntr <- isCentroided(x[ceiling(lx / 3)])
            if (is.na(cntr) || !cntr)
                warning("Your data appears to be not centroided! ",
                        "CentWave works best on data in centroid mode.")
        }
    }
    bpmapply(
        split(peaksData(x, columns = c("mz", "intensity"), f = factor(),
                        BPPARAM = SerialParam()), f),
        split(rtime(x), f),
        FUN = function(p, rt, prm, msl) {
            vals_per_spect <- vapply(p, nrow, integer(1), USE.NAMES = FALSE)
            p <- do.call(rbind, p)
            if (!length(p))
                return(NULL)            # not returning matrix because of rbind
            if (is.unsorted(rt))
                stop("Spectra are not ordered by retention time", .call = FALSE)
            do.call(
                .param_to_fun(prm),
                args = c(list(mz = p[, 1L], int = p[, 2L], scantime = rt,
                              valsPerSpect = vals_per_spect), as(prm, "list")))
        }, MoreArgs = list(prm = param, msl = msLevel), SIMPLIFY = FALSE,
        USE.NAMES = FALSE, BPPARAM = BPPARAM)
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

#' Create a `profMat` for spectra of each sample in the chunk of spectra.
#'
#' @param ... can also include `returnBreaks = TRUE`.
#'
#' @noRd
.mse_profmat_chunk <- function(x, method = "bin", step = 0.1,
                               baselevel = NULL, basespace = NULL,
                               mzrange. = NULL, msLevel = 1L, ...,
                               BPPARAM = bpparam()) {
    sidx <- unique(x$.SAMPLE_IDX)
    x <- filterMsLevel(x, msLevel = msLevel)
    if (length(x))
        f <- factor(x$.SAMPLE_IDX, levels = sidx)
    else f <- factor(integer(), levels = sidx)
    bplapply(
        split(Spectra::peaksData(x, columns = c("mz", "intensity"),
                                 f = factor(),
                                 BPPARAM = SerialParam()), f),
        FUN = .peaksdata_profmat, method = method, step = step,
        baselevel = baselevel, basespace = basespace, mzrange. = mzrange.,
        ..., BPPARAM = BPPARAM)
}

#' Calculate a profile matrix from a `Spectra` (should be from a single
#' file/sample and from a single MS level).
#'
#' @noRd
.peaksdata_profmat <- function(x, method = "bin", step = 0.1, baselevel = NULL,
                             basespace = NULL, mzrange. = NULL, ...) {
    pk_count <- vapply(x, nrow, integer(1), USE.NAMES = FALSE)
    empty_spectra <- which(!pk_count)
    if (length(empty_spectra))
        pk_count <- pk_count[-empty_spectra]
    x <- do.call(rbind, x)
    if (length(x)) {
        res <- .createProfileMatrix(mz = x[, 1], int = x[, 2],
                                    valsPerSpect = pk_count,
                                    method = method, step = step,
                                    baselevel = baselevel,
                                    basespace = basespace,
                                    mzrange. = mzrange., ...)
        if (length(empty_spectra)) {
            if (any(names(res) == "profMat"))
                res$profMat <- .insertColumn(res$profMat,
                                             empty_spectra, 0)
            else
                res <- .insertColumn(res, empty_spectra, 0)
        }
        res
    } else matrix(numeric(), nrow = 0, ncol = 0)
}

.mse_profmat_chunks <- function(x, msLevel = 1L, method = "bin", step = 0.1,
                               baselevel = NULL, basespace = NULL,
                               mzrange. = NULL, fileIndex = seq_along(x),
                               chunkSize = 1L, ..., BPPARAM = bpparam()) {
    if (!all(fileIndex %in% seq_along(x)))
        stop("fileIndex out of bounds", call. = FALSE)
    unlist(
        .mse_spectrapply_chunks(x[fileIndex], FUN = .mse_profmat_chunk,
                                method = method, step = step,
                                baselevel = baselevel, basespace = basespace,
                                mzrange. = mzrange., msLevel = msLevel, ...,
                                chunkSize = chunkSize, BPPARAM = BPPARAM),
        recursive = FALSE, use.names = FALSE)
}

.mse_obiwarp_chunks <- function(x, param, msLevel = 1L, chunkSize = 1L,
                                BPPARAM = bpparam()) {
    message("value ", param@rtimeDifferenceThreshold)
    rt_raw <- split(rtime(x), fromFile(x))
    subset_idx <- subset(param)
    if (length(subset_idx))
        x <- x[subset_idx]

    if (!length(centerSample(param)))
        centerSample(param) <- floor(median(seq_along(x)))
    ref_idx <- centerSample(param)
    if (!(ref_idx %in% seq_along(x)))
        stop("'centerSample' needs to be an integer between 1 and ", length(x))
    ref_sps <- filterMsLevel(spectra(x[ref_idx]), msLevel = msLevel)
    ref_pm <- .peaksdata_profmat(peaksData(ref_sps, f = factor()),
                                 method = "bin", step = binSize(param),
                                 returnBreaks = TRUE)
    res <- unlist(.mse_spectrapply_chunks(
        x, FUN = function(z, ref, ref_pm, param, msLevel, BPPARAM) {
            z <- setBackend(
                selectSpectraVariables(z, c("rtime", "msLevel", ".SAMPLE_IDX",
                                            "dataStorage", "scanIndex",
                                            "mz", "intensity")),
                MsBackendMemory(), BPPARAM = SerialParam())
            bplapply(split(z, f = as.factor(z$.SAMPLE_IDX)),
                     FUN = .obiwarp_spectra, ref = ref, ref_pm = ref_pm,
                     param = param, msLevel = msLevel, BPPARAM = BPPARAM)
        }, ref = ref_sps, ref_pm = ref_pm, param = param, msLevel = msLevel,
        chunkSize = chunkSize, BPPARAM = BPPARAM),
        recursive = FALSE, use.names = FALSE)

    if (length(subset_idx)) {
        res[[ref_idx]] <- rt_raw[subset_idx][[ref_idx]]
        rt_adj <- vector("list", length(rt_raw))
        rt_adj[subset_idx] <- res
        res <- adjustRtimeSubset(rt_raw, rt_adj, subset = subset_idx,
                                 method = subsetAdjust(param))
    } else
        res[[ref_idx]] <- rt_raw[[ref_idx]]
    res
}

#' Performs alignment of other against ref and returns the adjusted retention
#' times (for ALL spectra, even in other MS levels). other and ref are expected
#' to be `Spectra` objects.
#'
#' @noRd
.obiwarp_spectra <- function(other, ref, ref_pm = list(), param, msLevel = 1L) {
    ## why is that not idea? well, we can't do that in chunks!
    rt_raw <- rtime(other)
    n_all <- length(rt_raw)
    rt_ms <- which(msLevel(other) == msLevel)
    ref <- filterMsLevel(ref, msLevel = msLevel)
    other <- filterMsLevel(other, msLevel = msLevel)
    if (!(length(ref) & length(other)))
        stop("No spectra with MS level ", msLevel, " present")
    if (!length(ref_pm))
        ref_pm <- .peaksdata_profmat(peaksData(ref, f = factor()),
                                     method = "bin",
                                     step = binSize(param),
                                     returnBreaks = TRUE)
    other_pm <- .peaksdata_profmat(peaksData(other, f = factor()),
                                   method = "bin",
                                   step = binSize(param),
                                   returnBreaks = TRUE)
    adj <- .obiwarp_bare(rtime(ref), rtime(other), ref_pr = ref_pm,
                         other_pr = other_pm, param = param)
    n_adj <- length(adj)
    if (length(n_all) != n_adj) {
        ## Have to adjust rts for MS levels other than msLevel
        adj_fun <- approxfun(x = rtime(other), y = adj)
        rt_adj <- adj_fun(rt_raw)
        tmp_rt <- rtime(other[1L])
        idx_below <- which(rt_raw < tmp_rt)
        if (length(idx_below))
            rt_adj[idx_below] <- rt_raw[idx_below] + adj[1L] - tmp_rt
        tmp_rt <- rtime(other[n_adj])
        idx_above <- which(rt_raw > tmp_rt)
        if (length(idx_above))
            rt_adj[idx_above] <- rt_raw[idx_above] + adj[n_adj] - tmp_rt
        rt_adj
    } else adj
}

#' This function extracts a chromatogram for the provided rt, m/z ranges, MS
#' level and isolationWindow from the `MsExperiment`. Parameter
#' `isolationWindow` ensures that, for MS2 spectra, not simply all MS2 spectra
#' are used for the chromatogram, but only those with matching
#' `isolationWindowTargetMz` (and hence the same set of ions). Chromatograms
#' for MS1 will not need `isolationWindow`.
#'
#' @importFrom MsExperiment sampleData
#'
#' @param msLevel MS level. Can be of length 1 or equal to the number of rows
#'     of rt.
#'
#' @param isolationWindow additional variable eventually subsetting spectra,
#'     e.g. for SWATH data or similar to ensure chromatograms for MS level 2
#'     data are extracted from MS2 spectra of the correct (same) isolation
#'     window. Set to `NULL` (the default) if there are not isolation windows
#'     defined. Has to be the same lengths than there are chromatograms to be
#'     extracted.
#'
#' @note
#'
#' This function will also pass the `isolationWindowTargetMz` spectra variable
#' to the downstream function to ensure chromatograms from MS2 levels are
#' extracted from the isolation window matching `isolationWindow`. If one of
#' the two variables is not provided (or just `NA`) no chromatogram for MS2
#' data will be extracted!
#'
#' @noRd
.mse_chromatogram <- function(x, rt = matrix(nrow = 0, ncol = 2),
                              mz = matrix(nrow = 0, ncol = 2),
                              aggregationFun = "sum", msLevel = 1L,
                              isolationWindow = NULL,
                              chunkSize = 2L, progressbar = TRUE,
                              BPPARAM = bpparam()) {
    if (!nrow(rt))
        rt <- matrix(c(-Inf, Inf), ncol = 2)
    if (!nrow(mz))
        mz <- matrix(c(-Inf, Inf), ncol = 2)
    if (is.matrix(rt) && ncol(rt) != 2)
        stop("'rt' is expected to be a two-column matrix", call. = FALSE)
    if (is.matrix(mz) && ncol(mz) != 2)
        stop("'mz' is expected to be a two-column matrix", call. = FALSE)
    pks <- cbind(mz, rt)
    npks <- nrow(pks)
    if (length(msLevel) != npks)
        msLevel <- rep(msLevel[1L], npks)
    if (!length(isolationWindow))
        isolationWindow <- rep(1L, npks)
    if (length(isolationWindow) && length(isolationWindow) != npks)
        stop("Length of 'isolationWindow' (if provided) should match the ",
             "number of chromatograms to extract.")
    colnames(pks) <- c("mzmin", "mzmax", "rtmin", "rtmax")
    res <- .mse_spectrapply_chunks(
        x, FUN = function(z, pks, msl, afun, BPPARAM) {
            sidx <- unique(z$.SAMPLE_IDX)
            z <- filterMsLevel(z, msLevel = msLevel)
            rtr <- range(pks[, c("rtmin", "rtmax")], na.rm = TRUE)
            if (all(is.finite(rtr)))
                z <- filterRt(z, rt = rtr)
            lz <- length(z)
            if (lz)
                f <- factor(z$.SAMPLE_IDX, levels = sidx)
            else f <- factor(integer(), levels = sidx)
            bpmapply(
                split(Spectra::peaksData(z, columns = c("mz", "intensity"),
                                         f = factor(),
                                         BPPARAM = SerialParam()), f),
                split(rtime(z), f),
                split(msLevel(z), f),
                sidx,
                split(isolationWindowTargetMz(z), f),
                FUN = .chromatograms_for_peaks,
                MoreArgs = list(pks = pks, pks_msl = msl,
                                pks_tmz = isolationWindow,
                                aggregationFun = afun),
                SIMPLIFY = FALSE, USE.NAMES = FALSE, BPPARAM = BPPARAM)
        }, pks = pks, msl = msLevel, afun = aggregationFun,
        chunkSize = chunkSize, progressbar = progressbar, BPPARAM = BPPARAM)
    res <- as(do.call(cbind, unlist(res, recursive = FALSE, use.names = FALSE)),
              "MChromatograms")
    fd <- annotatedDataFrameFrom(res, byrow = TRUE)
    fd$mzmin <- mz[, 1]
    fd$mzmax <- mz[, 2]
    fd$rtmin <- rt[, 1]
    fd$rtmax <- rt[, 2]
    res@featureData <- fd
    rownames(res@.Data) <- rownames(fd)
    res@phenoData <- AnnotatedDataFrame(as.data.frame(sampleData(x)))
    colnames(res@.Data) <- rownames(pData(res))
    res
}

#' Split an `MsExperiment` by a spectra variable keeping sample to spectra
#' mapping.
#'
#' @author Johannes Rainer
#'
#' @noRd
.mse_split_spectra_variable <- function(x, f = msLevel(spectra(x))) {
    ls <- length(spectra(x))
    have_links <- length(x@sampleDataLinks[["spectra"]]) > 0
    if (have_links)
        x@spectra$._SPECTRA_IDX <- seq_len(ls)
    spl <- split(x@spectra, f)
    lapply(spl, function(z) {
        tmp <- x
        tmp@spectra <- z
        if (have_links)
            tmp <- .update_sample_data_links_spectra(tmp)
        svs <- unique(c(spectraVariables(z), "mz", "intensity"))
        tmp@spectra <- selectSpectraVariables(z, svs[svs != "._SPECTRA_IDX"])
        tmp
    })
}

#' Update the sampleDataLinks for a subsetted `@spectra` slot within an
#' `MsExperiment`. This function requires the presence of a spectra variable
#' `"._SPECTRA_IDX"` in `@spectra`. Note that this function  **only** updates
#' the `@sampleDataLinks[["spectra"]]` matrix but does **not** update or
#' subset the `@sampleData`.
#'
#' @noRd
.update_sample_data_links_spectra <- function(x) {
    sdl <- x@sampleDataLinks[["spectra"]]
    idx <- match(sdl[, 2L], x@spectra$._SPECTRA_IDX)
    keep <- !is.na(idx)
    sdl <- sdl[keep, , drop = FALSE]
    sdl[, 2L] <- idx[keep]
    x@sampleDataLinks[["spectra"]] <- sdl
    x
}
