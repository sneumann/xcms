library(MsExperiment)
mse <- MsExperiment()
fls <- normalizePath(faahko_3_files)
df <- data.frame(mzML_file = basename(fls),
                 dataOrigin = fls,
                 sample = c("ko15", "ko16", "ko18"))

spectra(mse) <- Spectra::Spectra(fls)
sampleData(mse) <- DataFrame(df)
## Link samples to spectra.
mse <- linkSampleData(mse, with = "sampleData.dataOrigin = spectra.dataOrigin")

test_that(".param_to_fun works", {
    expect_equal(.param_to_fun(CentWaveParam()), "do_findChromPeaks_centWave")
    expect_error(.param_to_fun(DataFrame()), "No peak detection")
})

test_that(".mse_sample_apply works", {
    dummy <- MsExperiment()
    spectra(dummy) <- spectra(mse)

    res <- .mse_sample_apply(dummy, length, BPPARAM = SerialParam())
    expect_true(length(res) == 0)

    res <- .mse_sample_apply(mse, length, BPPARAM = SerialParam())
    expect_equal(res, list(`1` = 1L, `2` = 1L, `3` = 1L))

    res <- .mse_sample_apply(mse, function(z, msLevel) {
        length(filterMsLevel(spectra(z), msLevel))
    }, msLevel = 1L, BPPARAM = SerialParam())
    expect_equal(unlist(res, use.names = FALSE),
                 as.integer(table(fromFile(od_x))))
    res <- .mse_sample_apply(mse, function(z, msLevel) {
        length(filterMsLevel(spectra(z), msLevel))
    }, msLevel = 2L, BPPARAM = SerialParam())
    expect_equal(unlist(res), c(`1` = 0L, `2` = 0L, `3` = 0L))
})

test_that(".mse_sample_spectra_apply works", {
    dummy <- MsExperiment()
    expect_error(.mse_sample_spectra_apply(dummy), "No spectra")

    spectra(dummy) <- spectra(mse)
    expect_error(.mse_sample_spectra_apply(dummy, length), "No link")

    res <- .mse_sample_spectra_apply(mse, length, BPPARAM = SerialParam())
    expect_equal(res, list(`1` = 1278, `2` = 1278, `3` = 1278))

    res <- .mse_sample_spectra_apply(mse, function(z, msLevel) {
        length(filterMsLevel(z, msLevel))
    }, msLevel = 1L, BPPARAM = SerialParam())
    expect_equal(unlist(res, use.names = FALSE),
                 as.integer(table(fromFile(od_x))))
    res <- .mse_sample_spectra_apply(mse, function(z, msLevel) {
        length(filterMsLevel(z, msLevel))
    }, msLevel = 2L, BPPARAM = SerialParam())
    expect_equal(unlist(res), c(`1` = 0L, `2` = 0L, `3` = 0L))
})

test_that(".mse_find_chrom_peaks_sample works", {
    p <- CentWaveParam(noise = 10000, snthresh = 40, prefilter = c(3, 10000))
    res <- .mse_find_chrom_peaks_sample(spectra(mse[2L]), param = p)
    tmp <- chromPeaks(faahko_xod)
    tmp <- tmp[tmp[, "sample"] == 2, colnames(tmp) != "sample"]
    rownames(tmp) <- NULL
    expect_equal(res, tmp)

    res <- .mse_find_chrom_peaks_sample(spectra(mse[1L]), param = p,
                                        msLevel = 2L)
    expect_true(is.null(res))
})

test_that(".mse_find_chrom_peaks works", {
    p <- CentWaveParam(noise = 10000, snthresh = 40, prefilter = c(3, 10000))
    res <- .mse_find_chrom_peaks(mse, param = p)
    tmp <- chromPeaks(faahko_xod)
    rownames(tmp) <- NULL
    expect_equal(tmp, res)

    res <- .mse_find_chrom_peaks(mse, param = p, msLevel = 2L)
    expect_true(nrow(res) == 0)
    expect_equal(colnames(res), colnames(.empty_chrom_peaks()))
})

test_that(".mse_find_chrom_peaks_chunks works", {
    p <- CentWaveParam(noise = 10000, snthresh = 40, prefilter = c(3, 10000))
    res <- .mse_find_chrom_peaks_chunks(mse, param = p)
    tmp <- chromPeaks(faahko_xod)
    rownames(tmp) <- NULL
    expect_equal(tmp, res)

    res <- .mse_find_chrom_peaks_chunks(mse, param = p, msLevel = 2L)
    expect_true(nrow(res) == 0)
    expect_equal(res, .empty_chrom_peaks())
})

test_that(".mse_find_chrom_peaks_chunk works", {
    p <- CentWaveParam(noise = 10000, snthresh = 40, prefilter = c(3, 10000))
    sps <- spectra(mse[1:2])[mse[1:2]@sampleDataLinks[["spectra"]][, 2L]]
    sps$.SAMPLE_IDX <- mse[1:2]@sampleDataLinks[["spectra"]][, 1L]

    res <- .mse_find_chrom_peaks_chunk(sps, param = p)
    expect_true(is.list(res))
    expect_true(length(res) == 2)

    cp <- chromPeaks(faahko_xod)
    f <- cp[, "sample"]
    rownames(cp) <- NULL
    cpl <- split.data.frame(cp[, colnames(cp) != "sample"], f)
    expect_equal(cpl[[1L]], res[[1L]])
    expect_equal(cpl[[2L]], res[[2L]])

    res <- .mse_find_chrom_peaks_chunk(sps, param = p, msLevel = 2L)
    expect_true(is.list(res))
    expect_true(length(res) == 2)
    expect_true(is.null(res[[1L]]))
    expect_true(is.null(res[[2L]]))
})
