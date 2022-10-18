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

test_that(".mse_filter_spectra works", {
    ## Create a custom, small test object.
    fls <- normalizePath(faahko_3_files)
    df <- data.frame(mzML_file = basename(fls),
                     dataOrigin = fls,
                     sample = c("ko15", "ko16", "ko18"))
    a <- Spectra::Spectra(fls[1])
    b <- Spectra::Spectra(fls[2])
    c <- Spectra::Spectra(fls[3])

    ## Select first 10 spectra from a and c and last 10 from b and combine them
    sps <- c(a[1:10], b[(length(b)-9):length(b)], c[1:10])

    tst <- MsExperiment()
    spectra(tst) <- sps
    sampleData(tst) <- DataFrame(df)
    ## Link samples to spectra.
    tst <- linkSampleData(tst, with = "sampleData.dataOrigin = spectra.dataOrigin")

    res <- .mse_filter_spectra(tst, filterRt, rt = c(2502, 2505))
    expect_true(length(spectra(res)) == 4L)
    expect_equal(res@sampleDataLinks[["spectra"]],
                 cbind(c(1L, 1L, 3L, 3L), 1:4))

    ## Some artificial sample assignment.
    tst@sampleDataLinks[["spectra"]] <- cbind(
        c(1, 1, 1, 2, 2, 2, 3, 3, 3),
        c(2, 3, 4, 2, 3, 4, 2, 3, 4))
    res <- .mse_filter_spectra(tst, filterRt, rt = c(2502, 2505))
    ## Filtering will filter based on rt
    expect_true(length(spectra(res)) == 4L)
    ## Sample assignment should map the first 2 spectra to all 3 samples
    expect_equal(res@sampleDataLinks[["spectra"]],
                 cbind(c(1L, 1L, 2L, 2L, 3L, 3L), c(1L, 2L, 1L, 2L, 1L, 2L)))
})

test_that("filterRt,MsExperiment works", {
    res <- filterRt(mse, rt = c(2700, 2900))
    expect_true(all(rtime(spectra(res)) > 2700 & rtime(spectra(res)) < 2900))
    b <- spectra(mse[2])
    B <- spectra(res[2])
    expect_equal(rtime(B), rtime(filterRt(b, rt = c(2700, 2900))))

    res <- filterRt(mse, rt = c(2700, 2900), msLevel = 2L)
    expect_equal(rtime(spectra(res)), rtime(spectra(mse)))
})

test_that("filterFile,XcmsExperiment works", {
    res <- filterFile(mse)
    expect_s4_class(res, "MsExperiment")
    expect_true(length(res) == 0)

    res <- filterFile(mse, 2)
    expect_equal(res, mse[2])
    res <- filterFile(mse, c(3, 1))
    expect_equal(res, mse[c(1, 3)])
})
