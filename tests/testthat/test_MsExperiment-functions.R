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

test_that(".mse_spectrapply_chunks works", {
    expect_error(.mse_spectrapply_chunks(MsExperiment), "spectra")

    myident <- function(z, ...) {z}
    res <- .mse_spectrapply_chunks(mse, FUN = myident)
    expect_true(is.list(res))
    expect_true(length(res) == 3)
    expect_equal(rtime(res[[1L]]), rtime(spectra(mse[1L])))
    expect_equal(rtime(res[[2L]]), rtime(spectra(mse[2L])))
    expect_equal(rtime(res[[3L]]), rtime(spectra(mse[3L])))

    res <- .mse_spectrapply_chunks(mse, FUN = myident, chunkSize = 2)
    expect_true(is.list(res))
    expect_true(length(res) == 2)
    expect_equal(rtime(res[[1L]]), c(rtime(spectra(mse[1L])),
                                     rtime(spectra(mse[2L]))))
    expect_equal(rtime(res[[2L]]), rtime(spectra(mse[3L])))
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

test_that(".mse_check_spectra_sample_mapping works", {
    expect_true(length(.mse_check_spectra_sample_mapping(mse)) == 0)

    tmp <- mse
    tmp@sampleDataLinks[["spectra"]] <-
        mse@sampleDataLinks[["spectra"]][1:100, ]
    expect_error(.mse_check_spectra_sample_mapping(tmp), "assigned to a sample")

    tmp@sampleDataLinks[["spectra"]] <- mse@sampleDataLinks[["spectra"]]
    tmp@sampleDataLinks[["spectra"]][3, ] <- c(2L, 2L)
    expect_error(.mse_check_spectra_sample_mapping(tmp), "single sample")
})

test_that(".mse_profmat_chunk works", {
    tmp <- mse[1]
    ref <- profMat(faahko_od, fileIndex = 1)

    sps <- spectra(tmp)
    sps$.SAMPLE_IDX <- 1L
    res <- .mse_profmat_chunk(sps)
    expect_equal(unname(res), ref)

    ref <- profMat(faahko_od, fileIndex = 1:2, step = 2, returnBreaks = TRUE)
    tmp <- mse[1:2]
    sps <- spectra(tmp)
    sps$.SAMPLE_IDX <- tmp@sampleDataLinks[["spectra"]][, 1L]
    res <- .mse_profmat_chunk(sps, step = 2, returnBreaks = TRUE)
    expect_equal(unname(res), ref)
    expect_true(all(names(res[[1L]]) == c("profMat", "breaks")))

    res <- .mse_profmat_chunk(sps, step = 2, msLevel = 2)
    expect_equal(length(res), 2)
    expect_true(nrow(res[[1L]]) == 0)
    expect_true(nrow(res[[2L]]) == 0)
})

test_that(".mse_profmat_chunks works", {
    expect_error(.mse_profmat_chunks(mse, fileIndex = 5), "bounds")
    expect_error(.mse_profmat_chunks(mse, fileIndex = 1:5), "bounds")

    ref <- profMat(faahko_od, fileIndex = 3)
    res <- .mse_profmat_chunks(mse, fileIndex = 3)
    expect_equal(ref, res)

    ref <- profMat(faahko_od, returnBreaks = TRUE, step = 4)
    res <- .mse_profmat_chunks(mse, chunkSize = 2L, step = 4,
                               returnBreaks = TRUE)
    expect_equal(res, ref)
    expect_equal(names(res[[1L]]), c("profMat", "breaks"))

    res <- .mse_profmat_chunks(mse, chunkSize = 3L, msLevel = 2L)
    expect_true(length(res) == 3)

    ## Testing the method.
    res <- profMat(mse, chunkSize = 2L, step = 4, returnBreaks = TRUE)
    expect_equal(res, ref)
})

test_that(".obiwarp_spectra works", {
    p <- ObiwarpParam(binSize = 3.4, centerSample = 1L)
    ref <- split(adjustRtime(faahko_od, param = p), fromFile(faahko_od))

    a <- spectra(mse[1L])
    b <- spectra(mse[2L])
    res <- xcms:::.obiwarp_spectra(b, a, param = p)
    expect_equal(res, unname(ref[[2L]]))
    res <- xcms:::.obiwarp_spectra(spectra(mse[3L]), a, param = p)
    expect_equal(res, unname(ref[[3L]]))

    ## Test with different MS levels.
    res_sub <- xcms:::.obiwarp_spectra(b[-seq(1, length(b), by = 3)],
                                   a[-seq(1, length(a), by = 3)],
                                   param = p)
    a$msLevel[seq(1, length(a), by = 3)] <- 2L
    b$msLevel[seq(1, length(b), by = 3)] <- 2L
    res2 <- xcms:::.obiwarp_spectra(b, a, param = p)
    expect_equal(length(res2), length(b))
    expect_equal(res_sub, res2[-seq(1, length(b), by = 3)])
    expect_true(cor(res, res2) > 0.999)
})

test_that(".mse_obiwarp_chunks works", {
    p <- ObiwarpParam(binSize = 50, centerSample = 1L)
    ref <- split(adjustRtime(faahko_od, param = p), fromFile(faahko_od))

    expect_error(.mse_obiwarp_chunks(mse, ObiwarpParam(centerSample = 6)),
                 "integer between 1 and 3")

    res <- xcms:::.mse_obiwarp_chunks(mse, p)
    expect_true(is.list(res))
    expect_equal(length(res), length(mse))
    expect_equal(unname(ref[[1L]]), res[[1L]])
    expect_equal(unname(ref[[2L]]), res[[2L]])
    expect_equal(unname(ref[[3L]]), res[[3L]])

    ## Subset alignment...
    p <- ObiwarpParam(binSize = 30, centerSample = 1L, subset = c(1, 3))
    ref <- split(adjustRtime(faahko_od, param = p), fromFile(faahko_od))

    res <- xcms:::.mse_obiwarp_chunks(mse, p, chunkSize = 2L)
    expect_equal(unname(ref[[1L]]), res[[1L]])
    expect_equal(unname(ref[[2L]]), res[[2L]])
    expect_equal(unname(ref[[3L]]), res[[3L]])

    expect_error(.mse_obiwarp_chunks(mse, p, msLevel = 2), "MS level")
})

test_that("readMsExperiment works", {
    expect_error(a <- readMsExperiment(), "'files'")
    expect_error(a <- readMsExperiment("a"), "not found")

    a <- readMsExperiment(faahko_3_files[1:2])
    expect_s4_class(a, "MsExperiment")
    expect_true(length(a) == 2)
    expect_true(nrow(sampleData(a)) == 2)

    df <- data.frame(sidx = 1:3, other_ann = c("a", "b", "c"))
    expect_error(readMsExperiment(faahko_3_files[1:2], df), "of files")
    a <- readMsExperiment(faahko_3_files[1:2], df[1:2, ])
    expect_s4_class(a, "MsExperiment")
    expect_true(length(a) == 2)
    expect_true(nrow(sampleData(a)) == 2)
    expect_equal(sampleData(a)$other_ann, c("a", "b"))
})
