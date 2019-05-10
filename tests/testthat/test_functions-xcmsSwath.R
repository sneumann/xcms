## This should be added to the msdata package!
fl <- "/Users/jo/Projects/git/michaelwitting/metabolomics2018/data/PestMix1_SWATH.mzML"
if (file.exists(fl)) {
    swth <- readMSData(fl, mode = "onDisk")
    cwp <- CentWaveParam(snthresh = 5, noise = 100, ppm = 10,
                         peakwidth = c(3, 20))
    swth <- findChromPeaks(swth, param = cwp)
    swth <- findChromPeaksIsolationWindow(swth, param = cwp)
} else swth <- NULL


test_that(".which_mz_in_range works", {
    mz <- 3.4
    lowerMz <- c(1, 3, 8, 12)
    upperMz <- c(3, 7, 11, 15)
    res <- .which_mz_in_range(mz, lowerMz, upperMz)
    expect_equal(res, 2L)
    res <- .which_mz_in_range(c(3, 3.4, 9), lowerMz, upperMz)
    expect_equal(res, list(c(1L, 2L), 2L, 3L))
})

test_that(".which_chrom_peak_overlap_rt works", {
    pks <- cbind(rtmin = c(1, 2, 3, 4, 5), rtmax = c(2, 3, 4, 5, 6))
    res <- .which_chrom_peak_overlap_rt(c(rtmin = 1.1, rtmax = 2.2), pks)
    expect_equal(res, c(1L, 2L))
    res <- .which_chrom_peak_overlap_rt(c(rtmin = 3.1, rtmax = 3.2), pks)
    expect_equal(res, 3L)
})

test_that(".which_chrom_peak_diff_rt works", {
    pks <- cbind(rt = c(1, 2, 3, 4, 5), rtmax = c(2, 3, 4, 5, 6))
    res <- .which_chrom_peak_diff_rt(c(rt = 1.1), pks, diffRt = 2)
    expect_equal(res, c(1L, 2L, 3L))
    res <- .which_chrom_peak_diff_rt(c(rt = 3.1), pks, diffRt = 0.101)
    expect_equal(res, 3L)
})

test_that(".reconstruct_ms2_for_chrom_peak works", {
    if (!is.null(swth)) {
        # Fluopicolide, exact mass = 381.965430576, [M+H]+ = 382.972706
        pk <- chromPeaks(swth, mz = 382.972706, ppm = 10)
        res <- xcms:::.reconstruct_ms2_for_chrom_peak(pk, swth, fromFile = 7L,
                                                      expandRt = 3, diffRt = 2,
                                                      minCor = 0.8)
        expect_true(is(res, "Spectra"))
        expect_equal(unname(fromFile(res)), 7L)
        expect_equal(length(mcols(res)$ms2_peak_id[[1]]), 14)
        ## chr1 <- chromatogram(swth, mz = pk[, c("mzmin", "mzmax")],
        ##                      rt = pk[, c("rtmin", "rtmax")])
        ## pks <- chromPeaks(swth)[mcols(res)$ms2_peak_id[[1]], ]
        ## chr2 <- chromatogram(filterIsolationWindow(swth, mz = pk[, "mz"]),
        ##                      mz = pks[, c("mzmin", "mzmax")],
        ##                      rt = pks[, c("rtmin", "rtmax")], msLevel = 2L)
        ## plot(chr1, ylim = c(0, 2500))
        ## for (i in 1:nrow(chr2)) {
        ##     points(rtime(chr2[i, 1]), chr2[i, 1]@intensity, col = "#0000ff80",
        ##            type = "l")
        ## }
        res <- xcms:::.reconstruct_ms2_for_chrom_peak(pk, swth, fromFile = 7L,
                                                      expandRt = 3, diffRt = 2,
                                                      minCor = 1)
        expect_true(is(res, "Spectra"))
        expect_equal(unname(fromFile(res)), 7L)
        expect_equal(length(mcols(res)$ms2_peak_id[[1]]), 0)
    }
})

test_that(".reconstruct_ms2_for_peaks_file works", {
    ## No MS1 level peaks: expect empty.
})

test_that("reconstructChromPeakSpectra works", {
})
