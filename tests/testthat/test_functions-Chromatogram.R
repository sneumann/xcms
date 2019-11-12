test_that(".align_chromatogram_approx works", {
    chr1 <- Chromatogram(rtime = c(1, 2, 3, 4, 5, 6, 7, 8),
                         intensity = c(5, 9, 3, 1, 4, 3, 6, 9))
    chr2 <- Chromatogram(rtime = c(3, 4, 6), intensity = c(3, 1, 3))
    res <- .align_chromatogram_approx(chr2, chr1)
    expect_equal(length(chr1), length(res))
    expect_equal(rtime(res), rtime(chr1))
    expect_equal(intensity(res), c(NA, NA, 3, 1, 2, 3, NA, NA))

    res <- .align_chromatogram_approx(chr1, chr2)
    expect_equal(length(chr2), length(res))
    expect_equal(rtime(res), rtime(chr2))
    expect_equal(intensity(res), c(3, 1, 3))

    ## Not perfectly matching rtimes:
    chr1 <- Chromatogram(rtime = c(1.1, 2.1, 3.1, 4.1, 5.1),
                         intensity = c(1, 2, 3, 2, 1))
    chr2 <- Chromatogram(rtime = c(2, 3), intensity = c(3, 5))
    res <- .align_chromatogram_approx(chr2, chr1)
    expect_equal(rtime(res), rtime(chr1))
    expect_equal(intensity(res), c(NA, 3.2, NA, NA, NA))

    res <- .align_chromatogram_approx(chr1, chr2)
    expect_equal(rtime(res), rtime(chr2))
    expect_equal(intensity(res), c(1.9, 2.9))
})

test_that(".align_chromatogram_match_rtime works", {
    chr1 <- Chromatogram(rtime = c(1, 2, 3, 4, 5, 6, 7, 8),
                         intensity = c(5, 9, 3, 1, 4, 3, 6, 9))
    chr2 <- Chromatogram(rtime = c(3, 4, 6), intensity = c(3, 1, 3))
    res <- .align_chromatogram_match_rtime(chr2, chr1)
    expect_equal(length(chr1), length(res))
    expect_equal(rtime(res), rtime(chr1))
    expect_equal(intensity(res), c(NA, NA, 3, 1, NA, 3, NA, NA))

    res <- .align_chromatogram_match_rtime(chr1, chr2)
    expect_equal(rtime(res), rtime(chr2))
    expect_equal(intensity(res), c(3, 4, 6))

    ## Not perfectly matching rtimes:
    chr1 <- Chromatogram(rtime = c(1.1, 2.1, 3.1, 4.1, 5.1),
                         intensity = c(1, 2, 3, 2, 1))
    chr2 <- Chromatogram(rtime = c(2, 3), intensity = c(3, 5))
    res <- .align_chromatogram_match_rtime(chr2, chr1)
    expect_equal(rtime(res), rtime(chr1))
    expect_equal(intensity(res), c(NA, 3, 5, NA, NA))

    res <- .align_chromatogram_match_rtime(chr1, chr2)
    expect_equal(rtime(res), rtime(chr2))
    expect_equal(intensity(res), c(2, 3))
})

test_that(".align_chromatogram works", {
    chr1 <- Chromatogram(rtime = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
                         intensity = c(3, 5, 14, 30, 24, 6, 2, 1, 1, 0))
    chr2 <- Chromatogram(rtime = c(2.5, 3.42, 4.5, 5.43, 6.5),
                         intensity = c(5, 12, 15, 11, 5))
    res <- .align_chromatogram(chr1, chr2)
    expect_equal(res, .align_chromatogram_match_rtime(chr1, chr2))
    res <- .align_chromatogram(chr1, chr2, method = "approx")
    expect_equal(res, .align_chromatogram_approx(chr1, chr2))

    expect_error(.align_chromatogram(chr1, chr2, method = "other"),
                 "should be one of")
})

test_that("align,Chromatogram,Chromatograms works", {
    library(testthat)
    chr1 <- Chromatogram(rtime = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
                         intensity = c(3, 5, 14, 30, 24, 6, 2, 1, 1, 0))
    chr2 <- Chromatogram(rtime = c(2.5, 3.42, 4.5, 5.43, 6.5),
                         intensity = c(5, 12, 15, 11, 5))
    chr3 <- Chromatogram(rtime = c(2.3, 3.2, 4.3, 5.2),
                         intensity = c(8, 9, 19, 8))
    chr4 <- Chromatogram(rtime = c(7, 8, 9, 10),
                         intensity = c(3, 5, 7, 8))

    chrs <- Chromatograms(list(chr2, chr3, chr4))
    res <- align(chrs, chr1)
    expect_equal(dimnames(res), dimnames(chrs))
    expect_equal(res[1, 1], align(chr2, chr1))
    expect_equal(res[2, 1], align(chr3, chr1))
    expect_equal(res[3, 1], align(chr4, chr1))

    chrs <- Chromatograms(list(chr2, chr3, chr4, chr2), nrow = 2)
    res <- align(chrs, chr1, method = "approx")
    expect_equal(dimnames(res), dimnames(chrs))
    expect_equal(res[1, 1], align(chr2, chr1, method = "approx"))
    expect_equal(res[2, 1], align(chr3, chr1, method = "approx"))
    expect_equal(res[1, 2], align(chr4, chr1, method = "approx"))
    expect_equal(res[2, 2], align(chr2, chr1, method = "approx"))
})

test_that(".correlate_chromatogram works", {
    set.seed(112)
    ## create
    rtime1 <- seq(5,25, by = 1)
    intensity1 <- dnorm(rtime1, mean=14, sd=1.0)*200

    rtime2 <- seq(5.1, 25.1, by = 1)
    intensity2 <- dnorm(rtime2, mean=14, sd=1.0)*500

    ## bogus chromatograms
    ch1 <- new("Chromatogram",
               rtime = rtime1,
               intensity = intensity1)

    ch2 <- new("Chromatogram",
               rtime = rtime2,
               intensity = intensity2)

    ## check that correlation with NA values fails
    expect_equal(xcms:::.correlate_chromatogram(ch1, ch2),
                 cor(intensity1, intensity2))
    expect_equal(xcms:::.correlate_chromatogram(ch1, ch2),
                 xcms:::.correlate_chromatogram(ch2, ch1))

    expect_equal(
        xcms:::.correlate_chromatogram(ch1, ch2, align = "approx"),
        xcms:::.correlate_chromatogram(ch2, ch1, align = "approx"),
        tolerance = 0.0001
    )

    ch2 <- filterRt(ch2, rt = c(5.1, 23))
    res <- xcms:::.correlate_chromatogram(ch1, ch2)
    expect_equal(res, cor(intensity(ch2), intensity(ch1)[1:length(ch2)]))

    res <- xcms:::.correlate_chromatogram(ch2, ch1, align = "approx",
                                          use = "everything")
    expect_equal(res, NA_real_)
})

test_that(".chrom_merge_neighboring_peaks works", {
    ints <- c(0.5, 1, 1, 3, 6, 9, 12, 13, 11, 6, 5, 3, 1, 1, 1.5, 1, 4, 6,
              8, 9, 8, 6, 3, 2, 1.3, 1, 0.7, 0.5, 1, 1, 1, 0.5, 3, 5, 8,
              12, 10, 9, 6, 3, 2, 1, 1, 1)
    rts <- 1:length(ints)
    chr <- Chromatogram(rts, ints)
    cwp <- CentWaveParam(snthresh = 0, prefilter = c(1, 1), peakwidth = c(1, 4))
    xchr <- findChromPeaks(chr, param = cwp)
    res <- .chrom_merge_neighboring_peaks(chr, chromPeaks(xchr))
    rownames(res) <- NULL
    expect_equal(res, chromPeaks(xchr))
    res <- .chrom_merge_neighboring_peaks(chr, chromPeaks(xchr), diffRt = 5)
    rownames(res) <- NULL
    expect_equal(res, chromPeaks(xchr))

    ints <- c(0.5, 1, 1, 3, 6, 9, 12, 13, 11, 7, 6, 5.5, 5.2, 5, 5.2, 5.4,
              5.7, 6, 8, 9, 8, 6, 3, 2, 1.3, 1, 0.7, 0.5, 1, 1, 1, 0.5, 3,
              5, 8, 12, 10, 9, 6, 3, 2, 1, 1, 1)
    rts <- 1:length(ints) + rnorm(length(ints), sd = 0.05)
    chr <- Chromatogram(rts, ints)
    xchr <- findChromPeaks(chr, param = cwp)
    pks <- chromPeaks(xchr)
    res <- .chrom_merge_neighboring_peaks(chr, pks, diffRt = 5, minProp = 0.5)
    expect_true(nrow(res) == 2)
    expect_equal(rownames(res), c(NA_character_, "3"))
    expect_equal(res[1, "rtmin"], unname(pks[1, "rtmin"]))
    expect_equal(res[1, "rtmax"], unname(pks[2, "rtmax"]))
    res <- .chrom_merge_neighboring_peaks(chr, pks, diffRt = 10, minProp = 0.01)
    expect_true(nrow(res) == 1)
    expect_equal(res[1, "rtmin"], unname(pks[1, "rtmin"]))
    expect_equal(res[1, "rtmax"], unname(pks[3, "rtmax"]))
    ## Check "into" calculation.
    pks <- rbind(pks[-c(1, 2), ],
                 c(18, pks[2, "rtmin"], pks[2, "rtmin"] + 4, NA_real_,
                   NA_real_, 3, 3),
                 c(20, pks[2, "rtmax"] - 5, pks[2, "rtmax"], NA_real_,
                   NA_real_, 9, 8))
    res <- .chrom_merge_neighboring_peaks(chr, pks, diffRt = 5, minProp = 0.75)
    expect_equal(unname(res[1, "into"]), unname(chromPeaks(xchr)[2, "into"]))
})
