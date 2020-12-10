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
    expect_equal(.correlate_chromatogram(ch1, ch2),
                 cor(intensity1, intensity2))
    expect_equal(.correlate_chromatogram(ch1, ch2),
                 .correlate_chromatogram(ch2, ch1))

    expect_equal(
        .correlate_chromatogram(ch1, ch2, align = "approx"),
        .correlate_chromatogram(ch2, ch1, align = "approx"),
        tolerance = 0.0001
    )

    ch2 <- filterRt(ch2, rt = c(5.1, 23))
    res <- .correlate_chromatogram(ch1, ch2)
    expect_equal(res, cor(intensity(ch2), intensity(ch1)[1:length(ch2)]))

    res <- .correlate_chromatogram(ch2, ch1, align = "approx",
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
    pkd <- chromPeakData(xchr)
    pkd$index <- 1:nrow(pkd)
    res <- .chrom_merge_neighboring_peaks(chr, chromPeaks(xchr), pkd)
    expect_true(all(names(res) %in% c("chromPeaks", "chromPeakData")))
    rownames(res$chromPeaks) <- NULL
    expect_equal(res$chromPeaks, chromPeaks(xchr))
    res <- .chrom_merge_neighboring_peaks(chr, chromPeaks(xchr), pkd,
                                          diffRt = 5)
    rownames(res$chromPeaks) <- NULL
    expect_equal(res$chromPeaks, chromPeaks(xchr))
    res <- .chrom_merge_neighboring_peaks(
        chr, chromPeaks(xchr)[1, , drop = FALSE],
        pkd[1, ],
        diffRt = 5)
    rownames(res$chromPeaks) <- NULL
    expect_equal(res$chromPeaks, chromPeaks(xchr)[1, , drop = FALSE])

    set.seed(123)
    ints <- c(0.5, 1, 1, 3, 6, 9, 12, 13, 11, 7, 6, 5.5, 5.2, 5, 5.2, 5.4,
              5.7, 6, 8, 9, 8, 6, 3, 2, 1.3, 1, 0.7, 0.5, 1, 1, 1, 0.5, 3,
              5, 8, 14, 10, 9, 6, 3, 2, 1, 1, 1)
    rts <- 1:length(ints) + rnorm(length(ints), sd = 0.05)
    chr <- Chromatogram(rts, ints)
    xchr <- findChromPeaks(chr, param = cwp)
    pks <- chromPeaks(xchr)
    pkd <- chromPeakData(xchr)
    pkd$index <- 1:nrow(pkd)
    res <- .chrom_merge_neighboring_peaks(
        chr, pks, pkd, diffRt = 5, minProp = 0.5)
    expect_true(nrow(res$chromPeaks) == 2)
    expect_equal(rownames(res$chromPeaks), c(NA_character_, "3"))
    expect_equal(res$chromPeaks[1, "rtmin"], unname(pks[1, "rtmin"]))
    expect_equal(res$chromPeaks[1, "rtmax"], unname(pks[2, "rtmax"]))
    expect_equal(nrow(res$chromPeaks), nrow(res$chromPeakData))
    expect_equal(res$chromPeakData$index, c(1L, 3L))

    res <- .chrom_merge_neighboring_peaks(
        chr, pks, pkd, diffRt = 10, minProp = 0.01)
    expect_true(nrow(res$chromPeaks) == 1)
    expect_equal(res$chromPeaks[1, "rtmin"], unname(pks[1, "rtmin"]))
    expect_equal(res$chromPeaks[1, "rtmax"], unname(pks[3, "rtmax"]))
    expect_equal(res$chromPeakData$index, 3L)

    ## Check "into" calculation.
    pks <- rbind(pks[-c(1, 2), ],
                 c(18, pks[2, "rtmin"], pks[2, "rtmin"] + 4, NA_real_,
                   NA_real_, 3, 3),
                 c(20, pks[2, "rtmax"] - 5, pks[2, "rtmax"], NA_real_,
                   NA_real_, 9, 8))
    res <- .chrom_merge_neighboring_peaks(
        chr, pks, pkd, diffRt = 5, minProp = 0.75)
    expect_equal(unname(res$chromPeaks[1, "into"]),
                 unname(chromPeaks(xchr)[2, "into"]))
})
