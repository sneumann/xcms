test_that("chromPeaks,XChromatograms works", {
    pks1 <- matrix(c(3, 2, 4, 339.2, 343, NA), nrow = 1,
                   dimnames = list(NULL, xcms:::.CHROMPEAKS_REQ_NAMES))
    chr1 <- XChromatogram(
        rtime = 1:8, intensity = c(3, 24.2, 343, 32, 3.3, 5, 2, 9),
        chromPeaks = pks1)
    chr2 <- XChromatogram(rtime = 1:4, intensity = c(45, 3, 34, 2))
    pks3 <- matrix(c(3, 2, 4, 145, 54, NA), nrow = 1,
                   dimnames = list(NULL, xcms:::.CHROMPEAKS_REQ_NAMES))
    chr3 <- XChromatogram(
        rtime = 1:7, intensity = c(12, 34, 54, 34, 23, 2, NA),
        chromPeaks = pks3)
    chr4 <- XChromatogram(rtime = 1:3, intensity = c(3, 4, 1))
    chr5 <- XChromatogram(rtime = 1:6, intensity = c(3, 4, 6, 7, 2, 4))
    pks6 <- matrix(c(2, 2, 3, 108, 65, NA, 3, 5, 7, 123, 4, NA),
                   nrow = 2, byrow = TRUE,
                   dimnames = list(NULL, xcms:::.CHROMPEAKS_REQ_NAMES))
    chr6 <- XChromatogram(
        rtime = 2:5, intensity = c(3, 65, 43, 12),
        chromPeaks = pks6)
    xchrs <- XChromatograms(list(chr1, chr2, chr3, chr4, chr5, chr6), nrow = 2)

    res <- chromPeaks(xchrs)
    expect_equal(res[1, 1:6, drop = FALSE], pks1)
    expect_equal(res[2, 1:6, drop = FALSE], pks3)
    expect_equal(res[3:4, 1:6, drop = FALSE], pks6)
    expect_equal(res[, "row"], c(1, 1, 2, 2))
    expect_equal(res[, "column"], c(1, 2, 3, 3))

    xchrs <- XChromatograms(list(chr4, chr5, chr2))
    res <- chromPeaks(xchrs)
    expect_true(nrow(res) == 0)
    expect_equal(colnames(res), c("rt", "rtmin", "rtmax", "into", "maxo", "sn",
                                  "row", "column"))

    xchrs <- XChromatograms(list(chr2, chr4, chr5, chr1))
    res <- chromPeaks(xchrs)
    expect_equal(res[, 1:6, drop = FALSE], pks1)
    expect_equal(unname(res[, "row"]), c(4))
    expect_equal(unname(res[, "column"]), c(1))

    xchrs <- XChromatograms(list(chr2, chr4, chr5, chr1), ncol = 2)
    res <- chromPeaks(xchrs)
    expect_equal(res[, 1:6, drop = FALSE], pks1)
    expect_equal(unname(res[, "row"]), c(2))
    expect_equal(unname(res[, "column"]), c(2))
})

test_that("filterRt, filterMz for XChromatograms works", {
    chr1 <- XChromatogram()
    chr2 <- XChromatogram(rtime = 1:4, intensity = c(45, 3, 34, 2))
    pks3 <- matrix(c(3, 2, 4, 145, 54, NA), nrow = 1,
                   dimnames = list(NULL, xcms:::.CHROMPEAKS_REQ_NAMES))
    chr3 <- XChromatogram(
        rtime = 1:7, intensity = c(12, 34, 54, 34, 23, 2, NA),
        chromPeaks = pks3)
    chr4 <- XChromatogram(rtime = 1:3, intensity = c(3, 4, 1))

    xchrs <- XChromatograms()
    expect_equal(xchrs, filterRt(xchrs))
    expect_equal(xchrs, filterMz(xchrs))

    xchrs <- XChromatograms(list(chr1, chr2))
    res <- filterMz(xchrs)
    expect_equal(xchrs, res)
    res <- filterRt(xchrs, rt = c(3, 4))
    expect_equal(res[1, 1], chr1)
    expect_equal(rtime(res[2, 1]), c(3, 4))
    expect_equal(intensity(res[2, 1]), c(34, 2))

    xchrs <- XChromatograms(list(chr3, chr4), nrow = 1)
    res <- filterMz(xchrs)
    expect_equal(xchrs, res)
    res <- filterRt(xchrs, rt = c(6, 7))
    expect_true(length(rtime(res[1, 2])) == 0)
    expect_equal(rtime(res[1, 1]), 6:7)
    expect_equal(intensity(res[1, 1]), c(2, NA))
    expect_true(nrow(chromPeaks(res)) == 0)

    res <- filterRt(xchrs, rt = c(3, 6))
    expect_equal(rtime(res[1, 1]), 3:6)
    expect_equal(intensity(res[1, 1]), c(54, 34, 23, 2))
    expect_true(nrow(chromPeaks(res)) == 1)
    expect_equal(chromPeaks(res[1, 1]), pks3)
})

test_that("plot,XChromatogram works", {
    mzr <- matrix(c(335, 335, 344, 344), ncol = 2, byrow = TRUE)
    rtr <- matrix(c(2700, 2900, 2600, 2750), ncol = 2, byrow = TRUE)

    ## Full rt range.
    xchr_rt <- chromatogram(xod_chr, mz = mzr)

    cols_smple <- c("#ff000060", "#00ff0060", "#0000ff60")
    plot(as(xchr_rt[1, ], "Chromatograms"), col = cols_smple)
    ## Plotting the data is fine, just as above.
    ## Then we have to loop over each chromatogram...
    x <- xchr_rt[1, ]
    pks <- chromPeaks(x)
    .add_chromatogram_peaks(x, pks, type = "rectangle", bg = rep(NA, nrow(pks)),
                            col = cols_smple[pks[, "sample"]])
    .add_chromatogram_peaks(x, pks, type = "point", bg = rep(NA, nrow(pks)),
                            col = cols_smple[pks[, "sample"]], pch = 15)
    .add_chromatogram_peaks(x, pks, type = "polygon",
                            bg = cols_smple[pks[, "sample"]],
                            col = cols_smple[pks[, "sample"]])
    x <- xchr_rt[2, ]
    pks <- chromPeaks(x)
    plot(as(x, "Chromatograms"), col = cols_smple, lwd = 2)
    .add_chromatogram_peaks(x, pks, type = "rectangle", bg = rep(NA, nrow(pks)),
                            col = cols_smple[pks[, "sample"]])
    .add_chromatogram_peaks(x, pks, type = "point", bg = rep(NA, nrow(pks)),
                            col = cols_smple[pks[, "sample"]], pch = 15)
    .add_chromatogram_peaks(x, pks, type = "polygon",
                            bg = cols_smple[pks[, "sample"]],
                            col = cols_smple[pks[, "sample"]])

    plot(xchr_rt, peakCol = cols_smple[chromPeaks(xchr_rt)[, "sample"]],
         peakBg = cols_smple[chromPeaks(xchr_rt)[, "sample"]], xlab = "RT")
    xsub <- xchr_rt[2, ]
    ## Use one color per peak
    library(RColorBrewer)
    cls <- paste0(brewer.pal(nrow(chromPeaks(xsub)), "Dark2"), 40)
    plot(xsub, peakBg = cls)

    ## Narrow on rt.
    xchr <- chromatogram(xod_chr, mz = mzr, rt = rtr)
    plot(xchr)
})

test_that("processHistory,XChromatograms works", {
    mzr <- matrix(c(335, 335, 344, 344), ncol = 2, byrow = TRUE)
    xchrs <- chromatogram(xod_xgrg, mz = mzr, rt = c(2600, 3600))
    res <- processHistory(xchrs)
    expect_equal(length(res), 4)
    res <- processHistory(xchrs, type = .PROCSTEP.PEAK.DETECTION)
    expect_equal(length(res), 1)

    xchrs2 <- findChromPeaks(xchrs, param = CentWaveParam())
    expect_equal(length(processHistory(xchrs2)), 5)
    expect_equal(length(processHistory(xchrs2,
                                       type = .PROCSTEP.PEAK.DETECTION)), 2)
})

test_that("groupChromPeaks,XChromatograms,PeakDensityParam works", {
    mzr <- matrix(c(335, 335, 344, 344), ncol = 2, byrow = TRUE)
    xchrs <- chromatogram(xod_xgrg, mz = mzr, rt = c(2600, 3600))
    param <- PeakDensityParam(sampleGroups = c(1, 1, 1))
    res <- groupChromPeaks(xchrs, param = param)
    expect_true(hasFeatures(res))
    expect_true(nrow(res@featureDefinitions) == 1)

    param <- PeakDensityParam(sampleGroups = c(1, 2, 3))
    res <- groupChromPeaks(xchrs, param = param)
    expect_true(nrow(res@featureDefinitions) == 2)
    expect_true(length(processHistory(res)) == 5)

    res <- groupChromPeaks(res,
                           param = PeakDensityParam(sampleGroups = c(1, 1, 1)))
    expect_true(length(processHistory(res)) == 5)
    expect_true(nrow(res@featureDefinitions) == 1)
})

test_that("dropFeatureDefinitions,XChromatograms works", {
    mzr <- matrix(c(335, 335, 344, 344), ncol = 2, byrow = TRUE)
    xchrs <- chromatogram(xod_xgrg, mz = mzr, rt = c(2600, 3600))
    param <- PeakDensityParam(sampleGroups = c(1, 1, 1))
    res <- groupChromPeaks(xchrs, param = param)
    expect_true(hasFeatures(res))
    expect_true(nrow(res@featureDefinitions) == 1)
    expect_true(length(res@.processHistory) == 5)
    res <- dropFeatureDefinitions(res)
    expect_false(hasFeatures(res))
    expect_true(length(res@.processHistory) == 4)
    expect_equal(res, xchrs)
})
