test_that("findChromPeaks,Chromatogram,MatchedFilterParam works", {
    od <- filterFile(faahko_od, file = 1)
    mzr <- c(272.1, 272.2)

    chr <- chromatogram(od, mz = mzr)[1, 1]

    mfp <- MatchedFilterParam()
    res <- findChromPeaks(chr, mfp)
    expect_true(is.matrix(res))
    expect_true(nrow(res) == 2)
    mfp <- MatchedFilterParam(fwhm = 60)
    res_2 <- findChromPeaks(chr, mfp)
    expect_true(all(res_2[, "rtmin"] < res[, "rtmin"]))
    expect_true(all(res_2[, "rtmax"] > res[, "rtmax"]))
    expect_true(all(xcms:::.CPEAKS_CHROMPEAKS_REQ_NAMES[-(1:2)] %in%
                    colnames(res)))
})

test_that("findChromPeaks,Chromatogram,CentWaveParam works", {
    od <- filterFile(faahko_od, file = 1)
    mzr <- c(272.1, 272.2)

    chr <- chromatogram(od, mz = mzr)[1, 1]

    cwp <- CentWaveParam()
    res <- findChromPeaks(chr, cwp)
    expect_true(is.matrix(res))
    expect_true(nrow(res) == 2)
    cwp <- CentWaveParam(peakwidth = c(10, 60))
    res_2 <- findChromPeaks(chr, cwp)
    expect_true(nrow(res_2) > nrow(res))
    cwp <- CentWaveParam(snthresh = 5)
    res_2 <- findChromPeaks(chr, cwp)
    expect_true(nrow(res_2) > nrow(res))
    expect_true(all(xcms:::.CPEAKS_CHROMPEAKS_REQ_NAMES[-(1:2)] %in%
                    colnames(res)))
    expect_equal(res[, "rtmin"], res_2[1:2, "rtmin"])
    expect_equal(res[, "rtmax"], res_2[1:2, "rtmax"])
    expect_equal(res[, "rt"], res_2[1:2, "rt"])
    expect_equal(res[, "maxo"], res_2[1:2, "maxo"])
    expect_equal(res[, "into"], res_2[1:2, "into"])
})

