test_that("findChromPeaks,Chromatogram,MatchedFilterParam works", {
    skip_on_os(os = "windows", arch = "i386")

    od <- filterFile(faahko_od, file = 1)
    mzr <- c(272.1, 272.2)

    chr <- chromatogram(od, mz = mzr)[1, 1]

    mfp <- MatchedFilterParam()
    res <- findChromPeaks(chr, mfp)
    expect_true(is(res, "XChromatogram"))
    expect_true(nrow(chromPeaks(res)) == 2)
    mfp <- MatchedFilterParam(fwhm = 60)
    res_2 <- findChromPeaks(chr, mfp)
    expect_true(all(chromPeaks(res_2)[, "rtmin"] < chromPeaks(res)[, "rtmin"]))
    expect_true(all(chromPeaks(res_2)[, "rtmax"] > chromPeaks(res)[, "rtmax"]))
    expect_true(validObject(res))
})

test_that("findChromPeaks,Chromatogram,CentWaveParam works", {
    skip_on_os(os = "windows", arch = "i386")

    od <- filterFile(faahko_od, file = 1)
    mzr <- c(272.1, 272.2)

    chr <- chromatogram(od, mz = mzr)[1, 1]

    cwp <- CentWaveParam()
    res <- findChromPeaks(chr, cwp)
    expect_true(is(res, "XChromatogram"))
    expect_true(nrow(chromPeaks(res)) == 2)
    cwp <- CentWaveParam(peakwidth = c(10, 60))
    res_2 <- findChromPeaks(chr, cwp)
    expect_true(nrow(chromPeaks(res_2)) > nrow(chromPeaks(res)))
    cwp <- CentWaveParam(snthresh = 5)
    res_2 <- findChromPeaks(chr, cwp)
    expect_true(nrow(chromPeaks(res_2)) > nrow(chromPeaks(res)))
    res <- chromPeaks(res)
    res_2 <- chromPeaks(res_2)
    expect_equal(res[, "rtmin"], res_2[1:2, "rtmin"])
    expect_equal(res[, "rtmax"], res_2[1:2, "rtmax"])
    expect_equal(res[, "rt"], res_2[1:2, "rt"])
    expect_equal(res[, "maxo"], res_2[1:2, "maxo"])
    expect_equal(res[, "into"], res_2[1:2, "into"])
})

test_that("removeIntensity,Chromatogram works", {
    skip_on_os(os = "windows", arch = "i386")

    chr <- Chromatogram(rtime = c(1, 2, 3, 4, 5, 6, 7),
                        intensity = c(NA_real_, 13, 16, 22, 34, 15, 6))
    res <- removeIntensity(chr)
    expect_equal(chr, res)
    res <- removeIntensity(chr, threshold = 20)
    expect_equal(intensity(res), c(NA_real_, NA_real_, NA_real_, 22, 34,
                                   NA_real_, NA_real_))
})
