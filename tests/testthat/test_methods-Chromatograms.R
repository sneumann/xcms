test_that("findChromPeaks,Chromatograms works", {
    mzr <- matrix(c(335, 335, 344, 344), ncol = 2, byrow = TRUE)
    chrs <- chromatogram(od_x, mz = mzr)
    res <- findChromPeaks(chrs, param = CentWaveParam())
    expect_true(is(res, "XChromatograms"))
    expect_equal(intensity(res[1, 2]), intensity(chrs[1, 2]))
    expect_equal(intensity(res[2, 3]), intensity(chrs[2, 3]))
    expect_true(length(res@.processHistory) == 1)
    expect_warning(res_2 <- findChromPeaks(chrs, param = CentWaveParam(sn = 50)))
    expect_true(nrow(chromPeaks(res)) > nrow(chromPeaks(res_2)))

    ## MatchedFilter
    res_m <- findChromPeaks(chrs, param = MatchedFilterParam())
    expect_true(is(res_m, "XChromatograms"))
    expect_true(nrow(chromPeaks(res_m)) < nrow(chromPeaks(res)))

    ## on a XChromatograms
    res_3 <- findChromPeaks(res_2, param = CentWaveParam())
    expect_true(length(res_3@.processHistory) == 1)
    expect_equal(chromPeaks(res), chromPeaks(res_3))
})
