test_that("findChromPeaks,Chromatograms works", {
    mzr <- matrix(c(335, 335, 344, 344), ncol = 2, byrow = TRUE)
    chrs <- chromatogram(od_x, mz = mzr)
    res <- findChromPeaks(chrs, param = CentWaveParam())
    expect_true(is(res, "XChromatograms"))
    expect_equal(intensity(res[1, 2]), intensity(chrs[1, 2]))
    expect_equal(intensity(res[2, 3]), intensity(chrs[2, 3]))
    expect_true(length(res@.processHistory) == 1)
    res_2 <- findChromPeaks(chrs, param = CentWaveParam(sn = 50))
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

test_that(".correlate_chromatograms_self works", {
    chr1 <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
                         intensity = c(5, 29, 50, NA, 100, 12, 3, 4, 1, 3))
    chr2 <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
                         intensity = c(80, 50, 20, 10, 9, 4, 3, 4, 1, 3))
    chr3 <- Chromatogram(rtime = 3:9 + rnorm(7, sd = 0.3),
                         intensity = c(53, 80, 130, 15, 5, 3, 2))
    res <- .correlate_chromatograms_self(list(chr1, chr2, chr3))

    expect_equal(nrow(res), ncol(res))
    expect_equal(nrow(res), 3)
    expect_equal(res[c(1, 5, 9)], c(1, 1, 1))
    expect_true(res[1, 3] > 0.9)
    expect_true(res[1, 2] < 0.5)

    chrs <- Chromatograms(list(chr1, chr2, chr3))
    expect_equal(.correlate_chromatograms_self(chrs), res)

    chrs <- Chromatograms(list(chr1, chr2, chr3, chr1), ncol = 2)
    expect_error(.correlate_chromatograms_self(chrs), "single column")
})

test_that(".correlate_chromatograms works", {
    chr1 <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
                         intensity = c(5, 29, 50, NA, 100, 12, 3, 4, 1, 3))
    chr2 <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
                         intensity = c(80, 50, 20, 10, 9, 4, 3, 4, 1, 3))
    chr3 <- Chromatogram(rtime = 3:9 + rnorm(7, sd = 0.3),
                         intensity = c(53, 80, 130, 15, 5, 3, 2))
    res <- .correlate_chromatograms(list(chr1), list(chr2, chr3))
    expect_equal(ncol(res), 2)
    expect_equal(nrow(res), 1)
    expect_true(res[1, 1] < 0.5)
    expect_true(res[1, 2] > 0.9)

    expect_equal(res, .correlate_chromatograms(list(chr1), list(chr2, chr3),
                                               full = FALSE))

    res <- .correlate_chromatograms(list(chr1, chr2, chr3),
                                    list(chr1, chr2, chr3))
    expect_equal(res[1, 2], correlate(chr1, chr2))
    expect_equal(res[1, 3], correlate(chr1, chr3))
    expect_equal(res[2, 1], correlate(chr2, chr1))
    expect_equal(res[2, 3], correlate(chr2, chr3))
    expect_equal(res[3, 1], correlate(chr3, chr1))

    res <- .correlate_chromatograms(list(chr1, chr2, chr3),
                                    list(chr1, chr2, chr3),
                                    full = FALSE)
    expect_equal(res[1, 2], correlate(chr1, chr2))
    expect_equal(res[1, 3], correlate(chr1, chr3))
    expect_equal(res[2, 1], NA_real_)
    expect_equal(res[2, 3], correlate(chr2, chr3))
    expect_equal(res[3, 1], NA_real_)

    res <- .correlate_chromatograms(list(chr1, chr2, chr3),
                                    list(chr1, chr2),
                                    full = FALSE)

    chrs <- Chromatograms(list(chr1, chr2, chr3, chr1), ncol = 2)
    expect_error(.correlate_chromatograms(chrs, list(chr1, chr2)), "single column")
    expect_error(.correlate_chromatograms(list(chr1, chr2), chrs), "single column")
})

test_that("correlate,Chromatograms works", {
    chr1 <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
                         intensity = c(5, 29, 50, NA, 100, 12, 3, 4, 1, 3))
    chr2 <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
                         intensity = c(80, 50, 20, 10, 9, 4, 3, 4, 1, 3))
    chr3 <- Chromatogram(rtime = 3:9 + rnorm(7, sd = 0.3),
                         intensity = c(53, 80, 130, 15, 5, 3, 2))
    chrs <- Chromatograms(list(chr1, chr2, chr3))

    res <- correlate(chrs)
    expect_true(nrow(res) == 3)
    expect_true(ncol(res) == 3)
    expect_true(res[1, 3] > 0.9)
    expect_true(res[1, 2] < 0.5)

    res_2 <- correlate(chrs, chrs)
    expect_equal(res_2, res)

    res <- correlate(chrs, full = FALSE)
    expect_true(is.na(res[2, 1]))
    expect_true(is.na(res[3, 1]))

})

test_that("removeIntensity,Chromatograms works", {
    chr1 <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
                         intensity = c(5, 29, 50, NA, 100, 12, 3, 4, 1, 3))
    chr2 <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
                         intensity = c(80, 50, 20, 10, 9, 4, 3, 4, 1, 3))
    chr3 <- Chromatogram(rtime = 3:9 + rnorm(7, sd = 0.3),
                         intensity = c(53, 80, 130, 15, 5, 3, 2))
    chrs <- Chromatograms(list(chr1, chr2, chr3))

    res <- removeIntensity(chrs)
    expect_equal(res, chrs)

    res <- removeIntensity(chrs, threshold = 20)
    expect_equal(intensity(res[1, 1]), c(NA_real_, 29, 50, NA_real_, 100,
                                         NA_real_, NA_real_, NA_real_, NA_real_,
                                         NA_real_))
    expect_equal(intensity(res[3, 1]), c(53, 80, 130, NA_real_, NA_real_,
                                         NA_real_, NA_real_))

    chrs <- Chromatograms(list(chr1, chr2, chr2, chr3), ncol = 2)
    res <- removeIntensity(chrs, threshold = 20)
    expect_equal(intensity(res[2, 2]), c(53, 80, 130, NA_real_, NA_real_,
                                         NA_real_, NA_real_))
})
