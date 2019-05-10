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
