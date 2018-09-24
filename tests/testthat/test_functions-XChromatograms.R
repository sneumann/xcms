test_that("XChromatograms class works", {
    library(testthat)
    library(xcms)

    c1 <- Chromatogram(rtime = c(1, 2, 3, 4), intensity = c(5, 3, 5, 7))
    c2 <- Chromatogram(rtime = c(1, 2, 3, 4), intensity = c(7, 9, 0, 4))
    cs <- Chromatograms(list(c1, c2), ncol = 2)

    tst <- as(cs, "XChromatograms")
    expect_true(validObject(tst))

    dummymat <- matrix(ncol = 6, nrow = 2, 1:12)
    tst@chromPeaks <- matrix(list(dummymat, matrix(ncol = 6, nrow = 0)),
                             ncol = 2)
    expect_error(validObject(tst))
    colnames(dummymat) <- c("rt", "rtmin", "rtmax", "into", "maxo", "sn")
    tst@chromPeaks[1, 1][[1]] <- dummymat
    colnames(tst@chromPeaks[1, 2][[1]]) <- c("rt", "rtmin", "rtmax",
                                             "into", "maxo", "sn")
    expect_true(validObject(tst))
    tst@chromPeaks <- matrix(list(dummymat, dummymat), ncol = 1)
    expect_error(validObject(tst))
    tst@chromPeaks <- matrix(list(dummymat, dummymat), ncol = 2)
    expect_true(validObject(tst))

    ## constructor
    tst <- XChromatograms(list(c1), chromPeaks = list(dummymat))
    expect_true(validObject(tst))
    expect_error(XChromatograms(list(c1, c2), chromPeaks = list(dummymat)))
})
