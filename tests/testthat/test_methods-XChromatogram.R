test_that("show,XChromatogram works", {
    show(XChromatogram())
})

test_that("chromPeaks,chromPeaks<-,XChromatogram works", {
    chr <- Chromatogram(rtime = 1:10,
                        intensity = c(4, 12, 18, 24, 23, 18, 15, 3, 2, 5))
    xchr <- XChromatogram(chr)
    expect_true(nrow(chromPeaks(xchr)) == 0)
    pks <- matrix(nrow = 4, ncol = 6)
    colnames(pks) <- xcms:::.CHROMPEAKS_REQ_NAMES
    pks[1, ] <- c(4, 2, 8, 24, NA, NA)
    pks[2, ] <- c(3, 2, 7, 24, NA, NA)
    pks[3, ] <- c(9, 7, 10, 2, NA, NA)
    pks[4, ] <- c(8, 5, 10, 3, NA, NA)
    expect_error(chromPeaks(xchr) <- 4)
    chromPeaks(xchr) <- pks
    expect_equal(chromPeaks(xchr), pks)

    expect_true(nrow(chromPeaks(xchr, rt = c(20, 30))) == 0)
    expect_equal(chromPeaks(xchr, rt = c(2, 7)), pks)
    expect_equal(chromPeaks(xchr, rt = c(2, 7), type = "apex_within"),
                 pks[1:2, ])
    expect_equal(chromPeaks(xchr, rt = c(2, 7), type = "within"),
                 pks[2, , drop = FALSE])
})
