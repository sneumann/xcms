test_that(".validXChromatogram works", {
    xc <- new("XChromatogram")
    expect_true(.validXChromatogram(xc))
    xc@chromPeaks <- matrix(1:10, ncol = 5)
    expect_true(is.character(.validXChromatogram(xc)))
    mat <- matrix("A", ncol = 6, nrow = 2)
    xc@chromPeaks <- mat
    expect_true(is.character(.validXChromatogram(xc)))
    mat <- matrix(ncol = 6, nrow = 2)
    colnames(mat) <- xcms:::.CHROMPEAKS_REQ_NAMES
    mat[, "rtmin"] <- c(3, 3)
    mat[, "rtmax"] <- c(3, 4)
    xc@chromPeaks <- mat
    expect_true(.validXChromatogram(xc))
    xc@chromPeaks[, "rtmin"] <- c(4, 3)
    expect_true(is.character(.validXChromatogram(xc)))
})

test_that("XChromatogram works", {
    chr <- Chromatogram(rtime = 1:10, intensity = 1:10)
    xc <- XChromatogram(chr)
    expect_true(nrow(xc@chromPeaks) == 0)
    expect_equal(rtime(xc), 1:10)

    xc <- XChromatogram()
    expect_true(nrow(xc@chromPeaks) == 0)
    expect_error(XChromatogram(chr, 4))

    pks <- matrix(nrow = 2, ncol = length(.CHROMPEAKS_REQ_NAMES),
                  dimnames = list(character(), .CHROMPEAKS_REQ_NAMES))
    pks[, "rtmin"] <- c(2, 4)
    pks[, "rtmax"] <- c(3, 5)
    xc <- XChromatogram(chr, pks)
    expect_equal(xc@chromPeaks[, "rtmin"], c(2, 4))
})
