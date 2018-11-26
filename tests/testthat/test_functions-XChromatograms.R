test_that("XChromatograms, as, validator, hasChromPeaks work", {
    chr1 <- Chromatogram(rtime = 1:8,
                         intensity = c(3, 24.2, 343, 32, 3.3, 5, 2, 9))
    chr2 <- Chromatogram(rtime = 1:4, intensity = c(45, 3, 34, 2))
    chr3 <- Chromatogram(rtime = 1:7, intensity = c(12, 34, 54, 34, 23, 2, NA))
    chr4 <- Chromatogram(rtime = 1:3, intensity = c(3, 4, 1))
    chr5 <- Chromatogram(rtime = 1:6, intensity = c(3, 4, 6, 7, 2, 4))
    chr6 <- Chromatogram(rtime = 2:5, intensity = c(3, 65, 43, 12))
    chrs <- Chromatograms(list(chr1, chr2, chr3, chr4, chr5, chr6), nrow = 2)

    expect_error(new("XChromatograms", matrix(list(chr1, chr2), nrow = 1)))

    res <- as(chrs, "XChromatograms")
    expect_true(validObject(res))
    expect_true(is(res, "XChromatograms"))

    colnames(chrs) <- c("A", "B", "C")
    res <- as(chrs, "XChromatograms")
    expect_equal(colnames(res), colnames(chrs))

    xchrs <- XChromatograms(list(chr1, chr2, chr3), ncol = 3)
    expect_equal(ncol(xchrs), 3)
    expect_true(is(xchrs, "XChromatograms"))
    expect_true(all(vapply(xchrs, is, logical(1), "XChromatogram")))

    pks1 <- matrix(c(3, 2, 4, 339.2, 343, NA), nrow = 1,
                   dimnames = list(NULL, xcms:::.CHROMPEAKS_REQ_NAMES))
    pks3 <- matrix(c(3, 2, 4, 145, 54, NA), nrow = 1,
                   dimnames = list(NULL, xcms:::.CHROMPEAKS_REQ_NAMES))
    pks6 <- matrix(c(2, 2, 3, 108, 65, NA), nrow = 1,
                   dimnames = list(NULL, xcms:::.CHROMPEAKS_REQ_NAMES))
    ## With peak matrix.
    xchrs1 <- XChromatograms(list(chr1, chr2, chr3, chr4, chr5, chr6), ncol = 3,
                             chromPeaks = list(pks1, NULL, pks3, NULL, NULL,
                                               pks6))
    expect_true(is(xchrs1, "XChromatograms"))
    expect_equal(unname(hasChromPeaks(xchrs1)[1, ]), c(TRUE, TRUE, FALSE))
    expect_equal(unname(hasChromPeaks(xchrs1)[2, ]), c(FALSE, FALSE, TRUE))

    xchrs1 <- XChromatograms(list(chr1, chr2, chr3, chr4, chr5, chr6), ncol = 3,
                             chromPeaks = list(pks1, NULL, pks3, NULL, NULL,
                                               pks6), byrow = TRUE)
    expect_true(is(xchrs1, "XChromatograms"))
    expect_equal(unname(hasChromPeaks(xchrs1)[1, ]), c(TRUE, FALSE, TRUE))
    expect_equal(unname(hasChromPeaks(xchrs1)[2, ]), c(FALSE, FALSE, TRUE))
    expect_equal(intensity(chr3), intensity(xchrs1[1, 3]))

    ## With XChromatogram objects
    xchr1 <- as(chr1, "XChromatogram")
    xchr3 <- as(chr3, "XChromatogram")
    xchr4 <- as(chr4, "XChromatogram")
    xchr6 <- as(chr6, "XChromatogram")
    chromPeaks(xchr1) <- pks1
    chromPeaks(xchr3) <- pks3
    chromPeaks(xchr6) <- pks6

    xchrs2 <- XChromatograms(list(xchr1, xchr4, xchr3, xchr6), ncol = 2)
    expect_equal(unname(hasChromPeaks(xchrs2)[1, ]), c(TRUE, TRUE))
    expect_equal(unname(hasChromPeaks(xchrs2)[2, ]), c(FALSE, TRUE))
    expect_equal(chromPeaks(xchrs2[1, 2]), chromPeaks(xchr3))
})
