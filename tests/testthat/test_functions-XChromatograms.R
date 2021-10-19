test_that("XChromatograms, as, validator, hasChromPeaks work", {
    skip_on_os(os = "windows", arch = "i386")

    chr1 <- Chromatogram(rtime = 1:8,
                         intensity = c(3, 24.2, 343, 32, 3.3, 5, 2, 9))
    chr2 <- Chromatogram(rtime = 1:4, intensity = c(45, 3, 34, 2))
    chr3 <- Chromatogram(rtime = 1:7, intensity = c(12, 34, 54, 34, 23, 2, NA))
    chr4 <- Chromatogram(rtime = 1:3, intensity = c(3, 4, 1))
    chr5 <- Chromatogram(rtime = 1:6, intensity = c(3, 4, 6, 7, 2, 4))
    chr6 <- Chromatogram(rtime = 2:5, intensity = c(3, 65, 43, 12))
    chrs <- MChromatograms(list(chr1, chr2, chr3, chr4, chr5, chr6), nrow = 2)

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
                   dimnames = list(NULL, .CHROMPEAKS_REQ_NAMES))
    pks3 <- matrix(c(3, 2, 4, 145, 54, NA), nrow = 1,
                   dimnames = list(NULL, .CHROMPEAKS_REQ_NAMES))
    pks6 <- matrix(c(2, 2, 3, 108, 65, NA), nrow = 1,
                   dimnames = list(NULL, .CHROMPEAKS_REQ_NAMES))
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

test_that(".subset_chrom_peaks_xchromatograms works", {
    skip_on_os(os = "windows", arch = "i386")

    ## Matrix is:    with elements
    ## A B C D       2 1 3 1
    ## E F G H       1 4 0 2
    ## I J K L       3 2 1 3

    testm <- data.frame(el = c("A", "A", "E", "I", "I", "I",
                               "B", "F", "F", "F", "F", "J", "J",
                               "C", "C", "C", "K",
                               "D", "H", "H", "L", "L", "L"),
                        row = c(1, 1, 2, 3, 3, 3,
                                1, 2, 2, 2, 2, 3, 3,
                                1, 1, 1, 3,
                                1, 2, 2, 3, 3, 3),
                        column = c(1, 1, 1, 1, 1, 1,
                                   2, 2, 2, 2, 2, 2, 2,
                                   3, 3, 3, 3,
                                   4, 4, 4, 4, 4, 4),
                        stringsAsFactors = FALSE)

    res <- .subset_chrom_peaks_xchromatograms(testm, i = 2:3, j = 2:4)
    expect_equal(res$el,  c("F", "F", "F", "F", "H", "H",
                            "J", "J", "K", "L", "L", "L"))
    expect_equal(res$row, c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2))
    expect_equal(res$column, c(1, 1, 1, 1, 3, 3, 1, 1, 2, 3, 3, 3))

    res <- .subset_chrom_peaks_xchromatograms(testm, i = c(3, 1), j = 2:4)
    expect_equal(res$el, c("J", "J", "K", "L", "L", "L",
                           "B", "C", "C", "C", "D"))
    expect_equal(res$row, c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2))
    expect_equal(res$column, c(1, 1, 2, 3, 3, 3, 1, 2, 2, 2, 3))

    res <- .subset_chrom_peaks_xchromatograms(testm, i = 2, j = c(3, 2))
    expect_equal(res$el, c("F", "F", "F", "F"))
    expect_equal(res$row, c(1, 1, 1, 1))
    expect_equal(res$column, c(2, 2, 2, 2))

    res <- .subset_chrom_peaks_xchromatograms(testm, i = 2, j = c(3, 2, 4))
    expect_equal(res$el, c("F", "F", "F", "F", "H", "H"))
    expect_equal(res$row, c(1, 1, 1, 1, 1, 1))
    expect_equal(res$column, c(2, 2, 2, 2, 3, 3))

    mzr <- matrix(c(335, 335, 344, 344), ncol = 2, byrow = TRUE)
    chrs <- chromatogram(od_x, mz = mzr)
    chrs <- findChromPeaks(chrs, param = CentWaveParam())
    prm <- PeakDensityParam(sampleGroups = c(1, 1, 1))
    chrs <- groupChromPeaks(chrs, param = prm)

    pks <- chromPeaks(chrs)
    rownames(pks) <- letters[1:nrow(pks)]
    res <- .subset_chrom_peaks_xchromatograms(pks, j = c(3, 1, 2))
    expect_equal(rownames(res), c("f", "g", "a", "b", "c", "d", "e",
                                  "l", "m", "h", "i", "j", "k"))
    res <- .subset_chrom_peaks_xchromatograms(pks, i = c(2, 1), j = c(2, 1, 3))
    expect_equal(rownames(res), c("j", "k", "h", "i", "l", "m",
                                  "e", "a", "b", "c", "d", "f", "g"))
})

test_that(".subset_features_on_chrom_peaks works", {
    skip_on_os(os = "windows", arch = "i386")

    chrs <- as(od_chrs, "XChromatograms")
    chrs <- findChromPeaks(chrs, param = CentWaveParam())
    prm <- PeakDensityParam(sampleGroups = c(1, 1, 1))
    chrs <- groupChromPeaks(chrs, param = prm)

    fts <- featureDefinitions(chrs)
    pks <- chromPeaks(chrs)
    res <- .subset_features_on_chrom_peaks(fts, pks, pks)
    expect_equal(fts, res)

    pks_sub <- pks[1:2, ]
    res <- .subset_features_on_chrom_peaks(fts, pks, pks_sub)
    expect_equal(rownames(res), "FT1")
    expect_equal(res$peakidx, list(1))

    pks_sub <- pks[c(8, 10, 12), ]
    res <- .subset_features_on_chrom_peaks(fts, pks, pks_sub)
    expect_equal(rownames(res), "FT3")
    expect_equal(res$peakidx, list(1:3))

    pks_sub <- pks[2, , drop = FALSE]
    res <- .subset_features_on_chrom_peaks(fts, pks, pks_sub)
    expect_true(nrow(res) == 0)

    pks_sub <- rbind(pks, pks[pks[, "row"] == 1, ])
    pks_sub[(nrow(pks)+1):nrow(pks_sub), "row"] <- 3
    fts <- rbind(fts, fts[fts$row == 1, ])
    fts$row[5:6] <- 3
    res <- .subset_features_on_chrom_peaks(fts, pks, pks_sub)
    expect_true(nrow(res) == 6)
    expect_equal(pks_sub[res$peakidx[[1]], "into"],
                 pks_sub[res$peakidx[[5]], "into"])
    expect_equal(pks_sub[res$peakidx[[2]], "into"],
                 pks_sub[res$peakidx[[6]], "into"])
})

test_that(".plot_chrom_peak_density works", {
    skip_on_os(os = "windows", arch = "i386")

    chrs <- as(od_chrs, "XChromatograms")
    chrs <- findChromPeaks(chrs, param = CentWaveParam())

    pks_1 <- chromPeaks(chrs[1, ])
    prm <- PeakDensityParam(sampleGroups = c(1, 1, 1))
    .plot_chrom_peak_density(pks_1, param = prm, lwd = 2)
    expect_warning(.plot_chrom_peak_density(pks_1, param = prm,
                                            peakCol = c(1, 2)))
    ## An individual color for each point.
    .plot_chrom_peak_density(pks_1, param = prm, peakCol = 1:nrow(pks_1),
                             peakPch = 16)

    chrs <- groupChromPeaks(chrs, param = prm)
    pks_1 <- chromPeaks(chrs[1, ])
    fts_1 <- featureDefinitions(chrs[1, ])
    .plot_chrom_peak_density(pks_1, fts = fts_1, param = prm,
                             peakCol = 1:nrow(pks_1), peakPch = 16,
                             simulate = FALSE)
})
