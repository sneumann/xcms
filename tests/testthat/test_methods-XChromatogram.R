test_that("show,XChromatogram works", {
    skip_on_os(os = "windows", arch = "i386")

    show(XChromatogram())
})

test_that("chromPeaks and chromPeakData for XChromatogram work", {
    skip_on_os(os = "windows", arch = "i386")

    chr <- Chromatogram(rtime = 1:10,
                        intensity = c(4, 12, 18, 24, 23, 18, 15, 3, 2, 5))
    xchr <- as(chr, "XChromatogram")
    expect_true(nrow(chromPeaks(xchr)) == 0)
    pks <- matrix(nrow = 4, ncol = 6)
    colnames(pks) <- .CHROMPEAKS_REQ_NAMES
    pks[1, ] <- c(4, 2, 8, 24, NA, NA)
    pks[2, ] <- c(3, 2, 7, 24, NA, NA)
    pks[3, ] <- c(9, 7, 10, 2, NA, NA)
    pks[4, ] <- c(8, 5, 10, 3, NA, NA)
    expect_error(chromPeaks(xchr) <- 4)
    chromPeaks(xchr) <- pks
    expect_equal(chromPeaks(xchr), pks)
    expect_equal(nrow(chromPeakData(xchr)), nrow(pks))
    expect_equal(chromPeakData(xchr)$is_filled, rep(FALSE, nrow(pks)))
    expect_equal(chromPeakData(xchr)$ms_level, rep(1L, nrow(pks)))
    chromPeakData(xchr)$id <- "a"
    expect_equal(chromPeakData(xchr)$id, rep("a", nrow(pks)))
    expect_error(chromPeakData(xchr) <- 4)
    expect_error(chromPeakData(xchr) <- DataFrame())
    expect_error(chromPeakData(xchr) <- DataFrame(id = letters[1:nrow(pks)]))

    expect_true(nrow(chromPeaks(xchr, rt = c(20, 30))) == 0)
    expect_equal(chromPeaks(xchr, rt = c(2, 7)), pks)
    expect_equal(chromPeaks(xchr, rt = c(2, 7), type = "apex_within"),
                 pks[1:2, ])
    expect_equal(chromPeaks(xchr, rt = c(2, 7), type = "within"),
                 pks[2, , drop = FALSE])

    expect_equal(chromPeaks(xchr, mz = 123), pks)

    ## with m/z
    pks <- cbind(pks, mz = c(123, 332, 332, 432))
    pks <- cbind(pks, mzmin = c(122.9, 331.9, 331.8, 431.9))
    pks <- cbind(pks, mzmax = c(123.1, 332.1, 332.1, 432.2))
    chromPeaks(xchr) <- pks
    expect_true(nrow(chromPeaks(xchr, mz = 23)) == 0)
    expect_equal(chromPeaks(xchr, mz = 123, type = "apex_within"),
                 pks[1, , drop = FALSE])
    expect_equal(chromPeaks(xchr, mz = 331.89, ppm = 100), pks[2:3, ])
    expect_equal(chromPeaks(xchr, mz = 331.89, ppm = 10), pks[3, , drop = FALSE])

    ## with msLevel
    res <- chromPeaks(xchr, msLevel = 2L)
    expect_true(nrow(res) == 0)
    chromPeakData(xchr)$ms_level <- c(1L, 2L, 3L, 2L)
    res <- chromPeaks(xchr, msLevel = 2L)
    expect_equal(res, pks[c(2, 4), ])
})

test_that("plot,XChromatogram works", {
    skip_on_os(os = "windows", arch = "i386")

    chr <- Chromatogram(rtime = 1:10,
                        intensity = c(4, 12, 18, 24, 23, 18, 15, 3, 2, 5))
    xchr <- as(chr, "XChromatogram")
    pks <- matrix(nrow = 4, ncol = 6)
    colnames(pks) <- .CHROMPEAKS_REQ_NAMES
    pks[1, ] <- c(4, 2, 8, 24, 24, NA)
    pks[2, ] <- c(3, 2, 7, 24, 18, NA)
    pks[3, ] <- c(9, 7, 10, 2, 2, NA)
    pks[4, ] <- c(8, 5, 10, 3, 3, NA)
    chromPeaks(xchr) <- pks
    plot(xchr)
    plot(xchr, peakType = "point")
    plot(xchr, peakType = "rectangle", col = "red", lwd = 3)
    plot(xchr, peakType = "polygon", col = "red",
         peakCol = c("#00000020", "#ff000020", "#00ff0020", "#0000ff20"))
    plot(xchr, peakType = "polygon", col = "red",
         peakBg = c("#00000020", "#ff000020", "#00ff0020", "#0000ff20"))
    ## highlight the 3rd peak with an apex at 9
    plot(xchr, peakType = "point", peakCol = c("#00000040", "#00000040",
                                               "#ff000060", "#00000040"))
    plot(xchr, peakType = "rectangle", peakCol = c("#00000040", "#00000040",
                                                   "#ff000060", "#00000040"))
    plot(xchr, peakType = "polygon", peakCol = c("#00000040", "#00000040",
                                                 "#ff000060", "#00000040"))
})

test_that("filterMz,filterRt,XChromatogram work", {
    skip_on_os(os = "windows", arch = "i386")

    chr <- Chromatogram(rtime = 1:10,
                        intensity = c(4, 12, 18, 24, 23, 18, 15, 3, 2, 5))
    xchr <- as(chr, "XChromatogram")
    pks <- matrix(nrow = 4, ncol = 6)
    colnames(pks) <- .CHROMPEAKS_REQ_NAMES
    pks[1, ] <- c(4, 2, 8, 24, 24, NA)
    pks[2, ] <- c(3, 2, 7, 24, 18, NA)
    pks[3, ] <- c(9, 7, 10, 2, 2, NA)
    pks[4, ] <- c(8, 5, 10, 3, 3, NA)
    chromPeaks(xchr) <- pks

    expect_equal(filterRt(xchr), xchr)
    res <- filterRt(xchr, rt = c(3, 7))
    expect_equal(rtime(res), 3:7)
    expect_equal(intensity(res), intensity(xchr)[3:7])
    expect_equal(chromPeaks(res), pks[1:2, ])

    pks <- cbind(pks, mz = c(123, 124, 232, 234))
    chromPeaks(xchr) <- pks
    expect_equal(xchr, filterMz(xchr))
    res <- filterMz(xchr, mz = c(2, 3))
    expect_true(nrow(chromPeaks(res)) == 0)
    res <- filterMz(xchr, mz = c(200, 233))
    expect_equal(chromPeaks(res), pks[3, , drop = FALSE])
})

test_that("hasChromPeaks,XChromatogram works", {
    skip_on_os(os = "windows", arch = "i386")

    chr <- Chromatogram(rtime = 1:10,
                        intensity = c(4, 12, 18, 24, 23, 18, 15, 3, 2, 5))
    xchr <- as(chr, "XChromatogram")
    pks <- matrix(nrow = 4, ncol = 6)
    colnames(pks) <- .CHROMPEAKS_REQ_NAMES
    pks[1, ] <- c(4, 2, 8, 24, 24, NA)
    pks[2, ] <- c(3, 2, 7, 24, 18, NA)
    pks[3, ] <- c(9, 7, 10, 2, 2, NA)
    pks[4, ] <- c(8, 5, 10, 3, 3, NA)
    chromPeaks(xchr) <- pks
    expect_true(hasChromPeaks(xchr))

    xchr <- as(chr, "XChromatogram")
    expect_false(hasChromPeaks(xchr))
})

test_that("removeIntensity,XChromatogram(s) works", {
    skip_on_os(os = "windows", arch = "i386")

    set.seed(123)
    chr1 <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
                         intensity = c(5, 29, 50, NA, 100, 12, 3, 4, 1, 3))
    xchr <- as(chr1, "XChromatogram")
    expect_warning(res <- removeIntensity(xchr, which = "outside_chromPeak"))
    expect_equal(res, xchr)

    xchrs <- featureChromatograms(xod_xgrg, features = c("FT07", "FT13"),
                                  expandRt = 20)
    res <- removeIntensity(xchrs[[1]], which = "outside_chromPeak")
    expect_true(all(is.na(intensity(res)[rtime(res) < chromPeaks(res)[, "rtmin"]])))
    expect_true(all(is.na(intensity(res)[rtime(res) > chromPeaks(res)[, "rtmax"]])))
    expect_true(!all(is.na(intensity(xchrs[1, 1])[rtime(res) < chromPeaks(res)[, "rtmin"]])))
    expect_true(!all(is.na(intensity(xchrs[1, 1])[rtime(res) > chromPeaks(res)[, "rtmax"]])))

    res <- removeIntensity(xchrs, which = "outside_chromPeak")
    expect_true(all(is.na(intensity(res[1, 1])[rtime(res[1, 1]) < chromPeaks(res[1, 1])[, "rtmin"]])))
    expect_true(all(is.na(intensity(res[1, 1])[rtime(res[1, 1]) > chromPeaks(res[1, 1])[, "rtmax"]])))
    expect_true(all(is.na(intensity(res[2, 2])[rtime(res[2, 2]) < chromPeaks(res[2, 2])[, "rtmin"]])))
    expect_true(all(is.na(intensity(res[2, 2])[rtime(res[2, 2]) > chromPeaks(res[2, 2])[, "rtmax"]])))
    expect_true(all(is.na(intensity(res[2, 3])[rtime(res[2, 3]) < chromPeaks(res[2, 3])[, "rtmin"]])))
    expect_true(all(is.na(intensity(res[2, 3])[rtime(res[2, 3]) > chromPeaks(res[2, 3])[, "rtmax"]])))

    ## with two peaks per chromatogram.
    mzr <- 462.2 + c(-0.04, 0.04)
    chr <- chromatogram(od_x, mz = mzr)
    chr <- findChromPeaks(chr, CentWaveParam())
    res <- removeIntensity(chr, which = "outside_chromPeak")

    ## And if we had an XChromatograms
    chr2 <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
                         intensity = c(80, 50, 20, 10, 9, 4, 3, 4, 1, 3))
    chr3 <- Chromatogram(rtime = 3:9 + rnorm(7, sd = 0.3),
                         intensity = c(53, 80, 130, 15, 5, 3, 2))
    chrs <- MChromatograms(list(chr1, chr2, chr3))
    chrs <- as(chrs, "XChromatograms")

    res <- removeIntensity(chrs)
    expect_true(is(res, "XChromatograms"))
    expect_equal(res, chrs)

    res <- removeIntensity(chrs, threshold = 20)
    expect_equal(intensity(res[1, 1]), c(NA_real_, 29, 50, NA_real_, 100,
                                         NA_real_, NA_real_, NA_real_, NA_real_,
                                         NA_real_))
    expect_equal(intensity(res[3, 1]), c(53, 80, 130, NA_real_, NA_real_,
                                         NA_real_, NA_real_))
})

test_that("filterChromPeaks,XChromatogram works", {
    skip_on_os(os = "windows", arch = "i386")

    chrs <- chromatogram(xod_x, mz = rbind(305.1 + c(-0.01, 0.01),
                                           462.2 + c(-0.04, 0.04)))
    res <- filterChromPeaks(chrs[1, 1], n = 1L, by = "keepTop")
    expect_true(nrow(chromPeaks(res)) == 1L)
})

test_that("transformIntensity,XChromatogram works", {
    skip_on_os(os = "windows", arch = "i386")

    chrs <- chromatogram(xod_x, mz = rbind(305.1 + c(-0.01, 0.01),
                                           462.2 + c(-0.04, 0.04)))

    chr <- chrs[1L, 1L]
    res <- transformIntensity(chr)

    expect_equal(intensity(res), intensity(chr))
    expect_equal(chromPeaks(res)[, "into"], chromPeaks(chr)[, "into"])
    expect_equal(chromPeaks(res)[, "maxo"], chromPeaks(chr)[, "maxo"])

    res <- transformIntensity(chr, log2)
    expect_equal(intensity(res), log2(intensity(chr)))
    expect_equal(chromPeaks(res)[, "into"], log2(chromPeaks(chr)[, "into"]))
    expect_equal(chromPeaks(res)[, "maxo"], log2(chromPeaks(chr)[, "maxo"]))

    chr <- XChromatogram()
    res <- transformIntensity(chr)
    expect_equal(intensity(res), intensity(chr))
    expect_equal(chromPeaks(res)[, "into"], chromPeaks(chr)[, "into"])
    expect_equal(chromPeaks(res)[, "maxo"], chromPeaks(chr)[, "maxo"])

    res <- transformIntensity(chr, log2)
    expect_equal(intensity(res), log2(intensity(chr)))
    expect_equal(chromPeaks(res)[, "into"], log2(chromPeaks(chr)[, "into"]))
    expect_equal(chromPeaks(res)[, "maxo"], log2(chromPeaks(chr)[, "maxo"]))
})
