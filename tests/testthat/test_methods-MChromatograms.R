test_that("findChromPeaks,MChromatograms works", {
    skip_on_os(os = "windows", arch = "i386")

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

test_that("correlate,MChromatograms works", {
    skip_on_os(os = "windows", arch = "i386")

    set.seed(123)
    chr1 <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
                         intensity = c(5, 29, 50, NA, 100, 12, 3, 4, 1, 3))
    chr2 <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
                         intensity = c(80, 50, 20, 10, 9, 4, 3, 4, 1, 3))
    chr3 <- Chromatogram(rtime = 3:9 + rnorm(7, sd = 0.3),
                         intensity = c(53, 80, 130, 15, 5, 3, 2))
    chrs <- MChromatograms(list(chr1, chr2, chr3))

    res <- correlate(chrs)
    expect_true(nrow(res) == 3)
    expect_true(ncol(res) == 3)
    expect_true(res[1, 3] > 0.9)
    expect_true(res[1, 2] < 0.5)

    res_2 <- correlate(chrs, chrs)
    expect_equal(res_2, res)

    res <- correlate(chrs)
    expect_equal(res[2, 1], res[1, 2])
    expect_equal(res[3, 1], res[1, 3])

})

test_that("removeIntensity,MChromatograms works", {
    skip_on_os(os = "windows", arch = "i386")

    set.seed(123)
    chr1 <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
                         intensity = c(5, 29, 50, NA, 100, 12, 3, 4, 1, 3))
    chr2 <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
                         intensity = c(80, 50, 20, 10, 9, 4, 3, 4, 1, 3))
    chr3 <- Chromatogram(rtime = 3:9 + rnorm(7, sd = 0.3),
                         intensity = c(53, 80, 130, 15, 5, 3, 2))
    chrs <- MChromatograms(list(chr1, chr2, chr3))

    res <- removeIntensity(chrs)
    expect_equal(res, chrs)

    res <- removeIntensity(chrs, threshold = 20)
    expect_equal(intensity(res[1, 1]), c(NA_real_, 29, 50, NA_real_, 100,
                                         NA_real_, NA_real_, NA_real_, NA_real_,
                                         NA_real_))
    expect_equal(intensity(res[3, 1]), c(53, 80, 130, NA_real_, NA_real_,
                                         NA_real_, NA_real_))

    chrs <- MChromatograms(list(chr1, chr2, chr2, chr3), ncol = 2)
    res <- removeIntensity(chrs, threshold = 20)
    expect_equal(intensity(res[2, 2]), c(53, 80, 130, NA_real_, NA_real_,
                                         NA_real_, NA_real_))
})

test_that("filterColumnsIntensityAbove,MChromatograms works", {
    skip_on_os(os = "windows", arch = "i386")

    set.seed(123)
    chr1 <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
                         intensity = c(5, 29, 50, NA, 100, 12, 3, 4, 1, 3))
    chr2 <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
                         intensity = c(80, 50, 20, 10, 9, 4, 3, 4, 1, 3))
    chr3 <- Chromatogram(rtime = 3:9 + rnorm(7, sd = 0.3),
                         intensity = c(53, 80, 130, 15, 5, 3, 2))
    chrs <- MChromatograms(list(chr1, chr2, chr3, chr1, chr2, chr3), ncol = 3)

    expect_error(filterColumnsIntensityAbove(chrs, threshold = c(1.1, 1.4)),
                 "should be")
    expect_error(filterColumnsIntensityAbove(chrs, threshold = TRUE),
                 "should be")

    res <- filterColumnsIntensityAbove(chrs)
    expect_equal(res, chrs)

    res <- filterColumnsIntensityAbove(chrs, threshold = 90)
    expect_equal(res, chrs)

    res <- filterColumnsIntensityAbove(chrs, threshold = 90, which = "all")
    expect_equal(res, chrs[, 2])

    res <- filterColumnsIntensityAbove(chrs, threshold = 200, which = "any",
                                       value = "tic")
    expect_equal(res, chrs)

    res <- filterColumnsIntensityAbove(chrs, threshold = 200, which = "all",
                                       value = "tic")
    expect_equal(res, chrs[, 2])
})

test_that("filterChromatogramsKeepTop,MChromatograms works", {
    skip_on_os(os = "windows", arch = "i386")

    set.seed(123)
    chr1 <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
                         intensity = c(5, 29, 50, NA, 100, 12, 3, 4, 1, 3))
    chr2 <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
                         intensity = c(80, 50, 20, 10, 9, 4, 3, 4, 1, 3))
    chr3 <- Chromatogram(rtime = 3:9 + rnorm(7, sd = 0.3),
                         intensity = c(53, 80, 130, 15, 5, 3, 2))
    chr4 <- Chromatogram(rtime = 1:10,
                         intensity = c(NA, NA, 4, NA, NA, 9, NA, 10, 9, 1))
    chr5 <- Chromatogram(rtime = 1:4, intensity = c(345, 5554, 323, 2000))
    chr6 <- Chromatogram(rtime = 1:3, intensity = c(400, 244, 133))


    chrs <- MChromatograms(list(chr1, chr2, chr3, chr4, chr5, chr6), ncol = 3)

    expect_error(filterColumnsKeepTop(chrs, n = c(1, 2)), "of length 1")
    expect_error(filterColumnsKeepTop(chrs, n = "b"), "of length 1")
    expect_error(filterColumnsKeepTop(chrs, n = 10), "number of columns")

    res <- filterColumnsKeepTop(chrs, n = 1)
    expect_equal(res, chrs[, 3])

    res <- filterColumnsKeepTop(chrs, n = 2)
    expect_equal(res[, 1], chrs[, 1])
    expect_equal(res[, 2], chrs[, 3])

    res <- filterColumnsKeepTop(chrs, n = 2, aggregationFun = max)
    expect_equal(res[, 1], chrs[, 2])
    expect_equal(res[, 2], chrs[, 3])

    res <- filterColumnsKeepTop(chrs, n = 0)
    expect_true(ncol(res) == 0)
    expect_true(nrow(res) == 2)

    res <- filterColumnsKeepTop(chrs, n = 1, sortBy = "tic")
    expect_equal(res, chrs[, 3])

    res <- filterColumnsKeepTop(chrs, n = 1, aggregationFun = mean)
    expect_equal(res, chrs[, 3])
})

test_that("normalize,MChromatograms works", {
    skip_on_os(os = "windows", arch = "i386")

    set.seed(123)
    chr1 <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
                         intensity = c(5, 29, 50, NA, 100, 12, 3, 4, 1, 3))
    chr2 <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
                         intensity = c(80, 50, 20, 10, 9, 4, 3, 4, 1, 3))
    chr3 <- Chromatogram(rtime = 3:9 + rnorm(7, sd = 0.3),
                         intensity = c(53, 80, 130, 15, 5, 3, 2))
    chr4 <- Chromatogram(rtime = 1:10,
                         intensity = c(NA, NA, 4, NA, NA, 9, NA, 10, 9, 1))
    chrs <- MChromatograms(list(chr1, chr2, chr3, chr4), ncol = 2)
    res <- normalize(chrs)

    expect_true(ncol(res) == ncol(chrs))
    expect_true(nrow(res) == nrow(chrs))

    expect_equal(intensity(res[1, 2]) * max(intensity(chrs[1, 2]), na.rm = TRUE),
                 intensity(chrs[1, 2]))
})

test_that(".plot_xchromatograms_overlay works", {
    skip_on_os(os = "windows", arch = "i386")

    set.seed(123)
    chr1 <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
                         intensity = c(5, 29, 50, NA, 100, 12, 3, 4, 1, 3))
    chr2 <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
                         intensity = c(80, 50, 20, 10, 9, 4, 3, 4, 1, 3))
    chrs <- MChromatograms(list(chr1, chr2), ncol = 1)
    .plot_xchromatograms_overlay(chrs)
    .plot_xchromatograms_overlay(chrs, xlim = c(-10, 20),
                                 ylim = c(0, 150), fill = "red")
    .plot_xchromatograms_overlay(chrs, xlim = c(-10, 20), ylim = c(0, 150),
                                 yoffset = 10, fill = c("red", "blue"))

})

test_that("plotChromatogramsOverlay,MChromatograms,XChromatograms work", {
    skip_on_os(os = "windows", arch = "i386")

    data(xdata)
    dirname(xdata) <- c(rep(system.file("cdf", "KO", package = "faahKO"), 4),
                        rep(system.file("cdf", "WT", package = "faahKO"), 4))
    fts <- c("FT097", "FT163", "FT165")
    xdata <- filterFile(xdata, file = 1:2, keepFeatures = TRUE)
    chrs <- featureChromatograms(xdata, features = fts)

    plotChromatogramsOverlay(chrs)
    plotChromatogramsOverlay(chrs, transform = log10)
    plotChromatogramsOverlay(chrs, peakType = "rectangle", peakBg = NA)
    plotChromatogramsOverlay(chrs, peakType = "rectangle", peakBg = NA,
                             transform = log2)
    plotChromatogramsOverlay(
        chrs, peakType = "rectangle", peakBg = NA, yoffset = 100000,
        fill = c("#ff000040", "#00ff0040", "#0000ff40"))

    res <- plotChromatogramsOverlay(chrs, stacked = 0.5, bty = "n")
    expect_equal(length(res), ncol(chrs))
    res <- plotChromatogramsOverlay(chrs, stacked = 0.5, bty = "n",
                                    transform = log2)
    res <- plotChromatogramsOverlay(chrs, stacked = 0.1, bty = "n")

    plotChromatogramsOverlay(chrs[1, ])

    chr <- chrs[, 1]
    plotChromatogramsOverlay(chr, peakBg = c("red", "blue"))
    plotChromatogramsOverlay(chr, peakBg = c("blue", "red"))

    chrs <- as(chrs, "MChromatograms")
    plotChromatogramsOverlay(chrs)
    plotChromatogramsOverlay(chrs, transform = log2)
})
