test_that("chromPeaks,XChromatograms works", {
    pks1 <- matrix(c(3, 2, 4, 339.2, 343, NA), nrow = 1,
                   dimnames = list(NULL, xcms:::.CHROMPEAKS_REQ_NAMES))
    chr1 <- XChromatogram(
        rtime = 1:8, intensity = c(3, 24.2, 343, 32, 3.3, 5, 2, 9),
        chromPeaks = pks1)
    chr2 <- XChromatogram(rtime = 1:4, intensity = c(45, 3, 34, 2))
    pks3 <- matrix(c(3, 2, 4, 145, 54, NA), nrow = 1,
                   dimnames = list(NULL, xcms:::.CHROMPEAKS_REQ_NAMES))
    chr3 <- XChromatogram(
        rtime = 1:7, intensity = c(12, 34, 54, 34, 23, 2, NA),
        chromPeaks = pks3)
    chr4 <- XChromatogram(rtime = 1:3, intensity = c(3, 4, 1))
    chr5 <- XChromatogram(rtime = 1:6, intensity = c(3, 4, 6, 7, 2, 4))
    pks6 <- matrix(c(2, 2, 3, 108, 65, NA, 3, 5, 7, 123, 4, NA),
                   nrow = 2, byrow = TRUE,
                   dimnames = list(NULL, xcms:::.CHROMPEAKS_REQ_NAMES))
    chr6 <- XChromatogram(
        rtime = 2:5, intensity = c(3, 65, 43, 12),
        chromPeaks = pks6)
    xchrs <- XChromatograms(list(chr1, chr2, chr3, chr4, chr5, chr6), nrow = 2)

    res <- chromPeaks(xchrs)
    expect_equal(res[1, 1:6, drop = FALSE], pks1)
    expect_equal(res[2, 1:6, drop = FALSE], pks3)
    expect_equal(res[3:4, 1:6, drop = FALSE], pks6)
    expect_equal(res[, "row"], c(1, 1, 2, 2))
    expect_equal(res[, "column"], c(1, 2, 3, 3))

    xchrs <- XChromatograms(list(chr4, chr5, chr2))
    res <- chromPeaks(xchrs)
    expect_true(nrow(res) == 0)
    expect_equal(colnames(res), c("rt", "rtmin", "rtmax", "into", "maxo", "sn",
                                  "row", "column"))

    xchrs <- XChromatograms(list(chr2, chr4, chr5, chr1))
    res <- chromPeaks(xchrs)
    expect_equal(res[, 1:6, drop = FALSE], pks1)
    expect_equal(unname(res[, "row"]), c(4))
    expect_equal(unname(res[, "column"]), c(1))

    xchrs <- XChromatograms(list(chr2, chr4, chr5, chr1), ncol = 2)
    res <- chromPeaks(xchrs)
    expect_equal(res[, 1:6, drop = FALSE], pks1)
    expect_equal(unname(res[, "row"]), c(2))
    expect_equal(unname(res[, "column"]), c(2))
})

test_that("filterRt, filterMz for XChromatograms works", {
    chr1 <- XChromatogram()
    chr2 <- XChromatogram(rtime = 1:4, intensity = c(45, 3, 34, 2))
    pks3 <- matrix(c(3, 2, 4, 145, 54, NA), nrow = 1,
                   dimnames = list(NULL, xcms:::.CHROMPEAKS_REQ_NAMES))
    chr3 <- XChromatogram(
        rtime = 1:7, intensity = c(12, 34, 54, 34, 23, 2, NA),
        chromPeaks = pks3)
    chr4 <- XChromatogram(rtime = 1:3, intensity = c(3, 4, 1))

    xchrs <- XChromatograms()
    expect_equal(xchrs, filterRt(xchrs))
    expect_equal(xchrs, filterMz(xchrs))

    xchrs <- XChromatograms(list(chr1, chr2))
    res <- filterMz(xchrs)
    expect_equal(xchrs, res)
    res <- filterRt(xchrs, rt = c(3, 4))
    expect_equal(res[1, 1], chr1)
    expect_equal(rtime(res[2, 1]), c(3, 4))
    expect_equal(intensity(res[2, 1]), c(34, 2))

    xchrs <- XChromatograms(list(chr3, chr4), nrow = 1)
    res <- filterMz(xchrs)
    expect_equal(xchrs, res)
    res <- filterRt(xchrs, rt = c(6, 7))
    expect_true(length(rtime(res[1, 2])) == 0)
    expect_equal(rtime(res[1, 1]), 6:7)
    expect_equal(intensity(res[1, 1]), c(2, NA))
    expect_true(nrow(chromPeaks(res)) == 0)

    res <- filterRt(xchrs, rt = c(3, 6))
    expect_equal(rtime(res[1, 1]), 3:6)
    expect_equal(intensity(res[1, 1]), c(54, 34, 23, 2))
    expect_true(nrow(chromPeaks(res)) == 1)
    expect_equal(chromPeaks(res[1, 1]), pks3)
})

test_that("plot,XChromatogram works", {
    mzr <- matrix(c(335, 335, 344, 344), ncol = 2, byrow = TRUE)
    rtr <- matrix(c(2700, 2900, 2600, 2750), ncol = 2, byrow = TRUE)

    ## Full rt range.
    xchr_rt <- chromatogram(xod_chr, mz = mzr)

    cols_smple <- c("#ff000060", "#00ff0060", "#0000ff60")
    plot(as(xchr_rt[1, ], "Chromatograms"), col = cols_smple)
    ## Plotting the data is fine, just as above.
    ## Then we have to loop over each chromatogram...
    x <- xchr_rt[1, ]
    pks <- chromPeaks(x)
    .add_chromatogram_peaks(x, pks, type = "rectangle", bg = rep(NA, nrow(pks)),
                            col = cols_smple[pks[, "sample"]])
    .add_chromatogram_peaks(x, pks, type = "point", bg = rep(NA, nrow(pks)),
                            col = cols_smple[pks[, "sample"]], pch = 15)
    .add_chromatogram_peaks(x, pks, type = "polygon",
                            bg = cols_smple[pks[, "sample"]],
                            col = cols_smple[pks[, "sample"]])
    x <- xchr_rt[2, ]
    pks <- chromPeaks(x)
    plot(as(x, "Chromatograms"), col = cols_smple, lwd = 2)
    .add_chromatogram_peaks(x, pks, type = "rectangle", bg = rep(NA, nrow(pks)),
                            col = cols_smple[pks[, "sample"]])
    .add_chromatogram_peaks(x, pks, type = "point", bg = rep(NA, nrow(pks)),
                            col = cols_smple[pks[, "sample"]], pch = 15)
    .add_chromatogram_peaks(x, pks, type = "polygon",
                            bg = cols_smple[pks[, "sample"]],
                            col = cols_smple[pks[, "sample"]])

    plot(xchr_rt, peakCol = cols_smple[chromPeaks(xchr_rt)[, "sample"]],
         peakBg = cols_smple[chromPeaks(xchr_rt)[, "sample"]], xlab = "RT")
    xsub <- xchr_rt[2, ]
    ## Use one color per peak
    library(RColorBrewer)
    cls <- paste0(brewer.pal(nrow(chromPeaks(xsub)), "Dark2"), 40)
    plot(xsub, peakBg = cls)

    ## Narrow on rt.
    xchr <- chromatogram(xod_chr, mz = mzr, rt = rtr)
    plot(xchr)
})

test_that("processHistory,XChromatograms works", {
    mzr <- matrix(c(335, 335, 344, 344), ncol = 2, byrow = TRUE)
    xchrs <- chromatogram(xod_xgrg, mz = mzr, rt = c(2600, 3600))
    res <- processHistory(xchrs)
    expect_equal(length(res), 4)
    res <- processHistory(xchrs, type = .PROCSTEP.PEAK.DETECTION)
    expect_equal(length(res), 1)

    xchrs2 <- findChromPeaks(xchrs, param = CentWaveParam())
    expect_equal(length(processHistory(xchrs2)), 5)
    expect_equal(length(processHistory(xchrs2,
                                       type = .PROCSTEP.PEAK.DETECTION)), 2)
})

test_that("groupChromPeaks,XChromatograms,PeakDensityParam works", {
    mzr <- matrix(c(335, 335, 344, 344), ncol = 2, byrow = TRUE)
    xchrs <- chromatogram(xod_xgrg, mz = mzr, rt = c(2600, 3600))
    param <- PeakDensityParam(sampleGroups = c(1, 1, 1))
    res <- groupChromPeaks(xchrs, param = param)
    expect_true(hasFeatures(res))
    expect_true(nrow(res@featureDefinitions) == 1)

    param <- PeakDensityParam(sampleGroups = c(1, 2, 3))
    res <- groupChromPeaks(xchrs, param = param)
    expect_true(nrow(res@featureDefinitions) == 2)
    expect_true(length(processHistory(res)) == 5)

    res <- groupChromPeaks(res,
                           param = PeakDensityParam(sampleGroups = c(1, 1, 1)))
    expect_true(length(processHistory(res)) == 5)
    expect_true(nrow(res@featureDefinitions) == 1)

    ## The same on artificial data.
    chrs <- as(xchrs, "Chromatograms")
    res <- findChromPeaks(chrs, param = CentWaveParam(snthresh = 1))
    expect_equal(nrow(chromPeaks(res)), 31)
    res <- groupChromPeaks(res, param = param)
    expect_equal(nrow(featureDefinitions(res)), 8)
})

test_that("dropFeatureDefinitions,XChromatograms works", {
    mzr <- matrix(c(335, 335, 344, 344), ncol = 2, byrow = TRUE)
    xchrs <- chromatogram(xod_xgrg, mz = mzr, rt = c(2600, 3600))
    param <- PeakDensityParam(sampleGroups = c(1, 1, 1))
    res <- groupChromPeaks(xchrs, param = param)
    expect_true(hasFeatures(res))
    expect_true(nrow(res@featureDefinitions) == 1)
    expect_true(length(res@.processHistory) == 5)
    res <- dropFeatureDefinitions(res)
    expect_false(hasFeatures(res))
    expect_true(length(res@.processHistory) == 4)
    expect_equal(res, xchrs)
})

test_that("featureDefinitions,XChromatograms works", {
    mzr <- matrix(c(335, 335, 344, 344), ncol = 2, byrow = TRUE)
    xchrs <- chromatogram(xod_xgrg, mz = mzr, rt = c(2600, 3600))
    param <- PeakDensityParam(sampleGroups = c(1, 1, 1))
    res <- groupChromPeaks(xchrs, param = param)
    expect_true(hasFeatures(res))
    expect_true(nrow(featureDefinitions(res)) == 1)

    xchrs <- findChromPeaks(xchrs, param = CentWaveParam())
    res <- groupChromPeaks(xchrs, param = param)
    expect_true(nrow(featureDefinitions(res)) == 3)
    fts <- featureDefinitions(res, rt = c(2500, 2800))
    expect_equal(nrow(fts), 2)

    res <- dropFeatureDefinitions(res)
    expect_equal(nrow(featureDefinitions(res)), 0)
})

test_that(".subset_chrom_peaks_xchromatograms works", {
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

test_that("[,XChromatograms works", {
    mzr <- matrix(c(335, 335, 344, 344), ncol = 2, byrow = TRUE)
    chrs <- chromatogram(od_x, mz = mzr)
    chrs <- findChromPeaks(chrs, param = CentWaveParam())
    prm <- PeakDensityParam(sampleGroups = c(1, 1, 1))
    chrs <- groupChromPeaks(chrs, param = prm)

    pks <- chromPeaks(chrs)
    fts <- featureDefinitions(chrs)

    res <- chrs[2, 2]
    expect_true(is(res, "XChromatogram"))
    expect_equal(chromPeaks(res), pks[pks[, "row"] == 2 &
                                      pks[, "column"] == 2,
                                      colnames(chromPeaks(res))])

    res <- chrs[2, 2:3]
    expect_true(is(res, "XChromatograms"))
    expect_true(ncol(res) == 2)
    expect_equal(res[1, 1], chrs[2, 2])
    expect_equal(res[1, 2], chrs[2, 3])
    pks_tmp <- pks[pks[, "row"] == 2 & pks[, "column"] %in% 2:3, ]
    pks_tmp <- pks_tmp[order(pks_tmp[, "column"], pks_tmp[, "row"]), ]
    expect_equal(chromPeaks(res)[, 1:6], pks_tmp[, 1:6])
    expect_equal(rownames(featureDefinitions(res)), c("FT3", "FT4"))
    res_fts <- featureDefinitions(res)
    expect_equal(res_fts$peakidx, list(c(1, 3), c(2, 4)))

    res <- chrs[c(2, 1), 1]
    expect_true(is(res, "XChromatograms"))
    expect_true(ncol(res) == 1)
    expect_true(nrow(res) == 2)
    expect_equal(res[1, 1], chrs[2, 1])
    expect_equal(res[2, 1], chrs[1, 1])
    expect_equal(res$sampleNames, factor("ko15.CDF"))
    expect_equal(fData(res), fData(chrs)[c(2, 1), ])
    pks_tmp <- pks[pks[, "row"] %in% c(1, 2) & pks[, "column"] == 1, ]
    pks_tmp <- pks_tmp[c(5, 6, 1, 2, 3, 4), ]
    expect_equal(chromPeaks(res)[, 1:6], pks_tmp[, 1:6])
    res_fts <- featureDefinitions(res)
    expect_equal(rownames(res_fts), c("FT3", "FT4", "FT1", "FT2"))
    expect_equal(res_fts$peakidx, list(c(1), c(2), c(3), c(6)))

    res <- chrs[c(2, 1), c(1, 3)]
    expect_true(is(res, "XChromatograms"))
    expect_true(ncol(res) == 2)
    expect_true(nrow(res) == 2)
    expect_equal(res[1, 1], chrs[2, 1])
    expect_equal(res[2, 1], chrs[1, 1])
    expect_equal(res[1, 2], chrs[2, 3])
    expect_equal(res[2, 2], chrs[1, 3])
    expect_equal(res$sampleNames, factor(c("ko15.CDF", "ko18.CDF")))
    expect_equal(fData(res), fData(chrs)[c(2, 1), ])
    pks_tmp <- pks[c(8, 9, 12, 13, 1, 2, 3, 4, 6, 7), ]
    expect_equal(chromPeaks(res)[, 1:6], pks_tmp[, 1:6])
    res_fts <- featureDefinitions(res)
    expect_equal(rownames(res_fts), c("FT3", "FT4", "FT1", "FT2"))
    expect_equal(res_fts$peakidx, list(c(1, 3), c(2, 4), c(5, 9), c(8, 10)))
})

test_that("featureValues,XChromatograms works", {
    mzr <- matrix(c(335, 335, 344, 344), ncol = 2, byrow = TRUE)
    chrs <- chromatogram(od_x, mz = mzr)
    chrs <- as(chrs, "XChromatograms")
    expect_error(featureValues(chrs))
    chrs <- findChromPeaks(chrs, param = CentWaveParam())
    expect_error(featureValues(chrs))
    prm <- PeakDensityParam(sampleGroups = c(1, 1, 1))
    chrs <- groupChromPeaks(chrs, param = prm)

    vls <- featureValues(chrs)
    expect_equal(colnames(vls), colnames(chrs))
    expect_equal(rownames(vls), rownames(featureDefinitions(chrs)))
    exp_mat <- matrix(c(1, 5, 6,
                        4, NA, 7,
                        8, 10, 12,
                        9, 11, 13), byrow = TRUE, ncol = 3,
                      dimnames = list(rownames(featureDefinitions(chrs)),
                                      colnames(chrs)))
    expect_equal(exp_mat, vls)
    ## into.
    vls <- featureValues(chrs, value = "into")
    vls_exp <- matrix(chromPeaks(chrs)[exp_mat, "into"], ncol = 3, byrow = FALSE,
                      dimnames = dimnames(exp_mat))
    expect_equal(vls, vls_exp)

    vls <- featureValues(chrs, value = "into", missing = 13)
    vls_exp[is.na(vls_exp)] <- 13
    expect_equal(vls, vls_exp)

    ## After subsetting/re-ordering.
    chrs_sub <- chrs[2, c(3, 1, 2)]
    vls_sub <- featureValues(chrs_sub, value = "into")
    vls <- featureValues(chrs, value = "into")
    expect_equal(colnames(vls_sub), colnames(chrs_sub))
    expect_equal(rownames(vls_sub), rownames(featureDefinitions(chrs_sub)))
    expect_equal(vls_sub, vls[3:4, c(3, 1, 2)])

    chrs_sub <- chrs[c(2, 1), 2]
    vls_sub <- featureValues(chrs_sub, value = "into")
    expect_equal(colnames(vls_sub), colnames(chrs_sub))
    expect_equal(rownames(vls_sub), rownames(featureDefinitions(chrs_sub)))
    expect_equal(vls_sub, vls[c(3, 4, 1), 2, drop = FALSE])
})
