test_that("chromPeaks, chromPeakData for XChromatograms work", {
    skip_on_os(os = "windows", arch = "i386")

    pks1 <- matrix(c(3, 2, 4, 339.2, 343, NA), nrow = 1,
                   dimnames = list(NULL, .CHROMPEAKS_REQ_NAMES))
    chr1 <- XChromatogram(
        rtime = 1:8, intensity = c(3, 24.2, 343, 32, 3.3, 5, 2, 9),
        chromPeaks = pks1)
    chr2 <- XChromatogram(rtime = 1:4, intensity = c(45, 3, 34, 2))
    pks3 <- matrix(c(3, 2, 4, 145, 54, NA), nrow = 1,
                   dimnames = list(NULL, .CHROMPEAKS_REQ_NAMES))
    chr3 <- XChromatogram(
        rtime = 1:7, intensity = c(12, 34, 54, 34, 23, 2, NA),
        chromPeaks = pks3)
    chr4 <- XChromatogram(rtime = 1:3, intensity = c(3, 4, 1))
    chr5 <- XChromatogram(rtime = 1:6, intensity = c(3, 4, 6, 7, 2, 4))
    pks6 <- matrix(c(2, 2, 3, 108, 65, NA, 3, 5, 7, 123, 4, NA),
                   nrow = 2, byrow = TRUE,
                   dimnames = list(NULL, .CHROMPEAKS_REQ_NAMES))
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
    resd <- chromPeakData(xchrs)
    expect_equal(resd$row, c(1, 1, 2, 2))
    expect_equal(resd$column, c(1, 2, 3, 3))
    expect_equal(colnames(resd), c("ms_level", "is_filled", "row", "column"))

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
    resd <- chromPeakData(xchrs)
    expect_equal(resd$column, 1)

    xchrs <- XChromatograms(list(chr2, chr4, chr5, chr1), ncol = 2)
    res <- chromPeaks(xchrs)
    expect_equal(res[, 1:6, drop = FALSE], pks1)
    expect_equal(unname(res[, "row"]), c(2))
    expect_equal(unname(res[, "column"]), c(2))
})

test_that("filterRt, filterMz for XChromatograms works", {
    skip_on_os(os = "windows", arch = "i386")

    chr1 <- XChromatogram()
    chr2 <- XChromatogram(rtime = 1:4, intensity = c(45, 3, 34, 2))
    pks3 <- matrix(c(3, 2, 4, 145, 54, NA), nrow = 1,
                   dimnames = list(NULL, .CHROMPEAKS_REQ_NAMES))
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

    chrs <- as(od_chrs, "XChromatograms")
    chrs <- findChromPeaks(chrs, param = CentWaveParam())
    prm <- PeakDensityParam(sampleGroups = c(1, 1, 1))
    chrs <- groupChromPeaks(chrs, param = prm)
    ## Filter on rt on the one above
    rtr <- c(2500, 3000)
    res <- filterRt(chrs, rt = rtr)
    pks_all <- chromPeaks(chrs)
    excl <- !(pks_all[, "rtmin"] < rtr[2] & pks_all[, "rtmax"] > rtr[1])
    pks <- chromPeaks(res)
    expect_true(all(pks[, "rtmin"] < rtr[2] & pks[, "rtmax"] > rtr[1]))
    expect_equal(pks, chromPeaks(chrs)[!excl, ])
    expect_equal(rownames(featureDefinitions(res)), c("FT1", "FT3", "FT4"))
    expect_equal(featureDefinitions(res)$peakidx,
                 list(c(1, 4, 5), c(6, 8, 10), c(7, 9, 11)))
    rtr <- c(2500, 2600)
    res <- filterRt(chrs, rtr)
    expect_equal(featureDefinitions(res)$peakidx, list(1:3))
    expect_equal(rownames(featureDefinitions(res)), "FT3")

    ## Filter on mz for a chrs extracted from a real object.
    mzr <- matrix(c(335, 335, 344, 344), ncol = 2, byrow = TRUE)
    chrs <- chromatogram(xod_xgrg, mz = mzr)
    expect_equal(chromPeakData(chrs)$ms_level, rep(1L, 4))
    expect_equal(chromPeakData(chrs)$is_filled, rep(FALSE, 4))
    expect_equal(chromPeakData(chrs)$row, c(1, 2, 2, 2))
    expect_equal(chromPeakData(chrs)$column, c(2, 1, 2, 3))

    res <- filterMz(chrs, mz = 335)
    expect_equal(nrow(featureDefinitions(res)), 0)
    expect_equal(nrow(chromPeaks(res)), 1)
    expect_equal(nrow(chromPeakData(res)), 1)

    res <- filterMz(chrs, mz = 344)
    expect_equal(nrow(featureDefinitions(res)), 1)
    expect_equal(nrow(chromPeaks(res)), 3)
    expect_equal(nrow(chromPeakData(res)), 3)

    res <- filterMz(chrs, mz = 444)
    expect_equal(nrow(featureDefinitions(res)), 0)
    expect_equal(nrow(chromPeaks(res)), 0)
    expect_equal(nrow(chromPeakData(res)), 0)
})

test_that("plot,XChromatogram works", {
    skip_on_os(os = "windows", arch = "i386")

    mzr <- matrix(c(335, 335, 344, 344), ncol = 2, byrow = TRUE)
    rtr <- matrix(c(2700, 2900, 2600, 2750), ncol = 2, byrow = TRUE)

    ## Full rt range.
    xchr_rt <- chromatogram(xod_chr, mz = mzr)

    cols_smple <- c("#ff000060", "#00ff0060", "#0000ff60")
    plot(as(xchr_rt[1, ], "MChromatograms"), col = cols_smple)
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
    plot(as(x, "MChromatograms"), col = cols_smple, lwd = 2)
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
    skip_on_os(os = "windows", arch = "i386")

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
    skip_on_os(os = "windows", arch = "i386")

    mzr <- matrix(c(335, 335, 344, 344), ncol = 2, byrow = TRUE)
    xchrs <- chromatogram(xod_xgrg, mz = mzr, rt = c(2600, 3600))
    param <- PeakDensityParam(sampleGroups = c(1, 1, 1))
    res <- groupChromPeaks(xchrs, param = param)
    expect_true(hasFeatures(res))
    expect_true(nrow(res@featureDefinitions) == 1)
    expect_equal(processHistory(xchrs)[1:3], processHistory(res)[1:3])
    expect_true(processHistory(xchrs)[[4]]@date !=
                processHistory(res)[[4]]@date)

    param <- PeakDensityParam(sampleGroups = c(1, 2, 3))
    res <- groupChromPeaks(xchrs, param = param)
    expect_true(nrow(res@featureDefinitions) == 2)
    expect_true(length(processHistory(res)) == 4)
    expect_equal(processHistory(xchrs)[1:3], processHistory(res)[1:3])
    expect_true(processHistory(xchrs)[[4]]@date !=
                processHistory(res)[[4]]@date)

    res <- groupChromPeaks(res,
                           param = PeakDensityParam(sampleGroups = c(1, 1, 1)))
    expect_true(length(processHistory(res)) == 4)
    expect_true(nrow(res@featureDefinitions) == 1)

    ## The same on artificial data.
    chrs <- as(xchrs, "MChromatograms")
    res <- findChromPeaks(chrs, param = CentWaveParam(snthresh = 1))
    expect_equal(nrow(chromPeaks(res)), 31)
    res <- groupChromPeaks(res, param = param)
    expect_equal(nrow(featureDefinitions(res)), 8)
})

test_that("dropFeatureDefinitions,XChromatograms works", {
    skip_on_os(os = "windows", arch = "i386")

    mzr <- matrix(c(335, 335, 344, 344), ncol = 2, byrow = TRUE)
    xchrs <- chromatogram(xod_xgrg, mz = mzr, rt = c(2600, 3600))
    param <- PeakDensityParam(sampleGroups = c(1, 1, 1))
    res <- groupChromPeaks(xchrs, param = param)
    expect_true(hasFeatures(res))
    expect_true(nrow(res@featureDefinitions) == 1)
    expect_true(length(res@.processHistory) == 4)
    res <- dropFeatureDefinitions(res)
    expect_false(hasFeatures(res))
    expect_true(length(res@.processHistory) == 3)
    expect_equal(chromPeaks(res), chromPeaks(xchrs))
    expect_equal(processHistory(res)[1:3], processHistory(xchrs)[1:3])
})

test_that("featureDefinitions,XChromatograms works", {
    skip_on_os(os = "windows", arch = "i386")

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

test_that("[,XChromatograms works", {
    skip_on_os(os = "windows", arch = "i386")

    chrs <- findChromPeaks(od_chrs, param = CentWaveParam())
    prm <- PeakDensityParam(sampleGroups = c(1, 1, 1))
    chrs <- groupChromPeaks(chrs, param = prm)

    pks <- chromPeaks(chrs)
    fts <- featureDefinitions(chrs)

    res <- chrs[2, 2]
    expect_true(is(res, "XChromatogram"))
    expect_equal(chromPeaks(res), pks[pks[, "row"] == 2 &
                                      pks[, "column"] == 2,
                                      colnames(chromPeaks(res))])
    res <- chrs[2, 2, drop = FALSE]
    expect_true(is(res, "XChromatograms"))
    expect_equal(chromPeaks(res)[, 1:7],
                 pks[pks[, "row"] == 2 & pks[, "column"] == 2, 1:7])

    res <- chrs[2, 2:3]
    expect_true(is(res, "XChromatograms"))
    expect_true(ncol(res) == 2)
    res_2 <- chrs[2, 2:3, drop = TRUE]
    expect_equal(res, res_2)
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
    expect_equal(as.character(res$sampleNames), "ko15.CDF")
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
    expect_equal(as.character(res$sampleNames), c("ko15.CDF", "ko18.CDF"))
    expect_equal(fData(res), fData(chrs)[c(2, 1), ])
    pks_tmp <- pks[c(8, 9, 12, 13, 1, 2, 3, 4, 6, 7), ]
    expect_equal(chromPeaks(res)[, 1:6], pks_tmp[, 1:6])
    res_fts <- featureDefinitions(res)
    expect_equal(rownames(res_fts), c("FT3", "FT4", "FT1", "FT2"))
    expect_equal(res_fts$peakidx, list(c(1, 3), c(2, 4), c(5, 9), c(8, 10)))

    ## Data from an XCMSnExp.
    mzm <- rbind(305.1 + c(-0.01, 0.01), 496.2 + c(-0.01, 0.01))
    xchr <- chromatogram(xod_xgrg, mz = mzm)
    pks <- chromPeaks(xchr)
    fts <- featureDefinitions(xchr)
    res <- xchr[2:1, ]
    pks_sub <- chromPeaks(res)
    fts_sub <- featureDefinitions(res)
    expect_equal(pks_sub[pks_sub[, "row"] == 1, "into"],
                 pks[pks[, "row"] == 2, "into"])
    expect_equal(pks_sub[fts_sub$peakidx[[1]], "into"],
                 pks[fts$peakidx[[4]], "into"])

    mzm <- rbind(mzm, mzm[1, ])
    xchr <- chromatogram(xod_xgrg, mz = mzm)
    pks_2 <- chromPeaks(xchr)
    expect_equal(pks_2[pks_2[, "row"] == 3, "into"],
                 pks_2[pks_2[, "row"] == 1, "into"])
    expect_error(xchr[c(2, 1, 1, 2), ], "rownames of object")

    res <- xchr[3:2, ]
    pks_2 <- chromPeaks(res)
    fts_2 <- featureDefinitions(res)
    expect_equal(pks_2, pks)
    expect_equal(fts_2, fts)

    res <- xchr[3:2, 2]
    expect_equal(chromPeaks(res)[, "into"],
                 pks[pks[, "row"] == 2 & pks[, "column"] == 2, "into"])
    expect_equal(featureDefinitions(res)[, "rtmed"],
                 fts[fts$row == 2, "rtmed"])

    expect_equal(unname(featureValues(res)[, 1]),
                 unname(chromPeaks(res)[, "into"]))
})

test_that("featureValues,XChromatograms works", {
    skip_on_os(os = "windows", arch = "i386")

    chrs <- as(od_chrs, "XChromatograms")
    expect_error(featureValues(chrs))
    chrs <- findChromPeaks(chrs, param = CentWaveParam())
    expect_error(featureValues(chrs))
    prm <- PeakDensityParam(sampleGroups = c(1, 1, 1))
    chrs <- groupChromPeaks(chrs, param = prm)

    vls <- featureValues(chrs, value = "index")
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

test_that("plotChromPeakDensity,XChromatograms works", {
    skip_on_os(os = "windows", arch = "i386")

    chrs <- as(od_chrs, "XChromatograms")
    chrs <- findChromPeaks(chrs, param = CentWaveParam())
    prm <- PeakDensityParam(sampleGroups = c(1, 1, 1))
    chrs <- groupChromPeaks(chrs, param = prm)
    expect_error(plotChromPeakDensity(chrs))
    frst <- chrs[1, ]
    library(RColorBrewer)
    plotChromPeakDensity(frst, peakBg = paste0(brewer.pal(7, "Set1"), 60),
                         peakPch = 16, peakCol = paste0(brewer.pal(7, "Set1")))
    plotChromPeakDensity(frst, peakBg = paste0(brewer.pal(7, "Set1"), 60),
                         peakPch = 16, peakCol = paste0(brewer.pal(7, "Set1")),
                         simulate = FALSE)
    frst <- dropFeatureDefinitions(frst)
    expect_error(plotChromPeakDensity(frst))
    plotChromPeakDensity(frst, param = prm)
    scnd <- chrs[2, ]
    plotChromPeakDensity(scnd)
    plotChromPeakDensity(scnd[, c(1, 3)])
})

test_that("dropFilledChromPeaks,XChromatogram and XChromatograms work", {
    skip_on_os(os = "windows", arch = "i386")

    ## With filled-in data
    mzr <- matrix(c(335, 335, 344, 344), ncol = 2, byrow = TRUE)
    rtr <- matrix(c(2700, 2900, 2600, 2750), ncol = 2, byrow = TRUE)
    ## group
    xod_tmp <- groupChromPeaks(
        xod_xgr, param = PeakDensityParam(sampleGroups = rep(1, 3),
                                          minFraction = 0.25))
    xod_tmpf <- fillChromPeaks(
        xod_tmp, param = FillChromPeaksParam(fixedRt = 30))
    expect_true(.hasFilledPeaks(xod_tmpf))
    xchr <- chromatogram(xod_tmpf, rt = rtr, mz = mzr)
    expect_false(any(.hasFilledPeaks(xchr)))
    ch <- dropFilledChromPeaks(xchr[1, 1])
    expect_equal(ch, xchr[1, 1])
    ch <- dropFilledChromPeaks(xchr[1, 2])
    expect_equal(ch, xchr[1, 2])
    res <- dropFilledChromPeaks(xchr)
    expect_true(length(res@.processHistory) < length(xchr@.processHistory))
    res@.processHistory <- list()
    xchr@.processHistory <- list()
    expect_equal(res, xchr)

    xchrf <- chromatogram(xod_tmpf, rt = rtr, mz = mzr, filled = TRUE)
    res <- hasFilledChromPeaks(xchrf)
    expect_true(is.matrix(res))
    expect_true(nrow(res) == 2)
    expect_true(ncol(res) == 3)
    expect_true(all(res[1, ] == c(TRUE, FALSE, TRUE)))
    expect_true(all(res[2, ] == FALSE))
    expect_equal(nrow(chromPeaks(xchrf)), 6)
    expect_equal(chromPeakData(xchrf)$is_filled, c(TRUE, FALSE, TRUE,
                                                   FALSE, FALSE, FALSE))
    ch <- dropFilledChromPeaks(xchr[1, 1])
    expect_equal(ch, xchr[1, 1])
    ch <- dropFilledChromPeaks(xchr[1, 2])
    expect_equal(ch, xchrf[1, 2])
    res <- dropFilledChromPeaks(xchrf)
    expect_true(all(hasFilledChromPeaks(res) == FALSE))
    expect_equal(nrow(chromPeaks(res)), 4)
    expect_equal(chromPeakData(res)$is_filled, c(FALSE, FALSE, FALSE, FALSE))
    expect_true(length(res@.processHistory) < length(xchrf@.processHistory))
    expect_equal(chromPeaks(res), chromPeaks(xchr))
    expect_equal(featureDefinitions(res), featureDefinitions(xchr))
})

test_that("refineChromPeaks,XChromatograms,MergeNeighboringPeaksParam works", {
    skip_on_os(os = "windows", arch = "i386")

    mzr <- 305.1 + c(-0.01, 0.01)
    chr <- chromatogram(filterFile(xod_x, 1), mz = mzr)
    res <- refineChromPeaks(chr, MergeNeighboringPeaksParam())
    expect_equal(chromPeaks(res), chromPeaks(chr))
    expect_true(length(processHistory(res)) > length(processHistory(chr)))

    res <- refineChromPeaks(chr, MergeNeighboringPeaksParam(expandRt = 3))
    expect_true(nrow(chromPeaks(res)) < nrow(chromPeaks(chr)))
    expect_true(sum(chromPeakData(res)$merged) == 1)

    ## With multiple files:
    chr <- chromatogram(xod_x, mz = mzr)
    res <- refineChromPeaks(chr, MergeNeighboringPeaksParam(expandRt = 5,
                                                            minProp = 0))
    expect_true(sum(chromPeakData(res)$merged) == 2)
    expect_true(validObject(res))
    res <- refineChromPeaks(chr, MergeNeighboringPeaksParam(expandRt = 5))
    expect_true(sum(chromPeakData(res)$merged) == 1)
    expect_true(validObject(res))

    ## Doing peak detection from scratch
    mzr <- 462.2 + c(-0.04, 0.04)
    chr <- chromatogram(od_x, mz = mzr)
    chr <- findChromPeaks(chr, CentWaveParam())
    res <- refineChromPeaks(chr, MergeNeighboringPeaksParam(minProp = 0.5,
                                                            expandRt = 20))
    expect_true(nrow(chromPeaks(res)) < nrow(chromPeaks(chr)))
    expect_true(sum(chromPeakData(res)$merged) == 3)
})

test_that("filterColumnsIntensityAbove,XChromatograms works", {
    skip_on_os(os = "windows", arch = "i386")

    mzr <- rbind(305.1 + c(-0.01, 0.01),
                 462.2 + c(-0.04, 0.04))
    chrs <- chromatogram(xod_x, mz = mzr)

    res <- filterColumnsIntensityAbove(chrs)
    expect_true(is(res, "XChromatograms"))
    expect_equal(res, chrs)

    res <- filterColumnsIntensityAbove(chrs, threshold = 20000, value = "maxo")
    expect_true(is(res, "XChromatograms"))
    expect_equal(res[, 1], chrs[, 1])
    expect_equal(res[, 2], chrs[, 3])

    res <- filterColumnsIntensityAbove(chrs, threshold = 20000, value = "maxo",
                                       which = "all")
    expect_equal(res[, 1], chrs[, 1])
    expect_true(ncol(res) == 1)
})

test_that("filterColumnsKeepTop,XChromatograms works", {
    skip_on_os(os = "windows", arch = "i386")

    mzr <- rbind(305.1 + c(-0.01, 0.01),
                 462.2 + c(-0.04, 0.04))
    chrs <- chromatogram(xod_x, mz = mzr)

    res <- filterColumnsKeepTop(chrs, n = 3)
    expect_true(is(res, "XChromatograms"))
    expect_equal(res, chrs)

    res <- filterColumnsKeepTop(chrs, n = 1)
    expect_equal(res, chrs[, 3])

    res <- filterColumnsKeepTop(chrs, n = 1, sortBy = "maxo")
    expect_equal(res, chrs[, 3])
})

test_that("filterChromPeaks,XChromatograms works", {
    skip_on_os(os = "windows", arch = "i386")

    chrs <- chromatogram(xod_x, mz = rbind(305.1 + c(-0.01, 0.01),
                                           462.2 + c(-0.04, 0.04)))
    res <- filterChromPeaks(chrs, n = 2L)
    expect_equal(nrow(chromPeaks(res[1, 1])), 2L)
    expect_equal(nrow(chromPeaks(res[1, 2])), 0L)
    expect_equal(nrow(chromPeaks(res[1, 3])), 2L)

    chrs <- chromatogram(xod_xg, mz = rbind(305.1 + c(-0.01, 0.01),
                                            462.2 + c(-0.04, 0.04)))
    res <- filterChromPeaks(chrs, n = 2L)
    a <- featureValues(res)
    b <- featureValues(chrs)
    expect_equal(a, b[2:3, ])
})

test_that("transformIntensity,XChromatogram works", {
    skip_on_os(os = "windows", arch = "i386")

    chrs <- chromatogram(xod_x, mz = rbind(305.1 + c(-0.01, 0.01),
                                           462.2 + c(-0.04, 0.04)))

    res <- transformIntensity(chrs)

    expect_equal(intensity(res[1, 2]), intensity(chrs[1, 2]))
    expect_equal(chromPeaks(res[1, 2])[, "into"],
                 chromPeaks(chrs[1, 2])[, "into"])
    expect_equal(chromPeaks(res[1, 2])[, "maxo"],
                 chromPeaks(chrs[1, 2])[, "maxo"])

    res <- transformIntensity(chrs, log2)
    expect_equal(intensity(res[1, 2]), log2(intensity(chrs[1, 2])))
    expect_equal(chromPeaks(res[1, 2])[, "into"],
                 log2(chromPeaks(chrs[1, 2])[, "into"]))
    expect_equal(chromPeaks(res[1, 2])[, "maxo"],
                 log2(chromPeaks(chrs[1, 2])[, "maxo"]))
})
