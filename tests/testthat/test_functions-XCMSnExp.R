test_that("XCMSnExp new works", {
    ## Basic contructor.
    x <- new("XCMSnExp")
    x@.processHistory <- list("a")
    expect_error(validObject(x))
    x@.processHistory <- list(ProcessHistory())
    expect_true(validObject(x))
    xod <- as(faahko_od, "XCMSnExp")
    expect_true(validObject(xod))
    ## MsFeatureData error: environment is locked
    expect_error(xod@msFeatureData$chromPeaks <- 3)
    expect_error(xod@msFeatureData$bla <- 4)
    expect_true(validObject(xod))
})

test_that("adjustRtimePeakGroups works", {
    pkGrp <- adjustRtimePeakGroups(xod_xg,
                                   param = PeakGroupsParam(minFraction = 1))
    expect_equal(colnames(pkGrp), basename(fileNames(xod_xg)))
    ## No NAs allowed across samples:
    isNa <- apply(pkGrp, MARGIN = 1, function(z) sum(is.na(z)))
    expect_true(all(isNa == 0))
    pkGrp <- adjustRtimePeakGroups(
        xod_xg, param = PeakGroupsParam(minFraction = 0.5))
    isNa <- apply(pkGrp, MARGIN = 1, function(z) sum(is.na(z)))
    expect_true(max(isNa) == 1)

    ## Test adjustRtime adjusting also MS level > 1.
    ## Artificially changing the MS level of some spectra.
    xod_mod <- xod_xg
    ## Select the spectra for MS level 2:
    idx_ms2 <- c(300:500, 300:500 + 1277, 300:500 + 2554)
    xod_mod@featureData$msLevel[idx_ms2] <- 2
    xod_mod_adj <- adjustRtime(xod_mod,
                               param = PeakGroupsParam(span = 0.4))
    ## rtime of the MS level 2 spectra are expected to be adjusted too
    expect_equal(rtime(xod_xgr), rtime(xod_mod_adj))
    expect_true(all(rtime(xod_mod)[idx_ms2] != rtime(xod_mod_adj)[idx_ms2]))

    res <- adjustRtimePeakGroups(xod_xg,
                                 param = PeakGroupsParam(subset = c(1, 3)))
    expect_equal(colnames(res), basename(fileNames(xod_xg)[c(1, 3)]))
})

test_that("plotAdjustedRtime works", {
    plotAdjustedRtime(xod_xgr, ylim = c(-20, 40))
    plotAdjustedRtime(xod_xgrg)
    expect_warning(plotAdjustedRtime(xod_x))
    expect_warning(plotAdjustedRtime(xod_xg))
})

test_that("plotChromPeakDensity works", {
    mzr <- c(305.05, 305.15)
    .plotChromPeakDensity(xod_x, mz = mzr)
    plotChromPeakDensity(xod_x, mz = mzr)

    ## Use the full range.
    .plotChromPeakDensity(xod_x)
    plotChromPeakDensity(xod_x)

    .plotChromPeakDensity(xod_x, mz = c(0, 1))
    plotChromPeakDensity(xod_x, mz = c(0, 1))
    .plotChromPeakDensity(xod_x, mz = c(300, 310), pch = 16,
                          xlim = c(2500, 4000))
    plotChromPeakDensity(xod_x, mz = c(300, 310), pch = 16,
                         xlim = c(2500, 4000))

    expect_error(plotChromPeakDensity(xod_x, mz = c(0, 1), type = "dunno"))
    expect_error(.plotChromPeakDensity(xod_x, mz = c(0, 1), type = "dunno"))
})

test_that("plotChromPeaks works", {
    ## Plot the full range.
    plotChromPeaks(xod_x)

    ## mz range
    plotChromPeaks(xod_x, ylim = c(453, 455))
    plotChromPeaks(xod_x, ylim = c(453.2, 453.201), xlim = c(2500, 3500))
})

test_that("plotChromPeakImage works", {
    plotChromPeakImage(xod_x, binSize = 30, log = FALSE)
    ## Check that it works if no peaks were found in one sample.
    tmp <- xod_x
    pks <- chromPeaks(tmp)
    pks <- pks[pks[, "sample"] != 1, ]
    chromPeaks(tmp) <- pks
    plotChromPeakImage(tmp, binSize = 30, log = FALSE)
    plotChromPeakImage(tmp, binSize = 20, log = FALSE)
    plotChromPeakImage(tmp, binSize = 10, log = FALSE, col = topo.colors(64))
})

test_that("applyAdjustedRtime works", {
    expect_error(applyAdjustedRtime(faahko_od))
    expect_equal(applyAdjustedRtime(faahko_xod), faahko_xod)
    ## Now really replacing the stuff.
    tmp <- applyAdjustedRtime(xod_r)
    expect_true(!hasAdjustedRtime(tmp))
    expect_equal(rtime(tmp), adjustedRtime(xod_r))
    expect_true(length(processHistory(tmp)) == 1)
    tmp <- applyAdjustedRtime(xod_xgrg)
    expect_true(!hasAdjustedRtime(tmp))
    expect_equal(rtime(tmp), adjustedRtime(xod_xgr))
    expect_equal(tmp, dropAdjustedRtime(tmp))
})

test_that(".concatenate_XCMSnExp works", {
    od1 <- readMSData(faahko_3_files[1], mode = "onDisk")
    od2 <- readMSData(faahko_3_files[2], mode = "onDisk")
    od3 <- readMSData(faahko_3_files[3], mode = "onDisk")
    xod1 <- findChromPeaks(od1, param = CentWaveParam(noise = 10000,
                                                      snthresh = 40))
    xod2 <- findChromPeaks(od2, param = CentWaveParam(noise = 10000,
                                                      snthresh = 40))
    xod3 <- findChromPeaks(od3, param = CentWaveParam(noise = 10000,
                                                      snthresh = 40))
    res <- .concatenate_XCMSnExp(xod1, xod2, xod3)
    expect_equal(pData(res), pData(faahko_xod))
    expect_equal(fData(res), fData(faahko_xod))
    expect_equal(chromPeaks(res), chromPeaks(faahko_xod))
    res2 <- c(xod1, xod2, xod3)
    expect_equal(pData(res), pData(res2))
    expect_equal(fData(res), fData(res2))
    expect_equal(chromPeaks(res), chromPeaks(res2))
})

test_that("filterFeatureDefinitions works", {
    tmp <- xod_xgrg
    expect_error(filterFeatureDefinitions("a"))
    expect_error(filterFeatureDefinitions(xod_xgr, 1:3))
    expect_error(filterFeatureDefinitions(tmp, c("FT01", "other")))
    expect_error(filterFeatureDefinitions(tmp, 1:4000))
    expect_error(filterFeatureDefinitions(tmp, features = c(1, 3.4)))
    tmp <- filterFeatureDefinitions(tmp, features = 11:30)
    expect_equal(nrow(featureDefinitions(tmp)), 20)
    expect_equal(rownames(featureDefinitions(tmp)), paste0("FT", 11:30))
    ## Check if we have the correct process history
    ph <- processHistory(tmp)[[length(processHistory(tmp))]]
    expect_true(is(processParam(ph), "GenericParam"))
    expect_true(processParam(ph)@fun == "filterFeatureDefinitions")
    tmp <- dropFeatureDefinitions(tmp)
    expect_true(length(processHistory(tmp)) == 3)
    ph <- processHistory(tmp)[[length(processHistory(tmp))]]
    expect_true(!is(processParam(ph), "GenericParam"))
})

test_that("featureSummary works", {
    expect_error(featureSummary(1:3))
    expect_error(featureSummary(xod_xgrg, group = 1:5))
    expect_error(featureSummary(xod_xgr))

    res <- featureSummary(xod_xgrg)
    expect_equal(colnames(res), c("count", "perc", "multi_count",
                                  "multi_perc", "rsd"))
    expect_equal(rownames(res), rownames(featureDefinitions(xod_xgrg)))
    ## Check that RSD calculation is correct.
    rsds <- apply(featureValues(xod_xgrg, value = "into", method = "maxint"),
                  MARGIN = 1, function(z) {
                      sd(z, na.rm = TRUE) / mean(z, na.rm = TRUE)
                  })
    expect_equal(rsds, res[, "rsd"])

    res <- featureSummary(xod_xgrg, group = c(2, 1, 1))
    expect_equal(colnames(res), c("count", "perc", "multi_count", "multi_perc",
                                  "rsd", "2_count", "2_perc", "2_multi_count",
                                  "2_multi_perc", "2_rsd", "1_count", "1_perc",
                                  "1_multi_count", "1_multi_perc", "1_rsd"))
    expect_equal(rownames(res), rownames(featureDefinitions(xod_xgrg)))
    expect_equal(res[, "rsd"], rsds)
    ## Sum of individual has to match overall numbers.
    expect_equal(res[, "count"], rowSums(res[, c("1_count", "2_count")]))
    expect_equal(res[, "multi_count"],
                 rowSums(res[, c("1_multi_count", "2_multi_count")]))
    rsds <- apply(featureValues(xod_xgrg, value = "into",
                                method = "maxint")[, 2:3],
                  MARGIN = 1, function(z) {
                      sd(z, na.rm = TRUE) / mean(z, na.rm = TRUE)
                  })
    expect_equal(rsds, res[, "1_rsd"])
    ## Columns except rsd are not allowed to have NAs
    expect_true(sum(is.na(res[, -grep(colnames(res), pattern = "rsd")])) == 0)

    res <- featureSummary(xod_xgrg, perSampleCounts = TRUE)
    expect_equal(colnames(res), c("count", "perc", "multi_count", "multi_perc",
                                  "rsd",
                                  basename(fileNames(xod_xgrg))))
})

test_that("overlappingFeatures works", {
    ## Errors
    expect_error(overlappingFeatures())
    expect_error(overlappingFeatures(4))
    expect_error(overlappingFeatures(xod_x))

    expect_true(length(overlappingFeatures(xod_xg)) == 0)
    res <- overlappingFeatures(xod_xg, expandRt = 60)
    expect_equal(res, list(c(5, 6), c(31, 32)))
})

test_that("exportMetaboAnalyst works", {
    expect_error(exportMetaboAnalyst(xod_x))
    expect_error(exportMetaboAnalyst(4))
    expect_error(exportMetaboAnalyst(xod_xg))
    expect_error(exportMetaboAnalyst(xod_xg, label = "group"))
    expect_error(exportMetaboAnalyst(xod_xg, label = c("a", "a")))
    res <- exportMetaboAnalyst(xod_xg, label = c("a", "a", "a"))

    expect_true(is.matrix(res))
    expect_equal(unname(res[1, ]), sampleNames(xod_x))
    expect_equal(unname(res[2, ]), c("a", "a", "a"))
    tmp <- xod_xg
    tmp$group <- c("a", "a", "a")
    res2 <- exportMetaboAnalyst(tmp, label = "group")
    expect_equal(res, res2)

    fl <- tempfile()
    exportMetaboAnalyst(tmp, file = fl, label = "group")
    res3 <- read.table(fl, sep = ",", row.names = 1, as.is = TRUE)
    colnames(res3) <- colnames(res)
    expect_equal(as.matrix(res3), res)

    res <- exportMetaboAnalyst(xod_xg, label = c("a", "a", "a"),
                               groupnames = TRUE)
    expect_equal(rownames(res), c("Sample", "Label", groupnames(xod_xg)))
})

test_that("chromPeakSpectra works", {
    ## For now we don't have MS1/MS2 data, so we have to stick to errors etc.
    expect_error(ms2_spectra_for_peaks(xod_x, method = "other"))
    expect_error(res <- chromPeakSpectra(od_x))
    expect_warning(res <- chromPeakSpectra(xod_x, return.type = "list"))
    expect_true(length(res) == nrow(chromPeaks(xod_x)))
    expect_equal(names(res), rownames(chromPeaks(xod_x)))
    expect_warning(res <- chromPeakSpectra(xod_x, return.type = "Spectra"))
    expect_true(is(res, "Spectra"))
    expect_true(length(res) == 0)
    expect_warning(res <- chromPeakSpectra(xod_x, msLevel = 1L))
    expect_true(length(res) == 0)
})

test_that("featureSpectra works", {
    ## For now we don't have MS1/MS2 data, so we have to stick to errors etc.
    expect_error(ms2_spectra_for_features(xod_x, method = "other"))
    expect_error(res <- featureSpectra(xod_x))
    expect_warning(res <- featureSpectra(xod_xg, return.type = "list"))
    expect_true(length(res) == nrow(featureDefinitions(xod_xg)))
    expect_equal(names(res), rownames(featureDefinitions(xod_xg)))
    expect_warning(res <- featureSpectra(xod_xg, return.type = "Spectra"))
    expect_true(is(res, "Spectra"))
    expect_true(length(res) == 0)
    expect_warning(res <- featureSpectra(xod_xg, msLevel = 1L))
    expect_true(length(res) == 0)
})

test_that("featureChromatograms works", {
    expect_error(featureChromatograms(xod_x))
    chrs <- featureChromatograms(xod_xgrg)
    expect_equal(nrow(chrs), nrow(featureDefinitions(xod_xgrg)))
    expect_equal(ncol(chrs), length(fileNames(xod_xgrg)))
    chrs_ext <- featureChromatograms(xod_xgrg, expandRt = 2)
    rts <- do.call(rbind, lapply(chrs, function(z) range(rtime(z))))
    rts_ext <- do.call(rbind, lapply(chrs_ext, function(z) range(rtime(z))))
    expect_true(all(rts[, 1] > rts_ext[, 1]))
    expect_true(all(rts[, 2] < rts_ext[, 2]))

    res_2 <- featureChromatograms(xod_xgrg, features = c(1, 5))
    expect_equal(chrs[1, ], res_2[1, ])
    expect_equal(chrs[5, 1], res_2[2, 1])
    expect_equal(chrs[5, 2], res_2[2, 2])
    expect_equal(chrs[5, 3], res_2[2, 3])

    res_3 <- featureChromatograms(xod_xgrg, features = c("FT01", "FT05"))
    expect_equal(res_2, res_3)

    res <- featureChromatograms(xod_xgrg, features = character())
    expect_true(nrow(res) == 0)
    expect_error(featureChromatograms(xod_xgrg,
                                      features = c(TRUE, FALSE, FALSE)))
    expect_error(featureChromatograms(xod_xgrg, features = c(100000, 1000002)))
    expect_error(featureChromatograms(xod_xgrg, features = c("a", "FT02")))

    ## Test with filled-in peaks.
    xod_tmp <- groupChromPeaks(
        xod_xgr, param = PeakDensityParam(sampleGroups = rep(1, 3),
                                          minFraction = 0.25))
    xod_tmpf <- fillChromPeaks(
        xod_tmp, param = FillChromPeaksParam(fixedRt = 30))
    fts <- c("FT036", "FT042")
    fchrs <- featureChromatograms(xod_tmp, features = fts, filled = TRUE)
    fchrsf <- featureChromatograms(xod_tmpf, features = fts, filled = TRUE)
    expect_equal(nrow(chromPeaks(fchrs)), 4)
    expect_equal(nrow(chromPeaks(fchrsf)), 6)
    expect_equal(chromPeakData(fchrsf)$is_filled, c(TRUE, FALSE, TRUE, FALSE,
                                                    FALSE, FALSE))
    expect_equal(featureDefinitions(fchrs)$peakidx, list(1, c(2, 3, 4)))
    expect_equal(featureDefinitions(fchrsf)$peakidx, list(c(2, 1, 3), 4:6))
    fchrsf2 <- featureChromatograms(xod_tmpf, features = fts, filled = FALSE)
    expect_equal(chromPeaks(fchrsf2), chromPeaks(fchrs))
    expect_equal(featureDefinitions(fchrsf2), featureDefinitions(fchrs))
    expect_equal(featureValues(fchrsf2), featureValues(fchrs))
})

test_that("highlightChromPeaks works", {
    mzr <- c(279, 279)
    rtr <- c(2700, 2850)
    chr <- chromatogram(xod_xgrg, mz = mzr, rt = rtr)
    plot(chr)
    pks <- chromPeaks(xod_xgrg, mz = mzr, rt = rtr, type = "apex_within")
    highlightChromPeaks(xod_xgrg, mz = mzr, rt = rtr)
    highlightChromPeaks(xod_xgrg, mz = mzr, rt = rtr, type = "polygon",
                        col = c("#ff000020", "#00ff0020", "#0000ff20"))
    expect_error(highlightChromPeaks(xod_xgrg, peakIds = c("a", "b")))

    highlightChromPeaks(xod_xgrg, mz = mzr, rt = c(rtr[1] - 20, rtr[2] + 20),
                        type  = "rect", col = c("#00000040"))
    plot(chr)
    highlightChromPeaks(xod_xgrg, mz = mzr, rt = c(rtr[1] - 30, rtr[2] + 20),
                        type = "polygon",
                        col = c("#ff000020", "#00ff0020", "#0000ff20"))
})

test_that(".swath_collect_chrom_peaks works", {
    ## This should be added to the msdata package!
    fl <- "/Users/jo/Projects/git/michaelwitting/metabolomics2018/data/PestMix1_SWATH.mzML"
    if (file.exists(fl)) {
        obj <- as(readMSData(fl, mode = "onDisk"), "XCMSnExp")
        msf <- new("MsFeatureData")
        msf@.xData <- .copy_env(obj@msFeatureData)
        cwp <- CentWaveParam(snthresh = 5, noise = 100, ppm = 10,
                             peakwidth = c(3, 30))
        x <- lapply(split(obj, f = isolationWindowTargetMz(obj)),
                    findChromPeaks, msLevel = 2L, param = cwp)
        res <- .swath_collect_chrom_peaks(x, msf, fileNames(obj))

        x_mod <- x
        chromPeaks(x_mod[[2]]) <- chromPeaks(x_mod[[2]])[integer(), ]
        msf <- new("MsFeatureData")
        msf@.xData <- .copy_env(obj@msFeatureData)
        res_mod <- .swath_collect_chrom_peaks(x_mod, msf, fileNames(obj))

        a <- chromPeaks(res_mod)
        b <- chromPeaks(res)[chromPeakData(res)$isolationWindowTargetMZ != 208.95, ]
        expect_equal(unname(a), unname(b))

        obj <- findChromPeaks(obj, param = cwp)
        msf <- new("MsFeatureData")
        msf@.xData <- .copy_env(obj@msFeatureData)
        x <- lapply(split(obj, f = isolationWindowTargetMz(obj)),
                    findChromPeaks, msLevel = 2L, param = cwp)
        res_2 <- .swath_collect_chrom_peaks(x, msf, fileNames(obj))
        expect_true(nrow(chromPeaks(res_2)) > nrow(chromPeaks(res)))
        expect_equal(chromPeaks(obj),
                     chromPeaks(res_2)[1:nrow(chromPeaks(obj)), ])
        expect_equal(nrow(chromPeaks(res_2)),
                     nrow(chromPeaks(obj)) + nrow(chromPeaks(res)))

        expect_equal(colnames(chromPeakData(res_2)),
                     c("ms_level", "is_filled", "isolationWindowTargetMZ",
                       "isolationWindowLowerMz",
                       "isolationWindowUpperMz"))

        expect_true(hasChromPeaks(obj))
        msf <- new("MsFeatureData")
        msf@.xData <- .copy_env(obj@msFeatureData)
        cwp <- CentWaveParam(snthresh = 200, noise = 1000, ppm = 10,
                             peakwidth = c(3, 30))
        x <- lapply(split(obj, f = isolationWindowTargetMz(obj)),
                    findChromPeaks, msLevel = 2L, param = cwp)
        res_3 <- .swath_collect_chrom_peaks(x, msf, fileNames(obj))
        ## First two isolation windows do not have any peaks.
        target_mz <- unique(isolationWindowTargetMz(obj))
        target_mz <- target_mz[!is.na(target_mz)]

        expect_equal(
            sort(intersect(res$chromPeakData$isolationWindowTargetMZ, target_mz)),
            sort(target_mz))
        expect_equal(
            sort(intersect(res_3$chromPeakData$isolationWindowTargetMZ, target_mz)),
            sort(target_mz)[-1])
        expect_equal(chromPeaks(obj),
                     chromPeaks(res_3)[1:nrow(chromPeaks(obj)), ])

        obj <- dropChromPeaks(obj)
        expect_false(hasChromPeaks(obj))
        msf <- new("MsFeatureData")
        msf@.xData <- .copy_env(obj@msFeatureData)
        x <- lapply(split(obj, f = isolationWindowTargetMz(obj)),
                    findChromPeaks, msLevel = 2L, param = cwp)
        res_4 <- .swath_collect_chrom_peaks(x, msf, fileNames(obj))

        ## No chromPeaks found:
        cwp <- CentWaveParam(snthresh = 10000, noise = 1e6)
        x <- lapply(split(obj, f = isolationWindowTargetMz(obj)),
                    findChromPeaks, msLevel = 2L, param = cwp)
        res_5 <- .swath_collect_chrom_peaks(x, msf, fileNames(obj))
        expect_equal(res_5, msf)
    }
})

test_that("findChromPeaksIsolationWindow works", {
    ## This should be added to the msdata package!
    fl <- "/Users/jo/Projects/git/michaelwitting/metabolomics2018/data/PestMix1_SWATH.mzML"
    if (file.exists(fl)) {
        ## OnDiskMSnExp
        obj <- readMSData(fl, mode = "onDisk")
        cwp <- CentWaveParam(snthresh = 5, noise = 100, ppm = 10,
                             peakwidth = c(3, 30))
        res <- findChromPeaksIsolationWindow(obj, param = cwp)
        expect_true(is(res, "XCMSnExp"))
        expect_equal(length(processHistory(res)), 1)
        expect_true(all(c("isolationWindow", "isolationWindowTargetMZ") %in%
                        colnames(chromPeakData(res))))
        expect_true(all(chromPeakData(res)$ms_level == 2L))

        ## Add to existing peaks
        obj <- findChromPeaks(obj, param = cwp)
        res_2 <- findChromPeaksIsolationWindow(obj, param = cwp)
        expect_equal(chromPeaks(res_2)[1:nrow(chromPeaks(obj)), ],
                     chromPeaks(obj))
        expect_true(length(processHistory(res_2)) == 2)

    }
    ## no isolation window/add isolation window
    expect_error(findChromPeaksIsolationWindow(od_x), "are NA")
    tmp <- od_x
    cwp <- CentWaveParam(noise = 10000, snthresh = 40)
    fData(tmp)$my_win <- 1
    res_3 <- findChromPeaksIsolationWindow(
        tmp, param = cwp, isolationWindow = fData(tmp)$my_win, msLevel = 1L)
    expect_equal(chromPeaks(xod_x), chromPeaks(res_3))
    expect_true(all(chromPeakData(res_3)$isolationWindow == 1))

    res_4 <- findChromPeaksIsolationWindow(
        xod_x, param = cwp, isolationWindow = rep(1, length(xod_x)),
        msLevel = 1L)
    expect_equal(chromPeaks(res_4)[1:nrow(chromPeaks(xod_x)), ],
                 chromPeaks(xod_x))
})
