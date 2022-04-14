test_that("XCMSnExp new works", {
    skip_on_os(os = "windows", arch = "i386")

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
    skip_on_os(os = "windows", arch = "i386")

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
    skip_on_os(os = "windows", arch = "i386")

    plotAdjustedRtime(xod_xgr, ylim = c(-20, 40))
    plotAdjustedRtime(xod_xgrg)
    expect_warning(plotAdjustedRtime(xod_x))
    expect_warning(plotAdjustedRtime(xod_xg))
})

test_that("plotChromPeakDensity works", {
    skip_on_os(os = "windows", arch = "i386")

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
    skip_on_os(os = "windows", arch = "i386")

    ## Plot the full range.
    plotChromPeaks(xod_x)

    ## mz range
    plotChromPeaks(xod_x, ylim = c(453, 455))
    plotChromPeaks(xod_x, ylim = c(453.2, 453.201), xlim = c(2500, 3500))
})

test_that("plotChromPeakImage works", {
    skip_on_os(os = "windows", arch = "i386")

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
    skip_on_os(os = "windows", arch = "i386")

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
    skip_on_os(os = "windows", arch = "i386")

    xod1 <- filterFile(faahko_xod, 1)
    xod2 <- filterFile(faahko_xod, 2)
    xod3 <- filterFile(faahko_xod, 3)
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
    skip_on_os(os = "windows", arch = "i386")

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
    skip_on_os(os = "windows", arch = "i386")

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
    skip_on_os(os = "windows", arch = "i386")

    ## Errors
    expect_error(overlappingFeatures())
    expect_error(overlappingFeatures(4))
    expect_error(overlappingFeatures(xod_x))

    expect_true(length(overlappingFeatures(xod_xg)) == 0)
    res <- overlappingFeatures(xod_xg, expandRt = 60)
    expect_equal(res, list(c(5, 6), c(31, 32)))
})

test_that("exportMetaboAnalyst works", {
    skip_on_os(os = "windows", arch = "i386")

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
    skip_on_os(os = "windows", arch = "i386")

    ## For now we don't have MS1/MS2 data, so we have to stick to errors etc.
    expect_error(ms2_mspectrum_for_all_peaks(xod_x, method = "other"))
    expect_error(res <- chromPeakSpectra(od_x))
    expect_warning(res <- chromPeakSpectra(xod_x, return.type = "list"))
    expect_true(length(res) == nrow(chromPeaks(xod_x)))
    expect_equal(names(res), rownames(chromPeaks(xod_x)))
    expect_warning(res <- chromPeakSpectra(xod_x, return.type = "MSpectra"))
    expect_true(is(res, "MSpectra"))
    expect_true(length(res) == 0)
    expect_warning(res <- chromPeakSpectra(xod_x, msLevel = 1L))
    expect_true(length(res) == 0)

    dta <- pest_dda

    ## ms2_mspectrum_for_peaks_from_file
    pks <- chromPeaks(dta)
    pks[, "sample"] <- 5
    res_all <- ms2_mspectrum_for_peaks_from_file(dta, pks)
    expect_equal(length(res_all), nrow(pks))
    expect_equal(names(res_all), rownames(pks))
    expect_true(any(lengths(res_all) > 1))
    tmp <- unlist(res_all)
    expect_true(all(vapply(tmp, fromFile, integer(1)) == 5L))

    res_sub <- ms2_mspectrum_for_peaks_from_file(dta, pks, method = "closest_rt")
    expect_true(all(lengths(res_sub) <= 1))
    pks[, "mz"] <- NA
    res_na <- ms2_mspectrum_for_peaks_from_file(dta, pks)
    expect_true(all(lengths(res_na) == 0))

    ## ms2_mspectrum_for_all_peaks
    res_all <- ms2_mspectrum_for_all_peaks(dta)
    expect_equal(rownames(chromPeaks(dta)), names(res_all))

    ## With subset.
    subs <- sample(1:nrow(chromPeaks(dta)), 20)
    res_subs <- ms2_mspectrum_for_all_peaks(dta, subset = subs)
    expect_true(all(lengths(res_subs[-subs]) == 0))

    ## With Spectra
    if (requireNamespace("Spectra", quietly = TRUE)) {
        res <- chromPeakSpectra(pest_dda, msLevel = 1L, return.type = "Spectra",
                                method = "closest_rt")
        expect_true(is(res, "Spectra"))
        expect_equal(rtime(res), unname(chromPeaks(pest_dda)[, "rt"]))

        res <- chromPeakSpectra(pest_dda, msLevel = 2L, return.type = "List")
        expect_true(is(res, "List"))
        expect_true(length(res) == nrow(chromPeaks(pest_dda)))
    }
})

test_that("featureSpectra works", {
    skip_on_os(os = "windows", arch = "i386")

    ## For now we don't have MS1/MS2 data, so we have to stick to errors etc.
    expect_error(ms2_mspectrum_for_features(xod_x, method = "other"))
    expect_error(res <- featureSpectra(xod_x))
    expect_warning(res <- featureSpectra(xod_xg, return.type = "list"))
    expect_true(length(res) == nrow(featureDefinitions(xod_xg)))
    expect_equal(names(res), rownames(featureDefinitions(xod_xg)))
    expect_warning(res <- featureSpectra(xod_xg, return.type = "MSpectra"))
    expect_true(is(res, "MSpectra"))
    expect_true(length(res) == 0)
    expect_warning(res <- featureSpectra(xod_xg, msLevel = 1L))
    expect_true(length(res) == 0)

    res <- featureSpectra(xod_xg, method = "closest_rt", msLevel = 1L,
                          return.type = "List")
    expect_equal(length(res), nrow(featureDefinitions(xod_xg)))
    for (i in seq_along(res)) {
        expect_true(is(res[[i]], "Spectra"))
    }

    res2 <- featureSpectra(xod_xg, msLevel = 2L, return.type = "Spectra")
    expect_true(length(res2) == 0)
})

test_that("featureChromatograms works", {
    skip_on_os(os = "windows", arch = "i386")

    expect_error(featureChromatograms(xod_x))

    fts <- rownames(featureDefinitions(xod_xgrg))
    chrs <- featureChromatograms(xod_xgrg, features = fts[c(1, 2, 1)])
    expect_true(ncol(chrs) == 3)
    expect_equal(featureDefinitions(chrs)$row, 1:3)
    expect_equal(featureValues(chrs),
                 featureValues(xod_xgrg)[fts[c(1, 2, 1)], ])
    chrs_ext <- featureChromatograms(xod_xgrg, expandRt = 2,
                                     features = fts[c(1, 2, 1)])
    rts <- do.call(rbind, lapply(chrs, function(z) range(rtime(z))))
    rts_ext <- do.call(rbind, lapply(chrs_ext, function(z) range(rtime(z))))
    expect_true(all(rts[, 1] > rts_ext[, 1]))
    expect_true(all(rts[, 2] < rts_ext[, 2]))

    expect_warning(res_n <- featureChromatograms(xod_xgrg, expandRt = 2,
                                                 features = fts[c(1, 2, 1)],
                                                 n = 1))
    expect_true(ncol(res_n) == 1)
    fvals <- featureValues(xod_xgrg, value = "maxo", method = "maxint",
                           intensity = "maxo")[fts[c(1, 2, 1)], ]
    fvals_sum <- apply(fvals, MARGIN = 2, sum, na.rm = TRUE)
    expect_true(colnames(res_n) == "ko15.CDF")
    expect_equal(chromPeaks(chrs[, 1])[, 1:11], chromPeaks(res_n)[, 1:11])
    expect_equal(featureValues(chrs[, 1]), featureValues(res_n))

    res_n <- featureChromatograms(xod_xgrg, expandRt = 2,
                                  features = fts[c(1, 2, 1)], n = 2)
    expect_true(ncol(res_n) == 2)
    expect_equal(featureValues(res_n), featureValues(chrs[, c(1, 3)]))

    res_2 <- featureChromatograms(xod_xgrg, features = c(1, 5))
    expect_equal(featureDefinitions(res_2)$row, 1:nrow(res_2))
    res_1 <- featureChromatograms(xod_xgrg, features = fts[c(1, 5)])
    expect_equal(res_1[1, ], res_2[1, ])
    expect_equal(res_1[2, 1], res_2[2, 1])
    expect_equal(res_1[2, 2], res_2[2, 2])
    expect_equal(res_1[2, 3], res_2[2, 3])

    res_3 <- featureChromatograms(xod_xgrg, features = c("FT01", "FT05"))
    expect_equal(res_2, res_3)

    res <- featureChromatograms(xod_xgrg, features = character())
    expect_true(nrow(res) == 0)
    expect_error(featureChromatograms(xod_xgrg,
                                      features = c(TRUE, FALSE, FALSE)))
    expect_error(featureChromatograms(xod_xgrg, features = c(100000, 1000002)))
    expect_error(featureChromatograms(xod_xgrg, features = c("a", "FT02")))

    ## expandMz
    res_4 <- featureChromatograms(xod_xgrg, features = c("FT01", "FT05"),
                                  expandMz = 2)
    expect_equal(mz(res_4)[1, ], mz(res_3)[1, ] + c(-2, 2))
    expect_true(sum(is.na(intensity(res_4[1, 2]))) <
                sum(is.na(intensity(res_3[1, 2]))))

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
    expect_equal(featureDefinitions(fchrsf)$peakidx, list(1:3, 4:6))
    fchrsf2 <- featureChromatograms(xod_tmpf, features = fts, filled = FALSE)
    expect_equal(chromPeaks(fchrsf2), chromPeaks(fchrs))
    expect_equal(featureDefinitions(fchrsf2), featureDefinitions(fchrs))
    expect_equal(featureValues(fchrsf2), featureValues(fchrs))
})

test_that("highlightChromPeaks works", {
    skip_on_os(os = "windows", arch = "i386")

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
    skip_on_os(os = "windows", arch = "i386")

    obj <- dropChromPeaks(pest_swth)
    msf <- new("MsFeatureData")
    ## msf@.xData <- .copy_env(obj@msFeatureData)
    cwp <- CentWaveParam(snthresh = 5, noise = 100, ppm = 10,
                         peakwidth = c(3, 30), prefilter = c(3, 1000))
    x <- lapply(split(obj, f = isolationWindowTargetMz(obj)),
                findChromPeaks, msLevel = 2L, param = cwp)
    res <- .swath_collect_chrom_peaks(x, msf, fileNames(obj))
    expect_equal(names(res), c("chromPeakData", "chromPeaks"))

    x_mod <- x
    x_mod[[2]] <- dropChromPeaks(x_mod[[2]])
    chromPeaks(x_mod[[2]]) <- chromPeaks(x_mod[[3]])[integer(), ]
    msf <- new("MsFeatureData")
    ## msf@.xData <- .copy_env(obj@msFeatureData)
    res_mod <- .swath_collect_chrom_peaks(x_mod, msf, fileNames(obj))

    a <- chromPeaks(res_mod)
    b <- chromPeaks(res)[chromPeakData(res)$isolationWindowTargetMZ != 208.95, ]
    expect_equal(unname(a), unname(b))

    ## obj <- findChromPeaks(obj, param = cwp)
    ## msf <- new("MsFeatureData")
    ## msf@.xData <- .copy_env(obj@msFeatureData)
    ## x <- lapply(split(obj, f = isolationWindowTargetMz(obj)),
    ##             findChromPeaks, msLevel = 2L, param = cwp)
    ## res_2 <- .swath_collect_chrom_peaks(x, msf, fileNames(obj))
    ## expect_true(nrow(chromPeaks(res_2)) > nrow(chromPeaks(res)))
    ## expect_equal(chromPeaks(obj),
    ##              chromPeaks(res_2)[1:nrow(chromPeaks(obj)), ])
    ## expect_equal(nrow(chromPeaks(res_2)),
    ##              nrow(chromPeaks(obj)) + nrow(chromPeaks(res)))

    ## expect_equal(colnames(chromPeakData(res_2)),
    ##              c("ms_level", "is_filled", "isolationWindow",
    ##                "isolationWindowTargetMZ",
    ##                "isolationWindowLowerMz",
    ##                "isolationWindowUpperMz"))

    ## expect_true(hasChromPeaks(obj))
    ## msf <- new("MsFeatureData")
    ## msf@.xData <- xcms:::.copy_env(obj@msFeatureData)
    ## cwp <- CentWaveParam(snthresh = 200, noise = 1000, ppm = 10,
    ##                      peakwidth = c(3, 30))
    ## x <- lapply(split(obj, f = isolationWindowTargetMz(obj)),
    ##             findChromPeaks, msLevel = 2L, param = cwp)
    ## res_3 <- xcms:::.swath_collect_chrom_peaks(x, msf, fileNames(obj))
    ## ## First two isolation windows do not have any peaks.
    ## target_mz <- unique(isolationWindowTargetMz(obj))
    ## target_mz <- target_mz[!is.na(target_mz)]

    ## expect_equal(
    ##     sort(intersect(res_3$chromPeakData$isolationWindowTargetMZ, target_mz)),
    ##     sort(target_mz)[-1])
    ## expect_equal(chromPeaks(obj),
    ##              chromPeaks(res_3)[1:nrow(chromPeaks(obj)), ])

    ## No chromPeaks found:
    cwp <- CentWaveParam(snthresh = 10000, noise = 1e6,
                         prefilter = c(4, 10000))
    x <- lapply(split(obj, f = isolationWindowTargetMz(obj)),
                findChromPeaks, msLevel = 2L, param = cwp)
    res_5 <- .swath_collect_chrom_peaks(x, msf, fileNames(obj))
    expect_equal(res_5, msf)
})

test_that("findChromPeaksIsolationWindow works", {
    skip_on_os(os = "windows", arch = "i386")

    obj <- filterRt(as(pest_swth, "OnDiskMSnExp"), rt = c(0, 300))
    cwp <- CentWaveParam(snthresh = 5, noise = 100, ppm = 10,
                         peakwidth = c(3, 30), prefilter = c(2, 1000))
    res <- findChromPeaksIsolationWindow(obj, param = cwp)
    expect_true(is(res, "XCMSnExp"))
    expect_equal(length(processHistory(res)), 1)
    expect_true(all(c("isolationWindow", "isolationWindowTargetMZ") %in%
                    colnames(chromPeakData(res))))
    expect_true(all(chromPeakData(res)$ms_level == 2L))

    ## Add to existing peaks
    obj <- findChromPeaks(obj, param = cwp)
    res_2 <- findChromPeaksIsolationWindow(obj, param = cwp)
    expect_equal(chromPeaks(res_2)[1:nrow(chromPeaks(obj)), , drop = FALSE],
                 chromPeaks(obj))
    expect_true(length(processHistory(res_2)) == 2)

    ## no isolation window/add isolation window
    expect_error(findChromPeaksIsolationWindow(od_x), "are NA")
    tmp <- od_x
    cwp <- CentWaveParam(noise = 10000, snthresh = 40, prefilter = c(3, 10000))
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

test_that("reconstructChromPeakSpectra works", {
    skip_on_os(os = "windows", arch = "i386")

    res <- reconstructChromPeakSpectra(
        pest_swth, peakId = rownames(chromPeaks(pest_swth))[5:6])
    expect_true(length(res) == 2)
    expect_true(length(intensity(res[[2]])) == 2)

    ## errors
    expect_error(reconstructChromPeakSpectra(od_x), "object with")

    ## peakId
    res_3 <- reconstructChromPeakSpectra(pest_swth, peakId = c("CP06"))
    expect_identical(intensity(res_3), intensity(res[2]))

    expect_warning(res <- reconstructChromPeakSpectra(
                       pest_swth, peakId = c("CP06", "other")))
    expect_identical(res_3, res)
    expect_error(reconstructChromPeakSpectra(pest_swth, peakId = c("a", "b")),
                 "None of the provided")
})

test_that(".plot_XIC works", {
    skip_on_os(os = "windows", arch = "i386")

    mzr <- c(453, 453.5)
    rtr <- c(2400, 2700)

    tmp <- filterMz(filterRt(xod_x, rtr), mzr)
    .plot_XIC(tmp)

    mzr <- c(301.9, 302.1)
    rtr <- c(2500, 2650)
    tmp <- filterMz(filterRt(xod_x, rtr), mzr)
    .plot_XIC(tmp, peakCol = "#ff0000", lwd = 10)
})

test_that(".group_overlapping_peaks works", {
    skip_on_os(os = "windows", arch = "i386")

    mzmin <- c(123.3, 123.35, 123.5, 341, 342.1, 343.2, 564, 564.3)
    mzmax <- c(123.4, 123.5, 124, 342, 343, 344, 564.1, 566)
    pks <- cbind(mzmin, mzmax)
    rownames(pks) <- letters[1:nrow(pks)]

    res <- .group_overlapping_peaks(pks)
    expect_true(is.list(res))
    expect_true(all(lengths(res) > 0))
    expect_equal(res[[1]], c("a", "b", "c"))

    res <- .group_overlapping_peaks(pks, expand = 0.05)
    expect_true(length(res) == 5)
    expect_equal(res[[1]], c("a", "b", "c"))
    expect_equal(res[[2]], c("d", "e"))

    res <- .group_overlapping_peaks(pks, expand = 0.1)
    expect_true(length(res) == 3)
    expect_equal(res[[1]], c("a", "b", "c"))
    expect_equal(res[[2]], c("d", "e", "f"))
    expect_equal(res[[3]], c("g", "h"))
})

test_that(".merge_neighboring_peaks works", {
    skip_on_os(os = "windows", arch = "i386")

    xod_x1 <- filterFile(xod_x, 1L)
    res <- .merge_neighboring_peaks(xod_x1, expandRt = 4)
    expect_true(is.list(res))
    expect_true(is.matrix(res$chromPeaks))
    expect_true(is(res$chromPeakData, "DataFrame"))
    expect_true(nrow(res$chromPeakData) == nrow(res$chromPeaks))
    expect_true(nrow(res$chromPeaks) < nrow(chromPeaks(xod_x1)))

    mz_groups <- .group_overlapping_peaks(chromPeaks(xod_x1), ppm = 10)
    mz_groups <- mz_groups[lengths(mz_groups) > 1]

    ## mz of 305.1: nice example of a split peak.
    tmp <- chromPeaks(xod_x1)[mz_groups[[1]], ]
    mzr <- range(tmp[, c("mzmin", "mzmax")])
    chr <- chromatogram(xod_x1, mz = mzr)
    ## plot(chr)
    pks <- res$chromPeaks
    pks <- pks[pks[, "mzmin"] >= mzr[1] & pks[, "mzmax"] <= mzr[2], ]
    expect_true(nrow(pks) == 2)
    expect_true(nrow(pks) < nrow(chromPeaks(xod_x1, mz = mzr)))
    ## rect(pks[, "rtmin"], 0, pks[, "rtmax"], pks[, "maxo"], border = "red")

    ## mz of 462.2:
    tmp <- chromPeaks(xod_x1)[mz_groups[[4]], ]
    mzr <- range(tmp[, c("mzmin", "mzmax")])
    chr <- chromatogram(xod_x1, mz = mzr)
    ## plot(chr)
    pks <- res$chromPeaks
    res_mzr <- pks[pks[, "mzmin"] >= mzr[1] & pks[, "mzmax"] <= mzr[2], , drop = FALSE]
    ## rect(res_mzr[, "rtmin"], 0, res_mzr[, "rtmax"], res_mzr[, "maxo"], border = "red")
    expect_true(nrow(res_mzr) == 1)

    ## mz of 496.2: two peaks that DON'T get merged (and that's OK).
    tmp <- chromPeaks(xod_x1)[mz_groups[[5]], ]
    mzr <- range(tmp[, c("mzmin", "mzmax")])
    mzr <- mzr + c(-0.01, 0.01)
    chr <- chromatogram(xod_x1, mz = mzr)
    pks <- res$chromPeaks
    pks <- pks[pks[, "mzmin"] >= mzr[1] & pks[, "mzmax"] <= mzr[2], ]
    ## plot(chr)
    ## rect(res_mzr[, "rtmin"], 0, res_mzr[, "rtmax"], res_mzr[, "maxo"], border = "red")
    expect_true(nrow(pks) == 2)
    expect_true(nrow(pks) == nrow(chromPeaks(xod_x1, mz = mzr)))
    expect_equal(rownames(pks), rownames(chromPeaks(chr)))
})

test_that(".XCMSnExp2SummarizedExperiment works", {
    skip_on_os(os = "windows", arch = "i386")

    expect_error(.XCMSnExp2SummarizedExperiment(xod_x), "No correspondence")

    res <- .XCMSnExp2SummarizedExperiment(xod_xgrg)
    expect_equal(SummarizedExperiment::assay(res), featureValues(xod_xgrg))

    res <- .XCMSnExp2SummarizedExperiment(xod_xgrg, value = "maxo")
    expect_equal(SummarizedExperiment::assay(res),
                 featureValues(xod_xgrg, value = "maxo"))

    res <- quantify(xod_xgrg, value = "intb")
    expect_equal(SummarizedExperiment::assay(res),
                 featureValues(xod_xgrg, value = "intb"))
})

test_that(".features_ms_region works", {
    skip_on_os(os = "windows", arch = "i386")

    res <- .features_ms_region(xod_xgrg, msLevel = 1L)
    expect_equal(nrow(res), nrow(featureDefinitions(xod_xgrg)))
    expect_equal(colnames(res), c("mzmin", "mzmax", "rtmin", "rtmax"))
    expect_true(all(res[, "mzmin"] <= res[, "mzmax"]))
    expect_true(all(res[, "rtmin"] < res[, "rtmax"]))

    expect_error(.features_ms_region(xod_xgrg, msLevel = 1L,
                                     features = c("a", "b")), "not available")
})

test_that(".which_peaks_above_threshold works", {
    skip_on_os(os = "windows", arch = "i386")

    xsub <- filterRt(filterFile(xod_x, 1L), rt = c(2500, 3500))
    pks <- chromPeaks(xsub)
    res <- .chrom_peaks_above_threshold(xsub, threshold = 100, nValues = 4)
    expect_equal(res, rep(TRUE, nrow(pks)))
    res <- .chrom_peaks_above_threshold(xsub, threshold = 50000, nValues = 1)
    expect_equal(res, unname(pks[, "maxo"] >= 50000))

    res <- .chrom_peaks_above_threshold(xsub, threshold = 50000, nValues = 1,
                                        msLevel = 2L)
    expect_equal(res, rep(TRUE, nrow(pks)))
})

test_that("manualChromPeaks works", {
    skip_on_os(os = "windows", arch = "i386")

    cp <- cbind(mzmin = c(453, 301.9, 100),
                mzmax = c(453.5, 302.1, 102),
                rtmin = c(2400, 2500, 2460),
                rtmax = c(2700, 2650, 2500))
    ## Errors
    expect_error(manualChromPeaks(xod_x, msLevel = 1:2), "can only add")
    expect_error(manualChromPeaks(1:2), "either an OnDiskMSnExp")
    expect_error(manualChromPeaks(xod_x, 1:2), "lacks one or more of the")
    expect_error(manualChromPeaks(xod_x, cp, samples = 10), "out of bounds")
    ## With an XCMSnExp
    res <- manualChromPeaks(xod_x, cp)
    expect_true(nrow(chromPeaks(res)) > nrow(chromPeaks(xod_x)))
    expect_equal(chromPeaks(res)[!is.na(chromPeaks(res)[, "intb"]), ],
                 chromPeaks(xod_x))
    ## With an OnDiskMSnExp
    res2 <- manualChromPeaks(od_x, cp)
    expect_true(is(res2, "XCMSnExp"))
    expect_true(hasChromPeaks(res2))

    res3 <- manualChromPeaks(od_x, cp, samples = 2)
    expect_true(all(chromPeaks(res3)[, "sample"] == 2))
})

test_that(".spectra_for_peaks works", {
    skip_on_os(os = "windows", arch = "i386")

    if (requireNamespace("Spectra", quietly = TRUE)) {
        res_all <- .spectra_for_peaks(pest_dda, method = "all")
        expect_true(length(res_all) == nrow(chromPeaks(pest_dda)))
        res_1 <- .spectra_for_peaks(pest_dda, method = "closest_rt")
        expect_true(all(lengths(res_1) < 2))

        res <- .spectra_for_peaks(pest_dda, msLevel = 3)
        expect_true(all(lengths(res) == 0))

        res <- .spectra_for_peaks(pest_dda, msLevel = 1L, method = "closest_rt")
        res <- do.call(c, unname(res))
        expect_equal(unname(rtime(res)), unname(chromPeaks(pest_dda)[, "rt"]))

        expect_warning(res <- .spectra_for_peaks(pest_dda, msLevel = 1L,
                                                 method = "signal"), "Changing")

        res_56 <- .spectra_for_peaks(pest_dda, method = "all", peaks = c(5, 6))
        expect_equal(res_56[[1]], res_all[[5]])
        expect_equal(res_56[[2]], res_all[[6]])
    }
})

test_that(".spectra_for_features works", {
    skip_on_os(os = "windows", arch = "i386")

    if (requireNamespace("Spectra", quietly = TRUE)) {
        res <- .spectra_for_features(xod_xgrg, method = "closest_rt",
                                     msLevel = 2L)
        expect_true(length(res) ==
                    nrow(featureDefinitions(xod_xgrg)))
        expect_true(all(lengths(res) == 0))

        res <- .spectra_for_features(xod_xgrg, method = "closest_rt",
                                     msLevel = 1L)
        expect_true(all(vapply(res, is, logical(1), "Spectra")))
        fds <- featureDefinitions(xod_xgrg)
        for (i in seq_len(nrow(fds))) {
            expect_true(all(res[[i]]$feature_id == rownames(fds)[i]))
        }

        ## with subset
        idx <- c(1, 400)
        expect_error(.spectra_for_features(
            xod_xg, msLevel = 1L, features = idx), "out of bounds")
        res_all <- .spectra_for_features(xod_xg, msLevel = 1L)
        res_sub <- .spectra_for_features(xod_xg, msLevel = 1L,
                                         features = c(5, 12, 45))
        res_sub2 <- .spectra_for_features(
            xod_xg, msLevel = 1L,
            features = rownames(featureDefinitions(xod_xg))[c(5, 12, 45)])
        expect_equal(length(res_sub), 3)
        expect_equal(rtime(res_sub[[1L]]), rtime(res_all[[5L]]))
        expect_equal(rtime(res_sub[[2L]]), rtime(res_all[[12L]]))
        expect_equal(rtime(res_sub[[3L]]), rtime(res_all[[45L]]))

        expect_equal(length(res_sub), length(res_sub2))
        expect_equal(rtime(res_sub[[1L]]), rtime(res_sub2[[1L]]))
        expect_equal(rtime(res_sub[[2L]]), rtime(res_sub2[[2L]]))
        expect_equal(rtime(res_sub[[3L]]), rtime(res_sub2[[3L]]))
    }
})

test_that("manualFeatures works", {
    skip_on_os(os = "windows", arch = "i386")

    idx <- list(1:4, c(4, "a"))
    expect_error(manualFeatures(od_x, idx), "XCMSnExp")
    ## Add features to an XCMSnExp without features.
    expect_error(manualFeatures(xod_x, idx), "out of bounds")
    idx <- list(1:4, c(5, 500, 500))
    expect_error(manualFeatures(xod_x, idx), "out of bounds")
    idx <- list(1:5, c(6, 34, 234))
    res <- manualFeatures(xod_x, idx)
    expect_true(hasFeatures(res))
    expect_true(nrow(featureDefinitions(res)) == 2)
    expect_equal(featureDefinitions(res)$peakidx, idx)
    ## Append features to an XCMSnExp.
    idx <- featureDefinitions(xod_xg)$peakidx[c(3, 5, 7)]
    res <- manualFeatures(xod_xg, idx)
    nfd <- nrow(featureDefinitions(xod_xg))
    expect_true(nrow(featureDefinitions(res)) == nfd + 3)
    expect_equal(featureDefinitions(res)[nfd + 1, "mzmed"],
                 featureDefinitions(xod_xg)[3, "mzmed"])
    expect_equal(featureDefinitions(res)[nfd + 2, "rtmin"],
                 featureDefinitions(xod_xg)[5, "rtmin"])
})
