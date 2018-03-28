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
})

test_that("plotAdjustedRtime works", {
    plotAdjustedRtime(xod_xgr)
    plotAdjustedRtime(xod_xgrg)
    expect_warning(plotAdjustedRtime(xod_x))
    expect_warning(plotAdjustedRtime(xod_xg))
})

test_that("plotChromPeakDensity works", {
    mzr <- c(305.05, 305.15)
    plotChromPeakDensity(xod_x, mz = mzr)

    ## Use the full range.
    plotChromPeakDensity(xod_x)
    
    plotChromPeakDensity(xod_x, mz = c(0, 1))
    plotChromPeakDensity(xod_x, mz = c(300, 310), pch = 16, xlim = c(2500, 4000))
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
