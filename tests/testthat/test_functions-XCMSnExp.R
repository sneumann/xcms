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

