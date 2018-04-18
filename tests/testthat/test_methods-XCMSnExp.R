test_that("XCMSnExp, XCMSnExp works", {
    rts <- rtime(faahko_od)
    rts_2 <- rtime(od_x)
    expect_equal(rts, rts_2)
    ## Test with bySample.
    rts_3 <- rtime(xod_x, bySample = TRUE)
    expect_equal(rts_3, split(rts, f = fromFile(faahko_od)))
    ## Check if rtimes are correctly ordered for bySample
    rts_4 <- rtime(filterFile(faahko_od, file = 2))
    expect_equal(rts_4, rts_3[[2]])
    rts_4 <- rtime(filterFile(faahko_od, file = 3))
    expect_equal(rts_4, rts_3[[3]])
    ## Compare with the values we get from an xcmsSet:
    rtx <- faahko_xs@rt$raw
    expect_equal(unlist(rtx, use.names = FALSE),
                 unlist(rtime(faahko_xod, bySample = TRUE), use.names = FALSE))
})

test_that("mz,XCMSnExp works", {
    mzs <- mz(faahko_od)
    ## The check below has to work, since we're calling the mz,OnDiskMSnExp.
    ## mzs_2 <- mz(od_x)
    ## expect_equal(mzs, mzs_2)
    mzs_2 <- mz(xod_x, bySample = TRUE)
    tmp <- split(mzs, fromFile(faahko_od))
    expect_equal(lapply(tmp, unlist, use.names = FALSE), mzs_2)
    ## Check if mz are correctly ordered for bySample
    mzs_3 <- mz(filterFile(faahko_od, file = 2))
    expect_equal(unlist(mzs_3, use.names = FALSE), mzs_2[[2]])
})

test_that("intensity,XCMSnExp works", {
    ints <- intensity(faahko_od)
    ## The check below has to work, since we're calling the intensity,OnDiskMSnExp.
    ## ints_2 <- intensity(od_x)
    ## expect_equal(ints, ints_2)
    ints_2 <- intensity(xod_x, bySample = TRUE)
    tmp <- split(ints, fromFile(faahko_od))
    expect_equal(lapply(tmp, unlist, use.names = FALSE), ints_2)
    ## Check if mz are correctly ordered for bySample
    ints_3 <- intensity(filterFile(faahko_od, file = 2))
    expect_equal(unlist(ints_3, use.names = FALSE), ints_2[[2]])
})

test_that("spectra,XCMSnExp works", {
    xod <- as(faahko_od, "XCMSnExp")
    res <- spectra(xod)
    res_2 <- spectra(xod, bySample = TRUE)
    expect_equal(split(res, fromFile(xod)), res_2)
    ## xod_x
    tmp <- filterRt(xod_x, rt = c(2700, 2900))
    res <- spectra(tmp)
    rts <- unlist(lapply(res, rtime))
    expect_equal(rts, rtime(tmp))
    ## Check with adjusted retention times.
    tmp2 <- filterRt(xod_xgr, rt = c(2700, 2900))
    res2 <- spectra(tmp2)
    rts2 <- unlist(lapply(res2, rtime))
    expect_equal(rts2, rtime(tmp2))
    ## Now do it on one file:
    tmp <- filterFile(xod_x, file = 2)
    res <- spectra(tmp)
    expect_equal(rtime(tmp), unlist(lapply(res, rtime)))
    tmp2 <- filterFile(xod_xgr, file = 2, keepAdjustedRtime = TRUE)
    res2 <- spectra(tmp2)
    expect_equal(rtime(tmp2), unlist(lapply(res2, rtime)))
    expect_true(sum(unlist(lapply(res2, rtime)) ==
                    unlist(lapply(res, rtime))) < length(rtime(tmp)) / 4)
    res3 <- spectra(tmp2, adjusted = FALSE)
    expect_equal(res, res3)    
    ## adjusted rt
    tmp <- filterFile(xod_xgr, file = 2, keepAdjustedRtime = TRUE)
    expect_true(hasAdjustedRtime(tmp))
    expect_true(is.character(all.equal(rtime(tmp, adjusted = FALSE),
                                       adjustedRtime(tmp))))
    res <- spectra(tmp)
    expect_equal(rtime(tmp), unlist(lapply(res, rtime)))
    res <- unlist(spectrapply(tmp, FUN = function(x) {rtime(x)}))
    expect_equal(res, adjustedRtime(tmp))
})

test_that("XCMSnExp accessors work", {
    ## Filling with data...
    xod <- as(faahko_od, "XCMSnExp")
    ## peaks 
    expect_true(!hasChromPeaks(xod))
    chromPeaks(xod) <- chromPeaks(xod_x)
    expect_true(hasChromPeaks(xod))
    expect_equal(chromPeaks(xod), chromPeaks(xod_x))
    expect_error(chromPeaks(xod) <- 4)
    tmp <- chromPeaks(xod, bySample = TRUE)
    expect_true(length(tmp) == length(fileNames(xod)))
    tmp <- do.call(rbind, tmp)
    rownames(tmp) <- NULL
    expect_equal(tmp, chromPeaks(xod))
    ## chromPeaks with rt
    all_pks <- chromPeaks(xod_x)
    pks <- chromPeaks(xod_x, rt = c(2000, 2600), type = "within")
    expect_true(nrow(pks) < nrow(all_pks))
    expect_true(all(pks[, "rtmin"] >= 2000 & pks[, "rtmax"] <= 2600))
    pks <- chromPeaks(xod_x, rt = c(2000, 2600), bySample = TRUE,
                      type = "within")
    expect_true(nrow(pks[[2]]) == 0)
    pks <- chromPeaks(xod_x, rt = c(2000, 2600), type = "any")
    expect_true(all(pks[, "rtmax"] >= 2000 & pks[, "rtmin"] <= 2600))
    pks <- chromPeaks(xod_x, rt = c(2000, 2200))
    expect_true(nrow(pks) == 0)
    pks <- chromPeaks(xod_x, rt = c(2000, 2200), bySample = TRUE)
    expect_true(all(lengths(pks) == 0))
    ## chromPeaks with mz
    pks <- chromPeaks(xod_x, mz = c(280, 281), type = "within")
    expect_true(all(pks[, "mzmin"] >= 280 & pks[, "mzmax"] <= 281))
    pks <- chromPeaks(xod_x, mz = c(280, 281), bySample = TRUE, type = "within")
    expect_true(nrow(pks[[1]]) == 0)
    expect_true(nrow(pks[[3]]) == 0)
    expect_true(nrow(pks[[2]]) == 1)
    pks <- chromPeaks(xod_x, mz = c(280, 300), bySample = FALSE, type = "within")
    expect_true(all(pks[, "mzmin"] >= 280 & pks[, "mzmax"] <= 300))
    pks <- chromPeaks(xod_x, mz = c(280, 300), bySample = FALSE, type = "any")
    expect_true(all(pks[, "mzmax"] >= 280 & pks[, "mzmin"] <= 300))
    pks <- chromPeaks(xod_x, mz = c(200, 210), bySample = FALSE)
    expect_true(nrow(pks) == 0)
    pks <- chromPeaks(xod_x, mz = c(200, 210), bySample = TRUE)
    expect_true(all(lengths(pks) == 0))
    ## chromPeaks with both
    pks <- chromPeaks(xod_x, mz = c(280, 300), rt = c(3000, 3300),
                      type = "within")
    expect_true(all(pks[, "mzmin"] >= 280 & pks[, "mzmax"] <= 300))
    expect_true(all(pks[, "rtmin"] >= 3000 & pks[, "rtmax"] <= 3300))
    pks <- chromPeaks(xod_x, mz = c(280, 300), rt = c(3000, 3300),
                      type = "any")
    expect_true(all(pks[, "mzmax"] >= 280 & pks[, "mzmin"] <= 300))
    expect_true(all(pks[, "rtmax"] >= 3000 & pks[, "rtmin"] <= 3300))
    ## Wrong assignments.
    pks <- chromPeaks(xod_x)
    pks[1, "sample"] <- 40
    expect_error(chromPeaks(xod) <- pks)
    ## featureDefinitions
    expect_true(!hasFeatures(xod))
    library(S4Vectors)
    fd <- DataFrame(faahko_xsg@groups)
    fd$peakidx <- faahko_xsg@groupidx
    featureDefinitions(xod) <- fd
    expect_true(hasChromPeaks(xod))
    expect_true(hasFeatures(xod))
    expect_equal(featureDefinitions(xod), fd)
    ## featureDefinitions with mz and/or rt range:
    feat_def <- featureDefinitions(xod_xg)
    ## Within
    mzr <- c(300, 330)
    keep_mz <- feat_def$mzmin > mzr[1] & feat_def$mzmax < mzr[2]
    expect_equal(featureDefinitions(xod_xg, mz = mzr, type = "within"),
                 feat_def[keep_mz, ])
    rtr <- c(3000, 3800)
    keep_rt <- feat_def$rtmin > rtr[1] & feat_def$rtmax < rtr[2]
    expect_equal(featureDefinitions(xod_xg, rt = rtr, type = "within"),
                 feat_def[keep_rt, ])
    expect_equal(featureDefinitions(xod_xg, rt = rtr, mz = mzr, type = "within"),
                 feat_def[keep_rt & keep_mz, ])
    ## Any
    mzr <- range(featureDefinitions(xod_xg)[2, "mzmed"])
    keep_mz <- feat_def$mzmax >= mzr[1] & feat_def$mzmin <= mzr[2]
    expect_equal(featureDefinitions(xod_xg, mz = mzr, type = "any"),
                 feat_def[keep_mz, , drop = FALSE])
    rtr <- range(3420.006)
    keep_rt <- feat_def$rtmax >= rtr[1] & feat_def$rtmin <= rtr[2]
    expect_true(nrow(featureDefinitions(xod_xg, rt = rtr, type = "within")) !=
                nrow(featureDefinitions(xod_xg, rt = rtr, type = "any")))
    expect_equal(featureDefinitions(xod_xg, rt = rtr, type = "any"),
                 feat_def[keep_rt, , drop = FALSE])
    expect_equal(featureDefinitions(xod_xg, rt = rtr, mz = mzr, type = "any"),
                 feat_def[keep_rt & keep_mz, , drop = FALSE])
    
    ## adjustedRtime
    expect_true(!hasAdjustedRtime(xod))
    expect_true(hasAdjustedRtime(xod_r))
    suppressWarnings(expect_equal(adjustedRtime(xod), NULL))
    expect_equal(rtime(xod_r, adjusted = FALSE), rtime(xod))
    expect_equal(rtime(xod_r), adjustedRtime(xod_r))
    expect_true(is.character(all.equal(rtime(xod_r), rtime(xod_r, adjusted = FALSE))))
    ## Indirect test that the ordering of the adjusted retention times matches
    ## ordering of rtime.
    ## From MSnbase version >= 2.3.9 values are ordered first by file then by
    ## spectrum.
    if (grepl("^F", names(rtime(xod_r)[1]))) {
        rts_by_sample <- adjustedRtime(xod_r, bySample = TRUE)
        rts <- adjustedRtime(xod_r)
        expect_equal(unname(rts_by_sample[[2]]),
                     unname(rts[grep(names(rts), pattern = "F2")]))
        expect_equal(unname(unlist(rts_by_sample)),
                     unname(rts))
    }
    xod2 <- xod_r
    ## Wrong assignments.
    expect_error(adjustedRtime(xod2) <- faahko_xsg@rt$corrected[1:2])
    ## bracket subset
    tmp <- xod2[1]
    expect_true(length(tmp[[1]]) == 1)
    expect_true(length(xod2[[1]]) == 1)
})

test_that("findChromPeaks,XCMSnExp works", {
    ## Call findChromPeaks on an XCMSnExp
    tmp <- findChromPeaks(xod_x, param = CentWaveParam(noise = 10000,
                                                       snthresh = 40))
    expect_equal(chromPeaks(tmp), chromPeaks(xod_x))
    ## Check that it works also on adjusted retention times:
    tmp <- findChromPeaks(xod_r, param = CentWaveParam(noise = 10000,
                                                       snthresh = 40))
    expect_true(hasAdjustedRtime(tmp))
    expect_equal(
        length(processHistory(tmp, type = .PROCSTEP.RTIME.CORRECTION)),1)
    expect_true(sum(chromPeaks(tmp)[, "rt"] != chromPeaks(xod_x)[, "rt"]) >
                ncol(chromPeaks(tmp)))
    tmp_sub <- filterFile(xod_r, file = 1, keepAdjustedRtime = TRUE)
    expect_equal(rtime(tmp_sub, adjusted = TRUE),
                 rtime(xod_r, bySample = TRUE, adjusted = TRUE)[[1]])
    spctr <- spectra(tmp_sub)
    mz_values <- lapply(spctr, mz)
    int_values <- unlist(lapply(spctr, intensity))
    res_2 <- do_findChromPeaks_centWave(mz = unlist(mz_values),
                                        int = int_values,
                                        scantime = rtime(tmp_sub,
                                                         adjusted = TRUE),
                                        valsPerSpect = lengths(mz_values),
                                        noise = 10000, snthresh = 40)
    pks <- chromPeaks(tmp)
    pks <- pks[pks[, "sample"] == 1, colnames(res_2)]
    expect_equal(res_2, pks)
    ## Second try:
    tmp <- findChromPeaks(xod_xgrg, param = CentWaveParam(noise = 10000,
                                                          snthresh = 40))
    expect_true(hasAdjustedRtime(tmp))
    expect_equal(
        length(processHistory(tmp, type = .PROCSTEP.RTIME.CORRECTION)),1)
    expect_true(sum(chromPeaks(tmp)[, "rt"] != chromPeaks(xod_x)[, "rt"]) >
                ncol(chromPeaks(tmp)))
    tmp_sub <- filterFile(xod_xgrg, file = 3, keepAdjustedRtime = TRUE)
    expect_equal(unname(rtime(tmp_sub, adjusted = TRUE)),
                 unname(rtime(xod_xgrg, bySample = TRUE, adjusted = TRUE)[[3]]))
    spctr <- spectra(tmp_sub)
    mz_values <- lapply(spctr, mz)
    int_values <- unlist(lapply(spctr, intensity))
    res_2 <- do_findChromPeaks_centWave(mz = unlist(mz_values),
                                        int = int_values,
                                        scantime = rtime(tmp_sub,
                                                         adjusted = TRUE),
                                        valsPerSpect = lengths(mz_values),
                                        noise = 10000, snthresh = 40)
    pks <- chromPeaks(tmp)
    pks <- pks[pks[, "sample"] == 3, colnames(res_2)]
    expect_equal(res_2, pks)
})

test_that("processHistory,XCMSnExp works", {
    ph <- ProcessHistory(fileIndex. = 2, info. = "For file 2")
    ph_2 <- ProcessHistory(fileIndex. = 1:2, info. = "For files 1 to 2")
    xod <- as(od_x, "XCMSnExp")
    xod@.processHistory <- list(ph, ph_2)
    expect_equal(processHistory(xod), list(ph, ph_2))
    expect_equal(processHistory(xod, fileIndex = 2), list(ph, ph_2))
    expect_equal(processHistory(xod, fileIndex = 1), list(ph_2))
    expect_error(processHistory(xod, fileIndex = 5))

    ph_3 <- XProcessHistory(fileIndex = 1, param = CentWaveParam())
    xod <- addProcessHistory(xod, ph_3)
    expect_equal(length(processHistory(xod)), 3)
    expect_equal(processHistory(xod)[[3]], ph_3)
    expect_equal(processHistory(xod, fileIndex = 1), list(ph_2, ph_3))
    expect_true(validObject(xod))
})

test_that("XCMSnExp droppers work", {
    ## How are the drop functions expected to work?
    type_feat_det <- .PROCSTEP.PEAK.DETECTION
    type_feat_algn <- .PROCSTEP.PEAK.GROUPING
    type_rt_adj <- .PROCSTEP.RTIME.CORRECTION
    ## Perform alignment.
    ## xod_xg <- groupChromPeaks(xod_x, param = PeakDensityParam())
    expect_true(hasFeatures(xod_xg))
    expect_true(hasChromPeaks(xod_x))
    expect_true(hasChromPeaks(xod_xg))
    expect_true(!hasAdjustedRtime(xod_xg))
    expect_true(length(processHistory(xod_xg, type = type_feat_algn)) == 1)
    ## Retention time adjustment.
    ## xod_xgr <- adjustRtime(xod_xg, param = PeakGroupsParam(span = 1))
    expect_true(hasChromPeaks(xod_xgr))
    expect_true(length(processHistory(xod_xgr, type = type_feat_det)) == 1)
    expect_true(!hasFeatures(xod_xgr))  ## These should have been removed
    expect_true(length(processHistory(xod_xgr, type = type_feat_algn)) == 1)
    expect_true(hasAdjustedRtime(xod_xgr))
    expect_true(length(processHistory(xod_xgr, type = type_rt_adj)) == 1)
    ## Most of the retention times are different
    expect_true(sum(chromPeaks(xod_xgr)[, "rt"] != chromPeaks(xod_x)[, "rt"]) >
                nrow(chromPeaks(xod_x)) / 2)
    expect_true(sum(rtime(xod_xgr) == rtime(xod_xg)) < length(rtime(xod_xg) / 2))
    ## Alignment after retention time adjustment.
    ## xod_xgrg <- groupChromPeaks(xod_xgr, param = PeakDensityParam())
    expect_true(hasChromPeaks(xod_xgrg))
    expect_equal(chromPeaks(xod_xgrg), chromPeaks(xod_xgr))
    expect_true(hasAdjustedRtime(xod_xgrg))
    expect_equal(rtime(xod_xgrg), rtime(xod_xgr))
    expect_equal(rtime(xod_xgrg, adjusted = FALSE), rtime(od_x))
    expect_true(length(processHistory(xod_xgr, type = type_feat_algn)) == 1)
    expect_true(hasFeatures(xod_xgrg))
    expect_true(length(processHistory(xod_xgrg, type = type_feat_algn)) == 2)
    
    ## 1) dropDetectedFeatures: delete all process history steps and all data.
    res <- dropChromPeaks(xod_x)
    expect_true(!hasChromPeaks(res))
    expect_true(length(processHistory(res, type = type_feat_det)) == 0)
    expect_true(!hasFeatures(res))
    expect_true(length(processHistory(res, type = type_feat_algn)) == 0)
    expect_true(!hasAdjustedRtime(res))
    expect_true(length(processHistory(res, type = type_rt_adj)) == 0)
    expect_equal(rtime(res), rtime(od_x))
    ##
    res <- dropChromPeaks(xod_xg)
    expect_true(!hasChromPeaks(res))
    expect_true(length(processHistory(res, type = type_feat_det)) == 0)
    expect_true(!hasFeatures(res))
    expect_true(length(processHistory(res, type = type_feat_algn)) == 0)
    expect_true(!hasAdjustedRtime(res))
    expect_true(length(processHistory(res, type = type_rt_adj)) == 0)
    expect_equal(rtime(res), rtime(od_x))
    ##
    res <- dropChromPeaks(xod_xgr)
    expect_true(!hasChromPeaks(res))
    expect_true(length(processHistory(res, type = type_feat_det)) == 0)
    expect_true(!hasFeatures(res))
    expect_true(length(processHistory(res, type = type_feat_algn)) == 0)
    expect_true(!hasAdjustedRtime(res))
    expect_true(length(processHistory(res, type = type_rt_adj)) == 0)
    expect_equal(rtime(res), rtime(od_x))
    res <- dropChromPeaks(xod_xgr, keepAdjustedRtime = TRUE)
    expect_true(hasAdjustedRtime(res))
    expect_true(length(processHistory(res, type = type_rt_adj)) == 1)
    expect_equal(rtime(res), rtime(xod_xgr))
    expect_true(!hasChromPeaks(res))
    expect_true(length(processHistory(res, type = type_feat_det)) == 0)
    expect_true(!hasFeatures(res))
    expect_true(length(processHistory(res, type = type_feat_algn)) == 0)
    ##
    res <- dropChromPeaks(xod_xgrg)
    expect_true(!hasChromPeaks(res))
    expect_true(length(processHistory(res, type = type_feat_det)) == 0)
    expect_true(!hasFeatures(res))
    expect_true(length(processHistory(res, type = type_feat_algn)) == 0)
    expect_true(!hasAdjustedRtime(res))
    expect_true(length(processHistory(res, type = type_rt_adj)) == 0)
    expect_equal(rtime(res), rtime(od_x))
    res <- dropChromPeaks(xod_xgrg, keepAdjustedRtime = TRUE)
    expect_true(hasAdjustedRtime(res))
    expect_true(length(processHistory(res, type = type_rt_adj)) == 1)
    expect_equal(rtime(res), rtime(xod_xgr))
    expect_true(!hasChromPeaks(res))
    expect_true(length(processHistory(res, type = type_feat_det)) == 0)
    expect_true(!hasFeatures(res))
    expect_true(length(processHistory(res, type = type_feat_algn)) == 0)
    
    ## 2) dropFeatureDefinitions:
    ##    a) drop the feature groups and the latest related process history
    ##    b) if retention time correction was performed AFTER the latest feature
    ##       grouping, drop also the retention time correction and all related
    ##       process histories.
    res <- dropFeatureDefinitions(xod_xg)
    expect_equal(res, xod_x)
    expect_true(hasChromPeaks(res))
    expect_true(length(processHistory(res, type = type_feat_det)) == 1)
    expect_true(!hasFeatures(res))
    expect_true(length(processHistory(res, type = type_feat_algn)) == 0)
    expect_true(!hasAdjustedRtime(res))
    expect_true(length(processHistory(res, type = type_rt_adj)) == 0)
    expect_equal(rtime(res), rtime(od_x))
    ## No feature groups - so there is nothing that this function does here.
    res <- dropFeatureDefinitions(xod_xgr)
    expect_equal(res, xod_xgr)
    expect_true(hasChromPeaks(res))
    expect_true(length(processHistory(res, type = type_feat_det)) == 1)
    expect_true(!hasFeatures(res))
    expect_true(length(processHistory(res, type = type_feat_algn)) == 1)
    expect_true(hasAdjustedRtime(res))
    expect_true(length(processHistory(res, type = type_rt_adj)) == 1)
    ## Remove the latest ones.
    res <- dropFeatureDefinitions(xod_xgrg)
    expect_equal(res, xod_xgr)
    expect_true(hasChromPeaks(res))
    expect_true(length(processHistory(res, type = type_feat_det)) == 1)
    expect_true(!hasFeatures(res))
    expect_true(length(processHistory(res, type = type_feat_algn)) == 1)
    expect_true(hasAdjustedRtime(res))
    expect_true(length(processHistory(res, type = type_rt_adj)) == 1)
    expect_equal(rtime(res, adjusted = FALSE), rtime(od_x))
    expect_equal(rtime(res, adjusted = TRUE), rtime(xod_xgr))

    ## 3) dropAdjustedRtime:
    ##    a) drop the retention time adjustment and related process histories
    ##    b) if grouping has been performed AFTER retention time correction,
    ##       drop the feature alignment and all related process histories.
    ##    c) if grouping has been performed BEFORE retention time correction,
    ##       do nothing.
    res <- dropAdjustedRtime(xod_xg)
    expect_equal(res, xod_xg)
    ## This drops also the process history for alignment.
    res <- dropAdjustedRtime(xod_xgr)
    expect_true(hasChromPeaks(res))
    expect_true(length(processHistory(res, type = type_feat_det)) == 1)
    expect_true(!hasFeatures(res))
    expect_true(length(processHistory(res, type = type_feat_algn)) == 0)
    expect_true(!hasAdjustedRtime(res))
    expect_true(length(processHistory(res, type = type_rt_adj)) == 0)
    expect_equal(chromPeaks(res), chromPeaks(xod_x))
    expect_equal(res, xod_x)
    expect_equal(rtime(res), rtime(xod_x))
    expect_equal(rtime(res), rtime(xod_xgr, adjusted = FALSE))
    ## This drops also the feature alignment performed later.
    res <- dropAdjustedRtime(xod_xgrg)
    expect_true(hasChromPeaks(res))
    expect_true(length(processHistory(res, type = type_feat_det)) == 1)
    expect_true(!hasFeatures(res))
    expect_true(length(processHistory(res, type = type_feat_algn)) == 0)
    expect_true(!hasAdjustedRtime(res))
    expect_true(length(processHistory(res, type = type_rt_adj)) == 0)
    expect_equal(chromPeaks(res), chromPeaks(xod_x))
    expect_equal(res, xod_x)
    expect_equal(rtime(res), rtime(xod_xgrg, adjusted = FALSE))
})

test_that("XCMSnExp inherited methods work", {
    ## [
    tmp_1 <- faahko_od[1:10]
    expect_warning(tmp_2 <- xod_x[1:10])
    expect_true(length(processHistory(tmp_2)) == 0)
    expect_true(!hasChromPeaks(tmp_2))
    tmp_1@processingData <- new("MSnProcess")
    tmp_2@processingData <- new("MSnProcess")
    expect_equal(tmp_1, as(tmp_2, "OnDiskMSnExp"))
    expect_error(xod_r[1, 1])
    idxs <- c(1432, 1621, 2492, 3001, 3013)
    tmp <- xod_r[idxs]
    expect_true(length(tmp) == length(idxs))
    expect_equal(mz(xod_r)[idxs], mz(tmp))
    expect_true(hasAdjustedRtime(xod_r) != hasAdjustedRtime(tmp))
    ## keeping adjusted retention times:
    tmp <- xod_r[idxs, keepAdjustedRtime = TRUE]
    expect_true(hasAdjustedRtime(tmp))
    expect_equal(rtime(xod_r)[idxs], rtime(tmp))
    ## Same with object containing also peaks and features
    expect_warning(tmp <- xod_xgrg[idxs])
    expect_true(!hasAdjustedRtime(tmp))
    expect_true(!hasChromPeaks(tmp))
    expect_true(!hasFeatures(tmp))
    expect_warning(tmp <- xod_xgrg[idxs, keepAdjusted = TRUE])
    expect_true(hasAdjustedRtime(tmp))
    expect_equal(rtime(xod_xgrg)[idxs], rtime(tmp))
    expect_true(length(processHistory(tmp)) == 1)
    
    ## [[
    spct <- xod_x[[13]]
    expect_true(is(spct, "Spectrum1"))
    expect_equal(rtime(spct), unname(rtime(xod_x)[13]))
    expect_equal(mz(spct), mz(xod_x)[[13]])
    ## Have to ensure that, if x has adjusted retention times, that these are
    ## reported in the Spectrum.
    spct <- xod_r[[13]]
    expect_equal(rtime(spct), unname(rtime(xod_r, adjusted = TRUE)[13]))
    expect_true(rtime(spct) != rtime(xod_r, adjusted = FALSE)[13])
    
    ## bin
    tmp_1 <- bin(faahko_od)
    expect_warning(tmp_2 <- bin(xod_x))
    expect_true(length(processHistory(tmp_2)) == 0)
    expect_true(!hasChromPeaks(tmp_2))
    tmp_1@processingData <- new("MSnProcess")
    tmp_2@processingData <- new("MSnProcess")
    expect_equal(tmp_1, as(tmp_2, "OnDiskMSnExp"))
    ## clean
    tmp_1 <- clean(faahko_od)
    expect_warning(tmp_2 <- clean(xod_x))
    expect_true(length(processHistory(tmp_2)) == 0)
    expect_true(!hasChromPeaks(tmp_2))
    tmp_1@processingData <- new("MSnProcess")
    tmp_2@processingData <- new("MSnProcess")
    expect_equal(tmp_1, as(tmp_2, "OnDiskMSnExp"))
    ## filterAcquisitionNum
    tmp_1 <- filterAcquisitionNum(faahko_od)
    expect_warning(tmp_2 <- filterAcquisitionNum(xod_x))
    expect_true(length(tmp_2[[1]]) > 0)
    expect_true(length(processHistory(tmp_2)) == 0)
    expect_true(!hasChromPeaks(tmp_2))
    tmp_1@processingData <- new("MSnProcess")
    tmp_2@processingData <- new("MSnProcess")
    expect_equal(tmp_1, as(tmp_2, "OnDiskMSnExp"))
    ## filterMsLevel
    tmp_1 <- filterMsLevel(faahko_od)
    tmp_2 <- filterMsLevel(xod_x)
    expect_equal(tmp_2, xod_x)
    expect_equal(length(filterMsLevel(xod_x, msLevel = 2)), 0)
    ## If we've got adjusted retention times, keep them.
    expect_warning(tmp_1 <- filterMsLevel(xod_xgr, msLevel = 1))
    expect_true(hasAdjustedRtime(tmp_1))
    expect_equal(rtime(tmp_1), rtime(xod_xgr)) # adjusted rt present
    expect_warning(tmp_1 <- filterMsLevel(xod_xgrg, msLevel = 1,
                                          keepAdjustedRtime = FALSE))
    expect_true(!hasAdjustedRtime(tmp_1))
    expect_equal(rtime(tmp_1), rtime(xod_xgr, adjusted = FALSE))
    ## normalize
    tmp_1 <- normalize(faahko_od)
    expect_warning(tmp_2 <- normalize(xod_x))
    expect_true(length(processHistory(tmp_2)) == 0)
    expect_true(!hasChromPeaks(tmp_2))
    tmp_1@processingData <- new("MSnProcess")
    tmp_2@processingData <- new("MSnProcess")
    expect_equal(tmp_1, as(tmp_2, "OnDiskMSnExp"))
    ## pickPeaks
    tmp_1 <- pickPeaks(faahko_od)
    expect_warning(tmp_2 <- pickPeaks(xod_x))
    expect_true(length(processHistory(tmp_2)) == 0)
    expect_true(!hasChromPeaks(tmp_2))
    tmp_1@processingData <- new("MSnProcess")
    tmp_2@processingData <- new("MSnProcess")
    expect_equal(tmp_1, as(tmp_2, "OnDiskMSnExp"))
    ## removePeaks
    tmp_1 <- removePeaks(faahko_od)
    expect_warning(tmp_2 <- removePeaks(xod_x))
    expect_true(length(processHistory(tmp_2)) == 0)
    expect_true(!hasChromPeaks(tmp_2))
    tmp_1@processingData <- new("MSnProcess")
    tmp_2@processingData <- new("MSnProcess")
    expect_equal(tmp_1, as(tmp_2, "OnDiskMSnExp"))
    ## smooth
    tmp_1 <- smooth(faahko_od)
    expect_warning(tmp_2 <- smooth(xod_x))
    expect_true(length(processHistory(tmp_2)) == 0)
    expect_true(!hasChromPeaks(tmp_2))
    tmp_1@processingData <- new("MSnProcess")
    tmp_2@processingData <- new("MSnProcess")
    expect_equal(tmp_1, as(tmp_2, "OnDiskMSnExp"))
})

test_that("filterFile,XCMSnExp works", {
    ## filterFile
    tmp <- filterFile(xod_x, file = 2)
    expect_error(tmp@msFeatureData$bla <- 3)
    expect_true(!hasAdjustedRtime(tmp))
    expect_true(!hasFeatures(tmp))
    expect_true(all(chromPeaks(tmp)[, "sample"] == 1))
    expect_equal(chromPeaks(tmp)[, -(ncol(chromPeaks(tmp)) - 1)],
                 chromPeaks(xod_x)[chromPeaks(xod_x)[, "sample"] == 2,
                                   -(ncol(chromPeaks(xod_x)) - 1)])
    expect_equal(fileIndex(processHistory(tmp)[[1]]), 1)
    ## check with other index.
    tmp <- filterFile(xod_x, file = c(1, 3))
    expect_true(length(tmp[[1]]) == 1)
    expect_true(!hasAdjustedRtime(tmp))
    expect_true(!hasFeatures(tmp))
    expect_true(all(chromPeaks(tmp)[, "sample"] %in% c(1, 2)))
    a <- chromPeaks(tmp)
    b <- chromPeaks(xod_x)
    expect_equal(a[, -(ncol(a) - 1)],
                 b[b[, "sample"] %in% c(1, 3), -(ncol(b) - 1)])
    expect_equal(fileIndex(processHistory(tmp)[[1]]), c(1, 2))

    ## Errors
    expect_error(filterFile(xod_x, file = 5))
    expect_error(filterFile(xod_x, file = 1:5))

    ## Little mockup to check correctness of Process history.
    od_2 <- xod_x
    od_2 <- addProcessHistory(
        od_2,
        ProcessHistory(
            type = .PROCSTEP.RTIME.CORRECTION))
    od_2 <- addProcessHistory(
        od_2,
        ProcessHistory(type = .PROCSTEP.UNKNOWN,
                       fileIndex = 2,
                       info. = "I should be here"))
    od_2 <- addProcessHistory(
        od_2,
        ProcessHistory(type = .PROCSTEP.UNKNOWN,
                       fileIndex = 1, info. = "EEEEEE"))

    tmp <- filterFile(od_2, file = 2)
    ph <- processHistory(tmp)
    expect_true(length(ph) == 2)
    expect_equal(processType(ph[[2]]), .PROCSTEP.UNKNOWN)
    b <- unlist(lapply(ph, function(z) {
        processInfo(z) == "I should be here"
    }))
    expect_true(any(b))
    b <- unlist(lapply(ph, function(z) {
        processInfo(z) == "EEEEEE"
    }))
    expect_true(!any(b))
    ## Do filterFile on xod_xg
    res <- filterFile(xod_xg, file = 2)
    expect_true(hasChromPeaks(res))
    expect_true(!hasAdjustedRtime(res))
    expect_true(!hasFeatures(res))
    tmp <- chromPeaks(xod_xg)
    expect_equal(chromPeaks(res)[, -(ncol(tmp) - 1)],
                 tmp[tmp[, "sample"] == 2, -(ncol(tmp) - 1)])
    expect_equal(rtime(res), rtime(xod_xg, bySample = TRUE)[[2]])
    ## Do filterFile on xod_xgr
    ## Should remove adjusted rts and revert the original peak rts.
    res <- filterFile(xod_xgr, file = 2)
    expect_true(hasChromPeaks(res))
    tmp <- chromPeaks(xod_xg)
    expect_equal(chromPeaks(res)[, -(ncol(tmp) - 1)],
                 tmp[tmp[, "sample"] == 2, -(ncol(tmp) - 1)])
    expect_equal(rtime(res), rtime(xod_xg, bySample = TRUE)[[2]])
    expect_true(!hasAdjustedRtime(res))
    expect_true(!hasFeatures(res))
    expect_true(length(processHistory(res)) == 1)
    expect_equal(processType(processHistory(res)[[1]]), "Peak detection")
    ## The same but keep the adjusted retention times.
    res <- filterFile(xod_xgr, file = 2, keepAdjustedRtime = TRUE)
    expect_true(hasChromPeaks(res))
    tmp <- chromPeaks(xod_xgr)
    expect_equal(chromPeaks(res)[, -(ncol(tmp) - 1)],
                 tmp[tmp[, "sample"] == 2, -(ncol(tmp) - 1)])
    ## has to be different from the ones in xod_x
    tmp <- chromPeaks(xod_x)
    expect_true(sum(chromPeaks(res)[, "rt"] == tmp[tmp[, "sample"] == 2, "rt"]) <
                nrow(tmp) / 4)
    expect_equal(rtime(res), rtime(xod_xgr, bySample = TRUE)[[2]])
    expect_equal(adjustedRtime(res), adjustedRtime(xod_xgr, bySample = TRUE)[[2]])
    expect_true(hasAdjustedRtime(res))
    expect_true(!hasFeatures(res))
    expect_true(length(processHistory(res)) == 3)
    expect_equal(processType(processHistory(res)[[1]]), "Peak detection")
    expect_equal(processType(processHistory(res)[[2]]), "Peak grouping")
    expect_equal(processType(processHistory(res)[[3]]), "Retention time correction")
    ## Do filterFile on xod_xgrg
    res <- filterFile(xod_xgrg, file = c(1, 3))
    expect_true(hasChromPeaks(res))
    tmp <- chromPeaks(xod_x)
    expect_equal(chromPeaks(res)[, -(ncol(tmp) - 1)],
                 tmp[tmp[, "sample"] %in% c(1, 3), -(ncol(tmp) - 1)])
    expect_equal(unname(rtime(res, bySample = TRUE)),
                 unname(rtime(xod_xg, bySample = TRUE)[c(1, 3)]))
    expect_true(!hasAdjustedRtime(res))
    expect_true(!hasFeatures(res))
    expect_true(length(processHistory(res)) == 1)
    expect_equal(processType(processHistory(res)[[1]]), "Peak detection")
    ## keep adjusted rtime
    res <- filterFile(xod_xgrg, file = c(1, 3), keepAdjustedRtime = TRUE)
    expect_true(hasChromPeaks(res))
    tmp <- chromPeaks(xod_xgr)
    expect_equal(chromPeaks(res)[, -(ncol(tmp) - 1)],
                 tmp[tmp[, "sample"] %in% c(1, 3), -(ncol(tmp) - 1)])
    ## has to be different from the ones in xod_x
    tmp <- chromPeaks(xod_x)
    expect_true(sum(chromPeaks(res)[, "rt"] == tmp[tmp[, "sample"] %in% c(1, 3), "rt"]) <
                nrow(tmp) / 4)
    expect_equal(rtime(res, bySample = TRUE),
                 rtime(xod_xgr, bySample = TRUE)[c(1, 3)])
    expect_equal(adjustedRtime(res, bySample = TRUE),
                 adjustedRtime(xod_xgr, bySample = TRUE)[c(1, 3)])
    expect_true(hasAdjustedRtime(res))
    expect_true(!hasFeatures(res))
    expect_true(length(processHistory(res)) == 3)
    expect_equal(processType(processHistory(res)[[1]]), "Peak detection")
    expect_equal(processType(processHistory(res)[[2]]), "Peak grouping")
    expect_equal(processType(processHistory(res)[[3]]), "Retention time correction")
})

test_that("filterMz,XCMSnExp works", {
    ## subset on xod_x
    res <- filterMz(xod_x, mz = c(300, 400))
    expect_true(length(res[[1]]) == 1)
    expect_true(length(res@spectraProcessingQueue) == 1)
    expect_true(hasChromPeaks(res))
    expect_true(all(chromPeaks(res)[, "mz"] >= 300 & chromPeaks(res)[, "mz"] <= 400))
    expect_true(nrow(chromPeaks(res)) < nrow(chromPeaks(xod_x)))
    idx <- which(chromPeaks(xod_x)[, "mzmin"] >= 300 &
                 chromPeaks(xod_x)[, "mzmax"] <= 400)
    expect_equal(chromPeaks(res), chromPeaks(xod_x)[idx, ])
    expect_true(!hasAdjustedRtime(res))
    expect_true(!hasFeatures(res))
    ## subset on xod_xg
    res <- filterMz(xod_xg, mz = c(300, 400))
    expect_true(validObject(res))
    expect_true(length(res[[1]]) == 1)
    expect_true(length(res@spectraProcessingQueue) == 1)
    expect_true(hasChromPeaks(res))
    expect_true(all(chromPeaks(res)[, "mz"] >= 300 & chromPeaks(res)[, "mz"] <= 400))
    expect_true(nrow(chromPeaks(res)) < nrow(chromPeaks(xod_x)))
    idx <- which(chromPeaks(xod_xg)[, "mzmin"] >= 300 &
                 chromPeaks(xod_xg)[, "mzmax"] <= 400)
    expect_equal(chromPeaks(res), chromPeaks(xod_xg)[idx, ])
    expect_true(!hasAdjustedRtime(res))
    expect_true(hasFeatures(res))
    expect_true(nrow(featureDefinitions(res)) < nrow(featureDefinitions(xod_xg)))
    expect_true(all(featureDefinitions(res)[, "mzmed"] >= 300 &
                    featureDefinitions(res)[, "mzmed"] <= 400))
    expect_true(all(featureDefinitions(res)[, "mzmin"] >= 300 &
                    featureDefinitions(res)[, "mzmin"] <= 400))
    expect_true(all(featureDefinitions(res)[, "mzmax"] >= 300 &
                    featureDefinitions(res)[, "mzmax"] <= 400))
    ## subset on xod_xgr
    ## o keep chromPeaks
    ## o keep adjusted rtime
    res <- filterMz(xod_xgr, mz = c(300, 400))
    expect_true(validObject(res))
    expect_true(length(res[[1]]) == 1)
    expect_true(length(res@spectraProcessingQueue) == 1)
    expect_true(hasChromPeaks(res))
    expect_true(all(chromPeaks(res)[, "mz"] >= 300 & chromPeaks(res)[, "mz"] <= 400))
    expect_true(nrow(chromPeaks(res)) < nrow(chromPeaks(xod_x)))
    idx <- which(chromPeaks(xod_xgr)[, "mzmin"] >= 300 &
                 chromPeaks(xod_xgr)[, "mzmax"] <= 400)
    expect_equal(chromPeaks(res), chromPeaks(xod_xgr)[idx, ])
    expect_true(hasAdjustedRtime(res))
    expect_equal(adjustedRtime(res), adjustedRtime(xod_xgr))
    expect_true(!hasFeatures(res))
    expect_true(length(processHistory(res)) == 3)
    ## subset xod_xgrg
    res <- filterMz(xod_xgrg, mz = c(300, 400))
    expect_true(validObject(res))
    expect_true(length(res[[1]]) == 1)
    expect_true(length(res@spectraProcessingQueue) == 1)
    expect_true(hasChromPeaks(res))
    expect_true(all(chromPeaks(res)[, "mz"] >= 300 & chromPeaks(res)[, "mz"] <= 400))
    expect_true(nrow(chromPeaks(res)) < nrow(chromPeaks(xod_xgrg)))
    idx <- which(chromPeaks(xod_xgrg)[, "mzmin"] >= 300 &
                 chromPeaks(xod_xgrg)[, "mzmax"] <= 400)
    expect_equal(chromPeaks(res), chromPeaks(xod_xgrg)[idx, ])
    expect_true(hasAdjustedRtime(res))
    expect_equal(adjustedRtime(res), adjustedRtime(xod_xgrg))
    expect_true(hasFeatures(res))
    expect_true(nrow(featureDefinitions(res)) < nrow(featureDefinitions(xod_xgrg)))
    expect_true(all(featureDefinitions(res)[, "mzmed"] >= 300 &
                    featureDefinitions(res)[, "mzmed"] <= 400))
    expect_true(all(featureDefinitions(res)[, "mzmin"] >= 300 &
                    featureDefinitions(res)[, "mzmin"] <= 400))
    expect_true(all(featureDefinitions(res)[, "mzmax"] >= 300 &
                    featureDefinitions(res)[, "mzmax"] <= 400))
    ## With groups - no groups within this range
    mzr <- c(120, 130)
    res <- filterMz(xod_xg, mz = mzr)
    expect_true(!hasFeatures(res))
    expect_true(hasChromPeaks(res))
    expect_true(all(chromPeaks(res)[, "mz"] >= 120 & chromPeaks(res)[, "mz"] <= 130))
    res <- filterMz(xod_xgrg, mz = mzr)
    expect_true(!hasFeatures(res))
    expect_true(hasChromPeaks(res))
    expect_true(all(chromPeaks(res)[, "mz"] >= 120 & chromPeaks(res)[, "mz"] <= 130))    
})

test_that("filterRt,XCMSnExp works", {
    ## xod_x
    res <- filterRt(xod_x, rt = c(2700, 2900))
    ## Check if the object is OK:
    expect_equal(pData(res), pData(xod_x))
    spct <- spectra(res)
    expect_true(length(spct) > 0)    
    ## MsFeatureData has to be locked!
    expect_error(res@msFeatureData$bla <- 3)
    ## Retention time has to be within the range.
    expect_true(all(rtime(res) >= 2700 & rtime(res) <= 2900))
    rtm <- unlist(lapply(spct, rtime))
    expect_true(all(rtm >= 2700 & rtm <= 2900))
    ## peaks have to be within the range.
    expect_true(all(chromPeaks(res)[, "rt"] >= 2700 &
                    chromPeaks(res)[, "rt"] <= 2900))
    are_within <- chromPeaks(xod_x)[, "rt"] >= 2700 &
        chromPeaks(xod_x)[, "rt"] <= 2900
    expect_equal(chromPeaks(res), chromPeaks(xod_x)[are_within,])
    ## Have a feature detection process history.
    expect_equal(length(processHistory(res)), 1)
    expect_equal(processType(processHistory(res)[[1]]),
                 .PROCSTEP.PEAK.DETECTION)
    ## filter such that we keep some spectra but no chromPeaks:
    res <- filterRt(xod_x, rt = c(4200, 4400))
    expect_true(all(rtime(res) >= 4200 & rtime(res) <= 4400))
    expect_true(!hasChromPeaks(res))
    expect_true(length(processHistory(res)) == 0)
    ## No rt
    res <- filterRt(xod_x, rt = c(10, 20))
    expect_true(length(res) == 0)

    ## xod_xg
    ## o keep also the feature groups that are within the window.
    res <- filterRt(xod_xg, rt = c(2700, 2900))
    expect_true(all(rtime(res) >= 2700 & rtime(res) <= 2900))
    expect_equal(hasChromPeaks(res), hasChromPeaks(xod_xg))
    expect_true(all(chromPeaks(res)[, "rt"] >= 2700 &
                    chromPeaks(res)[, "rt"] <= 2900))
    are_within <- chromPeaks(xod_x)[, "rt"] >= 2700 &
        chromPeaks(xod_x)[, "rt"] <= 2900
    expect_equal(chromPeaks(res), chromPeaks(xod_xg)[are_within,])
    expect_true(!hasAdjustedRtime(res))
    expect_true(hasFeatures(res))
    expect_true(all(featureDefinitions(res)$rtmed >= 2700 &
                                           featureDefinitions(res)$rtmed <= 2900))
    expect_true(nrow(featureDefinitions(res)) < nrow(featureDefinitions(xod_xg)))
    expect_true(length(processHistory(res)) == 2)
    expect_true(length(processHistory(res, type = "Peak detection")) == 1)
    expect_true(length(processHistory(res, type = "Peak grouping")) == 1)
    ## All feature idx have to match.
    expect_true(all(unlist(featureDefinitions(res)$peakidx) %in%
                    1:nrow(chromPeaks(res))))
    ## Filter such that we don't have any chromPeaks.
    res <- filterRt(xod_xg, rt = c(4200, 4400))
    expect_true(all(rtime(res) >= 4200 & rtime(res) <= 4400))
    expect_true(!hasChromPeaks(res))
    expect_true(!hasFeatures(res))
    expect_true(length(processHistory(res)) == 0)
    ## No rt
    res <- filterRt(xod_xg, rt = c(10, 20))
    expect_true(length(res) == 0)

    ## xod_xgr
    res <- filterRt(xod_xgr, rt = c(2700, 2900))
    expect_true(all(rtime(res) >= 2700 & rtime(res) <= 2900))
    expect_equal(hasChromPeaks(res), hasChromPeaks(xod_xg))
    expect_true(all(chromPeaks(res)[, "rt"] >= 2700 &
                    chromPeaks(res)[, "rt"] <= 2900))
    are_within <- chromPeaks(xod_xgr)[, "rt"] >= 2700 &
        chromPeaks(xod_xgr)[, "rt"] <= 2900
    expect_equal(chromPeaks(res), chromPeaks(xod_xgr)[are_within,])
    expect_true(hasAdjustedRtime(res))
    expect_true(all(adjustedRtime(res) >= 2700 & adjustedRtime(res) <= 2900))
    expect_true(!all(rtime(res, adjusted = FALSE) >= 2700 &
                     rtime(res, adjusted = FALSE) <= 2900))
    expect_true(!hasFeatures(res))
    expect_true(length(processHistory(res, type = "Peak detection")) == 1)
    expect_true(length(processHistory(res, type = "Peak grouping")) == 1)
    expect_true(length(processHistory(res, type = "Retention time correction")) == 1)
    ## Filter such that we don't have any chromPeaks.
    res <- filterRt(xod_xgr, rt = c(4200, 4400), adjusted = TRUE)
    expect_true(hasAdjustedRtime(res))
    expect_true(all(adjustedRtime(res) >= 4200 & adjustedRtime(res) <= 4400))
    expect_true(all(rtime(res) >= 4200 & rtime(res) <= 4400))
    expect_true(!all(rtime(res, adjusted = FALSE) >= 4200 &
                     rtime(res, adjusted = FALSE) <= 4400))
    expect_true(!hasChromPeaks(res))
    expect_true(!hasFeatures(res))
    expect_true(length(processHistory(res)) == 1)
    expect_true(length(processHistory(res, type = "Retention time correction")) == 1)
    ## No rt
    res <- filterRt(xod_xgr, rt = c(10, 20))
    expect_true(length(res) == 0)
    ## filter using raw rt
    res <- filterRt(xod_xgr, rt = c(2700, 2900), adjusted = FALSE)
    expect_true(!all(rtime(res) >= 2700 & rtime(res) <= 2900))
    expect_equal(hasChromPeaks(res), hasChromPeaks(xod_xg))
    expect_true(all(chromPeaks(res)[, "rt"] >= 2700 &
                    chromPeaks(res)[, "rt"] <= 2900))
    are_within <- chromPeaks(xod_xgr)[, "rt"] >= 2700 &
        chromPeaks(xod_xgr)[, "rt"] <= 2900
    expect_equal(chromPeaks(res), chromPeaks(xod_xgr)[are_within,])
    expect_true(hasAdjustedRtime(res))
    expect_true(!all(adjustedRtime(res) >= 2700 & adjustedRtime(res) <= 2900))
    expect_true(all(rtime(res, adjusted = FALSE) >= 2700 &
                    rtime(res, adjusted = FALSE) <= 2900))
    expect_true(!hasFeatures(res))

    ## xod_xgrg
    res <- filterRt(xod_xgrg, rt = c(2700, 2900))
    expect_true(all(rtime(res) >= 2700 & rtime(res) <= 2900))
    expect_equal(hasChromPeaks(res), hasChromPeaks(xod_xg))
    expect_true(all(chromPeaks(res)[, "rt"] >= 2700 &
                    chromPeaks(res)[, "rt"] <= 2900))
    are_within <- chromPeaks(xod_xgrg)[, "rt"] >= 2700 &
        chromPeaks(xod_xgr)[, "rt"] <= 2900
    expect_equal(chromPeaks(res), chromPeaks(xod_xgrg)[are_within,])
    expect_true(hasAdjustedRtime(res))
    expect_true(all(adjustedRtime(res) >= 2700 & adjustedRtime(res) <= 2900))
    expect_true(!all(rtime(res, adjusted = FALSE) >= 2700 &
                     rtime(res, adjusted = FALSE) <= 2900))
    expect_true(length(processHistory(res, type = "Peak detection")) == 1)
    expect_true(length(processHistory(res, type = "Peak grouping")) == 2)
    expect_true(length(processHistory(res, type = "Retention time correction")) == 1)
    expect_true(hasFeatures(res))
    expect_true(all(featureDefinitions(res)$rtmed >= 2700 &
                                           featureDefinitions(res)$rtmed <= 2900))
    ## Filter such that we don't have any chromPeaks.
    res <- filterRt(xod_xgrg, rt = c(4200, 4400), adjusted = TRUE)
    expect_true(hasAdjustedRtime(res))
    expect_true(all(adjustedRtime(res) >= 4200 & adjustedRtime(res) <= 4400))
    expect_true(all(rtime(res) >= 4200 & rtime(res) <= 4400))
    expect_true(!all(rtime(res, adjusted = FALSE) >= 4200 &
                     rtime(res, adjusted = FALSE) <= 4400))
    expect_true(!hasChromPeaks(res))
    expect_true(!hasFeatures(res))
    expect_true(length(processHistory(res)) == 1)
    expect_true(length(processHistory(res, type = "Retention time correction")) == 1)
    ## No rt
    res <- filterRt(xod_xgrg, rt = c(10, 20))
    expect_true(length(res) == 0)
    ## filter using raw rt
    res <- filterRt(xod_xgrg, rt = c(2700, 2900), adjusted = FALSE)
    expect_true(!all(rtime(res) >= 2700 & rtime(res) <= 2900))
    expect_equal(hasChromPeaks(res), hasChromPeaks(xod_xg))
    expect_true(all(chromPeaks(res)[, "rt"] >= 2700 &
                    chromPeaks(res)[, "rt"] <= 2900))
    are_within <- chromPeaks(xod_xgrg)[, "rt"] >= 2700 &
        chromPeaks(xod_xgrg)[, "rt"] <= 2900
    expect_equal(chromPeaks(res), chromPeaks(xod_xgrg)[are_within,])
    expect_true(hasAdjustedRtime(res))
    expect_true(!all(adjustedRtime(res) >= 2700 & adjustedRtime(res) <= 2900))
    expect_true(all(rtime(res, adjusted = FALSE) >= 2700 &
                    rtime(res, adjusted = FALSE) <= 2900))
    expect_true(hasFeatures(res))
    expect_true(all(featureDefinitions(res)$rtmed >= 2700 &
                                           featureDefinitions(res)$rtmed <= 2900))
})

## Test the coercion method.
test_that("as,XCMSnExp,xcmsSet works", {
    od_x <- faahko_xod
    res <- .XCMSnExp2xcmsSet(od_x)
    res <- as(od_x, "xcmsSet")
    ## Results should be the same as in xs.
    expect_equal(res@peaks, chromPeaks(od_x))
    expect_equal(res@.processHistory, processHistory(od_x))
    expect_equal(phenoData(res), pData(od_x))
    expect_equal(filepaths(res), fileNames(od_x))
    expect_equal(res@rt$raw, res@rt$corrected)
    expect_equal(res@rt$raw, rtime(od_x, bySample = TRUE))
    expect_equal(profMethod(res), "bin")
    expect_equal(profStep(res), 0.1)
    ## Can we further process this?
    sampclass(res) <- rep("K", 3)
    res <- group.density(res, minfrac = 0.5)
    ## res <- fillPeaks(res)

    ## Add groups.
    od_2 <- groupChromPeaks(
        od_x,
        param = PeakDensityParam(sampleGroups =rep(1, length(fileNames(od_x)))))
    expect_equal(unname(featureDefinitions(od_2)$peakidx), groupidx(res))

    ## rt correction
    od_3 <- adjustRtime(od_2, param = PeakGroupsParam(minFraction = 1,
                                                      span = 0.4))
    ## With groups.
    res <- as(od_2, "xcmsSet")
    ftDef <- featureDefinitions(od_2)[, -ncol(featureDefinitions(od_2))]
    ftDef <- S4Vectors::as.matrix(ftDef)
    rownames(ftDef) <- NULL
    expect_equal(res@groups, ftDef)
    expect_equal(res@groupidx, unname(featureDefinitions(od_2)$peakidx))

    ## With adjusted retention time.
    res_2 <- retcor.peakgroups(res, missing = 0, span = 0.4)
    res <- as(od_3, "xcmsSet")
    expect_true(any(unlist(res@rt$raw) != unlist(res@rt$corrected)))
    expect_equal(res@rt$corrected, res_2@rt$corrected)
    expect_equal(chromPeaks(od_3), peaks(res))
    expect_equal(peaks(res_2), peaks(res))
    
    ## Test with different binning methods:
    ## o binlin
    mfp <- MatchedFilterParam(impute = "lin", binSize = 3)
    od_2 <- od_x
    processParam(od_2@.processHistory[[1]]) <- mfp
    res <- as(od_2, "xcmsSet")
    expect_equal(profStep(res), 3)
    expect_equal(profMethod(res), "binlin")
    ## o binlinbase
    mfp <- MatchedFilterParam(impute = "linbase", binSize = 2)
    processParam(od_2@.processHistory[[1]]) <- mfp
    expect_warning(res <- as(od_2, "xcmsSet"))
    expect_equal(profStep(res), 2)
    expect_equal(profMethod(res), "binlinbase")
})

test_that("chromatogram,XCMSnExp works", {
    ## Have: od_x: OnDiskMSNnExp
    ## xod_x: XCMSnExp, with detected chromPeaks.
    ## xod_xg: with feature groups.
    ## xod_xgr: with adjusted retention times (no feature groups)
    ## xod_xgrg: adjusted rt and feature groups.

    ## XCMSnExp: TIC - can NOT compare with the reported TIC, as that is
    ## different! Eventually some background adjustment performed?
    ## BPC - CDF don't habe a BPC.
    rtr <- c(2600, 2700)
    tmp_obj <- filterFile(xod_x, file = c(1, 2))
    res <- chromatogram(tmp_obj, aggregationFun = "max", rt = rtr)
    expect_true(all(rtime(res[1, 1]) >= rtr[1]))
    expect_true(all(rtime(res[1, 1]) <= rtr[2]))
    expect_true(all(rtime(res[1, 2]) >= rtr[1]))
    expect_true(all(rtime(res[1, 2]) <= rtr[2]))
    tmp <- filterRt(filterFile(xod_x, file = 2), rt = rtr)
    expect_equal(rtime(tmp), rtime(res[1, 2]))
    ints <- spectrapply(tmp, function(z) return(max(intensity(z))))
    expect_equal(unlist(ints), intensity(res[1, 2]))
    ## Check names
    expect_equal(names(rtime(res[1, 1])), names(intensity(res[1, 1])))
    ## Assure we get the same with an OnDiskMSnExp and grouped XCMSnExp
    res_2 <- chromatogram(filterFile(od_x, file = c(1, 2)),
                          aggregationFun = "max", rt = rtr)
    expect_equal(res, res_2)
    res_3 <- chromatogram(filterFile(xod_xg, file = c(1, 2)),
                          aggregationFun = "max", rt = rtr)
    expect_equal(res, res_3)
    
    ## XCMSnExp: with mzrange and rtrange:
    mzr <- c(120, 130)
    tmp <- filterMz(xod_xg, mz = mzr)
    expect_warning(fts <- featureDefinitions(tmp))
    expect_true(is.null(fts))
    tmp <- filterRt(xod_xg, rt = rtr)
    featureDefinitions(tmp)
    res_2 <- chromatogram(xod_xg, rt = rtr, mz = mzr)
    ##

    ## XCMSnExp with adjusted rtime
    ## SEE runit.Chromatogram.R
})

test_that("signal integration is correct", {
    ## Testing the signal integration of peaks.
    ## For centWave
    tmp <- xod_xgrg
    rtr <- chromPeaks(tmp)[1, c("rtmin", "rtmax")]
    mzr <- chromPeaks(tmp)[1, c("mzmin", "mzmax")]
    chr <- chromatogram(tmp, rt = rtr, mz = mzr)
    pkInt <- sum(intensity(chr[1, 1]) *
                 ((rtr[2] - rtr[1]) / (length(chr[1, 1]) - 1)))
    expect_equal(pkInt, unname(chromPeaks(tmp)[1, "into"]))

    tmp <- filterFile(xod_xgrg, file = 2)
    idxs <- sample(1:nrow(chromPeaks(tmp)), 5)
    ## Now, for i = 20, for 6 rt I got an NA. Should I remove these measurements?
    ## idxs <- 1:nrow(chromPeaks(tmp))
    for (i in idxs) {
        rtr <- chromPeaks(tmp)[i, c("rtmin", "rtmax")]
        mzr <- chromPeaks(tmp)[i, c("mzmin", "mzmax")]
        chr <- chromatogram(tmp, rt = rtr, mz = mzr)[1, 1]
        ints <- intensity(chr)
        pkI <- sum(ints, na.rm = TRUE) * ((rtr[2] - rtr[1]) / (length(ints) - 1))
        ## cat(" ", chromPeaks(tmp)[i, "into"], " - ", pkI, "\n")
        expect_equal(unname(pkI), unname(chromPeaks(tmp)[i, "into"]))
    }
    ## pkI2 <- .getPeakInt2(tmp, chromPeaks(tmp)[idxs, , drop = FALSE])
    ## expect_equal(unname(pkI2), unname(chromPeaks(tmp)[idxs, "into"]))
    
    ## Now for matchedfilter.
    tmp <- findChromPeaks(filterFile(od_x, 2), param = MatchedFilterParam())
    rtr <- chromPeaks(tmp)[1, c("rtmin", "rtmax")]
    mzr <- chromPeaks(tmp)[1, c("mzmin", "mzmax")]
    chr <- chromatogram(tmp, rt = rtr, mz = mzr)
    pkInt <- sum(intensity(chr[1, 1]) *
                 ((rtr[2] - rtr[1]) / (length(chr[1, 1]) - 1)))
    chromPeaks(tmp)[1, "into"]
    expect_equal(pkInt, unname(chromPeaks(tmp)[1, "into"]))
    idxs <- sample(1:nrow(chromPeaks(tmp)), 5)
    ## idxs <- 1:nrow(chromPeaks(tmp))
    for (i in idxs) {
        rtr <- chromPeaks(tmp)[i, c("rtmin", "rtmax")]
        mzr <- chromPeaks(tmp)[i, c("mzmin", "mzmax")]
        chr <- chromatogram(tmp, rt = rtr, mz = mzr)[1, 1]
        ints <- intensity(chr)
        pkI <- sum(ints, na.rm = TRUE) * ((rtr[2] - rtr[1]) / (length(ints) - 1))
        ## cat(" ", chromPeaks(tmp)[i, "into"], " - ", pkI, "\n")
        expect_equal(unname(pkI), unname(chromPeaks(tmp)[i, "into"]))
    }
    ## pkI2 <- .getPeakInt2(tmp, chromPeaks(tmp)[idxs, , drop = FALSE])
    ## expect_equal(unname(pkI2), unname(chromPeaks(tmp)[idxs, "into"]))

    ## ## matchedFilter with wide mz bins.
    ## ## For matchedFilter I will have to do this on the profile matrix!
    ## tmp <- findChromPeaks(filterFile(od_x, 2),
    ##                       param = MatchedFilterParam(binSize = 2))
    ## idxs <- 1:nrow(chromPeaks(tmp))
    ## pkI2 <- .getPeakInt2(tmp, chromPeaks(tmp)[idxs, , drop = FALSE])
    ## expect_equal(unname(pkI2), unname(chromPeaks(tmp)[idxs, "into"]))
})

test_that("featureValues,XCMSnExp works", {
    fdp <- PeakDensityParam(sampleGroups = faahko_xs$class)
    od_x <- groupChromPeaks(xod_x, param = fdp)
    xs <- group(faahko_xs, method = "density")
    fvs <- featureValues(od_x, value = "into")
    expect_equal(rownames(fvs), rownames(featureDefinitions(od_x)))
    rownames(fvs) <- NULL
    colnames(fvs) <- NULL
    gvs <- groupval(xs, value = "into")
    rownames(gvs) <- NULL
    colnames(gvs) <- NULL
    expect_equal(fvs, gvs)
})

test_that("peakIndex,XCMSnExp works", {
    pkI <- .peakIndex(xod_xg)
    expect_equal(names(pkI), rownames(featureDefinitions(xod_xg)))
    expect_equal(unname(pkI), featureDefinitions(xod_xg)$peakidx)
})

test_that("MS1 MS2 data works on XCMSnExp", {
    ## That's to test stuff for issues #208 and related (also issue #214).
    ## Set every other spectra in the original files to MS2.

    ## OnDiskMSnExp: od_x
    od_mod <- faahko_od
    fDat <- fData(od_mod)
    idx_1 <- which(fDat$fileIdx == 1)
    idx_1 <- idx_1[rep(c(TRUE, FALSE), length.out = length(idx_1))]
    idx_1 <- sort(unique(c(idx_1, tail(which(fDat$fileIdx == 1)))))
    fDat[idx_1, "msLevel"] <- 2
    idx_1 <- which(fDat$fileIdx == 2)
    idx_1 <- idx_1[rep(c(TRUE, FALSE), length.out = length(idx_1))]
    idx_1 <- sort(unique(c(idx_1, tail(which(fDat$fileIdx == 2)))))
    fDat[idx_1, "msLevel"] <- 2
    idx_1 <- which(fDat$fileIdx == 3)
    idx_1 <- idx_1[rep(c(TRUE, FALSE), length.out = length(idx_1))]
    idx_1 <- sort(unique(c(idx_1, tail(which(fDat$fileIdx == 3)))))
    fDat[idx_1, "msLevel"] <- 2
    fData(od_mod) <- fDat

    res <- adjustRtime(od_mod, param = ObiwarpParam())
    res_2 <- adjustRtime(filterMsLevel(od_mod, msLevel = 1),
                         param = ObiwarpParam())
    ## Expect:
    ## - adjusted rtime of any other spectrum is identical to the
    ##   ones performed on the data sub set.
    expect_equal(res[msLevel(od_mod) == 1], res_2)
    ## - difference between raw and adjusted rtime at the end and beginning are
    ##   constant.
    res_by_file <- split(res, fromFile(od_mod))
    raw_by_file <- split(rtime(od_mod), fromFile(od_mod))
    expect_true(raw_by_file[[1]][1] != res_by_file[[1]][1])
    expect_true(raw_by_file[[2]][1] != res_by_file[[2]][1])
    expect_true(raw_by_file[[3]][1] != res_by_file[[3]][1])
    diffs <- tail(res_by_file[[1]]) - tail(raw_by_file[[1]])
    expect_equal(unname(diff(diffs)), rep(0, 5))
    diffs <- tail(res_by_file[[2]]) - tail(raw_by_file[[2]])
    expect_equal(unname(diff(diffs)), rep(0, 5))
    diffs <- tail(res_by_file[[3]]) - tail(raw_by_file[[3]])
    expect_equal(unname(diff(diffs)), rep(0, 5))    
    ## - adjusted rtime of the MS level 2 are in interpolated between rts of
    ##   MS level 2.
    ## rtime for 3 should be interpolated between 2 and 4:
    adj_fun <- approxfun(x = raw_by_file[[1]][c(2, 4)],
                         y = res_by_file[[1]][c(2, 4)])
    expect_equal(adj_fun(raw_by_file[[1]][3]), unname(res_by_file[[1]][3]))
    adj_fun <- approxfun(x = raw_by_file[[2]][c(2, 4)],
                         y = res_by_file[[2]][c(2, 4)])
    expect_equal(adj_fun(raw_by_file[[2]][3]), unname(res_by_file[[2]][3]))
    adj_fun <- approxfun(x = raw_by_file[[3]][c(2, 4)],
                         y = res_by_file[[3]][c(2, 4)])
    expect_equal(adj_fun(raw_by_file[[3]][3]), unname(res_by_file[[3]][3]))

    ## XCMSnExp: xod_x, repeat the stuff above
    xod_mod <- xod_x
    fDat <- fData(xod_mod)
    idx_1 <- which(fDat$fileIdx == 1)
    idx_1 <- idx_1[rep(c(TRUE, FALSE), length.out = length(idx_1))]
    idx_1 <- sort(unique(c(idx_1, tail(which(fDat$fileIdx == 1)))))
    fDat[idx_1, "msLevel"] <- 2
    idx_1 <- which(fDat$fileIdx == 2)
    idx_1 <- idx_1[rep(c(TRUE, FALSE), length.out = length(idx_1))]
    idx_1 <- sort(unique(c(idx_1, tail(which(fDat$fileIdx == 2)))))
    fDat[idx_1, "msLevel"] <- 2
    idx_1 <- which(fDat$fileIdx == 3)
    idx_1 <- idx_1[rep(c(TRUE, FALSE), length.out = length(idx_1))]
    idx_1 <- sort(unique(c(idx_1, tail(which(fDat$fileIdx == 3)))))
    fDat[idx_1, "msLevel"] <- 2
    fData(xod_mod) <- fDat

    res <- adjustRtime(xod_mod, param = ObiwarpParam())
    suppressWarnings(res_2 <- adjustRtime(filterMsLevel(xod_mod, msLevel = 1),
                                          param = ObiwarpParam()))
    ## Expect:
    ## - adjusted rtime of any other spectrum is identical to the
    ##   ones performed on the data sub set.
    expect_equal(rtime(res, adjusted = TRUE)[msLevel(res) == 1],
                 rtime(res_2, adjusted = TRUE))
    ## - difference between raw and adjusted rtime at the end and beginning are
    ##   constant.
    res_by_file <- rtime(res, bySample = TRUE)
    raw_by_file <- rtime(xod_mod, bySample = TRUE)
    expect_true(raw_by_file[[1]][1] != res_by_file[[1]][1])
    expect_true(raw_by_file[[2]][1] != res_by_file[[2]][1])
    expect_true(raw_by_file[[3]][1] != res_by_file[[3]][1])
    diffs <- tail(res_by_file[[1]]) - tail(raw_by_file[[1]])
    expect_equal(unname(diff(diffs)), rep(0, 5))
    diffs <- tail(res_by_file[[2]]) - tail(raw_by_file[[2]])
    expect_equal(unname(diff(diffs)), rep(0, 5))
    diffs <- tail(res_by_file[[3]]) - tail(raw_by_file[[3]])
    expect_equal(unname(diff(diffs)), rep(0, 5))    
    ## - adjusted rtime of the MS level 2 are in interpolated between rts of
    ##   MS level 2.
    ## rtime for 3 should be interpolated between 2 and 4:
    adj_fun <- approxfun(x = raw_by_file[[1]][c(2, 4)],
                         y = res_by_file[[1]][c(2, 4)])
    expect_equal(adj_fun(raw_by_file[[1]][3]), unname(res_by_file[[1]][3]))
    adj_fun <- approxfun(x = raw_by_file[[2]][c(2, 4)],
                         y = res_by_file[[2]][c(2, 4)])
    expect_equal(adj_fun(raw_by_file[[2]][3]), unname(res_by_file[[2]][3]))
    adj_fun <- approxfun(x = raw_by_file[[3]][c(2, 4)],
                         y = res_by_file[[3]][c(2, 4)])
    expect_equal(adj_fun(raw_by_file[[3]][3]), unname(res_by_file[[3]][3]))
})

test_that("extractMsData,XCMSnExp works", {
    ## All the data
    ## all <- extractMsData(od_x)
    ## expect_equal(length(all), length(fileNames(od_x)))
    ## rts <- split(rtime(od_x), f = fromFile(od_x))
    ## expect_equal(lengths(rts), unlist(lapply(all, nrow)))
    ## On an OnDiskMSnExp with only mz
    mzr <- c(300, 302)
    res <- extractMsData(filterFile(od_x, 1:2), mz = mzr)
    expect_equal(length(res), 2)
    expect_true(all(res[[1]][, "mz"] >= mzr[1] & res[[1]][, "mz"] <= mzr[2]))
    expect_true(all(res[[2]][, "mz"] >= mzr[1] & res[[2]][, "mz"] <= mzr[2]))
    ## On an OnDiskMSnExp with only rt
    rtr <- c(2500, 2800)
    res <- extractMsData(filterFile(od_x, 1:2), rt = rtr)
    expect_true(all(res[[1]][, "rt"] >= rtr[1] & res[[1]][, "rt"] <= rtr[2]))
    expect_true(all(res[[2]][, "rt"] >= rtr[1] & res[[2]][, "rt"] <= rtr[2]))
    ## LLLLL TODO Continue here, and then add example to the extractMsData
    ## help page.
    ## On an OnDiskMSnExp with mz and rt
    res <- extractMsData(filterFile(od_x, 1:2), rt = rtr, mz = mzr)
    expect_true(all(res[[1]][, "rt"] >= rtr[1] & res[[1]][, "rt"] <= rtr[2]))
    expect_true(all(res[[2]][, "rt"] >= rtr[1] & res[[2]][, "rt"] <= rtr[2]))
    expect_true(all(res[[1]][, "mz"] >= mzr[1] & res[[1]][, "mz"] <= mzr[2]))
    expect_true(all(res[[2]][, "mz"] >= mzr[1] & res[[2]][, "mz"] <= mzr[2]))
    
    ## XCMSnExp, xod_xgr
    ## with adjusted retention times
    tmp <- filterFile(xod_xgr, 1:2, keepAdjustedRtime = TRUE)
    expect_true(hasAdjustedRtime(tmp))
    res <- extractMsData(tmp, rt = rtr, mz = mzr)
    mzs <- mz(tmp)
    rts <- rtime(tmp, bySample = TRUE, adjusted = TRUE)
    expect_true(all(res[[1]][, "rt"] >= rtr[1] & res[[1]][, "rt"] <= rtr[2]))
    expect_true(all(res[[2]][, "rt"] >= rtr[1] & res[[2]][, "rt"] <= rtr[2]))
    expect_true(all(res[[1]][, "mz"] >= mzr[1] & res[[1]][, "mz"] <= mzr[2]))
    expect_true(all(res[[2]][, "mz"] >= mzr[1] & res[[2]][, "mz"] <= mzr[2]))
    tmp_rts <- rts[[1]]
    tmp_rts <- tmp_rts[tmp_rts >= rtr[1] & tmp_rts <= rtr[2]]
    res_rts <- res[[1]][, 1]
    expect_equal(unique(res_rts), unname(tmp_rts))
    tmp_rts <- rts[[2]]
    tmp_rts <- tmp_rts[tmp_rts >= rtr[1] & tmp_rts <= rtr[2]]
    res_rts <- res[[2]][, 1]
    expect_equal(unique(res_rts), unname(tmp_rts))
    
    ## without adjusted retention times
    res_2 <- extractMsData(filterFile(xod_xgr, 1:2), adjustedRtime = FALSE,
                           rt = rtr, mz = mzr)
    expect_true(all(res_2[[1]][, "rt"] >= rtr[1] & res_2[[1]][, "rt"] <= rtr[2]))
    expect_true(all(res_2[[2]][, "rt"] >= rtr[1] & res_2[[2]][, "rt"] <= rtr[2]))
    expect_true(all(res_2[[1]][, "mz"] >= mzr[1] & res_2[[1]][, "mz"] <= mzr[2]))
    expect_true(all(res_2[[2]][, "mz"] >= mzr[1] & res_2[[2]][, "mz"] <= mzr[2]))
    ## expect_true(nrow(res[[1]]) != nrow(res_2[[1]]))
    ## expect_true(nrow(res[[2]]) != nrow(res_2[[2]]))

    ## rt and mzr out of range.
    res <- extractMsData(od_x, rt = c(6000, 6300), mz = c(0, 3))
    expect_equal(length(res), 3)
    expect_true(all(unlist(lapply(res, FUN = nrow)) == 0))
    res <- extractMsData(od_x, rt = c(6000, 6300))
    expect_equal(length(res), 3)
    expect_true(all(unlist(lapply(res, FUN = nrow)) == 0))
    ## res <- extractMsData(od_x, mz = c(0, 3))
    ## expect_equal(length(res), 3)
    ## expect_true(all(unlist(lapply(res, FUN = nrow)) == 0))
})

test_that("spectrapply and spectra,XCMSnExp work", {
    ## With adjusted retention time
    tmp <- filterFile(xod_r, file = 3, keepAdjustedRtime = TRUE)
    expect_true(hasAdjustedRtime(tmp))
    expect_true(is.character(all.equal(rtime(tmp, adjusted = FALSE),
                                       rtime(tmp, adjusted = TRUE))))
    sps <- spectra(tmp)
    expect_equal(unlist(lapply(sps, rtime)), rtime(tmp, adjusted = TRUE))
    expect_equal(unlist(spectrapply(tmp, FUN = function(x) rtime(x))),
                 rtime(tmp, adjusted = TRUE))
    sps_2 <- spectra(tmp, adjusted = FALSE)
    expect_equal(unlist(lapply(sps_2, rtime)), rtime(tmp, adjusted = FALSE))
    ## without adjusted retention time
    tmp <- filterFile(xod_x, file = 3)
    expect_true(!hasAdjustedRtime(tmp))
    sps <- spectra(tmp)
    expect_equal(unlist(lapply(sps, rtime)), rtime(tmp))
    expect_equal(unlist(lapply(sps, mz)), unlist(mz(tmp)))
    expect_equal(unlist(spectrapply(tmp, FUN = function(x) rtime(x))),
                 rtime(tmp))
})

test_that("processHistory,XCMSnExp works", {
    type_peak_det <- .PROCSTEP.PEAK.DETECTION
    type_align <- .PROCSTEP.RTIME.CORRECTION
    type_corr <- .PROCSTEP.PEAK.GROUPING
    expect_true(length(processHistory(xod_x, type = type_corr)) == 0)
    ph <- processHistory(xod_x, type = type_peak_det)
    expect_equal(as.character(class(processParam(ph[[1]]))), "CentWaveParam")

    ph <- processHistory(xod_xgrg)
    expect_true(length(ph) == 4)
    ph <- processHistory(xod_xgrg, msLevel = 1L)
    expect_true(length(ph) == 2)
    expect_equal(as.character(class(processParam(ph[[1]]))), "CentWaveParam")
})

test_that("split,XCMSnExp works", {
    xod <- as(faahko_od, "XCMSnExp")
    tmp <- split(xod_xgr, f = fromFile(xod_xgr))
    ## Split by file.
    expect_equal(spectra(tmp[[1]]), spectra(filterFile(xod, file = 1)))
    expect_equal(spectra(tmp[[3]]), spectra(filterFile(xod, file = 3)))
    ## Split by acquisitionNum.
    tmp <- filterRt(xod_xgr, rt = c(2500, 2700))
    expect_true(hasChromPeaks(tmp))
    expect_true(hasAdjustedRtime(tmp))
    tmp_2 <- split(tmp, f = acquisitionNum(tmp))
    expect_true(all(acquisitionNum(tmp_2[[1]]) == acquisitionNum(tmp)[1]))
    expect_true(all(acquisitionNum(tmp_2[[14]]) == acquisitionNum(tmp)[14]))
    ## with keepAdjustedRtime
    tmp <- split(xod_xgr, f = fromFile(xod_xgr), keepAdjustedRtime = TRUE)
    tmp_1 <- filterFile(xod_xgr, file = 1, keepAdjustedRtime = TRUE)
    expect_true(hasAdjustedRtime(tmp_1))
    expect_equal(rtime(tmp[[1]]), rtime(tmp_1))
    tmp_2 <- filterFile(xod_xgr, file = 2, keepAdjustedRtime = TRUE)
    expect_true(hasAdjustedRtime(tmp_2))
    expect_equal(rtime(tmp[[2]]), rtime(tmp_2))
    tmp_3 <- filterFile(xod_xgr, file = 3, keepAdjustedRtime = TRUE)
    expect_true(hasAdjustedRtime(tmp_3))
    expect_equal(rtime(tmp[[3]]), rtime(tmp_3))
    expect_true(!all(rtime(tmp[[3]]) == rtime(tmp[[3]], adjusted = FALSE)))
})

test_that("groupnames,XCMSnExp works", {
    gn <- groupnames(xod_xgrg)
    expect_error(groupnames(xod_x))
})

test_that("calibrate,XCMSnExp works", {
    do_plot <- FALSE
    tmp <- filterFile(faahko_xod, file = 1)
    
    ## Check shift calibration.
    mzs <- chromPeaks(tmp)[c(3, 6, 7, 13, 17, 32, 45)]
    mzs_shift <- mzs + 0.0001
    prm <- CalibrantMassParam(mz = mzs_shift, method = "shift")
    res <- calibrate(tmp, prm)
    expect_true(isCalibrated(res))
    expect_equal(chromPeaks(tmp)[, -1], chromPeaks(res)[, -1])
    expect_equal(chromPeaks(tmp)[, 1] + 0.0001, chromPeaks(res)[, 1])
    diffs <- chromPeaks(res)[, "mz"] - chromPeaks(tmp)[, "mz"]
    X <- chromPeaks(res)[, "mz"]
    if (do_plot)
        plot(X, diffs)

    ## Check linear.
    mzs_lin <- mzs + 0.00005 + mzs * 0.000002
    max_dif <- max(mzs_lin - mzs)
    prm <- CalibrantMassParam(mz = mzs_lin, method = "linear", mzabs = max_dif)
    res <- calibrate(tmp, prm)
    expect_true(isCalibrated(res))
    diffs <- chromPeaks(res)[, "mz"] - chromPeaks(tmp)[, "mz"]
    X <- chromPeaks(res)[, "mz"]
    if (do_plot)
        plot(X, diffs)
    res_lm <- lm(diffs ~ X)
    expect_equal(unname(coefficients(res_lm)[1]), 0.00005, tolerance = 1e-5)
    expect_equal(unname(coefficients(res_lm)[2]), 0.000002, tolerance = 1e-5)

    ## edgeshift
    prm <- CalibrantMassParam(mz = mzs_lin, method = "edgeshift",
                              mzabs = max_dif)    
    res <- calibrate(tmp, prm)
    expect_true(isCalibrated(res))
    diffs <- chromPeaks(res)[, "mz"] - chromPeaks(tmp)[, "mz"]
    X <- chromPeaks(res)[, "mz"]
    if (do_plot)
        plot(X, diffs)
    mz_sorted <- chromPeaks(tmp)[, "mz"]
    ## Diff has to be constant before and after the linear range.
    lower_idx <- which(chromPeaks(tmp)[, "mz"] < min(mzs))
    expect_true(all(diffs[lower_idx] == diffs[lower_idx][1]))
    upper_idx <- which(chromPeaks(tmp)[, "mz"] > max(mzs))
    expect_true(all(diffs[upper_idx] == diffs[upper_idx][1]))
    lin_idx <- 1:length(diffs)
    lin_idx <- lin_idx[!(lin_idx %in% lower_idx)]
    lin_idx <- lin_idx[!(lin_idx %in% upper_idx)]
    lin_mod <- lm(diffs[lin_idx] ~ X[lin_idx])
    expect_equal(unname(coefficients(lin_mod)[1]), 0.00005, tolerance = 1e-5)
    expect_equal(unname(coefficients(lin_mod)[2]), 0.000002, tolerance = 1e-5)
    
    ## Test with a single mass, fall back to shift.
    prm <- CalibrantMassParam(mz = mzs_lin[1], method = "edgeshift",
                              mzabs = max_dif)    
    expect_warning(res <- calibrate(tmp, prm))
    diffs <- chromPeaks(res)[, "mz"] - chromPeaks(tmp)[, "mz"]
    min_diff <- min(abs(chromPeaks(tmp)[, "mz"] - mzs_lin[1]))
    expect_equal(diffs, rep(min_diff, length(diffs)))

    ## Check errors.
    expect_error(calibrate(tmp, 4))
    expect_error(calibrate(tmp, CalibrantMassParam(mz = list(mzs, mzs))))    
})

test_that("adjustRtime,peakGroups works", {
    xod <- faahko_xod
    xs <- faahko_xs
    ## Group these
    xsg <- group(xs)
    xodg <- groupChromPeaks(xod,
                            param = PeakDensityParam(sampleGroups = xs$class))
    expect_equal(peaks(xsg), chromPeaks(xodg)[, colnames(peaks(xsg))])
    expect_equal(xsg@groupidx, featureDefinitions(xodg)$peakidx)
    expect_true(length(processHistory(xodg,
                                      type = .PROCSTEP.PEAK.DETECTION)) == 1)
    expect_true(length(processHistory(xodg,
                                      type = .PROCSTEP.PEAK.GROUPING)) == 1)
    ## Now do the retention time correction
    xsr <- retcor(xsg, method = "peakgroups", missing = 0, span = 0.3)
    ## minFr <- (length(fileNames(xod)) - 1) / length(fileNames(xod))
    p <- PeakGroupsParam(minFraction = 1, span = 0.3)
    xodr <- adjustRtime(xodg, param = p)
    ## Check that we've got process histories.
    expect_true(validObject(xodr))
    expect_true(hasChromPeaks(xodr))
    expect_true(!hasFeatures(xodr))
    ## But we would like to keep the related process history step:
    expect_true(hasAdjustedRtime(xodr))
    expect_true(hasFeatures(xodg))
    ## We want to keep the process history step of the feature alignment!
    expect_true(length(processHistory(xodr,
                                      type = .PROCSTEP.PEAK.GROUPING)) == 1)
    expect_true(length(processHistory(xodr,
                                      type = .PROCSTEP.RTIME.CORRECTION)) == 1)
    ## Different from original:
    expect_true(sum(chromPeaks(xod)[, "rt"] != chromPeaks(xodr)[, "rt"]) > 200)
    expect_true(sum(chromPeaks(xod)[, "rtmin"] != chromPeaks(xodr)[, "rtmin"]) > 200)
    expect_true(sum(chromPeaks(xod)[, "rtmax"] != chromPeaks(xodr)[, "rtmax"]) > 200)
    ## between xcmsSet and XCMSnExp
    expect_equal(chromPeaks(xodr)[, colnames(peaks(xsr))], peaks(xsr))
    ## To compare the adjusted retention time we have to extract it by sample!
    ## Otherwise the ordering will not be the same, as rtime is ordered by
    ## retention time, but @rt$raw by sample.
    expect_equal(unlist(adjustedRtime(xodr, bySample = TRUE), use.names = FALSE),
                 unlist(xsr@rt$corrected, use.names = FALSE))
    ## Just to ensure - are the raw rt the same?
    expect_equal(unlist(rtime(xod, bySample = TRUE), use.names = FALSE),
                 unlist(xs@rt$raw, use.names = FALSE))
    ## Check that we get the same by supplying the peakGroupsMatrix.
    pgm <- adjustRtimePeakGroups(xodg, param = p)
    p_2 <- p
    minFraction(p_2) <- 0.5
    extraPeaks(p_2) <- 20
    peakGroupsMatrix(p_2) <- pgm
    xodr_2 <- adjustRtime(xodg, param = p_2)
    expect_equal(adjustedRtime(xodr), adjustedRtime(xodr_2))
    expect_equal(chromPeaks(xodr), chromPeaks(xodr_2))
    p_got <- processParam(
        processHistory(xodr, type = .PROCSTEP.RTIME.CORRECTION)[[1]])
    peakGroupsMatrix(p_got) <- matrix(ncol = 0, nrow = 0)
    expect_equal(p_got, p)
    expect_equal(processParam(
        processHistory(xodr_2, type = .PROCSTEP.RTIME.CORRECTION)[[1]]),
        p_2)
    ## Doing an additional grouping
    xodrg <- groupChromPeaks(xodr, param = PeakDensityParam(sampleGroups =
                                                                xs$class))
    expect_true(length(processHistory(xodrg,
                                      type = .PROCSTEP.PEAK.GROUPING)) == 2)
    expect_true(hasAdjustedRtime(xodrg))
    expect_true(hasFeatures(xodrg))
    xsrg <- group(xsr)
    expect_equal(xsrg@groupidx, featureDefinitions(xodrg)$peakidx)
    
    ## Mod settings:
    xsr <- retcor(xsg, method = "peakgroups", missing = 0, span = 1)
    xodr <- adjustRtime(xodg, param = PeakGroupsParam(minFraction = 1,
                                                      span = 1))
    expect_equal(chromPeaks(xodr)[, colnames(peaks(xsr))], peaks(xsr))
    expect_equal(unlist(adjustedRtime(xodr, bySample = TRUE), use.names = FALSE),
                 unlist(xsr@rt$corrected, use.names = FALSE))

    xsr <- retcor(xsg, method = "peakgroups", missing = 0, span = 1,
                  smooth = "linear")
    xodr <- adjustRtime(xodg, param = PeakGroupsParam(minFraction = 1,
                                                      span = 1,
                                                      smooth = "linear"))
    expect_equal(chromPeaks(xodr)[, colnames(peaks(xsr))], peaks(xsr))
    expect_equal(unlist(adjustedRtime(xodr, bySample = TRUE), use.names = FALSE),
                 unlist(xsr@rt$corrected, use.names = FALSE))

    xsr <- retcor(xsg, method = "peakgroups", missing = 0, span = 1,
                  family = "symmetric")
    xodr <- adjustRtime(xodg, param = PeakGroupsParam(minFraction = 1,
                                                      span = 1,
                                                      family = "symmetric"))
    expect_equal(chromPeaks(xodr)[, colnames(peaks(xsr))], peaks(xsr))
    expect_equal(unlist(adjustedRtime(xodr, bySample = TRUE), use.names = FALSE),
                 unlist(xsr@rt$corrected, use.names = FALSE))
    ## Dropping results.
    tmp <- dropAdjustedRtime(xodr)
    expect_equal(tmp, xod)    
})

test_that("findChromPeaks,MSWParam works", {
    od <- microtofq_od
    ## Restrict to first spectrum
    od1 <- od[1]
    sp1 <- od[[1]]
    res_1 <- do_findPeaks_MSW(mz = mz(sp1), int = intensity(sp1))
    mp <- MSWParam()
    expect_error(findChromPeaks(od1, param = mp, msLevel = 2))
    res_2 <- findChromPeaks(od1, param = mp)
    expect_equal(res_1, chromPeaks(res_2)[, colnames(res_1), drop = FALSE])
    ## Changing settings.
    snthresh(mp) <- 1
    nearbyPeak(mp) <- FALSE
    res_1 <- do_findPeaks_MSW(mz = mz(sp1), int = intensity(sp1),
                              snthresh = 1, nearbyPeak = FALSE)
    res_2 <- findChromPeaks(od1, param = mp, return.type = "list")
    expect_equal(res_1, res_2[[1]][, colnames(res_1)])
    peakThr(mp) <- 200
    res_1 <- do_findPeaks_MSW(mz = mz(sp1), int = intensity(sp1),
                              snthresh = 1, nearbyPeak = FALSE,
                              peakThr = 200)
    res_2 <- findChromPeaks(od1, param = mp, return.type = "list")
    expect_equal(res_1, res_2[[1]][, colnames(res_1)])
    addParams(mp) <- list(forder = 2)
    res_3 <- do_findPeaks_MSW(mz = mz(sp1), int = intensity(sp1),
                              snthresh = 1, nearbyPeak = FALSE,
                              peakThr = 200, forder = 2)
    res_4 <- findChromPeaks(od1, param = mp, return.type = "list")
    expect_equal(res_3, res_4[[1]][, colnames(res_3)])
    addParams(mp) <- list(forder = 2, dorder = 1)
    res_3 <- do_findPeaks_MSW(mz = mz(sp1), int = intensity(sp1),
                              snthresh = 1, nearbyPeak = FALSE,
                              peakThr = 200, forder = 2, dorder = 1)
    res_4 <- findChromPeaks(od1, param = mp, return.type = "list")
    expect_equal(res_3, res_4[[1]][, colnames(res_3)])
    ## Compare old vs new:
    expect_equal(chromPeaks(fticr_xod)[, -ncol(chromPeaks(fticr_xod))],
                 peaks(fticr_xs))
})

test_that("featureValues,XCMSnExp works", {
    od_x <- faahko_xod
    xs <- faahko_xs
    
    p <- PeakDensityParam(sampleGroups = xs$class)
    od_x <- groupChromPeaks(od_x, param = p)

    xs <- group(xs, method = "density")

    expect_equal(unname(groupval(xs, value = "into")),
                 unname(featureValues(od_x, value = "into")))
    expect_equal(unname(groupval(xs, method = "maxint", value = "into")),
                 unname(featureValues(od_x, method = "maxint", value = "into")))
    ## Checking errors
    expect_error(featureValues(od_x, value = "bla"))    
})

test_that("groupChromPeaks,XCMSnExp,PeakDensityParam works", {
    od_x <- faahko_xod
    xs <- faahko_xs
    ## Check error if no features were found. issue #273
    pdp <- PeakDensityParam(sampleGroups = xs$class, minSamples = 30)
    expect_error(groupChromPeaks(od_x, param = pdp), "Unable to group any chromatographic peaks.")
    
    fdp <- PeakDensityParam(sampleGroups = xs$class)
    od_x <- groupChromPeaks(od_x, param = fdp)
    xs <- group(xs, method = "density")
    expect_equal(xs@groupidx, featureDefinitions(od_x)$peakidx)
    fg <- featureDefinitions(od_x)
    fg <- S4Vectors::as.matrix(fg[, -ncol(fg)])
    rownames(fg) <- NULL
    expect_equal(xs@groups, fg)
    expect_true(length(processHistory(od_x)) == 2)
    ph <- processHistory(od_x, type = .PROCSTEP.PEAK.GROUPING)[[1]]
    expect_equal(processParam(ph), fdp)
    expect_equal(rownames(featureDefinitions(od_x)),
                 .featureIDs(nrow(featureDefinitions(od_x))))
    
    fdp2 <- PeakDensityParam(sampleGroups = xs$class, binSize = 2,
                             minFraction = 0.8)
    od_x <- groupChromPeaks(od_x, param = fdp2)
    xs <- group(xs, method = "density", minfrac = 0.8, mzwid = 2)
    expect_equal(xs@groupidx, featureDefinitions(od_x)$peakidx)
    fg <- featureDefinitions(od_x)
    fg <- S4Vectors::as.matrix(fg[, -ncol(fg)])
    rownames(fg) <- NULL
    expect_equal(xs@groups, fg)
    expect_true(length(processHistory(od_x)) == 2)
    ph <- processHistory(od_x, type = .PROCSTEP.PEAK.GROUPING)[[1]]
    expect_equal(processParam(ph), fdp2)
    expect_equal(rownames(featureDefinitions(od_x)),
                 .featureIDs(nrow(featureDefinitions(od_x))))
})

test_that("groupPeaks,XCMSnExp,MzClustParam works", {
    p <- MzClustParam(sampleGroups = sampclass(fticr_xs))
    fticr_xod2 <- groupChromPeaks(fticr_xod, param = p)
    fticr_xs2 <- group(fticr_xs, method = "mzClust")
    expect_equal(fticr_xs2@groupidx, featureDefinitions(fticr_xod2)$peakidx)
    fg <- featureDefinitions(fticr_xod2)
    fg <- S4Vectors::as.matrix(fg[, -ncol(fg)])
    rownames(fg) <- NULL
    expect_equal(fticr_xs2@groups, fg)
    expect_true(length(processHistory(fticr_xod2)) == 2)
    ph <- processHistory(fticr_xod2,
                         type = .PROCSTEP.PEAK.GROUPING)[[1]]
    expect_equal(processParam(ph), p)
    expect_equal(rownames(featureDefinitions(fticr_xod2)),
                 .featureIDs(nrow(featureDefinitions(fticr_xod2))))
    p2 <- MzClustParam(sampleGroups = fticr_xs$class, absMz = 1,
                       minFraction = 0.8)
    fticr_xod2 <- groupChromPeaks(fticr_xod, param = p2)
    fticr_xs2 <- group(fticr_xs, method = "mzClust", minfrac = 0.8, mzabs = 1)
    expect_equal(fticr_xs2@groupidx, featureDefinitions(fticr_xod2)$peakidx)
    fg <- featureDefinitions(fticr_xod2)
    fg <- S4Vectors::as.matrix(fg[, -ncol(fg)])
    rownames(fg) <- NULL
    expect_equal(fticr_xs2@groups, fg)
    expect_true(length(processHistory(fticr_xod2)) == 2)
    ph <- processHistory(fticr_xod2,
                         type = .PROCSTEP.PEAK.GROUPING)[[1]]
    expect_equal(processParam(ph), p2)    
    expect_equal(rownames(featureDefinitions(fticr_xod2)),
                 .featureIDs(nrow(featureDefinitions(fticr_xod2))))
})

test_that("groupChromPeaks,XCMSnExp,NearestPeaksParam works", {
    od_x <- faahko_xod
    xs <- faahko_xs
    p <- NearestPeaksParam(sampleGroups = xs$class)
    od_x <- groupChromPeaks(od_x, param = p)
    xs <- group(xs, method = "nearest")
    expect_equal(xs@groupidx, featureDefinitions(od_x)$peakidx)
    fg <- featureDefinitions(od_x)
    fg <- S4Vectors::as.matrix(fg[, -ncol(fg)])
    rownames(fg) <- NULL
    expect_equal(xs@groups, fg)
    expect_true(length(processHistory(od_x)) == 2)
    ph <- processHistory(od_x, type = .PROCSTEP.PEAK.GROUPING)[[1]]
    expect_equal(processParam(ph), p)
    expect_equal(rownames(featureDefinitions(od_x)),
                 .featureIDs(nrow(featureDefinitions(od_x))))
    fdp2 <- NearestPeaksParam(sampleGroups = xs$class, kNN = 3)
    od_x <- groupChromPeaks(od_x, param = fdp2)
    xs <- group(xs, method = "nearest", kNN = 3)
    expect_equal(xs@groupidx, featureDefinitions(od_x)$peakidx)
    fg <- featureDefinitions(od_x)
    fg <- S4Vectors::as.matrix(fg[, -ncol(fg)])
    rownames(fg) <- NULL
    expect_equal(xs@groups, fg)
    expect_true(length(processHistory(od_x)) == 2)
    ph <- processHistory(od_x, type = .PROCSTEP.PEAK.GROUPING)[[1]]
    expect_equal(processParam(ph), fdp2)    
    expect_equal(rownames(featureDefinitions(od_x)),
                 .featureIDs(nrow(featureDefinitions(od_x))))
})

test_that("fillChromPeaks,XCMSnExp works", {
    ## No adjusted retention times
    expect_true(!.hasFilledPeaks(xod_xg))
    res <- fillChromPeaks(xod_xg)
    expect_true(.hasFilledPeaks(res))
    ph <- processHistory(res, type = .PROCSTEP.PEAK.FILLING)
    expect_true(length(ph) == 1)
    expect_equal(ph[[1]]@param, FillChromPeaksParam())
    ## Check parameter filled in featureValues (issue #157)
    expect_equal(featureValues(res, filled = FALSE), featureValues(xod_xg))
    
    ## Check if the signal corresponds to what we expect for some peaks.
    fp <- chromPeaks(res)
    fp <- fp[fp[, "is_filled"] == 1, ]
    idxs <- sample(1:nrow(fp), 5)
    for (i in idxs) {
        cfp <- fp[i, , drop = FALSE]
        tmp <- filterFile(xod_xg, file = cfp[1, "sample"])
        chr <- chromatogram(tmp, rt = cfp[1, c("rtmin", "rtmax")],
                            mz = cfp[1, c("mzmin", "mzmax")])[1, 1]
        into <- sum(intensity(chr), na.rm = TRUE) *
            (cfp[1, "rtmax"] - cfp[1, "rtmin"]) / (length(chr) - 1)
        expect_equal(unname(into), unname(cfp[1, "into"]))
    }
    ## Plot the data for some...
    if (FALSE) {
        pk_idx <- featureValues(res)[1, ]
        pks <- chromPeaks(res)[pk_idx, ]
        rtr <- c(min(pks[, "rtmin"]), max(pks[, "rtmax"]))
        rtr[1] <- rtr[1] - 10
        rtr[2] <- rtr[2] + 10
        chrs <- chromatogram(res, rt = rtr, mz = c(min(pks[, "mzmin"]),
                                                   max(pks[, "mzmax"])))[1, ]
        plot(3, 3, pch = NA, xlim = range(lapply(chrs, rtime), na.rm = TRUE),
             ylim = range(lapply(chrs, intensity), na.rm = TRUE), xlab = "rt",
             ylab = "int")
        for (i in 1:length(chrs)) {
            points(rtime(chrs[[i]]), intensity(chrs[[i]]), type = "l",
                   col = ifelse(pks[i, "is_filled"], yes = "red", no = "black"))
            abline(v = pks[i, c("rtmin", "rtmax")],
                   col = ifelse(pks[i, "is_filled"], yes = "red", no = "black"))
        }
    }

    ## Check if the results are similar that we get with findChromPeaks
    for (i in 1:length(fileNames(xod_xg))) {
        fnd_pks <- chromPeaks(xod_xg)[chromPeaks(xod_xg)[, "sample"] == i, ]
        prm <- processHistory(tmp, type ="Peak detection")[[1]]@param
        ## Extract the data for these using the internal function.
        fld_pks <- .getChromPeakData(filterFile(xod_xg, i),
                                     peakArea = fnd_pks,
                                     sample_idx = i,
                                     cn = colnames(fnd_pks))
        ## rt
        expect_true(cor(fnd_pks[, "rt"], fld_pks[, "rt"]) > 0.99)
        ## mz
        expect_true(cor(fnd_pks[, "mz"], fld_pks[, "mz"]) > 0.99)
        expect_equal(fnd_pks[, "mz"], fld_pks[, "mz"])
        ## into
        expect_true(cor(fnd_pks[, "into"], fld_pks[, "into"]) > 0.99)
        expect_equal(fnd_pks[, "into"], fld_pks[, "into"])
        ## expect_equal(fnd_pks[, "into"], fld_pks[, "into"])    
        ## maxo
        expect_equal(fnd_pks[, "maxo"], fld_pks[, "maxo"])    
        expect_equal(fnd_pks[, "maxo"], fld_pks[, "maxo"])
    }
    
    ## Check for the NAs if there is really no signal
    gv <- featureValues(res)
    feat_i <- which(is.na(gv[, 1]))
    tmp <- chromPeaks(res)[featureDefinitions(res)$peakidx[[feat_i]],
                           c("rtmin", "rtmax", "mzmin", "mzmax")]
    ## Get the intensities for the first one.
    pkArea <- apply(tmp, median, MARGIN = 2)
    chr <- chromatogram(res, rt = pkArea[1:2], mz = pkArea[3:4])[1, ]
    expect_true(all(unlist(lapply(chr, function(z) is.na(intensity(z))))))
    ## Get also the spectra:
    spctr <- spectra(filterRt(filterFile(xod_xg, file = 1), rt = pkArea[1:2]))
    mzs <- unlist(lapply(spctr, mz))
    ## No spectra for the fiven mz:
    expect_equal(sum(mzs >= pkArea[3] & mzs <= pkArea[4]), 0)
    
    ## Check increasing the expandRt and expandMz to see whether we get rid of
    ## the NA.
    res_2 <- fillChromPeaks(xod_xg, param = FillChromPeaksParam(expandMz = 1))
    ## Check if the mzrange is now indeed broader for the integrated ones.
    fp <- chromPeaks(res)
    fp <- fp[fp[, "is_filled"] == 1, ]
    fp2 <- chromPeaks(res_2)
    fp2 <- fp2[fp2[, "is_filled"] == 1, ]
    expect_equal(fp2[, "mzmax"] - fp2[, "mzmin"],
                 2 * (fp[, "mzmax"] - fp[, "mzmin"]))
    
    res_2 <- fillChromPeaks(xod_xg, param = FillChromPeaksParam(expandRt = 1))
    ## Check if the mzrange is now indeed broader for the integrated ones.
    fp <- chromPeaks(res)
    fp <- fp[fp[, "is_filled"] == 1, ]
    fp2 <- chromPeaks(res_2)
    fp2 <- fp2[fp2[, "is_filled"] == 1, ]
    expect_equal(fp2[, "rtmax"] - fp2[, "rtmin"],
                 2 * (fp[, "rtmax"] - fp[, "rtmin"]))
    ## Check using ppm
    res_2 <- fillChromPeaks(xod_xg, param = FillChromPeaksParam(ppm = 40,
                                                                expandMz = 5,
                                                                expandRt = 2))
    expect_true(all(!is.na(rowSums(featureValues(res_2)))))
    ## Drop them.
    res_rem <- dropFilledChromPeaks(res)
    expect_true(!.hasFilledPeaks(res_rem))
    expect_equal(res_rem, xod_xg)
    ## Drop feature definitions from res -> also filled peaks should be dropped.
    res_rem <- dropFeatureDefinitions(res)
    expect_true(!.hasFilledPeaks(res_rem))
    expect_true(!any(chromPeaks(res_rem)[, "is_filled"] == 1))
    expect_equal(res_rem, xod_x)
    
    ## With adjusted rtime.
    res_2 <- fillChromPeaks(xod_xgrg)
    ## Check if the signal corresponds to what we expect for some peaks.
    fp <- chromPeaks(res_2)
    fp <- fp[fp[, "is_filled"] == 1, ]
    ## These have to be different from before!
    fp_raw <- chromPeaks(res)
    fp_raw <- fp_raw[fp_raw[, "is_filled"] == 1, ]
    expect_true(all(fp_raw[, "rt"] != fp[, "rt"]))
    expect_true(all(fp_raw[, "rtmin"] != fp[, "rtmin"]))
    expect_true(all(fp_raw[, "rtmax"] != fp[, "rtmax"]))
    expect_equal(fp_raw[, "mz"], fp[, "mz"])
    expect_equal(fp_raw[, "mzmin"], fp[, "mzmin"])
    expect_equal(fp_raw[, "mzmax"], fp[, "mzmax"])
    ## Values are expected to be different, but still correlated!
    expect_true(all(fp_raw[, "into"] != fp[, "into"]))
    expect_true(cor(fp_raw[, "into"], fp[, "into"]) > 0.99)
    ## Check if we can get the same data using the provided range.
    ## Use the .rawMat function
    first <- filterFile(xod_xgrg, file = 1, keepAdjustedRtime = TRUE)
    spctr <- spectra(first)
    mzs <- lapply(spctr, mz)
    vps <- lengths(mzs)
    ints <- unlist(lapply(spctr, intensity), use.names = FALSE)
    mzs <- unlist(mzs, use.names = FALSE)
    rtim <- rtime(first)
    idx <- which(fp[, "sample"] == 1)
    for (i in idx) {
        mtx <- .rawMat(mz = mzs, int = ints, scantime = rtim,
                       valsPerSpect = vps,
                       rtrange = fp[i, c("rtmin", "rtmax")],
                       mzrange = fp[i, c("mzmin", "mzmax")])
        into <- sum(mtx[, 3], na.rm = TRUE) *
            ((fp[i, "rtmax"] - fp[i, "rtmin"]) /
             (sum(rtim >= fp[i, "rtmin"] & rtim <= fp[i, "rtmax"]) - 1))
        expect_equal(unname(into), unname(fp[i, "into"]))
    }
    ## Drop them.
    res_rem <- dropFilledChromPeaks(res_2)
    expect_true(!.hasFilledPeaks(res_rem))
    expect_equal(res_rem, xod_xgrg)
})

test_that("fillChromPeaks,XCMSnExp with MSW works", {
    p <- MzClustParam()
    fticr_xodg <- groupChromPeaks(fticr_xod, param = p)
    expect_error(res <- fillChromPeaks(fticr_xod))
    res <- fillChromPeaks(fticr_xodg)
    
    ## Got a signal for all of em.
    expect_true(!any(is.na(featureValues(res))))
    ## 1) Compare with what I get for xcmsSet.
    tmp_x <- fticr_xs
    tmp_x <- group(tmp_x, method = "mzClust")
    tmp_x <- fillPeaks(tmp_x, method = "MSW")
    ## Compare
    expect_equal(unname(groupval(tmp_x)), unname(featureValues(res)))
    expect_equal(unname(groupval(tmp_x, value = "maxo")),
                 unname(featureValues(res, value = "maxo")))
    expect_equal(unname(groupval(tmp_x, value = "into")),
                 unname(featureValues(res, value = "into")))
    expect_equal(unname(groupval(tmp_x, value = "mz")),
                 unname(featureValues(res, value = "mz")))
    expect_equal(unname(groupval(tmp_x, value = "mzmin")),
                 unname(featureValues(res, value = "mzmin")))
    expect_equal(unname(groupval(tmp_x, value = "mzmax")),
                 unname(featureValues(res, value = "mzmax")))
    ## OK
    ## 2) Check if the fillChromPeaks returns same/similar data than the
    ##    findChromPeaks does:
    fdef <- featureDefinitions(fticr_xodg)
    pkArea <- do.call(
        rbind,
        lapply(
            fdef$peakidx, function(z) {
                tmp <- chromPeaks(fticr_xodg)[z, c("rtmin", "rtmax",
                                                   "mzmin", "mzmax"),
                                              drop = FALSE]
                pa <- c(median(tmp[, 1]), median(tmp[, 2]),
                        median(tmp[, 3]), median(tmp[, 4]))
                return(pa)
            }
        ))
    colnames(pkArea) <- c("rtmin", "rtmax", "mzmin", "mzmax")
    pkArea <- cbind(group_idx = 1:nrow(pkArea), pkArea,
                    mzmed = fdef$mzmed)
    ## Get peak data for all peaks in the first file
    allPks <- .getMSWPeakData(filterFile(fticr_xodg, file = 1),
                              peakArea = pkArea,
                              sample_idx = 1,
                              cn = colnames(chromPeaks(fticr_xodg)))
    curP <- chromPeaks(res)[chromPeaks(res)[, "sample"] == 1, ]
    curP <- curP[order(curP[, "mz"]), ]
    expect_equal(allPks[, "mz"], curP[, "mz"])
    expect_equal(allPks[, "maxo"], curP[, "maxo"])
    expect_true(cor(allPks[, "into"], curP[, "into"]) > 0.99) ## Not exactly the
    ## same but highly similar.
})

test_that("fillChromPeaks,XCMSnExp with matchedFilter works", {
    tmp <- findChromPeaks(faahko_od, param = MatchedFilterParam())
    sg <- rep(1, length(fileNames(tmp)))
    tmp <- groupChromPeaks(tmp, param = PeakDensityParam(sampleGroups = sg))

    tmp_filled <- fillChromPeaks(tmp)
    expect_true(sum(is.na(featureValues(tmp_filled))) <
                sum(is.na(featureValues(tmp))))
    nas <- is.na(featureValues(tmp)[, 1]) | is.na(featureValues(tmp)[, 2])
    expect_true(cor(featureValues(tmp, value = "into")[!nas, 1],
                    featureValues(tmp, value = "into")[!nas, 2]) > 0.97)
    expect_true(cor(featureValues(tmp_filled, value = "into")[, 1],
                    featureValues(tmp_filled, value = "into")[, 2],
                    use = "complete.obs") > 0.97)

    ## Check signal generation for already found peaks.
    for (i in 1:length(fileNames(tmp))) {
        fnd_pks <- chromPeaks(tmp)[chromPeaks(tmp)[, "sample"] == i, ]
        prm <- processHistory(tmp, type ="Peak detection")[[1]]@param
        ## Extract the data for these using the internal function.
        fld_pks <- .getChromPeakData_matchedFilter(filterFile(tmp, i),
                                                   peakArea = fnd_pks,
                                                   sample_idx = i,
                                                   param = prm,
                                                   cn = colnames(fnd_pks))
        ## rt can not be the same, since for fillChromPeaks it is the rt of the
        ## maximum signal and for findChromPeaks it is the rt of the apex of the
        ## filtered/fitted peak.
        expect_true(cor(fnd_pks[, "rt"], fld_pks[, "rt"]) > 0.99)
        ## mz: also not the same; most likely due to slightly different binning.
        diffs <- fnd_pks[, "mz"] - fld_pks[, "mz"]
        expect_true(max(diffs) < 1e-4)
        ## into
        expect_equal(fnd_pks[, "into"], fld_pks[, "into"])    
        ## maxo
        expect_equal(fnd_pks[, "maxo"], fld_pks[, "maxo"])    
    }

    ## modify fillChromPeaks settings.
    tmp_fld_2 <- fillChromPeaks(
        tmp, param = FillChromPeaksParam(ppm = 40, expandRt = 1))
    expect_true(sum(is.na(featureValues(tmp_filled))) <
                sum(is.na(featureValues(tmp))))
    expect_true(sum(is.na(featureValues(tmp_fld_2))) <
                sum(is.na(featureValues(tmp_filled))))
    nas <- is.na(featureValues(tmp)[, 1]) | is.na(featureValues(tmp)[, 2])
    expect_true(cor(featureValues(tmp_fld_2, value = "into")[, 1],
                    featureValues(tmp_fld_2, value = "into")[, 2],
                    use = "complete.obs") > 0.97)
})
