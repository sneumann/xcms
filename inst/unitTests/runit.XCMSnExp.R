## tests related to the new XCMSnExp object.

## od_x <- faahko_od
## xod_x <- faahko_xod
## xod_xg <- groupChromPeaks(xod_x, param = PeakDensityParam())
## xod_xgr <- adjustRtime(xod_xg, param = PeakGroupsParam(span = 0.4))
## xod_xgrg <- groupChromPeaks(xod_xgr, param = PeakDensityParam())

xs <- faahko_xs

xs_2 <- group(xs)
suppressWarnings(
    xs_2 <- retcor(xs_2)
)
xs_2 <- group(xs_2)
od_fa <- faahko_od

## Just checking if we can create an empty objects - got strange cases in which
## a newly created object was not empty.
.checkCreationOfEmptyObject <- function() {
    x <- new("XCMSnExp")
    checkTrue(!hasAdjustedRtime(x))
    checkTrue(!hasFeatures(x))
    checkTrue(!hasChromPeaks(x))
}

test_XCMSnExp_class <- function() {
    .checkCreationOfEmptyObject()
    ## Basic contructor.
    x <- new("XCMSnExp")
    x@.processHistory <- list("a")
    checkException(validObject(x))
    x@.processHistory <- list(xcms:::ProcessHistory())
    checkTrue(validObject(x))
    xod <- as(od_fa, "XCMSnExp")
    checkTrue(validObject(xod))
    ## MsFeatureData error: environment is locked
    checkException(xod@msFeatureData$chromPeaks <- 3)
    checkException(xod@msFeatureData$bla <- 4)
    checkTrue(validObject(xod))

    .checkCreationOfEmptyObject()
    ## xod@msFeatureData$chromPeaks <- xs_2@peaks
    ## checkTrue(validObject(xod))
    ## xod@msFeatureData$chromPeaks[1, "sample"] <- 40
    ## checkException(validObject(xod))
    ## xod@msFeatureData$chromPeaks[1, "sample"] <- 3
    ## xod@msFeatureData$adjustedRtime <- xs_2@rt$corrected
    ## checkTrue(validObject(xod))
    ## xod@msFeatureData$adjustedRtime[[2]] <- 1:4
    ## checkException(validObject(xod))
    ## xod@msFeatureData$adjustedRtime <- xs_2@rt$corrected[1:2]
    ## checkException(validObject(xod))
    ## That code above resulted in the strange effect that newly created objects
    ## were not empty.
}

test_XCMSnExp_rtime <- function() {
    rts <- rtime(faahko_od)
    rts_2 <- rtime(od_x)
    checkEquals(rts, rts_2)
    ## Test with bySample.
    rts_3 <- rtime(xod_x, bySample = TRUE)
    checkEquals(rts_3, split(rts, f = fromFile(faahko_od)))
    ## Check if rtimes are correctly ordered for bySample
    rts_4 <- rtime(filterFile(faahko_od, file = 2))
    checkEquals(rts_4, rts_3[[2]])
    rts_4 <- rtime(filterFile(faahko_od, file = 3))
    checkEquals(rts_4, rts_3[[3]])
    ## Compare with the values we get from an xcmsSet:
    rtx <- faahko_xs@rt$raw
    checkEquals(unlist(rtx, use.names = FALSE),
                unlist(rtime(faahko_xod, bySample = TRUE), use.names = FALSE))
}

test_XCMSnExp_mz <- function() {
    mzs <- mz(faahko_od)
    ## The check below has to work, since we're calling the mz,OnDiskMSnExp.
    ## mzs_2 <- mz(od_x)
    ## checkEquals(mzs, mzs_2)
    mzs_2 <- mz(xod_x, bySample = TRUE)
    tmp <- split(mzs, fromFile(faahko_od))
    checkEquals(lapply(tmp, unlist, use.names = FALSE), mzs_2)
    ## Check if mz are correctly ordered for bySample
    mzs_3 <- mz(filterFile(faahko_od, file = 2))
    checkEquals(unlist(mzs_3, use.names = FALSE), mzs_2[[2]])
}

test_XCMSnExp_intensity <- function() {
    ints <- intensity(faahko_od)
    ## The check below has to work, since we're calling the intensity,OnDiskMSnExp.
    ## ints_2 <- intensity(od_x)
    ## checkEquals(ints, ints_2)
    ints_2 <- intensity(xod_x, bySample = TRUE)
    tmp <- split(ints, fromFile(faahko_od))
    checkEquals(lapply(tmp, unlist, use.names = FALSE), ints_2)
    ## Check if mz are correctly ordered for bySample
    ints_3 <- intensity(filterFile(faahko_od, file = 2))
    checkEquals(unlist(ints_3, use.names = FALSE), ints_2[[2]])
}

test_XCMSnExp_spectra <- function() {
    xod <- as(faahko_od, "XCMSnExp")
    res <- spectra(xod)
    res_2 <- spectra(xod, bySample = TRUE)
    checkEquals(split(res, fromFile(xod)), res_2)
    ## xod_x
    tmp <- filterRt(xod_x, rt = c(2700, 2900))
    res <- spectra(tmp)
    rts <- unlist(lapply(res, rtime))
    checkEquals(rts, rtime(tmp))
    ## Check with adjusted retention times.
    tmp2 <- filterRt(xod_xgr, rt = c(2700, 2900))
    res2 <- spectra(tmp2)
    rts2 <- unlist(lapply(res2, rtime))
    checkEquals(rts2, rtime(tmp2))
    ## Now do it on one file:
    tmp <- filterFile(xod_x, file = 2)
    res <- spectra(tmp)
    checkEquals(rtime(tmp), unlist(lapply(res, rtime)))
    tmp2 <- filterFile(xod_xgr, file = 2, keepAdjustedRtime = TRUE)
    res2 <- spectra(tmp2)
    checkEquals(rtime(tmp2), unlist(lapply(res2, rtime)))
    checkTrue(sum(unlist(lapply(res2, rtime)) ==
                  unlist(lapply(res, rtime))) < length(rtime(tmp)) / 4)
    res3 <- spectra(tmp2, adjusted = FALSE)
    checkEquals(res, res3)    
    ## adjusted rt
    tmp <- filterFile(xod_xgr, file = 2, keepAdjustedRtime = TRUE)
    checkTrue(hasAdjustedRtime(tmp))
    checkTrue(is.character(all.equal(rtime(tmp, adjusted = FALSE),
                                     adjustedRtime(tmp))))
    res <- spectra(tmp)
    checkEquals(rtime(tmp), unlist(lapply(res, rtime)))
    res <- unlist(spectrapply(tmp, FUN = function(x) {rtime(x)}))
    checkEquals(res, adjustedRtime(tmp))
}

test_XCMSnExp_class_accessors <- function() {
    .checkCreationOfEmptyObject()
    ## Filling with data...
    xod <- as(od_fa, "XCMSnExp")
    ## peaks 
    checkTrue(!hasChromPeaks(xod))
    chromPeaks(xod) <- xs_2@peaks
    checkTrue(hasChromPeaks(xod))
    checkEquals(chromPeaks(xod), xs_2@peaks)
    checkException(chromPeaks(xod) <- 4)
    tmp <- chromPeaks(xod, bySample = TRUE)
    checkTrue(length(tmp) == length(fileNames(xod)))
    tmp <- do.call(rbind, tmp)
    rownames(tmp) <- NULL
    checkEquals(tmp, chromPeaks(xod))
    ## chromPeaks with rt
    all_pks <- chromPeaks(xod_x)
    pks <- chromPeaks(xod_x, rt = c(2000, 2600), type = "within")
    checkTrue(nrow(pks) < nrow(all_pks))
    checkTrue(all(pks[, "rtmin"] >= 2000 & pks[, "rtmax"] <= 2600))
    pks <- chromPeaks(xod_x, rt = c(2000, 2600), bySample = TRUE,
                      type = "within")
    checkTrue(nrow(pks[[2]]) == 0)
    pks <- chromPeaks(xod_x, rt = c(2000, 2600), type = "any")
    checkTrue(all(pks[, "rtmax"] >= 2000 & pks[, "rtmin"] <= 2600))
    pks <- chromPeaks(xod_x, rt = c(2000, 2200))
    checkTrue(nrow(pks) == 0)
    pks <- chromPeaks(xod_x, rt = c(2000, 2200), bySample = TRUE)
    checkTrue(all(lengths(pks) == 0))
    ## chromPeaks with mz
    pks <- chromPeaks(xod_x, mz = c(280, 281), type = "within")
    checkTrue(all(pks[, "mzmin"] >= 280 & pks[, "mzmax"] <= 281))
    pks <- chromPeaks(xod_x, mz = c(280, 281), bySample = TRUE, type = "within")
    checkTrue(nrow(pks[[1]]) == 0)
    checkTrue(nrow(pks[[3]]) == 0)
    checkTrue(nrow(pks[[2]]) == 1)
    pks <- chromPeaks(xod_x, mz = c(280, 300), bySample = FALSE, type = "within")
    checkTrue(all(pks[, "mzmin"] >= 280 & pks[, "mzmax"] <= 300))
    pks <- chromPeaks(xod_x, mz = c(280, 300), bySample = FALSE, type = "any")
    checkTrue(all(pks[, "mzmax"] >= 280 & pks[, "mzmin"] <= 300))
    pks <- chromPeaks(xod_x, mz = c(200, 210), bySample = FALSE)
    checkTrue(nrow(pks) == 0)
    pks <- chromPeaks(xod_x, mz = c(200, 210), bySample = TRUE)
    checkTrue(all(lengths(pks) == 0))
    ## chromPeaks with both
    pks <- chromPeaks(xod_x, mz = c(280, 300), rt = c(3000, 3300),
                      type = "within")
    checkTrue(all(pks[, "mzmin"] >= 280 & pks[, "mzmax"] <= 300))
    checkTrue(all(pks[, "rtmin"] >= 3000 & pks[, "rtmax"] <= 3300))
    pks <- chromPeaks(xod_x, mz = c(280, 300), rt = c(3000, 3300),
                      type = "any")
    checkTrue(all(pks[, "mzmax"] >= 280 & pks[, "mzmin"] <= 300))
    checkTrue(all(pks[, "rtmax"] >= 3000 & pks[, "rtmin"] <= 3300))
    ## Wrong assignments.
    pks <- xs_2@peaks
    pks[1, "sample"] <- 40
    checkException(chromPeaks(xod) <- pks)
    ## featureDefinitions
    checkTrue(!hasFeatures(xod))
    library(S4Vectors)
    fd <- DataFrame(xs_2@groups)
    fd$peakidx <- xs_2@groupidx
    featureDefinitions(xod) <- fd
    checkTrue(hasChromPeaks(xod))
    checkTrue(hasFeatures(xod))
    checkEquals(featureDefinitions(xod), fd)
    ## featureDefinitions with mz and/or rt range:
    feat_def <- featureDefinitions(xod_xg)
    ## Within
    mzr <- c(300, 330)
    keep_mz <- feat_def$mzmin > mzr[1] & feat_def$mzmax < mzr[2]
    checkEquals(featureDefinitions(xod_xg, mz = mzr, type = "within"),
                feat_def[keep_mz, ])
    rtr <- c(3000, 3800)
    keep_rt <- feat_def$rtmin > rtr[1] & feat_def$rtmax < rtr[2]
    checkEquals(featureDefinitions(xod_xg, rt = rtr, type = "within"),
                feat_def[keep_rt, ])
    checkEquals(featureDefinitions(xod_xg, rt = rtr, mz = mzr, type = "within"),
                feat_def[keep_rt & keep_mz, ])
    ## Any
    mzr <- range(featureDefinitions(xod_xg)[2, "mzmed"])
    keep_mz <- feat_def$mzmax >= mzr[1] & feat_def$mzmin <= mzr[2]
    checkEquals(featureDefinitions(xod_xg, mz = mzr, type = "any"),
                feat_def[keep_mz, , drop = FALSE])
    rtr <- range(3420.006)
    keep_rt <- feat_def$rtmax >= rtr[1] & feat_def$rtmin <= rtr[2]
    checkTrue(nrow(featureDefinitions(xod_xg, rt = rtr, type = "within")) !=
              nrow(featureDefinitions(xod_xg, rt = rtr, type = "any")))
    checkEquals(featureDefinitions(xod_xg, rt = rtr, type = "any"),
                feat_def[keep_rt, , drop = FALSE])
    checkEquals(featureDefinitions(xod_xg, rt = rtr, mz = mzr, type = "any"),
                feat_def[keep_rt & keep_mz, , drop = FALSE])
    
    ## adjustedRtime
    checkTrue(!hasAdjustedRtime(xod))
    checkTrue(hasAdjustedRtime(xod_r))
    suppressWarnings(checkEquals(adjustedRtime(xod), NULL))
    checkEquals(rtime(xod_r, adjusted = FALSE), rtime(xod))
    checkEquals(rtime(xod_r), adjustedRtime(xod_r))
    checkTrue(is.character(all.equal(rtime(xod_r), rtime(xod_r, adjusted = FALSE))))
    ## Indirect test that the ordering of the adjusted retention times matches
    ## ordering of rtime.
    ## From MSnbase version >= 2.3.9 values are ordered first by file then by
    ## spectrum.
    if (grepl("^F", names(rtime(xod_r)[1]))) {
        rts_by_sample <- adjustedRtime(xod_r, bySample = TRUE)
        rts <- adjustedRtime(xod_r)
        checkEquals(unname(rts_by_sample[[2]]),
                    unname(rts[grep(names(rts), pattern = "F2")]))
        checkEquals(unname(unlist(rts_by_sample)),
                    unname(rts))
    }
    xod2 <- xod_r
    ## Wrong assignments.
    checkException(adjustedRtime(xod2) <- xs_2@rt$corrected[1:2])
    ## bracket subset
    tmp <- xod2[1]
    checkTrue(length(tmp[[1]]) == 1)
    checkTrue(length(xod2[[1]]) == 1)
    .checkCreationOfEmptyObject()
}

test_XCMSnExp_findChromPeaks <- function() {
    ## Call findChromPeaks on an XCMSnExp
    tmp <- findChromPeaks(xod_x, param = CentWaveParam(noise = 10000,
                                                       snthresh = 40))
    checkEquals(chromPeaks(tmp), chromPeaks(xod_x))
    ## Check that it works also on adjusted retention times:
    tmp <- findChromPeaks(xod_r, param = CentWaveParam(noise = 10000,
                                                       snthresh = 40))
    checkTrue(hasAdjustedRtime(tmp))
    checkEquals(
        length(processHistory(tmp, type = xcms:::.PROCSTEP.RTIME.CORRECTION)),1)
    checkTrue(sum(chromPeaks(tmp)[, "rt"] != chromPeaks(xod_x)[, "rt"]) >
              ncol(chromPeaks(tmp)))
    tmp_sub <- filterFile(xod_r, file = 1, keepAdjustedRtime = TRUE)
    checkEquals(rtime(tmp_sub, adjusted = TRUE),
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
    checkEquals(res_2, pks)
    ## Second try:
    tmp <- findChromPeaks(xod_xgrg, param = CentWaveParam(noise = 10000,
                                                          snthresh = 40))
    checkTrue(hasAdjustedRtime(tmp))
    checkEquals(
        length(processHistory(tmp, type = xcms:::.PROCSTEP.RTIME.CORRECTION)),1)
    checkTrue(sum(chromPeaks(tmp)[, "rt"] != chromPeaks(xod_x)[, "rt"]) >
              ncol(chromPeaks(tmp)))
    tmp_sub <- filterFile(xod_xgrg, file = 3, keepAdjustedRtime = TRUE)
    checkEquals(unname(rtime(tmp_sub, adjusted = TRUE)),
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
    checkEquals(res_2, pks)
}


test_XCMSnExp_processHistory <- function() {
    ph <- xcms:::ProcessHistory(fileIndex. = 2, info. = "For file 2")
    ph_2 <- xcms:::ProcessHistory(fileIndex. = 1:2, info. = "For files 1 to 2")
    xod <- as(od_fa, "XCMSnExp")
    xod@.processHistory <- list(ph, ph_2)
    checkEquals(processHistory(xod), list(ph, ph_2))
    checkEquals(processHistory(xod, fileIndex = 2), list(ph, ph_2))
    checkEquals(processHistory(xod, fileIndex = 1), list(ph_2))
    checkException(processHistory(xod, fileIndex = 5))

    ph_3 <- xcms:::XProcessHistory(fileIndex = 1, param = CentWaveParam())
    xod <- xcms:::addProcessHistory(xod, ph_3)
    checkEquals(length(processHistory(xod)), 3)
    checkEquals(processHistory(xod)[[3]], ph_3)
    checkEquals(processHistory(xod, fileIndex = 1), list(ph_2, ph_3))
    checkTrue(validObject(xod))
}

test_XCMSnExp_droppers <- function() {
    ## How are the drop functions expected to work?
    .checkCreationOfEmptyObject()
    type_feat_det <- xcms:::.PROCSTEP.PEAK.DETECTION
    type_feat_algn <- xcms:::.PROCSTEP.PEAK.GROUPING
    type_rt_adj <- xcms:::.PROCSTEP.RTIME.CORRECTION
    ## Perform alignment.
    ## xod_xg <- groupChromPeaks(xod_x, param = PeakDensityParam())
    checkTrue(hasFeatures(xod_xg))
    checkTrue(hasChromPeaks(xod_x))
    checkTrue(hasChromPeaks(xod_xg))
    checkTrue(!hasAdjustedRtime(xod_xg))
    checkTrue(length(processHistory(xod_xg, type = type_feat_algn)) == 1)
    ## Retention time adjustment.
    ## xod_xgr <- adjustRtime(xod_xg, param = PeakGroupsParam(span = 1))
    checkTrue(hasChromPeaks(xod_xgr))
    checkTrue(length(processHistory(xod_xgr, type = type_feat_det)) == 1)
    checkTrue(!hasFeatures(xod_xgr))  ## These should have been removed
    checkTrue(length(processHistory(xod_xgr, type = type_feat_algn)) == 1)
    checkTrue(hasAdjustedRtime(xod_xgr))
    checkTrue(length(processHistory(xod_xgr, type = type_rt_adj)) == 1)
    ## Most of the retention times are different
    checkTrue(sum(chromPeaks(xod_xgr)[, "rt"] != chromPeaks(xod_x)[, "rt"]) >
              nrow(chromPeaks(xod_x)) / 2)
    checkTrue(sum(rtime(xod_xgr) == rtime(xod_xg)) < length(rtime(xod_xg) / 2))
    ## Alignment after retention time adjustment.
    ## xod_xgrg <- groupChromPeaks(xod_xgr, param = PeakDensityParam())
    checkTrue(hasChromPeaks(xod_xgrg))
    checkEquals(chromPeaks(xod_xgrg), chromPeaks(xod_xgr))
    checkTrue(hasAdjustedRtime(xod_xgrg))
    checkEquals(rtime(xod_xgrg), rtime(xod_xgr))
    checkEquals(rtime(xod_xgrg, adjusted = FALSE), rtime(od_x))
    checkTrue(length(processHistory(xod_xgr, type = type_feat_algn)) == 1)
    checkTrue(hasFeatures(xod_xgrg))
    checkTrue(length(processHistory(xod_xgrg, type = type_feat_algn)) == 2)
       
    ## 1) dropDetectedFeatures: delete all process history steps and all data.
    res <- dropChromPeaks(xod_x)
    checkTrue(!hasChromPeaks(res))
    checkTrue(length(processHistory(res, type = type_feat_det)) == 0)
    checkTrue(!hasFeatures(res))
    checkTrue(length(processHistory(res, type = type_feat_algn)) == 0)
    checkTrue(!hasAdjustedRtime(res))
    checkTrue(length(processHistory(res, type = type_rt_adj)) == 0)
    checkEquals(rtime(res), rtime(od_x))
    ##
    res <- dropChromPeaks(xod_xg)
    checkTrue(!hasChromPeaks(res))
    checkTrue(length(processHistory(res, type = type_feat_det)) == 0)
    checkTrue(!hasFeatures(res))
    checkTrue(length(processHistory(res, type = type_feat_algn)) == 0)
    checkTrue(!hasAdjustedRtime(res))
    checkTrue(length(processHistory(res, type = type_rt_adj)) == 0)
    checkEquals(rtime(res), rtime(od_x))
    ##
    res <- dropChromPeaks(xod_xgr)
    checkTrue(!hasChromPeaks(res))
    checkTrue(length(processHistory(res, type = type_feat_det)) == 0)
    checkTrue(!hasFeatures(res))
    checkTrue(length(processHistory(res, type = type_feat_algn)) == 0)
    checkTrue(!hasAdjustedRtime(res))
    checkTrue(length(processHistory(res, type = type_rt_adj)) == 0)
    checkEquals(rtime(res), rtime(od_x))
    res <- dropChromPeaks(xod_xgr, keepAdjustedRtime = TRUE)
    checkTrue(hasAdjustedRtime(res))
    checkTrue(length(processHistory(res, type = type_rt_adj)) == 1)
    checkEquals(rtime(res), rtime(xod_xgr))
    checkTrue(!hasChromPeaks(res))
    checkTrue(length(processHistory(res, type = type_feat_det)) == 0)
    checkTrue(!hasFeatures(res))
    checkTrue(length(processHistory(res, type = type_feat_algn)) == 0)
    ##
    res <- dropChromPeaks(xod_xgrg)
    checkTrue(!hasChromPeaks(res))
    checkTrue(length(processHistory(res, type = type_feat_det)) == 0)
    checkTrue(!hasFeatures(res))
    checkTrue(length(processHistory(res, type = type_feat_algn)) == 0)
    checkTrue(!hasAdjustedRtime(res))
    checkTrue(length(processHistory(res, type = type_rt_adj)) == 0)
    checkEquals(rtime(res), rtime(od_x))
    res <- dropChromPeaks(xod_xgrg, keepAdjustedRtime = TRUE)
    checkTrue(hasAdjustedRtime(res))
    checkTrue(length(processHistory(res, type = type_rt_adj)) == 1)
    checkEquals(rtime(res), rtime(xod_xgr))
    checkTrue(!hasChromPeaks(res))
    checkTrue(length(processHistory(res, type = type_feat_det)) == 0)
    checkTrue(!hasFeatures(res))
    checkTrue(length(processHistory(res, type = type_feat_algn)) == 0)
    
    ## 2) dropFeatureDefinitions:
    ##    a) drop the feature groups and the latest related process history
    ##    b) if retention time correction was performed AFTER the latest feature
    ##       grouping, drop also the retention time correction and all related
    ##       process histories.
    res <- dropFeatureDefinitions(xod_xg)
    checkEquals(res, xod_x)
    checkTrue(hasChromPeaks(res))
    checkTrue(length(processHistory(res, type = type_feat_det)) == 1)
    checkTrue(!hasFeatures(res))
    checkTrue(length(processHistory(res, type = type_feat_algn)) == 0)
    checkTrue(!hasAdjustedRtime(res))
    checkTrue(length(processHistory(res, type = type_rt_adj)) == 0)
    checkEquals(rtime(res), rtime(od_x))
    ## No feature groups - so there is nothing that this function does here.
    res <- dropFeatureDefinitions(xod_xgr)
    checkEquals(res, xod_xgr)
    checkTrue(hasChromPeaks(res))
    checkTrue(length(processHistory(res, type = type_feat_det)) == 1)
    checkTrue(!hasFeatures(res))
    checkTrue(length(processHistory(res, type = type_feat_algn)) == 1)
    checkTrue(hasAdjustedRtime(res))
    checkTrue(length(processHistory(res, type = type_rt_adj)) == 1)
    ## Remove the latest ones.
    res <- dropFeatureDefinitions(xod_xgrg)
    checkEquals(res, xod_xgr)
    checkTrue(hasChromPeaks(res))
    checkTrue(length(processHistory(res, type = type_feat_det)) == 1)
    checkTrue(!hasFeatures(res))
    checkTrue(length(processHistory(res, type = type_feat_algn)) == 1)
    checkTrue(hasAdjustedRtime(res))
    checkTrue(length(processHistory(res, type = type_rt_adj)) == 1)
    checkEquals(rtime(res, adjusted = FALSE), rtime(od_x))
    checkEquals(rtime(res, adjusted = TRUE), rtime(xod_xgr))

    ## 3) dropAdjustedRtime:
    ##    a) drop the retention time adjustment and related process histories
    ##    b) if grouping has been performed AFTER retention time correction,
    ##       drop the feature alignment and all related process histories.
    ##    c) if grouping has been performed BEFORE retention time correction,
    ##       do nothing.
    res <- dropAdjustedRtime(xod_xg)
    checkEquals(res, xod_xg)
    ## This drops also the process history for alignment.
    res <- dropAdjustedRtime(xod_xgr)
    checkTrue(hasChromPeaks(res))
    checkTrue(length(processHistory(res, type = type_feat_det)) == 1)
    checkTrue(!hasFeatures(res))
    checkTrue(length(processHistory(res, type = type_feat_algn)) == 0)
    checkTrue(!hasAdjustedRtime(res))
    checkTrue(length(processHistory(res, type = type_rt_adj)) == 0)
    checkEquals(chromPeaks(res), chromPeaks(xod_x))
    checkEquals(res, xod_x)
    checkEquals(rtime(res), rtime(xod_x))
    checkEquals(rtime(res), rtime(xod_xgr, adjusted = FALSE))
    ## This drops also the feature alignment performed later.
    res <- dropAdjustedRtime(xod_xgrg)
    checkTrue(hasChromPeaks(res))
    checkTrue(length(processHistory(res, type = type_feat_det)) == 1)
    checkTrue(!hasFeatures(res))
    checkTrue(length(processHistory(res, type = type_feat_algn)) == 0)
    checkTrue(!hasAdjustedRtime(res))
    checkTrue(length(processHistory(res, type = type_rt_adj)) == 0)
    checkEquals(chromPeaks(res), chromPeaks(xod_x))
    checkEquals(res, xod_x)
    checkEquals(rtime(res), rtime(xod_xgrg, adjusted = FALSE))
    
    .checkCreationOfEmptyObject()
}


## Testing all inherited methods for XCMSnExp classes that perform data
## manipulations - for all of these we have to ensure that eventual preprocessing
## results are removed.
test_XCMSnExp_inherited_methods <- function() {
    .checkCreationOfEmptyObject()
    ## [
    tmp_1 <- od_fa[1:10]
    suppressWarnings(
        tmp_2 <- xod_x[1:10]
    )
    checkTrue(length(processHistory(tmp_2)) == 0)
    checkTrue(!hasChromPeaks(tmp_2))
    tmp_1@processingData <- new("MSnProcess")
    tmp_2@processingData <- new("MSnProcess")
    checkEquals(tmp_1, as(tmp_2, "OnDiskMSnExp"))
    checkException(xod_r[1, 1])
    idxs <- c(1432, 1621, 2492, 3001, 3013)
    tmp <- xod_r[idxs]
    checkTrue(length(tmp) == length(idxs))
    checkEquals(mz(xod_r)[idxs], mz(tmp))
    checkTrue(hasAdjustedRtime(xod_r) != hasAdjustedRtime(tmp))
    ## keeping adjusted retention times:
    tmp <- xod_r[idxs, keepAdjustedRtime = TRUE]
    checkTrue(hasAdjustedRtime(tmp))
    checkEquals(rtime(xod_r)[idxs], rtime(tmp))
    ## Same with object containing also peaks and features
    tmp <- xod_xgrg[idxs]
    checkTrue(!hasAdjustedRtime(tmp))
    checkTrue(!hasChromPeaks(tmp))
    checkTrue(!hasFeatures(tmp))
    tmp <- xod_xgrg[idxs, keepAdjusted = TRUE]
    checkTrue(hasAdjustedRtime(tmp))
    checkEquals(rtime(xod_xgrg)[idxs], rtime(tmp))
    checkTrue(length(processHistory(tmp)) == 1)
    
    ## [[
    spct <- xod_x[[13]]
    checkTrue(is(spct, "Spectrum1"))
    checkEquals(rtime(spct), unname(rtime(xod_x)[13]))
    checkEquals(mz(spct), mz(xod_x)[[13]])
    ## Have to ensure that, if x has adjusted retention times, that these are
    ## reported in the Spectrum.
    spct <- xod_r[[13]]
    checkEquals(rtime(spct), unname(rtime(xod_r, adjusted = TRUE)[13]))
    checkTrue(rtime(spct) != rtime(xod_r, adjusted = FALSE)[13])
    
    ## bin
    tmp_1 <- bin(od_fa)
    suppressWarnings(
        tmp_2 <- bin(xod_x)
    )
    checkTrue(length(processHistory(tmp_2)) == 0)
    checkTrue(!hasChromPeaks(tmp_2))
    tmp_1@processingData <- new("MSnProcess")
    tmp_2@processingData <- new("MSnProcess")
    checkEquals(tmp_1, as(tmp_2, "OnDiskMSnExp"))
    ## clean
    tmp_1 <- clean(od_fa)
    suppressWarnings(
        tmp_2 <- clean(xod_x)
    )
    checkTrue(length(processHistory(tmp_2)) == 0)
    checkTrue(!hasChromPeaks(tmp_2))
    tmp_1@processingData <- new("MSnProcess")
    tmp_2@processingData <- new("MSnProcess")
    checkEquals(tmp_1, as(tmp_2, "OnDiskMSnExp"))
    ## filterAcquisitionNum
    tmp_1 <- filterAcquisitionNum(od_fa)
    suppressWarnings(
        tmp_2 <- filterAcquisitionNum(xod_x)
    )
    checkTrue(length(tmp_2[[1]]) > 0)
    checkTrue(length(processHistory(tmp_2)) == 0)
    checkTrue(!hasChromPeaks(tmp_2))
    tmp_1@processingData <- new("MSnProcess")
    tmp_2@processingData <- new("MSnProcess")
    checkEquals(tmp_1, as(tmp_2, "OnDiskMSnExp"))
    ## filterMsLevel
    tmp_1 <- filterMsLevel(od_fa)
    suppressWarnings(
        tmp_2 <- filterMsLevel(xod_x)
    )
    checkEquals(tmp_2, xod_x)
    suppressWarnings(
        checkEquals(length(filterMsLevel(xod_x, msLevel = 2)), 0)
    )
    ## If we've got adjusted retention times, keep them.
    suppressWarnings(tmp_1 <- filterMsLevel(xod_xgr, msLevel = 1))
    checkTrue(hasAdjustedRtime(tmp_1))
    checkEquals(rtime(tmp_1), rtime(xod_xgr)) # adjusted rt present
    suppressWarnings(
        tmp_1 <- filterMsLevel(xod_xgrg, msLevel = 1, keepAdjustedRtime = FALSE)
    )
    checkTrue(!hasAdjustedRtime(tmp_1))
    checkEquals(rtime(tmp_1), rtime(xod_xgr, adjusted = FALSE))
    ## normalize
    tmp_1 <- normalize(od_fa)
    suppressWarnings(
        tmp_2 <- normalize(xod_x)
    )
    checkTrue(length(processHistory(tmp_2)) == 0)
    checkTrue(!hasChromPeaks(tmp_2))
    tmp_1@processingData <- new("MSnProcess")
    tmp_2@processingData <- new("MSnProcess")
    checkEquals(tmp_1, as(tmp_2, "OnDiskMSnExp"))
    ## pickPeaks
    tmp_1 <- pickPeaks(od_fa)
    suppressWarnings(
        tmp_2 <- pickPeaks(xod_x)
    )
    checkTrue(length(processHistory(tmp_2)) == 0)
    checkTrue(!hasChromPeaks(tmp_2))
    tmp_1@processingData <- new("MSnProcess")
    tmp_2@processingData <- new("MSnProcess")
    checkEquals(tmp_1, as(tmp_2, "OnDiskMSnExp"))
    ## removePeaks
    tmp_1 <- removePeaks(od_fa)
    suppressWarnings(
        tmp_2 <- removePeaks(xod_x)
    )
    checkTrue(length(processHistory(tmp_2)) == 0)
    checkTrue(!hasChromPeaks(tmp_2))
    tmp_1@processingData <- new("MSnProcess")
    tmp_2@processingData <- new("MSnProcess")
    checkEquals(tmp_1, as(tmp_2, "OnDiskMSnExp"))
    ## smooth
    tmp_1 <- smooth(od_fa)
    suppressWarnings(
        tmp_2 <- smooth(xod_x)
    )
    checkTrue(length(processHistory(tmp_2)) == 0)
    checkTrue(!hasChromPeaks(tmp_2))
    tmp_1@processingData <- new("MSnProcess")
    tmp_2@processingData <- new("MSnProcess")
    checkEquals(tmp_1, as(tmp_2, "OnDiskMSnExp"))
    .checkCreationOfEmptyObject()
}

## Test XCMSnExp filter methods.
test_XCMSnExp_filterFile <- function() {
    ## filterFile
    tmp <- filterFile(xod_x, file = 2)
    checkException(tmp@msFeatureData$bla <- 3)
    checkTrue(!hasAdjustedRtime(tmp))
    checkTrue(!hasFeatures(tmp))
    checkTrue(all(chromPeaks(tmp)[, "sample"] == 1))
    checkEquals(chromPeaks(tmp)[, -(ncol(chromPeaks(tmp)) - 1)],
                chromPeaks(xod_x)[chromPeaks(xod_x)[, "sample"] == 2,
                                  -(ncol(chromPeaks(xod_x)) - 1)])
    checkEquals(fileIndex(processHistory(tmp)[[1]]), 1)
    ## check with other index.
    tmp <- filterFile(xod_x, file = c(1, 3))
    checkTrue(length(tmp[[1]]) == 1)
    checkTrue(!hasAdjustedRtime(tmp))
    checkTrue(!hasFeatures(tmp))
    checkTrue(all(chromPeaks(tmp)[, "sample"] %in% c(1, 2)))
    a <- chromPeaks(tmp)
    b <- chromPeaks(xod_x)
    checkEquals(a[, -(ncol(a) - 1)],
                b[b[, "sample"] %in% c(1, 3), -(ncol(b) - 1)])
    checkEquals(fileIndex(processHistory(tmp)[[1]]), c(1, 2))

    ## Errors
    checkException(filterFile(xod_x, file = 5))
    checkException(filterFile(xod_x, file = 1:5))

    ## Little mockup to check correctness of Process history.
    od_2 <- xod_x
    od_2 <- xcms:::addProcessHistory(
                       od_2,
                       xcms:::ProcessHistory(
                                  type = xcms:::.PROCSTEP.RTIME.CORRECTION))
    od_2 <- xcms:::addProcessHistory(
                       od_2,
                       xcms:::ProcessHistory(type = xcms:::.PROCSTEP.UNKNOWN,
                                             fileIndex = 2,
                                             info. = "I should be here"))
    od_2 <- xcms:::addProcessHistory(
                       od_2,
                       xcms:::ProcessHistory(type = xcms:::.PROCSTEP.UNKNOWN,
                                             fileIndex = 1, info. = "EEEEEE"))

    tmp <- filterFile(od_2, file = 2)
    ph <- processHistory(tmp)
    checkTrue(length(ph) == 2)
    checkEquals(processType(ph[[2]]), xcms:::.PROCSTEP.UNKNOWN)
    b <- unlist(lapply(ph, function(z) {
        processInfo(z) == "I should be here"
    }))
    checkTrue(any(b))
    b <- unlist(lapply(ph, function(z) {
        processInfo(z) == "EEEEEE"
    }))
    checkTrue(!any(b))
    ## Do filterFile on xod_xg
    res <- filterFile(xod_xg, file = 2)
    checkTrue(hasChromPeaks(res))
    checkTrue(!hasAdjustedRtime(res))
    checkTrue(!hasFeatures(res))
    tmp <- chromPeaks(xod_xg)
    checkEquals(chromPeaks(res)[, -(ncol(tmp) - 1)],
                tmp[tmp[, "sample"] == 2, -(ncol(tmp) - 1)])
    checkEquals(rtime(res), rtime(xod_xg, bySample = TRUE)[[2]])
    ## Do filterFile on xod_xgr
    ## Should remove adjusted rts and revert the original peak rts.
    res <- filterFile(xod_xgr, file = 2)
    checkTrue(hasChromPeaks(res))
    tmp <- chromPeaks(xod_xg)
    checkEquals(chromPeaks(res)[, -(ncol(tmp) - 1)],
                tmp[tmp[, "sample"] == 2, -(ncol(tmp) - 1)])
    checkEquals(rtime(res), rtime(xod_xg, bySample = TRUE)[[2]])
    checkTrue(!hasAdjustedRtime(res))
    checkTrue(!hasFeatures(res))
    checkTrue(length(processHistory(res)) == 1)
    checkEquals(processType(processHistory(res)[[1]]), "Peak detection")
    ## The same but keep the adjusted retention times.
    res <- filterFile(xod_xgr, file = 2, keepAdjustedRtime = TRUE)
    checkTrue(hasChromPeaks(res))
    tmp <- chromPeaks(xod_xgr)
    checkEquals(chromPeaks(res)[, -(ncol(tmp) - 1)],
                tmp[tmp[, "sample"] == 2, -(ncol(tmp) - 1)])
    ## has to be different from the ones in xod_x
    tmp <- chromPeaks(xod_x)
    checkTrue(sum(chromPeaks(res)[, "rt"] == tmp[tmp[, "sample"] == 2, "rt"]) <
              nrow(tmp) / 4)
    checkEquals(rtime(res), rtime(xod_xgr, bySample = TRUE)[[2]])
    checkEquals(adjustedRtime(res), adjustedRtime(xod_xgr, bySample = TRUE)[[2]])
    checkTrue(hasAdjustedRtime(res))
    checkTrue(!hasFeatures(res))
    checkTrue(length(processHistory(res)) == 3)
    checkEquals(processType(processHistory(res)[[1]]), "Peak detection")
    checkEquals(processType(processHistory(res)[[2]]), "Peak grouping")
    checkEquals(processType(processHistory(res)[[3]]), "Retention time correction")
    ## Do filterFile on xod_xgrg
    res <- filterFile(xod_xgrg, file = c(1, 3))
    checkTrue(hasChromPeaks(res))
    tmp <- chromPeaks(xod_x)
    checkEquals(chromPeaks(res)[, -(ncol(tmp) - 1)],
                tmp[tmp[, "sample"] %in% c(1, 3), -(ncol(tmp) - 1)])
    checkEquals(unname(rtime(res, bySample = TRUE)),
                unname(rtime(xod_xg, bySample = TRUE)[c(1, 3)]))
    checkTrue(!hasAdjustedRtime(res))
    checkTrue(!hasFeatures(res))
    checkTrue(length(processHistory(res)) == 1)
    checkEquals(processType(processHistory(res)[[1]]), "Peak detection")
    ## keep adjusted rtime
    res <- filterFile(xod_xgrg, file = c(1, 3), keepAdjustedRtime = TRUE)
    checkTrue(hasChromPeaks(res))
    tmp <- chromPeaks(xod_xgr)
    checkEquals(chromPeaks(res)[, -(ncol(tmp) - 1)],
                tmp[tmp[, "sample"] %in% c(1, 3), -(ncol(tmp) - 1)])
    ## has to be different from the ones in xod_x
    tmp <- chromPeaks(xod_x)
    checkTrue(sum(chromPeaks(res)[, "rt"] == tmp[tmp[, "sample"] %in% c(1, 3), "rt"]) <
              nrow(tmp) / 4)
    checkEquals(rtime(res, bySample = TRUE),
                rtime(xod_xgr, bySample = TRUE)[c(1, 3)])
    checkEquals(adjustedRtime(res, bySample = TRUE),
                adjustedRtime(xod_xgr, bySample = TRUE)[c(1, 3)])
    checkTrue(hasAdjustedRtime(res))
    checkTrue(!hasFeatures(res))
    checkTrue(length(processHistory(res)) == 3)
    checkEquals(processType(processHistory(res)[[1]]), "Peak detection")
    checkEquals(processType(processHistory(res)[[2]]), "Peak grouping")
    checkEquals(processType(processHistory(res)[[3]]), "Retention time correction")
}

test_XCMSnExp_filterMz <- function() {
    ## od_x2 <- od_x
    ## new_e <- xcms:::.copy_env(od_x2@msFeatureData)
    ## xcms:::adjustedRtime(od_x2) <- xs_2@rt$corrected
    ## library(S4Vectors)
    ## fd <- DataFrame(xs_2@groups)
    ## fd$peakidx <- xs_2@groupidx
    ## featureDefinitions(od_x2) <- fd

    ## subset on xod_x
    res <- filterMz(xod_x, mz = c(300, 400))
    suppressWarnings(
        checkTrue(length(res[[1]]) == 1)
    )
    checkTrue(length(res@spectraProcessingQueue) == 1)
    checkTrue(hasChromPeaks(res))
    checkTrue(all(chromPeaks(res)[, "mz"] >= 300 & chromPeaks(res)[, "mz"] <= 400))
    checkTrue(nrow(chromPeaks(res)) < nrow(chromPeaks(xod_x)))
    idx <- which(chromPeaks(xod_x)[, "mzmin"] >= 300 &
                 chromPeaks(xod_x)[, "mzmax"] <= 400)
    checkEquals(chromPeaks(res), chromPeaks(xod_x)[idx, ])
    checkTrue(!hasAdjustedRtime(res))
    checkTrue(!hasFeatures(res))
    ## subset on xod_xg
    res <- filterMz(xod_xg, mz = c(300, 400))
    checkTrue(validObject(res))
    suppressWarnings(
        checkTrue(length(res[[1]]) == 1)
    )
    checkTrue(length(res@spectraProcessingQueue) == 1)
    checkTrue(hasChromPeaks(res))
    checkTrue(all(chromPeaks(res)[, "mz"] >= 300 & chromPeaks(res)[, "mz"] <= 400))
    checkTrue(nrow(chromPeaks(res)) < nrow(chromPeaks(xod_x)))
    idx <- which(chromPeaks(xod_xg)[, "mzmin"] >= 300 &
                 chromPeaks(xod_xg)[, "mzmax"] <= 400)
    checkEquals(chromPeaks(res), chromPeaks(xod_xg)[idx, ])
    checkTrue(!hasAdjustedRtime(res))
    checkTrue(hasFeatures(res))
    checkTrue(nrow(featureDefinitions(res)) < nrow(featureDefinitions(xod_xg)))
    checkTrue(all(featureDefinitions(res)[, "mzmed"] >= 300 &
                  featureDefinitions(res)[, "mzmed"] <= 400))
    checkTrue(all(featureDefinitions(res)[, "mzmin"] >= 300 &
                  featureDefinitions(res)[, "mzmin"] <= 400))
    checkTrue(all(featureDefinitions(res)[, "mzmax"] >= 300 &
                  featureDefinitions(res)[, "mzmax"] <= 400))
    ## subset on xod_xgr
    ## o keep chromPeaks
    ## o keep adjusted rtime
    res <- filterMz(xod_xgr, mz = c(300, 400))
    checkTrue(validObject(res))
    suppressWarnings(
        checkTrue(length(res[[1]]) == 1)
    )
    checkTrue(length(res@spectraProcessingQueue) == 1)
    checkTrue(hasChromPeaks(res))
    checkTrue(all(chromPeaks(res)[, "mz"] >= 300 & chromPeaks(res)[, "mz"] <= 400))
    checkTrue(nrow(chromPeaks(res)) < nrow(chromPeaks(xod_x)))
    idx <- which(chromPeaks(xod_xgr)[, "mzmin"] >= 300 &
                 chromPeaks(xod_xgr)[, "mzmax"] <= 400)
    checkEquals(chromPeaks(res), chromPeaks(xod_xgr)[idx, ])
    checkTrue(hasAdjustedRtime(res))
    checkEquals(adjustedRtime(res), adjustedRtime(xod_xgr))
    checkTrue(!hasFeatures(res))
    checkTrue(length(processHistory(res)) == 3)
    ## subset xod_xgrg
    res <- filterMz(xod_xgrg, mz = c(300, 400))
    checkTrue(validObject(res))
    suppressWarnings(
        checkTrue(length(res[[1]]) == 1)
    )
    checkTrue(length(res@spectraProcessingQueue) == 1)
    checkTrue(hasChromPeaks(res))
    checkTrue(all(chromPeaks(res)[, "mz"] >= 300 & chromPeaks(res)[, "mz"] <= 400))
    checkTrue(nrow(chromPeaks(res)) < nrow(chromPeaks(xod_xgrg)))
    idx <- which(chromPeaks(xod_xgrg)[, "mzmin"] >= 300 &
                 chromPeaks(xod_xgrg)[, "mzmax"] <= 400)
    checkEquals(chromPeaks(res), chromPeaks(xod_xgrg)[idx, ])
    checkTrue(hasAdjustedRtime(res))
    checkEquals(adjustedRtime(res), adjustedRtime(xod_xgrg))
    checkTrue(hasFeatures(res))
    checkTrue(nrow(featureDefinitions(res)) < nrow(featureDefinitions(xod_xgrg)))
    checkTrue(all(featureDefinitions(res)[, "mzmed"] >= 300 &
                  featureDefinitions(res)[, "mzmed"] <= 400))
    checkTrue(all(featureDefinitions(res)[, "mzmin"] >= 300 &
                  featureDefinitions(res)[, "mzmin"] <= 400))
    checkTrue(all(featureDefinitions(res)[, "mzmax"] >= 300 &
                  featureDefinitions(res)[, "mzmax"] <= 400))
    ## With groups - no groups within this range
    mzr <- c(120, 130)
    res <- filterMz(xod_xg, mz = mzr)
    checkTrue(!hasFeatures(res))
    checkTrue(hasChromPeaks(res))
    checkTrue(all(chromPeaks(res)[, "mz"] >= 120 & chromPeaks(res)[, "mz"] <= 130))
    res <- filterMz(xod_xgrg, mz = mzr)
    checkTrue(!hasFeatures(res))
    checkTrue(hasChromPeaks(res))
    checkTrue(all(chromPeaks(res)[, "mz"] >= 120 & chromPeaks(res)[, "mz"] <= 130))    
}

test_XCMSnExp_filterRt <- function() {

    ## xod_x
    res <- filterRt(xod_x, rt = c(2700, 2900))
    ## Check if the object is OK:
    checkEquals(pData(res), pData(xod_x))
    spct <- spectra(res)
    checkTrue(length(spct) > 0)    
    ## MsFeatureData has to be locked!
    checkException(res@msFeatureData$bla <- 3)
    ## Retention time has to be within the range.
    checkTrue(all(rtime(res) >= 2700 & rtime(res) <= 2900))
    rtm <- unlist(lapply(spct, rtime))
    checkTrue(all(rtm >= 2700 & rtm <= 2900))
    ## peaks have to be within the range.
    checkTrue(all(chromPeaks(res)[, "rt"] >= 2700 &
                  chromPeaks(res)[, "rt"] <= 2900))
    are_within <- chromPeaks(xod_x)[, "rt"] >= 2700 &
        chromPeaks(xod_x)[, "rt"] <= 2900
    checkEquals(chromPeaks(res), chromPeaks(xod_x)[are_within,])
    ## Have a feature detection process history.
    checkEquals(length(processHistory(res)), 1)
    checkEquals(processType(processHistory(res)[[1]]),
                xcms:::.PROCSTEP.PEAK.DETECTION)
    ## filter such that we keep some spectra but no chromPeaks:
    res <- filterRt(xod_x, rt = c(4200, 4400))
    checkTrue(all(rtime(res) >= 4200 & rtime(res) <= 4400))
    checkTrue(!hasChromPeaks(res))
    checkTrue(length(processHistory(res)) == 0)
    ## No rt
    res <- filterRt(xod_x, rt = c(10, 20))
    checkTrue(length(res) == 0)

    ## xod_xg
    ## o keep also the feature groups that are within the window.
    res <- filterRt(xod_xg, rt = c(2700, 2900))
    checkTrue(all(rtime(res) >= 2700 & rtime(res) <= 2900))
    checkEquals(hasChromPeaks(res), hasChromPeaks(xod_xg))
    checkTrue(all(chromPeaks(res)[, "rt"] >= 2700 &
                  chromPeaks(res)[, "rt"] <= 2900))
    are_within <- chromPeaks(xod_x)[, "rt"] >= 2700 &
        chromPeaks(xod_x)[, "rt"] <= 2900
    checkEquals(chromPeaks(res), chromPeaks(xod_xg)[are_within,])
    checkTrue(!hasAdjustedRtime(res))
    checkTrue(hasFeatures(res))
    checkTrue(all(featureDefinitions(res)$rtmed >= 2700 &
                                    featureDefinitions(res)$rtmed <= 2900))
    checkTrue(nrow(featureDefinitions(res)) < nrow(featureDefinitions(xod_xg)))
    checkTrue(length(processHistory(res)) == 2)
    checkTrue(length(processHistory(res, type = "Peak detection")) == 1)
    checkTrue(length(processHistory(res, type = "Peak grouping")) == 1)
    ## All feature idx have to match.
    checkTrue(all(unlist(featureDefinitions(res)$peakidx) %in%
                  1:nrow(chromPeaks(res))))
    ## Filter such that we don't have any chromPeaks.
    res <- filterRt(xod_xg, rt = c(4200, 4400))
    checkTrue(all(rtime(res) >= 4200 & rtime(res) <= 4400))
    checkTrue(!hasChromPeaks(res))
    checkTrue(!hasFeatures(res))
    checkTrue(length(processHistory(res)) == 0)
    ## No rt
    res <- filterRt(xod_xg, rt = c(10, 20))
    checkTrue(length(res) == 0)

    ## xod_xgr
    res <- filterRt(xod_xgr, rt = c(2700, 2900))
    checkTrue(all(rtime(res) >= 2700 & rtime(res) <= 2900))
    checkEquals(hasChromPeaks(res), hasChromPeaks(xod_xg))
    checkTrue(all(chromPeaks(res)[, "rt"] >= 2700 &
                  chromPeaks(res)[, "rt"] <= 2900))
    are_within <- chromPeaks(xod_xgr)[, "rt"] >= 2700 &
        chromPeaks(xod_xgr)[, "rt"] <= 2900
    checkEquals(chromPeaks(res), chromPeaks(xod_xgr)[are_within,])
    checkTrue(hasAdjustedRtime(res))
    checkTrue(all(adjustedRtime(res) >= 2700 & adjustedRtime(res) <= 2900))
    checkTrue(!all(rtime(res, adjusted = FALSE) >= 2700 &
                   rtime(res, adjusted = FALSE) <= 2900))
    checkTrue(!hasFeatures(res))
    checkTrue(length(processHistory(res, type = "Peak detection")) == 1)
    checkTrue(length(processHistory(res, type = "Peak grouping")) == 1)
    checkTrue(length(processHistory(res, type = "Retention time correction")) == 1)
    ## Filter such that we don't have any chromPeaks.
    res <- filterRt(xod_xgr, rt = c(4200, 4400), adjusted = TRUE)
    checkTrue(hasAdjustedRtime(res))
    checkTrue(all(adjustedRtime(res) >= 4200 & adjustedRtime(res) <= 4400))
    checkTrue(all(rtime(res) >= 4200 & rtime(res) <= 4400))
    checkTrue(!all(rtime(res, adjusted = FALSE) >= 4200 &
                   rtime(res, adjusted = FALSE) <= 4400))
    checkTrue(!hasChromPeaks(res))
    checkTrue(!hasFeatures(res))
    checkTrue(length(processHistory(res)) == 1)
    checkTrue(length(processHistory(res, type = "Retention time correction")) == 1)
    ## No rt
    res <- filterRt(xod_xgr, rt = c(10, 20))
    checkTrue(length(res) == 0)
    ## filter using raw rt
    res <- filterRt(xod_xgr, rt = c(2700, 2900), adjusted = FALSE)
    checkTrue(!all(rtime(res) >= 2700 & rtime(res) <= 2900))
    checkEquals(hasChromPeaks(res), hasChromPeaks(xod_xg))
    checkTrue(all(chromPeaks(res)[, "rt"] >= 2700 &
                  chromPeaks(res)[, "rt"] <= 2900))
    are_within <- chromPeaks(xod_xgr)[, "rt"] >= 2700 &
        chromPeaks(xod_xgr)[, "rt"] <= 2900
    checkEquals(chromPeaks(res), chromPeaks(xod_xgr)[are_within,])
    checkTrue(hasAdjustedRtime(res))
    checkTrue(!all(adjustedRtime(res) >= 2700 & adjustedRtime(res) <= 2900))
    checkTrue(all(rtime(res, adjusted = FALSE) >= 2700 &
                  rtime(res, adjusted = FALSE) <= 2900))
    checkTrue(!hasFeatures(res))

    ## xod_xgrg
    res <- filterRt(xod_xgrg, rt = c(2700, 2900))
    checkTrue(all(rtime(res) >= 2700 & rtime(res) <= 2900))
    checkEquals(hasChromPeaks(res), hasChromPeaks(xod_xg))
    checkTrue(all(chromPeaks(res)[, "rt"] >= 2700 &
                  chromPeaks(res)[, "rt"] <= 2900))
    are_within <- chromPeaks(xod_xgrg)[, "rt"] >= 2700 &
        chromPeaks(xod_xgr)[, "rt"] <= 2900
    checkEquals(chromPeaks(res), chromPeaks(xod_xgrg)[are_within,])
    checkTrue(hasAdjustedRtime(res))
    checkTrue(all(adjustedRtime(res) >= 2700 & adjustedRtime(res) <= 2900))
    checkTrue(!all(rtime(res, adjusted = FALSE) >= 2700 &
                   rtime(res, adjusted = FALSE) <= 2900))
    checkTrue(length(processHistory(res, type = "Peak detection")) == 1)
    checkTrue(length(processHistory(res, type = "Peak grouping")) == 2)
    checkTrue(length(processHistory(res, type = "Retention time correction")) == 1)
    checkTrue(hasFeatures(res))
    checkTrue(all(featureDefinitions(res)$rtmed >= 2700 &
                                    featureDefinitions(res)$rtmed <= 2900))
    ## Filter such that we don't have any chromPeaks.
    res <- filterRt(xod_xgrg, rt = c(4200, 4400), adjusted = TRUE)
    checkTrue(hasAdjustedRtime(res))
    checkTrue(all(adjustedRtime(res) >= 4200 & adjustedRtime(res) <= 4400))
    checkTrue(all(rtime(res) >= 4200 & rtime(res) <= 4400))
    checkTrue(!all(rtime(res, adjusted = FALSE) >= 4200 &
                   rtime(res, adjusted = FALSE) <= 4400))
    checkTrue(!hasChromPeaks(res))
    checkTrue(!hasFeatures(res))
    checkTrue(length(processHistory(res)) == 1)
    checkTrue(length(processHistory(res, type = "Retention time correction")) == 1)
    ## No rt
    res <- filterRt(xod_xgrg, rt = c(10, 20))
    checkTrue(length(res) == 0)
    ## filter using raw rt
    res <- filterRt(xod_xgrg, rt = c(2700, 2900), adjusted = FALSE)
    checkTrue(!all(rtime(res) >= 2700 & rtime(res) <= 2900))
    checkEquals(hasChromPeaks(res), hasChromPeaks(xod_xg))
    checkTrue(all(chromPeaks(res)[, "rt"] >= 2700 &
                  chromPeaks(res)[, "rt"] <= 2900))
    are_within <- chromPeaks(xod_xgrg)[, "rt"] >= 2700 &
        chromPeaks(xod_xgrg)[, "rt"] <= 2900
    checkEquals(chromPeaks(res), chromPeaks(xod_xgrg)[are_within,])
    checkTrue(hasAdjustedRtime(res))
    checkTrue(!all(adjustedRtime(res) >= 2700 & adjustedRtime(res) <= 2900))
    checkTrue(all(rtime(res, adjusted = FALSE) >= 2700 &
                  rtime(res, adjusted = FALSE) <= 2900))
    checkTrue(hasFeatures(res))
    checkTrue(all(featureDefinitions(res)$rtmed >= 2700 &
                                    featureDefinitions(res)$rtmed <= 2900))
}

## Test the coercion method.
test_as_XCMSnExp_xcmsSet <- function() {
    od_x <- faahko_xod
    res <- xcms:::.XCMSnExp2xcmsSet(od_x)
    res <- as(od_x, "xcmsSet")
    ## Results should be the same as in xs.
    checkEquals(res@peaks, chromPeaks(od_x))
    checkEquals(res@.processHistory, processHistory(od_x))
    checkEquals(phenoData(res), pData(od_x))
    checkEquals(filepaths(res), fileNames(od_x))
    checkEquals(res@rt$raw, res@rt$corrected)
    checkEquals(res@rt$raw, rtime(od_x, bySample = TRUE))
    checkEquals(profMethod(res), "bin")
    checkEquals(profStep(res), 0.1)
    ## Can we further process this?
    sampclass(res) <- rep("K", 3)
    res <- group.density(res, minfrac = 0.5)
    ## res <- fillPeaks(res)

    ## Add groups.
    od_2 <- groupChromPeaks(
        od_x,
        param = PeakDensityParam(sampleGroups =rep(1, length(fileNames(od_x)))))
    checkEquals(unname(featureDefinitions(od_2)$peakidx), groupidx(res))

    ## rt correction
    od_3 <- adjustRtime(od_2, param = PeakGroupsParam(minFraction = 1,
                                                         span = 0.4))
    ## With groups.
    res <- as(od_2, "xcmsSet")
    ftDef <- featureDefinitions(od_2)[, -ncol(featureDefinitions(od_2))]
    ftDef <- S4Vectors::as.matrix(ftDef)
    rownames(ftDef) <- NULL
    checkEquals(res@groups, ftDef)
    checkEquals(res@groupidx, unname(featureDefinitions(od_2)$peakidx))

    ## With adjusted retention time.
    res_2 <- retcor.peakgroups(res, missing = 0, span = 0.4)
    res <- as(od_3, "xcmsSet")
    checkTrue(any(unlist(res@rt$raw) != unlist(res@rt$corrected)))
    checkEquals(res@rt$corrected, res_2@rt$corrected)
    checkEquals(chromPeaks(od_3), peaks(res))
    checkEquals(peaks(res_2), peaks(res))
    
    ## Test with different binning methods:
    ## o binlin
    mfp <- MatchedFilterParam(impute = "lin", binSize = 3)
    od_2 <- od_x
    xcms:::processParam(od_2@.processHistory[[1]]) <- mfp
    res <- as(od_2, "xcmsSet")
    checkEquals(profStep(res), 3)
    checkEquals(profMethod(res), "binlin")
    ## o binlinbase
    mfp <- MatchedFilterParam(impute = "linbase", binSize = 2)
    xcms:::processParam(od_2@.processHistory[[1]]) <- mfp
    suppressWarnings(
        res <- as(od_2, "xcmsSet")
    )
    checkEquals(profStep(res), 2)
    checkEquals(profMethod(res), "binlinbase")
}

test_MsFeatureData_class_validation <- function() {
    fd <- new("MsFeatureData")
    library(S4Vectors)
    ## Check error for wrong elements.
    fd$a <- 5
    checkTrue(!is.logical(xcms:::validateMsFeatureData(fd)))
    rm("a", envir = fd)
    ## Check chromPeaks
    fd$chromPeaks <- 4
    checkTrue(!is.logical(xcms:::validateMsFeatureData(fd)))
    fdm <- matrix(ncol = 3, nrow = 5)
    colnames(fdm) <- c("a", "b", "sample")
    fd$chromPeaks <- fdm
    checkTrue(!is.logical(xcms:::validateMsFeatureData(fd)))
    rm("chromPeaks", envir = fd)
    ## featureDefinitions
    fd$chromPeaks <- xs_2@peaks
    fd$featureDefinitions <- 4
    checkTrue(!is.logical(xcms:::validateMsFeatureData(fd)))
    fg <- DataFrame(fdm)
    fd$featureDefinitions <- fg
    checkTrue(!is.logical(xcms:::validateMsFeatureData(fd)))
    fg <- DataFrame(xs_2@groups)
    fg$peakidx <- xs_2@groupidx
    fg_2 <- fg
    fg_2$mzmin <- "a"
    fd$featureDefinitions <- fg_2
    checkTrue(!is.logical(xcms:::validateMsFeatureData(fd)))
    fg_2 <- fg
    fg_2$peakidx[[1]] <- c(50000, 3)
    fd$featureDefinitions <- fg_2
    checkTrue(!is.logical(xcms:::validateMsFeatureData(fd)))
    ## adjustedRtime
    fd$featureDefinitions <- fg
    fd$adjustedRtime <- 4
    checkTrue(!is.logical(xcms:::validateMsFeatureData(fd)))
    fd$adjustedRtime <- list(1:5, "b")
    checkTrue(!is.logical(xcms:::validateMsFeatureData(fd)))
    ## Now check that we pass if we put all correct data into the object:
    fd <- new("MsFeatureData")
    fd$chromPeaks <- xs_2@peaks
    checkTrue(length(xcms:::validateMsFeatureData(fd)) == 0)
    fd$adjustedRtime <- xs_2@rt$corrected
    checkTrue(length(xcms:::validateMsFeatureData(fd)) == 0)
    fg <- DataFrame(xs_2@groups)
    fg$peakidx <- xs_2@groupidx
    checkTrue(length(xcms:::validateMsFeatureData(fd)) == 0)
}

test_MsFeatureData_class_accessors <- function() {
    fd <- new("MsFeatureData")
    library(S4Vectors)
    checkTrue(!hasChromPeaks(fd))
    checkTrue(!hasAdjustedRtime(fd))
    checkTrue(!hasFeatures(fd))
    suppressWarnings(checkEquals(chromPeaks(fd), NULL))
    suppressWarnings(checkEquals(featureDefinitions(fd), NULL))
    suppressWarnings(checkEquals(adjustedRtime(fd), NULL))
    ## chromPeaks
    chromPeaks(fd) <- xs_2@peaks
    checkTrue(hasChromPeaks(fd))
    checkEquals(chromPeaks(fd), xs_2@peaks)
    ## featureDefinitions
    fg <- DataFrame(xs_2@groups)
    fg$peakidx <- xs_2@groupidx
    featureDefinitions(fd) <- fg
    checkTrue(hasFeatures(fd))
    checkEquals(featureDefinitions(fd), fg)
    ## adjustedRtime
    adjustedRtime(fd) <- xs_2@rt$corrected
    checkTrue(hasAdjustedRtime(fd))
    checkEquals(adjustedRtime(fd), xs_2@rt$corrected)
}


## Test extraction of chromatograms.
test_chromatogram <- function() {
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
    checkTrue(all(rtime(res[1, 1]) >= rtr[1]))
    checkTrue(all(rtime(res[1, 1]) <= rtr[2]))
    checkTrue(all(rtime(res[1, 2]) >= rtr[1]))
    checkTrue(all(rtime(res[1, 2]) <= rtr[2]))
    tmp <- filterRt(filterFile(xod_x, file = 2), rt = rtr)
    checkEquals(rtime(tmp), rtime(res[1, 2]))
    ints <- spectrapply(tmp, function(z) return(max(intensity(z))))
    checkEquals(unlist(ints), intensity(res[1, 2]))
    ## Check names
    checkEquals(names(rtime(res[1, 1])), names(intensity(res[1, 1])))
    ## Assure we get the same with an OnDiskMSnExp and grouped XCMSnExp
    res_2 <- chromatogram(filterFile(od_x, file = c(1, 2)),
                          aggregationFun = "max", rt = rtr)
    checkEquals(res, res_2)
    res_3 <- chromatogram(filterFile(xod_xg, file = c(1, 2)),
                          aggregationFun = "max", rt = rtr)
    checkEquals(res, res_3)
    
    ## XCMSnExp: with mzrange and rtrange:
    mzr <- c(120, 130)
    tmp <- filterMz(xod_xg, mz = mzr)
    featureDefinitions(tmp)
    tmp <- filterRt(xod_xg, rt = rtr)
    featureDefinitions(tmp)
    res_2 <- chromatogram(xod_xg, rt = rtr, mz = mzr)
    ##

    ## XCMSnExp with adjusted rtime
    ## SEE runit.Chromatogram.R
}

test_signal_integration <- function() {
    ## Testing the signal integration of peaks.
    ## For centWave
    tmp <- xod_xgrg
    rtr <- chromPeaks(tmp)[1, c("rtmin", "rtmax")]
    mzr <- chromPeaks(tmp)[1, c("mzmin", "mzmax")]
    chr <- chromatogram(tmp, rt = rtr, mz = mzr)
    pkInt <- sum(intensity(chr[1, 1]) *
                 ((rtr[2] - rtr[1]) / (length(chr[1, 1]) - 1)))
    checkEquals(pkInt, unname(chromPeaks(tmp)[1, "into"]))

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
        checkEquals(unname(pkI), unname(chromPeaks(tmp)[i, "into"]))
    }
    ## pkI2 <- xcms:::.getPeakInt2(tmp, chromPeaks(tmp)[idxs, , drop = FALSE])
    ## checkEquals(unname(pkI2), unname(chromPeaks(tmp)[idxs, "into"]))
    
    ## Now for matchedfilter.
    tmp <- findChromPeaks(filterFile(od_x, 2), param = MatchedFilterParam())
    rtr <- chromPeaks(tmp)[1, c("rtmin", "rtmax")]
    mzr <- chromPeaks(tmp)[1, c("mzmin", "mzmax")]
    chr <- chromatogram(tmp, rt = rtr, mz = mzr)
    pkInt <- sum(intensity(chr[1, 1]) *
                 ((rtr[2] - rtr[1]) / (length(chr[1, 1]) - 1)))
    chromPeaks(tmp)[1, "into"]
    checkEquals(pkInt, unname(chromPeaks(tmp)[1, "into"]))
    idxs <- sample(1:nrow(chromPeaks(tmp)), 5)
    ## idxs <- 1:nrow(chromPeaks(tmp))
    for (i in idxs) {
        rtr <- chromPeaks(tmp)[i, c("rtmin", "rtmax")]
        mzr <- chromPeaks(tmp)[i, c("mzmin", "mzmax")]
        chr <- chromatogram(tmp, rt = rtr, mz = mzr)[1, 1]
        ints <- intensity(chr)
        pkI <- sum(ints, na.rm = TRUE) * ((rtr[2] - rtr[1]) / (length(ints) - 1))
        ## cat(" ", chromPeaks(tmp)[i, "into"], " - ", pkI, "\n")
        checkEquals(unname(pkI), unname(chromPeaks(tmp)[i, "into"]))
    }
    ## pkI2 <- xcms:::.getPeakInt2(tmp, chromPeaks(tmp)[idxs, , drop = FALSE])
    ## checkEquals(unname(pkI2), unname(chromPeaks(tmp)[idxs, "into"]))

    ## ## matchedFilter with wide mz bins.
    ## ## For matchedFilter I will have to do this on the profile matrix!
    ## tmp <- findChromPeaks(filterFile(od_x, 2),
    ##                       param = MatchedFilterParam(binSize = 2))
    ## idxs <- 1:nrow(chromPeaks(tmp))
    ## pkI2 <- xcms:::.getPeakInt2(tmp, chromPeaks(tmp)[idxs, , drop = FALSE])
    ## checkEquals(unname(pkI2), unname(chromPeaks(tmp)[idxs, "into"]))
}

## Test the featureValues method.
test_featureValues <- function() {
    od_x <- faahko_xod
    xs <- faahko_xs

    fdp <- PeakDensityParam(sampleGroups = xs$class)
    od_x <- groupChromPeaks(od_x, param = fdp)
    xs <- group(xs, method = "density")

    fvs <- featureValues(od_x, value = "into")
    checkEquals(rownames(fvs), rownames(featureDefinitions(od_x)))
    rownames(fvs) <- NULL
    colnames(fvs) <- NULL
    gvs <- groupval(xs, value = "into")
    rownames(gvs) <- NULL
    colnames(gvs) <- NULL
    checkEquals(fvs, gvs)
}

## Test internal helpers.
test_peakIndex <- function() {
    pkI <- xcms:::.peakIndex(xod_xg)
    checkEquals(names(pkI), rownames(featureDefinitions(xod_xg)))
    checkEquals(unname(pkI), featureDefinitions(xod_xg)$peakidx)
}

test_adjustRtimePeakGroups <- function() {
    pkGrp <- xcms:::adjustRtimePeakGroups(xod_xg,
                                          param = PeakGroupsParam(minFraction = 1))
    checkEquals(colnames(pkGrp), basename(fileNames(xod_xg)))
    ## No NAs allowed across samples:
    isNa <- apply(pkGrp, MARGIN = 1, function(z) sum(is.na(z)))
    checkTrue(all(isNa == 0))
    pkGrp <- xcms:::adjustRtimePeakGroups(
                        xod_xg, param = PeakGroupsParam(minFraction = 0.5))
    isNa <- apply(pkGrp, MARGIN = 1, function(z) sum(is.na(z)))
    checkTrue(max(isNa) == 1)

    ## Test adjustRtime adjusting also MS level > 1.
    ## Artificially changing the MS level of some spectra.
    xod_mod <- xod_xg
    ## Select the spectra for MS level 2:
    idx_ms2 <- c(300:500, 300:500 + 1277, 300:500 + 2554)
    xod_mod@featureData$msLevel[idx_ms2] <- 2
    xod_mod_adj <- adjustRtime(xod_mod,
                               param = PeakGroupsParam(span = 0.4))
    ## rtime of the MS level 2 spectra are expected to be adjusted too
    checkEquals(rtime(xod_xgr), rtime(xod_mod_adj))
    checkTrue(all(rtime(xod_mod)[idx_ms2] != rtime(xod_mod_adj)[idx_ms2]))
}

test_MS1_MS2_data <- function() {
    ## That's to test stuff for issues #208 and related (also issue #214).
    ## Set every other spectra in the original files to MS2.

    ## OnDiskMSnExp: od_x
    od_mod <- od_x
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
    checkEquals(res[msLevel(od_mod) == 1], res_2)
    ## - difference between raw and adjusted rtime at the end and beginning are
    ##   constant.
    res_by_file <- split(res, fromFile(od_mod))
    raw_by_file <- split(rtime(od_mod), fromFile(od_mod))
    checkTrue(raw_by_file[[1]][1] != res_by_file[[1]][1])
    checkTrue(raw_by_file[[2]][1] != res_by_file[[2]][1])
    checkTrue(raw_by_file[[3]][1] != res_by_file[[3]][1])
    diffs <- tail(res_by_file[[1]]) - tail(raw_by_file[[1]])
    checkEquals(unname(diff(diffs)), rep(0, 5))
    diffs <- tail(res_by_file[[2]]) - tail(raw_by_file[[2]])
    checkEquals(unname(diff(diffs)), rep(0, 5))
    diffs <- tail(res_by_file[[3]]) - tail(raw_by_file[[3]])
    checkEquals(unname(diff(diffs)), rep(0, 5))    
    ## - adjusted rtime of the MS level 2 are in interpolated between rts of
    ##   MS level 2.
    ## rtime for 3 should be interpolated between 2 and 4:
    adj_fun <- approxfun(x = raw_by_file[[1]][c(2, 4)],
                         y = res_by_file[[1]][c(2, 4)])
    checkEquals(adj_fun(raw_by_file[[1]][3]), unname(res_by_file[[1]][3]))
    adj_fun <- approxfun(x = raw_by_file[[2]][c(2, 4)],
                         y = res_by_file[[2]][c(2, 4)])
    checkEquals(adj_fun(raw_by_file[[2]][3]), unname(res_by_file[[2]][3]))
    adj_fun <- approxfun(x = raw_by_file[[3]][c(2, 4)],
                         y = res_by_file[[3]][c(2, 4)])
    checkEquals(adj_fun(raw_by_file[[3]][3]), unname(res_by_file[[3]][3]))

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
    checkEquals(rtime(res, adjusted = TRUE)[msLevel(res) == 1],
                rtime(res_2, adjusted = TRUE))
    ## - difference between raw and adjusted rtime at the end and beginning are
    ##   constant.
    res_by_file <- rtime(res, bySample = TRUE)
    raw_by_file <- rtime(xod_mod, bySample = TRUE)
    checkTrue(raw_by_file[[1]][1] != res_by_file[[1]][1])
    checkTrue(raw_by_file[[2]][1] != res_by_file[[2]][1])
    checkTrue(raw_by_file[[3]][1] != res_by_file[[3]][1])
    diffs <- tail(res_by_file[[1]]) - tail(raw_by_file[[1]])
    checkEquals(unname(diff(diffs)), rep(0, 5))
    diffs <- tail(res_by_file[[2]]) - tail(raw_by_file[[2]])
    checkEquals(unname(diff(diffs)), rep(0, 5))
    diffs <- tail(res_by_file[[3]]) - tail(raw_by_file[[3]])
    checkEquals(unname(diff(diffs)), rep(0, 5))    
    ## - adjusted rtime of the MS level 2 are in interpolated between rts of
    ##   MS level 2.
    ## rtime for 3 should be interpolated between 2 and 4:
    adj_fun <- approxfun(x = raw_by_file[[1]][c(2, 4)],
                         y = res_by_file[[1]][c(2, 4)])
    checkEquals(adj_fun(raw_by_file[[1]][3]), unname(res_by_file[[1]][3]))
    adj_fun <- approxfun(x = raw_by_file[[2]][c(2, 4)],
                         y = res_by_file[[2]][c(2, 4)])
    checkEquals(adj_fun(raw_by_file[[2]][3]), unname(res_by_file[[2]][3]))
    adj_fun <- approxfun(x = raw_by_file[[3]][c(2, 4)],
                         y = res_by_file[[3]][c(2, 4)])
    checkEquals(adj_fun(raw_by_file[[3]][3]), unname(res_by_file[[3]][3]))
}

test_extractMsData <- function() {
    ## All the data
    ## all <- extractMsData(od_x)
    ## checkEquals(length(all), length(fileNames(od_x)))
    ## rts <- split(rtime(od_x), f = fromFile(od_x))
    ## checkEquals(lengths(rts), unlist(lapply(all, nrow)))
    ## On an OnDiskMSnExp with only mz
    mzr <- c(300, 302)
    res <- extractMsData(filterFile(od_x, 1:2), mz = mzr)
    checkEquals(length(res), 2)
    checkTrue(all(res[[1]][, "mz"] >= mzr[1] & res[[1]][, "mz"] <= mzr[2]))
    checkTrue(all(res[[2]][, "mz"] >= mzr[1] & res[[2]][, "mz"] <= mzr[2]))
    ## On an OnDiskMSnExp with only rt
    rtr <- c(2500, 2800)
    res <- extractMsData(filterFile(od_x, 1:2), rt = rtr)
    checkTrue(all(res[[1]][, "rt"] >= rtr[1] & res[[1]][, "rt"] <= rtr[2]))
    checkTrue(all(res[[2]][, "rt"] >= rtr[1] & res[[2]][, "rt"] <= rtr[2]))
    ## LLLLL TODO Continue here, and then add example to the extractMsData
    ## help page.
    ## On an OnDiskMSnExp with mz and rt
    res <- extractMsData(filterFile(od_x, 1:2), rt = rtr, mz = mzr)
    checkTrue(all(res[[1]][, "rt"] >= rtr[1] & res[[1]][, "rt"] <= rtr[2]))
    checkTrue(all(res[[2]][, "rt"] >= rtr[1] & res[[2]][, "rt"] <= rtr[2]))
    checkTrue(all(res[[1]][, "mz"] >= mzr[1] & res[[1]][, "mz"] <= mzr[2]))
    checkTrue(all(res[[2]][, "mz"] >= mzr[1] & res[[2]][, "mz"] <= mzr[2]))
    
    ## XCMSnExp, xod_xgr
    ## with adjusted retention times
    tmp <- filterFile(xod_xgr, 1:2, keepAdjustedRtime = TRUE)
    checkTrue(hasAdjustedRtime(tmp))
    res <- extractMsData(tmp, rt = rtr, mz = mzr)
    mzs <- mz(tmp)
    rts <- rtime(tmp, bySample = TRUE, adjusted = TRUE)
    checkTrue(all(res[[1]][, "rt"] >= rtr[1] & res[[1]][, "rt"] <= rtr[2]))
    checkTrue(all(res[[2]][, "rt"] >= rtr[1] & res[[2]][, "rt"] <= rtr[2]))
    checkTrue(all(res[[1]][, "mz"] >= mzr[1] & res[[1]][, "mz"] <= mzr[2]))
    checkTrue(all(res[[2]][, "mz"] >= mzr[1] & res[[2]][, "mz"] <= mzr[2]))
    tmp_rts <- rts[[1]]
    tmp_rts <- tmp_rts[tmp_rts >= rtr[1] & tmp_rts <= rtr[2]]
    res_rts <- res[[1]][, 1]
    checkEquals(unique(res_rts), unname(tmp_rts))
    tmp_rts <- rts[[2]]
    tmp_rts <- tmp_rts[tmp_rts >= rtr[1] & tmp_rts <= rtr[2]]
    res_rts <- res[[2]][, 1]
    checkEquals(unique(res_rts), unname(tmp_rts))
    
    ## without adjusted retention times
    res_2 <- extractMsData(filterFile(xod_xgr, 1:2), adjustedRtime = FALSE,
                           rt = rtr, mz = mzr)
    checkTrue(all(res_2[[1]][, "rt"] >= rtr[1] & res_2[[1]][, "rt"] <= rtr[2]))
    checkTrue(all(res_2[[2]][, "rt"] >= rtr[1] & res_2[[2]][, "rt"] <= rtr[2]))
    checkTrue(all(res_2[[1]][, "mz"] >= mzr[1] & res_2[[1]][, "mz"] <= mzr[2]))
    checkTrue(all(res_2[[2]][, "mz"] >= mzr[1] & res_2[[2]][, "mz"] <= mzr[2]))
    ## checkTrue(nrow(res[[1]]) != nrow(res_2[[1]]))
    ## checkTrue(nrow(res[[2]]) != nrow(res_2[[2]]))

    ## rt and mzr out of range.
    res <- extractMsData(od_x, rt = c(6000, 6300), mz = c(0, 3))
    checkEquals(length(res), 3)
    checkTrue(all(unlist(lapply(res, FUN = nrow)) == 0))
    res <- extractMsData(od_x, rt = c(6000, 6300))
    checkEquals(length(res), 3)
    checkTrue(all(unlist(lapply(res, FUN = nrow)) == 0))
    ## res <- extractMsData(od_x, mz = c(0, 3))
    ## checkEquals(length(res), 3)
    ## checkTrue(all(unlist(lapply(res, FUN = nrow)) == 0))
}

test_spectrapply_spectra <- function() {
    ## With adjusted retention time
    tmp <- filterFile(xod_r, file = 3, keepAdjustedRtime = TRUE)
    checkTrue(hasAdjustedRtime(tmp))
    checkTrue(is.character(all.equal(rtime(tmp, adjusted = FALSE),
                                     rtime(tmp, adjusted = TRUE))))
    sps <- spectra(tmp)
    checkEquals(unlist(lapply(sps, rtime)), rtime(tmp, adjusted = TRUE))
    checkEquals(unlist(spectrapply(tmp, FUN = function(x) rtime(x))),
                rtime(tmp, adjusted = TRUE))
    sps_2 <- spectra(tmp, adjusted = FALSE)
    checkEquals(unlist(lapply(sps_2, rtime)), rtime(tmp, adjusted = FALSE))
    ## without adjusted retention time
    tmp <- filterFile(xod_x, file = 3)
    checkTrue(!hasAdjustedRtime(tmp))
    sps <- spectra(tmp)
    checkEquals(unlist(lapply(sps, rtime)), rtime(tmp))
    checkEquals(unlist(lapply(sps, mz)), unlist(mz(tmp)))
    checkEquals(unlist(spectrapply(tmp, FUN = function(x) rtime(x))),
                rtime(tmp))
}

test_processHistory <- function() {
    type_peak_det <- xcms:::.PROCSTEP.PEAK.DETECTION
    type_align <- xcms:::.PROCSTEP.RTIME.CORRECTION
    type_corr <- xcms:::.PROCSTEP.PEAK.GROUPING
    checkTrue(length(processHistory(xod_x, type = type_corr)) == 0)
    ph <- processHistory(xod_x, type = type_peak_det)
    checkEquals(as.character(class(processParam(ph[[1]]))), "CentWaveParam")

    ph <- processHistory(xod_xgrg)
    checkTrue(length(ph) == 4)
    ph <- processHistory(xod_xgrg, msLevel = 1L)
    checkTrue(length(ph) == 2)
    checkEquals(as.character(class(processParam(ph[[1]]))), "CentWaveParam")
}

test_split <- function() {
    xod <- as(od_x, "XCMSnExp")
    tmp <- split(xod_xgr, f = fromFile(xod_xgr))
    ## Split by file.
    checkEquals(spectra(tmp[[1]]), spectra(filterFile(xod, file = 1)))
    checkEquals(spectra(tmp[[3]]), spectra(filterFile(xod, file = 3)))
    ## Split by acquisitionNum.
    tmp <- filterRt(xod_xgr, rt = c(2500, 2700))
    checkTrue(hasChromPeaks(tmp))
    checkTrue(hasAdjustedRtime(tmp))
    tmp_2 <- split(tmp, f = acquisitionNum(tmp))
    checkTrue(all(acquisitionNum(tmp_2[[1]]) == acquisitionNum(tmp)[1]))
    checkTrue(all(acquisitionNum(tmp_2[[14]]) == acquisitionNum(tmp)[14]))
    ## with keepAdjustedRtime
    tmp <- split(xod_xgr, f = fromFile(xod_xgr), keepAdjustedRtime = TRUE)
    tmp_1 <- filterFile(xod_xgr, file = 1, keepAdjustedRtime = TRUE)
    checkTrue(hasAdjustedRtime(tmp_1))
    checkEquals(rtime(tmp[[1]]), rtime(tmp_1))
    tmp_2 <- filterFile(xod_xgr, file = 2, keepAdjustedRtime = TRUE)
    checkTrue(hasAdjustedRtime(tmp_2))
    checkEquals(rtime(tmp[[2]]), rtime(tmp_2))
    tmp_3 <- filterFile(xod_xgr, file = 3, keepAdjustedRtime = TRUE)
    checkTrue(hasAdjustedRtime(tmp_3))
    checkEquals(rtime(tmp[[3]]), rtime(tmp_3))
    checkTrue(!all(rtime(tmp[[3]]) == rtime(tmp[[3]], adjusted = FALSE)))
}


############################################################
## Test getEIC alternatives.
dontrun_getEIC_alternatives <- function() {
    library(RUnit)
    library(xcms)

    fls <- c(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
             system.file('cdf/KO/ko16.CDF', package = "faahKO"),
             system.file('cdf/KO/ko18.CDF', package = "faahKO"))
    od <- readMSData(fls, mode = "inDisk")
    cwp <- CentWaveParam(noise = 10000, snthresh = 40)
    od_x <- findChromPeaks(od, param = cwp)

    ## ## with this one we get 3 spectras back, one in each file.
    ## rtr <- c(2787, 2788)    
    ## res <- filterRt(od_x, rt = rtr)

    ## ## -----------
    ## ## That's to test .extractChromatogram
    ## mzr <- c(279, 279)
    ## chrs <- extractChromatograms(od_x, mzrange = mzr)
    ## ##   input parameter
    ## x <- od_x
    ## rm(rtrange)
    ## rm(mzrange)
    ## mzrange <- mzr
    ## aggregationFun <- "sum"
    ## ##   function call
    ## ## -----------
    
    ## od_xg <- groupChromPeaks(od_x, param = PeakDensityParam())
    ## od_xgr <- adjustRtime(od_xg, param = PeakGroupsParam(span = 0.4))

    ## rtr <- as.matrix(featureDefinitions(od_xg)[1:5, c("rtmin", "rtmax")])
    ## mzr <- as.matrix(featureDefinitions(od_xg)[1:5, c("mzmin", "mzmax")])

    ## system.time(
    ##     res1 <- xcms:::.extractMsData(od, rtrange = rtr[1, ], mzrange = mzr[1, ])
    ## )
    ## system.time(
    ##     res2 <- xcms:::.extractMsData(od_xgr, rtrange = rtr[1, ], mzrange = mzr[1, ])
    ## )
    ## system.time(
    ##     res1 <- xcms:::.sliceApply(od, rtrange = rtr[1, ], mzrange = mzr[1, ])
    ## )
    ## system.time(
    ##     res1 <- xcms:::.sliceApply(od_xgr, rtrange = rtr[1, ], mzrange = mzr[1, ])
    ## )
    
    ## library(profvis)
    ## profvis(res <- xcms:::.extractMsData(od, rtrange = rtr[1, ], mzrange = mzr[1, ]))
    

    ## Compare with getEIC
    xs <- as(od_x, "xcmsSet")
    sampclass(xs) <- rep("KO", 3)
    xs_2 <- group(xs)
    suppressWarnings(
        xs_2 <- retcor(xs_2)
    )
    xs_2 <- group(xs_2)

    rtr <- groups(xs_2)[1:5, c("rtmin", "rtmax")]
    mzr <- groups(xs_2)[1:5, c("mzmin", "mzmax")]

    ##

    register(SerialParam())
    od <- as(od_x, "OnDiskMSnExp")
    ## Get all of em.
    chrs <- xcms:::.extractMultipleChromatograms(od, rt = rtr, mz = mzr)
    for (i in 1:nrow(rtr)) {
        chrs1 <- extractChromatograms(od_x, rt = rtr[i, ], mz = mzr[i, ])
        checkEquals(unname(chrs1), unname(chrs[[i]]))
    }

    library(microbenchmark)
    microbenchmark(xcms:::.extractChromatogram(od, rt = rtr[1, ], mz = mzr[1, ]),
                   xcms:::.extractMultipleChromatograms(od, rt = rtr[1, , drop = FALSE],
                                                        mz = mzr[1, , drop = FALSE]),
                   times = 10)

    library(profvis)
    profvis(xcms:::.extractMultipleChromatograms(od, rt = rtr[1, , drop = FALSE],
                                                 mz = mzr[1, , drop = FALSE]))
    
    ## Extract the EIC:
    system.time(
        eic <- getEIC(xs_2, rtrange = rtr, mzrange = mzr, rt = "raw")
    ) ## 5.5 sec
    system.time(
        eic2 <- xcms:::.extractMultipleChromatograms(od, rt = rtr, mz = mzr)
    ) ## 0.13 sec

    

    
    ## Now try to do the same using MSnbase stuff.
    system.time(
        res <- xcms:::.extractMsData(od, rtrange = rtr[1, ], mzrange = mzr[1, ])
    ) ## 0.7 sec.
    system.time(
        res <- xcms:::.extractMsData(od_x, rtrange = rtr[1, ], mzrange = mzr[1, ])
    ) ## 0.74 sec.
    
    
    rts <- rtime(od_x)
    idx <- apply(rtr, MARGIN = 1, function(z) {
        which(rts >= z[1] & rts <= z[2])
    })
    ## Subset the od_x
    od_ss <- od[sort(unique(unlist(idx)))]
    system.time(
        scts <- spectra(od_ss)
    ) ## 0.7 secs.

    ## Tests.
    rtrange <- rtr[3, ]
    mzrange <- mzr[3, ]
    system.time(
        tmp <- extractMsData(od, rtrange = rtr[1, ], mzrange = mzr[1, ])
    )

    ## For all of em:
    system.time(
        tmp <- mapply(rtrange = split(rtr, f = seq_len(nrow(rtr))),
                      mzrange = split(mzr, f = seq_len(nrow(mzr))),
                      FUN = extractMsData, MoreArgs = list(x = od),
                      SIMPLIFY = FALSE, USE.NAMES = FALSE)
    ) ## 3.589
    ## That's fine as long as we don't extract too many of em.

    ############################################################
    ## Alternative: do it by file.
    ## IF it's an XCMSnExp: coerce to OnDiskMSnExp by replacing the rtime with
    ## the adjusted rtime.
    ## 1) Subset the od selecting all spectra that fall into the rt ranges.
    ##    keep_logical <- have_rt >= rt[1] & have_rt <= rt[2]
    ##    tmp <- as(object, "OnDiskMSnExp")[base::which(keep_logical)]
    ## 2) Call a spectrapply, passing the matrix of rts and mzs.
    ##    (Load the spectra (without any filtering now).)
    ## 3) spectrapply function loops over the rtrange and mzrange:
    ##    - select all spectra that are within the range.
    ##    - lapply on those, apply filterMz with the current mz range.
    ##    - return a list of Chromatogram classes.
    ## 4) We get a list of list of Chromatogram objects. [[files]][[ranges]].
    ##    Rearrange the lists: [[ranges]][[files]].
    ## Could also put that into a DataFrame... [ranges, files]
    


    ## For a single one:
    ## BE AWARE: we have use the fileNames matching to the original file names to
    ## match to the CORRECT indices.
    system.time(
        dfs <- spectrapply(tmp, as.data.frame)
    )

    
    ## mz outside:
    mzrange <- c(600, 601)
    tmp <- filterMz(filterRt(od, rt = rtrange), mz = mzrange)
    dfs <- spectrapply(tmp, as.data.frame)
    ## Funny, it returns something - 0.
}


