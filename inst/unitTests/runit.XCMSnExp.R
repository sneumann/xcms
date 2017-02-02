## tests related to the new XCMSnExp object.
library(RUnit)

od_x <- faahko_xod
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
    checkTrue(!hasAlignedFeatures(x))
    checkTrue(!hasDetectedFeatures(x))
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
    checkException(xod@msFeatureData$features <- 3)
    checkException(xod@msFeatureData$bla <- 4)
    checkTrue(validObject(xod))

    .checkCreationOfEmptyObject()
    ## xod@msFeatureData$features <- xs_2@peaks
    ## checkTrue(validObject(xod))
    ## xod@msFeatureData$features[1, "sample"] <- 40
    ## checkException(validObject(xod))
    ## xod@msFeatureData$features[1, "sample"] <- 3
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
    rts_3 <- rtime(od_x, bySample = TRUE)
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
    mzs_2 <- mz(od_x, bySample = TRUE)
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
    ints_2 <- intensity(od_x, bySample = TRUE)
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
}

test_XCMSnExp_class_accessors <- function() {
    .checkCreationOfEmptyObject()
    ## Filling with data...
    xod <- as(od_fa, "XCMSnExp")
    ## features
    checkTrue(!hasDetectedFeatures(xod))
    features(xod) <- xs_2@peaks
    checkTrue(hasDetectedFeatures(xod))
    checkEquals(features(xod), xs_2@peaks)
    checkException(features(xod) <- 4)
    tmp <- features(xod, bySample = TRUE)
    checkTrue(length(tmp) == length(fileNames(xod)))
    tmp <- do.call(rbind, tmp)
    rownames(tmp) <- NULL
    checkEquals(tmp, features(xod))
    ## Wrong assignments.
    pks <- xs_2@peaks
    pks[1, "sample"] <- 40
    checkException(features(xod) <- pks)
    ## featureGroups
    checkTrue(!hasAlignedFeatures(xod))
    library(S4Vectors)
    fd <- DataFrame(xs_2@groups)
    fd$featureidx <- xs_2@groupidx
    featureGroups(xod) <- fd
    checkTrue(hasDetectedFeatures(xod))
    checkTrue(hasAlignedFeatures(xod))
    checkEquals(featureGroups(xod), fd)
    ## adjustedRtime
    checkTrue(!hasAdjustedRtime(xod))
    adjustedRtime(xod) <- xs_2@rt$corrected
    checkTrue(hasAdjustedRtime(xod))
    checkTrue(hasDetectedFeatures(xod))
    checkTrue(hasAlignedFeatures(xod))
    checkEquals(adjustedRtime(xod, bySample = TRUE), xs_2@rt$corrected)
    ## rtime should now also return adjusted retention times
    checkEquals(rtime(xod), adjustedRtime(xod))
    checkEquals(rtime(xod, adjusted = FALSE), rtime(as(xod, "OnDiskMSnExp")))
    checkEquals(rtime(xod, adjusted = TRUE), adjustedRtime(xod))
    ## Indirect test that the ordering of the adjusted retention times matches
    ## ordering of rtime.
    tmp <- unlist(adjustedRtime(xod, bySample = TRUE))
    tmp_diff <- tmp - rtime(xod)
    tmp_diff_2 <- adjustedRtime(xod, bySample = FALSE) - rtime(xod)
    checkTrue(max(tmp_diff) > max(tmp_diff_2))
    checkEquals(names(adjustedRtime(xod)), names(rtime(xod)))
    ## Wrong assignments.
    checkException(adjustedRtime(xod) <- xs_2@rt$corrected[1:2])
    .checkCreationOfEmptyObject()
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
    type_feat_det <- xcms:::.PROCSTEP.FEATURE.DETECTION
    type_feat_algn <- xcms:::.PROCSTEP.FEATURE.ALIGNMENT
    type_rt_adj <- xcms:::.PROCSTEP.RTIME.CORRECTION
    ## Perform alignment.
    od_xg <- groupFeatures(od_x, param = FeatureDensityParam())
    checkTrue(hasAlignedFeatures(od_xg))
    checkTrue(length(processHistory(od_xg, type = type_feat_algn)) == 1)
    ## Retention time adjustment.
    od_xgr <- adjustRtime(od_xg, param = FeatureGroupsParam(span = 1))
    checkTrue(hasDetectedFeatures(od_xgr))
    checkTrue(length(processHistory(od_xgr, type = type_feat_det)) == 1)
    checkTrue(!hasAlignedFeatures(od_xgr))  ## These should have been removed
    checkTrue(length(processHistory(od_xgr, type = type_feat_algn)) == 1)
    checkTrue(hasAdjustedRtime(od_xgr))
    checkTrue(length(processHistory(od_xgr, type = type_rt_adj)) == 1)
    ## Most of the retention times are different
    checkTrue(sum(features(od_xgr)[, "rt"] != features(od_x)[, "rt"]) >
              nrow(features(od_x)) / 2)
    ## Alignment after retention time adjustment.
    od_xgrg <- groupFeatures(od_xgr, param = FeatureDensityParam())
    checkTrue(hasDetectedFeatures(od_xgrg))
    checkEquals(features(od_xgrg), features(od_xgr))
    checkTrue(hasAdjustedRtime(od_xgrg))
    checkTrue(length(processHistory(od_xgr, type = type_feat_algn)) == 1)
    checkTrue(hasAlignedFeatures(od_xgrg))
    checkTrue(length(processHistory(od_xgrg, type = type_feat_algn)) == 2)
       
    ## 1) dropDetectedFeatures: delete all process history steps and all data.
    res <- dropFeatures(od_x)
    checkTrue(!hasDetectedFeatures(res))
    checkTrue(length(processHistory(res, type = type_feat_det)) == 0)
    checkTrue(!hasAlignedFeatures(res))
    checkTrue(length(processHistory(res, type = type_feat_algn)) == 0)
    checkTrue(!hasAdjustedRtime(res))
    checkTrue(length(processHistory(res, type = type_rt_adj)) == 0)
    ##
    res <- dropFeatures(od_xg)
    checkTrue(!hasDetectedFeatures(res))
    checkTrue(length(processHistory(res, type = type_feat_det)) == 0)
    checkTrue(!hasAlignedFeatures(res))
    checkTrue(length(processHistory(res, type = type_feat_algn)) == 0)
    checkTrue(!hasAdjustedRtime(res))
    checkTrue(length(processHistory(res, type = type_rt_adj)) == 0)
    ##
    res <- dropFeatures(od_xgr)
    checkTrue(!hasDetectedFeatures(res))
    checkTrue(length(processHistory(res, type = type_feat_det)) == 0)
    checkTrue(!hasAlignedFeatures(res))
    checkTrue(length(processHistory(res, type = type_feat_algn)) == 0)
    checkTrue(!hasAdjustedRtime(res))
    checkTrue(length(processHistory(res, type = type_rt_adj)) == 0)
    ##
    res <- dropFeatures(od_xgrg)
    checkTrue(!hasDetectedFeatures(res))
    checkTrue(length(processHistory(res, type = type_feat_det)) == 0)
    checkTrue(!hasAlignedFeatures(res))
    checkTrue(length(processHistory(res, type = type_feat_algn)) == 0)
    checkTrue(!hasAdjustedRtime(res))
    checkTrue(length(processHistory(res, type = type_rt_adj)) == 0)
    
    ## 2) dropFeatureGroups:
    ##    a) drop the feature groups and the latest related process history
    ##    b) if retention time correction was performed AFTER the latest feature
    ##       grouping, drop also the retention time correction and all related
    ##       process histories.
    res <- dropFeatureGroups(od_xg)
    checkEquals(res, od_x)
    checkTrue(hasDetectedFeatures(res))
    checkTrue(length(processHistory(res, type = type_feat_det)) == 1)
    checkTrue(!hasAlignedFeatures(res))
    checkTrue(length(processHistory(res, type = type_feat_algn)) == 0)
    checkTrue(!hasAdjustedRtime(res))
    checkTrue(length(processHistory(res, type = type_rt_adj)) == 0)
    ## No feature groups - so there is nothing that this function does here.
    res <- dropFeatureGroups(od_xgr)
    checkEquals(res, od_xgr)
    checkTrue(hasDetectedFeatures(res))
    checkTrue(length(processHistory(res, type = type_feat_det)) == 1)
    checkTrue(!hasAlignedFeatures(res))
    checkTrue(length(processHistory(res, type = type_feat_algn)) == 1)
    checkTrue(hasAdjustedRtime(res))
    checkTrue(length(processHistory(res, type = type_rt_adj)) == 1)
    ## Remove the latest ones.
    res <- dropFeatureGroups(od_xgrg)
    checkEquals(res, od_xgr)
    checkTrue(hasDetectedFeatures(res))
    checkTrue(length(processHistory(res, type = type_feat_det)) == 1)
    checkTrue(!hasAlignedFeatures(res))
    checkTrue(length(processHistory(res, type = type_feat_algn)) == 1)
    checkTrue(hasAdjustedRtime(res))
    checkTrue(length(processHistory(res, type = type_rt_adj)) == 1)

    ## 3) dropAdjustedRtime:
    ##    a) drop the retention time adjustment and related process histories
    ##    b) if grouping has been performed AFTER retention time correction,
    ##       drop the feature alignment and all related process histories.
    ##    c) if grouping has been performed BEFORE retention time correction,
    ##       do nothing.
    res <- dropAdjustedRtime(od_xg)
    checkEquals(res, od_xg)
    ## This drops also the process history for alignment.
    res <- dropAdjustedRtime(od_xgr)
    checkTrue(hasDetectedFeatures(res))
    checkTrue(length(processHistory(res, type = type_feat_det)) == 1)
    checkTrue(!hasAlignedFeatures(res))
    checkTrue(length(processHistory(res, type = type_feat_algn)) == 0)
    checkTrue(!hasAdjustedRtime(res))
    checkTrue(length(processHistory(res, type = type_rt_adj)) == 0)
    checkEquals(features(res), features(od_x))
    checkEquals(res, od_x)
    checkEquals(rtime(res), rtime(od_x))
    checkEquals(rtime(res), rtime(od_xgr, adjusted = FALSE))
    ## This drops also the feature alignment performed later.
    res <- dropAdjustedRtime(od_xgrg)
    checkTrue(hasDetectedFeatures(res))
    checkTrue(length(processHistory(res, type = type_feat_det)) == 1)
    checkTrue(!hasAlignedFeatures(res))
    checkTrue(length(processHistory(res, type = type_feat_algn)) == 0)
    checkTrue(!hasAdjustedRtime(res))
    checkTrue(length(processHistory(res, type = type_rt_adj)) == 0)
    checkEquals(features(res), features(od_x))
    checkEquals(res, od_x)
    checkEquals(rtime(res), rtime(od_xgrg, adjusted = FALSE))
    
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
        tmp_2 <- od_x[1:10]
    )
    checkTrue(length(processHistory(tmp_2)) == 0)
    checkTrue(!hasDetectedFeatures(tmp_2))
    tmp_1@processingData <- new("MSnProcess")
    tmp_2@processingData <- new("MSnProcess")
    checkEquals(tmp_1, as(tmp_2, "OnDiskMSnExp"))
    ## bin
    tmp_1 <- bin(od_fa)
    suppressWarnings(
        tmp_2 <- bin(od_x)
    )
    checkTrue(length(processHistory(tmp_2)) == 0)
    checkTrue(!hasDetectedFeatures(tmp_2))
    tmp_1@processingData <- new("MSnProcess")
    tmp_2@processingData <- new("MSnProcess")
    checkEquals(tmp_1, as(tmp_2, "OnDiskMSnExp"))
    ## clean
    tmp_1 <- clean(od_fa)
    suppressWarnings(
        tmp_2 <- clean(od_x)
    )
    checkTrue(length(processHistory(tmp_2)) == 0)
    checkTrue(!hasDetectedFeatures(tmp_2))
    tmp_1@processingData <- new("MSnProcess")
    tmp_2@processingData <- new("MSnProcess")
    checkEquals(tmp_1, as(tmp_2, "OnDiskMSnExp"))
    ## filterAcquisitionNum
    tmp_1 <- filterAcquisitionNum(od_fa)
    suppressWarnings(
        tmp_2 <- filterAcquisitionNum(od_x)
    )
    checkTrue(length(processHistory(tmp_2)) == 0)
    checkTrue(!hasDetectedFeatures(tmp_2))
    tmp_1@processingData <- new("MSnProcess")
    tmp_2@processingData <- new("MSnProcess")
    checkEquals(tmp_1, as(tmp_2, "OnDiskMSnExp"))
    ## filterMsLevel
    tmp_1 <- filterMsLevel(od_fa)
    suppressWarnings(
        tmp_2 <- filterMsLevel(od_x)
    )
    checkTrue(length(processHistory(tmp_2)) == 0)
    checkTrue(!hasDetectedFeatures(tmp_2))
    tmp_1@processingData <- new("MSnProcess")
    tmp_2@processingData <- new("MSnProcess")
    checkEquals(tmp_1, as(tmp_2, "OnDiskMSnExp"))
    ## normalize
    tmp_1 <- normalize(od_fa)
    suppressWarnings(
        tmp_2 <- normalize(od_x)
    )
    checkTrue(length(processHistory(tmp_2)) == 0)
    checkTrue(!hasDetectedFeatures(tmp_2))
    tmp_1@processingData <- new("MSnProcess")
    tmp_2@processingData <- new("MSnProcess")
    checkEquals(tmp_1, as(tmp_2, "OnDiskMSnExp"))
    ## pickPeaks
    tmp_1 <- pickPeaks(od_fa)
    suppressWarnings(
        tmp_2 <- pickPeaks(od_x)
    )
    checkTrue(length(processHistory(tmp_2)) == 0)
    checkTrue(!hasDetectedFeatures(tmp_2))
    tmp_1@processingData <- new("MSnProcess")
    tmp_2@processingData <- new("MSnProcess")
    checkEquals(tmp_1, as(tmp_2, "OnDiskMSnExp"))
    ## removePeaks
    tmp_1 <- removePeaks(od_fa)
    suppressWarnings(
        tmp_2 <- removePeaks(od_x)
    )
    checkTrue(length(processHistory(tmp_2)) == 0)
    checkTrue(!hasDetectedFeatures(tmp_2))
    tmp_1@processingData <- new("MSnProcess")
    tmp_2@processingData <- new("MSnProcess")
    checkEquals(tmp_1, as(tmp_2, "OnDiskMSnExp"))
    ## smooth
    tmp_1 <- smooth(od_fa)
    suppressWarnings(
        tmp_2 <- smooth(od_x)
    )
    checkTrue(length(processHistory(tmp_2)) == 0)
    checkTrue(!hasDetectedFeatures(tmp_2))
    tmp_1@processingData <- new("MSnProcess")
    tmp_2@processingData <- new("MSnProcess")
    checkEquals(tmp_1, as(tmp_2, "OnDiskMSnExp"))
    .checkCreationOfEmptyObject()
}

## Test XCMSnExp filter methods.
test_XCMSnExp_filterFile <- function() {
    ## filterFile
    tmp <- filterFile(od_x, file = 2)
    checkException(tmp@msFeatureData$bla <- 3)
    checkTrue(!hasAdjustedRtime(tmp))
    checkTrue(!hasAlignedFeatures(tmp))
    checkTrue(all(features(tmp)[, "sample"] == 1))
    checkEquals(features(tmp)[, -ncol(features(tmp))],
                features(od_x)[features(od_x)[, "sample"] == 2,
                               -ncol(features(od_x))])
    checkEquals(fileIndex(processHistory(tmp)[[1]]), 1)
    ## check with other index.
    tmp <- filterFile(od_x, file = c(1, 3))
    checkTrue(!hasAdjustedRtime(tmp))
    checkTrue(!hasAlignedFeatures(tmp))
    checkTrue(all(features(tmp)[, "sample"] %in% c(1, 2)))
    a <- features(tmp)
    b <- features(od_x)
    checkEquals(a[, -ncol(a)], b[b[, "sample"] %in% c(1, 3), -ncol(b)])
    checkEquals(fileIndex(processHistory(tmp)[[1]]), c(1, 2))

    ## Errors
    checkException(filterFile(od_x, file = 5))
    checkException(filterFile(od_x, file = 1:5))

    ## Little mockup to check correctness of Process history.
    od_2 <- od_x
    od_2 <- xcms:::addProcessHistory(od_2,
                                     xcms:::ProcessHistory(type = xcms:::.PROCSTEP.RTIME.CORRECTION))
    od_2 <- xcms:::addProcessHistory(od_2, xcms:::ProcessHistory(type = xcms:::.PROCSTEP.UNKNOWN, fileIndex = 2, info. = "I should be here"))
    od_2 <- xcms:::addProcessHistory(od_2, xcms:::ProcessHistory(type = xcms:::.PROCSTEP.UNKNOWN, fileIndex = 1, info. = "EEEEEE"))

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
}

test_XCMSnExp_filterMz <- function() {
    od_x2 <- od_x
    new_e <- xcms:::.copy_env(od_x2@msFeatureData)
    xcms:::adjustedRtime(od_x2) <- xs_2@rt$corrected
    library(S4Vectors)
    fd <- DataFrame(xs_2@groups)
    fd$featureidx <- xs_2@groupidx
    featureGroups(od_x2) <- fd

    ## Subset
    tmp <- filterMz(od_x2, mz = c(300, 400))
    checkException(tmp@msFeatureData$bla <- 3)
    checkTrue(length(tmp@spectraProcessingQueue) == 1)
    checkTrue(all(features(tmp)[, "mz"] >= 300 & features(tmp)[, "mz"] <= 400))
    checkTrue(validObject(tmp@msFeatureData))
    checkTrue(all(featureGroups(tmp)$mzmed >= 300 &
                                    featureGroups(tmp)$mzmed <= 400))
    checkEquals(adjustedRtime(tmp), adjustedRtime(od_x2))
    checkTrue(nrow(features(tmp)) < nrow(features(od_x2)))
    checkTrue(hasDetectedFeatures(tmp))
    checkTrue(hasAlignedFeatures(tmp))
    ## second
    tmp <- filterMz(od_x, mz = c(300, 400))
    checkTrue(hasDetectedFeatures(tmp))
    checkTrue(!hasAlignedFeatures(tmp))
    checkTrue(all(features(tmp)[, "mz"] >= 300 & features(tmp)[, "mz"] <= 400))
    checkTrue(validObject(tmp@msFeatureData))
}

test_XCMSnExp_filterRt <- function() {
    od_x2 <- od_x

    ## Testing with only feature data present:
    res <- filterRt(od_x2, rt = c(2700, 2900))
    ## MsFeatureData has to be locked!
    checkException(res@msFeatureData$bla <- 3)
    ## Retention time has to be within the range.
    checkTrue(all(rtime(res) >= 2700 & rtime(res) <= 2900))
    ## features have to be within the range.
    checkTrue(all(features(res)[, "rt"] >= 2700 &
                  features(res)[, "rt"] <= 2900))
    ## features have to match the subsetted ones.
    are_within <- features(od_x2)[, "rt"] >= 2700 &
        features(od_x2)[, "rt"] <= 2900
    checkEquals(features(res), features(od_x2)[are_within,])
    ## Have a feature detection process history.
    checkEquals(processType(processHistory(res)[[1]]),
                xcms:::.PROCSTEP.FEATURE.DETECTION)
    ## filter such that we keep some spectra but no features:
    res <- filterRt(od_x2, rt = c(4200, 4400))
    checkTrue(all(rtime(res) >= 4200 & rtime(res) <= 4400))
    checkTrue(!hasDetectedFeatures(res))
    checkTrue(length(processHistory(res)) == 0)
    ## No rt
    res <- filterRt(od_x2, rt = c(10, 20))
    checkTrue(length(res) == 0)

    ## With feature groups:
    od_xg <- groupFeatures(od_x, param = FeatureDensityParam())
    res <- filterRt(od_xg, rt = c(2700, 2900))
    checkEquals(hasDetectedFeatures(res), hasDetectedFeatures(od_xg))
    checkEquals(hasAlignedFeatures(res), hasAlignedFeatures(od_xg))
    checkTrue(length(processHistory(res, type = "Feature detection")) == 1)
    checkTrue(length(processHistory(res, type = "Feature alignment")) == 1)
    ## Retention time has to be within the range.
    checkTrue(all(rtime(res) >= 2700 & rtime(res) <= 2900))
    ## features have to be within the range.
    checkTrue(all(features(res)[, "rt"] >= 2700 &
                  features(res)[, "rt"] <= 2900))
    ## feature groups have to be within range:
    checkTrue(all(featureGroups(res)$rtmed >= 2700 &
                                    featureGroups(res)$rtmed <= 2900))
    ## All feature idx have to match.
    checkTrue(all(unlist(featureGroups(res)$featureidx) %in%
                  1:nrow(features(res))))
    ## Filter such that we don't have any features.
    res <- filterRt(od_xg, rt = c(4200, 4400))
    checkTrue(all(rtime(res) >= 4200 & rtime(res) <= 4400))
    checkTrue(!hasDetectedFeatures(res))
    checkTrue(!hasAlignedFeatures(res))
    checkTrue(length(processHistory(res)) == 0)
    ## No rt
    res <- filterRt(od_xg, rt = c(10, 20))
    checkTrue(length(res) == 0)

    ## With adjusted retention time.
    od_xgr <- adjustRtime(od_xg, param = FeatureGroupsParam(span = 0.5))
    res <- filterRt(od_xgr, rt = c(2700, 2900))
    checkEquals(hasDetectedFeatures(res), hasDetectedFeatures(od_xgr))
    checkEquals(hasAlignedFeatures(res), hasAlignedFeatures(od_xgr))
    checkEquals(hasAdjustedRtime(res), hasAdjustedRtime(od_xgr))
    checkTrue(length(processHistory(res, type = "Feature detection")) == 1)
    checkTrue(length(processHistory(res, type = "Feature alignment")) == 1)
    checkTrue(length(processHistory(res, type = "Retention time correction")) == 1)
    ## Retention time has to be within the range. By default we're supposed to
    ## filter by adjusted rtime.
    checkTrue(all(rtime(res) >= 2700 & rtime(res) <= 2900))
    checkTrue(!all(rtime(res, adjusted = FALSE) >= 2700 &
                   rtime(res, adjusted = FALSE) <= 2900))
    ## Check if rtime is what we expect it to be:
    keep_em <- rtime(od_xgr) >= 2700 & rtime(od_xgr) <= 2900
    checkEquals(rtime(res), rtime(od_xgr)[keep_em])
    checkEquals(rtime(res, adjusted = FALSE),
                rtime(od_xgr, adjusted = FALSE)[keep_em])
    checkEquals(adjustedRtime(res), adjustedRtime(od_xgr)[keep_em])
    ## Check features.
    keep_em <- features(od_xgr)[, "rt"] >= 2700 & features(od_xgr)[, "rt"] <= 2900
    checkEquals(features(od_xgr)[keep_em, ], features(res))
    ## Filtering for rt with no features drops also the adjusted rt.
    res <- filterRt(od_xgr, rt = c(4200, 4400))
    checkTrue(!hasDetectedFeatures(res))
    checkTrue(!hasAdjustedRtime(res))
    checkTrue(length(processHistory(res)) == 0)
    ## Enforce filtering on adjusted:
    checkEquals(filterRt(od_xgr, rt = c(2700, 2900)),
                filterRt(od_xgr, rt = c(2700, 2900), adjusted = TRUE))
    ## Filter using raw retention times
    res <- filterRt(od_xgr, rt = c(2700, 2900), adjusted = FALSE)
    checkTrue(!all(rtime(res) >= 2700 & rtime(res) <= 2900))
    checkTrue(all(rtime(res, adjusted = FALSE) >= 2700 &
                  rtime(res, adjusted = FALSE) <= 2900))
    ## Features - can not really check here.
    checkTrue(all(features(res)[, "rt"] >= 2700 &
                  features(res)[, "rt"] <= 2900))
    
    ## Adjusted retention time AND grouping
    od_xgrg <- groupFeatures(od_xgr, param = FeatureDensityParam())
    res <- filterRt(od_xgrg, rt = c(2700, 2900))
    checkEquals(hasDetectedFeatures(res), hasDetectedFeatures(od_xgrg))
    checkEquals(hasAlignedFeatures(res), hasAlignedFeatures(od_xgrg))
    checkEquals(hasAdjustedRtime(res), hasAdjustedRtime(od_xgrg))
    checkTrue(length(processHistory(res, type = "Feature detection")) == 1)
    checkTrue(length(processHistory(res, type = "Feature alignment")) == 2)
    checkTrue(length(processHistory(res, type = "Retention time correction")) == 1)
    ## Retention time has to be within the range. By default we're supposed to
    ## filter by adjusted rtime.
    checkTrue(all(rtime(res) >= 2700 & rtime(res) <= 2900))
    checkTrue(!all(rtime(res, adjusted = FALSE) >= 2700 &
                   rtime(res, adjusted = FALSE) <= 2900))
    ## Check if rtime is what we expect it to be:
    keep_em <- rtime(od_xgrg) >= 2700 & rtime(od_xgrg) <= 2900
    checkEquals(rtime(res), rtime(od_xgrg)[keep_em])
    checkEquals(rtime(res, adjusted = FALSE),
                rtime(od_xgrg, adjusted = FALSE)[keep_em])
    checkEquals(adjustedRtime(res), adjustedRtime(od_xgrg)[keep_em])
    ## Check features.
    keep_em <- features(od_xgrg)[, "rt"] >= 2700 & features(od_xgrg)[, "rt"] <= 2900
    checkEquals(features(od_xgrg)[keep_em, ], features(res))
    ## Feature groups.
    checkTrue(all(featureGroups(res)$rtmed >= 2700 &
                                    featureGroups(res)$rtmed <= 2900))
    validObject(res)
    ## Filtering for rt with no features drops also the adjusted rt.
    res <- filterRt(od_xgrg, rt = c(4200, 4400))
    checkTrue(!hasDetectedFeatures(res))
    checkTrue(!hasAdjustedRtime(res))
    checkTrue(!hasAlignedFeatures(res))
    checkTrue(length(processHistory(res)) == 0)
    ## Enforce filtering on adjusted:
    checkEquals(filterRt(od_xgrg, rt = c(2700, 2900)),
                filterRt(od_xgrg, rt = c(2700, 2900), adjusted = TRUE))
    ## Filter using raw retention times
    res <- filterRt(od_xgrg, rt = c(2700, 2900), adjusted = FALSE)
    checkTrue(!all(rtime(res) >= 2700 & rtime(res) <= 2900))
    checkTrue(all(rtime(res, adjusted = FALSE) >= 2700 &
                  rtime(res, adjusted = FALSE) <= 2900))
    ## Features - can not really check here.
    checkTrue(all(features(res)[, "rt"] >= 2700 &
                  features(res)[, "rt"] <= 2900))

}

## Test the coercion method.
test_as_XCMSnExp_xcmsSet <- function() {
    od_x <- faahko_xod
    res <- xcms:::.XCMSnExp2xcmsSet(od_x)
    res <- as(od_x, "xcmsSet")
    ## Results should be the same as in xs.
    checkEquals(res@peaks, features(od_x))
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
    od_2 <- groupFeatures(od_x, param = FeatureDensityParam())
    checkEquals(featureGroups(od_2)$featureidx, groupidx(res))

    ## rt correction
    od_3 <- adjustRtime(od_2, param = FeatureGroupsParam(minFraction = 1,
                                                         span = 0.4))
    ## With groups.
    res <- as(od_2, "xcmsSet")
    checkEquals(res@groups,
                as.matrix(featureGroups(od_2)[, -ncol(featureGroups(od_2))]))
    checkEquals(res@groupidx, featureGroups(od_2)$featureidx)

    ## With adjusted retention time.
    res_2 <- retcor.peakgroups(res, missing = 0, span = 0.4)
    res <- as(od_3, "xcmsSet")
    checkTrue(any(unlist(res@rt$raw) != unlist(res@rt$corrected)))
    checkEquals(res@rt$corrected, res_2@rt$corrected)
    checkEquals(features(od_3), peaks(res))
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
    ## Check features
    fd$features <- 4
    checkTrue(!is.logical(xcms:::validateMsFeatureData(fd)))
    fdm <- matrix(ncol = 3, nrow = 5)
    colnames(fdm) <- c("a", "b", "sample")
    fd$features <- fdm
    checkTrue(!is.logical(xcms:::validateMsFeatureData(fd)))
    rm("features", envir = fd)
    ## featureGroups
    fd$features <- xs_2@peaks
    fd$featureGroups <- 4
    checkTrue(!is.logical(xcms:::validateMsFeatureData(fd)))
    fg <- DataFrame(fdm)
    fd$featureGroups <- fg
    checkTrue(!is.logical(xcms:::validateMsFeatureData(fd)))
    fg <- DataFrame(xs_2@groups)
    fg$featureidx <- xs_2@groupidx
    fg_2 <- fg
    fg_2$mzmin <- "a"
    fd$featureGroups <- fg_2
    checkTrue(!is.logical(xcms:::validateMsFeatureData(fd)))
    fg_2 <- fg
    fg_2$featureidx[[1]] <- c(50000, 3)
    fd$featureGroups <- fg_2
    checkTrue(!is.logical(xcms:::validateMsFeatureData(fd)))
    ## adjustedRtime
    fd$featureGroups <- fg
    fd$adjustedRtime <- 4
    checkTrue(!is.logical(xcms:::validateMsFeatureData(fd)))
    fd$adjustedRtime <- list(1:5, "b")
    checkTrue(!is.logical(xcms:::validateMsFeatureData(fd)))
    ## Now check that we pass if we put all correct data into the object:
    fd <- new("MsFeatureData")
    fd$features <- xs_2@peaks
    checkTrue(xcms:::validateMsFeatureData(fd))
    fd$adjustedRtime <- xs_2@rt$corrected
    checkTrue(xcms:::validateMsFeatureData(fd))
    fg <- DataFrame(xs_2@groups)
    fg$featureidx <- xs_2@groupidx
    checkTrue(xcms:::validateMsFeatureData(fd))
}

test_MsFeatureData_class_accessors <- function() {
    fd <- new("MsFeatureData")
    library(S4Vectors)
    checkTrue(!hasDetectedFeatures(fd))
    checkTrue(!hasAdjustedRtime(fd))
    checkTrue(!hasAlignedFeatures(fd))
    suppressWarnings(checkEquals(features(fd), NULL))
    suppressWarnings(checkEquals(featureGroups(fd), NULL))
    suppressWarnings(checkEquals(adjustedRtime(fd), NULL))
    ## features
    features(fd) <- xs_2@peaks
    checkTrue(hasDetectedFeatures(fd))
    checkEquals(features(fd), xs_2@peaks)
    ## featureGroups
    fg <- DataFrame(xs_2@groups)
    fg$featureidx <- xs_2@groupidx
    featureGroups(fd) <- fg
    checkTrue(hasAlignedFeatures(fd))
    checkEquals(featureGroups(fd), fg)
    ## adjustedRtime
    adjustedRtime(fd) <- xs_2@rt$corrected
    checkTrue(hasAdjustedRtime(fd))
    checkEquals(adjustedRtime(fd), xs_2@rt$corrected)
}

############################################################
## Test getEIC alternatives.
dontrun_getEIC_alternatives <- function() {
    library(xcms)

    fls <- c(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
             system.file('cdf/KO/ko16.CDF', package = "faahKO"),
             system.file('cdf/KO/ko18.CDF', package = "faahKO"))
    od <- readMSData2(fls)
    cwp <- CentWaveParam(noise = 10000, snthresh = 40)
    od_x <- detectFeatures(od, param = cwp)
    xs <- as(od_x, "xcmsSet")
    sampclass(xs) <- rep("KO", 3)
    xs_2 <- group(xs)
    suppressWarnings(
        xs_2 <- retcor(xs_2)
    )
    xs_2 <- group(xs_2)

    rtr <- groups(xs_2)[1:5, c("rtmin", "rtmax")]
    mzr <- groups(xs_2)[1:5, c("mzmin", "mzmax")]

    ## Extract the EIC:
    system.time(
        eic <- getEIC(xs_2, rtrange = rtr, mzrange = mzr, rt = "raw")
    ) ## 3.7sec

    ## Now try to do the same using MSnbase stuff.
    rts <- rtime(od_x)
    idx <- apply(rtr, MARGIN = 1, function(z) {
        which(rts >= z[1] & rts <= z[2])
    })
    ## Subset the od_x
    od_ss <- od[sort(unique(unlist(idx)))]
    system.time(
        scts <- spectra(od_ss)
    ) ## 0.7 secs.

    ## This is somewhat similar to the getEIC, just that it extracts for each
    ## mz/rt range pair a data.frame with rt, mz, intensity per sample.
    ## This version works on a single rt/mz range pair at a time.
    ## CHECK:
    ## 1) mz range outside.
    ## 2) rt range outside.
    extractMsData <- function(x, rtrange, mzrange) {
        ## Subset the OnDiskMSnExp
        fns <- fileNames(x)
        tmp <- filterMz(filterRt(x, rt = rtrange), mz = mzrange)
        fromF <- match(fileNames(tmp), fns)
        ## Now extract mz-intensity pairs from each spectrum.
        ## system.time(
        ## suppressWarnings(
        ##     dfs <- spectrapply(tmp, as.data.frame)
        ## )
        ## ) ## 0.73sec
        ## system.time(
        suppressWarnings(
            dfs <- spectrapply(tmp, function(z) {
                if (peaksCount(z))
                    return(data.frame(rt = rep_len(rtime(z), length(z@mz)),
                                      as.data.frame(z)))
                else
                    return(data.frame(rt = numeric(), mz = numeric(),
                                      i = integer()))
            })
        )
        ## I should check if returning just the spectra is not faster!
        ## ) ## 0.701
        ## dfs[] <- mapply(FUN = function(y, z) {
        ##     return(cbind(rt = rep.int(z, nrow(y)), y))
        ## }, y = dfs, z = rtime(tmp), SIMPLIFY = FALSE, USE.NAMES = FALSE)
        res <- vector(mode = "list", length = length(fns))
        res[fromF] <- split(dfs, f = fromFile(tmp))
        return(lapply(res, do.call, what = rbind))
    }

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
    ## 1) Subset the od selecting all spectra that fall into the rt ranges.
    ## 2) Work on that subsetted od: load into memory.
    ## 3) loop over the rtrange and mzrange.


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


