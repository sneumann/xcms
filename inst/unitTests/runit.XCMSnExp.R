## tests related to the new XCMSnExp object.
library(RUnit)

cwp <- CentWaveParam(noise = 10000, snthresh = 40)
## od <- filterRt(od, rt = c(3000, 4000))
od_x <- detectFeatures(faahko_od, param = cwp)
xs <- xcmsSet(faahko_3_files, profparam = list(step = 0), method = "centWave",
              noise = 10000, snthresh = 40)
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
    checkEquals(adjustedRtime(xod), xs_2@rt$corrected)
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
    ##checkTrue(FALSE)
    .checkCreationOfEmptyObject()

    res <- dropFeatures(od_x)
    checkTrue(!hasDetectedFeatures(res))

    ## Add grouping.
    od_2 <- od_x
    od_2 <- xcms:::addProcessHistory(od_2,
                                     xcms:::ProcessHistory(fileIndex. = 1:3,
                                                           type = xcms:::.PROCSTEP.FEATURE.ALIGNMENT))
    library(S4Vectors)
    fd <- DataFrame(xs_2@groups)
    fd$featureidx <- xs_2@groupidx
    featureGroups(od_2) <- fd
    ## Add retention time adjustment.
    od_3 <- od_2
    od_3 <- xcms:::addProcessHistory(od_3,
                                     xcms:::ProcessHistory(fileIndex. = 1:3,
                                                           type = xcms:::.PROCSTEP.RTIME.CORRECTION))
    adjustedRtime(od_3) <- xs_2@rt$corrected
    ## and grouping
    od_4 <- od_3
    od_4 <- xcms:::addProcessHistory(od_4,
                                     xcms:::ProcessHistory(fileIndex. = 1:3,
                                                           type = xcms:::.PROCSTEP.FEATURE.ALIGNMENT))
    ## Tests
    ## drop features:
    res <- dropFeatures(od_2)
    checkTrue(!hasDetectedFeatures(res))
    checkTrue(!hasAlignedFeatures(res))
    checkTrue(!hasAdjustedRtime(res))
    types <- unlist(lapply(processHistory(res), processType))
    checkTrue(!any(types == xcms:::.PROCSTEP.FEATURE.DETECTION))
    checkTrue(!any(types == xcms:::.PROCSTEP.FEATURE.ALIGNMENT))
    checkTrue(!any(types == xcms:::.PROCSTEP.RTIME.CORRECTION))

    res <- dropFeatures(od_3)
    checkTrue(!hasDetectedFeatures(res))
    checkTrue(!hasAlignedFeatures(res))
    checkTrue(!hasAdjustedRtime(res))
    types <- unlist(lapply(processHistory(res), processType))
    checkTrue(!any(types == xcms:::.PROCSTEP.FEATURE.DETECTION))
    checkTrue(!any(types == xcms:::.PROCSTEP.FEATURE.ALIGNMENT))
    checkTrue(!any(types == xcms:::.PROCSTEP.RTIME.CORRECTION))

    res <- dropFeatures(od_4)
    checkTrue(!hasDetectedFeatures(res))
    checkTrue(!hasAlignedFeatures(res))
    checkTrue(!hasAdjustedRtime(res))
    types <- unlist(lapply(processHistory(res), processType))
    checkTrue(!any(types == xcms:::.PROCSTEP.FEATURE.DETECTION))
    checkTrue(!any(types == xcms:::.PROCSTEP.FEATURE.ALIGNMENT))
    checkTrue(!any(types == xcms:::.PROCSTEP.RTIME.CORRECTION))

    ## Drop feature groups
    ## adjusted retention times are always dropped.
    res <- dropFeatureGroups(od_x)
    checkTrue(!hasAlignedFeatures(res))
    checkTrue(!hasAdjustedRtime(res))
    checkTrue(hasDetectedFeatures(res))
    types <- unlist(lapply(processHistory(res), processType))
    checkTrue(any(types == xcms:::.PROCSTEP.FEATURE.DETECTION))
    checkTrue(!any(types == xcms:::.PROCSTEP.FEATURE.ALIGNMENT))
    checkTrue(!any(types == xcms:::.PROCSTEP.RTIME.CORRECTION))

    res <- dropFeatureGroups(od_2)
    checkTrue(!hasAlignedFeatures(res))
    checkTrue(!hasAdjustedRtime(res))
    checkTrue(hasDetectedFeatures(res))
    types <- unlist(lapply(processHistory(res), processType))
    checkTrue(any(types == xcms:::.PROCSTEP.FEATURE.DETECTION))
    checkTrue(!any(types == xcms:::.PROCSTEP.FEATURE.ALIGNMENT))
    checkTrue(!any(types == xcms:::.PROCSTEP.RTIME.CORRECTION))

    res <- dropFeatureGroups(od_3)
    checkTrue(!hasAlignedFeatures(res))
    checkTrue(!hasAdjustedRtime(res))
    checkTrue(hasDetectedFeatures(res))
    types <- unlist(lapply(processHistory(res), processType))
    checkTrue(any(types == xcms:::.PROCSTEP.FEATURE.DETECTION))
    checkTrue(!any(types == xcms:::.PROCSTEP.FEATURE.ALIGNMENT))
    checkTrue(!any(types == xcms:::.PROCSTEP.RTIME.CORRECTION))

    res <- dropFeatureGroups(od_4)
    checkTrue(!hasAlignedFeatures(res))
    checkTrue(!hasAdjustedRtime(res))
    checkTrue(hasDetectedFeatures(res))
    types <- unlist(lapply(processHistory(res), processType))
    checkTrue(any(types == xcms:::.PROCSTEP.FEATURE.DETECTION))
    checkTrue(!any(types == xcms:::.PROCSTEP.FEATURE.ALIGNMENT))
    checkTrue(!any(types == xcms:::.PROCSTEP.RTIME.CORRECTION))

    ## Drop aligned rtime
    ## Feature alignments are only dropped if they were performed after the
    ## retention time adjustments.
    res <- dropAdjustedRtime(od_x)
    checkTrue(!hasAlignedFeatures(res))
    checkTrue(!hasAdjustedRtime(res))
    checkTrue(hasDetectedFeatures(res))
    types <- unlist(lapply(processHistory(res), processType))
    checkTrue(any(types == xcms:::.PROCSTEP.FEATURE.DETECTION))
    checkTrue(!any(types == xcms:::.PROCSTEP.FEATURE.ALIGNMENT))
    checkTrue(!any(types == xcms:::.PROCSTEP.RTIME.CORRECTION))

    res <- dropAdjustedRtime(od_2)
    checkTrue(hasAlignedFeatures(res))
    checkTrue(!hasAdjustedRtime(res))
    checkTrue(hasDetectedFeatures(res))
    types <- unlist(lapply(processHistory(res), processType))
    checkTrue(any(types == xcms:::.PROCSTEP.FEATURE.DETECTION))
    checkTrue(any(types == xcms:::.PROCSTEP.FEATURE.ALIGNMENT))
    checkTrue(!any(types == xcms:::.PROCSTEP.RTIME.CORRECTION))

    ## Now dropping the adjusted retention time, but not the feature alignment.
    res <- dropAdjustedRtime(od_3)
    checkTrue(hasAlignedFeatures(res))
    checkTrue(!hasAdjustedRtime(res))
    checkTrue(hasDetectedFeatures(res))
    checkTrue(hasAdjustedRtime(od_3))
    types <- unlist(lapply(processHistory(res), processType))
    checkTrue(any(types == xcms:::.PROCSTEP.FEATURE.DETECTION))
    checkTrue(any(types == xcms:::.PROCSTEP.FEATURE.ALIGNMENT))
    checkTrue(!any(types == xcms:::.PROCSTEP.RTIME.CORRECTION))

    ## Drop adjusted retention time AND feature alignment.
    res <- dropAdjustedRtime(od_4)
    checkTrue(!hasAlignedFeatures(res))
    checkTrue(!hasAdjustedRtime(res))
    checkTrue(hasDetectedFeatures(res))
    checkTrue(hasAdjustedRtime(od_4))
    types <- unlist(lapply(processHistory(res), processType))
    checkTrue(any(types == xcms:::.PROCSTEP.FEATURE.DETECTION))
    checkTrue(!any(types == xcms:::.PROCSTEP.FEATURE.ALIGNMENT))
    checkTrue(!any(types == xcms:::.PROCSTEP.RTIME.CORRECTION))
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

dontrun_test_XCMSnExp <- function() {
    ## MsFeatureData

    ## Some checks how to deal and best represent the data from an xcmsSet

    ## Cast from OnDiskMSnExp.
}

