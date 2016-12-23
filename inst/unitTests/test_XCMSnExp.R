## tests related to the new XCMSnExp object.
library(RUnit)

library(faahKO)
fs <- c(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
        system.file('cdf/KO/ko16.CDF', package = "faahKO"),
        system.file('cdf/KO/ko18.CDF', package = "faahKO"))
## library(msdata)
## f <- msdata::proteomics(full.names = TRUE, pattern = "TMT_Erwinia")
## mzf <- c(system.file("microtofq/MM14.mzML", package = "msdata"),
##          system.file("microtofq/MM8.mzML", package = "msdata"))
od <- readMSData2(fs)
cwp <- CentWaveParam(noise = 10000, snthresh = 40)
od <- filterRt(od, rt = c(3000, 4000))
od_x <- detectFeatures(od, param = cwp)
xs <- xcmsSet(fs, profparam = list(step = 0), method = "centWave",
              noise = 10000, snthresh = 40)
xs_2 <- group(xs)
suppressWarnings(
    xs_2 <- retcor(xs_2)
)
xs_2 <- group(xs_2)

test_XCMSnExp_class <- function() {
    ## Basic contructor.
    xs <- new("XCMSnExp")
    xs@.processHistory <- list("a")
    checkException(validObject(xs))
    xs@.processHistory <- list(xcms:::ProcessHistory())
    checkTrue(validObject(xs))
    xod <- as(od, "XCMSnExp")
    checkTrue(validObject(xod))
    ## MsFeatureData error
    xod@msFeatureData$features <- 3
    checkException(validObject(xod))
    xod@msFeatureData$features <- xs_2@peaks
    checkTrue(validObject(xod))
    xod@msFeatureData$features[1, "sample"] <- 40
    checkException(validObject(xod))
    xod@msFeatureData$features[1, "sample"] <- 4
    xod@msFeatureData$adjustedRtime <- xs_2@rt$corrected
    checkTrue(validObject(xod))
    xod@msFeatureData$adjustedRtime[[2]] <- 1:4
    checkException(validObject(xod))
    xod@msFeatureData$adjustedRtime <- xs_2@rt$corrected[1:3]
    checkException(validObject(xod))
}

test_XCMSnExp_class_accessors <- function() {
    xs <- new("XCMSnExp")
    checkTrue(!hasAdjustedRtime(xs))
    checkTrue(!hasAlignedFeatures(xs))
    checkTrue(!hasDetectedFeatures(xs))
    ## Filling with data...
    xod <- as(od, "XCMSnExp")
    ## features
    checkTrue(!hasDetectedFeatures(xod))
    features(xod) <- xs_2@peaks
    checkTrue(hasDetectedFeatures(xod))
    checkEquals(features(xod), xs_2@peaks)
    checkException(features(xod) <- 4)
    ## featureGroups
    checkTrue(!hasAlignedFeatures(xod))
    library(S4Vectors)
    fd <- DataFrame(xs_2@groups)
    fd$featureidx <- xs_2@groupidx
    featureGroups(xod) <- fd
    checkTrue(hasAlignedFeatures(xod))
    checkEquals(featureGroups(xod), fd)
    ## adjustedRtime
    checkTrue(!hasAdjustedRtime(xod))
    adjustedRtime(xod) <- xs_2@rt$corrected
    checkTrue(hasAdjustedRtime(xod))
    checkEquals(adjustedRtime(xod), xs_2@rt$corrected)
    ## rtime
    rtm <- rtime(xod)
    checkEquals(rtm, split(rtime(od), fromFile(od)))
}

test_XCMSnExp_processHistory <- function() {
    ph <- xcms:::ProcessHistory(fileIndex. = 2, info. = "For file 2")
    ph_2 <- xcms:::ProcessHistory(fileIndex. = 1:2, info. = "For files 1 to 2")
    xod <- as(od, "XCMSnExp")
    xod@.processHistory <- list(ph, ph_2)
    checkEquals(processHistory(xod), list(ph, ph_2))
    checkEquals(processHistory(xod, fileIndex = 2), list(ph, ph_2))
    checkEquals(processHistory(xod, fileIndex = 1), list(ph))
    checkException(processHistory(xod, fileIndex = 5))

    ph_3 <- xcms:::XProcessHistory(fileIndex = 1, param = CentWaveParam())
    xod <- xcms:::addProcessHistory(xod, ph_3)
    checkEquals(length(processHistory(xod)), 3)
    checkEquals(processHistory(xod)[[3]], ph_3)
    checkEquals(processHistory(xod, fileIndex = 1), list(ph_2, ph_3))
}

## Testing all inherited methods for XCMSnExp classes that perform data
## manipulations - for all of these we have to ensure that eventual preprocessing
## results are removed.
test_XCMSnExp_inherited_methods <- function() {
    ## [
    tmp_1 <- od[1:10]
    suppressWarnings(
        tmp_2 <- od_x[1:10]
    )
    checkTrue(length(processHistory(tmp_2)) == 0)
    checkTrue(!hasDetectedFeatures(tmp_2))
    tmp_1@processingData <- new("MSnProcess")
    tmp_2@processingData <- new("MSnProcess")
    checkEquals(tmp_1, as(tmp_2, "OnDiskMSnExp"))
    ## bin
    tmp_1 <- bin(od)
    suppressWarnings(
        tmp_2 <- bin(od_x)
    )
    checkTrue(length(processHistory(tmp_2)) == 0)
    checkTrue(!hasDetectedFeatures(tmp_2))
    tmp_1@processingData <- new("MSnProcess")
    tmp_2@processingData <- new("MSnProcess")
    checkEquals(tmp_1, as(tmp_2, "OnDiskMSnExp"))
    ## clean
    tmp_1 <- clean(od)
    suppressWarnings(
        tmp_2 <- clean(od_x)
    )
    checkTrue(length(processHistory(tmp_2)) == 0)
    checkTrue(!hasDetectedFeatures(tmp_2))
    tmp_1@processingData <- new("MSnProcess")
    tmp_2@processingData <- new("MSnProcess")
    checkEquals(tmp_1, as(tmp_2, "OnDiskMSnExp"))
    ## filterAcquisitionNum
    tmp_1 <- filterAcquisitionNum(od)
    suppressWarnings(
        tmp_2 <- filterAcquisitionNum(od_x)
    )
    checkTrue(length(processHistory(tmp_2)) == 0)
    checkTrue(!hasDetectedFeatures(tmp_2))
    tmp_1@processingData <- new("MSnProcess")
    tmp_2@processingData <- new("MSnProcess")
    checkEquals(tmp_1, as(tmp_2, "OnDiskMSnExp"))
    ## filterMsLevel
    tmp_1 <- filterMsLevel(od)
    suppressWarnings(
        tmp_2 <- filterMsLevel(od_x)
    )
    checkTrue(length(processHistory(tmp_2)) == 0)
    checkTrue(!hasDetectedFeatures(tmp_2))
    tmp_1@processingData <- new("MSnProcess")
    tmp_2@processingData <- new("MSnProcess")
    checkEquals(tmp_1, as(tmp_2, "OnDiskMSnExp"))
    ## normalize
    tmp_1 <- normalize(od)
    suppressWarnings(
        tmp_2 <- normalize(od_x)
    )
    checkTrue(length(processHistory(tmp_2)) == 0)
    checkTrue(!hasDetectedFeatures(tmp_2))
    tmp_1@processingData <- new("MSnProcess")
    tmp_2@processingData <- new("MSnProcess")
    checkEquals(tmp_1, as(tmp_2, "OnDiskMSnExp"))
    ## pickPeaks
    tmp_1 <- pickPeaks(od)
    suppressWarnings(
        tmp_2 <- pickPeaks(od_x)
    )
    checkTrue(length(processHistory(tmp_2)) == 0)
    checkTrue(!hasDetectedFeatures(tmp_2))
    tmp_1@processingData <- new("MSnProcess")
    tmp_2@processingData <- new("MSnProcess")
    checkEquals(tmp_1, as(tmp_2, "OnDiskMSnExp"))
    ## removePeaks
    tmp_1 <- removePeaks(od)
    suppressWarnings(
        tmp_2 <- removePeaks(od_x)
    )
    checkTrue(length(processHistory(tmp_2)) == 0)
    checkTrue(!hasDetectedFeatures(tmp_2))
    tmp_1@processingData <- new("MSnProcess")
    tmp_2@processingData <- new("MSnProcess")
    checkEquals(tmp_1, as(tmp_2, "OnDiskMSnExp"))
    ## smooth
    tmp_1 <- smooth(od)
    suppressWarnings(
        tmp_2 <- smooth(od_x)
    )
    checkTrue(length(processHistory(tmp_2)) == 0)
    checkTrue(!hasDetectedFeatures(tmp_2))
    tmp_1@processingData <- new("MSnProcess")
    tmp_2@processingData <- new("MSnProcess")
    checkEquals(tmp_1, as(tmp_2, "OnDiskMSnExp"))
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

