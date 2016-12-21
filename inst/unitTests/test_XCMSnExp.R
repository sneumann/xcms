## tests related to the new XCMSnExp object.
library(RUnit)

library(faahKO)
fs <- c(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
        system.file('cdf/KO/ko16.CDF', package = "faahKO"),
        system.file('cdf/KO/ko18.CDF', package = "faahKO"),
        system.file('cdf/KO/ko19.CDF', package = "faahKO"))
library(msdata)
f <- msdata::proteomics(full.names = TRUE, pattern = "TMT_Erwinia")
mzf <- c(system.file("microtofq/MM14.mzML", package = "msdata"),
         system.file("microtofq/MM8.mzML", package = "msdata"))
library(MSnbase)
od <- readMSData2(fs)
xs <- xcmsSet(fs, profparam = list(step = 1))
xs_2 <- group(xs)
xs_2 <- retcor(xs_2)
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

