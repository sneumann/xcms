test_that("MsFeatureData class validation works", {
    fd <- new("MsFeatureData")
    ## Check error for wrong elements.
    fd$a <- 5
    expect_true(!is.logical(xcms:::validateMsFeatureData(fd)))
    rm("a", envir = fd)
    ## Check chromPeaks
    fd$chromPeaks <- 4
    expect_true(!is.logical(xcms:::validateMsFeatureData(fd)))
    fdm <- matrix(ncol = 3, nrow = 5)
    colnames(fdm) <- c("a", "b", "sample")
    fd$chromPeaks <- fdm
    expect_true(!is.logical(xcms:::validateMsFeatureData(fd)))
    rm("chromPeaks", envir = fd)
    ## featureDefinitions
    fd$chromPeaks <- chromPeaks(xod_x)
    fd$featureDefinitions <- 4
    expect_true(!is.logical(xcms:::validateMsFeatureData(fd)))
    fg <- featureDefinitions(xod_xgrg)
    fd$featureDefinitions <- fg[, 1:8]
    expect_true(!is.logical(xcms:::validateMsFeatureData(fd)))
    fg_2 <- fg
    fg_2$mzmin <- "a"
    fd$featureDefinitions <- fg_2
    expect_true(!is.logical(xcms:::validateMsFeatureData(fd)))
    fg_2 <- fg
    fg_2$peakidx[[1]] <- c(50000, 3)
    fd$featureDefinitions <- fg_2
    expect_true(!is.logical(xcms:::validateMsFeatureData(fd)))
    ## adjustedRtime
    fd$featureDefinitions <- fg
    fd$adjustedRtime <- 4
    expect_true(!is.logical(xcms:::validateMsFeatureData(fd)))
    fd$adjustedRtime <- list(1:5, "b")
    expect_true(!is.logical(xcms:::validateMsFeatureData(fd)))
    ## Now check that we pass if we put all correct data into the object:
    fd <- new("MsFeatureData")
    fd$chromPeaks <- chromPeaks(xod_xgrg)
    expect_true(length(xcms:::validateMsFeatureData(fd)) == 0)
    fd$adjustedRtime <- xod_xgrg@msFeatureData$adjustedRtime
    expect_true(length(xcms:::validateMsFeatureData(fd)) == 0)
    fd$featureDefinitions <- featureDefinitions(xod_xgrg)
    expect_true(length(xcms:::validateMsFeatureData(fd)) == 0)
})

test_that("MsFeatureData class_accessors work", {
    fd <- new("MsFeatureData")
    expect_true(!hasChromPeaks(fd))
    expect_true(!hasAdjustedRtime(fd))
    expect_true(!hasFeatures(fd))
    expect_warning(expect_equal(chromPeaks(fd), NULL))
    expect_warning(expect_equal(featureDefinitions(fd), NULL))
    expect_warning(expect_equal(adjustedRtime(fd), NULL))
    ## chromPeaks
    chromPeaks(fd) <- chromPeaks(xod_xgrg)
    expect_true(hasChromPeaks(fd))
    expect_equal(chromPeaks(fd), chromPeaks(xod_xgrg))
    ## featureDefinitions
    featureDefinitions(fd) <- featureDefinitions(xod_xgrg)
    expect_true(hasFeatures(fd))
    expect_equal(featureDefinitions(fd), featureDefinitions(xod_xgrg))
    ## adjustedRtime
    expect_error(adjustedRtime(fd) <- adjustedRtime(xod_xgrg))
    adjustedRtime(fd) <- xod_xgrg@msFeatureData$adjustedRtime
    expect_true(hasAdjustedRtime(fd))
    expect_equal(adjustedRtime(fd), xod_xgrg@msFeatureData$adjustedRtime)
})

