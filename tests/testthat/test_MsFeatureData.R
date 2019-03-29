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
    ## chromPeaks
    cp <- matrix(nrow = 3, ncol = length(xcms:::.REQ_PEAKS_COLS))
    colnames(cp) <- xcms:::.REQ_PEAKS_COLS
    expect_true(is.character(xcms:::.validChromPeaksMatrix(cp)))
    cp[3] <- 3.4
    expect_true(xcms:::.validChromPeaksMatrix(cp))
    rownames(cp) <- letters[1:nrow(cp)]
    chromPeaks(fd) <- cp
    expect_true(validObject(fd))
    ## featureDefinitions
    xcms:::chromPeakData(fd) <- 5
    expect_error(validObject(fd), "is supposed to be a 'DataFrame'")
    fdef <- DataFrame(ms_level = rep(1L, nrow(cp)), row.names = rownames(cp))
    xcms:::chromPeakData(fd) <- fdef
    expect_true(validObject(fd))
    xcms:::chromPeakData(fd)$ms_level <- "a"
    expect_error(validObject(fd), "column 'ms_level' should contain")
    xcms:::chromPeakData(fd)$ms_level <- 1L
    rownames(xcms:::chromPeakData(fd)) <- 1:nrow(cp)
    expect_error(validObject(fd), "rownames differ")
    rm("chromPeaks", envir = fd)
    expect_error(validObject(fd), "'chromPeakData' present but 'chromPeaks'")
    rm("chromPeakData", envir = fd)

    ## Additional tests.
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
