test_that("MsFeatureData class validation works", {
    skip_on_os(os = "windows", arch = "i386")

    fd <- new("MsFeatureData")
    ## Check error for wrong elements.
    fd$a <- 5
    expect_true(!is.logical(validateMsFeatureData(fd)))
    rm("a", envir = fd)
    ## Check chromPeaks
    fd$chromPeaks <- 4
    expect_true(!is.logical(validateMsFeatureData(fd)))
    fdm <- matrix(ncol = 3, nrow = 5)
    colnames(fdm) <- c("a", "b", "sample")
    fd$chromPeaks <- fdm
    expect_true(!is.logical(validateMsFeatureData(fd)))
    rm("chromPeaks", envir = fd)
    ## chromPeaks
    cp <- matrix(nrow = 3, ncol = length(.REQ_PEAKS_COLS))
    colnames(cp) <- .REQ_PEAKS_COLS
    expect_true(is.character(.validChromPeaksMatrix(cp)))
    cp[3] <- 3.4
    expect_true(.validChromPeaksMatrix(cp))
    rownames(cp) <- letters[1:nrow(cp)]
    chromPeaks(fd) <- cp
    expect_true(validObject(fd))
    ## featureDefinitions
    chromPeakData(fd) <- 5
    expect_error(validObject(fd), "is supposed to be a 'DataFrame'")
    fdef <- DataFrame(ms_level = rep(1L, nrow(cp)),
                      is_filled = rep(FALSE, nrow(cp)),
                      row.names = rownames(cp))
    chromPeakData(fd) <- fdef
    expect_true(validObject(fd))
    chromPeakData(fd)$ms_level <- "a"
    expect_error(validObject(fd), "column 'ms_level' should contain")
    chromPeakData(fd)$ms_level <- 1L
    rownames(chromPeakData(fd)) <- 1:nrow(cp)
    expect_error(validObject(fd), "rownames differ")
    rm("chromPeaks", envir = fd)
    expect_error(validObject(fd), "'chromPeakData' present but 'chromPeaks'")
    rm("chromPeakData", envir = fd)

    ## Additional tests.
    fd$chromPeaks <- chromPeaks(xod_x)
    fd$featureDefinitions <- 4
    expect_true(!is.logical(validateMsFeatureData(fd)))
    fg <- featureDefinitions(xod_xgrg)
    fd$featureDefinitions <- fg[, 1:8]
    expect_true(!is.logical(validateMsFeatureData(fd)))
    fg_2 <- fg
    fg_2$mzmin <- "a"
    fd$featureDefinitions <- fg_2
    expect_true(!is.logical(validateMsFeatureData(fd)))
    fg_2 <- fg
    fg_2$peakidx[[1]] <- c(50000, 3)
    fd$featureDefinitions <- fg_2
    expect_true(!is.logical(validateMsFeatureData(fd)))
    ## adjustedRtime
    fd$featureDefinitions <- fg
    fd$adjustedRtime <- 4
    expect_true(!is.logical(validateMsFeatureData(fd)))
    fd$adjustedRtime <- list(1:5, "b")
    expect_true(!is.logical(validateMsFeatureData(fd)))
    ## Now check that we pass if we put all correct data into the object:
    fd <- new("MsFeatureData")
    fd$chromPeaks <- chromPeaks(xod_xgrg)
    expect_true(length(validateMsFeatureData(fd)) == 0)
    fd$adjustedRtime <- xod_xgrg@msFeatureData$adjustedRtime
    expect_true(length(validateMsFeatureData(fd)) == 0)
    fd$featureDefinitions <- featureDefinitions(xod_xgrg)
    expect_true(length(validateMsFeatureData(fd)) == 0)
})

test_that("MsFeatureData class_accessors work", {
    skip_on_os(os = "windows", arch = "i386")

    fd <- new("MsFeatureData")
    expect_true(!hasChromPeaks(fd))
    expect_true(!hasAdjustedRtime(fd))
    expect_true(!hasFeatures(fd))
    expect_equal(chromPeaks(fd), NULL)
    expect_warning(expect_equal(featureDefinitions(fd), DataFrame()))
    expect_warning(expect_equal(adjustedRtime(fd), NULL))
    ## chromPeaks
    chromPeaks(fd) <- chromPeaks(xod_xgrg)
    chromPeakData(fd) <- chromPeakData(xod_xgrg)
    expect_true(hasChromPeaks(fd))
    expect_false(hasChromPeaks(fd, msLevel = 2L))
    expect_equal(chromPeaks(fd), chromPeaks(xod_xgrg))
    ## featureDefinitions
    featureDefinitions(fd) <- featureDefinitions(xod_xgrg)
    expect_true(hasFeatures(fd))
    expect_equal(featureDefinitions(fd), featureDefinitions(xod_xgrg))
    expect_false(hasFeatures(fd, msLevel = 2L))
    expect_true(nrow(featureDefinitions(fd, msLevel = 2L)) == 0)
    ## adjustedRtime
    adjustedRtime(fd) <- adjustedRtime(xod_xgrg)
    expect_error(validObject(fd))
    adjustedRtime(fd) <- xod_xgrg@msFeatureData$adjustedRtime
    expect_true(hasAdjustedRtime(fd))
    expect_equal(adjustedRtime(fd), xod_xgrg@msFeatureData$adjustedRtime)
})
