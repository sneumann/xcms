test_that(".obiwarp works", {

    od <- faahko_od
    xod <- faahko_xod
    prm <- ObiwarpParam(binSize = 1)

    raw_rt <- split(rtime(od), fromFile(od))
    ## And the OnDiskMSnExp implementation:
    res <- .obiwarp(od, param = prm)
    expect_true(all(raw_rt[[1]] != res[[1]]))
    expect_equal(res[[2]], unname(raw_rt[[2]]))
    res_2 <- adjustRtime(od, param = prm)
    res_3 <- adjustRtime(xod, param = prm)
    expect_equal(adjustedRtime(res_3), res_2)
    expect_equal(lapply(adjustedRtime(res_3, bySample = TRUE), unname), res)
    expect_equal(adjustedRtime(res_3), res_2)
    ## Check if peaks were corrected correctly
    expect_true(sum(chromPeaks(res_3)[, "rt"] == chromPeaks(xod)) <
                nrow(chromPeaks(res_3)) / 2)
    ## Dropping the adjusted rtime on these
    expect_true(hasAdjustedRtime(res_3))
    tmp <- dropAdjustedRtime(res_3)
    expect_equal(chromPeaks(tmp), chromPeaks(xod))

    ## File issue on that! retcor.obiwarp does use round for the adjustment of
    ## the peak!
    ## -> issue #122
    ## expect_equal(chromPeaks(res_3), peaks(xs_2))

    ## Manually specify center Sample
    centerSample(prm) <- 3
    res <- .obiwarp(od, param = prm)
    expect_equal(res[[3]], unname(raw_rt[[3]]))
    expect_true(all(res[[2]] != raw_rt[[2]]))

    ## With subset.
    prm <- ObiwarpParam(binSize = 1, subset = c(1, 3))
    res <- .obiwarp(od, param = prm)
    expect_equal(res[[1]], unname(raw_rt[[1]]))
    expect_true(all(res[[2]] != unname(raw_rt[[2]])))
    expect_true(all(res[[3]] != unname(raw_rt[[3]])))

    prm <- ObiwarpParam(binSize = 1, subset = c(2, 3))
    res <- .obiwarp(od, param = prm)
    expect_equal(res[[2]], unname(raw_rt[[2]]))
    expect_true(sum(res[[1]] == unname(raw_rt[[1]])) > 500)
    expect_true(all(res[[3]] != unname(raw_rt[[3]])))
})

test_that(".concatenate_OnDiskMSnExp works", {
    od1 <- readMSData(faahko_3_files[1], mode = "onDisk")
    od2 <- readMSData(faahko_3_files[2:3], mode = "onDisk")
    res <- .concatenate_OnDiskMSnExp(od1, od2)
    expect_equal(fileNames(faahko_od), fileNames(res))
    expect_equal(experimentData(faahko_od), experimentData(res))
    expect_equal(fData(faahko_od), fData(res))
    expect_equal(pData(faahko_od), pData(res))
    ## expect_equal(spectra(faahko_od), spectra(res))
    ## Checking data manipulations
    od1 <- smooth(od1)
    od1 <- pickPeaks(od1)
    od2 <- smooth(od2)
    expect_warning(res <- .concatenate_OnDiskMSnExp(od1, od2))
    expect_true(length(res@spectraProcessingQueue) == 0)
    expect_equal(fileNames(faahko_od), fileNames(res))
    expect_equal(experimentData(faahko_od), experimentData(res))
    expect_equal(pData(faahko_od), pData(res))

    od2 <- pickPeaks(od2)
    res <- .concatenate_OnDiskMSnExp(od1, od2)
    expect_true(length(res@spectraProcessingQueue) == 2)
    expect_equal(fileNames(faahko_od), fileNames(res))
    expect_equal(experimentData(faahko_od), experimentData(res))
    expect_equal(pData(faahko_od), pData(res))

    expect_equal(.concatenate_OnDiskMSnExp(od1), od1)
})

test_that(".split_by_file works", {
    a <- lapply(seq_along(fileNames(od_x)), filterFile, object = od_x)
    b <- .split_by_file(od_x)
    expect_equal(a[[2]][[19]], b[[2]][[19]])
    expect_equal(spectra(a[[3]]), spectra(b[[3]]))
})
