test_that(".obiwarp works", {

    xs <- faahko_xs
    od <- faahko_od
    xod <- faahko_xod
    ## Feature alignment on those:
    ## object <- findChromPeaks(faahko_od, param = CentWaveParam(noise = 10000,
    ##                                                           snthresh = 40))
    prm <- ObiwarpParam(binSize = 1)
    xs_2 <- retcor.obiwarp(xs, profStep = binSize(prm))
    expect_equal(xs_2@rt$raw[[2]], xs_2@rt$corrected[[2]])
    expect_true(sum(xs_2@rt$raw[[1]] != xs_2@rt$corrected[[1]]) > 500)
    expect_true(sum(xs_2@rt$raw[[3]] != xs_2@rt$corrected[[3]]) > 500)
    
    ## And the OnDiskMSnExp implementation:
    res <- .obiwarp(od, param = prm)
    expect_equal(xs_2@rt$corrected, res)
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
    xs_2 <- retcor.obiwarp(xs, profStep = binSize(prm), center = centerSample(prm))
    expect_equal(xs_2@rt$raw[[centerSample(prm)]],
                 xs_2@rt$corrected[[centerSample(prm)]])
    res <- .obiwarp(od, param = prm)
    expect_equal(xs_2@rt$corrected, res)
    ## change some settings
    gapInit(prm) <- 3.1
    gapExtend(prm) <- 0.9
    xs_2 <- retcor.obiwarp(xs, profStep = binSize(prm), gapInit = gapInit(prm),
                           center = centerSample(prm), gapExtend = gapExtend(prm))
    expect_equal(xs_2@rt$raw[[centerSample(prm)]],
                 xs_2@rt$corrected[[centerSample(prm)]])
    res <- .obiwarp(od, param = prm)
    expect_equal(xs_2@rt$corrected, res)
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


