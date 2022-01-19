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

    ## Manually specify center Sample
    centerSample(prm) <- 3
    xs_2 <- retcor.obiwarp(xs, profStep = binSize(prm), center = centerSample(prm))
    expect_equal(xs_2@rt$raw[[centerSample(prm)]],
                 xs_2@rt$corrected[[centerSample(prm)]])
    res <- .obiwarp(od, param = prm)
    expect_equal(res[[3]], unname(rtime(xod, bySample = TRUE)[[3]]))
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
