test_that("getPeakGroupsRtMatrix works", {
    param <- PeakGroupsParam()
    nSamples <- length(fileNames(xod_xg))
    pkGrp <- .getPeakGroupsRtMatrix(
        peaks = chromPeaks(xod_xg),
        peakIndex = .peakIndex(xod_xg),
        nSamples = nSamples,
        missingSample = nSamples - (nSamples * minFraction(param)),
        extraPeaks = extraPeaks(param)
    )
    ## expect_equal(colnames(pkGrp), colnames(chromPeaks(xod_xg)))
    fts <- featureDefinitions(xod_xg)[rownames(pkGrp), ]
    expect_true(all(pkGrp[, 1] >= fts$rtmin & pkGrp[, 1] <= fts$rtmax))
    expect_true(all(pkGrp[, 2] >= fts$rtmin & pkGrp[, 2] <= fts$rtmax))
    expect_true(all(pkGrp[, 3] >= fts$rtmin & pkGrp[, 3] <= fts$rtmax))
})

test_that("do_adjustRtime_peakGroups works", {
    xs <- faahko
    xsg <- group(xs)
    misSamp <- 1
    xsa <- retcor(xsg, method = "peakgroups", missing = misSamp)
    minFr <- (length(sampnames(xs)) - misSamp) / length(sampnames(xs))
    res <- do_adjustRtime_peakGroups(peaks = peaks(xs),
                                     peakIndex = xsg@groupidx,
                                     rtime = xsg@rt$raw,
                                     minFraction = minFr)
    expect_equal(xsa@rt$corrected, res)
    ## Change settings.
    misSamp <- 3
    xsa <- retcor(xsg, method = "peakgroups", missing = misSamp)
    minFr <- (length(sampnames(xs)) - misSamp) / length(sampnames(xs))
    res <- do_adjustRtime_peakGroups(peaks = peaks(xs),
                                     peakIndex = xsg@groupidx,
                                     rtime = xsg@rt$raw,
                                     minFraction = minFr)
    expect_equal(xsa@rt$corrected, res)
    misSamp <- 2
    xtr <- 2
    xsa <- retcor(xsg, method = "peakgroups", missing = misSamp, extra = xtr)
    minFr <- (length(sampnames(xs)) - misSamp) / length(sampnames(xs))
    res <- do_adjustRtime_peakGroups(peaks = peaks(xs),
                                     peakIndex = xsg@groupidx,
                                     rtime = xsg@rt$raw,
                                     minFraction = minFr, extraPeaks = xtr)
    expect_equal(xsa@rt$corrected, res)
    xsa <- retcor(xsg, method = "peakgroups", missing = misSamp, extra = xtr,
                  smooth = "linear")
    minFr <- (length(sampnames(xs)) - misSamp) / length(sampnames(xs))
    res <- do_adjustRtime_peakGroups(peaks = peaks(xs),
                                     peakIndex = xsg@groupidx,
                                     rtime = xsg@rt$raw,
                                     minFraction = minFr, extraPeaks = xtr,
                                     smooth = "linear")
    expect_equal(xsa@rt$corrected, res)
    xsa <- retcor(xsg, method = "peakgroups", missing = misSamp, extra = xtr,
                  family = "symmetric")
    minFr <- (length(sampnames(xs)) - misSamp) / length(sampnames(xs))
    res <- do_adjustRtime_peakGroups(peaks = peaks(xs),
                                     peakIndex = xsg@groupidx,
                                     rtime = xsg@rt$raw,
                                     minFraction = minFr, extraPeaks = xtr,
                                     family = "symmetric")
    expect_equal(xsa@rt$corrected, res)
    xsa <- retcor(xsg, method = "peakgroups", missing = misSamp, extra = xtr,
                  span = 1)
    minFr <- (length(sampnames(xs)) - misSamp) / length(sampnames(xs))
    res <- do_adjustRtime_peakGroups(peaks = peaks(xs),
                                     peakIndex = xsg@groupidx,
                                     rtime = xsg@rt$raw,
                                     minFraction = minFr, extraPeaks = xtr,
                                     span = 1)
    expect_equal(xsa@rt$corrected, res)
})

test_that("applyRtAdjustment works", {
    xs <- faahko
    ## group em.
    xsg <- group(xs)
    ## align em.
    xsa <- retcor(xsg, method = "peakgroups")
    pksAdj <- .applyRtAdjToChromPeaks(peaks(xsg),
                                      rtraw = xsa@rt$raw,
                                      rtadj = xsa@rt$corrected)
    expect_equal(pksAdj, peaks(xsa))
    ## Reset em.
    pksRaw <- .applyRtAdjToChromPeaks(pksAdj,
                                      rtraw = xsa@rt$corrected,
                                      rtadj = xsa@rt$raw)
    expect_equal(pksRaw, peaks(xsg))
})
