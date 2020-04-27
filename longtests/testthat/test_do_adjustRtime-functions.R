test_that("getPeakGroupsRtMatrix works", {
    param <- PeakGroupsParam()
    nSamples <- length(fileNames(xod_xg))
    pkGrp <- .getPeakGroupsRtMatrix(
        peaks = chromPeaks(xod_xg),
        peakIndex = .peakIndex(xod_xg),
        sampleIndex = seq_len(nSamples),
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
    expect_error(do_adjustRtime_peakGroups(peaks = peaks(xs),
                                           peakIndex = xsg@groupidx,
                                           rtime = xsg@rt$raw,
                                           minFraction = minFr,
                                           subset = "4"))
    expect_error(do_adjustRtime_peakGroups(peaks = peaks(xs),
                                           peakIndex = xsg@groupidx,
                                           rtime = xsg@rt$raw,
                                           minFraction = minFr,
                                           subset = 4L))
    expect_error(do_adjustRtime_peakGroups(peaks = peaks(xs),
                                           peakIndex = xsg@groupidx,
                                           rtime = xsg@rt$raw,
                                           minFraction = minFr,
                                           subset = c(1, 2, 5, 14)))

    res <- do_adjustRtime_peakGroups(peaks = peaks(xs),
                                     peakIndex = xsg@groupidx,
                                     rtime = xsg@rt$raw,
                                     minFraction = minFr)
    res_orig <- xcms:::do_adjustRtime_peakGroups_orig(peaks = peaks(xs),
                                     peakIndex = xsg@groupidx,
                                     rtime = xsg@rt$raw,
                                     minFraction = minFr)
    expect_equal(res, res_orig)
    expect_equal(xsa@rt$corrected, res)
    ## Use only a subset.
    res_sub <- do_adjustRtime_peakGroups(peaks = peaks(xs),
                                         peakIndex = xsg@groupidx,
                                         rtime = xsg@rt$raw,
                                         minFraction = minFr,
                                         subset = c(1, 3, 5, 7, 9, 11))
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

    rt_raw <- rtime(xod_xgr, adjusted = FALSE, bySample = TRUE)[[1]]
    rt_adj <- rtime(xod_xgr, bySample = TRUE)[[1]]

    rt_new <- xcms:::.applyRtAdjustment(rt_raw, rt_raw, rt_adj)
    expect_equal(unname(rt_new), unname(rt_adj))

    rt_new2 <- .applyRtAdjustment(rt_raw, rt_raw[200:1000], rt_adj[200:1000])

    ## Artificial examples.
    a_raw <- c(1, 2, 3, 5, 6, 7, 8, 10, 12, 13, 14, 16)
    a_adj <- a_raw + 2 # shift by 2
    b <- .applyRtAdjustment(a_raw, a_raw, a_adj)
    expect_equal(a_adj, b)
    b_2 <- .applyRtAdjustment(a_raw, a_raw[4:8], a_adj[4:8])
    expect_equal(b, b_2)

    a_adj <- a_raw - 2
    b <- .applyRtAdjustment(a_raw, a_raw, a_adj)
    expect_equal(a_adj, b)
    b_2 <- .applyRtAdjustment(a_raw, a_raw[4:8], a_adj[4:8])
    expect_equal(b, b_2)
})

test_that(".get_closest_index works", {
    expect_equal(.get_closest_index(2, c(1, 3, 5, 7)), 3)
    expect_equal(.get_closest_index(2, c(1, 3, 5, 7), method = "previous"), 1)
    expect_equal(.get_closest_index(2, c(1, 3, 5, 7), method = "closest"), 1)
    expect_equal(.get_closest_index(6, c(1, 3, 5)), 5)
    expect_equal(.get_closest_index(6, c(1, 3, 5), method = "previous"), 5)
    expect_equal(.get_closest_index(6, c(1, 3, 5), method = "closest"), 5)
    expect_equal(.get_closest_index(10, c(1, 3, 5)), 5)
    expect_equal(.get_closest_index(10, c(1, 3, 5), method = "previous"), 5)
    expect_equal(.get_closest_index(10, c(1, 3, 5), method = "closest"), 5)
    expect_equal(.get_closest_index(2, c(5, 7, 9)), 5)
    expect_equal(.get_closest_index(2, c(5, 7, 9), method = "previous"), 5)
    expect_equal(.get_closest_index(2, c(5, 7, 9), method = "closest"), 5)
    expect_equal(.get_closest_index(2, c(1, 5, 9)), 5)
    expect_equal(.get_closest_index(2, c(1, 5, 9), method = "previous"), 1)
    expect_equal(.get_closest_index(2, c(1, 5, 9), method = "closest"), 1)
    expect_equal(.get_closest_index(3, c(1, 5, 9)), 5)
    expect_equal(.get_closest_index(3, c(1, 5, 9), method = "previous"), 1)
    expect_equal(.get_closest_index(3, c(1, 5, 9), method = "closest"), 1)
    expect_equal(.get_closest_index(4, c(1, 5, 9)), 5)
    expect_equal(.get_closest_index(4, c(1, 5, 9), method = "previous"), 1)
    expect_equal(.get_closest_index(4, c(1, 5, 9), method = "closest"), 5)
    expect_equal(.get_closest_index(6, c(1, 5, 9)), 9)
    expect_equal(.get_closest_index(6, c(1, 5, 9), method = "previous"), 5)
    expect_equal(.get_closest_index(6, c(1, 5, 9), method = "closest"), 5)
    expect_equal(.get_closest_index(7, c(1, 5, 9)), 9)
    expect_equal(.get_closest_index(7, c(1, 5, 9), method = "previous"), 5)
    expect_equal(.get_closest_index(7, c(1, 5, 9), method = "closest"), 5)
    expect_equal(.get_closest_index(8, c(1, 5, 9)), 9)
    expect_equal(.get_closest_index(8, c(1, 5, 9), method = "previous"), 5)
    expect_equal(.get_closest_index(8, c(1, 5, 9), method = "closest"), 9)
})

test_that(".match_trim_vectors and index works", {
    res <- .match_trim_vectors(list(1:10, 3:10))
    expect_equal(res, list(3:10, 3:10))
    res <- .match_trim_vectors(list(3:10, 4:15))
    expect_equal(res, list(3:10, 4:11))
    res <- .match_trim_vectors(list(1:5, 1:20))
    expect_equal(res, list(1:5, 1:5))
    res <- .match_trim_vector_index(list(1:10, 3:10))
    expect_equal(res, list(3:10, 1:8))
    res <- .match_trim_vector_index(list(2:10, 2:8))
    expect_equal(res, list(1:7, 1:7))
})

test_that("adjustRtimeSubset works", {
    rt_raw <- rtime(xod_xgr, adjusted = FALSE, bySample = TRUE)
    rt_adj <- rtime(xod_xgr, adjusted = TRUE, bySample = TRUE)

    res <- adjustRtimeSubset(rt_raw, rt_adj, subset = c(1, 3),
                             method = "previous")
    expect_equal(res[[1]], rt_adj[[1]])
    expect_equal(res[[3]], rt_adj[[3]])
    expect_true(all(res[[2]] != rt_adj[[2]]))
    expect_equal(names(res[[2]]), names(rt_adj[[2]]))
    expect_equal(unname(res[[2]]), unname(rt_adj[[1]]))

    a <- res[[1]] - rt_raw[[1]]
    b <- res[[2]] - rt_raw[[2]]
    c <- res[[3]] - rt_raw[[3]]
    plot(res[[1]], a, type = "l", col = "#ff000040", lty = 2,
         ylim = range(a, b, c))
    points(res[[2]], b, type = "l", col = "#00ff0060", lty = 1)
    points(res[[3]], c, type = "l", col = "#0000ff40", lty = 2)

    res <- adjustRtimeSubset(rt_raw, rt_adj, subset = c(1, 3),
                             method = "average")
    expect_equal(res[[1]], rt_adj[[1]])
    expect_equal(res[[3]], rt_adj[[3]])
    expect_true(all(res[[2]] != rt_adj[[2]]))
    expect_true(all(res[[2]] != rt_adj[[1]]))
    expect_true(all(res[[2]] != rt_adj[[3]]))

    a <- res[[1]] - rt_raw[[1]]
    b <- res[[2]] - rt_raw[[2]]
    c <- res[[3]] - rt_raw[[3]]
    plot(res[[1]], a, type = "l", col = "#ff000040", lty = 2,
         ylim = range(a, b, c))
    points(res[[2]], b, type = "l", col = "#00ff0060", lty = 1)
    points(res[[3]], c, type = "l", col = "#0000ff40", lty = 2)
})
