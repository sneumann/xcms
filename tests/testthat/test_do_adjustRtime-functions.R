test_that("getPeakGroupsRtMatrix works", {
    skip_on_os(os = "windows", arch = "i386")

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
    skip_on_os(os = "windows", arch = "i386")

    xsa <- xod_xgr
    xsg <- xod_xg
    misSamp <- 1
    minFr <- (length(fileNames(xsa)) - misSamp) / length(fileNames(xsa))
    expect_error(do_adjustRtime_peakGroups(
        peaks = chromPeaks(xsg), peakIndex = featureDefinitions(xsg)$peakidx,
        rtime = rtime(xsg, bySample = TRUE), minFraction = minFr, subset = "4"),
        "expected to be an integer")
    expect_error(do_adjustRtime_peakGroups(
        peaks = chromPeaks(xsg), peakIndex = featureDefinitions(xsg)$peakidx,
        rtime = rtime(xsg, bySample = TRUE), minFraction = minFr,
        subset = 4L), "out of range")
    expect_error(do_adjustRtime_peakGroups(
        peaks = chromPeaks(xsg), peakIndex = featureDefinitions(xsg)$peakidx,
        rtime = rtime(xsg, bySample = TRUE), minFraction = minFr,
        subset = c(1, 2, 5, 14)), "out of range")

    ## Use only a subset.
    res_sub <- do_adjustRtime_peakGroups(
        peaks = chromPeaks(xsg), peakIndex = featureDefinitions(xsg)$peakidx,
        rtime = rtime(xsg, bySample = TRUE), minFraction = minFr,
        subset = c(1, 3))
})

test_that("applyRtAdjustment works", {
    skip_on_os(os = "windows", arch = "i386")

    xs <- faahko
    ## group em.
    ## xsg <- group(xs)
    ## ## align em.
    ## xsa <- retcor(xsg, method = "peakgroups")
    pksAdj <- .applyRtAdjToChromPeaks(chromPeaks(xod_xg),
                                      rtraw = rtime(xod_xg, bySample = TRUE),
                                      rtadj = rtime(xod_xgr, bySample = TRUE))
    expect_equal(pksAdj, chromPeaks(xod_xgr))
    ## Reset em.
    pksRaw <- .applyRtAdjToChromPeaks(pksAdj,
                                      rtraw = rtime(xod_xgr, bySample = TRUE),
                                      rtadj = rtime(xod_xg, bySample = TRUE))
    expect_equal(pksRaw, chromPeaks(xod_xg))

    rt_raw <- rtime(xod_xgr, adjusted = FALSE, bySample = TRUE)[[1]]
    rt_adj <- rtime(xod_xgr, bySample = TRUE)[[1]]

    rt_new <- .applyRtAdjustment(rt_raw, rt_raw, rt_adj)
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
    skip_on_os(os = "windows", arch = "i386")

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
    skip_on_os(os = "windows", arch = "i386")

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
    skip_on_os(os = "windows", arch = "i386")

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

test_that(".adjustRtime_peakGroupsMatrix works", {
    ## Expect the same results by running just this function with a pre-defined
    ## peakGroupsMatrix.
    ph <- processHistory(xod_xgr)
    pgm <- peakGroupsMatrix(ph[[length(ph)]]@param)
    rts <- split(rtime(xod_xg), fromFile(xod_xg))

    ## Errors
    expect_error(.adjustRtime_peakGroupsMatrix(rts), "peakGroupsMatrix")
    expect_error(.adjustRtime_peakGroupsMatrix(peakGroupsMatrix = pgm), "rtime")
    expect_error(.adjustRtime_peakGroupsMatrix(rts[[1L]], pgm), "list of")
    expect_error(.adjustRtime_peakGroupsMatrix(list(1:3, "a"), pgm), "numeric")
    expect_error(.adjustRtime_peakGroupsMatrix(rts, pgm, subset = "a"),
                 "integer")
    expect_error(.adjustRtime_peakGroupsMatrix(rts, pgm, subset = 1), "is 2")
    expect_error(.adjustRtime_peakGroupsMatrix(rts, pgm, subset = c(1, 10)),
                 "out of range")
    expect_error(.adjustRtime_peakGroupsMatrix(rts[1:2], pgm), "columns of")

    ## Results are identical
    res <- .adjustRtime_peakGroupsMatrix(rts, pgm, span = 0.4, smooth = "loess")
    expect_equal(res, split(rtime(xod_xgr), fromFile(xod_xgr)))

    ## Works with subset.
    tmp <- adjustRtime(xod_xg, PeakGroupsParam(subset = c(1, 3), span = 0.4))
    ph <- processHistory(tmp)
    pgm <- peakGroupsMatrix(ph[[length(ph)]]@param)

    res <- .adjustRtime_peakGroupsMatrix(rts, pgm, span = 0.4,
                                         smooth = "loess", subset = c(1, 3))
    expect_equal(res, split(rtime(tmp), fromFile(tmp)))

})

test_that(".define_peak_groups works", {
    expect_error(.define_peak_groups(), "peaks")
    expect_error(.define_peak_groups(
        chromPeaks(xod_xgrg)[1:10, ], featureDefinitions(xod_xgrg)$peakidx),
        "outside")
    expect_error(.define_peak_groups(
        chromPeaks(xod_xgrg), featureDefinitions(xod_xgrg)$peakidx,
        subset = "a"), "integer")
    expect_error(.define_peak_groups(
        chromPeaks(xod_xgrg), featureDefinitions(xod_xgrg)$peakidx,
        subset = c(1, 10), total_samples = 3), "out of range")
    expect_error(.define_peak_groups(
        chromPeaks(xod_xgrg), featureDefinitions(xod_xgrg)$peakidx,
        subset = c(1), total_samples = 3), "is 2")
    res <- .define_peak_groups(chromPeaks(
        xod_xg), featureDefinitions(xod_xg)$peakidx, subset = integer(),
        total_samples = 3)
    expect_true(is.matrix(res))
    expect_true(ncol(res) == 3)
    expect_true(is.numeric(res))

    ## Compare with results from full run.
    ph <- processHistory(xod_xgr)
    pgm <- peakGroupsMatrix(ph[[length(ph)]]@param)
    expect_equal(unname(res), unname(pgm))

    ## Run with subset.
    res <- .define_peak_groups(chromPeaks(
        xod_xg), featureDefinitions(xod_xg)$peakidx, subset = c(1, 3),
        total_samples = 3, minFraction = 0.9, extraPeaks = 1)
    expect_true(is.matrix(res))
    expect_true(ncol(res) == 2)
    expect_true(is.numeric(res))
    expect_true(!any(is.na(res)))

    pgm <- adjustRtimePeakGroups(
        xod_xg, param = PeakGroupsParam(span = 0.4, subset = c(1, 3),
                                        minFraction = 0.9, extraPeaks = 1))
    expect_equal(unname(res), unname(pgm))

    res <- xcms:::.define_peak_groups(chromPeaks(
        xod_xg), featureDefinitions(xod_xg)$peakidx, subset = c(1, 3),
        total_samples = 3, minFraction = 0.5, extraPeaks = 1)
    expect_true(is.matrix(res))
    expect_true(ncol(res) == 2)
    expect_true(is.numeric(res))
    expect_true(any(is.na(res)))
})

test_that(".getPeakGroupsRtMatrix works with subsets", {
    res <- .getPeakGroupsRtMatrix(chromPeaks(xod_xg),
                                  featureDefinitions(xod_xg)$peakidx,
                                  sampleIndex = c(1, 3),
                                  missingSample = 0,
                                  extraPeaks = 1L)
    expect_true(!any(is.na(res)))
})
