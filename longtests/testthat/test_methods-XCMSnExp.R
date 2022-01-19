test_that("XCMSnExp, XCMSnExp works", {
    rts <- rtime(faahko_od)
    rts_2 <- rtime(od_x)
    expect_equal(rts, rts_2)
    ## Test with bySample.
    rts_3 <- rtime(xod_x, bySample = TRUE)
    expect_equal(rts_3, split(rts, f = fromFile(faahko_od)))
    ## Check if rtimes are correctly ordered for bySample
    rts_4 <- rtime(filterFile(faahko_od, file = 2))
    expect_equal(rts_4, rts_3[[2]])
    rts_4 <- rtime(filterFile(faahko_od, file = 3))
    expect_equal(rts_4, rts_3[[3]])
    ## Compare with the values we get from an xcmsSet:
    rtx <- faahko_xs@rt$raw
    expect_equal(unlist(rtx, use.names = FALSE),
                 unlist(rtime(faahko_xod, bySample = TRUE), use.names = FALSE))
})

test_that("featureValues,XCMSnExp works", {
    fdp <- PeakDensityParam(sampleGroups = faahko_xs$class)
    od_x <- groupChromPeaks(xod_x, param = fdp)
    xs <- group(faahko_xs, method = "density")
    fvs <- featureValues(od_x, value = "into")
    expect_equal(rownames(fvs), rownames(featureDefinitions(od_x)))
    rownames(fvs) <- NULL
    colnames(fvs) <- NULL
    gvs <- groupval(xs, value = "into")
    rownames(gvs) <- NULL
    colnames(gvs) <- NULL
    expect_equal(fvs, gvs)

    ## Use the internal function
    res <- xcms:::.feature_values(chromPeaks(od_x), featureDefinitions(od_x),
                           value = "into", method = "medret",
                           intensity = "into",
                           colnames = basename(fileNames(od_x)))
    expect_equal(featureValues(od_x, value = "into"), res)
    res <- xcms:::.feature_values(chromPeaks(od_x), featureDefinitions(od_x),
                           value = "into", method = "sum",
                           intensity = "into",
                           colnames = basename(fileNames(od_x)))
    expect_equal(featureValues(od_x, value = "into", method = "sum"), res)


    fsum <- featureSummary(xod_xg)
    fv <- featureValues(xod_xg, method = "maxint", value = "into")
    ## For feature 3 we have 2 peaks in sample 3
    idx <- unlist(featureDefinitions(xod_xg)[3, "peakidx"])
    pks <- chromPeaks(xod_xg)[idx, ]
    expect_equal(max(pks[pks[, "sample"] == 3, "into"]), fv[3, 3])
    ## For feature 37 we have 2 peaks per sample
    idx <- unlist(featureDefinitions(xod_xg)[37, "peakidx"])
    pks <- chromPeaks(xod_xg)[idx, ]
    expect_equal(max(pks[pks[, "sample"] == 1, "into"]), fv[37, 1])
    expect_equal(max(pks[pks[, "sample"] == 2, "into"]), fv[37, 2])
    expect_equal(max(pks[pks[, "sample"] == 3, "into"]), fv[37, 3])

    ## method sum
    fv <- featureValues(xod_xg, method = "sum", value = "into")
    ## For feature 3 we have 2 peaks in sample 3
    idx <- unlist(featureDefinitions(xod_xg)[3, "peakidx"])
    pks <- chromPeaks(xod_xg)[idx, ]
    expect_equal(sum(pks[pks[, "sample"] == 3, "into"]), fv[3, 3])
    ## For feature 37 we have 2 peaks per sample
    idx <- unlist(featureDefinitions(xod_xg)[37, "peakidx"])
    pks <- chromPeaks(xod_xg)[idx, ]
    expect_equal(sum(pks[pks[, "sample"] == 1, "into"]), fv[37, 1])
    expect_equal(sum(pks[pks[, "sample"] == 2, "into"]), fv[37, 2])
    expect_equal(sum(pks[pks[, "sample"] == 3, "into"]), fv[37, 3])

    ## missing
    na_num <- sum(is.na(featureValues(od_x, value = "into")))
    res <- featureValues(od_x, value = "into", missing = 123)
    expect_equal(sum(res == 123), na_num)
    res <- featureValues(od_x, value = "into", missing = "rowmin_half")
    res_na <- featureValues(od_x, value = "into")
    is_na <- is.na(rowMeans(res_na))
    for (i in which(is_na)) {
        are_na <- is.na(res_na[i, ])
        expect_true(all(res[i, are_na] == min(res_na[i, ], na.rm = TRUE) / 2))
    }
    ## Check errors
    expect_error(featureValues(od_x, value = "into", missing = "b"))
    expect_error(featureValues(od_x, value = "into", missing = TRUE))

    ## feature values with MS level > 1
    expect_error(featureValues(xod_xg, msLevel = 2), "No feature definitions")
    ## Fake feature definitions for MS level 2
    cwp <- CentWaveParam(noise = 10000, snthresh = 40,
                         prefilter = c(3, 10000))
    tmp <- xod_xg
    fd <- new("MsFeatureData")
    fd@.xData <- xcms:::.copy_env(tmp@msFeatureData)
    chromPeakData(fd)$ms_level <- 2L
    fd$featureDefinitions$ms_level <- 2L
    lockEnvironment(fd, bindings = TRUE)
    tmp@msFeatureData <- fd
    expect_true(hasChromPeaks(tmp, msLevel = 2L))
    expect_true(hasFeatures(tmp, msLevel = 2L))
    expect_equal(featureValues(tmp, msLevel = 2L), featureValues(xod_xg))

    tmp <- findChromPeaks(tmp, add = TRUE, param = cwp)
    expect_equal(unname(chromPeaks(tmp, msLevel = 1L)[, "into"]),
                 unname(chromPeaks(tmp, msLevel = 2L)[, "into"]))
    ## correspondence
    pdp <- PeakDensityParam(sampleGroups = rep(1, 3))
    tmp <- groupChromPeaks(tmp, param = pdp, msLevel = 1L)
    tmp <- groupChromPeaks(tmp, param = pdp, msLevel = 2L, add = TRUE)
    expect_true(hasFeatures(tmp, msLevel = 1L))
    expect_true(hasFeatures(tmp, msLevel = 2L))

    all <- featureValues(tmp)
    ms1 <- featureValues(tmp, msLevel = 1L)
    ms2 <- featureValues(tmp, msLevel = 2L)
    expect_equal(all, rbind(ms1, ms2))
    rownames(ms1) <- rownames(ms2) <- NULL
    expect_equal(ms1, ms2)
})

test_that("adjustRtime,peakGroups works", {
    xod <- faahko_xod
    xs <- faahko_xs
    ## Group these
    xsg <- group(xs)
    xodg <- groupChromPeaks(xod,
                            param = PeakDensityParam(sampleGroups = xs$class))
    pks <- chromPeaks(xodg)
    rownames(pks) <- NULL
    expect_equal(peaks(xsg), pks[, colnames(peaks(xsg))])
    expect_equal(xsg@groupidx, featureDefinitions(xodg)$peakidx)
    expect_true(length(processHistory(xodg,
                                      type = .PROCSTEP.PEAK.DETECTION)) == 1)
    expect_true(length(processHistory(xodg,
                                      type = .PROCSTEP.PEAK.GROUPING)) == 1)
    ## Now do the retention time correction
    xsr <- retcor(xsg, method = "peakgroups", missing = 0, span = 0.3)
    ## minFr <- (length(fileNames(xod)) - 1) / length(fileNames(xod))
    p <- PeakGroupsParam(minFraction = 1, span = 0.3)
    xodr <- adjustRtime(xodg, param = p)
    ## Check that we've got process histories.
    expect_true(validObject(xodr))
    expect_true(hasChromPeaks(xodr))
    expect_true(!hasFeatures(xodr))
    ## But we would like to keep the related process history step:
    expect_true(hasAdjustedRtime(xodr))
    expect_true(hasFeatures(xodg))
    ## We want to keep the process history step of the feature alignment!
    expect_true(length(processHistory(xodr,
                                      type = .PROCSTEP.PEAK.GROUPING)) == 1)
    expect_true(length(processHistory(xodr,
                                      type = .PROCSTEP.RTIME.CORRECTION)) == 1)
    ## Different from original:
    expect_true(sum(chromPeaks(xod)[, "rt"] != chromPeaks(xodr)[, "rt"]) > 200)
    expect_true(sum(chromPeaks(xod)[, "rtmin"] != chromPeaks(xodr)[, "rtmin"]) > 200)
    expect_true(sum(chromPeaks(xod)[, "rtmax"] != chromPeaks(xodr)[, "rtmax"]) > 200)
    ## between xcmsSet and XCMSnExp
    pks <- chromPeaks(xodr)
    rownames(pks) <- NULL
    expect_equal(pks[, colnames(peaks(xsr))], peaks(xsr))
    ## To compare the adjusted retention time we have to extract it by sample!
    ## Otherwise the ordering will not be the same, as rtime is ordered by
    ## retention time, but @rt$raw by sample.
    expect_equal(unlist(adjustedRtime(xodr, bySample = TRUE), use.names = FALSE),
                 unlist(xsr@rt$corrected, use.names = FALSE))
    ## Just to ensure - are the raw rt the same?
    expect_equal(unlist(rtime(xod, bySample = TRUE), use.names = FALSE),
                 unlist(xs@rt$raw, use.names = FALSE))
    ## Check that we get the same by supplying the peakGroupsMatrix.
    pgm <- adjustRtimePeakGroups(xodg, param = p)
    p_2 <- p
    minFraction(p_2) <- 0.5
    extraPeaks(p_2) <- 20
    peakGroupsMatrix(p_2) <- pgm
    xodr_2 <- adjustRtime(xodg, param = p_2)
    expect_equal(adjustedRtime(xodr), adjustedRtime(xodr_2))
    expect_equal(chromPeaks(xodr), chromPeaks(xodr_2))
    p_got <- processParam(
        processHistory(xodr, type = .PROCSTEP.RTIME.CORRECTION)[[1]])
    peakGroupsMatrix(p_got) <- matrix(ncol = 0, nrow = 0)
    expect_equal(p_got, p)
    expect_equal(processParam(
        processHistory(xodr_2, type = .PROCSTEP.RTIME.CORRECTION)[[1]]),
        p_2)
    ## Doing an additional grouping
    xodrg <- groupChromPeaks(xodr, param = PeakDensityParam(sampleGroups =
                                                                xs$class))
    expect_true(length(processHistory(xodrg,
                                      type = .PROCSTEP.PEAK.GROUPING)) == 2)
    expect_true(hasAdjustedRtime(xodrg))
    expect_true(hasFeatures(xodrg))
    xsrg <- group(xsr)
    expect_equal(xsrg@groupidx, featureDefinitions(xodrg)$peakidx)

    ## Mod settings:
    xsr <- retcor(xsg, method = "peakgroups", missing = 0, span = 1)
    xodr <- adjustRtime(xodg, param = PeakGroupsParam(minFraction = 1,
                                                      span = 1))
    pks <- chromPeaks(xodr)
    rownames(pks) <- NULL
    expect_equal(pks[, colnames(peaks(xsr))], peaks(xsr))
    expect_equal(unlist(adjustedRtime(xodr, bySample = TRUE), use.names = FALSE),
                 unlist(xsr@rt$corrected, use.names = FALSE))

    xsr <- retcor(xsg, method = "peakgroups", missing = 0, span = 1,
                  smooth = "linear")
    xodr <- adjustRtime(xodg, param = PeakGroupsParam(minFraction = 1,
                                                      span = 1,
                                                      smooth = "linear"))
    pks <- chromPeaks(xodr)
    rownames(pks) <- NULL
    expect_equal(pks[, colnames(peaks(xsr))], peaks(xsr))
    expect_equal(unlist(adjustedRtime(xodr, bySample = TRUE), use.names = FALSE),
                 unlist(xsr@rt$corrected, use.names = FALSE))

    xsr <- retcor(xsg, method = "peakgroups", missing = 0, span = 1,
                  family = "symmetric")
    xodr <- adjustRtime(xodg, param = PeakGroupsParam(minFraction = 1,
                                                      span = 1,
                                                      family = "symmetric"))
    pks <- chromPeaks(xodr)
    rownames(pks) <- NULL
    expect_equal(pks[, colnames(peaks(xsr))], peaks(xsr))
    expect_equal(unlist(adjustedRtime(xodr, bySample = TRUE), use.names = FALSE),
                 unlist(xsr@rt$corrected, use.names = FALSE))
    ## Dropping results.
    tmp <- dropAdjustedRtime(xodr)
    expect_equal(tmp, xod)

    ## With subset.
    res_sub <- adjustRtime(
        xodg, param = PeakGroupsParam(subset = c(1, 3),
                                      subsetAdjust = "previous"))
    expect_true(all(rtime(res_sub, bySample = TRUE)[[1]] !=
                    rtime(xodg, bySample = TRUE)[[1]]))
    expect_true(all(rtime(res_sub, bySample = TRUE)[[2]] !=
                    rtime(xodg, bySample = TRUE)[[2]]))
    expect_true(all(rtime(res_sub, bySample = TRUE)[[3]] !=
                    rtime(xodg, bySample = TRUE)[[3]]))
    expect_equal(unname(rtime(res_sub, bySample = TRUE)[[1]]),
                 unname(rtime(res_sub, bySample = TRUE)[[2]]))
    expect_equal(rtime(res_sub, bySample = TRUE)[[2]],
                 xcms:::.applyRtAdjustment(rtime(xodg, bySample = TRUE)[[2]],
                                           rtime(xodg, bySample = TRUE)[[1]],
                                           rtime(res_sub, bySample = TRUE)[[1]]))
    res_sub <- adjustRtime(
        xodg, param = PeakGroupsParam(subset = c(1, 3),
                                      subsetAdjust = "average"))
    expect_true(all(rtime(res_sub, bySample = TRUE)[[1]] !=
                    rtime(xodg, bySample = TRUE)[[1]]))
    expect_true(all(rtime(res_sub, bySample = TRUE)[[2]] !=
                    rtime(xodg, bySample = TRUE)[[2]]))
    expect_true(all(rtime(res_sub, bySample = TRUE)[[3]] !=
                    rtime(xodg, bySample = TRUE)[[3]]))
    expect_true(all(rtime(res_sub, bySample = TRUE)[[1]] !=
                    rtime(res_sub, bySample = TRUE)[[2]]))
    tmp <- adjustRtime(xodg, param = PeakGroupsParam())

    ## With subsetAdjust = "average" and the left-out being at the end.
    res_sub <- adjustRtime(
        xodg, param = PeakGroupsParam(subset = 1:2, subsetAdjust = "average"))
    res_2 <- adjustRtime(
        xodg, param = PeakGroupsParam(subset = 1:2, subsetAdjust = "previous"))
    expect_equal(rtime(res_sub), rtime(res_2))
})

test_that("featureValues,XCMSnExp works as with groupval", {
    od_x <- faahko_xod
    xs <- faahko_xs

    p <- PeakDensityParam(sampleGroups = xs$class)
    od_x <- groupChromPeaks(od_x, param = p)

    xs <- group(xs, method = "density")

    expect_equal(unname(groupval(xs, value = "into")),
                 unname(featureValues(od_x, value = "into")))
    expect_equal(unname(groupval(xs, method = "maxint", value = "into")),
                 unname(featureValues(od_x, method = "maxint", value = "into")))
    ## Checking errors
    expect_error(featureValues(od_x, value = "bla"))
})

test_that("groupChromPeaks,XCMSnExp,PeakDensityParam works", {
    od_x <- faahko_xod
    xs <- faahko_xs
    ## Check error if no features were found. issue #273
    pdp <- PeakDensityParam(sampleGroups = xs$class, minSamples = 30)
    expect_warning(groupChromPeaks(od_x, param = pdp), "Unable to group any chromatographic peaks.")

    fdp <- PeakDensityParam(sampleGroups = xs$class)
    od_x <- groupChromPeaks(od_x, param = fdp)
    xs <- group(xs, method = "density")
    expect_equal(xs@groupidx, featureDefinitions(od_x)$peakidx)
    fg <- featureDefinitions(od_x)
    fg <- S4Vectors::as.matrix(fg[, !(colnames(fg) %in% c("peakidx", "ms_level"))])
    rownames(fg) <- NULL
    expect_equal(xs@groups, fg)
    expect_true(length(processHistory(od_x)) == 2)
    ph <- processHistory(od_x, type = .PROCSTEP.PEAK.GROUPING)[[1]]
    expect_equal(processParam(ph), fdp)
    expect_equal(rownames(featureDefinitions(od_x)),
                 .featureIDs(nrow(featureDefinitions(od_x))))

    fdp2 <- PeakDensityParam(sampleGroups = xs$class, binSize = 2,
                             minFraction = 0.8)
    od_x <- groupChromPeaks(od_x, param = fdp2)
    xs <- group(xs, method = "density", minfrac = 0.8, mzwid = 2)
    expect_equal(xs@groupidx, featureDefinitions(od_x)$peakidx)
    fg <- featureDefinitions(od_x)
    fg <- S4Vectors::as.matrix(fg[, !(colnames(fg) %in% c("peakidx", "ms_level"))])
    rownames(fg) <- NULL
    expect_equal(xs@groups, fg)
    expect_true(length(processHistory(od_x)) == 2)
    ph <- processHistory(od_x, type = .PROCSTEP.PEAK.GROUPING)[[1]]
    expect_equal(processParam(ph), fdp2)
    expect_equal(rownames(featureDefinitions(od_x)),
                 .featureIDs(nrow(featureDefinitions(od_x))))

    pdp <- PeakDensityParam(sampleGroups = xs$class)
    res <- groupChromPeaks(od_x, param = pdp)
    res_2 <- groupChromPeaks(res, param = pdp)
    expect_equal(featureDefinitions(res), featureDefinitions(res_2))
    res_2 <- groupChromPeaks(res, param = pdp, add = TRUE)
    expect_true(nrow(featureDefinitions(res_2)) ==
                2 * nrow(featureDefinitions(res)))
    nr <- nrow(featureDefinitions(res))
    expect_equal(featureDefinitions(res),
                 featureDefinitions(res_2)[1:nr, ])
    expect_equal(featureDefinitions(res)$mzmed,
                 featureDefinitions(res_2)$mzmed[(nr + 1):(2 * nr)])
    expect_equal(featureDefinitions(res)$peakidx,
                 featureDefinitions(res_2)$peakidx[(nr + 1):(2 * nr)])

    expect_error(groupChromPeaks(od_x, param = pdp, msLevel = 2), "MS level 2")
    expect_error(groupChromPeaks(od_x, param = pdp, msLevel = 1:4),
                 "one MS level at a time")
})

test_that("groupPeaks,XCMSnExp,MzClustParam works", {
    p <- MzClustParam(sampleGroups = sampclass(fticr_xs))
    fticr_xod2 <- groupChromPeaks(fticr_xod, param = p)
    fticr_xs2 <- group(fticr_xs, method = "mzClust")
    expect_equal(fticr_xs2@groupidx, featureDefinitions(fticr_xod2)$peakidx)
    fg <- featureDefinitions(fticr_xod2)
    fg <- S4Vectors::as.matrix(fg[, -ncol(fg)])
    rownames(fg) <- NULL
    expect_equal(fticr_xs2@groups, fg)
    expect_true(length(processHistory(fticr_xod2)) == 2)
    ph <- processHistory(fticr_xod2,
                         type = .PROCSTEP.PEAK.GROUPING)[[1]]
    expect_equal(processParam(ph), p)
    expect_equal(rownames(featureDefinitions(fticr_xod2)),
                 .featureIDs(nrow(featureDefinitions(fticr_xod2))))
    p2 <- MzClustParam(sampleGroups = fticr_xs$class, absMz = 1,
                       minFraction = 0.8)
    fticr_xod2 <- groupChromPeaks(fticr_xod, param = p2)
    fticr_xs2 <- group(fticr_xs, method = "mzClust", minfrac = 0.8, mzabs = 1)
    expect_equal(fticr_xs2@groupidx, featureDefinitions(fticr_xod2)$peakidx)
    fg <- featureDefinitions(fticr_xod2)
    fg <- S4Vectors::as.matrix(fg[, -ncol(fg)])
    rownames(fg) <- NULL
    expect_equal(fticr_xs2@groups, fg)
    expect_true(length(processHistory(fticr_xod2)) == 2)
    ph <- processHistory(fticr_xod2,
                         type = .PROCSTEP.PEAK.GROUPING)[[1]]
    expect_equal(processParam(ph), p2)
    expect_equal(rownames(featureDefinitions(fticr_xod2)),
                 .featureIDs(nrow(featureDefinitions(fticr_xod2))))
})

test_that("groupChromPeaks,XCMSnExp,NearestPeaksParam works", {
    od_x <- faahko_xod
    xs <- faahko_xs
    p <- NearestPeaksParam(sampleGroups = xs$class)
    od_x <- groupChromPeaks(od_x, param = p)
    xs <- group(xs, method = "nearest")
    expect_equal(xs@groupidx, featureDefinitions(od_x)$peakidx)
    fg <- featureDefinitions(od_x)
    fg <- S4Vectors::as.matrix(fg[, !(colnames(fg) %in% c("peakidx", "ms_level"))])
    rownames(fg) <- NULL
    expect_equal(xs@groups, fg)
    expect_true(length(processHistory(od_x)) == 2)
    ph <- processHistory(od_x, type = .PROCSTEP.PEAK.GROUPING)[[1]]
    expect_equal(processParam(ph), p)
    expect_equal(rownames(featureDefinitions(od_x)),
                 .featureIDs(nrow(featureDefinitions(od_x))))
    fdp2 <- NearestPeaksParam(sampleGroups = xs$class, kNN = 3)
    od_x <- groupChromPeaks(od_x, param = fdp2)
    xs <- group(xs, method = "nearest", kNN = 3)
    expect_equal(xs@groupidx, featureDefinitions(od_x)$peakidx)
    fg <- featureDefinitions(od_x)
    fg <- S4Vectors::as.matrix(fg[, !(colnames(fg) %in% c("peakidx", "ms_level"))])
    rownames(fg) <- NULL
    expect_equal(xs@groups, fg)
    expect_true(length(processHistory(od_x)) == 2)
    ph <- processHistory(od_x, type = .PROCSTEP.PEAK.GROUPING)[[1]]
    expect_equal(processParam(ph), fdp2)
    expect_equal(rownames(featureDefinitions(od_x)),
                 .featureIDs(nrow(featureDefinitions(od_x))))

    expect_error(groupChromPeaks(od_x, param = p, msLevel = 2), "MS level 2")
    expect_error(groupChromPeaks(od_x, param = p, msLevel = 1:3), " at a time")
    res <- groupChromPeaks(od_x, param = p)
    res_2 <- groupChromPeaks(res, param = p)
    expect_equal(featureDefinitions(res), featureDefinitions(res_2))
    res_2 <- groupChromPeaks(res, param = p, add = TRUE)
    expect_true(nrow(featureDefinitions(res_2)) ==
                2 * nrow(featureDefinitions(res)))
    nr <- nrow(featureDefinitions(res))
    expect_equal(featureDefinitions(res), featureDefinitions(res)[1:nr, ])
    expect_equal(featureDefinitions(res)$peakidx,
                 featureDefinitions(res_2)$peakidx[(nr + 1):(2 * nr)])
    expect_equal(featureDefinitions(res)$mzmed,
                 featureDefinitions(res_2)$mzmed[(nr + 1):(2 * nr)])
})

test_that("fillChromPeaks,XCMSnExp with MSW works", {
    p <- MzClustParam()
    fticr_xodg <- groupChromPeaks(fticr_xod, param = p)
    expect_error(res <- fillChromPeaks(fticr_xod))
    res <- fillChromPeaks(fticr_xodg)

    ## Got a signal for all of em.
    expect_true(!any(is.na(featureValues(res))))
    ## 1) Compare with what I get for xcmsSet.
    tmp_x <- fticr_xs
    tmp_x <- group(tmp_x, method = "mzClust")
    tmp_x <- fillPeaks(tmp_x, method = "MSW")
    ## Compare
    expect_equal(unname(groupval(tmp_x)),
                 unname(featureValues(res, value = "index")))
    expect_equal(unname(groupval(tmp_x, value = "maxo")),
                 unname(featureValues(res, value = "maxo")))
    expect_equal(unname(groupval(tmp_x, value = "into")),
                 unname(featureValues(res, value = "into")))
    expect_equal(unname(groupval(tmp_x, value = "mz")),
                 unname(featureValues(res, value = "mz")))
    expect_equal(unname(groupval(tmp_x, value = "mzmin")),
                 unname(featureValues(res, value = "mzmin")))
    expect_equal(unname(groupval(tmp_x, value = "mzmax")),
                 unname(featureValues(res, value = "mzmax")))
    ## OK
    ## 2) Check if the fillChromPeaks returns same/similar data than the
    ##    findChromPeaks does:
    fdef <- featureDefinitions(fticr_xodg)
    pkArea <- do.call(
        rbind,
        lapply(
            fdef$peakidx, function(z) {
                tmp <- chromPeaks(fticr_xodg)[z, c("rtmin", "rtmax",
                                                   "mzmin", "mzmax"),
                                              drop = FALSE]
                pa <- c(median(tmp[, 1]), median(tmp[, 2]),
                        median(tmp[, 3]), median(tmp[, 4]))
                return(pa)
            }
        ))
    colnames(pkArea) <- c("rtmin", "rtmax", "mzmin", "mzmax")
    pkArea <- cbind(group_idx = 1:nrow(pkArea), pkArea,
                    mzmed = fdef$mzmed)
    ## Get peak data for all peaks in the first file
    allPks <- .getMSWPeakData(filterFile(fticr_xodg, file = 1),
                              peakArea = pkArea,
                              sample_idx = 1,
                              cn = colnames(chromPeaks(fticr_xodg)))
    curP <- chromPeaks(res)[chromPeaks(res)[, "sample"] == 1, ]
    curP <- curP[order(curP[, "mz"]), ]
    expect_equal(unname(allPks[, "mz"]), unname(curP[, "mz"]))
    expect_equal(unname(allPks[, "maxo"]), unname(curP[, "maxo"]))
    expect_true(cor(allPks[, "into"], curP[, "into"]) > 0.99) ## Not exactly the
    ## same but highly similar.
})
