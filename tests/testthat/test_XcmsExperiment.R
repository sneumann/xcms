library(MsExperiment)
fls <- normalizePath(faahko_3_files)
df <- data.frame(mzML_file = basename(fls),
                 dataOrigin = fls,
                 sample = c("ko15", "ko16", "ko18"))
mse <- readMsExperiment(spectraFiles = fls, sampleData = df)
p <- CentWaveParam(noise = 10000, snthresh = 40, prefilter = c(3, 10000))
xmse <- findChromPeaks(mse, param = p)
pdp <- PeakDensityParam(sampleGroups = rep(1, 3))
xmseg <- groupChromPeaks(xmse, param = pdp, add = FALSE)

fl <- system.file("TripleTOF-SWATH", "PestMix1_SWATH.mzML", package = "msdata")
mse_dia <- readMsExperiment(fl)

test_that(".empty_chrom_peaks works", {
    res <- .empty_chrom_peaks()
    expect_true(nrow(res) == 0)
    expect_equal(colnames(res), c(.REQ_PEAKS_COLS, "maxo"))

    res <- .empty_chrom_peaks(sample = FALSE)
    expect_true(nrow(res) == 0)
    expect_true(!any(colnames(res) == "sample"))
})

test_that("XcmsExperiment validation works", {
    a <- new("XcmsExperiment")
    expect_true(validObject(a))

    a@chromPeakData <- data.frame(a = 1:3, b = 1:3)
    expect_error(validObject(a), "chromPeakData")

    a <- new("XcmsExperiment")
    a@chromPeaks <- cbind(a = 1, b = 2)
    expect_error(validObject(a), "chromPeaks")
})

## All chrom peak related functions:
## - findChromPeaks
## - hasChromPeaks
## - dropChromPeaks
## - chromPeaks
## - chromPeakData
test_that("findChromPeaks,MsExperiment et al works", {
    expect_error(findChromPeaks(MsExperiment(), param = p), "No spectra")

    a <- MsExperiment()
    spectra(a) <- spectra(mse)
    expect_error(findChromPeaks(a, param = p), "No link")

    res <- xmse
    expect_equal(res@chromPeaks, chromPeaks(faahko_xod))
    expect_equal(res@chromPeakData, as.data.frame(chromPeakData(faahko_xod)))
    expect_true(hasChromPeaks(res))

    ## chromPeaks
    expect_equal(chromPeaks(res), res@chromPeaks)
    cp <- chromPeaks(res, isFilledColumn = TRUE)
    expect_true(any(colnames(cp) == "is_filled"))
    cp <- chromPeaks(res, msLevel = 2)
    expect_true(is.matrix(cp))
    expect_true(nrow(cp) == 0)
    cp <- chromPeaks(res, rt = c(3000, 3500), type = "within")
    expect_true(all(cp[, "rt"] >= 3000 & cp[, "rt"] <= 3500))
    cp <- chromPeaks(res, mz = c(300, 310), type = "within")
    expect_true(all(cp[, "mz"] >= 300 & cp[, "mz"] <= 310))

    ## chromPeakData
    expect_equal(chromPeakData(res), DataFrame(res@chromPeakData))
    expect_true(is.data.frame(chromPeakData(xmse, return.type = "data.frame")))
    expect_s4_class(chromPeakData(xmse), "DataFrame")
    expect_true(nrow(chromPeakData(xmse, 2:3)) == 0)
    expect_true(is.integer(chromPeakData(res)$ms_level))

    ## dropChromPeaks
    rres <- dropChromPeaks(res)
    expect_true(length(rres@processHistory) == 0)
    expect_true(nrow(rres@chromPeaks) == 0)
    expect_false(hasChromPeaks(rres))

    res2 <- findChromPeaks(mse, param = p, msLevel = 2L)
    expect_true(nrow(res2@chromPeaks) == 0)
    expect_false(hasChromPeaks(res2))

    res2 <- findChromPeaks(res, param = p, msLevel = 2L, add = TRUE)
    expect_equal(res@chromPeaks, res2@chromPeaks)
    expect_equal(res@chromPeakData, res2@chromPeakData)
    expect_true(length(res2@processHistory) == 2)
    expect_true(is.integer(chromPeakData(res2)$ms_level))

    res2 <- findChromPeaks(res, param = p, msLevel = 2L, add = FALSE)
    expect_equal(nrow(res2@chromPeaks), 0)
    expect_true(length(res2@processHistory) == 1)
    expect_true(is.integer(chromPeakData(res2)$ms_level))

    res2 <- findChromPeaks(mse, param = p, chunkSize = -1)
    expect_equal(res@chromPeaks, res2@chromPeaks)
    expect_true(is.integer(chromPeakData(res2)$ms_level))

    expect_true(hasChromPeaks(res))
    expect_true(hasChromPeaks(res, msLevel = 1L))
    expect_true(hasChromPeaks(res, msLevel = 1:4))
    expect_false(hasChromPeaks(res, msLevel = 2))
})

## That's from XcmsExperiment-functions.R
test_that("subsetting,XcmsExperiment works", {
    expect_error(.subset_xcms_experiment(xmse, i = 1:4), "out of bounds")
    expect_error(.subset_xcms_experiment(xmse, i = c(1, 1, 2)), "Duplicated")

    res <- .subset_xcms_experiment(xmse, i = 2)
    expect_true(hasChromPeaks(res))
    expect_true(all(chromPeaks(res)[, "sample"] == 1L))
    expect_true(length(res) == 1)
    expect_equal(spectra(res), spectra(mse[2L]))
    cp <- chromPeaks(xmse)
    expect_equal(chromPeaks(res)[, colnames(chromPeaks(res)) != "sample"],
                 cp[cp[, "sample"] == 2L, colnames(cp) != "sample"])

    res <- .subset_xcms_experiment(xmse, i = c(3, 1))
    expect_true(hasChromPeaks(res))
    cp3 <- chromPeaks(res)[chromPeaks(res)[, "sample"] == 1, ]
    cp1 <- chromPeaks(res)[chromPeaks(res)[, "sample"] == 2, ]
    expect_equal(cp3[, colnames(cp3) != "sample"],
                 cp[cp[, "sample"] == 3, colnames(cp) != "sample"])
    expect_equal(cp1[, colnames(cp1) != "sample"],
                 cp[cp[, "sample"] == 1, colnames(cp) != "sample"])

    res <- .subset_xcms_experiment(xmse, i = 3, keepChromPeaks = FALSE)
    expect_false(hasChromPeaks(res))
    expect_true(length(res@processHistory) == 0)

    res <- .subset_xcms_experiment(xmse, i = 3, keepChromPeaks = FALSE,
                                   ignoreHistory = TRUE)
    expect_false(hasChromPeaks(res))
    expect_true(length(res@processHistory) == 1)

    expect_error(xmse[3, 4], "not supported")
    res <- xmse[c(3, 1)]
    expect_true(hasChromPeaks(res))
    cp3 <- chromPeaks(res)[chromPeaks(res)[, "sample"] == 1, ]
    cp1 <- chromPeaks(res)[chromPeaks(res)[, "sample"] == 2, ]
    expect_equal(cp3[, colnames(cp3) != "sample"],
                 cp[cp[, "sample"] == 3, colnames(cp) != "sample"])
    expect_equal(cp1[, colnames(cp1) != "sample"],
                 cp[cp[, "sample"] == 1, colnames(cp) != "sample"])

    ## peak grouping results.
    res <- .subset_xcms_experiment(xmseg, i = c(3, 1))
    expect_true(hasChromPeaks(res))
    expect_true(all(chromPeaks(res)[, "sample"] %in% 1:2))
    expect_false(hasFeatures(res))
    res <- .subset_xcms_experiment(xmseg, i = c(3, 1), keepFeatures = TRUE)
    expect_true(hasChromPeaks(res))
    expect_true(hasFeatures(res))
    expect_equal(featureValues(res), featureValues(xmseg)[, c(3, 1)])
})

test_that("filterRt,XcmsExperiment works", {
    res <- filterRt(xmse)
    expect_equal(res, xmse)
    expect_s4_class(res, "XcmsExperiment")

    res <- filterRt(xmse, rt = c(3000, 3500))
    expect_s4_class(res, "XcmsExperiment")
    expect_true(all(rtime(spectra(res)) >= 3000 &
                    rtime(spectra(res)) <= 3500))
    expect_true(all(chromPeaks(res)[, "rt"] >= 3000 &
                    chromPeaks(res)[, "rt"] <= 3500))
    ## Check error: define msLevel
    expect_warning(res_2 <- filterRt(xmse, rt = c(3000, 3500), msLevel = 2L),
                   "ignored")
    expect_equal(chromPeaks(res), chromPeaks(res_2))
    expect_equal(rtime(spectra(res)), rtime(spectra(res_2)))

    res <- filterRt(xmseg, rt = c(3000, 3500))
    expect_true(all(rtime(spectra(res)) >= 3000 &
                    rtime(spectra(res)) <= 3500))
    expect_true(all(chromPeaks(res)[, "rt"] >= 3000 &
                    chromPeaks(res)[, "rt"] <= 3500))
    expect_true(hasFeatures(res))
    expect_true(nrow(featureDefinitions(res)) < nrow(featureDefinitions(xmseg)))
    expect_true(validObject(res))
    fv <- featureValues(res)
    expect_equal(fv[, 1L], featureValues(xmseg)[rownames(fv), 1L])
    ## no match for other samples because some chrom peaks are out of rt range
    expect_true(all(unlist(featureDefinitions(res)$peakidx)) %in%
                seq_len(nrow(chromPeaks(res))))
})

test_that("filterFile,XcmsExperiment works", {
    res <- filterFile(xmse)
    expect_s4_class(res, "XcmsExperiment")
    expect_equal(res, xmse)

    res <- filterFile(xmse, 2)
    expect_equal(res, xmse[2])
    expect_true(hasChromPeaks(res))

    res <- filterFile(xmse, c(3, 1))
    expect_equal(res, xmse[c(1, 3)])

    res <- filterFile(xmseg, c(3, 1))
    expect_true(hasChromPeaks(res))
    expect_false(hasFeatures(res))

    res <- filterFile(xmseg, c(3, 1), keepFeatures = TRUE)
    expect_true(hasChromPeaks(res))
    expect_true(hasFeatures(res))
    expect_equal(featureValues(res), featureValues(xmseg)[, c(1, 3)])
})

test_that("adjustRtime,MsExperiment,XcmsExperiment,ObiwarpParam works", {
    op <- ObiwarpParam(binSize = 35.5)
    ref <- adjustRtime(faahko_xod, param = op)

    res <- adjustRtime(mse, param = op)
    expect_equal(unname(rtime(ref)), spectra(res)$rtime_adjusted)
    expect_s4_class(res, "XcmsExperiment")
    expect_true(length(res@processHistory) == 1L)
    expect_true(hasAdjustedRtime(res))
    expect_equal(rtime(res), unname(rtime(ref)))
    expect_equal(rtime(res, adjusted = FALSE),
                 unname(rtime(ref, adjusted = FALSE)))

    ## .plot_adjusted_rtime
    expect_warning(
        .plot_adjusted_rtime(rtime(res, adjusted = FALSE),
                             rtime(res), from_file = fromFile(res),
                             col = c("blue", "red")), "length"
        )

    ## applyAdjustedRtime
    res2 <- applyAdjustedRtime(res)
    expect_false(hasAdjustedRtime(res2))
    expect_equal(rtime(res2), rtime(res, adjusted = TRUE))
    res2 <- dropAdjustedRtime(res)
    expect_false(hasAdjustedRtime(res2))
    expect_true(length(res2@processHistory) == 0L)

    ## xcms object.
    res2 <- adjustRtime(xmse, param = op)
    expect_true(hasAdjustedRtime(res2))
    expect_equal(rtime(res, adjusted = FALSE), rtime(res2, adjusted = FALSE))
    expect_equal(rtime(res, adjusted = TRUE), rtime(res2, adjusted = TRUE))
    expect_equal(chromPeaks(ref), chromPeaks(res2))
    ## chrom peaks got adjusted too
    a <- chromPeaks(xmse)
    b <- chromPeaks(res2)
    expect_true(all(a[a[, "sample"] == 1L, "rt"] !=
                    b[b[, "sample"] == 1L, "rt"]))
    ## those of center sample are not changed
    expect_true(all(a[a[, "sample"] == 2L, "rt"] ==
                    b[b[, "sample"] == 2L, "rt"]))
    expect_true(length(res2@processHistory) == 2L)

    ## adjustedRtime
    expect_equal(adjustedRtime(res2), rtime(res2, adjusted = TRUE))

    ## Order: peak detection, alignment.
    ## dropAdjustedRtime:
    res3 <- dropAdjustedRtime(res2)
    expect_false(hasAdjustedRtime(res3))
    expect_true(hasChromPeaks(res3))
    ## chrom peak rt gets reverted
    expect_equal(chromPeaks(res3), chromPeaks(xmse))
    expect_true(length(res3@processHistory) == 1L)
    ref2 <- dropAdjustedRtime(ref)
    expect_equal(chromPeaks(res3), chromPeaks(ref2))
    ## dropChromPeaks
    res3 <- dropChromPeaks(res2)
    expect_false(hasChromPeaks(res3))
    expect_false(hasAdjustedRtime(res3))
    expect_true(length(res3@processHistory) == 0L)
    res3 <- dropChromPeaks(res2, keepAdjustedRtime = TRUE)
    expect_false(hasChromPeaks(res3))
    expect_true(hasAdjustedRtime(res3))
    expect_true(length(res3@processHistory) == 1L)
    expect_equal(rtime(res3, adjusted = TRUE), rtime(res2, adjusted = TRUE))

    ## Order: alignment, peak detection.
    res3 <- findChromPeaks(res, param = p)
    expect_true(hasChromPeaks(res3))
    expect_true(hasAdjustedRtime(res3))
    expect_true(length(res3@processHistory) == 2L)
    ## dropAdjustedRtime
    res4 <- dropAdjustedRtime(res3)
    expect_true(hasChromPeaks(res4))
    expect_false(hasAdjustedRtime(res4))
    expect_true(length(res4@processHistory) == 1L)
    ## chrom peak rt should be "reverted" (as with raw data) - but they
    ## are not identical because of the interpolation
    ## expect_equal(chromPeaks(res4), chromPeaks(xmse))
    ## a <- chromPeaks(res4)
    ## b <- chromPeaks(xmse)
    ## expect_true(all(a[a[, "sample"] == 1L, "rt"] ==
    ##                 b[b[, "sample"] == 1L, "rt"]))
    ## expect_true(all(a[a[, "sample"] == 2L, "rt"] ==
    ##                 b[b[, "sample"] == 2L, "rt"]))
    ## expect_true(all(a[a[, "sample"] == 3L, "rt"] ==
    ##                 b[b[, "sample"] == 3L, "rt"]))
    ## dropChromPeaks
    res4 <- dropChromPeaks(res3)
    expect_false(hasChromPeaks(res4))
    expect_true(hasAdjustedRtime(res4))
    expect_true(length(res4@processHistory) == 1L)

    ## With spectra that are NOT all associated to a sample.
    mse2 <- MsExperiment()
    sampleData(mse2) <- DataFrame(df)

    sps <- spectra(mse)
    tmp <- sps[1:10]
    tmp$dataOrigin <- "a"
    spectra(mse2) <- c(tmp, sps)
    mse2 <- linkSampleData(
        mse2, with = "sampleData.dataOrigin = spectra.dataOrigin")
    expect_error(adjustRtime(mse2, param = op), "to a sample")
})

test_that(".empty_feature_definitions works", {
    res <- .empty_feature_definitions()
    expect_true(is.data.frame(res))
    expect_true(nrow(res) == 0)
    expect_true(all(.REQ_PEAKG_COLS %in% colnames(res)))
})

## That's from XcmsExperiment-functions.R
test_that(".xmse_group_cpeaks works", {
    expect_error(.xmse_group_cpeaks(chromPeaks(xmse), p), "No correspondence")
    ## Just for PeakDensityParam.
    pdp <- PeakDensityParam(sampleGroups = rep(1, 3))
    cp <- chromPeaks(xmse, msLevel = 1L)
    res <- .xmse_group_cpeaks(cp, pdp)
    expect_true(is.data.frame(res))
    expect_true(all(.REQ_PEAKG_COLS %in% colnames(res)))
    expect_equal(res$mzmed, featureDefinitions(xod_xg)$mzmed)
    expect_equal(res$mzmin, featureDefinitions(xod_xg)$mzmin)
    expect_equal(res$mzmax, featureDefinitions(xod_xg)$mzmax)
    expect_equal(res$peakidx, featureDefinitions(xod_xg)$peakidx)

    res2 <- .xmse_group_cpeaks(cp, pdp, index = seq_len(nrow(cp)) + 13)
    idx <- lapply(res$peakidx, function(z) z + 13)
    expect_equal(idx, res2$peakidx)

    res <- .xmse_group_cpeaks(chromPeaks(xmse, msLevel = 2L), pdp)
    expect_true(all(.REQ_PEAKG_COLS %in% colnames(res)))

    ## NearestPeaksParam
    npp <- NearestPeaksParam(sampleGroups = c(1, 1))
    res <- .xmse_group_cpeaks(cp, npp)
    expect_true(is.data.frame(res))
    expect_true(all(.REQ_PEAKG_COLS %in% colnames(res)))
    expect_true(is.list(res$peakidx))

    ## MzClustParam
    cp <- chromPeaks(fticr_xod)
    mcp <- MzClustParam(sampleGroups = c(1, 1))
    res <- .xmse_group_cpeaks(cp, mcp)
    expect_true(is.data.frame(res))
    expect_true(is.list(res$peakidx))
})

test_that("groupChromPeaks,XcmsExperiment and related things work", {
    ## PeakDensityParam
    expect_false(hasFeatures(xmse))
    expect_false(hasFeatures(xmse, msLevel = 2L))
    pdp <- PeakDensityParam(sampleGroups = rep(1, 3))
    res <- groupChromPeaks(xmse, param = pdp, add = FALSE)
    expect_true(hasFeatures(res))
    expect_false(hasFeatures(res, msLevel = 2L))
    expect_equal(DataFrame(res@featureDefinitions), featureDefinitions(xod_xg))
    ## add FALSE
    res2 <- groupChromPeaks(res, param = pdp, add = FALSE)
    ## add TRUE
    res2 <- groupChromPeaks(res, param = pdp, add = TRUE)
    expect_true(length(res2@processHistory) > length(res@processHistory))
    expect_equal(res2@featureDefinitions[1:nrow(res@featureDefinitions), ],
                 res@featureDefinitions)
    a <- res2@featureDefinitions[1:nrow(res@featureDefinitions), ]
    b <- res2@featureDefinitions[(nrow(a) + 1):nrow(res2@featureDefinitions), ]
    expect_true(all(rownames(a) != rownames(b)))
    rownames(a) <- NULL
    rownames(b) <- NULL
    expect_equal(a, b)

    expect_error(groupChromPeaks(xmse, param = pdp, msLevel = 2L), "MS level 2")

    ## featureDefinitions
    expect_equal(featureDefinitions(xod_xg), DataFrame(featureDefinitions(res)))
    expect_equal(featureDefinitions(xod_xg, msLevel = 2L),
                 DataFrame(featureDefinitions(res, msLevel = 2L)))

    ## dropFeatureDefinitions
    res2 <- dropFeatureDefinitions(res)
    expect_false(hasFeatures(res2))
    expect_true(hasChromPeaks(res2))
    expect_true(length(res2@processHistory) == 1L)
    ## alignment before correspondence: keep adjusted rtime
    xmse2 <- adjustRtime(xmse, param = ObiwarpParam())
    res2 <- groupChromPeaks(xmse2, param = pdp)
    res3 <- dropFeatureDefinitions(res2)
    expect_false(hasFeatures(res3))
    expect_true(hasAdjustedRtime(res3))
    expect_true(length(res3@processHistory) == 2L)

    ## alignment after correspondence: drop adjusted rtime
    res2 <- adjustRtime(res, param = ObiwarpParam())
    expect_true(hasAdjustedRtime(res2))
    expect_true(hasFeatures(res2))
    res3 <- dropFeatureDefinitions(res2)
    expect_false(hasFeatures(res3))
    expect_false(hasAdjustedRtime(res3))
    expect_true(length(res3@processHistory) == 1L)
    expect_equal(chromPeaks(res3), chromPeaks(xmse))
    res3 <- dropFeatureDefinitions(res2, keepAdjustedRtime = TRUE)
    expect_false(hasFeatures(res3))
    expect_true(hasAdjustedRtime(res3))
    expect_true(length(res3@processHistory) == 2L)
    expect_equal(chromPeaks(res3), chromPeaks(res2))

    ## NearestPeaksParam
    npp <- NearestPeaksParam(sampleGroups = rep(1, 3), kNN = 3)
    ref <- groupChromPeaks(faahko_xod, param = npp)
    res <- groupChromPeaks(xmse, param = npp)
    expect_equal(featureDefinitions(ref), DataFrame(featureDefinitions(res)))
    expect_true(hasFeatures(res))
    expect_true(length(res@processHistory) == 2L)
})

test_that("adjustRtime,MsExperiment,PeakGroupsParam works", {
    a <- groupChromPeaks(xmse, param = PeakDensityParam(
                                   sampleGroups = c(1, 1, 1)))

    pgp <- PeakGroupsParam(span = 0.4)
    expect_false(hasAdjustedRtime(a))
    expect_true(hasFeatures(a))
    expect_error(adjustRtime(a, param = pgp, msLevel = 2L), "MS level 1")
    res <- adjustRtime(a, param = pgp)
    expect_true(hasAdjustedRtime(res))
    expect_false(hasFeatures(res))
    expect_equal(unname(rtime(xod_xgr)), unname(rtime(res)))
    expect_true(length(res@processHistory) == 3L)
})

test_that("findChromPeaks,XcmsExperiment,MatchedFilterParam works", {
    mfp <- MatchedFilterParam(binSize = 20, impute = "lin")
    ref <- findChromPeaks(faahko_od, param = mfp)
    res <- findChromPeaks(mse, param = mfp)
    expect_s4_class(res, "XcmsExperiment")
    expect_true(length(res@processHistory) == 1L)
    expect_equal(chromPeaks(res), chromPeaks(ref))
    expect_true(hasChromPeaks(res))

    res <- findChromPeaks(res, param = mfp, add = TRUE)
    expect_s4_class(res, "XcmsExperiment")
    expect_true(length(res@processHistory) == 2L)
    expect_equal(chromPeaks(res)[1:nrow(chromPeaks(ref)), ], chromPeaks(ref))
    expect_true(hasChromPeaks(res))

    res <- findChromPeaks(mse, param = mfp, msLevel = 2L)
    expect_s4_class(res, "XcmsExperiment")
    expect_true(length(res@processHistory) == 1L)
    expect_true(nrow(chromPeaks(res)) == 0L)
    expect_false(hasChromPeaks(res))
})

test_that("refineChromPeaks,XcmsExperiment,CleanPeaksParam works", {
    res <- refineChromPeaks(xmse, CleanPeaksParam(), msLevel = 2L)
    expect_equal(res, xmse)
    res <- refineChromPeaks(xmse, CleanPeaksParam(maxPeakwidth = 20))
    expect_true(length(res@processHistory) > length(xmse@processHistory))
    expect_true(nrow(chromPeaks(res)) < nrow(chromPeaks(xmse)))
    expect_equal(nrow(chromPeaks(res)), nrow(chromPeakData(res)))
    rtw <- chromPeaks(res)[, "rtmax"] - chromPeaks(res)[, "rtmin"]
    expect_true(all(rtw <= 20))
})

## That's from XcmsExperiment-functions.R
test_that(".merge_neighboring_peak_candidates works", {
    ## first file
    ref <- refineChromPeaks(filterFile(
        faahko_xod, 1L), MergeNeighboringPeaksParam(expandRt = 4))
    ref_pks <- chromPeaks(ref)[is.na(chromPeaks(ref)[, "intb"]), ]
    tmp <- xmse[1L, keepSampleIndex = TRUE, keepAdjustedRtime = TRUE]
    pd <- Spectra::peaksData(filterMsLevel(spectra(tmp), 1L))
    rt <- rtime(filterMsLevel(spectra(tmp), 1L))

    cand <- .define_merge_candidates(chromPeaks(tmp), expandRt = 4,
                                            expandMz = 0, ppm = 10)[[2L]]
    ## 1
    pks <- chromPeaks(tmp)[cand[[1L]], ]
    pkd <- chromPeakData(tmp)[cand[[1L]], ]
    res <- .merge_neighboring_peak_candidates(pd, rt, pks, pkd, diffRt = 8,
                                              ppm = 10, expandMz = 0)
    expect_equal(unname(res$chromPeaks), unname(ref_pks[1, , drop = FALSE]))

    ## 3
    pks <- chromPeaks(tmp)[cand[[3L]], ]
    pkd <- chromPeakData(tmp)[cand[[3L]], ]
    res <- .merge_neighboring_peak_candidates(pd, rt, pks, pkd, diffRt = 8,
                                              ppm = 10, expandMz = 0)
    expect_equal(unname(res$chromPeaks[1, ]),
                 unname(ref_pks[2, ]))
    expect_equal(unname(res$chromPeaks[3, ]),
                 unname(ref_pks[3, ]))

    ## Second file
    ref <- refineChromPeaks(filterFile(
        faahko_xod, 2L), MergeNeighboringPeaksParam(expandRt = 4))
    ref_pks <- chromPeaks(ref)[is.na(chromPeaks(ref)[, "intb"]), ]
    ref_pks["sample"] <- 2L
    tmp <- xmse[2L, keepSampleIndex = TRUE, keepAdjustedRtime = TRUE]
    pd <- Spectra::peaksData(filterMsLevel(spectra(tmp), 1L))
    rt <- rtime(filterMsLevel(spectra(tmp), 1L))

    cand <- .define_merge_candidates(chromPeaks(tmp), expandRt = 4,
                                            expandMz = 0, ppm = 10)[[2L]]
    pks <- chromPeaks(tmp)[cand[[1L]], ]
    pkd <- chromPeakData(tmp)[cand[[1L]], ]
    res <- .merge_neighboring_peak_candidates(pd, rt, pks, pkd, diffRt = 8,
                                              ppm = 10, expandMz = 0)
    expect_equal(names(res), c("chromPeaks", "chromPeakData"))
    expect_equal(nrow(res$chromPeaks), 1)
    expect_equal(res$chromPeaks[1, ], ref_pks)

    ## Not merging
    pks <- chromPeaks(tmp)[cand[[2L]], ]
    pkd <- chromPeakData(tmp)[cand[[2L]], ]
    res <- .merge_neighboring_peak_candidates(pd, rt, pks, pkd, diffRt = 8,
                                              ppm = 10, expandMz = 0)
    expect_equal(res$chromPeaks, pks)

    pks <- chromPeaks(tmp)[cand[[3L]], ]
    pkd <- chromPeakData(tmp)[cand[[3L]], ]
    res <- .merge_neighboring_peak_candidates(pd, rt, pks, pkd, diffRt = 8,
                                              ppm = 10, expandMz = 0)
    expect_equal(res$chromPeaks, pks)
})

## That's from XcmsExperiment-functions.R
test_that(".merge_neighboring_peaks2 works", {
    tmp1 <- filterFile(faahko_xod, 1L)
    tmp2 <- xmse[1L]
    x <- Spectra::peaksData(filterMsLevel(spectra(tmp2), 1L))
    pks <- chromPeaks(tmp2, msLevel = 1L)
    pkd <- chromPeaks(tmp2, msLevel = 1L)
    rt <- rtime(tmp2)[msLevel(spectra(tmp2)) == 1L]
    prm <- MergeNeighboringPeaksParam(expandRt = 6, expandMz = 1)

    ref <- refineChromPeaks(tmp1, prm)
    res <- .merge_neighboring_peaks2(x, pks, pkd, rt,
                                     expandRt = prm@expandRt,
                                     expandMz = prm@expandMz)
    expect_equal(unname(chromPeaks(ref)), unname(res$chromPeaks))

    prm <- MergeNeighboringPeaksParam(expandRt = 4, expandMz = 1)
    ref <- refineChromPeaks(tmp1, prm)
    res <- .merge_neighboring_peaks2(x, pks, pkd, rt,
                                     expandRt = prm@expandRt,
                                     expandMz = prm@expandMz)
    expect_equal(unname(chromPeaks(ref)), unname(res$chromPeaks))

    tmp1 <- filterFile(faahko_xod, 2L)
    tmp2 <- xmse[2L]
    x <- Spectra::peaksData(filterMsLevel(spectra(tmp2), 1L))
    pks <- chromPeaks(tmp2, msLevel = 1L)
    pkd <- chromPeaks(tmp2, msLevel = 1L)
    rt <- rtime(tmp2)[msLevel(spectra(tmp2)) == 1L]
    ref <- refineChromPeaks(tmp1, prm)
    res <- .merge_neighboring_peaks2(x, pks, pkd, rt,
                                     expandRt = prm@expandRt,
                                     expandMz = prm@expandMz)
    expect_equal(unname(chromPeaks(ref)), unname(res$chromPeaks))
})

## That's from XcmsExperiment-functions.R
test_that(".xmse_merge_neighboring_peaks etc works", {
    ref <- refineChromPeaks(faahko_xod, MergeNeighboringPeaksParam(
                                            expandRt = 6, expandMz = 1))
    a <- .xmse_merge_neighboring_peaks(xmse, expandRt = 6, expandMz = 1)
    expect_true(nrow(a$chromPeaks) == nrow(chromPeaks(ref)))
    expect_true(nrow(a$chromPeakData) == nrow(chromPeaks(ref)))
})

## That's from XcmsExperiment-functions.R
test_that(".xmse_apply_chunks works", {
    res <- .xmse_apply_chunks(xmse, FUN = identity, chunkSize = 2L)
    expect_true(length(res) == 2)
    expect_s4_class(res[[1L]], "XcmsExperiment")
    expect_s4_class(res[[2L]], "XcmsExperiment")
    expect_true(length(res[[1L]]) == 2)
    expect_true(length(res[[2L]]) == 1)
})

test_that("refineChromPeaks,XcmsExperiment,MergedChromPeaksParam works", {
    prm <- MergeNeighboringPeaksParam(expandRt = 4, ppm = 20)
    ref <- refineChromPeaks(faahko_xod, param = prm)

    res <- refineChromPeaks(xmse, param = prm)
    expect_true(validObject(res))
    expect_equal(unname(chromPeaks(res)), unname(chromPeaks(ref)))
    a <- chromPeakData(res)
    b <- chromPeakData(ref)
    rownames(a) <- NULL
    rownames(b) <- NULL
    expect_equal(a, b)
    expect_equal(rownames(chromPeaks(res)), rownames(chromPeakData(res)))
    expect_true(length(res@processHistory) > length(xmse@processHistory))

    prm <- MergeNeighboringPeaksParam(expandRt = 6, expandMz = 1)
    ref <- refineChromPeaks(faahko_xod, param = prm)

    res <- refineChromPeaks(xmse, param = prm)
    expect_true(validObject(res))
    expect_equal(unname(chromPeaks(res)), unname(chromPeaks(ref)))
    a <- chromPeakData(res)
    b <- chromPeakData(ref)
    rownames(a) <- NULL
    rownames(b) <- NULL
    expect_equal(a, b)
    expect_equal(rownames(chromPeaks(res)), rownames(chromPeakData(res)))
    expect_true(length(res@processHistory) > length(xmse@processHistory))
})

## That's from XcmsExperiment-functions.R
test_that(".xmse_filter_peaks_intensities works", {
    res <- .xmse_filter_peaks_intensities(xmse, nValues = 4, threshold = 0,
                                          msLevel = 1L)
    expect_true(is.logical(res))
    expect_true(length(res) == nrow(chromPeaks(xmse)))
    expect_true(all(res))

    res <- .xmse_filter_peaks_intensities(xmse, nValues = 20, threshold = 0,
                                          msLevel = 1L)
    expect_false(all(res))

    res <- .xmse_filter_peaks_intensities(xmse, nValues = 1, threshold = 50000)
    expect_equal(res, unname(chromPeaks(xmse)[, "maxo"] >= 50000))

    res <- .xmse_filter_peaks_intensities(xmse, nValues = 1, , msLeve = 2L)
    expect_true(length(res) == 0)
})

test_that("refineChromPeaks,XcmsExperiment,FilterIntensityParam works", {
    fip <- FilterIntensityParam(threshold = 13000, nValues = 3)
    ref <- refineChromPeaks(faahko_xod, fip)

    res <- refineChromPeaks(xmse, fip)
    expect_equal(chromPeaks(res), chromPeaks(ref))
    expect_true(nrow(chromPeaks(res)) < nrow(chromPeaks(xmse)))
    expect_true(length(res@processHistory) > length(xmse@processHistory))

    expect_warning(res <- refineChromPeaks(xmse, fip, msLevel = 3L), "level 3")
    expect_equal(chromPeaks(res), chromPeaks(xmse))

    expect_error(
        refineChromPeaks(xmse, FilterIntensityParam(1000, value = "other")),
        "not available")
    res <- refineChromPeaks(xmse, FilterIntensityParam(300000, nValues = 1,
                                                       value = "into"))
    expect_true(all(chromPeaks(res)[, "into"] >= 300000))
})

test_that("featureValues,XcmsExperiment works", {
    expect_error(featureValues(xmse), "feature definitions")

    expect_error(featureValues(xmseg, msLevel = 2L), "feature definitions")
    res <- featureValues(xmseg)
    expect_true(is.matrix(res))
    expect_equal(colnames(res), c("ko15.CDF", "ko16.CDF", "ko18.CDF"))
    expect_equal(rownames(res), rownames(featureDefinitions(xmseg)))
    res2 <- featureValues(xmseg, missing = 10)
    expect_true(all(res2[is.na(res)] == 10))

    expect_error(featureValues(xmseg, method = "sum", value = "index"),
                 "value is set to")
    expect_error(featureValues(xmseg, missing = "sum"), "or a numeric")
})

test_that("fillChromPeaks,XcmsExperiment,ChromPeakAreaParam works", {
    cpp <- ChromPeakAreaParam()
    ## pdp <- PeakDensityParam(sampleGroups = rep(1, 3))
    ## xmseg <- groupChromPeaks(xmse, param = pdp, add = FALSE)

    pal <- split.data.frame(chromPeaks(xmseg), chromPeaks(xmseg)[, "sample"])
    res <- .xmse_integrate_chrom_peaks(xmse, pal)
    expect_true(is.matrix(res))
    expect_equal(unname(res[, "into"]), unname(chromPeaks(xmse)[, "into"]))
    expect_equal(unname(res[, "maxo"]), unname(chromPeaks(xmse)[, "maxo"]))

    expect_error(fillChromPeaks(xmse, param = cpp), "MS level 1")
    expect_error(fillChromPeaks(xmseg, param = cpp, msLevel = 2L), "MS level 2")

    res <- fillChromPeaks(xmseg, param = cpp)
    expect_true(length(res@processHistory) > length(xmseg@processHistory))
    expect_true(nrow(chromPeaks(res)) > nrow(chromPeaks(xmseg)))
    expect_true(nrow(chromPeakData(res)) > nrow(chromPeakData(xmseg)))
    expect_true(sum(is.na(featureValues(res))) <
                sum(is.na(featureValues(xmseg))))
    expect_true(hasFilledChromPeaks(res))
    res <- dropFilledChromPeaks(res)
    expect_false(hasFilledChromPeaks(res))
    expect_equal(chromPeaks(res), chromPeaks(xmseg))
    expect_equal(featureDefinitions(res), featureDefinitions(xmseg))
    expect_true(length(res@processHistory) == length(xmseg@processHistory))

    ## With matched filter.
    mfp <- MatchedFilterParam(binSize = 0.2)
    tmp <- findChromPeaks(mse, mfp)
    tmp <- groupChromPeaks(tmp, pdp)
    res <- fillChromPeaks(tmp, cpp)
    expect_true(length(res@processHistory) > length(tmp@processHistory))
    expect_true(nrow(chromPeaks(res)) > nrow(chromPeaks(tmp)))
    expect_true(nrow(chromPeakData(res)) > nrow(chromPeakData(tmp)))
    expect_true(sum(is.na(featureValues(res))) <
                sum(is.na(featureValues(tmp))))
    expect_true(hasFilledChromPeaks(res))
})

## That's from XcmsExperiment-functions.R
test_that(".xmse_process_history works", {
    res <- .xmse_process_history(xmse)
    expect_equal(res, xmse@processHistory)
    res <- .xmse_process_history(xmse, msLevel = 1L)
    expect_equal(res, xmse@processHistory)
    res <- .xmse_process_history(xmse, msLevel = 2L)
    expect_equal(res, list())
    res <- .xmse_process_history(xmse, type = .PROCSTEP.PEAK.DETECTION)
    expect_equal(res, xmse@processHistory)
    res <- .xmse_process_history(xmse, type = .PROCSTEP.PEAK.DETECTION,
                                 msLevel = 2L)
    expect_equal(res, list())
    res <- .xmse_process_history(xmse, type = .PROCSTEP.FEATURE.GROUPING)
    expect_equal(res, list())
})

## That's from XcmsExperiment-functions.R
test_that(".chrom_peak_intensity_centWave works", {
    x <- Spectra::peaksData(spectra(xmse[2L]))
    rt <- rtime(spectra(xmse[2L]))
    pks <- chromPeaks(xmse)[chromPeaks(xmse)[, "sample"] == 2L, ]

    res <- .chrom_peak_intensity_centWave(x, rt, pks, sampleIndex = 2L,
                                          cn = colnames(pks))
    expect_equal(unname(res[, "mz"]), unname(pks[, "mz"]))
    ## expect_equal(res[, "rt"], unname(pks[, "rt"])) # that is different.
    expect_equal(unname(res[, "into"]), unname(pks[, "into"]))
    expect_equal(unname(res[, "maxo"]), unname(pks[, "maxo"]))
})

## That's from XcmsExperiment-functions.R
test_that(".chrom_peak_intensity_matchedFilter works", {
    x <- Spectra::peaksData(spectra(xmse[2L]))
    rt <- rtime(spectra(xmse[2L]))

    tmp <- findChromPeaks(mse[2L], param = MatchedFilterParam())
    pks <- chromPeaks(tmp)
    res <- .chrom_peak_intensity_matchedFilter(x, rt, pks, cn = colnames(pks),
                                               sampleIndex = 2L)
    ## expect_equal(res[, "rt"], pks[, "rt"]) # not the same: no gauss filter
    expect_equal(res[, "mz"], pks[, "mz"], tolerance = 0.0001)
    expect_equal(res[, "into"], pks[, "into"])
    expect_equal(res[, "maxo"], pks[, "maxo"])
})

## That's from XcmsExperiment-functions.R
test_that(".filter_chrom_peaks works", {
    res <- .filter_chrom_peaks(xmse, idx = c(4, 2, 34, 1))
    expect_s4_class(res, "XcmsExperiment")
    expect_true(nrow(chromPeaks(res)) == 4)
    expect_true(nrow(chromPeakData(res)) == 4)
    expect_equal(rownames(chromPeaks(res)),
                 c("CP004", "CP002", "CP034", "CP001"))
    ## with feature data.
    res <- .filter_chrom_peaks(xmseg, c(11, 199, 115, 205, 212))
    expect_true(hasFeatures(res))
    expect_equal(rownames(chromPeaks(res)), c("CP011", "CP199", "CP115",
                                              "CP205", "CP212"))
    expect_equal(featureDefinitions(res)$peakidx,
                 list(c(1L, 2L), c(3L, 4L), 5L))
})

## That's from XcmsExperiment-functions.R
test_that(".spectra_index_list works", {
    sp <- spectra(xmse[1L])
    pks <- chromPeaks(xmse[1L])

    res <- .spectra_index_list(sp, pks, msLevel = 1L)
    expect_true(is.list(res))
    expect_true(length(res) == nrow(pks))
    rt <- rtime(sp)
    for (i in seq_len(nrow(pks))) {
        expect_true(all(rt[res[[i]]] >= pks[i, "rtmin"]))
        expect_true(all(rt[res[[i]]] <= pks[i, "rtmax"]))
    }
    res <- .spectra_index_list(sp, cbind(rtmin = 10000, rtmax = 20000), 1L)
    expect_equal(res, list(integer()))

    res <- .spectra_index_list(sp, pks, msLevel = 3)
    expect_true(all(lengths(res) == 0))
})

## That's from XcmsExperiment-functions.R
test_that(".spectra_index_list_closest_rt works", {
    sp <- spectra(xmse[1L])
    pks <- chromPeaks(xmse[1L])

    res <- .spectra_index_list_closest_rt(sp, pks, msLevel = 1L)
    expect_true(is.list(res))
    expect_true(length(res) == nrow(pks))
    expect_true(all(lengths(res) == 1L))
    rt <- rtime(sp)
    for (i in seq_len(nrow(pks))) {
        expect_true(all(rt[res[[i]]] >= pks[i, "rtmin"]))
        expect_true(all(rt[res[[i]]] <= pks[i, "rtmax"]))
    }
    diffs <- rt[unlist(res)] - pks[, "rt"]
    expect_true(all(diffs == 0))
})

## That's from XcmsExperiment-functions.R
test_that(".spectra_index_list_closest_mz works", {
    sp <- spectra(xmse[1L])
    pks <- chromPeaks(xmse[1L])
    sp$precursorMz[65:75] <- pks[2, "mz"]

    res <- .spectra_index_list_closest_mz(sp, pks)
    expect_true(is.list(res))
    expect_true(length(res) == nrow(pks))
    expect_true(sum(lengths(res)) == 1L)
    expect_equal(res[[2]], 65L)
})

## That's from XcmsExperiment-functions.R
test_that(".mse_spectra_for_peaks works", {
    res <- .mse_spectra_for_peaks(xmse)
    expect_s4_class(res, "Spectra")
    expect_true(any(spectraVariables(res) == "peak_id"))
    expect_true(length(res) == 0)

    res <- .mse_spectra_for_peaks(xmse, msLevel = 1L, method = "closest_rt")
    expect_s4_class(res, "Spectra")
    expect_true(any(spectraVariables(res) == "peak_id"))
    expect_true(length(res) == nrow(chromPeaks(xmse)))

    res <- .mse_spectra_for_peaks(xmse, msLevel = 1L, method = "all",
                                  peaks = 220)
    expect_true(all(res$peak_id == "CP220"))

    ## Duplicates index?
    res <- .mse_spectra_for_peaks(xmse, msLevel = 1L, method = "closest_rt",
                                  peaks = c(3, 2, 3, 3, 1))
    expect_equal(res$peak_id, c("CP003", "CP002", "CP003", "CP003", "CP001"))
    expect_equal(rtime(res)[1], rtime(res)[3])
    expect_equal(mz(res)[1], mz(res)[3])
})

test_that("chromPeakSpectra works", {
    ## input errors
    expect_error(chromPeakSpectra(xmse, method = "other"), "'arg' should be")
    expect_error(chromPeakSpectra(xmse, return.type = "list"), "'arg' should")
    expect_error(chromPeakSpectra(xmse, peaks = "other"), "out of bounds")

    pks <- c("CP242", "CP007", "CP123")
    res <- chromPeakSpectra(xmse, peaks = pks)
    expect_s4_class(res, "Spectra")
    expect_equal(length(res), 0)
    res <- chromPeakSpectra(xmse, peaks = pks, msLevel = 1L,
                            return.type = "List")
    expect_s4_class(res, "List")
    expect_equal(names(res), pks)
    res <- chromPeakSpectra(xmse, peaks = pks, msLevel = 1L)
    expect_s4_class(res, "Spectra")
    expect_equal(unique(res$peak_id), pks)

    res2 <- chromPeakSpectra(xmse, msLevel = 1L, method = "closest_rt")
    expect_equal(length(res2), nrow(chromPeaks(xmse)))
    expect_equal(res2$peak_id, rownames(chromPeaks(xmse)))

    res2 <- chromPeakSpectra(xmse, msLevel = 1L, method = "largest_tic",
                             peaks = pks)
    expect_equal(length(res2), length(pks))
    expect_equal(res2$peak_id, pks)
    ic <- split(ionCount(res), factor(res$peak_id, levels = pks))
    idx <- vapply(ic, which.max, integer(1))
    expect_equal(rtime(res2[1L]), rtime(res[res$peak_id == pks[1L]])[idx[1L]])
    expect_equal(rtime(res2[2L]), rtime(res[res$peak_id == pks[2L]])[idx[2L]])
    expect_equal(rtime(res2[3L]), rtime(res[res$peak_id == pks[3L]])[idx[3L]])

    res2 <- chromPeakSpectra(xmse, msLevel = 1L, method = "largest_bpi",
                             peaks = pks, return.type = "List")
    expect_equal(length(res2), length(pks))
    expect_equal(names(res2), pks)
    expect_true(all(lengths(res2) == 1L))
    bpi <- split(max(intensity(res)), factor(res$peak_id, levels = pks))
    idx <- vapply(bpi, which.max, integer(1))
    expect_equal(rtime(res2[[1L]]), rtime(res[res$peak_id == pks[1L]])[idx[1L]])
    expect_equal(rtime(res2[[2L]]), rtime(res[res$peak_id == pks[2L]])[idx[2L]])
    expect_equal(rtime(res2[[3L]]), rtime(res[res$peak_id == pks[3L]])[idx[3L]])

    ## DDA data
    fl <- system.file("TripleTOF-SWATH/PestMix1_DDA.mzML", package = "msdata")
    tmp <- readMsExperiment(fl)
    tmp <- filterRt(tmp, c(200, 400))
    tmp <- findChromPeaks(tmp, CentWaveParam(peakwidth = c(5, 15),
                                             prefilter = c(5, 1000)))
    res <- chromPeakSpectra(tmp, return.type = "List")
    expect_equal(length(res), nrow(chromPeaks(tmp)))
    expect_true(all(precursorMz(res[[1L]]) >= chromPeaks(tmp)[1, "mzmin"] &
                    precursorMz(res[[1L]]) <= chromPeaks(tmp)[1, "mzmax"]))
    expect_equal(rtime(chromPeakSpectra(tmp, peaks = c("CP7", "CP1", "CP3"))),
                 rtime(chromPeakSpectra(tmp, peaks = c(7, 1, 3))))
})

test_that("manualChromPeaks,XcmsExperiment works", {
    pks <- chromPeaks(xmse)[chromPeaks(xmse)[, "sample"] == 2L, ]
    res <- manualChromPeaks(xmse, chromPeaks = matrix(numeric()))
    expect_equal(chromPeaks(res), chromPeaks(xmse))
    tmp <- as(mse, "XcmsExperiment")
    expect_error(manualChromPeaks(xmse, msLevel = c(1L, 2L)), "one MS")
    expect_error(manualChromPeaks(xmse, pks[, c("rt", "mz")]), "required")
    expect_error(manualChromPeaks(xmse, pks, samples = 2:5), "out of bounds")

    res <- manualChromPeaks(tmp, pks)
    expect_true(hasChromPeaks(res))
    pks_2 <- chromPeaks(res)[chromPeaks(res)[, "sample"] == 2L, ]
    expect_equal(nrow(pks), nrow(pks_2))
    expect_equal(unname(pks[, c("mz", "into", "maxo")]),
                 unname(pks_2[, c("mz", "into", "maxo")]))

    res2 <- manualChromPeaks(tmp, pks, samples = 2)
    expect_equal(unname(chromPeaks(res2)), unname(pks_2))
})

test_that("filterChromPeaks,XcmsExperiment works", {
    res <- filterChromPeaks(xmse, keep = 4:9)
    expect_equal(nrow(chromPeaks(res)), 6)
    expect_equal(chromPeaks(res), chromPeaks(xmse)[4:9, ])

    ## With features
    res <- filterChromPeaks(xmseg, keep = 13:60)
    expect_equal(chromPeaks(res), chromPeaks(xmseg)[13:60, ])

    a <- featureValues(res)
    b <- featureValues(xmseg)
    expect_equal(a[, 1], b[rownames(a), 1])
    expect_true(all(unlist(featureDefinitions(res)$peakidx) %in%
                    seq_len(nrow(chromPeaks(res)))))
})

## That's from XcmsExperiment-functions.R
test_that(".manual_feature_definitions works", {
    idx <- list(c(13, 15, 220), c(45, 46, 100, 200))
    expect_error(
        .manual_feature_definitions(chromPeaks(xmse),
                                    list(c(1, 2, 1000), c(1, 2))),
        "out of bounds")
    res <- .manual_feature_definitions(chromPeaks(xmse), idx)
    expect_true(is.data.frame(res))
    expect_true(nrow(res) == 2)
    expect_equal(colnames(res), c("mzmed", "mzmin", "mzmax", "rtmed",
                                  "rtmin", "rtmax", "npeaks", "peakidx"))
    expect_equal(res$npeaks, lengths(idx))
})

test_that("manualFeatures,XcmsExperiment works", {
    idx <- list(c(3), c(1, 5, 9), c(34, 121, 247))
    expect_error(manualFeatures(as(mse, "XcmsExperiment"), idx), "present")
    res <- manualFeatures(xmse)
    expect_equal(res, xmse)
    expect_error(manualFeatures(xmse, idx, msLevel = 1:23), "at a time")
    res <- manualFeatures(xmse, idx)
    expect_true(hasFeatures(res))
    expect_true(nrow(featureDefinitions(res)) == 3)
    expect_equal(featureDefinitions(res)$peakidx, idx)
    expect_equal(unname(featureValues(res)[3, ]),
                 unname(chromPeaks(res)[idx[[3]], "into"]))

    res2 <- manualFeatures(xmseg, idx)
    expect_true(nrow(featureDefinitions(res2)) ==
                (nrow(featureDefinitions(xmseg)) + 3))
    expect_equal(
        featureDefinitions(res2)[seq_len(nrow(featureDefinitions(xmseg))), ],
        featureDefinitions(xmseg))
})

test_that("filterFeatureDefinitions works", {
    expect_error(filterFeatureDefinitions(xmse, 1:3), "No feature")
    res <- filterFeatureDefinitions(xmseg)
    expect_equal(res, xmseg)

    expect_error(filterFeatureDefinitions(xmseg, c("a", "b")), "bounds")
    res <- filterFeatureDefinitions(xmseg, 1:10)
    expect_true(hasFeatures(res))
    expect_equal(featureDefinitions(res), featureDefinitions(xmseg)[1:10, ])
})

test_that("featureSpectra works", {
    expect_error(featureSpectra(xmse), "No feature definitions")
    res_all <- featureSpectra(xmseg, msLevel = 1L)
    expect_s4_class(res_all, "Spectra")
    expect_true(all(rownames(featureDefinitions(xmseg)) %in%
                    res_all$feature_id))

    res_all <- featureSpectra(xmseg, msLevel = 1L, method = "closest_rt",
                              return.type = "List")
    expect_s4_class(res_all, "List")
    expect_equal(length(res_all), nrow(featureDefinitions(xmseg)))
    expect_equal(names(res_all), rownames(featureDefinitions(xmseg)))

    res <- featureSpectra(xmseg, msLevel = 1L, features = c(3, 1, 2, 1, 5),
                          method = "closest_rt", return.type = "List")
    expect_true(length(res) == 5)
    expect_equal(rtime(res[[1L]]), rtime(res_all[[3L]]))
    expect_equal(rtime(res[[2L]]), rtime(res_all[[1L]]))
    expect_equal(rtime(res[[3L]]), rtime(res_all[[2L]]))
    expect_equal(rtime(res[[4L]]), rtime(res_all[[1L]]))
    expect_equal(rtime(res[[5L]]), rtime(res_all[[5L]]))

    res_2 <- featureSpectra(xmseg, msLevel = 1L, features = c("FT03", "FT01"),
                            return.type = "List")
    expect_true(length(res[[1L]]) < length(res_2[[1L]]))
    expect_true(all(res[[1L]]$peak_id %in% res_2[[1L]]$peak_id))
    expect_equal(unique(res[[1L]]$feature_id), unique(res_2[[1L]]$feature_id))
    expect_true(length(res[[2L]]) < length(res_2[[2L]]))
    expect_true(all(res[[2L]]$peak_id %in% res_2[[2L]]$peak_id))
    expect_equal(unique(res[[2L]]$feature_id), unique(res_2[[2L]]$feature_id))
})

test_that("chromatogram,XcmsExperiment and .xmse_extract_chromatograms_old", {
    expect_error(chromatogram(xmse, adjustedRtime = FALSE), "unused")
    expect_warning(res <- chromatogram(xmse, include = "apex_within",
                                       return.type = "MChromatograms"),
                   "deprecated")
    expect_s4_class(res, "MChromatograms")
    expect_true(nrow(res) == 1L)
    ref <- chromatogram(faahko_od)
    expect_equal(intensity(res[1, 1]), unname(intensity(ref[1, 1])))

    rtr <- c(2600, 2700)
    mzr <- c(340, 400)
    res <- chromatogram(xmse, mz = mzr, rt = rtr)
    expect_s4_class(res, "XChromatograms")
    expect_true(nrow(res) == 1L)
    expect_true(nrow(chromPeaks(res)) > 0)
    expect_true(all(chromPeaks(res)[, "mz"] >= 340 &
                    chromPeaks(res)[, "mz"] <= 400))
    expect_true(all(chromPeaks(res[1, 1])[, "sample"] == 1L))
    expect_true(all(chromPeaks(res[1, 2])[, "sample"] == 2L))
    expect_true(all(chromPeaks(res[1, 3])[, "sample"] == 3L))
    ref <- chromatogram(xod_x, mz = mzr, rt = rtr)
    expect_equal(chromPeaks(res), chromPeaks(ref))

    ## Multiple rows.
    res <- .xmse_extract_chromatograms_old(
        xmse, mz = chromPeaks(xmse)[1:10, c("mzmin", "mzmax")],
        rt = chromPeaks(xmse)[1:10, c("rtmin", "rtmax")],
        msLevel = 1L, aggregationFun = "sum", chunkSize = 2L,
        isolationWindow = NULL,
        chromPeaks = "apex_within", BPPARAM = bpparam(),
        return.type = "XChromatograms")
    expect_true(nrow(res) == 10L)

    ## with features
    res <- .xmse_extract_chromatograms_old(
        xmseg, mz = chromPeaks(xmse)[1:5, c("mzmin", "mzmax")],
        rt = chromPeaks(xmse)[1:5, c("rtmin", "rtmax")], chunkSize = 2L,
        BPPARAM = bpparam(), msLevel = 1L, aggregationFun = "sum",
        isolationWindow = NULL,
        chromPeaks = "apex_within", return.type = "XChromatograms")
    expect_true(nrow(featureDefinitions(res)) == 3)
    expect_true(all(unlist(featureDefinitions(res)$peakidx) %in%
                    seq_len(nrow(chromPeaks(res)))))
    ref <- chromatogram(xod_xg,
                        mz = chromPeaks(xmse)[1:5, c("mzmin", "mzmax")],
                        rt = chromPeaks(xmse)[1:5, c("rtmin", "rtmax")])
    expect_equal(featureDefinitions(ref), featureDefinitions(res))

    ## MS2 data.
    res <- chromatogram(xmseg, msLevel = 2L,
                        mz = chromPeaks(xmse)[1:5, c("mzmin", "mzmax")],
                        rt = chromPeaks(xmse)[1:5, c("rtmin", "rtmax")])
    expect_true(validObject(res))
    expect_true(length(intensity(res[[1L]])) == 0)
    expect_true(length(intensity(res[[2L]])) == 0)
    expect_s4_class(res, "XChromatograms")
    expect_true(nrow(chromPeaks(res)) == 0)

    ## real MS2 data.
    res <- chromatogram(mse_dia, msLevel = 2L, mz = c(50, 300),
                        rt = c(100, 600))
    expect_true(length(intensity(res[[1L]])) == 0)
    res <- chromatogram(mse_dia, msLevel = 2L, mz = c(50, 300),
                        rt = c(100, 600), isolationWindowTargetMz = 270.85,
                        aggregationFun = "sum")
    expect_true(all(intensity(res[[1L]]) > 0))

    ## Defining only mz or rt.
    rtr <- c(2600, 2700)
    mzr <- c(340, 400)
    res <- chromatogram(xmse, mz = mzr)
    expect_s4_class(res, "XChromatograms")
    expect_true(nrow(res) == 1L)
    expect_true(nrow(chromPeaks(res)) > 0)
    expect_true(all(chromPeaks(res)[, "mz"] >= 340 &
                    chromPeaks(res)[, "mz"] <= 400))
    expect_true(all(chromPeaks(res[1, 1])[, "sample"] == 1L))
    expect_true(all(chromPeaks(res[1, 2])[, "sample"] == 2L))
    expect_true(all(chromPeaks(res[1, 3])[, "sample"] == 3L))
    rrt <- range(lapply(res, rtime))
    expect_true(rrt[1] < 2600)
    expect_true(rrt[2] > 4400)

    res <- chromatogram(xmse, rt = rtr)
    expect_s4_class(res, "XChromatograms")
    expect_true(nrow(res) == 1L)
    expect_true(nrow(chromPeaks(res)) > 0)
    expect_true(any(chromPeaks(res)[, "mz"] < 340 |
                    chromPeaks(res)[, "mz"] > 400))
    expect_true(all(chromPeaks(res[1, 1])[, "sample"] == 1L))
    expect_true(all(chromPeaks(res[1, 2])[, "sample"] == 2L))
    expect_true(all(chromPeaks(res[1, 3])[, "sample"] == 3L))
    rrt <- range(lapply(res, rtime))
    expect_true(rrt[1] >= 2600)
    expect_true(rrt[2] <= 2700)
})

test_that("featureChromatograms,XcmsExperiment works", {
    res <- featureChromatograms(xmseg, return.type = "XChromatograms")
    expect_s4_class(res, "XChromatograms")
    expect_true(validObject(res))
    expect_equal(nrow(res), nrow(featureDefinitions(xmseg)))
    expect_equal(rownames(featureDefinitions(res)),
                 rownames(featureDefinitions(xmseg)))
    expect_equal(featureValues(res), featureValues(xmseg))

    ref <- featureChromatograms(xod_xg)
    expect_true(validObject(ref))
    expect_equal(chromPeaks(ref), chromPeaks(res))
    expect_equal(featureDefinitions(ref)$peakidx,
                 featureDefinitions(res)$peakidx)
    expect_equal(featureDefinitions(res),
                 featureDefinitions(ref))

    expect_error(featureChromatograms(xmseg, features = "a"), "out of")

    ## Duplicated features
    res <- featureChromatograms(xmseg, features = c("FT12", "FT03", "FT12"))
    expect_s4_class(res, "XChromatograms")
    expect_true(validObject(res))
    expect_true(nrow(res) == 3L)
    expect_equal(rownames(featureDefinitions(res)), c("FT12", "FT03", "FT12.1"))
    expect_equal(chromPeaks(res)[featureDefinitions(res)$peakidx[[1L]], 1:11],
                 chromPeaks(res)[featureDefinitions(res)$peakidx[[3L]], 1:11])
    expect_equal(unname(chromPeaks(res)[, "row"]),
                 c(1, 1, 1, 2, 2, 2, 2, 3, 3, 3))
})

test_that("processHistory,XcmsExperiment works", {
    res <- processHistory(new("XcmsExperiment"))
    expect_equal(res, list())
    res <- processHistory(xmse)
    expect_true(length(res) == 1)

    res2 <- processHistory(xmseg, type = .PROCSTEP.PEAK.DETECTION)
    expect_equal(res, res2)

    res <- processHistory(xmseg)
    expect_true(length(res) == 2)

    res <- processHistory(xmseg, type = "other")
    expect_true(length(res) == 0)
})

test_that("filterMzRange,XcmsExperiment works", {
    mzr <- c(340, 350)
    res <- filterMzRange(xmse, mzr)
    expect_s4_class(res, "XcmsExperiment")
    mzs <- unlist(mz(spectra(res)))
    expect_true(all(mzs >= 340 & mzs <= 350))
    expect_true(all(chromPeaks(res)[, "mz"] >= 340 &
                    chromPeaks(res)[, "mz"] <= 350))

    expect_warning(res2 <- filterMzRange(xmse, mzr, msLevel. = 2L), "not")
    expect_s4_class(res, "XcmsExperiment")
    expect_equal(mz(spectra(res2[1L])), mz(spectra(xmse[1L])))
    expect_equal(chromPeaks(xmse), chromPeaks(res2))

    res <- filterMzRange(xmse)
    expect_equal(res, xmse)
})

test_that("featureSummary works for XcmsExperiment", {
    res <- featureSummary(xmseg)
    expect_true(nrow(res) == nrow(featureDefinitions(xmseg)))
})

test_that("quantify,XcmsExperiment works", {
    expect_error(quantify(xmse), "No correspondence")
    res <- quantify(xmseg, method = "sum")
    expect_s4_class(res, "SummarizedExperiment")
    fd <- featureDefinitions(xmseg)
    expect_equal(rownames(fd), rownames(SummarizedExperiment::rowData(res)))
    a <- SummarizedExperiment::rowData(res)
    b <- as(fd[, colnames(fd) != "peakidx"], "DataFrame")
    expect_equal(a, b)
    expect_equal(SummarizedExperiment::assay(res),
                 featureValues(xmseg, method = "sum"))
    a <- SummarizedExperiment::colData(res)
    b <- sampleData(xmseg)
    rownames(b) <- colnames(SummarizedExperiment::assay(res))
    expect_equal(a, b)
})

test_that("addProcessHistory,XcmsExperiment works", {
    tmp <- xmse
    expect_error(addProcessHistory(tmp, "A"), "ProcessHistory")
    ph <- xcms:::ProcessHistory()
    tmp <- addProcessHistory(tmp, ph)
    expect_true(length(processHistory(tmp)) == 2L)
    expect_equal(processHistory(tmp)[[2L]], ph)
})

test_that("findChromPeaksIsolationWindow, etc, MsExperiment works", {

    ############################################################################
    ## findChromPeaksIsolationWindow
    cwp <- CentWaveParam(snthresh = 5, noise = 100, ppm = 10,
                         peakwidth = c(3, 20), prefilter = c(3, 1000))

    expect_error(
        findChromPeaksIsolationWindow(mse_dia, cwp, isolationWindow = 3),
        "Length")

    res <- findChromPeaksIsolationWindow(mse_dia, cwp)
    expect_s4_class(res, "XcmsExperiment")
    expect_true(all(chromPeakData(res)$ms_level == 2L))
    expect_true(length(processHistory(res)) == 1L)

    a <- findChromPeaks(mse_dia, cwp)
    a <- findChromPeaksIsolationWindow(a, cwp)
    expect_true(nrow(chromPeaks(a)) > nrow(chromPeaks(res)))
    expect_true(length(processHistory(a)) == 2L)

    ## Compare against XCMSnExp results
    expect_equal(nrow(chromPeaks(pest_swth)), nrow(chromPeaks(a)))
    expect_equal(chromPeaks(pest_swth), chromPeaks(a))

    ############################################################################
    ## reconstructChromPeakSpectra
    expect_error(reconstructChromPeakSpectra(a, peakId = c("a", "b")))

    ref <- reconstructChromPeakSpectra(pest_swth)
    res <- reconstructChromPeakSpectra(a)
    ref@processing <- ""
    res@processing <- ""
    expect_equal(ref, res)

    ## Change some settings.
    ref2 <- reconstructChromPeakSpectra(pest_swth, expandRt = 2, diffRt = 4,
                                        minCor = 0.9, intensity = "into")
    res2 <- reconstructChromPeakSpectra(a, expandRt = 2, diffRt = 4,
                                        minCor = 0.9, intensity = "into")
    ref2@processing <- ""
    res2@processing <- ""
    expect_equal(ref2, res2)

    expect_false(all(lengths(res2) == lengths(res)))

    ############################################################################
    ## filterIsolationWindow
    res <- filterIsolationWindow(mse_dia)
    expect_equal(length(spectra(res)), length(spectra(mse_dia)))
    expect_equal(rtime(spectra(res)), rtime(spectra(mse_dia)))

    ## with an isolation window.
    res <- filterIsolationWindow(mse_dia, mz = 301)
    expect_true(length(spectra(res)) < length(spectra(mse_dia)))
    expect_true(all(isolationWindowLowerMz(res@spectra) < 301))
    expect_true(all(isolationWindowUpperMz(res@spectra) > 301))
    expect_true(all(res@sampleDataLinks[["spectra"]][, 2L] %in%
                    seq_along(res@spectra)))

    ## on an XcmsExperiment.
    res <- filterIsolationWindow(a)
    expect_equal(chromPeaks(res), chromPeaks(a))

    res <- filterIsolationWindow(a, mz = 301)
    expect_true(nrow(chromPeaks(res)) < nrow(chromPeaks(a)))
    expect_true(all(chromPeakData(res)$isolationWindowLowerMz < 301))
    expect_true(all(chromPeakData(res)$isolationWindowUpperMz > 301))
})

test_that("chromPeaksChromatograms,XcmsExperiment works", {
    expect_error(chromPeakChromatograms(xmse, peaks = 1:3), "expected to")

    chrs <- chromPeakChromatograms(xmse)

    ## Test providing peaks. Only those from one file.
    pks <- rownames(chromPeaks(xmse)[chromPeaks(xmse)[, "sample"] == 2, ])
    res <- chromPeakChromatograms(xmse, peaks = pks)
    expect_equal(rownames(res), pks)
    ref <- chrs[rownames(chrs) %in% pks, 1L]
    expect_equal(ref@.Data, res@.Data)
    expect_equal(fData(ref), fData(res))
    expect_equal(chromPeaks(ref), chromPeaks(res))

    ## Test providing peaks. different order.
    pks <- sample(rownames(chromPeaks(xmseg)), 10)
    res <- chromPeakChromatograms(xmseg, peaks = pks)
    ref <- chrs[match(rownames(res), rownames(chrs)), 1L]
    expect_equal(res@.Data, ref@.Data)
    expect_equal(fData(res), fData(ref))
    expect_equal(chromPeaks(res), chromPeaks(ref))

    ## Test on a SWATH data set: are MS1 and MS2 chrom peaks extracted
    ## correctly?
    cwp <- CentWaveParam(snthresh = 5, noise = 100, ppm = 10,
                         peakwidth = c(3, 20), prefilter = c(3, 1000))
    xmse_dia <- findChromPeaks(mse_dia, param = cwp)
    xmse_dia <- findChromPeaksIsolationWindow(xmse_dia, param = cwp)
    res <- chromPeakChromatograms(xmse_dia)
    ## To compare against what?
    ints <- vapply(res, function(z) sum(intensity(z), na.rm = TRUE), numeric(1))
    expect_true(cor(chromPeaks(res)[, "into"], ints) >= 0.97)
})

test_that("setAs,XcmsExperiment,xcmsSet works", {
    res <- as(xmseg, "xcmsSet")
    expect_s4_class(res, "xcmsSet")
    expect_equal(peaks(res), chromPeaks(xmseg))
})

test_that(".index_chrom_peaks works", {
    ## xmse
    res <- .index_chrom_peaks(xmse)
    expect_equal(res, seq_len(nrow(chromPeaks(xmse))))

    ## MS level
    res <- .index_chrom_peaks(xmse, msLevel = c(3, 4))
    expect_equal(res, integer())

    ## rt
    res <- .index_chrom_peaks(xmse, rt = c(2500, 2700),
                              type = "apex_within")
    rts <- chromPeaks(xmse)[res, "rt"]
    expect_true(all(rts >= 2500 & rts <= 2700))

    ## mz
    res <- .index_chrom_peaks(xmse, mz = c(400, 600), type = "apex_within")
    mzs <- chromPeaks(xmse)[res, "mz"]
    expect_true(all(mzs >= 400 & mzs <= 600))

    ## rt and mz
    res <- .index_chrom_peaks(xmse, mz = c(400, 600), type = "apex_within",
                              rt = c(2500, 2700))
    pks <- chromPeaks(xmse)[res, c("mz", "rt")]
    expect_true(all(pks[, "mz"] >= 400 & pks[, "mz"] <= 600 &
                    pks[, "rt"] >= 2500 & pks[, "rt"] <= 2700))
})
