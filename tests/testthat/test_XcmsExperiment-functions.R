test_that(".mse_add_chrom_peaks works", {
    a <- new("XcmsExperiment")
    pks <- chromPeaks(faahko_xod)
    pkd <- as.data.frame(chromPeakData(faahko_xod))

    res <- .mse_add_chrom_peaks(a, pks, pkd)
    expect_equal(res@chromPeakData, pkd)
    expect_equal(res@chromPeaks, pks)

    pkd$add_col <- "a"
    res <- .mse_add_chrom_peaks(res, pks, pkd)
    expect_equal(rownames(res@chromPeaks), rownames(res@chromPeakData))
    expect_equal(rownames(res@chromPeaks),
                 .featureIDs(nrow(res@chromPeaks), prefix = "CP"))
    expect_equal(colnames(res@chromPeakData),
                 c("ms_level", "is_filled", "add_col"))
})

test_that(".mse_valid_chrom_peaks works", {
    a <- cbind(mz = 1:3, mzmin = 1:3)
    expect_match(.mse_valid_chrom_peaks(a), "missing")
    expect_equal(.mse_valid_chrom_peaks(chromPeaks(faahko_xod)), NULL)
})

test_that(".mse_valid_chrom_peak_data works", {
    expect_equal(.mse_valid_chrom_peak_data(chromPeakData(faahko_xod)), NULL)
    a <- data.frame(msLevel = 1:3, other = "a")
    expect_match(.mse_valid_chrom_peak_data(a), "missing")
})

test_that(".mse_valid_feature_def", {
    expect_equal(.mse_valid_feature_def(.empty_feature_definitions()), NULL)
    expect_match(.mse_valid_feature_def(3), "required columns")
    expect_match(.mse_valid_feature_def(data.frame(a = 1:4, b = 1:4)),
                 "required columns")
})

test_that(".mse_same_rownames works", {
    expect_equal(.mse_same_rownames(chromPeaks(faahko_xod),
                                    chromPeakData(faahko_xod)), NULL)
    a <- cbind(a = 1:3, b = 1:3)
    rownames(a) <- c("a", "b", "c")
    expect_match(.mse_same_rownames(a, data.frame(ms_level = 1:3)),
                 "don't match")
})

test_that(".aggregate_intensities works", {
    pks <- cbind(mz = c(12.3, 12.4, 12.5, 13, 13.2, 13.3, 13.4, 13.5),
                 intensity = 1:8)
    res <- .aggregate_intensities(pks)
    expect_equal(res, sum(1:8))

    res <- .aggregate_intensities(pks, mzr = c(12.4, 13.2))
    expect_equal(res, sum(2:5))

    res <- .aggregate_intensities(pks, mzr = c(12.4, 13.2), INTFUN = max)
    expect_equal(res, 5)
})

test_that(".history2fill_fun works", {
    expect_equal(.history2fill_fun(), .chrom_peak_intensity_centWave)
    h <- xcms:::XProcessHistory()
    expect_equal(.history2fill_fun(list(h)), .chrom_peak_intensity_centWave)
    h@param <- MSWParam()
    expect_equal(.history2fill_fun(list(h)), .chrom_peak_intensity_msw)
    h@param <- MatchedFilterParam()
    expect_equal(.history2fill_fun(list(h)),
                 .chrom_peak_intensity_matchedFilter)
})

test_that(".pmat_filter_mz works", {
    a <- cbind(mz = 1:4, intensity = c(1.2, 3.4, 5.6, 8.9))
    res <- .pmat_filter_mz(a, c(5, 6))
    expect_true(nrow(res) == 0)
    expect_equal(colnames(res), colnames(a))

    res <- .pmat_filter_mz(a)
    expect_equal(a, res)

    res <- .pmat_filter_mz(a, c(2, 3))
    expect_equal(res, a[2:3, ])
})

test_that(".chrom_peak_intensity_msw  works", {
    library(MsExperiment)
    fticrf <- list.files(system.file("fticr-mzML", package = "msdata"),
                         recursive = TRUE, full.names = TRUE)
    fls <- normalizePath(fticrf)[1:2]
    mp <- MSWParam(scales = c(1, 7), peakThr = 80000, ampTh = 0.005,
                   SNR.method = "data.mean", winSize.noise = 500)

    tmp <- MsExperiment()
    df <- data.frame(mzML_file = basename(fls),
                     dataOrigin = fls, sample = c("a", "b"))
    sampleData(tmp) <- DataFrame(df)
    spectra(tmp) <- Spectra::Spectra(fls)
    tmp <- linkSampleData(
        tmp, with = "sampleData.dataOrigin = spectra.dataOrigin")
    tmp <- findChromPeaks(tmp, mp)

    pks <- chromPeaks(tmp[1L])
    cn <- colnames(pks)
    colnames(pks)[cn == "mz"] <- "mzmed"
    p <- Spectra::peaksData(spectra(tmp[1L]))
    rt <- rtime(spectra(tmp[1L]))
    res <- .chrom_peak_intensity_msw(p, rt, pks, sampleIndex = 1L, cn = cn)
    expect_equal(res[, "mz"], pks[, "mzmed"])
    ## expect_equal(res[, "into"], pks[, "into"]) # not the same...
    expect_equal(res[, "maxo"], pks[, "maxo"])

    ## also check that it works on `XcmsExperiment`.
    mzc <- MzClustParam(sampleGroups = c(1, 1))
    tmp <- groupChromPeaks(tmp, param = mzc)
    cpp <- ChromPeakAreaParam()
    res <- fillChromPeaks(tmp, cpp)
    expect_true(length(res@processHistory) > length(tmp@processHistory))
    expect_true(sum(is.na(featureValues(res))) < sum(is.na(featureValues(tmp))))

    ## Compare with XCMSnExp.
    ## ref <- readMSData(fticrf[1:2], msLevel. = 1, mode = "onDisk")
    ## ref <- findChromPeaks(ref, mp)
    ## ref <- groupChromPeaks(ref, param = mzc)
    ## ref <- fillChromPeaks(ref, cpp)
    ## expect_equal(chromPeaks(ref), chromPeaks(res))
})

test_that(".index_chrom_peaks works", {
    ## xod_xg
    ## xmse
    res <- .index_chrom_peaks(xmse)
    expect_equal(res, seq_len(nrow(chromPeaks(xmse))))
    res2 <- .index_chrom_peaks(xod_x)
    expect_equal(res, res2)

    ## MS level
    res <- .index_chrom_peaks(xmse, msLevel = c(3, 4))
    expect_equal(res, integer())
    res2 <- .index_chrom_peaks(xod_x, msLevel = c(3, 4))
    expect_equal(res, res2)

    ## rt
    res <- .index_chrom_peaks(xmse, rt = c(2500, 2700),
                              type = "apex_within")
    rts <- chromPeaks(xmse)[res, "rt"]
    expect_true(all(rts >= 2500 & rts <= 2700))
    expect_equal(res, .index_chrom_peaks(xod_x, rt = c(2500, 2700),
                                         type = "apex_within"))

    ## mz
    res <- .index_chrom_peaks(xmse, mz = c(400, 600), type = "apex_within")
    mzs <- chromPeaks(xmse)[res, "mz"]
    expect_true(all(mzs >= 400 & mzs <= 600))
    expect_equal(res, .index_chrom_peaks(xod_x, mz = c(400, 600),
                                         type = "apex_within"))

    ## rt and mz
    res <- .index_chrom_peaks(xmse, mz = c(400, 600), type = "apex_within",
                              rt = c(2500, 2700))
    pks <- chromPeaks(xmse)[res, c("mz", "rt")]
    expect_true(all(pks[, "mz"] >= 400 & pks[, "mz"] <= 600 &
                    pks[, "rt"] >= 2500 & pks[, "rt"] <= 2700))
    expect_equal(res, .index_chrom_peaks(xod_x, mz = c(400, 600),
                                         type = "apex_within",
                                         rt = c(2500, 2700)))
})

test_that(".chromPeaks works", {
    x <- new("XcmsExperiment")
    res <- .chromPeaks(x)
    expect_equal(res, x@chromPeaks)

    res <- .chromPeaks(xmse)
    expect_equal(res, xmse@chromPeaks)
    res <- .chromPeaks(xod_x)
    expect_equal(res, chromPeaks(xod_x@msFeatureData))
})

test_that(".chromPeakData works", {
    x <- new("XcmsExperiment")
    res <- .chromPeakData(x)
    expect_equal(res, x@chromPeakData)

    res <- .chromPeakData(xmse)
    expect_equal(res, xmse@chromPeakData)
    res <- .chromPeakData(xmse, msLevel = 2L)
    expect_equal(res, xmse@chromPeakData[integer(), ])
})

test_that(".features_ms_region works", {
    res <- .features_ms_region(
        xod_xgrg, features = rownames(featureDefinitions(xod_xgrg)))
    expect_equal(nrow(res), nrow(featureDefinitions(xod_xgrg)))
    expect_equal(colnames(res), c("mzmin", "mzmax", "rtmin", "rtmax"))
    expect_true(all(res[, "mzmin"] <= res[, "mzmax"]))
    expect_true(all(res[, "rtmin"] < res[, "rtmax"]))

    expect_error(.features_ms_region(xod_xgrg,
                                     features = c("a", "b")), "out of")

    res <- .features_ms_region(
        xmseg, features = rownames(featureDefinitions(xmseg)))
    expect_equal(rownames(res), rownames(featureDefinitions(xmseg)))
})
