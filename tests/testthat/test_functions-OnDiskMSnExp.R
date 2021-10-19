test_that(".obiwarp works", {
    skip_on_os(os = "windows", arch = "i386")


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
    skip_on_os(os = "windows", arch = "i386")

    od1 <- filterFile(od_x, 1)
    od2 <- filterFile(od_x, 2:3)
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

test_that(".split_by_file and .split_bu_file2 work", {
    skip_on_os(os = "windows", arch = "i386")

    a <- lapply(seq_along(fileNames(od_x)), filterFile, object = od_x)
    b <- .split_by_file(od_x, subsetFeatureData = FALSE)
    expect_equal(a[[2]][[19]], b[[2]][[19]])
    expect_equal(spectra(a[[3]]), spectra(b[[3]]))
    d <- .split_by_file2(od_x, subsetFeatureData = FALSE)
    expect_equal(b, d)

    b_2 <- .split_by_file(od_x, subsetFeatureData = TRUE)
    expect_true(ncol(fData(b_2[[1]])) < ncol(fData(b[[1]])))
    d_2 <- .split_by_file2(od_x, subsetFeatureData = TRUE)

    res <- .split_by_file(xod_xgr)
    expect_true(length(res) == 3)
    expect_equal(rtime(res[[2]]),
                 rtime(xod_xgr, bySample = TRUE, adjusted = TRUE)[[2]])
    expect_true(ncol(fData(res[[1]])) < ncol(fData(xod_xgr)))
    res_2 <- .split_by_file(xod_xgr, subsetFeatureData = TRUE)
    expect_equal(res, res_2)
    expect_error(.split_by_file(xod_xgr, msLevel. = 2), "No MS level")

    a <- filterFile(xod_xgrg, 2, keepAdjustedRtime = TRUE)
    b <- .split_by_file(xod_xgrg, to_class = "XCMSnExp")
    expect_equal(rtime(a), rtime(b[[2]]))
    expect_true(is(b[[1]], "XCMSnExp"))
    expect_equal(chromPeaks(a), chromPeaks(b[[2]]))
    expect_equal(chromPeakData(a), chromPeakData(b[[2]]))
    d <- .split_by_file2(xod_xgrg, to_class = "XCMSnExp",
                         subsetFeatureData = TRUE)
    expect_equal(b, d)
})

test_that(".estimate_prec_intensity works", {
    skip_on_os(os = "windows", arch = "i386")

    tmp <- filterRt(as(pest_dda, "OnDiskMSnExp"), c(300, 400))
    res <- .estimate_prec_intensity(tmp, method = "previous")
    expect_true(length(res) == length(tmp))
    expect_true(all(is.na(res[msLevel(tmp) == 1L])))

    res <- .estimate_prec_intensity(tmp, method = "interpolation")
    expect_true(length(res) == length(tmp))
    expect_true(all(is.na(res[msLevel(tmp) == 1L])))
})

test_that("estimatePrecursorIntensity,OnDiskMSnExp works", {
    skip_on_os(os = "windows", arch = "i386")

    tmp <- filterRt(pest_dda, c(300, 400))
    res <- estimatePrecursorIntensity(tmp, BPPARAM = SerialParam())
    expect_true(length(res) == length(tmp))
    expect_true(all(is.na(res[msLevel(tmp) == 1L])))

    expect_error(estimatePrecursorIntensity(od_x))
})

test_that(".OnDiskMSnExp2MsBackendMzR works", {
    skip_on_os(os = "windows", arch = "i386")

    if (requireNamespace("Spectra", quietly = TRUE)) {
        res <- .OnDiskMSnExp2MsBackendMzR(xod_x)
        expect_equal(Spectra::rtime(res), unname(rtime(xod_x)))
    }
})

test_that(".fData2MsBackendMzR works", {
    skip_on_os(os = "windows", arch = "i386")

    if (requireNamespace("Spectra", quietly = TRUE)) {
        res <- .fData2MsBackendMzR(fData(xod_x), fileNames(xod_x))
        expect_true(is(res, "MsBackendMzR"))
        expect_equal(Spectra::rtime(res), unname(rtime(xod_x)))
        tmp <- xod_x[c(1, 45, 113)]
        res <- .fData2MsBackendMzR(fData(tmp), fileNames(tmp))
        expect_equal(Spectra::rtime(res), unname(rtime(tmp)))
    }
})
