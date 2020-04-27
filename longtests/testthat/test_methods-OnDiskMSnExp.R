test_that("profMat,OnDiskMSnExp works", {
    ## Get it from all 3 files in one go.
    res <- profMat(faahko_od, step = 2)
    res_2 <- profMat(xcmsRaw(faahko_3_files[2], profstep = 0), step = 2)
    expect_equal(res_2, res[[2]])
    res_2 <- profMat(xcmsRaw(faahko_3_files[3], profstep = 0), step = 2)
    expect_equal(res_2, res[[3]])
    res_2 <- profMat(faahko_xod, step = 2)
    expect_equal(res, res_2)
    res <- profMat(faahko_od, step = 2, method = "binlin", fileIndex = 2)
    res_2 <- profMat(xcmsRaw(faahko_3_files[2], profstep = 0), step = 2,
                     method = "binlin")
    expect_equal(res_2, res[[1]])

    ## Simulating issue #312
    od_1 <- filterFile(microtofq_od, 1)
    od_1_clnd <- clean(removePeaks(od_1, t = 1800))
    res_clnd <- profMat(od_1_clnd)
})

test_that("findChromPeaks,OnDiskMSnExp,CentWaveParam works", {
    fs <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xr <- deepCopy(faahko_xr_1)
    onDisk <- filterFile(faahko_od, file = 1)
    ppm <- 40
    snthresh <- 40
    res_x <- findPeaks.centWave(xr, ppm = ppm, snthresh = snthresh,
                                noise = 100000)@.Data
    ## Bypass xcmsRaw
    xs <- xcmsSet(fs[1], profparam = list(profstep = 0), ppm = ppm,
                  snthresh = snthresh, method = "centWave",
                  noise = 100000)
    expect_equal(xs@peaks[, colnames(res_x)], res_x)
    ## OnDiskMSnExp
    ## onDisk <- readMSData(fs[1], msLevel. = 1, mode = "onDisk")
    cwp <- CentWaveParam(ppm = ppm, snthresh = snthresh, noise = 100000,
                         prefilter = c(3, 10000))
    res <- findChromPeaks(onDisk, param = cwp, return.type = "list")
    expect_equal(res[[1]], peaks(xs)@.Data)

    expect_error(findChromPeaks(onDisk, param = cwp, msLevel = 2))

    ## returning an xcmsSet
    res <- findChromPeaks(onDisk, param = cwp, return.type = "xcmsSet")
    pks <- peaks(res)
    rownames(pks) <- NULL
    expect_equal(pks[, colnames(peaks(xs))], peaks(xs))
    expect_true(is(res, "xcmsSet"))

    ## Return type XCMSnExp
    res <- findChromPeaks(onDisk, param = cwp)
    expect_true(hasChromPeaks(res))
    expect_true(!hasAdjustedRtime(res))
    expect_true(!hasFeatures(res))
    pks <- chromPeaks(res)
    rownames(pks) <- NULL
    expect_equal(peaks(xs)@.Data, pks[, !colnames(pks) %in% c("is_filled", "ms_level")])

    ## check that rownames are set
    expect_true(!is.null(rownames(chromPeaks(res))))
    expect_true(length(grep("CP", rownames(chromPeaks(res)))) ==
                nrow(chromPeaks(res)))
})

test_that("findChromPeaks,OnDiskMSnExp,CentWavePredIsoParam works", {
    fs <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xr <- deepCopy(faahko_xr_1)
    snth <- 20
    ns <- 2500
    snthIso <- 5
    res_x <- findPeaks.centWaveWithPredictedIsotopeROIs(xr, noise = ns,
                                                        snthresh = snth,
                                                        snthreshIsoROIs = snthIso)@.Data
    ## Bypass xcmsRaw
    xs <- xcmsSet(fs[1], profparam = list(profstep = 0), snthresh = snth,
                  method = "centWaveWithPredictedIsotopeROIs", noise = ns,
                  snthreshIsoROIs = snthIso)
    expect_equal(xs@peaks[, colnames(res_x)], res_x)
    ## OnDiskMSnExp
    onDisk <- readMSData(fs[1], msLevel. = 1, mode = "onDisk")
    cwp <- CentWavePredIsoParam(snthresh = snth, noise = ns,
                                snthreshIsoROIs = snthIso)
    res <- findChromPeaks(onDisk, param = cwp, return.type = "list")
    expect_equal(res[[1]], peaks(xs)@.Data)
    expect_error(findChromPeaks(onDisk, param = cwp, msLevel = 2))

    ## returning an xcmsSet
    res <- findChromPeaks(onDisk, param = cwp, return.type = "xcmsSet")
    pks <- peaks(res)
    rownames(pks) <- NULL
    expect_equal(pks[, colnames(peaks(xs))], peaks(xs))
    expect_true(is(res, "xcmsSet"))

    ## Return an XCMSnExp
    res <- findChromPeaks(onDisk, param = cwp)
    expect_true(hasChromPeaks(res))
    expect_true(!hasAdjustedRtime(res))
    expect_true(!hasFeatures(res))
    pks <- chromPeaks(res)
    rownames(pks) <- NULL
    expect_equal(peaks(xs)@.Data, pks[, colnames(peaks(xs)@.Data)])
})

test_that("findChromPeaks,OnDiskMSnExp,MassifquantParam works", {
    mzf <- system.file("microtofq/MM14.mzML", package = "msdata")
    mqp <- MassifquantParam(ppm = 20, criticalValue = 1.2)
    res <- xcmsSet(mzf[1], method = "massifquant", ppm = 20,
                   criticalValue = 1.2)
    ## onDisk
    onDisk <- readMSData(mzf[1], mode = "onDisk")
    res_o <- findChromPeaks(onDisk, param = mqp, return.type = "xcmsSet")
    expect_equal(unname(peaks(res_o)[, colnames(peaks(res))]),
                 unname(peaks(res)))
    expect_equal(unname(res_o@rt$raw[[1]]), unname(res@rt$raw[[1]]))

    expect_error(findChromPeaks(onDisk, param = mqp, msLevel = 2))
})

test_that("findChromPeaks,OnDiskMSnExp,MatchedFilterParam works", {
    fs <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    mfp <- MatchedFilterParam(binSize = 20, impute = "lin")
    res <- xcmsSet(fs[1], method = "matchedFilter", profmethod = "binlin",
                   step = binSize(mfp))
    ## onDisk
    ## onDisk <- readMSData(fs[1], mode = "onDisk")
    onDisk <- filterFile(faahko_od, file = 1)
    res_o <- findChromPeaks(onDisk, param = mfp, return.type = "xcmsSet")
    expect_equal(unname(peaks(res_o)[, colnames(peaks(res))]),
                 unname(peaks(res)))
    expect_equal(unname(res_o@rt$raw[[1]]), unname(res@rt$raw[[1]]))

    expect_error(findChromPeaks(onDisk, param = mfp, msLevel = 2))
})
