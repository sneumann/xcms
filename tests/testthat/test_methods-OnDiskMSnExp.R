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
})

test_that(
    "findChromPeaks,OnDiskMSnExp,CentWaveParam works with multiple MS levels", {
        msn_file <- system.file(
            package = "msdata", 
            "proteomics/MS3TMT10_01022016_32917-33481.mzML.gz")
        msn_file <- system.file(
            package = "msdata",
            "proteomics/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzML.gz")
        msn_data <- readMSData(msn_file, mode = "onDisk")
        msn_xdata <- findChromPeaks(pickPeaks(msn_data), param = CentWaveParam())
        expect_equal(msLevel(msn_data), msLevel(msn_xdata))
    })

test_that("findChromPeaks,OnDiskMSnExp,CentWaveParam variants", {
    ## Reproduce with msdata files:
    fl <- system.file("microtofq/MM14.mzML", package = "msdata")
    raw <- readMSData(fl, mode = "onDisk")
    options(originalCentWave = TRUE)
    tmp <- findChromPeaks(raw, param = CentWaveParam(peakwidth = c(2, 10)))
    ## ## Use the getPeakInt2 which uses the rawMat function.
    ## pkI2 <- xcms:::.getPeakInt2(tmp, chromPeaks(tmp))
    ## ## Use the getPeakInt3 which uses the getEIC C function.
    ## pkI3 <- xcms:::.getPeakInt3(tmp, chromPeaks(tmp))
    ## ## These fail for the original centWave code.
    ## expect_true(sum(pkI2 != chromPeaks(tmp)[, "into"]) > length(pkI2) / 2)
    ## ## expect_equal(unname(pkI2), unname(chromPeaks(tmp)[, "into"]))
    ## ## expect_equal(unname(pkI3), unname(chromPeaks(tmp)[, "into"]))
    ## expect_equal(pkI2, pkI3)
    ## Try with new implementation.
    options(originalCentWave = FALSE)
    tmp2 <- findChromPeaks(raw, param = CentWaveParam(peakwidth = c(2, 10)))
    ## Find different number of peaks:
    expect_true(nrow(chromPeaks(tmp2)) != nrow(chromPeaks(tmp)))
    ## Are the peaks similar?
    id_1 <- paste(chromPeaks(tmp)[, "mz"], chromPeaks(tmp)[, "rt"])
    id_2 <- paste(chromPeaks(tmp2)[, "mz"], chromPeaks(tmp2)[, "rt"])
    ## But all of the ones from the old are ALSO in the new one.
    expect_true(all(id_1 %in% id_2))
    ## Are the peaks the same?
    cp2 <- chromPeaks(tmp2)[id_2 %in% id_1, ]
    cn <- colnames(cp2)
    cn <- cn[!(cn %in% c("intb", "into", "rtmin", "rtmax"))]
    expect_equal(cp2[, cn], chromPeaks(tmp)[, cn])
    ## Are the values related?
    plot(cp2[, "into"], chromPeaks(tmp)[, "into"])   ## Very similar
    plot(cp2[, "intb"], chromPeaks(tmp)[, "intb"])   ## Very similar
    plot(cp2[, "rtmin"], chromPeaks(tmp)[, "rtmin"])   ## Very similar
    plot(cp2[, "rtmax"], chromPeaks(tmp)[, "rtmax"])   ## Very similar
    ## Use the getPeakInt3 which uses the getEIC C function.
    ## pkI2_2 <- xcms:::.getPeakInt2(tmp2, chromPeaks(tmp2))
    ## pkI3_2 <- xcms:::.getPeakInt3(tmp2, chromPeaks(tmp2))
    ## ## These fail for the original centWave code.
    ## expect_equal(unname(pkI2_2), unname(chromPeaks(tmp2)[, "into"]))
    ## expect_equal(unname(pkI3_2), unname(chromPeaks(tmp2)[, "into"]))
    ## expect_equal(pkI2_2, pkI3_2)

    
    ## The same for one of the test files; this works even with the original
    ## centWave code
    options(originalCentWave = TRUE)
    tmp <- filterFile(xod_xgrg, file = 3)
    ## ## Use the getPeakInt2 which uses the rawMat function.
    ## pkI2 <- xcms:::.getPeakInt2(tmp, chromPeaks(tmp))
    ## ## Use the getPeakInt3 which uses the getEIC C function.
    ## pkI3 <- xcms:::.getPeakInt3(tmp, chromPeaks(tmp))
    ## expect_equal(pkI2, pkI3)
    ## expect_equal(unname(pkI2), unname(chromPeaks(tmp)[, "into"]))
    ## expect_equal(unname(pkI3), unname(chromPeaks(tmp)[, "into"]))
    ## New modified centWave.
    options(originalCentWave = FALSE)
    tmp2 <- findChromPeaks(filterFile(faahko_od, file = 3),
                           CentWaveParam(noise = 10000, snthresh = 40))
    ## Even the identified peaks are identical!
    expect_equal(chromPeaks(tmp), chromPeaks(tmp2))
    ## Use the getPeakInt2 which uses the rawMat function.
    ## pkI2 <- xcms:::.getPeakInt2(tmp2, chromPeaks(tmp2))
    ## ## Use the getPeakInt3 which uses the getEIC C function.
    ## pkI3 <- xcms:::.getPeakInt3(tmp2, chromPeaks(tmp2))
    ## expect_equal(pkI2, pkI3)
    ## expect_equal(unname(pkI2), unname(chromPeaks(tmp2)[, "into"]))
    ## expect_equal(unname(pkI3), unname(chromPeaks(tmp2)[, "into"]))
    options(originalCentWave = TRUE)
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
    cwp <- CentWaveParam(ppm = ppm, snthresh = snthresh, noise = 100000)
    res <- findChromPeaks(onDisk, param = cwp, return.type = "list")
    expect_equal(res[[1]], peaks(xs)@.Data)

    expect_error(findChromPeaks(onDisk, param = cwp, msLevel = 2))

    ## returning an xcmsSet
    res <- findChromPeaks(onDisk, param = cwp, return.type = "xcmsSet")
    expect_equal(peaks(res)[, colnames(peaks(xs))], peaks(xs))
    expect_true(is(res, "xcmsSet"))

    ## Return type XCMSnExp
    res <- findChromPeaks(onDisk, param = cwp)
    expect_true(hasChromPeaks(res))
    expect_true(!hasAdjustedRtime(res))
    expect_true(!hasFeatures(res))
    expect_equal(peaks(xs)@.Data, chromPeaks(res)[, -ncol(chromPeaks(res))])
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
    expect_equal(peaks(res)[, colnames(peaks(xs))], peaks(xs))
    expect_true(is(res, "xcmsSet"))

    ## Return an XCMSnExp
    res <- findChromPeaks(onDisk, param = cwp)
    expect_true(hasChromPeaks(res))
    expect_true(!hasAdjustedRtime(res))
    expect_true(!hasFeatures(res))
    expect_equal(peaks(xs)@.Data, chromPeaks(res)[, colnames(peaks(xs)@.Data)])
})

test_that("findChromPeaks,OnDiskMSnExp,MassifquantParam works", {
    mzf <- system.file("microtofq/MM14.mzML", package = "msdata")
    mqp <- MassifquantParam(ppm = 20, criticalValue = 1.2)
    res <- xcmsSet(mzf[1], method = "massifquant", ppm = 20,
                   criticalValue = 1.2)
    ## onDisk
    onDisk <- readMSData(mzf[1], mode = "onDisk")
    res_o <- findChromPeaks(onDisk, param = mqp, return.type = "xcmsSet")
    expect_equal(peaks(res_o)[, colnames(peaks(res))], peaks(res))
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
    expect_equal(peaks(res_o)[, colnames(peaks(res))], peaks(res))
    expect_equal(unname(res_o@rt$raw[[1]]), unname(res@rt$raw[[1]]))

    expect_error(findChromPeaks(onDisk, param = mfp, msLevel = 2))
})
