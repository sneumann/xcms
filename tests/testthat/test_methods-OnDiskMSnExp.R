test_that("profMat,OnDiskMSnExp works", {
    skip_on_os(os = "windows", arch = "i386")

    ## Get it from all 3 files in one go.
    res <- profMat(filterRt(faahko_od, c(2500, 3000)), step = 2)
    res_2 <- profMat(filterRt(faahko_xod, c(2500, 3000)), step = 2)
    expect_equal(res, res_2)

    ## Simulating issue #312
    od_1 <- filterFile(microtofq_od, 1)
    od_1_clnd <- clean(removePeaks(od_1, t = 1800))
    res_clnd <- profMat(od_1_clnd)
})

test_that("findChromPeaks,OnDiskMSnExp,CentWaveParam variants", {
    skip_on_os(os = "windows", arch = "i386")

    ## Reproduce with msdata files:
    fl <- system.file("microtofq/MM14.mzML", package = "msdata")
    raw <- readMSData(fl, mode = "onDisk")
    options(originalCentWave = TRUE)
    tmp <- findChromPeaks(raw, param = CentWaveParam(peakwidth = c(2, 10),
                                                     prefilter = c(3, 500)))
    ## ## Use the getPeakInt2 which uses the rawMat function.
    ## pkI2 <- .getPeakInt2(tmp, chromPeaks(tmp))
    ## ## Use the getPeakInt3 which uses the getEIC C function.
    ## pkI3 <- .getPeakInt3(tmp, chromPeaks(tmp))
    ## ## These fail for the original centWave code.
    ## expect_true(sum(pkI2 != chromPeaks(tmp)[, "into"]) > length(pkI2) / 2)
    ## ## expect_equal(unname(pkI2), unname(chromPeaks(tmp)[, "into"]))
    ## ## expect_equal(unname(pkI3), unname(chromPeaks(tmp)[, "into"]))
    ## expect_equal(pkI2, pkI3)
    ## Try with new implementation.
    options(originalCentWave = FALSE)
    tmp2 <- findChromPeaks(raw, param = CentWaveParam(peakwidth = c(2, 10),
                                                      prefilter = c(3, 500)))
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
    pks <- chromPeaks(tmp)
    rownames(pks) <- NULL
    rownames(cp2) <- NULL
    expect_equal(cp2[, cn], pks[, cn])
    ## Are the values related?
    ## plot(cp2[, "into"], pks[, "into"])   ## Very similar
    ## plot(cp2[, "intb"], pks[, "intb"])   ## Very similar
    ## plot(cp2[, "rtmin"], pks[, "rtmin"])   ## Very similar
    ## plot(cp2[, "rtmax"], pks[, "rtmax"])   ## Very similar
    ## Use the getPeakInt3 which uses the getEIC C function.
    ## pkI2_2 <- .getPeakInt2(tmp2, chromPeaks(tmp2))
    ## pkI3_2 <- .getPeakInt3(tmp2, chromPeaks(tmp2))
    ## ## These fail for the original centWave code.
    ## expect_equal(unname(pkI2_2), unname(chromPeaks(tmp2)[, "into"]))
    ## expect_equal(unname(pkI3_2), unname(chromPeaks(tmp2)[, "into"]))
    ## expect_equal(pkI2_2, pkI3_2)


    ## The same for one of the test files; this works even with the original
    ## centWave code
    options(originalCentWave = TRUE)
    tmp <- filterFile(xod_xgrg, file = 3, keepAdjustedRtime = FALSE)
    ## ## Use the getPeakInt2 which uses the rawMat function.
    ## pkI2 <- .getPeakInt2(tmp, chromPeaks(tmp))
    ## ## Use the getPeakInt3 which uses the getEIC C function.
    ## pkI3 <- .getPeakInt3(tmp, chromPeaks(tmp))
    ## expect_equal(pkI2, pkI3)
    ## expect_equal(unname(pkI2), unname(chromPeaks(tmp)[, "into"]))
    ## expect_equal(unname(pkI3), unname(chromPeaks(tmp)[, "into"]))
    ## New modified centWave.
    options(originalCentWave = FALSE)
    tmp2 <- findChromPeaks(filterFile(faahko_od, file = 3),
                           CentWaveParam(noise = 10000, snthresh = 40,
                                         prefilter = c(3, 10000)))
    ## Even the identified peaks are identical!
    expect_equal(unname(chromPeaks(tmp)), unname(chromPeaks(tmp2)))
    ## Use the getPeakInt2 which uses the rawMat function.
    ## pkI2 <- .getPeakInt2(tmp2, chromPeaks(tmp2))
    ## ## Use the getPeakInt3 which uses the getEIC C function.
    ## pkI3 <- .getPeakInt3(tmp2, chromPeaks(tmp2))
    ## expect_equal(pkI2, pkI3)
    ## expect_equal(unname(pkI2), unname(chromPeaks(tmp2)[, "into"]))
    ## expect_equal(unname(pkI3), unname(chromPeaks(tmp2)[, "into"]))
    options(originalCentWave = TRUE)
})

test_that("findChromPeaks,OnDiskMSnExp,CentWaveParam works", {
    skip_on_os(os = "windows", arch = "i386")

    onDisk <- filterFile(faahko_od, file = 1)
    ppm <- 40
    snthresh <- 40
    cwp <- CentWaveParam(ppm = ppm, snthresh = snthresh, noise = 100000,
                         prefilter = c(3, 10000))
    res <- findChromPeaks(onDisk, param = cwp)
    expect_true(hasChromPeaks(res))
    expect_equal(nrow(chromPeaks(res)), 6)

    expect_error(findChromPeaks(onDisk, param = cwp, msLevel = 2))

    pks <- chromPeaks(res)

    ## check that rownames are set
    expect_true(!is.null(rownames(chromPeaks(res))))
    expect_true(length(grep("CP", rownames(chromPeaks(res)))) ==
                nrow(chromPeaks(res)))
})

test_that("findChromPeaks,OnDiskMSnExp,CentWavePredIsoParam works", {
    skip_on_os(os = "windows", arch = "i386")

    ## OnDiskMSnExp
    onDisk <- filterFile(faahko_od, file = 1)
    cwp <- CentWavePredIsoParam(snthresh = 20, noise = 2500,
                                snthreshIsoROIs = 5, prefilter = c(5, 10000))
    expect_error(findChromPeaks(onDisk, param = cwp, msLevel = 2))

    ## Return an XCMSnExp
    res <- findChromPeaks(onDisk, param = cwp)
    expect_true(hasChromPeaks(res))
    expect_true(!hasAdjustedRtime(res))
    expect_true(!hasFeatures(res))
    pks <- chromPeaks(res)
})

test_that("findChromPeaks,OnDiskMSnExp,MassifquantParam works", {
    skip_on_os(os = "windows", arch = "i386")

    onDisk <- filterFile(microtofq_od, 1)
    res_o <- findChromPeaks(onDisk, param = MassifquantParam(prefilter = c(5, 5000)))
    expect_true(hasChromPeaks(res_o))
    expect_equal(nrow(chromPeaks(res_o)), 15)

    expect_error(findChromPeaks(onDisk, param = mqp, msLevel = 2))
})

test_that("findChromPeaks,OnDiskMSnExp,MatchedFilterParam works", {
    skip_on_os(os = "windows", arch = "i386")

    mfp <- MatchedFilterParam(binSize = 20, impute = "lin")
    onDisk <- filterFile(faahko_od, file = 1)
    res_o <- findChromPeaks(onDisk, param = mfp)
    expect_true(hasChromPeaks(res_o))
    expect_equal(nrow(chromPeaks(res_o)), 54)

    expect_error(findChromPeaks(onDisk, param = mfp, msLevel = 2))
})

test_that("isolationWindowTargetMz,OnDiskMSnExp works", {
    skip_on_os(os = "windows", arch = "i386")

    res <- isolationWindowTargetMz(xod_x)
    expect_true(all(is.na(res)))
    expect_true(length(res) == length(xod_x))

    f <- proteomics(full.names = TRUE)[5]
    tmt <- readMSData(f, mode = "onDisk")
    res <- isolationWindowTargetMz(tmt)
    expect_true(!all(is.na(res)))
})
