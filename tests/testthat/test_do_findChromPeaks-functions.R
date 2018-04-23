test_that("do_findPeaks_MSW works", {
    first_file <- filterFile(fticr, file = 1)
    spctr <- spectra(first_file)
    expect_true(length(spctr) == 1)
    mzs <- unname(mz(spctr[[1]]))
    ints <- unname(intensity(spctr[[1]]))
    feats1 <- do_findPeaks_MSW(mz = mzs[10000:20000],
                               int = ints[10000:20000],
                               snthresh = 100)
    feats2 <- do_findPeaks_MSW(mz = mzs[10000:20000],
                               int = ints[10000:20000],
                               snthresh = 50)
    expect_true(nrow(feats2) > nrow(feats1))
})

test_that("do_findChromPeaks_centWave works", {
    ## xr <- xcmsRaw(fs[1], profstep = 0)
    ## We expect that changing a parameter has an influence on the result.
    xr <- deepCopy(faahko_xr_1)
    mzVals <- xr@env$mz
    intVals <- xr@env$intensity
    ## Define the values per spectrum:
    valsPerSpect <- diff(c(xr@scanindex, length(mzVals)))
    res1 <- do_findChromPeaks_centWave(mz = mzVals,
                                       int = intVals,
                                       scantime = xr@scantime,
                                       valsPerSpect,
                                       snthresh = 200,
                                       noise = 4000)
    ## Eventually disable the sleep option to improve speed!
    res2 <- do_findChromPeaks_centWave(mz = mzVals,
                                       int = intVals,
                                       scantime = xr@scantime,
                                       valsPerSpect,
                                       snthresh = 500,
                                       noise = 4000, sleep = 0.00)
    expect_true(nrow(res1) > nrow(res2))

    ## Check scanrange on findPeaks.centWave.
    res_1 <- findPeaks.centWave(xr, scanrange = c(90, 345), noise = 2000)
    xr <- xr[90:345]
    mzVals <- xr@env$mz
    intVals <- xr@env$intensity
    ## Define the values per spectrum:
    valsPerSpect <- diff(c(xr@scanindex, length(mzVals)))
    res_2 <- do_findChromPeaks_centWave(mz = mzVals, int = intVals,
                                        scantime = xr@scantime, valsPerSpect,
                                        noise = 2000)
    expect_equal(res_1@.Data, res_2)
})

test_that("do_findChromPeaks_centWaveWithPredIsoROIs works", {
    xr <- deepCopy(faahko_xr_1)
    mzVals <- xr@env$mz
    intVals <- xr@env$intensity
    ## initial centWave:
    valsPerSpect <- diff(c(xr@scanindex, length(mzVals)))
    feats_1 <- do_findChromPeaks_centWave(mz = mzVals, int = intVals,
                                          scantime = xr@scantime,
                                          valsPerSpect = valsPerSpect,
                                          noise = 1500, verboseColumns = TRUE)
    feats_2 <- do_findChromPeaks_addPredIsoROIs(mz = mzVals,
                                                int = intVals,
                                                scantime = xr@scantime,
                                                valsPerSpect = valsPerSpect,
                                                noise = 1500,
                                                peaks. = feats_1)
    all_f <- do_findChromPeaks_centWaveWithPredIsoROIs(mz = mzVals,
                                                       int = intVals,
                                                       scantime = xr@scantime,
                                                       valsPerSpect = valsPerSpect,
                                                       noise = 1500)
    ## Comparisons.
    expect_equal(all_f, feats_2)
    ## old_all <- xcms:::.centWaveWithPredictedIsotopeROIs(xr, noise = 1500)
    ## checkEquals(all_f, old_all@.Data)
})

test_that("do_findChromPeaks_massifquant works", {
    xr <- deepCopy(faahko_xr_1)
    res <- findPeaks.massifquant(xr, snthresh = 100)
    mz <- xr@env$mz
    int <- xr@env$intensity
    valsPerSpect <- diff(c(xr@scanindex, length(mz)))
    scantime <- xr@scantime
    res_2 <- do_findChromPeaks_massifquant(mz = mz, int = int,
                                           valsPerSpect = valsPerSpect,
                                           scantime = scantime)
    expect_equal(res@.Data, res_2)
    ## With centWave:
    res_3 <- do_findChromPeaks_massifquant(mz = mz, int = int,
                                           valsPerSpect = valsPerSpect,
                                           scantime = scantime, withWave = TRUE,
                                           snthresh = 100, noise = 4000)
    res_4 <- findPeaks.massifquant(xr, withWave = 1, snthresh = 100,
                                   noise = 4000)
    expect_equal(res_3, res_4@.Data)
    expect_true(nrow(res_3) < nrow(res_2))

    ## Subsetted data and scanrange:
    res_1 <- findPeaks.massifquant(xr, scanrange = c(90, 345))
    xsub <- xr[90:345]
    mz <- xsub@env$mz
    int <- xsub@env$intensity
    valsPerSpect <- diff(c(xsub@scanindex, length(mz)))
    scantime <- xsub@scantime
    res_2 <- do_findChromPeaks_massifquant(mz = mz, int = int,
                                           valsPerSpect = valsPerSpect,
                                           scantime = scantime)
    expect_identical(res_1@.Data, res_2)
})

test_that("do_findChromPeaks_matchedFilter works", {
    xr <- deepCopy(faahko_xr_1)
    ## We expect that changing a parameter has an influence on the result.
    mzVals <- xr@env$mz
    intVals <- xr@env$intensity
    ## Define the values per spectrum:
    valsPerSpect <- diff(c(xr@scanindex, length(mzVals)))
    res1 <- do_findChromPeaks_matchedFilter(mz = mzVals,
                                            int = intVals,
                                            scantime = xr@scantime,
                                            valsPerSpect,
                                            binSize = 10)
    res2 <- do_findChromPeaks_matchedFilter(mz = mzVals,
                                            int = intVals,
                                            scantime = xr@scantime,
                                            valsPerSpect,
                                            binSize = 10,
                                            snthresh = 100)
    expect_true(nrow(res1) > nrow(res2))
    res2 <- do_findChromPeaks_matchedFilter(mz = mzVals,
                                            int = intVals,
                                            scantime = xr@scantime,
                                            valsPerSpect,
                                            binSize = 20)
    expect_true(nrow(res1) > nrow(res2))
})

test_that("peaksWithMatchedFilter is working", {
    od <- filterFile(faahko_od, file = 1)
    od_mf <- findChromPeaks(od, param = MatchedFilterParam())

    chr <- chromatogram(od, mz = c(272.1, 272.3))[1, 1]
    pks <- peaksWithMatchedFilter(intensity(chr), rtime(chr))
    pks_mf <- chromPeaks(od_mf, mz = c(272.1, 272.3))
    expect_equal(pks[, "rt"], pks_mf[, "rt"])
    expect_equal(pks[, "rtmin"], pks_mf[, "rtmin"])
    expect_equal(pks[, "rtmax"], pks_mf[, "rtmax"])
    expect_equal(pks[, "intf"], pks_mf[, "intf"])
    expect_equal(pks[, "into"], pks_mf[, "into"])
    expect_equal(pks[, "maxf"], pks_mf[, "maxf"])
    expect_equal(pks[, "maxo"], pks_mf[, "maxo"])
    ## Errors and empty data.
    expect_error(peaksWithMatchedFilter())
    expect_error(peaksWithMatchedFilter(int = rnorm(10)))
    expect_error(peaksWithMatchedFilter(int = rnorm(10)), rt = 1:4)
    expect_true(nrow(peaksWithMatchedFilter(rep(NA, 10), rt = 1:10)) == 0)
})
