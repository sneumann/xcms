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
                                       noise = 4000,
                                       prefilter = c(3, 10000))
    ## Eventually disable the sleep option to improve speed!
    res2 <- do_findChromPeaks_centWave(mz = mzVals,
                                       int = intVals,
                                       scantime = xr@scantime,
                                       valsPerSpect,
                                       snthresh = 500,
                                       noise = 4000, sleep = 0.00,
                                       prefilter = c(3, 10000))
    expect_true(nrow(res1) > nrow(res2))

    ## Check scanrange on findPeaks.centWave.
    res_1 <- findPeaks.centWave(xr, scanrange = c(90, 345), noise = 2000,
                                prefilter = c(3, 10000))
    xr <- xr[90:345]
    mzVals <- xr@env$mz
    intVals <- xr@env$intensity
    ## Define the values per spectrum:
    valsPerSpect <- diff(c(xr@scanindex, length(mzVals)))
    res_2 <- do_findChromPeaks_centWave(mz = mzVals, int = intVals,
                                        scantime = xr@scantime, valsPerSpect,
                                        noise = 2000, prefilter = c(3, 10000))
    expect_equal(res_1@.Data, res_2)
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
