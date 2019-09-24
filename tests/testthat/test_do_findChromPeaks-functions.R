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
    mzVals <- faahko_xr_1@env$mz
    intVals <- faahko_xr_1@env$intensity
    ## initial centWave:
    valsPerSpect <- diff(c(faahko_xr_1@scanindex, length(mzVals)))
    feats_1 <- do_findChromPeaks_centWave(
        mz = mzVals, int = intVals, scantime = faahko_xr_1@scantime,
        valsPerSpect = valsPerSpect, noise = 1500, verboseColumns = TRUE)
    feats_2 <- do_findChromPeaks_addPredIsoROIs(
        mz = mzVals, int = intVals, scantime = faahko_xr_1@scantime,
        valsPerSpect = valsPerSpect, noise = 1500, peaks. = feats_1)
    expect_true(nrow(feats_1) < nrow(feats_2))
    all_f <- do_findChromPeaks_centWaveWithPredIsoROIs(
        mz = mzVals, int = intVals, scantime = faahko_xr_1@scantime,
        valsPerSpect = valsPerSpect, noise = 1500)
    expect_equal(all_f, feats_2)
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

    ## with empty spectra - simulating issue #325
    od_sub <- filterMz(od_x, mz = c(334.9, 344.1))
    sps <- spectra(filterFile(od_sub, 1))
    ## Add an artificial peak at m/z 0 if spectrum is empty
    sps <- lapply(sps, function(z) {
        if (!length(z@mz)) {
            z@mz <- 0.0
            z@intensity <- 0.0
        }
        z
    })
    mzs <- lapply(sps, mz)
    n_peaks <- lengths(mzs, FALSE)
    mzs <- unlist(mzs, use.names = FALSE)
    ints <- unlist(lapply(sps, intensity), use.names = FALSE)
    rtms <- vapply(sps, rtime, numeric(1))
    res3 <- do_findChromPeaks_matchedFilter(mz = mzs, int = ints,
                                            scantime = rtms,
                                            valsPerSpect = n_peaks)
    full_data <- findChromPeaks(filterFile(od_x, 1),
                                param = MatchedFilterParam())
    pks_full <- chromPeaks(full_data, mz = c(335, 344))
    rownames(pks_full) <- NULL
    rownames(res3) <- NULL
    expect_equal(res3, pks_full[, colnames(res3)])
    res4 <- findChromPeaks(filterMz(filterFile(od_x, 1), mz = c(334.9, 344.1)),
                           param = MatchedFilterParam())
    res4 <- chromPeaks(res4)
    rownames(res4) <- NULL
    expect_equal(res4, pks_full)
})

test_that("peaksWithMatchedFilter is working", {
    od <- filterFile(faahko_od, file = 1)
    od_mf <- findChromPeaks(od, param = MatchedFilterParam())

    chr <- chromatogram(od, mz = c(272.1, 272.3))[1, 1]
    pks <- peaksWithMatchedFilter(intensity(chr), rtime(chr))
    pks_mf <- chromPeaks(od_mf, mz = c(272.1, 272.3))
    rownames(pks_mf) <- NULL
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

test_that(".getRtROI works", {
    od <- filterFile(faahko_od, file = 1)
    expect_error(.getRtROI())
    expect_error(.getRtROI(1:3))
    expect_error(.getRtROI(1:3, 1:5))
    chr <- chromatogram(od, mz = c(272.1, 272.3))[1, 1]

    int <- intensity(chr)
    int[is.na(int)] <- 0
    rt <- rtime(chr)
    res <- .getRtROI(int, rt)
    expect_true(is.matrix(res))
    expect_true(ncol(res) == 3)
    res_2 <- .getRtROI(int, rt, noise = 400)
    expect_true(nrow(res) > nrow(res_2))
    res_3 <- .getRtROI(int, rt, noise = 400, prefilter = c(4, 500))
    expect_true(nrow(res_2) > nrow(res_3))
    res_4 <- .getRtROI(int, rt, noise = 400, prefilter = c(100, 500))
    expect_true(nrow(res_4) == 0)
    
    # Generate a nice-looking peak 
    # Values from table(cut(rnorm(20000), breaks = 40))
    model_peak <- c(3, 4, 4, 9, 26, 31, 65, 123, 196, 260, 404, 523, 743, 893, 
                    1188, 1329, 1505, 1540, 1705, 1592, 1535, 1371, 1255, 929, 790, 
                    652, 438, 336, 223, 138, 78, 50, 25, 15, 11, 6, 3, 0, 1, 1)
    model_single_peak <- c(numeric(80), model_peak, numeric(80))
    single_peak_scans <- seq_along(model_single_peak)+200
    single_peak_rois <- .getRtROI(model_single_peak, single_peak_scans)
    expect_true(is.matrix(single_peak_rois))
    expect_true(nrow(single_peak_rois)==1)
    
    model_triple_peak <- c(numeric(20), model_peak, numeric(100),
                           model_peak/5, numeric(100),
                           rev(model_peak)*2, numeric(20))
    triple_peak_scans <- seq_along(model_triple_peak)+200
    # Get ROIs for a chromatogram with 3 good peaks
    triple_peak_rois <- .getRtROI(model_triple_peak, triple_peak_scans)
    expect_true(nrow(triple_peak_rois)==3)
    
    # Get ROIs for a chromatogram with three peaks
    # One of which doesn't pass prefilter check
    skipped_peak_rois <- .getRtROI(model_triple_peak, triple_peak_scans, 
                                   prefilter = c(3, 500))
    expect_true(nrow(skipped_peak_rois)==2)
    
    # Get ROIs for a chromatogram with three peaks
    # None of which pass stringent prefilter check
    skipped_peak_rois <- .getRtROI(model_triple_peak, triple_peak_scans, 
                                   prefilter = c(3, 5000))
    expect_true(nrow(skipped_peak_rois)==0)
    
    # Get ROIs for a chromatogram with three peaks
    # One of which passes stringent prefilter check
    tall_peak_rois <- .getRtROI(model_triple_peak, triple_peak_scans, 
                                prefilter = c(9, 1500))
    expect_true(nrow(tall_peak_rois)==1)
})

test_that("peaksWithCentWave works", {
    od <- filterFile(faahko_od, file = 1)
    mzr <- c(272.1, 272.2)

    od_cw <- findChromPeaks(filterMz(od, mz = c(270, 300)),
                            param = CentWaveParam())

    chr <- chromatogram(od, mz = mzr)[1, 1]
    pks <- peaksWithCentWave(intensity(chr), rtime(chr))
    pks_cw <- chromPeaks(od_cw, mz = mzr)
    rownames(pks_cw) <- NULL
    expect_equal(pks[2, "rt"], pks_cw[, "rt"])
    expect_equal(pks[2, "rtmin"], pks_cw[, "rtmin"])
    expect_equal(pks[2, "rtmax"], pks_cw[, "rtmax"])
    expect_equal(pks[2, "into"], pks_cw[, "into"])

    cwp <- CentWaveParam(fitgauss = TRUE)
    pks <- peaksWithCentWave(intensity(chr), rtime(chr), fitgauss = TRUE)

    ## Check errors
    expect_error(peaksWithCentWave())
    expect_error(peaksWithCentWave(int = 1:3, rt = 1:5))
    expect_warning(res <- peaksWithCentWave(int = rep(NA, 20), rt = 1:20))
    expect_true(nrow(res) == 0)
})

test_that(".narrow_rt_boundaries works", {
    d <- c(0, 0, 1, 2, 1, 3, 4, 6, 4, 3, 2, 0, 1, 0, 2, 0)

    ## Full range
    lm <- c(1, length(d))
    res <- .narrow_rt_boundaries(lm, d)
    expect_equal(res, c(2, 16))
    res <- .narrow_rt_boundaries(lm, d, thresh = 2)
    expect_equal(res, c(3, 16))
    res <- .narrow_rt_boundaries(lm, d, thresh = 3)
    expect_equal(res, c(5, 11))

    ## Subset (reflecting the real situation).
    lm <- c(3, 9)
    res <- .narrow_rt_boundaries(lm, d)
    expect_equal(res, c(3, 9))
    res <- .narrow_rt_boundaries(lm, d, thresh = 2)
    expect_equal(res, c(3, 9))
    res <- .narrow_rt_boundaries(lm, d, thresh = 3)
    expect_equal(res, c(5, 9))

    lm <- c(3, 13)
    res <- .narrow_rt_boundaries(lm, d)
    expect_equal(res, c(3, 13))
    res <- .narrow_rt_boundaries(lm, d, thresh = 3)
    expect_equal(res, c(5, 11))

    ## That's the fix for issue #300
    expect_equal(.narrow_rt_boundaries(lm, d, thresh = 100), lm)
    expect_equal(.narrow_rt_boundaries(c(1, length(d)), d, thresh = 100),
                 c(1, length(d)))
})
