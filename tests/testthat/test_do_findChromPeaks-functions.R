test_that("do_findPeaks_MSW works", {
    skip_on_os(os = "windows", arch = "i386")

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
    skip_on_os(os = "windows", arch = "i386")

    ## We expect that changing a parameter has an influence on the result.
    mzVals <- mz(xod_x)
    intVals <- unlist(intensity(xod_x), use.names = FALSE)
    ## Define the values per spectrum:
    valsPerSpect <- lengths(mzVals)
    res1 <- do_findChromPeaks_centWave(mz = unlist(mzVals, use.names = FALSE),
                                       int = intVals,
                                       scantime = rtime(xod_x),
                                       valsPerSpect,
                                       snthresh = 100,
                                       noise = 4000,
                                       prefilter = c(3, 10000))
    ## Eventually disable the sleep option to improve speed!
    res2 <- do_findChromPeaks_centWave(mz = unlist(mzVals, use.names = FALSE),
                                       int = intVals,
                                       scantime = rtime(xod_x),
                                       valsPerSpect,
                                       snthresh = 500,
                                       noise = 4000, sleep = 0.00,
                                       prefilter = c(3, 10000))
    expect_true(nrow(res1) > nrow(res2))
})

test_that("beta calculation returns expected values", {
  expect_lt(.get_beta_values(1:10, zero.rm = FALSE)["best_cor"], 0.0001)
  expect_lt(.get_beta_values(1:10, zero.rm = FALSE)["beta_snr"], 2)
  expect_lt(.get_beta_values(1:10)["best_cor"], 0.0001)
  expect_lt(.get_beta_values(1:10)["beta_snr"], 2)

  ideal_beta <- dbeta(seq(0, 1, length.out=10), 5, 5)
  expect_gte(.get_beta_values(ideal_beta, zero.rm = FALSE)["best_cor"], 1)
  expect_gte(.get_beta_values(ideal_beta, zero.rm = FALSE)["beta_snr"], 16)
  expect_gte(.get_beta_values(ideal_beta)["best_cor"], 0.97)
  expect_gte(.get_beta_values(ideal_beta)["beta_snr"], 1)

  skew_beta <- dbeta(seq(0, 1, length.out=10), 3, 5)
  expect_gte(.get_beta_values(ideal_beta, zero.rm = FALSE)["best_cor"], 1)
  expect_gte(.get_beta_values(ideal_beta, zero.rm = FALSE)["beta_snr"], 16)
  expect_gte(.get_beta_values(ideal_beta)["best_cor"], 0.97)
  expect_gte(.get_beta_values(ideal_beta)["beta_snr"], 1)

  rightskew_beta <- dbeta(seq(0, 1, length.out=10), 7, 5)
  expect_gt(.get_beta_values(rightskew_beta, skews = c(3,5,7))["best_cor"], 0.95)

  noise_beta <- dbeta(seq(0, 1, length.out=21), 5, 5)*10+runif(21)
  expect_gt(.get_beta_values(noise_beta)["best_cor"], 0.9)

  expect_no_error(.get_beta_values(runif(1)))
  expect_no_error(.get_beta_values(runif(10)))
  expect_no_error(.get_beta_values(runif(100)))

  expect_length(.get_beta_values(1), 2)
  expect_true(is.na(.get_beta_values(1)["best_cor"]))
  expect_true(is.na(.get_beta_values(1)["beta_snr"]))
})

test_that("New beta columns perform as expected", {
  # faahko_xod comes from testthat.R
  # faahko_xod <- findChromPeaks(
  #   faahko_od, param = CentWaveParam(noise = 10000, snthresh = 40,
  #                                    prefilter = c(3, 10000)))
  # Same params as before but with verboseBetaColumns = TRUE
  beta_cwp <- CentWaveParam(noise = 10000, snthresh = 40,
                            prefilter = c(3, 10000),
                            verboseBetaColumns = TRUE)
  faahko_xod_beta <- findChromPeaks(faahko_od, beta_cwp)

  # Unit test - check that the new object contains expected columns
  expect_contains(colnames(chromPeaks(faahko_xod_beta)),
                  c("beta_cor", "beta_snr"))

  # Unit test - check that everything else in the object is the same
  orig_chrompeaks <- chromPeaks(faahko_xod)
  beta_chrompeaks <- chromPeaks(faahko_xod_beta)
  expect_identical(orig_chrompeaks, beta_chrompeaks[,colnames(orig_chrompeaks)])

  # Object will contain NAs because there are peaks <5 scans wide
  expect_true(any(is.na(beta_chrompeaks[,"beta_snr"])))
  beta_chrompeaks <- beta_chrompeaks[!is.na(beta_chrompeaks[,"beta_cor"]),]

  # Unit test - check that beta values make sense
  expect_true(all(beta_chrompeaks[,"beta_cor"]<=1))
  expect_true(all(beta_chrompeaks[,"beta_cor"]>=-1))
  expect_true(all(beta_chrompeaks[,"beta_snr"]>=0))

  # Unit test - finds a single good peak (beta_cor>0.8, beta_snr>7)
  # Skinny peak copied from below peaksWithCentWave tests
  skinny_peak <- c(9107, 3326, 9523, 3245, 3429, 9394, 1123, 935, 5128, 8576,
                   2711, 3427, 7294, 8109, 9288, 6997, 9756, 8034, 1317, 8866, 13877,
                   14854, 28296, 57101, 92209, 151797, 222386, 299402, 365045, 394255,
                   402680, 363996, 293985, 222989, 147007, 94947, 52924, 32438,
                   11511, 10836, 8046, 601, 889, 5917, 2690, 5381, 9901, 8494, 3349,
                   8283, 3410, 5935, 3332, 7041, 3284, 7478, 76, 3739, 2158, 5507)
  skinny_peak_rt <- seq_along(skinny_peak)+100
  cw_output_beta <- .centWave_orig(int = skinny_peak, scantime = skinny_peak_rt,
                                   mz=sort(rnorm(60)/1000+100),
                                   valsPerSpect = rep(1, length(skinny_peak)),
                                   peakwidth = c(20, 80), extendLengthMSW = TRUE,
                                   verboseBetaColumns = TRUE, snthresh = 0)
  expect_equal(nrow(cw_output_beta), 1)
  # Known values to ensure performance doesn't degrade unexpectedly
  expect_gt(cw_output_beta[,"beta_cor"], 0.8)
  expect_gt(cw_output_beta[,"beta_snr"], 7)

  # Unit test - finds a single noise peak (beta_cor < 0.5, beta_snr < 6)
  # set.seed(123)
  # noise_peak <- round(runif(100), 3)
  noise_peak <- c(0.288, 0.788, 0.409, 0.883, 0.94, 0.046, 0.528, 0.892, 0.551,
                  0.457, 0.957, 0.453, 0.678, 0.573, 0.103, 0.9, 0.246, 0.042,
                  0.328, 0.955, 0.89, 0.693, 0.641, 0.994, 0.656, 0.709, 0.544,
                  0.594, 0.289, 0.147, 0.963, 0.902, 0.691, 0.795, 0.025, 0.478,
                  0.758, 0.216, 0.318, 0.232, 0.143, 0.415, 0.414, 0.369, 0.152,
                  0.139, 0.233, 0.466, 0.266, 0.858, 0.046, 0.442, 0.799, 0.122,
                  0.561, 0.207, 0.128, 0.753, 0.895, 0.374, 0.665, 0.095, 0.384,
                  0.274, 0.815, 0.449, 0.81, 0.812, 0.794, 0.44, 0.754, 0.629,
                  0.71, 0.001, 0.475, 0.22, 0.38, 0.613, 0.352, 0.111, 0.244,
                  0.668,
                  0.418, 0.788, 0.103, 0.435, 0.985, 0.893, 0.886, 0.175, 0.131,
                  0.653, 0.344, 0.657, 0.32, 0.188, 0.782, 0.094, 0.467, 0.512)
  cw_output_beta <- .centWave_orig(int = noise_peak*100000,
                                   scantime = seq_along(noise_peak),
                                   mz=rep(530.1, length(noise_peak)),
                                   valsPerSpect = rep(1, length(noise_peak)),
                                   peakwidth = c(20, 80), extendLengthMSW = TRUE,
                                   verboseBetaColumns = TRUE, snthresh = 0)
  expect_equal(nrow(cw_output_beta), 1)
  # Known values to ensure performance doesn't degrade unexpectedly
  expect_lt(cw_output_beta[,"beta_cor"], 0.2)
  expect_lt(cw_output_beta[,"beta_snr"], 6)
})

test_that(".get_beta_values works with chromatographic and MS data", {
    vals <- c(2, 3, 4, 6, 8, 10, 7, 6, 4, 3, 2)
    res <- .get_beta_values(vals)
    ## intensity values with duplicated retention times
    vals <- c(2, 3, 4, 6, 2, 8, 10, 7, 6, 4, 3, 2)
    rts <- c(1, 2, 3, 4, 4, 5, 6, 7, 8, 9, 10)
    res_a <- .get_beta_values(vals)
    expect_true(res_a[1L] < res[1L])
    ## WARNING: does not work if we have duplicated retentio times!
    expect_warning(
        res_b <- .get_beta_values(vals, rts)
    )
    expect_true(is.na(res_b[1L]))
})

test_that("do_findChromPeaks_centWaveWithPredIsoROIs works", {
    skip_on_os(os = "windows", arch = "i386")

    tmp <- filterFile(od_x, 1)
    mzVals <- mz(tmp)
    intVals <- unlist(intensity(tmp), use.names = FALSE)
    ## initial centWave:
    valsPerSpect <- lengths(mzVals)
    feats_1 <- do_findChromPeaks_centWave(
        mz = unlist(mzVals, use.names = FALSE), int = intVals,
        prefilter = c(3, 5000),
        scantime = rtime(tmp),
        valsPerSpect = valsPerSpect, noise = 1500, verboseColumns = TRUE)
    feats_2 <- do_findChromPeaks_addPredIsoROIs(
        mz = unlist(mzVals, use.names = FALSE), int = intVals,
        scantime = rtime(tmp), valsPerSpect = valsPerSpect, noise = 1500,
        prefilter = c(3, 5000),
        peaks. = feats_1)
    expect_true(nrow(feats_1) < nrow(feats_2))
    all_f <- do_findChromPeaks_centWaveWithPredIsoROIs(
        mz = unlist(mzVals, use.names = FALSE), int = intVals,
        scantime = rtime(tmp), prefilter = c(3, 5000),
        valsPerSpect = valsPerSpect, noise = 1500)
    expect_equal(all_f, feats_2)
})

test_that("do_findChromPeaks_massifquant works", {
    skip_on_os(os = "windows", arch = "i386")

    tmp <- filterFile(od_x, 1)
    mz <- mz(tmp)
    valsPerSpect <- lengths(mz)
    mz <- unlist(mz, use.names = FALSE)
    int <- unlist(intensity(tmp), use.names = FALSE)
    rtime <- rtime(tmp)
    res_2 <- do_findChromPeaks_massifquant(mz = mz, int = int,
                                           valsPerSpect = valsPerSpect,
                                           scantime = rtime)
    expect_true(nrow(res_2) == 1542)
    ## With centWave:
    res_3 <- do_findChromPeaks_massifquant(mz = mz, int = int,
                                           valsPerSpect = valsPerSpect,
                                           scantime = rtime, withWave = TRUE,
                                           snthresh = 100, noise = 4000)
    expect_true(nrow(res_3) < nrow(res_2))
})

test_that("do_findChromPeaks_matchedFilter works", {
    skip_on_os(os = "windows", arch = "i386")

    tmp <- filterFile(od_x, 1)
    mzVals <- mz(tmp)
    valsPerSpect <- lengths(mzVals)
    mzVals <- unlist(mzVals, use.names = FALSE)
    intVals <- unlist(intensity(tmp), use.names = FALSE)
    res1 <- do_findChromPeaks_matchedFilter(mz = mzVals,
                                            int = intVals,
                                            scantime = rtime(tmp),
                                            valsPerSpect,
                                            binSize = 10)
    res2 <- do_findChromPeaks_matchedFilter(mz = mzVals,
                                            int = intVals,
                                            scantime = rtime(tmp),
                                            valsPerSpect,
                                            binSize = 10,
                                            snthresh = 100)
    expect_true(nrow(res1) > nrow(res2))
    res2 <- do_findChromPeaks_matchedFilter(mz = mzVals,
                                            int = intVals,
                                            scantime = rtime(tmp),
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
    skip_on_os(os = "windows", arch = "i386")

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
    skip_on_os(os = "windows", arch = "i386")

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
    skip_on_os(os = "windows", arch = "i386")

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

    #Testing for Github issue #445 after introducing new CWP method
    # Values from round((c(numeric(20), dnorm(seq(-3, 3, length.out = 20)),
    #             numeric(20))*100 +runif(60))*10000)
    skinny_peak <- c(9107, 3326, 9523, 3245, 3429, 9394, 1123, 935, 5128, 8576,
                     2711, 3427, 7294, 8109, 9288, 6997, 9756, 8034, 1317, 8866, 13877,
                     14854, 28296, 57101, 92209, 151797, 222386, 299402, 365045, 394255,
                     402680, 363996, 293985, 222989, 147007, 94947, 52924, 32438,
                     11511, 10836, 8046, 601, 889, 5917, 2690, 5381, 9901, 8494, 3349,
                     8283, 3410, 5935, 3332, 7041, 3284, 7478, 76, 3739, 2158, 5507)
    skinny_peak_rt <- seq_along(skinny_peak)+100
    pks <- peaksWithCentWave(skinny_peak, rt=skinny_peak_rt,
                             snthresh = 0, peakwidth = c(20, 50),
                             extendLengthMSW = TRUE)
    expect_true(nrow(pks)==1)

    # Reducing minimum peakwidth shouldn't affect peak detection
    pks_widerpeakwidth <- peaksWithCentWave(skinny_peak, rt=skinny_peak_rt,
                                            snthresh = 0, peakwidth = c(2, 50),
                                            extendLengthMSW = TRUE)
    expect_true(nrow(pks)==nrow(pks_widerpeakwidth))

    # Test a wider peak
    # Values from round((dnorm(seq(-3, 3, length.out = 60))*100+runif(60))*10000)
    wider_peak <- c(5000, 12043, 15344, 12748, 20730, 20781, 24673, 36956, 44600,
                    48596, 57698, 76937, 89422, 106482, 122977, 143989, 157769, 181563,
                    206296, 226309, 251067, 283592, 307523, 324212, 341520, 368568,
                    375716, 388428, 401694, 408352, 399415, 403964, 394144, 382952,
                    368333, 341668, 330255, 301146, 276234, 254643, 231601, 211038,
                    184239, 155817, 140996, 123284, 100121, 90280, 77303, 58708,
                    52817, 44003, 36068, 24637, 20688, 14162, 14836, 16603, 8341,
                    8307)
    wider_peak_rt <- seq_along(wider_peak)+100
    pks <- peaksWithCentWave(wider_peak, rt=wider_peak_rt,
                             snthresh = 0, peakwidth = c(20, 80),
                             extendLengthMSW = TRUE)
    expect_true(nrow(pks)==1)
    pks_widerpeakwidth <- peaksWithCentWave(skinny_peak, rt=skinny_peak_rt,
                                            snthresh = 0, peakwidth = c(2, 80),
                                            extendLengthMSW = TRUE)
    expect_true(nrow(pks)==nrow(pks_widerpeakwidth))


    ## Check errors
    expect_error(peaksWithCentWave())
    expect_error(peaksWithCentWave(int = 1:3, rt = 1:5))
    expect_warning(res <- peaksWithCentWave(int = rep(NA, 20), rt = 1:20))
    expect_true(nrow(res) == 0)
})

test_that(".narrow_rt_boundaries works", {
    skip_on_os(os = "windows", arch = "i386")

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
