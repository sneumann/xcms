test_that("profMat,xcmsRaw works", {
    fs <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    ## Create a new profile matrix:
    xr <- deepCopy(faahko_xr_1)
    expect_error(pm <- profMat(xr))
    xr_2 <- xcmsRaw(fs, profstep = 2)
    expect_equal(profMat(xr_2), profMat(xr, step = 2))
    xr_3 <- xcmsRaw(fs, profstep = 2, profmethod = "binlinbase",
                    profparam = list(baselevel = 666666))
    expect_equal(xr_3@env$profile, profMat(xr, step = 2, method = "binlinbase",
                                           baselevel = 666666))
    expect_equal(xr_3@env$profile, profMat(xr_3))
})

test_that("profStep<-,xcmsRaw works", {
    ## Profile matrix will be generated/replaced if the step parameter is > 0
    ## and differs from the one within the object.
    ## xr <- xcmsRaw(fs, profstep = 0)
    xr <- deepCopy(faahko_xr_1)
    xr_2 <- xr
    expect_true(length(xr_2@env$profile) == 0)
    profStep(xr_2) <- 2
    expect_true(length(xr_2@env$profile) > 0)
    expect_equal(profMat(xr_2), profMat(xr, step = 2))
    expect_equal(profMat(xr_2, step = 2), xr_2@env$profile)
})

test_that("profMethod<-,xcmsRaw works", {
    fs <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    ## Profile matrix will be generated/replaced if profMethod is changed.
    xr <- deepCopy(faahko_xr_1)
    xr_2 <- xr
    expect_true(length(xr_2@env$profile) == 0)
    ## Just setting profMethod doesn't help here
    profMethod(xr_2) <- "binlin"
    expect_equal(profMethod(xr_2), "binlin")
    expect_true(length(xr_2@env$profile) == 0)
    profStep(xr_2) <- 2
    xr_3 <- xcmsRaw(fs, profstep = 2, profmethod = "binlin")
    expect_equal(profMat(xr_3), profMat(xr_2))
    ## binlinbase
    xr_3@profparam <- list(baselevel = 666666)
    profMethod(xr_3) <- "binlinbase"
    expect_equal(xr_3@env$profile, profMat(xr_3, method = "binlinbase",
                                           baselevel = 666666))
    xr_4 <- xcmsRaw(fs, profstep = 2, profmethod = "binlinbase",
                    profparam = list(baselevel = 666666))
    expect_equal(xr_4@env$profile, xr_3@env$profile)
})

test_that("findPeaks.centWave,xcmsRaw ordering works", {
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xr <- xcmsRaw(file, profstep = 0)
    p1 <- findPeaks.centWave(xr, noise = 10000)
    scan <- getScan(xr,1)
    o <- order(scan[,"mz"], decreasing=TRUE)

    xr@env$mz[1:length(o)] <- scan[o, "mz"]
    xr@env$intensity[1:length(o)] <- scan[o, "intensity"]

    p2 <- findPeaks.centWave(xr, noise = 10000)
    expect_equal(p1, p2)
})

test_that("findPeaksCentWave,xcmsSet doesn't fail", {
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xraw <- xcmsRaw(file, profstep = 0)
    p <- findPeaks.centWave(xraw, fitgauss = TRUE, verbose = TRUE, sleep = 0.00,
                            noise = 10000)
    expect_equal(nrow(p), 272)
})

test_that("findPeaks.addPredictedIsotopeFeatures,xcmsRaw works", {
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xr <- xcmsRaw(file, profstep = 0)
    p1 <- findPeaks.centWave(xr, verbose.columns = TRUE, noise = 10000)
    p2 <- findPeaks.addPredictedIsotopeFeatures(
        object = xr, xcmsPeaks = p1, noise = 10000)
    expect_true(nrow(p1) < nrow(p2))
    ## Now the same with the new modified centWave:
    options(originalCentWave = FALSE)
    p1_2 <- findPeaks.centWave(xr, verbose.columns = TRUE, noise = 10000)
    nrow(p1)
    nrow(p1_2)
    expect_equal(p1, p1_2)
    p2_2 <- findPeaks.addPredictedIsotopeFeatures(
        object = xr, xcmsPeaks = p1_2, noise = 10000)
    options(originalCentWave = TRUE)    
})

test_that("findPeaks,xcmsRaw massifquant works", {
    file <- system.file('microtofq/MM14.mzdata', package = "msdata")
    xraw <- xcmsRaw(file, profstep = 0)
    p <- findPeaks(xraw, method = "massifquant")
    expect_equal(nrow(p), 114)
})

test_that("profEIC,xcmsRaw works", {
    ## file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    ## Compare the results with manual calculations on the profile matrix.
    step <- 1
    ## xraw <- xcmsRaw(file, profstep = 0)
    xraw <- deepCopy(faahko_xr_1)
    profmat <- profMat(xraw, step = step)
    mzr <- c(200, 201)
    rtr <- c(3000, 3500)
    res_1 <- profEIC(xraw, mzrange = mzr, rtrange = rtr, step = step)
    res_x <- getEIC(xraw, mzrange = mzr, rtrange = rtr, step = step)
    ## manual calculation:
    scn_idx <- which(xraw@scantime >= rtr[1] & xraw@scantime <= rtr[2])
    mass <- seq(floor(xraw@mzrange[1] / step) * step,
                ceiling(xraw@mzrange[2] / step) * step, by = step)
    mass_idx <- which(mass >= mzr[1] & mass <= mzr[2])
    max_mass <- apply(profmat[mass_idx, scn_idx], MARGIN = 2, max)
    rts <- xraw@scantime[scn_idx]
    expect_equal(res_1@eic[[1]][[1]][, 1], rts)
    expect_equal(res_1@eic[[1]][[1]][, 2], max_mass)
    expect_equal(res_1, res_x)

    ## OK, now with binlinbase
    profMethod(xraw) <- "binlinbase"
    xraw@profparam <- list(baselevel = 666666)
    profmat <- profMat(xraw, step = 1)
    res_2 <- profEIC(xraw, mzrange = mzr, rtrange = rtr, step = step)
    res_x <- getEIC(xraw, mzrange = mzr, rtrange = rtr, step = step)
    max_mass <- apply(profmat[mass_idx, scn_idx], MARGIN = 2, max)
    expect_equal(res_2@eic[[1]][[1]][, 1], rts)
    expect_equal(res_2@eic[[1]][[1]][, 2], max_mass)
    expect_true(sum(res_2@eic[[1]][[1]][, 2] == 666666) > 0)
    expect_equal(res_2, res_x)
    ## Results should be different to previous ones.
    expect_true(any(res_2@eic[[1]][[1]][, 2] != res_1@eic[[1]][[1]][, 2]))
    expect_true(all(res_2@eic[[1]][[1]][, 1] == res_1@eic[[1]][[1]][, 1]))
    ## Without pre-calculated profile matrix:
    profmat <- profMat(xraw, step = 0.2)
    res_3 <- profEIC(xraw, mzrange = mzr, rtrange = rtr, step = 0.2)
    res_x <- getEIC(xraw, mzrange = mzr, rtrange = rtr, step = 0.2)
    max_mass <- apply(profmat[mass_idx, scn_idx], MARGIN = 2, max)
    expect_equal(res_3@eic[[1]][[1]][, 1], rts)
    expect_equal(res_3@eic[[1]][[1]][, 2], max_mass)
    expect_equal(res_3, res_x)
    ## Results should be different to previous ones.
    expect_true(any(res_3@eic[[1]][[1]][, 2] != res_2@eic[[1]][[1]][, 2]))
    expect_true(all(res_3@eic[[1]][[1]][, 1] == res_2@eic[[1]][[1]][, 1]))

    ## Now with intlin
    profMethod(xraw) <- "intlin"
    profmat <- profMat(xraw, step = 1)
    res_3 <- profEIC(xraw, mzrange = mzr, rtrange = rtr, step = step)
    res_x <- getEIC(xraw, mzrange = mzr, rtrange = rtr, step = step)
    max_mass <- apply(profmat[mass_idx, scn_idx], MARGIN = 2, max)
    expect_equal(res_3@eic[[1]][[1]][, 1], rts)
    expect_equal(res_3@eic[[1]][[1]][, 2], max_mass)
    expect_equal(res_3, res_x)

    ## Now test with multiple mz ranges
    mzrm <- rbind(c(300, 302), mzr, c(205, 207), c(500, 507))
    profMethod(xraw) <- "bin"
    profmat <- profMat(xraw, step = 1)
    max_mass <- apply(profmat[mass_idx, scn_idx], MARGIN = 2, max)
    res_1 <- profEIC(xraw, mzrange = mzrm, rtrange = rtr, step = step)
    res_x <- getEIC(xraw, mzrange = mzrm, rtrange = rtr, step = step)
    expect_equal(res_1@eic[[1]][[2]][, 2], max_mass)
    ## 4th range:
    mass_idx_2 <- which(mass >= 500 & mass <= 507)
    max_mass <- apply(profmat[mass_idx_2, scn_idx], MARGIN = 2, max)
    expect_equal(res_1@eic[[1]][[4]][, 2], max_mass)
    expect_equal(res_1, res_x)

    ## Different step:
    step <- 0.2
    res_1 <- profEIC(xraw, mzrange = mzrm, rtrange = rtr, step = step)
    res_x <- getEIC(xraw, mzrange = mzrm, rtrange = rtr, step = step)
    mass <- seq(floor(xraw@mzrange[1] / step) * step,
                ceiling(xraw@mzrange[2] / step) * step, by = step)
    mass_idx <- which(mass >= mzr[1] & mass <= mzr[2])
    scn_idx <- which(xraw@scantime >= rtr[1] & xraw@scantime <= rtr[2])
    profmat <- profMat(xraw, step = step)
    max_mass <- apply(profmat[mass_idx, scn_idx], MARGIN = 2, max)
    expect_equal(res_1@eic[[1]][[2]][, 2], max_mass)
    expect_equal(res_1, res_x)
    ## 4th range:
    mass_idx_2 <- which(mass >= 500 & mass <= 507)
    max_mass <- apply(profmat[mass_idx_2, scn_idx], MARGIN = 2, max)
    expect_equal(res_1@eic[[1]][[4]][, 2], max_mass)

    ## Test with the whole mzrange.
    step <- 0.5
    profMethod(xraw) <- "bin"
    profmat <- profMat(xraw, step = step)
    res_1 <- profEIC(xraw, rtrange = rtr, step = step)
    res_x <- getEIC(xraw, rtrange = rtr, step = step)
    scn_idx <- which(xraw@scantime >= rtr[1] & xraw@scantime <= rtr[2])
    max_mass <- apply(profmat[, scn_idx], MARGIN = 2, max)
    expect_equal(res_1@eic[[1]][[1]][, 2], max_mass) ## Compare max per rt
    expect_equal(res_1@eic[[1]][[1]][, 1], xraw@scantime[scn_idx])
    expect_equal(res_1, res_x)

    ## Whole rtrange:
    res_1 <- profEIC(xraw, mzrange = mzr, step = step)
    res_x <- getEIC(xraw, mzrange = mzr, step = step)
    mass <- seq(floor(xraw@mzrange[1] / step) * step,
                ceiling(xraw@mzrange[2] / step) * step, by = step)
    mass_idx <- which(mass >= mzr[1] & mass <= mzr[2])
    max_mass <- apply(profmat[mass_idx, ], MARGIN = 2, max)
    expect_equal(res_1@eic[[1]][[1]][, 2], max_mass)
    expect_equal(res_1@eic[[1]][[1]][, 1], xraw@scantime)
    expect_equal(res_1, res_x)

    ## Everything.
    res_1 <- profEIC(xraw, step = step)
    res_x <- getEIC(xraw, step = step)
    max_mass <- apply(profmat, MARGIN = 2, max)
    expect_equal(res_1@eic[[1]][[1]][, 2], max_mass)
    expect_equal(res_1, res_x)

    ## Error checking:
    ## mz range outside
    mzr <- c(0, 14)
    expect_error(profEIC(xraw, mzrange = mzr))
    expect_error(getEIC(xraw, mzrange = mzr))
    ## rt range outside
    expect_error(profEIC(xraw, rtrange = c(6000, 6060)))
    expect_error(getEIC(xraw, rtrange = c(6000, 6060)))
    ## mzrange with more than 2 columns
    expect_error(profEIC(xraw, mzrange = c(1, 3, 5)))
    expect_error(getEIC(xraw, mzrange = c(1, 3, 5)))
    ## rtrange with more than 2 columns
    expect_error(profEIC(xraw, rtrange = c(2, 4, 6)))
    expect_error(getEIC(xraw, rtrange = c(2, 4, 6)))
})


test_that("plotEIC,xcmsRaw works", {
    ## file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    ## xraw <- xcmsRaw(file)
    xraw <- deepCopy(faahko_xr_1)
    e <- getEIC(xraw, rtrange=cbind(3000,3500), mzrange=cbind(200,201))
    plot(e)
})

test_that("findPeaks.centWave,xcmsRaw works with scanrange", {
    ## Without sub-setting!
    xraw <- deepCopy(faahko_xr_1)
    res_1 <- findPeaks.centWave(xraw, noise = 10000)
    ## res_2 <- xcms:::.findPeaks.centWave_orig(xraw, noise = 10000)
    ## checkIdentical(res_1, res_2)

    scnr <- c(90, 345)
    res_1 <- findPeaks.centWave(xraw, scanrange = scnr, noise = 5000)
    ## res_2 <- xcms:::.findPeaks.centWave_orig(xraw, scanrange = scnr,
    ##                                          noise = 5000)
    ## checkIdentical(res_1, res_2)

    ## Compare with do_
    xsub <- xraw[90:345]
    res_3 <- do_findChromPeaks_centWave(mz = xsub@env$mz,
                                        int = xsub@env$intensity,
                                        scantime = xsub@scantime,
                                        noise = 5000,
                                        valsPerSpect = diff(c(xsub@scanindex,
                                                              length(xsub@env$mz))))
    expect_identical(res_3, res_1@.Data)

    scnr <- c(1, 400)
    res_1 <- findPeaks.centWave(xraw, scanrange = scnr, noise = 5000)
    ## res_2 <- xcms:::.findPeaks.centWave_orig(xraw, scanrange = scnr,
    ##                                          noise = 5000)
    ## checkIdentical(res_1, res_2)
})

test_that("findPeaks.matchedFilter,xcmsRaw works with scanrange", {
    scnr <- c(90, 50000)
    xraw <- deepCopy(faahko_xr_1)
    res_1 <- findPeaks.matchedFilter(xraw, scanrange = scnr)
    ## suppressWarnings(
    ##     res_2 <- xcms:::findPeaks.matchedFilter_orig(xraw, scanrange = scnr)
    ## )
    expect_warning(xsub <- xraw[90:50000])
    res_3 <- findPeaks.matchedFilter(xsub)
    ## checkIdentical(res_1, res_2)
    expect_identical(res_1, res_3)
})

test_that("findPeaks.massifquant,xcmsRaw works with scanrange", {
    ## Compare passing the scanrange with performing the search on a
    ## pre-subsetted object.
    xraw <- deepCopy(faahko_xr_1)
    scnr <- c(90, 150)
    res_1 <- findPeaks.massifquant(xraw, scanrange = scnr)
    xsub <- xraw[scnr[1]:scnr[2]]
    res_2 <- findPeaks.massifquant(xsub)
    expect_identical(res_1, res_2)
    ## Same with "withWave"
    scnr <- c(90, 200)
    res_1 <- findPeaks.massifquant(xraw, scanrange = scnr, withWave = 1)
    xsub <- xraw[scnr[1]:scnr[2]]
    res_2 <- findPeaks.massifquant(xsub, withWave = 1)
    expect_identical(res_1, res_2)
})

test_that("[,xcmsRaw works", {
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xraw <- xcmsRaw(file, profstep = 0)
    ## Get scans 1:10
    xsub <- xraw[1:10]
    expect_identical(xraw@scantime[1:10], xsub@scantime)
    expect_identical(xraw@scanindex[1:10], xsub@scanindex[1:10])
    expect_identical(xraw@env$mz[1:xraw@scanindex[11]], xsub@env$mz)
    expect_identical(xraw@env$intensity[1:xraw@scanindex[11]], xsub@env$intensity)
    ## Check if mz is sorted:
    ## Check using logical
    bm <- rep(FALSE, length(xraw@scanindex))
    bm[1:10] <- TRUE
    xsub_2 <- xraw[bm]
    expect_equal(xsub, xsub_2)
    ## Get none.
    xempty <- xraw[rep(FALSE, length(xraw@scanindex))]
    expect_true(length(xempty@scanindex) == 0)
    ## Get the full one:
    xsub <- xraw[1:length(xraw@scanindex)]
    expect_equal(xsub@env$mz, xraw@env$mz)
    expect_equal(xsub@env$intensity, xraw@env$intensity)
    ## Get some scans in the middle somewhere
    i <- c(5, 99, 317)
    xsub <- xraw[i]
    expect_identical(xsub@scantime, xraw@scantime[i])
    ## scanindex:
    vps <- diff(c(xraw@scanindex, length(xraw@env$mz)))
    scnidx <- xcms:::valueCount2ScanIndex(vps[i])
    expect_equal(scnidx, xsub@scanindex)
    whichIdx <- c(((xraw@scanindex[5] + 1):xraw@scanindex[6]),
    ((xraw@scanindex[99] + 1):xraw@scanindex[100]),
    ((xraw@scanindex[317] + 1):xraw@scanindex[318]))
    expect_identical(xsub@env$mz, xraw@env$mz[whichIdx])
    expect_identical(xsub@env$intensity, xraw@env$intensity[whichIdx])
    ## Finally check that we get the same object by calling getXcmsRaw from an
    ## xcmsSet and by subsetting or using xcmsRaw with scanrange:
    xraw <- xcmsRaw(file, profstep = 20)
    xraw_sub <- xraw[5:100]
    expect_warning(xset <- xcmsSet(file, scanrange = c(5, 100), step = 20))
    xraw_xset <- getXcmsRaw(xset)
    expect_equal(scanrange(xraw_sub), scanrange(xraw_xset))
    ## Load the object using xcmsRaw and scanrange
    xraw_2 <- xcmsRaw(file, scanrange = c(5, 100), profstep = 20)
    expect_equal(scanrange(xraw_sub), scanrange(xraw_2))
    ## Compare objects
    expect_equal(xraw_sub, xraw_xset)
    expect_equal(xraw_sub, xraw_2)
    expect_equal(xraw_xset, xraw_2)
})

test_that("deepCopy,xcmsRaw works", {
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xraw <- xcmsRaw(file, profstep = 0)
    xrawCopy <- deepCopy(xraw)
    expect_true(all(xraw@env$mz == xrawCopy@env$mz))
    xrawCopy@env$mz <-     xrawCopy@env$mz+1
    expect_true(all(xraw@env$mz != xrawCopy@env$mz))
})
