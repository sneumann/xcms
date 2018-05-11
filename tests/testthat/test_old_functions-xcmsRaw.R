test_that("xcmsRaw on MS1 asking for MS2 doesn't fail", {
    filename <- system.file('microtofq/MM14.mzdata', package = "msdata")
    ## This file has no MS/MS data at all, but should not fail
    expect_warning(x1 <- xcmsRaw(filename, includeMSn=TRUE, profstep = 0))
})

test_that("xcmsRaw with multiple MS levels works", {
    filename <- system.file('iontrap/extracted.mzData', package = "msdata")
    x1 <- xcmsRaw(filename, includeMSn=TRUE, profstep = 0)
    expect_warning(x2 <- xcmsRaw(filename, includeMSn=TRUE, mslevel=2,
                                 profstep = 0))
    expect_warning(x3 <- xcmsRaw(filename, includeMSn=TRUE, mslevel=3,
                                 profstep = 0))

    expect_true(length(x1@env$msnMz) == length(x2@env$mz) + length(x3@env$mz))
    expect_true(all(x1@msnLevel[1:6]==2))
    expect_true(all(x1@msnScanindex[1:6] == x2@scanindex[1:6]))
    expect_equal(nrow(getMsnScan(x1, scan=1)), 278)

    ## This would fail, since mslevel=2 above seems to use split(),
    ## which does drop MSn information ?
    ## expect_equalNumeric(nrow(xcms:::getMsnScan(x2, scan=1)), 278)

})

test_that("msn2xcmsRaw works", {
    msnfile <- system.file("microtofq/MSMSpos20_6.mzML", package = "msdata")
    xrmsn <- xcmsRaw(msnfile, includeMSn=TRUE)
    xr <- msn2xcmsRaw(xrmsn)

    expect_equal(length(xr@env$mz), 3132)
    expect_equal(length(xr@env$intensity), 3132)
    ## In reality, it seems there are 1612 MS2 spectra in the file, just that
    ## 1121 have a peaksCount > 0
    expect_equal(length(xr@scantime), 1121)
})

test_that("xcmsRaw works with scanrange", {
    fs <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    ## Use xcmsRaw with scanrange and check if results are identical to
    ## later subsetting.
    x_1 <- xcmsRaw(fs, scanrange = c(90, 345), profstep = 0)
    x_2 <- xcmsRaw(fs, profstep = 0)
    x_2_sub <- x_2[90:345]
    expect_true(length(x_1@scantime) < length(x_2@scantime))
    expect_equal(x_1@scantime, x_2_sub@scantime)
    length(x_1@acquisitionNum)
    length(x_2@acquisitionNum)
    length(x_2_sub@acquisitionNum)
    expect_true(length(x_1@acquisitionNum) < length(x_2@acquisitionNum))
    expect_equal(x_1@acquisitionNum, x_2_sub@acquisitionNum)
    x_1@scanrange
    x_2@scanrange
    x_2_sub@scanrange
    expect_equal(x_1@scanrange, x_2_sub@scanrange)
    expect_equal(x_1, x_2_sub)
})

test_that("rawMat,xcmsRaw works", {
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xraw <- xcmsRaw(file, profstep = 0)
    ## That's the old test; pretty arbitrary
    ##rawmat <- rawMat(xraw, mzrange = c(200,300), rtrange = c(2500,3000))
    ##expect_equalNumeric(nrow(rawmat),31770)

    ## New evaluations:
    ## Check with manually generated content:
    scnt <- xraw@scantime
    scnrange <- range(which((scnt >= 2500) & (scnt <= 3000)))
    ## Just to ensure: are the scantimes within our range?
    expect_true(all((scnt[scnrange] >= 2500 & scnt[scnrange] <= 3000)))
    ## Now get the values falling into this scan range:
    ## scanindex records the start index for each scan/spectrum -1.
    scanStart <- xraw@scanindex[scnrange[1]] + 1
    ## To get the end index we add also the total number of scans and
    ## pick the start index of the next scan outside our scnrange
    scanEnd <- c(xraw@scanindex, length(xraw@env$mz))[scnrange[2] + 1]
    mzVals <- xraw@env$mz[scanStart:scanEnd]
    ## Check if that's what we get from the function.
    res_scnr <- rawMat(xraw, scanrange = scnrange)
    expect_equal(res_scnr[, "mz"], mzVals)
    ## Is that the same if we used the rtrange instead?
    res_rtr <- rawMat(xraw, rtrange = c(2500, 3000))
    expect_identical(res_rtr[, "mz"], mzVals)
    ## Hm. Actually we have a discrepancy here:
    expect_identical(res_scnr, res_rtr)
    ## Compare the retention time for the subset: get the retention time
    ## for each scan and "rep" it the number of values per scan.
    valsPerSpect <- diff(c(xraw@scanindex, length(xraw@env$mz)))
    rtimes <- rep(scnt[scnrange[1]:scnrange[2]],
                  valsPerSpect[scnrange[1]:scnrange[2]])
    expect_identical(res_scnr[, "time"], rtimes)
    ## This failing test was reported as issue #58 on github.

    ## Add also an mzrange:
    mzr <- c(300, 330)
    mzInRange <- xraw@env$mz >= mzr[1] & xraw@env$mz <= mzr[2]
    res_mzr <- rawMat(xraw, mzrange = mzr)
    expect_equal(res_mzr[, "mz"], xraw@env$mz[mzInRange])
    ## Combine with scnrange:
    scnInRange <- rep(FALSE, length(mzInRange))
    scnInRange[scanStart:scanEnd] <- TRUE
    res_mzr_scnr <- rawMat(xraw, scanrange = scnrange, mzrange = mzr)
    expect_identical(res_mzr_scnr[, "mz"], xraw@env$mz[mzInRange & scnInRange])
    ## And the retention time.
    allRtimes <- rep(xraw@scantime, diff(c(xraw@scanindex, length(xraw@env$mz))))
    expect_identical(res_mzr_scnr[, "time"], allRtimes[mzInRange & scnInRange])

    ## And a second scan range:
    scnrange <- c(13, 54)
    scanStart <- xraw@scanindex[scnrange[1]] + 1
    scanEnd <- c(xraw@scanindex, length(xraw@env$mz))[scnrange[2] + 1]
    mzVals <- xraw@env$mz[scanStart:scanEnd]
    res_scnr <- rawMat(xraw, scanrange = scnrange)
    expect_identical(res_scnr[, "mz"], mzVals)
    valsPerSpect <- diff(c(xraw@scanindex, length(xraw@env$mz)))
    rtimes <- rep(scnt[scnrange[1]:scnrange[2]],
                  valsPerSpect[scnrange[1]:scnrange[2]])
    expect_identical(res_scnr[, "time"], rtimes)

    ## Test .rawMat directly. This is to evaluate potential problems in
    ## .getChromPeakData/fillChromPeaks. We're calling .rawMat ONLY based
    ## on rtrange and mzrange. scanrange can/could cause a problem being
    ## -Inf, Inf.
    mz <- xraw@env$mz
    int <- xraw@env$intensity
    rt <- xraw@scantime
    valsPerSpect <- diff(c(xraw@scanindex, length(xraw@env$mz)))
    ## rtrange outside range of rt
    res <- xcms:::.rawMat(mz, int, rt, valsPerSpect, rtrange = c(12, 30))
    expect_true(is.matrix(res))
    expect_true(nrow(res) == 0)
    ## rtrange partially outside
    res <- xcms:::.rawMat(mz, int, rt, valsPerSpect, rtrange = c(12, 2502))
    expect_true(is.matrix(res))
    expect_true(nrow(res) > 0)
    ## mzrange outside range of mz
    res <- xcms:::.rawMat(mz, int, rt, valsPerSpect, mzrange = c(20, 30))
    expect_true(is.matrix(res))
    expect_true(nrow(res) == 0)
})

test_that("xcmsRaw with scanrange works", {
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xraw <- xcmsRaw(file, profstep = 0)
    expect_equal(length(xraw@scantime), 1278)
    xraw <- xcmsRaw(file, scanrange=c(1,1), profstep = 0)
    expect_equal(length(xraw@scantime), 1)
    xraw <- xcmsRaw(file, scanrange=c(1278,1278), profstep = 0)
    expect_equal(length(xraw@scantime), 1)
    xraw <- xcmsRaw(file, scanrange=c(100,200), profstep = 0)
    expect_equal(length(xraw@scantime), 101)
    xraw <- xcmsRaw(file, scanrange=c(100,199), profstep = 0)
    expect_equal(length(xraw@scantime), 100)
})

test_that("split.xcmsRaw works", {
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xraw <- xcmsRaw(file, profstep = 0)
    ## even
    xrl <- split(xraw, f = xraw@scanindex%%2)
    expect_equal(length(xrl), 2)
    expect_true(length(xraw@scanindex) == length(xrl[[1]]@scanindex)
              + length(xrl[[2]]@scanindex))
    expect_true(length(xraw@env$mz) == length(xrl[[1]]@env$mz)
              + length(xrl[[2]]@env$mz))
    ## odd
    xrl <- split(xraw, f = (xraw@scanindex+1)%%2)
    expect_equal(length(xrl), 2)
    expect_true(length(xraw@scanindex) == length(xrl[[1]]@scanindex) +
              length(xrl[[2]]@scanindex))
    expect_true(length(xraw@env$mz) == length(xrl[[1]]@env$mz) +
              length(xrl[[2]]@env$mz))
    ## first
    xrl <- split(xraw, f=c(1,rep(2,length(xraw@scanindex)-1)))
    expect_equal(length(xrl), 2)
    ## last
    xrl <- split(xraw, f=c(rep(1,length(xraw@scanindex)-1), 2))
    expect_equal(length(xrl), 2)
    ## none
    xrl <- split(xraw, f=rep(1,length(xraw@scanindex)))
    expect_equal(length(xrl), 1)
    xraw2 <- xrl[[1]]
    expect_true(length(xraw@scanindex) == length(xrl[[1]]@scanindex))
    expect_true(all(xraw@scanindex == xraw2@scanindex))
})
