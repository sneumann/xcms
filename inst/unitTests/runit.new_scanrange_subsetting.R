############################################################
## This test file is supposed to be removed once the subsetting
## of xcmsRaw objects in the findPeaks methods prior to feature
## detection has been tested and is working as expected.
## issue #62
## library(faahKO)
fs <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
## xraw <- xcmsRaw(fs, profstep = 0)
xraw <- deepCopy(faahko_xr_1)

test_scanrange_centWave <- function() {
    ## Without sub-setting
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
    checkIdentical(res_3, res_1@.Data)

    scnr <- c(1, 400)
    res_1 <- findPeaks.centWave(xraw, scanrange = scnr, noise = 5000)
    ## res_2 <- xcms:::.findPeaks.centWave_orig(xraw, scanrange = scnr,
    ##                                          noise = 5000)
    ## checkIdentical(res_1, res_2)
}

test_scanrange_matchedFilter <- function() {
    scnr <- c(90, 50000)
    res_1 <- findPeaks.matchedFilter(xraw, scanrange = scnr)
    ## suppressWarnings(
    ##     res_2 <- xcms:::findPeaks.matchedFilter_orig(xraw, scanrange = scnr)
    ## )
    suppressWarnings(
        xsub <- xraw[90:50000]
    )
    res_3 <- findPeaks.matchedFilter(xsub)
    ## checkIdentical(res_1, res_2)
    checkIdentical(res_1, res_3)
}

test_scanrange_massifquant <- function() {
    ## Compare passing the scanrange with performing the search on a
    ## pre-subsetted object.
    scnr <- c(90, 150)
    res_1 <- findPeaks.massifquant(xraw, scanrange = scnr)
    xsub <- xraw[scnr[1]:scnr[2]]
    res_2 <- findPeaks.massifquant(xsub)
    checkIdentical(res_1, res_2)
    ## Same with "withWave"
    scnr <- c(90, 200)
    res_1 <- findPeaks.massifquant(xraw, scanrange = scnr, withWave = 1)
    xsub <- xraw[scnr[1]:scnr[2]]
    res_2 <- findPeaks.massifquant(xsub, withWave = 1)
    checkIdentical(res_1, res_2)
}

test_scanrange_xcmsRaw <- function() {
    ## Use xcmsRaw with scanrange and check if results are identical to
    ## later subsetting.
    x_1 <- xcmsRaw(fs, scanrange = c(90, 345), profstep = 0)
    x_2 <- xcmsRaw(fs, profstep = 0)
    x_2_sub <- x_2[90:345]
    checkTrue(length(x_1@scantime) < length(x_2@scantime))
    checkEquals(x_1@scantime, x_2_sub@scantime)
    length(x_1@acquisitionNum)
    length(x_2@acquisitionNum)
    length(x_2_sub@acquisitionNum)
    checkTrue(length(x_1@acquisitionNum) < length(x_2@acquisitionNum))
    checkEquals(x_1@acquisitionNum, x_2_sub@acquisitionNum)
    x_1@scanrange
    x_2@scanrange
    x_2_sub@scanrange
    checkEquals(x_1@scanrange, x_2_sub@scanrange)

    checkEquals(x_1, x_2_sub)
}
