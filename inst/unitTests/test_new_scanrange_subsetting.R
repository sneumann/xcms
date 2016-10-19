############################################################
## This test file is supposed to be removed once the subsetting
## of xcmsRaw objects in the findPeaks methods prior to feature
## detection has been tested and is working as expected.
## issue #62
library(faahKO)
fs <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
xraw <- xcmsRaw(fs, profstep = 0)

test_scanrange_centWave <- function() {
    ## Without sub-setting
    res_1 <- findPeaks.centWave(xraw)
    res_2 <- xcms:::.findPeaks.centWave_scanrange(xraw)
    checkIdentical(res_1, res_2)

    scnr <- c(90, 345)
    res_1 <- findPeaks.centWave(xraw, scanrange = scnr)
    ## res_2 <- xcms:::.findPeaks.centWave_scanrange(xraw, scanrange = scnr)
    ## checkIdentical(res_1, res_2)
    ## Does NOT work: scalerange e.g. depends on scantime.
    xsub <- xraw[90:345]
    res_1 <- findPeaks.centWave(xsub)
    ## res_2 <- xcms:::.findPeaks.centWave_scanrange(xraw, scanrange = scnr)

    scnr <- c(1, 400)
    res_1 <- findPeaks.centWave(xraw, scanrange = scnr)
    res_2 <- xcms:::.findPeaks.centWave_scanrange(xraw, scanrange = scnr)
    ## checkIdentical(res_1, res_2)
    ## Does NOT work: scalerange e.g. depends on scantime.
}

test_scanrange_matchedFilter <- function() {
    scnr <- c(90, 50000)
    res_1 <- findPeaks.matchedFilter(xraw, scanrange = scnr)
    res_2 <- xcms:::findPeaks.matchedFilter_orig(xraw, scanrange = scnr)
    xsub <- xraw[90:50000]
    res_3 <- findPeaks.matchedFilter(xsub)
    checkIdentical(res_1, res_2)
    checkIdentical(res_1, res_3)
}

test_scanrange_massifquant <- function() {
    ## Compare passing the scanrange with performing the search on a
    ## pre-subsetted object.
    scnr <- c(90, 345)
    res_1 <- findPeaks.massifquant(xraw, scanrange = scnr)
    xsub <- xraw[scnr[1]:scnr[2]]
    res_2 <- findPeaks.massifquant(xsub)
    checkIdentical(res_1, res_2)
    ## Same with "withWave"
    res_1 <- findPeaks.massifquant(xraw, scanrange = scnr, withWave = 1)
    res_2 <- findPeaks.massifquant(xsub, withWave = 1)
    checkIdentical(res_1, res_2)
}
