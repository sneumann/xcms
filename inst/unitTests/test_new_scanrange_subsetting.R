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
    res_2 <- xcms:::.findPeaks.centWave_scanrange(xraw, scanrange = scnr)
    ## checkIdentical(res_1, res_2)
    ## Does NOT work: scalerange e.g. depends on scantime.

    scnr <- c(1, 400)
    res_1 <- findPeaks.centWave(xraw, scanrange = scnr)
    res_2 <- xcms:::.findPeaks.centWave_scanrange(xraw, scanrange = scnr)
    ## checkIdentical(res_1, res_2)
    ## Does NOT work: scalerange e.g. depends on scantime.
}

test_scanrange_matchedFilter <- function() {
    ## Without sub-setting
    res_1 <- findPeaks.matchedFilter(xraw)
    res_2 <- xcms:::.findPeaks.matchedFilter_scanrange(xraw)
    checkIdentical(res_1, res_2)

    scnr <- c(90, 345)
    res_1 <- findPeaks.matchedFilter(xraw, scanrange = scnr)
    res_2 <- xcms:::.findPeaks.matchedFilter_scanrange(xraw, scanrange = scnr)
    ## checkIdentical(res_1, res_2)
    ## WHY??
}

test_scanrange_massifquant <- function() {
}
