############################################################
## do_detectFeatures_MSW tests

library(msdata)
mzf <- c(system.file("microtofq/MM14.mzML", package = "msdata"),
         system.file("microtofq/MM8.mzML", package = "msdata"))

xraw <- xcmsRaw(mzf[1], profstep = 0)

test_do_detectFeatures_MSW <- function() {
    feats1 <- xcms:::do_detectFeatures_MSW(xraw@env$intensity,
                                           xraw@env$mz,
                                           snthresh = 100)
    feats2 <- xcms:::do_detectFeatures_MSW(xraw@env$intensity,
                                           xraw@env$mz,
                                           snthresh = 50)
    checkTrue(nrow(feats2) > nrow(feats1))
}

test_detectFeatures_MSW <- function() {
    library(MSnbase)
    od <- readMSData2(mzf)
    ## Restrict to first spectrum
    od1 <- od[1]
    sp1 <- od[[1]]
    res_1 <- do_detectFeatures_MSW(mz = mz(sp1), int = intensity(sp1))
    mp <- MSWParam()
    res_2 <- detectFeatures(od1, param = mp)
    checkEquals(res_1, res_2)
    ## Changing settings.
    snthresh(mp) <- 1
    nearbyPeak(mp) <- FALSE
    res_1 <- do_detectFeatures_MSW(mz = mz(sp1), int = intensity(sp1),
                                   snthresh = 1, nearbyPeak = FALSE)
    res_2 <- detectFeatures(od1, param = mp)
    checkEquals(res_1, res_2[[1]][, colnames(res_1)])
    peakThr(mp) <- 200
    res_1 <- do_detectFeatures_MSW(mz = mz(sp1), int = intensity(sp1),
                                   snthresh = 1, nearbyPeak = FALSE,
                                   peakThr = 200)
    res_2 <- detectFeatures(od1, param = mp)
    checkEquals(res_1, res_2[[1]][, colnames(res_1)])
    addParams(mp) <- list(forder = 2)
    res_3 <- do_detectFeatures_MSW(mz = mz(sp1), int = intensity(sp1),
                                   snthresh = 1, nearbyPeak = FALSE,
                                   peakThr = 200, forder = 2)
    res_4 <- detectFeatures(od1, param = mp)
    checkEquals(res_3, res_4[[1]][, colnames(res_3)])
    addParams(mp) <- list(forder = 2, dorder = 1)
    res_3 <- do_detectFeatures_MSW(mz = mz(sp1), int = intensity(sp1),
                                   snthresh = 1, nearbyPeak = FALSE,
                                   peakThr = 200, forder = 2, dorder = 1)
    res_4 <- detectFeatures(od1, param = mp)
    checkEquals(res_3, res_4[[1]][, colnames(res_3)])
}


############################################################
## Test the implementation of the "do" function
dontrun_test_do_detectFeatures_MSW_impl <- function() {

    library(xcms)
    library(RUnit)

    library(faahKO)
    fs <- c(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
            system.file('cdf/KO/ko16.CDF', package = "faahKO"),
            system.file('cdf/KO/ko18.CDF', package = "faahKO"),
            system.file('cdf/KO/ko19.CDF', package = "faahKO"))
    i <- 1
    for (i in 1:length(fs)) {
        cat("============================================================\n")
        cat("|   file", i, ":", fs[i], "\n")
        cat("------------------------------------------------------------")

        xr <- xcmsRaw(fs[i])
        mz <- xr@env$mz
        int <- xr@env$intensity

        ## Default
        a <- findPeaks.MSW(xr)
        b <- xcms:::.MSW_orig(int, mz)
        checkEquals(a@.Data, b)

        ## Slightly modified one:
        c <- xcms:::.MSW(int, mz)
        checkEquals(a@.Data, c)
    }
}

dontrun_benchmark_MSW <- function() {
    ## Actually, we don't expect much of a speed improvement here.
    library(microbenchmark)
    microbenchmark(xcms:::.MSW(xraw@env$intensity, xraw@env$mz),
                   xcms:::.MSW_orig(xraw@env$intensity, xraw@env$mz),
                   times = 4)
}
