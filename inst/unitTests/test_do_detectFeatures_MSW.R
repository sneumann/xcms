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
