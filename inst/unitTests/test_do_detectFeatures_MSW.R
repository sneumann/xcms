############################################################
## do_detectFeatures_MSW tests

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
        checkEquals(a, b)
    }
}
