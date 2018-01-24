## Comparing the old profStepPad<- with the new one.
## library(faahKO)
## fs <- c(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
##         system.file('cdf/KO/ko16.CDF', package = "faahKO"),
##         system.file('cdf/KO/ko18.CDF', package = "faahKO"),
##         system.file('cdf/KO/ko19.CDF', package = "faahKO"))

## library(msdata)
## mzf <- c(system.file("microtofq/MM14.mzML", package = "msdata"),
##          system.file("microtofq/MM8.mzML", package = "msdata"))
## f <- msdata::proteomics(full.names = TRUE, pattern = "TMT_Erwinia")

dontrun_test_profStepPad <- function() {
    xr_1 <- deepCopy(faahko_xr_1)
    xr_2 <- deepCopy(faahko_xr_1)

    steps <- c(0.1, 0.3, 0.34, 1, 2, 2.1, 2.13)
    for (i in 1:length(steps)) {
        xcms:::profStepPad(xr_1) <- steps[i]
        xcms:::profStepPadOld(xr_2) <- steps[i]
        checkEquals(xr_1@env$profile, xr_2@env$profile)
    }
    ## I get differences here due to the difference between the old profBin
    ## method and the new (correct) implementation.
}
