############################################################
## XCMSnSet related tests.
library(RUnit)

library(faahKO)
fs <- c(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
        system.file('cdf/KO/ko16.CDF', package = "faahKO"),
        system.file('cdf/KO/ko18.CDF', package = "faahKO"),
        system.file('cdf/KO/ko19.CDF', package = "faahKO"))
library(msdata)
f <- msdata::proteomics(full.names = TRUE, pattern = "TMT_Erwinia")
mzf <- c(system.file("microtofq/MM14.mzML", package = "msdata"),
         system.file("microtofq/MM8.mzML", package = "msdata"))

test_XCMSnSet_class <- function() {
    ## Basic contructor.
    xs <- new("XCMSnSet")
    xs@.processHistory <- list("a")
    checkException(validObject(xs))
    xs@.processHistory <- list(xcms:::ProcessHistory())
    checkTrue(validObject(xs))
}


