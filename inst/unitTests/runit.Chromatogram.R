## Unit tests related to the Chromatogram class.
library(xcms)
library(RUnit)

test_Chromatogram_class <- function() {
    ch <- new("Chromatogram")
    ch@mz <- 3
    checkException(validObject(ch))
    int <- rnorm(100, mean = 200, sd = 2)
    rt <- rnorm(100, mean = 300, sd = 3)
    ## check exceptions:
    checkException(xcms:::Chromatogram(intensity = int))
    checkException(xcms:::Chromatogram(intensity = int, rtime = rt))
    rt <- sort(rt)
    ch <- xcms:::Chromatogram(intensity = int, rtime = rt)
    checkEquals(rtime(ch), rt)
    checkEquals(intensity(ch), int)
    checkException(xcms:::Chromatogram(aggregationFun = "other"))
    ch@aggregationFun <- "max"
    checkTrue(validObject(ch))
    checkEquals(aggregationFun(ch), "max")
    ch@aggregationFun <- "sum"
    checkTrue(validObject(ch))
    checkEquals(aggregationFun(ch), "sum")
    ch@aggregationFun <- "mean"
    checkTrue(validObject(ch))
    checkEquals(aggregationFun(ch), "mean")
    ch@aggregationFun <- "min"
    checkTrue(validObject(ch))
    checkEquals(aggregationFun(ch), "min")
    ch@fromFile <- 3L
    checkTrue(validObject(ch))
    checkEquals(fromFile(ch), 3L)
    checkEquals(length(ch), length(rt))
    ## as.data.frame
    df <- as.data.frame(ch)
    checkEquals(df, data.frame(rtime = rt, intensity = int))
    ch <- xcms:::Chromatogram(mz = c(1, 3))
    checkEquals(ch@mz, c(1, 3))
    checkEquals(mzrange(ch), c(1, 3))
    checkEquals(mzrange(ch, filter = TRUE), c(0, 0))
    ch <- xcms:::Chromatogram(filterMz = c(1, 3))
    checkEquals(ch@filterMz, c(1, 3))
    checkEquals(mzrange(ch, filter = TRUE), c(1, 3))
    checkEquals(mzrange(ch, filter = FALSE), c(0, 0))
}
