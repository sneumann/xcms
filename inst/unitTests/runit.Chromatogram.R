## Unit tests related to the Chromatogram class.
library(xcms)
library(RUnit)

test_Chromatogram_class <- function() {
    ch <- new("Chromatogram")
    ch@mzrange <- 3
    checkException(validObject(ch))
    int <- rnorm(100, mean = 200, sd = 2)
    rt <- rnorm(100, mean = 300, sd = 3)
    ## check exceptions:
    checkException(Chromatogram(intensity = int))
    checkException(Chromatogram(intensity = int, rtime = rt))
    rt <- sort(rt)
    ch <- Chromatogram(intensity = int, rtime = rt)
    checkEquals(rtime(ch), rt)
    checkEquals(intensity(ch), int)
    checkException(Chromatogram(aggregationFun = "other"))
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
    ch <- Chromatogram(mzrange = c(1, 3))
    checkEquals(ch@mzrange, c(1, 3))
    checkEquals(mzrange(ch), c(1, 3))
    checkEquals(mzrange(ch, filter = TRUE), c(0, 0))
    ch <- Chromatogram(filterMzrange = c(1, 3))
    checkEquals(ch@filterMzrange, c(1, 3))
    checkEquals(mzrange(ch, filter = TRUE), c(1, 3))
    checkEquals(mzrange(ch, filter = FALSE), c(0, 0))
}
