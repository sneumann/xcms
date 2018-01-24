test.xcmsRawScanrange.none <- function() {
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xraw <- xcmsRaw(file, profstep = 0)
    checkEquals(length(xraw@scantime), 1278)
}
test.xcmsRawScanrange.first <- function() {
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xraw <- xcmsRaw(file, scanrange=c(1,1), profstep = 0)
    checkEquals(length(xraw@scantime), 1)
}
test.xcmsRawScanrange.last <- function() {
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xraw <- xcmsRaw(file, scanrange=c(1278,1278), profstep = 0)
    checkEquals(length(xraw@scantime), 1)
}
test.xcmsRawScanrange.100 <- function() {
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xraw <- xcmsRaw(file, scanrange=c(100,200), profstep = 0)
    checkEquals(length(xraw@scantime), 100)
}
test.xcmsRawScanrange.100 <- function() {
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xraw <- xcmsRaw(file, scanrange=c(100,199), profstep = 0)
    checkEquals(length(xraw@scantime), 100)
}
