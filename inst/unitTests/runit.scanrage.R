test.xcmsRawScanrange.none <- function() {
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xraw <- xcmsRaw(file)
    checkEquals(length(xraw@scantime), 1278)
}
test.xcmsRawScanrange.first <- function() {
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xraw <- xcmsRaw(file, scanrange=c(1,1))
    checkEquals(length(xraw@scantime), 1)
}
test.xcmsRawScanrange.last <- function() {
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xraw <- xcmsRaw(file, scanrange=c(1278,1278))
    checkEquals(length(xraw@scantime), 1)
}
test.xcmsRawScanrange.100 <- function() {
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xraw <- xcmsRaw(file, scanrange=c(100,200))
    checkEquals(length(xraw@scantime), 100)
}
test.xcmsRawScanrange.100 <- function() {
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xraw <- xcmsRaw(file, scanrange=c(100,199))
    checkEquals(length(xraw@scantime), 100)
}
