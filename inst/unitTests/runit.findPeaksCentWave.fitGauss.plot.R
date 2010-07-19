
test.findPeaksCentWave <- function() {
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xraw <- xcmsRaw(file)
    p <- findPeaks.centWave(xraw, fitgauss=T, verbose=T, sleep=0.001)
}
