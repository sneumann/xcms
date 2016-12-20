
test.findPeaksCentWave <- function() {
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xraw <- xcmsRaw(file, profstep = 0)
    p <- findPeaks.centWave(xraw, fitgauss = TRUE, verbose = TRUE, sleep=0.001,
                            noise = 10000)
}
