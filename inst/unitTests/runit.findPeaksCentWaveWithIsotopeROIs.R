
test.addPredictedIsotopeFeatures <- function() {
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xr <- xcmsRaw(file, profstep = 0)
    p1 <- findPeaks.centWave(xr, verbose.columns = TRUE, noise = 10000)

    p2 <- findPeaks.addPredictedIsotopeFeatures(
        object = xr, xcmsPeaks = p1, noise = 10000)

    checkTrue(nrow(p1) < nrow(p2))
    ## Now the same with the new modified centWave:
    options(originalCentWave = FALSE)
    p1_2 <- findPeaks.centWave(xr, verbose.columns = TRUE, noise = 10000)
    nrow(p1)
    nrow(p1_2)
    checkEquals(p1, p1_2)
    
    p2_2 <- findPeaks.addPredictedIsotopeFeatures(
        object = xr, xcmsPeaks = p1_2, noise = 10000)

    options(originalCentWave = TRUE)
    
}
