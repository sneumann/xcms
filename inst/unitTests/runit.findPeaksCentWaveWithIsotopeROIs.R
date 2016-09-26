
test.addPredictedIsotopeFeatures <- function() {
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xr <- xcmsRaw(file)
    p1 <- findPeaks.centWave(xr, verbose.columns=TRUE)
    
    p2 <- findPeaks.addPredictedIsotopeFeatures(
      object = xr, xcmsPeaks=p1
    )
    
    checkTrue(nrow(p1)<nrow(p2))
}
