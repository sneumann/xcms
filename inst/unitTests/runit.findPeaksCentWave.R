
test.findPeaksCentWaveResort <- function() {
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xr <- xcmsRaw(file)
    p1 <- findPeaks.centWave(xr)

    scan <- getScan(xr,1)
    o <- order(scan[,"mz"], decreasing=TRUE)

    xr@env$mz[1:length(o)] <- scan[o, "mz"] 
    xr@env$intensity[1:length(o)] <- scan[o, "intensity"] 

    p2 <- findPeaks.centWave(xr)

    checkTrue(all(p1==p2))
    
}

