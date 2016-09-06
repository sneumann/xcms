
test.findPeaksCentWaveWithIsotopeROIs <- function() {
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xr <- xcmsRaw(file)
    p1 <- findPeaks.centWave(xr)
    
    addROIsParams <- list()
    addROIsParams$addNewROIs <- TRUE
    addROIsParams$addNewIsotopeROIs <- TRUE
    addROIsParams$addNewAdductROIs  <- FALSE
    addROIsParams$snthreshOfGeneratedROIs <- 6.25
    addROIsParams$maxcharge  <- 3
    addROIsParams$maxiso     <- 5
    addROIsParams$mzIntervalExtension <- TRUE
    
    p2 <- findPeaks.centWaveWithPredictedIsotopeROIs(
      xr, 
      ROI.list=list(), 
      addROIsParams = addROIsParams
    )
    
    checkTrue(nrow(p1)<nrow(p2))
    
}

