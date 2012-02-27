
test.findPeaksMassIfQuant <- function() {
    file <- system.file('microtofq/MM14.mzdata', package = "msdata")
    xraw <- xcmsRaw(file)
    p <- findPeaks(xraw, method="massifquant")
}
