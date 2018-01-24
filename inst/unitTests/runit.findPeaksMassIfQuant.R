
test.findPeaksMassIfQuant <- function() {
    file <- system.file('microtofq/MM14.mzdata', package = "msdata")
    xraw <- xcmsRaw(file, profstep = 0)
    p <- findPeaks(xraw, method = "massifquant")
}
