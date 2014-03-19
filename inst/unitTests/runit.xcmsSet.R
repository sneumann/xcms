test.findPeaksCentWave <- function() {
    file <- system.file('microtofq/MM14.mzdata', package = "msdata")

    xset1 <-  xcmsSet(files=file, method="centWave", peakwidth=c(5,12))
    xset2 <-  xcmsSet(files=file, method="centWave", peakwidth=c(5,12), scanrange=c(1,112))
    xset3 <-  xcmsSet(files=file, method="centWave", peakwidth=c(5,12), scanrange=c(1,80))

    checkTrue (nrow((peaks(xset1)@.Data)) == nrow((peaks(xset2)@.Data)))
    checkTrue (nrow((peaks(xset1)@.Data))  > nrow((peaks(xset3)@.Data)))

}

test.findPeaksMatchedFilter <- function() {
    file <- system.file('microtofq/MM14.mzdata', package = "msdata")

    xset1 <- xcmsSet(files=file, method="matchedFilter", fwhm=10)
    xset2 <- xcmsSet(files=file, method="matchedFilter", fwhm=10, scanrange=c(1,112))
    xset3 <- xcmsSet(files=file, method="matchedFilter", fwhm=10, scanrange=c(1,80))

    checkTrue (nrow((peaks(xset1)@.Data)) == nrow((peaks(xset2)@.Data)))
    checkTrue (nrow((peaks(xset1)@.Data))  > nrow((peaks(xset3)@.Data)))

}

test.xcmsSetParallel <- function() {
    file <- system.file('microtofq/MM14.mzdata', package = "msdata")

    xset1 <-  xcmsSet(files=file, method="centWave", peakwidth=c(5,12), scanrange=c(1,80))
    xset2 <-  xcmsSet(files=file, method="centWave", peakwidth=c(5,12), scanrange=c(1,80), nSlaves=2)

    checkTrue (nrow((peaks(xset1)@.Data)) == nrow((peaks(xset2)@.Data)))

}
