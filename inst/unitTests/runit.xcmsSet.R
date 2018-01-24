test.findPeaksCentWave <- function() {
    file <- system.file('microtofq/MM14.mzdata', package = "msdata")

    xset1 <-  xcmsSet(files=file, method="centWave", peakwidth=c(5,12),
                      profparam = list(step = 0))
    xset2 <-  xcmsSet(files=file, method="centWave", peakwidth=c(5,12),
                      profparam = list(step = 0), scanrange=c(1,112))
    xset3 <-  xcmsSet(files=file, method="centWave", peakwidth=c(5,12),
                      profparam = list(step = 0), scanrange=c(1,80))

    checkTrue (nrow((peaks(xset1)@.Data)) == nrow((peaks(xset2)@.Data)))
    checkTrue (nrow((peaks(xset1)@.Data))  > nrow((peaks(xset3)@.Data)))

}

test.findPeaksMatchedFilter <- function() {
    file <- system.file('microtofq/MM14.mzdata', package = "msdata")

    xset1 <- xcmsSet(files=file, method="matchedFilter", fwhm=10)
    xset2 <- xcmsSet(files=file, method="matchedFilter", fwhm=10,
                     scanrange=c(1,112))
    xset3 <- xcmsSet(files=file, method="matchedFilter", fwhm=10,
                     scanrange=c(1,80))

    checkTrue (nrow((peaks(xset1)@.Data)) == nrow((peaks(xset2)@.Data)))
    checkTrue (nrow((peaks(xset1)@.Data))  > nrow((peaks(xset3)@.Data)))

}

test.xcmsSetParallel <- function() {
    file <- system.file('microtofq/MM14.mzdata', package = "msdata")

    xset1 <-  xcmsSet(files=file, method="centWave", peakwidth=c(5,12),
                      scanrange=c(1,80), profparam = list(step = 0))
    xset2 <-  xcmsSet(files=file, method="centWave", peakwidth=c(5,12),
                      scanrange=c(1,80), profparam = list(step = 0))
    ## parallel disabled: , nSlaves=2)

    checkTrue (nrow((peaks(xset1)@.Data)) == nrow((peaks(xset2)@.Data)))

}

test.phenoDataFromPaths <- function() {
    base_dir <- system.file("cdf", package = "faahKO")
    cdf_files <- list.files(base_dir, recursive = TRUE, full.names = TRUE)
    pd <- phenoDataFromPaths(cdf_files)
    checkTrue(colnames(pd) == "class")
    checkEquals(levels(pd$class), c("KO", "WT"))
}
