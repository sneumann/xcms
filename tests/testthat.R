library(testthat)
library(xcms)
<<<<<<< Updated upstream
library(faahKO)
library(msdata)
=======
library(MSnbase)
library(MsDataHub)
library(BiocParallel)
prm <- SerialParam()
>>>>>>> Stashed changes

if (.Platform$OS.type == "unix") {
    prm <- MulticoreParam(3)
} else {
    # prm <- SnowParam(3)
    prm <- SerialParam()
}
register(bpstart(prm))

## Create some objects we can re-use in different tests:
faahko_3_files <- c(ko15.CDF(),
                     ko16.CDF(),
                     ko18.CDF())

faahko_od <- readMSData(faahko_3_files, mode = "onDisk")
faahko_xod <- findChromPeaks(
    faahko_od, param = CentWaveParam(noise = 10000, snthresh = 40,
                                     prefilter = c(3, 10000)))
od_x <- faahko_od
mzr <- matrix(c(335, 335, 344, 344), ncol = 2, byrow = TRUE)
od_chrs <- chromatogram(od_x, mz = mzr)
xod_x <- faahko_xod
pdp <- PeakDensityParam(sampleGroups = rep(1, 3))
xod_xg <- groupChromPeaks(xod_x, param = pdp)
xod_xgr <- adjustRtime(xod_xg, param = PeakGroupsParam(span = 0.4))
xod_xgrg <- groupChromPeaks(xod_xgr, param = pdp)
xod_r <- adjustRtime(as(od_x, "XCMSnExp"), param = ObiwarpParam())

xod_chr <- findChromPeaks(filterMz(filterRt(od_x, rt = c(2500, 3500)),
                                   mz = c(334.9, 344.1)),
                          param = CentWaveParam())

microtofq_fs <- c(MM14.mzML(),
                  MM8.mzML())

microtofq_od <- readMSData(microtofq_fs, mode = "onDisk")

## Direct injection data:
fticrf <- c(HAM004_641fE_14.11.07..Exp1.extracted.mzML(),
            HAM004_641fE_14.11.07..Exp2.extracted.mzML())

fticr <- readMSData(fticrf[1:2], msLevel. = 1, mode = "onDisk")
fticr_xod <- findChromPeaks(fticr, MSWParam(scales = c(1, 7),
                                            peakThr = 80000, ampTh = 0.005,
                                            SNR.method = "data.mean",
                                            winSize.noise = 500))
## Pesticide data
fl <- PestMix1_SWATH.mzML()

pest_swth <- readMSData(fl, mode = "onDisk")
cwp <- CentWaveParam(snthresh = 5, noise = 100, ppm = 10,
                     peakwidth = c(3, 20), prefilter = c(3, 1000))
pest_swth <- findChromPeaks(pest_swth, param = cwp)
pest_swth <- findChromPeaksIsolationWindow(pest_swth, param = cwp)

fl <- PestMix1_DDA.mzML()
pest_dda <- readMSData(fl, mode = "onDisk")
pest_dda <- findChromPeaks(pest_dda, param = cwp)

## Sciex test data.
## fl <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)
## sciex_data <- readMSData(fl, mode = "onDisk")
## sciex_data <- pickPeaks(sciex_data)

test_check("xcms")

bpstop(prm)

