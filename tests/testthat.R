library(testthat)
library(xcms)
library(faahKO)
library(MSnbase)
library(msdata)
library(BiocParallel)
prm <- SerialParam()

register(SerialParam())

## Create some objects we can re-use in different tests:
faahko_3_files <- c(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
                    system.file('cdf/KO/ko16.CDF', package = "faahKO"),
                    system.file('cdf/KO/ko18.CDF', package = "faahKO"))

cwp <- CentWaveParam(noise = 10000, snthresh = 40, prefilter = c(3, 10000))
faahko_od <- readMSData(faahko_3_files, mode = "onDisk")
faahko_xod <- findChromPeaks(faahko_od, param = cwp, BPPARAM = SerialParam())
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

microtofq_fs <- c(system.file("microtofq/MM14.mzML", package = "msdata"),
                  system.file("microtofq/MM8.mzML", package = "msdata"))
microtofq_od <- readMSData(microtofq_fs, mode = "onDisk")

## Direct injection data:
fticrf <- list.files(system.file("fticr-mzML", package = "msdata"),
                     recursive = TRUE, full.names = TRUE)
fticr <- readMSData(fticrf[1:2], msLevel. = 1, mode = "onDisk")
fticr_xod <- findChromPeaks(fticr, MSWParam(scales = c(1, 7),
                                            peakThr = 80000, ampTh = 0.005,
                                            SNR.method = "data.mean",
                                            winSize.noise = 500))

## Pesticide data
fl <- system.file("TripleTOF-SWATH", "PestMix1_SWATH.mzML", package = "msdata")
pest_swth <- readMSData(fl, mode = "onDisk")
cwp2 <- CentWaveParam(snthresh = 5, noise = 100, ppm = 10,
                      peakwidth = c(3, 20), prefilter = c(3, 1000))
pest_swth <- findChromPeaks(pest_swth, param = cwp2)
pest_swth <- findChromPeaksIsolationWindow(pest_swth, param = cwp2)

fl <- system.file("TripleTOF-SWATH", "PestMix1_DDA.mzML", package = "msdata")
pest_dda <- readMSData(fl, mode = "onDisk")
pest_dda <- findChromPeaks(pest_dda, param = cwp2)

## Sciex test data.
## fl <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)
## sciex_data <- readMSData(fl, mode = "onDisk")
## sciex_data <- pickPeaks(sciex_data)

library(MsExperiment)
fls <- normalizePath(faahko_3_files)
df <- data.frame(mzML_file = basename(fls),
                 dataOrigin = fls,
                 sample = c("ko15", "ko16", "ko18"))
mse <- readMsExperiment(spectraFiles = fls, sampleData = df)
xmse <- findChromPeaks(mse, param = cwp)
expect_true(length(processHistory(xmse)) == 1L)
pdp <- PeakDensityParam(sampleGroups = rep(1, 3))
xmseg <- groupChromPeaks(xmse, param = pdp, add = FALSE)
expect_true(length(processHistory(xmseg)) == 2L)

## Data for LamaParama checks
ref <- loadXcmsData("xmse")
f <- sampleData(ref)$sample_type
f[f == "QC"] <- NA
ref <- filterFeatures(ref, PercentMissingFilter(threshold = 0, f = f))
ref_mz_rt <- featureDefinitions(ref)[, c("mzmed","rtmed")]
tst <- loadXcmsData("faahko_sub2")

test_check("xcms")
