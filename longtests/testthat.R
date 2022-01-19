library(testthat)
library(xcms)
library(faahKO)
library(msdata)

attr(faahko, "filepaths") <- sapply(
    as.list(basename(attr(faahko, "filepaths"))),
    function(x) system.file("cdf", if (length(grep("ko",x)) > 0) "KO" else "WT",
                            x, package = "faahKO"))
if (.Platform$OS.type == "unix") {
    prm <- MulticoreParam(3)
} else {
    prm <- SnowParam(3)
}
register(bpstart(prm))
## register(SerialParam())

## Create some objects we can re-use in different tests:
faahko_3_files <- c(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
                    system.file('cdf/KO/ko16.CDF', package = "faahKO"),
                    system.file('cdf/KO/ko18.CDF', package = "faahKO"))

## An xcmsRaw for the first file:
faahko_xr_1 <- xcmsRaw(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
                       profstep = 0)
faahko_od <- readMSData(faahko_3_files, mode = "onDisk")
faahko_xod <- findChromPeaks(
    faahko_od, param = CentWaveParam(noise = 10000, snthresh = 40,
                                     prefilter = c(3, 10000)))
faahko_xs <- xcmsSet(faahko_3_files, profparam = list(step = 0),
                     method = "centWave", noise = 10000, snthresh = 40,
                     prefilter = c(3, 10000))
faahko_xsg <- group(faahko_xs)
## Doing also the retention time correction etc
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

faahko_grouped_filled <- fillPeaks(group(faahko))
faahko_grouped_retcor_filled <-
    fillPeaks(group(retcor(group(updateObject(faahko)))))

microtofq_fs <- c(system.file("microtofq/MM14.mzML", package = "msdata"),
                  system.file("microtofq/MM8.mzML", package = "msdata"))
microtofq_xr <- xcmsRaw(microtofq_fs[1], profstep = 0)
microtofq_od <- readMSData(microtofq_fs, mode = "onDisk")

## Direct injection data:
fticrf <- list.files(system.file("fticr", package = "msdata"),
                     recursive = TRUE, full.names = TRUE)
fticr <- readMSData(fticrf[1:2], msLevel. = 1, mode = "onDisk")
fticr_xod <- findChromPeaks(fticr, MSWParam(scales = c(1, 7),
                                            peakThr = 80000, ampTh = 0.005,
                                            SNR.method = "data.mean",
                                            winSize.noise = 500))
fticr_xs <- xcmsSet(method="MSW", files=fticrf[1:2], scales=c(1,7),
                    SNR.method='data.mean' , winSize.noise=500,
                    peakThr=80000,  amp.Th=0.005)

fs <- c(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
        system.file('cdf/KO/ko16.CDF', package = "faahKO"),
        system.file('cdf/KO/ko18.CDF', package = "faahKO"),
        system.file('cdf/KO/ko19.CDF', package = "faahKO"))
xs_1 <- xcmsSet(fs, profparam = list(step = 0), method = "centWave",
                noise = 10000, snthresh = 50, prefilter = c(3, 10000))

test_check("xcms")
