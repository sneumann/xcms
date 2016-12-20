############################################################
## do_detectFeatures_massifquant tests

library(faahKO)
fs <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
xraw <- xcmsRaw(fs, profstep = 0)

library(msdata)
f <- msdata::proteomics(full.names = TRUE, pattern = "TMT_Erwinia")
mzf <- c(system.file("microtofq/MM14.mzML", package = "msdata"),
         system.file("microtofq/MM8.mzML", package = "msdata"))

############################################################
## Simple test comparing results from various massifquant runs and
## centWave analyses.
test_do_detectFeatures_massifquant <- function() {
    res <- findPeaks.massifquant(xraw, snthresh = 100)
    mz <- xraw@env$mz
    int <- xraw@env$intensity
    valsPerSpect <- diff(c(xraw@scanindex, length(mz)))
    scantime <- xraw@scantime
    res_2 <- do_detectFeatures_massifquant(mz = mz, int = int,
                                           valsPerSpect = valsPerSpect,
                                           scantime = scantime)
    checkEquals(res@.Data, res_2)
    ## With centWave:
    res_3 <- do_detectFeatures_massifquant(mz = mz, int = int,
                                           valsPerSpect = valsPerSpect,
                                           scantime = scantime, withWave = TRUE,
                                           snthresh = 100, noise = 4000)
    res_4 <- findPeaks.massifquant(xraw, withWave = 1, snthresh = 100,
                                   noise = 4000)
    checkEquals(res_3, res_4@.Data)
    checkTrue(nrow(res_3) < nrow(res_2))

    ## Subsetted data and scanrange:
    res_1 <- findPeaks.massifquant(xraw, scanrange = c(90, 345))
    xsub <- xraw[90:345]
    mz <- xsub@env$mz
    int <- xsub@env$intensity
    valsPerSpect <- diff(c(xsub@scanindex, length(mz)))
    scantime <- xsub@scantime
    res_2 <- do_detectFeatures_massifquant(mz = mz, int = int,
                                           valsPerSpect = valsPerSpect,
                                           scantime = scantime)
    checkIdentical(res_1@.Data, res_2)
}

## Evaluate the featureDetection method using massifquant on MSnExp and
## OnDiskMSnExp objects.
test_featureDetection_massifquant <- function() {
    library(MSnbase)
    mqp <- MassifquantParam(ppm = 20, criticalValue = 1.2)
    res <- xcmsSet(mzf[1], method = "massifquant", ppm = 20, criticalValue = 1.2)
    ## onDisk
    onDisk <- readMSData2(mzf[1])
    res_o <- detectFeatures(onDisk, param = mqp, return.type = "xcmsSet")
    checkEquals(peaks(res_o), peaks(res))
    checkEquals(res_o@rt$raw, res@rt$raw, checkNames = FALSE)
    ## inMem
    inMem <- readMSData(mzf[1], msLevel. = 1)
    res_i <- detectFeatures(inMem, param = mqp, return.type = "xcmsSet")
    checkEquals(peaks(res_i), peaks(res))
    checkEquals(res_i@rt$raw, res@rt$raw, checkNames = FALSE)
}


############################################################
## Test the implementation of the "do" function, i.e. whether
## the results are the same between versions and implementations.
dontrun_test_do_detectFeatures_massifquant_impl <- function() {

    library(xcms)
    library(RUnit)

    library(faahKO)
    fs <- c(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
            system.file('cdf/KO/ko16.CDF', package = "faahKO"),
            system.file('cdf/KO/ko18.CDF', package = "faahKO"),
            system.file('cdf/KO/ko19.CDF', package = "faahKO"))
    i <- 1
    for (i in 1:length(fs)) {
        cat("============================================================\n")
        cat("|   file", i, ":", fs[i], "\n")
        cat("------------------------------------------------------------")

        xr <- xcmsRaw(fs[i])
        mz <- xr@env$mz
        int <- xr@env$intensity
        valsPerSpect <- diff(c(xr@scanindex, length(mz)))
        scantime <- xr@scantime

        ## Default
        ppm <- 10
        peakwidth <- c(20, 50)
        snthresh <- 10
        criticalValue <- 1.125
        consecMissedLimit <- 2
        unions <- 1
        withWave <- 0
        ## Now, this method calls do_detectFeatures_massifquant.
        a <- findPeaks.massifquant(xr, ppm = ppm, peakwidth = peakwidth,
                                   snthresh = snthresh,
                                   criticalValue = criticalValue,
                                   consecMissedLimit = consecMissedLimit,
                                   unions = unions, withWave = withWave)
        b <- xcms:::findPeaks.massifquant_orig(xr, ppm = ppm, peakwidth = peakwidth,
                                               snthresh = snthresh,
                                               criticalValue = criticalValue,
                                               consecMissedLimit = consecMissedLimit,
                                               unions = unions, withWave = withWave)
        checkEquals(a@.Data, b)
        ## LLL: compare the _orig method
        ## 1) check if scanrange works between both.
        ## 2) compare the orig_ with the do
        d <- do_detectFeatures_massifquant(mz, int, scantime = scantime,
                                           valsPerSpect = valsPerSpect,
                                           ppm = ppm, peakwidth = peakwidth,
                                           snthresh = snthresh,
                                           criticalValue = criticalValue,
                                           consecMissedLimit = consecMissedLimit,
                                           unions = unions,
                                           withWave = as.logical(withWave))
        checkEquals(a@.Data, d)


        ## Modify settings...
        ## Slightly modified one:
        ppm <- 20
        peakwidth <- c(20, 50)
        snthresh <- 10
        criticalValue <- 1.125
        consecMissedLimit <- 2
        unions <- 1
        withWave <- 0
        a <- findPeaks.massifquant(xr, ppm = ppm, peakwidth = peakwidth,
                                   snthresh = snthresh,
                                   criticalValue = criticalValue,
                                   consecMissedLimit = consecMissedLimit,
                                   unions = unions, withWave = withWave)
        b <- xcms:::findPeaks.massifquant_orig(xr, ppm = ppm, peakwidth = peakwidth,
                                               snthresh = snthresh,
                                               criticalValue = criticalValue,
                                               consecMissedLimit = consecMissedLimit,
                                               unions = unions, withWave = withWave)
        checkEquals(a@.Data, b)
        d <- do_detectFeatures_massifquant(mz, int, scantime = scantime,
                                 valsPerSpect = valsPerSpect,
                                 ppm = ppm, peakwidth = peakwidth,
                                 snthresh = snthresh,
                                 criticalValue = criticalValue,
                                 consecMissedLimit = consecMissedLimit,
                                 unions = unions,
                                 withWave = as.logical(withWave))
        checkEquals(a@.Data, d)

        ppm <- 20
        peakwidth <- c(20, 50)
        snthresh <- 10
        criticalValue <- 1.125
        consecMissedLimit <- 2
        unions <- 1
        withWave <- 1
        a <- findPeaks.massifquant(xr, ppm = ppm, peakwidth = peakwidth,
                                   snthresh = snthresh,
                                   criticalValue = criticalValue,
                                   consecMissedLimit = consecMissedLimit,
                                   unions = unions, withWave = withWave)
        b <- xcms:::findPeaks.massifquant_orig(xr, ppm = ppm, peakwidth = peakwidth,
                                               snthresh = snthresh,
                                               criticalValue = criticalValue,
                                               consecMissedLimit = consecMissedLimit,
                                               unions = unions, withWave = withWave)
        checkEquals(a, b)
        d <- do_detectFeatures_massifquant(mz, int, scantime = scantime,
                                           valsPerSpect = valsPerSpect,
                                           ppm = ppm, peakwidth = peakwidth,
                                           snthresh = snthresh,
                                           criticalValue = criticalValue,
                                           consecMissedLimit = consecMissedLimit,
                                           unions = unions, withWave = withWave)
        checkEquals(a@.Data, d)

        ppm <- 20
        peakwidth <- c(20, 50)
        snthresh <- 10
        criticalValue <- 1.75
        consecMissedLimit <- 2
        unions <- 1
        withWave <- 1
        a <- findPeaks.massifquant(xr, ppm = ppm, peakwidth = peakwidth,
                                   snthresh = snthresh,
                                   criticalValue = criticalValue,
                                   consecMissedLimit = consecMissedLimit,
                                   unions = unions, withWave = withWave)
        b <- xcms:::findPeaks.massifquant_orig(xr, ppm = ppm, peakwidth = peakwidth,
                                               snthresh = snthresh,
                                               criticalValue = criticalValue,
                                               consecMissedLimit = consecMissedLimit,
                                               unions = unions, withWave = withWave)
        checkEquals(a, b)
        d <- do_detectFeatures_massifquant(mz, int, scantime = scantime,
                                           valsPerSpect = valsPerSpect,
                                           ppm = ppm, peakwidth = peakwidth,
                                           snthresh = snthresh,
                                           criticalValue = criticalValue,
                                           consecMissedLimit = consecMissedLimit,
                                           unions = unions, withWave = withWave)
        checkEquals(a@.Data, d)

        ppm <- 20
        peakwidth <- c(20, 50)
        snthresh <- 299
        criticalValue <- 1.75
        consecMissedLimit <- 2
        unions <- 1
        withWave <- 0
        a <- findPeaks.massifquant(xr, ppm = ppm, peakwidth = peakwidth,
                                   snthresh = snthresh,
                                   criticalValue = criticalValue,
                                   consecMissedLimit = consecMissedLimit,
                                   unions = unions, withWave = withWave)
        b <- xcms:::findPeaks.massifquant_orig(xr, ppm = ppm, peakwidth = peakwidth,
                                               snthresh = snthresh,
                                               criticalValue = criticalValue,
                                               consecMissedLimit = consecMissedLimit,
                                               unions = unions, withWave = withWave)
        checkEquals(a@.Data, b)
        d <- do_detectFeatures_massifquant(mz, int, scantime = scantime,
                                           valsPerSpect = valsPerSpect,
                                           ppm = ppm, peakwidth = peakwidth,
                                           snthresh = snthresh,
                                           criticalValue = criticalValue,
                                           consecMissedLimit = consecMissedLimit,
                                           unions = unions, withWave = withWave)
        checkEquals(a@.Data, d)

        ppm <- 20
        peakwidth <- c(20, 50)
        snthresh <- 299
        criticalValue <- 1.75
        consecMissedLimit <- 4
        unions <- 1
        withWave <- 0
        a <- findPeaks.massifquant(xr, ppm = ppm, peakwidth = peakwidth,
                                   snthresh = snthresh,
                                   criticalValue = criticalValue,
                                   consecMissedLimit = consecMissedLimit,
                                   unions = unions, withWave = withWave)
        b <- xcms:::findPeaks.massifquant_orig(xr, ppm = ppm, peakwidth = peakwidth,
                                               snthresh = snthresh,
                                               criticalValue = criticalValue,
                                               consecMissedLimit = consecMissedLimit,
                                               unions = unions, withWave = withWave)
        checkEquals(a@.Data, b)
        d <- do_detectFeatures_massifquant(mz, int, scantime = scantime,
                                           valsPerSpect = valsPerSpect,
                                           ppm = ppm, peakwidth = peakwidth,
                                           snthresh = snthresh,
                                           criticalValue = criticalValue,
                                           consecMissedLimit = consecMissedLimit,
                                           unions = unions, withWave = withWave)
        checkEquals(a@.Data, d)

        ppm <- 20
        peakwidth <- c(20, 50)
        snthresh <- 299
        criticalValue <- 1.75
        consecMissedLimit <- 4
        unions <- 2
        withWave <- 1
        a <- findPeaks.massifquant(xr, ppm = ppm, peakwidth = peakwidth,
                                   snthresh = snthresh,
                                   criticalValue = criticalValue,
                                   consecMissedLimit = consecMissedLimit,
                                   unions = unions, withWave = withWave)
        b <- xcms:::findPeaks.massifquant_orig(xr, ppm = ppm, peakwidth = peakwidth,
                                               snthresh = snthresh,
                                               criticalValue = criticalValue,
                                               consecMissedLimit = consecMissedLimit,
                                               unions = unions, withWave = withWave)
        checkEquals(a, b)
        d <- do_detectFeatures_massifquant(mz, int, scantime = scantime,
                                           valsPerSpect = valsPerSpect,
                                           ppm = ppm, peakwidth = peakwidth,
                                           snthresh = snthresh,
                                           criticalValue = criticalValue,
                                           consecMissedLimit = consecMissedLimit,
                                           unions = unions, withWave = withWave)
        checkEquals(a@.Data, d)

        ## Specifying scanrange:
        scnr <- c(300, 800)
        a <- findPeaks.massifquant(xr, ppm = ppm, peakwidth = peakwidth,
                                   snthresh = snthresh,
                                   criticalValue = criticalValue,
                                   consecMissedLimit = consecMissedLimit,
                                   unions = unions, withWave = withWave,
                                   scanrange = scnr)
        ## That fails because scanrange is not considered in the Kalman filter
        ## based peak detection. This has been reported as issue #61
        checkException(
        b <- xcms:::findPeaks.massifquant_orig(xr, ppm = ppm, peakwidth = peakwidth,
                                               snthresh = snthresh,
                                               criticalValue = criticalValue,
                                               consecMissedLimit = consecMissedLimit,
                                               unions = unions, withWave = withWave,
                                               scanrange = scnr)
        )
        ## without withWave
        withWave <- 0
        scnr <- c(400, 410)
        ## ??? scanrange ignored???
        a <- findPeaks.massifquant(xr, ppm = ppm, peakwidth = peakwidth,
                                   snthresh = snthresh,
                                   criticalValue = criticalValue,
                                   consecMissedLimit = consecMissedLimit,
                                   unions = unions, withWave = withWave,
                                   scanrange = scnr)
        b <- xcms:::findPeaks.massifquant_orig(xr, ppm = ppm, peakwidth = peakwidth,
                                               snthresh = snthresh,
                                               criticalValue = criticalValue,
                                               consecMissedLimit = consecMissedLimit,
                                               unions = unions, withWave = withWave,
                                               scanrange = scnr)
        checkEquals(a@.Data, b)
    }
}

############################################################
## Benchmarking.
dontrun_benchmark_massifquant <- function() {
    library(microbenchmark)
    mz <- xraw@env$mz
    int <- xraw@env$intensity
    valsPerSpect <- diff(c(xraw@scanindex, length(mz)))
    scantime <- xraw@scantime

    microbenchmark(findPeaks.massifquant(xraw),
                   xcms:::.massifquant_orig(mz, int, scantime = scantime,
                                            valsPerSpect))
    microbenchmark(xcms:::.massifquant(mz, int, scantime = scantime,
                                       valsPerSpect),
                   xcms:::.massifquant_orig(mz, int, scantime = scantime,
                                            valsPerSpect))
    ## No difference;
}
