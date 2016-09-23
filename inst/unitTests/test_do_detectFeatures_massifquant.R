############################################################
## do_detectFeatures_massifquant tests

library(msdata)
mzf <- c(system.file("microtofq/MM14.mzML", package = "msdata"),
         system.file("microtofq/MM8.mzML", package = "msdata"))

xraw <- xcmsRaw(mzf[1], profstep = 0)


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
        a <- findPeaks.massifquant(xr, ppm = ppm, peakwidth = peakwidth,
                                   snthresh = snthresh,
                                   criticalValue = criticalValue,
                                   consecMissedLimit = consecMissedLimit,
                                   unions = unions, withWave = withWave)
        ## b <- xcms:::.massifquant_orig(mz, int, scantime = scantime,
        ##                               valsPerSpect = valsPerSpect,
        ##                               ppm = ppm, peakwidth = peakwidth,
        ##                               snthresh = snthresh,
        ##                               criticalValue = criticalValue,
        ##                               consecMissedLimit = consecMissedLimit,
        ##                               unions = unions, withWave = withWave)
        ## checkEquals(a@.Data, b)
        d <- xcms:::.massifquant(mz, int, scantime = scantime,
                                 valsPerSpect = valsPerSpect,
                                 ppm = ppm, peakwidth = peakwidth,
                                 snthresh = snthresh,
                                 criticalValue = criticalValue,
                                 consecMissedLimit = consecMissedLimit,
                                 unions = unions, withWave = withWave)
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
        ## b <- xcms:::.massifquant_orig(mz, int, scantime, valsPerSpect,
        ##                               ppm = ppm, peakwidth = peakwidth,
        ##                               snthresh = snthresh,
        ##                               criticalValue = criticalValue,
        ##                               consecMissedLimit = consecMissedLimit,
        ##                               unions = unions, withWave = withWave)
        ## checkEquals(a@.Data, b)
        d <- xcms:::.massifquant(mz, int, scantime = scantime,
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
        criticalValue <- 1.125
        consecMissedLimit <- 2
        unions <- 1
        withWave <- 1
        a <- findPeaks.massifquant(xr, ppm = ppm, peakwidth = peakwidth,
                                   snthresh = snthresh,
                                   criticalValue = criticalValue,
                                   consecMissedLimit = consecMissedLimit,
                                   unions = unions, withWave = withWave)
        ## b <- xcms:::.massifquant_orig(mz, int, scantime, valsPerSpect,
        ##                               ppm = ppm, peakwidth = peakwidth,
        ##                               snthresh = snthresh,
        ##                               criticalValue = criticalValue,
        ##                               consecMissedLimit = consecMissedLimit,
        ##                               unions = unions, withWave = withWave)
        ## checkEquals(a@.Data, b)
        d <- xcms:::.massifquant(mz, int, scantime = scantime,
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
        ## b <- xcms:::.massifquant_orig(mz, int, scantime, valsPerSpect,
        ##                               ppm = ppm, peakwidth = peakwidth,
        ##                               snthresh = snthresh,
        ##                               criticalValue = criticalValue,
        ##                               consecMissedLimit = consecMissedLimit,
        ##                               unions = unions, withWave = withWave)
        ## checkEquals(a@.Data, b)
        d <- xcms:::.massifquant(mz, int, scantime = scantime,
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
        ## b <- xcms:::.massifquant_orig(mz, int, scantime, valsPerSpect,
        ##                               ppm = ppm, peakwidth = peakwidth,
        ##                               snthresh = snthresh,
        ##                               criticalValue = criticalValue,
        ##                               consecMissedLimit = consecMissedLimit,
        ##                               unions = unions, withWave = withWave)
        ## checkEquals(a@.Data, b)
        d <- xcms:::.massifquant(mz, int, scantime = scantime,
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
        ## b <- xcms:::.massifquant_orig(mz, int, scantime, valsPerSpect,
        ##                               ppm = ppm, peakwidth = peakwidth,
        ##                               snthresh = snthresh,
        ##                               criticalValue = criticalValue,
        ##                               consecMissedLimit = consecMissedLimit,
        ##                               unions = unions, withWave = withWave)
        ## checkEquals(a@.Data, b)
        d <- xcms:::.massifquant(mz, int, scantime = scantime,
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
        ## b <- xcms:::.massifquant_orig(mz, int, scantime, valsPerSpect,
        ##                               ppm = ppm, peakwidth = peakwidth,
        ##                               snthresh = snthresh,
        ##                               criticalValue = criticalValue,
        ##                               consecMissedLimit = consecMissedLimit,
        ##                               unions = unions, withWave = withWave)
        ## checkEquals(a@.Data, b)
        d <- xcms:::.massifquant(mz, int, scantime = scantime,
                                 valsPerSpect = valsPerSpect,
                                 ppm = ppm, peakwidth = peakwidth,
                                 snthresh = snthresh,
                                 criticalValue = criticalValue,
                                 consecMissedLimit = consecMissedLimit,
                                 unions = unions, withWave = withWave)
        checkEquals(a@.Data, d)
    }
}

############################################################
## Benchmarking.
dontrun_benchmar_massifquant <- function() {
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
