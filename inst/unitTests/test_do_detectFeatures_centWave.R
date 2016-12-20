## Test detectFeatures centWave

library(faahKO)
fs <- c(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
        system.file('cdf/KO/ko16.CDF', package = "faahKO"),
        system.file('cdf/KO/ko18.CDF', package = "faahKO"),
        system.file('cdf/KO/ko19.CDF', package = "faahKO"))

library(msdata)
mzf <- c(system.file("microtofq/MM14.mzML", package = "msdata"),
         system.file("microtofq/MM8.mzML", package = "msdata"))
f <- msdata::proteomics(full.names = TRUE, pattern = "TMT_Erwinia")


test_do_detectFeatures_centWave <- function() {
    xr <- xcmsRaw(fs[1], profstep = 0)
    ## We expect that changing a parameter has an influence on the result.
    mzVals <- xr@env$mz
    intVals <- xr@env$intensity
    ## Define the values per spectrum:
    valsPerSpect <- diff(c(xr@scanindex, length(mzVals)))
    res1 <- do_detectFeatures_centWave(mz = mzVals,
                                       int = intVals,
                                       scantime = xr@scantime,
                                       valsPerSpect,
                                       snthresh = 200,
                                       noise = 4000)
    res2 <- do_detectFeatures_centWave(mz = mzVals,
                                       int = intVals,
                                       scantime = xr@scantime,
                                       valsPerSpect,
                                       snthresh = 500,
                                       noise = 4000)
    checkTrue(nrow(res1) > nrow(res2))

    ## Check scanrange on findPeaks.centWave.
    res_1 <- findPeaks.centWave(xr, scanrange = c(90, 345), noise = 2000)
    xr <- xr[90:345]
    mzVals <- xr@env$mz
    intVals <- xr@env$intensity
    ## Define the values per spectrum:
    valsPerSpect <- diff(c(xr@scanindex, length(mzVals)))
    res_2 <- do_detectFeatures_centWave(mz = mzVals, int = intVals,
                                        scantime = xr@scantime, valsPerSpect,
                                        noise = 2000)
    checkEquals(res_1@.Data, res_2)
}

## Evaluate the featureDetection method using the centWave method on
## OnDiskMSnExp and on MSnExp objects.
test_featureDetection_centWave <- function() {
    ## Control
    library(MSnbase)
    xr <- xcmsRaw(f[1], profstep = 0)
    ppm <- 40
    snthresh <- 40
    suppressWarnings(
        res_x <- findPeaks.centWave(xr, ppm = ppm, snthresh = snthresh,
                                    noise = 100000)@.Data
    )
    ## Bypass xcmsRaw
    suppressWarnings(
        xs <- xcmsSet(f[1], profparam = list(profstep = 0), ppm = ppm,
                      snthresh = snthresh, method = "centWave",
                      noise = 100000)
    )
    checkEquals(xs@peaks[, colnames(res_x)], res_x)
    ## OnDiskMSnExp
    onDisk <- readMSData2(f[1], msLevel. = 1)
    cwp <- CentWaveParam(ppm = ppm, snthresh = snthresh, noise = 100000)
    suppressWarnings(
        res <- detectFeatures(onDisk, param = cwp, return.type = "list")
    )
    checkEquals(res[[1]], peaks(xs)@.Data)

    ## MSnExp
    inMem <- readMSData(f[1], msLevel. = 1)
    suppressWarnings(
        res_2 <- detectFeatures(inMem, param = cwp, return.type = "list")
    )
    checkEquals(res_2[[1]], peaks(xs)@.Data)

    ## returning an xcmsSet
    suppressWarnings(
        res <- detectFeatures(onDisk, param = cwp, return.type = "xcmsSet")
    )
    checkEquals(peaks(res), peaks(xs))
    suppressWarnings(
        res <- detectFeatures(inMem, param = cwp, return.type = "xcmsSet")
    )
    checkEquals(peaks(res), peaks(xs))
}

dontrun_test_benchmark_centWaves <- function() {
    library(msdata)
    f <- msdata::proteomics(full.names = TRUE, pattern = "TMT_Erwinia")
    library(microbenchmark)
    library(MSnbase)
    library(xcms)
    ##
    xr <- xcmsRaw(f[1], profstep = 0)
    ppm <- 40
    snthresh <- 40

    cwp <- CentWaveParam(ppm = ppm, snthresh = snthresh)
    onDisk <- readMSData2(f[1], msLevel. = 1)
    register(SerialParam())
    system.time(
        tmp <- detectFeatures(onDisk, param = cwp)
    ) ## 9.7sec
    system.time(
        tmp <- detectFeatures(onDisk, param = cwp, return.type = "xcmsSet")
    ) ## 12sec
    system.time(
        tmp <- xcmsSet(f[1], profparam = list(profstep = 0), ppm = ppm,
                       snthresh = snthresh, method = "centWave")
    ) ## 11.99sec

    inMem <- readMSData(f[1], msLevel. = 1)
    register(SerialParam())

    ## detectFeatures,MSnExp and findPeaks.centWave should be about similar.
    microbenchmark(findPeaks.centWave(xr, ppm = ppm, snthresh = snthresh),
                   detectFeatures(inMem, param = cwp), times = 3)
    ## findPeaks.centWave is about 1 second faster.

    ## detectFeatures,OnDiskMSnExp and xcmsSet should be about similar.
    microbenchmark(xcmsSet(f[1], profparam = list(profstep = 0), ppm = ppm,
                           snthresh = snthresh, method = "centWave"),
                   detectFeatures(onDisk, param = cwp),
                   detectFeatures(inMem, param = cwp),
                   times = 3)
}



############################################################
## This is only relevant during development of the do_ function
## to evaluate that results are identical.
dontrun_test_do_detectFeatures_centWave_impl <- function() {

    for (i in 1:length(fs)) {
        ppm = 25
        peakwidth = c(20, 50)
        snthresh = 10
        prefilter = c(3, 100)
        mzCenterFun = "wMean"
        integrate = 1
        mzdiff = -0.001
        fitgauss = FALSE
        noise = 0
        verboseColumns = FALSE

        xr <- xcmsRaw(fs[i])

        ## Default settings
        .runAndCompare(xr, ppm, peakwidth, snthresh, prefilter, mzCenterFun,
                       integrate, mzdiff, fitgauss, noise, verboseColumns)
        ## xcms: 14.6 sec
        ## do_ : 13 sec

        ppm <- 10
        .runAndCompare(xr, ppm, peakwidth, snthresh, prefilter, mzCenterFun,
                   integrate, mzdiff, fitgauss, noise, verboseColumns)
        ## xcms: 15 sec
        ## do_ : 13.3 sec

        peakwidth <- c(3, 30)
        .runAndCompare(xr, ppm, peakwidth, snthresh, prefilter, mzCenterFun,
                       integrate, mzdiff, fitgauss, noise, verboseColumns)
        ## xcms: 11.4 sec
        ## do_ :  9.5 sec

        snthresh <- 15
        .runAndCompare(xr, ppm, peakwidth, snthresh, prefilter, mzCenterFun,
                       integrate, mzdiff, fitgauss, noise, verboseColumns)
        ## xcms: 10.6 sec
        ## do_ :  8.8 sec

        fitgauss <- TRUE
        .runAndCompare(xr, ppm, peakwidth, snthresh, prefilter, mzCenterFun,
                       integrate, mzdiff, fitgauss, noise, verboseColumns)
        ## xcms: 12.5 sec
        ## do_ : 10.7 sec

        verboseColumns <- TRUE
        .runAndCompare(xr, ppm, peakwidth, snthresh, prefilter, mzCenterFun,
                       integrate, mzdiff, fitgauss, noise, verboseColumns)
        ## xcms: 12.2 sec
        ## do_ : 10.6 sec
    }
}

## That's to compare the functions in version 1.49.7.
.runAndCompare <- function(xr, ppm, peakwidth, snthresh, prefilter, mzCenterFun,
                           integrate, mzdiff, fitgauss, noise, verboseColumns) {
    require(RUnit)
    mz <- xr@env$mz
    int <- xr@env$intensity
    scantime <- xr@scantime
    scanindex <- xr@scanindex
    a <- system.time(
        ## That's the method called inside do_...
        xrDo <- xcms:::.centWave_orig(mz = mz, int = int, scantime = scantime,
                                      valsPerSpect = diff(c(scanindex, length(mz))),
                                      ppm = ppm, peakwidth = peakwidth,
                                      snthresh = snthresh,
                                      prefilter = prefilter,
                                      mzCenterFun = mzCenterFun,
                                      integrate = integrate,
                                      mzdiff = mzdiff,
                                      fitgauss = fitgauss,
                                      noise = noise,
                                      verboseColumns = verboseColumns)
    ) ## 12.7
    ## Run the original centWave code on xcmsRaw:
    b <- system.time(
        xrPeaks <- xcms:::.findPeaks.centWave_orig(xr,
                                                   ppm = ppm,
                                                   peakwidth = peakwidth,
                                                   snthresh = snthresh,
                                                   prefilter = prefilter,
                                                   mzCenterFun = mzCenterFun,
                                                   integrate = integrate,
                                                   mzdiff = mzdiff,
                                                   fitgauss = fitgauss,
                                                   noise = noise,
                                                   verbose.columns = verboseColumns)
    )  ## 15.4
    ## Compare.
    cat("DO: ", a, "\n")
    cat("XCMS: ", b, "\n")
    if (!checkEquals(new("xcmsPeaks", xrDo), xrPeaks))
        stop("do_ and xcms yield different results!")
}


## Some speed tests.
.otherTest <- function() {
    Testv <- c(2, 4.2, 34.1, 34.5, 6.4, 6.3, 1.2)
    RforM <- matrix(nrow = 0, ncol = length(Testv))
    system.time(
        for(i in 1:5000){
            RforM <- rbind(RforM, Testv)
        }
    ) ## 1.27
    ## with append to list.
    RforL <- vector("list", 0)
    system.time(
        for(i in 1:5000){
            RforL <- c(RforL, Testv)
        }
    ) ## 1.12
    system.time(
        RapplyL <- lapply(1:5000, function(z) {return(Testv)})
    ) ## 0.003
    RM <- matrix(nrow=5000, ncol = length(Testv))
    system.time(
        for (i in 1:5000) {
            RM[i, ] <- Testv
        }
    ) ## 0.006

    ## Compare adding to list instead of adding to existing. [[]]
    RexL <- vector("list", 5000)
    system.time(
        for (i in 1:5000){
            RexL[[i]] <- Testv
        }
    ) ## 0.005
    ## Dynamically...
    RexL <- list()
    system.time(
        for (i in 1:5000){
            RexL[[i]] <- Testv
        }
    ) ## 0.005

}
