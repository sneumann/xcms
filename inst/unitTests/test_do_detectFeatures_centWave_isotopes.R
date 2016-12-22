## Tests related to the findPeaks.centWaveWithAddIsotopeROIs.
library(faahKO)
fs <- c(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
        system.file('cdf/KO/ko16.CDF', package = "faahKO"),
        system.file('cdf/KO/ko18.CDF', package = "faahKO"),
        system.file('cdf/KO/ko19.CDF', package = "faahKO"))
xr <- xcmsRaw(fs[1], profstep = 0)
mzVals <- xr@env$mz
intVals <- xr@env$intensity
## f <- msdata::proteomics(full.names = TRUE, pattern = "TMT_Erwinia")

test_do_detectFeatures_centWaveWithPredIsoROIs <- function() {
    ## initial centWave:
    valsPerSpect <- diff(c(xr@scanindex, length(mzVals)))
    feats_1 <- do_detectFeatures_centWave(mz = mzVals, int = intVals,
                                          scantime = xr@scantime,
                                          valsPerSpect = valsPerSpect,
                                          noise = 1500, verboseColumns = TRUE)
    feats_2 <- do_detectFeatures_addPredIsoROIs(mz = mzVals,
                                                int = intVals,
                                                scantime = xr@scantime,
                                                valsPerSpect = valsPerSpect,
                                                noise = 1500,
                                                features. = feats_1)
    all_f <- do_detectFeatures_centWaveWithPredIsoROIs(mz = mzVals,
                                                       int = intVals,
                                                       scantime = xr@scantime,
                                                       valsPerSpect = valsPerSpect,
                                                       noise = 1500)
    ## Comparisons.
    checkEquals(all_f, feats_2)
    old_all <- xcms:::.centWaveWithPredictedIsotopeROIs(xr, noise = 1500)
    checkEquals(all_f, old_all@.Data)
}

## Evaluate the featureDetection method using the centWaveWithPreIsoROIs method
## on OnDiskMSnExp and on MSnExp objects.
test_detectFeatures_centWaveWithPredIsoROIs <- function() {
    ## Control
    library(MSnbase)
    ##ppm <- 40
    snth <- 20
    ns <- 2500
    snthIso <- 5
    res_x <- findPeaks.centWaveWithPredictedIsotopeROIs(xr, noise = ns,
                                                        snthresh = snth,
                                                        snthreshIsoROIs = snthIso)@.Data
    ## Bypass xcmsRaw
    xs <- xcmsSet(fs[1], profparam = list(profstep = 0), snthresh = snth,
                  method = "centWaveWithPredictedIsotopeROIs", noise = ns,
                  snthreshIsoROIs = snthIso)
    checkEquals(xs@peaks[, colnames(res_x)], res_x)
    ## OnDiskMSnExp
    onDisk <- readMSData2(fs[1], msLevel. = 1)
    cwp <- CentWavePredIsoParam(snthresh = snth, noise = ns,
                                snthreshIsoROIs = snthIso)
    res <- detectFeatures(onDisk, param = cwp, return.type = "list")
    checkEquals(res[[1]], peaks(xs)@.Data)

    ## ## MSnExp
    ## inMem <- readMSData(fs[1], msLevel. = 1)
    ## res_2 <- detectFeatures(inMem, param = cwp, return.type = "list")
    ## checkEquals(res_2[[1]], peaks(xs)@.Data)

    ## returning an xcmsSet
    res <- detectFeatures(onDisk, param = cwp, return.type = "xcmsSet")
    checkEquals(peaks(res), peaks(xs))
    ## res <- detectFeatures(inMem, param = cwp, return.type = "xcmsSet")
    ## checkEquals(peaks(res), peaks(xs))

    ## Return an XCMSnExp
    res <- detectFeatures(onDisk, param = cwp)
    checkTrue(hasDetectedFeatures(res))
    checkTrue(!hasAdjustedRtime(res))
    checkTrue(!hasAlignedFeatures(res))
    checkEquals(peaks(xs)@.Data, features(res))

    ## Check on the full data.
    ## xs <- xcmsSet(fs, profparam = list(profstep = 0), snthresh = snth,
    ##               method = "centWaveWithPredictedIsotopeROIs", noise = ns,
    ##               snthreshIsoROIs = snthIso)
    ## onDisk <- readMSData2(fs, msLevel. = 1)
    ## res <- detectFeatures(onDisk, param = cwp)
    ## checkEquals(features(res), peaks(xs)@.Data)
}


## That is to test and evaluate the original code with the do_ code.
dontrun_test_impl_centWave_add <- function() {
    ## Using the do functions:
    valsPerSpect <- diff(c(xr@scanindex, length(mzVals)))
    do_1 <- do_detectFeatures_centWave(mz = mzVals, int = intVals,
                                       scantime = xr@scantime,
                                       valsPerSpect = valsPerSpect,
                                       verboseColumns = TRUE)
    do_2 <- do_detectFeatures_addPredIsoROIs(mz = mzVals, int = intVals,
                                             scantime = xr@scantime,
                                             valsPerSpect = valsPerSpect,
                                             features. = do_1)
    do_3 <- do_detectFeatures_centWaveWithPredIsoROIs(mz = mzVals, int = intVals,
                                                      scantime = xr@scantime,
                                                      valsPerSpect = valsPerSpect)
    checkEquals(do_2, do_3)
    checkTrue(nrow(do_1) < nrow(do_2))
    ## findPeaks using the do function.
    fp_1 <- findPeaks.centWave(xr, verbose.columns = TRUE)
    fp_2 <- findPeaks.addPredictedIsotopeFeatures(xr, xcmsPeaks = fp_1)
    fp_3 <- findPeaks.centWaveWithPredictedIsotopeROIs(xr)
    checkEquals(fp_2, fp_3)
    checkEquals(fp_2@.Data, do_2)
    ## Compare with original code:
    xs_2 <- xcms:::.addPredictedIsotopeFeatures(xr, xcmsPeaks = fp_1)
    xs_3 <- xcms:::.centWaveWithPredictedIsotopeROIs(xr)
    checkEquals(xs_2, xs_3)
    checkEquals(fp_2, xs_2)

    ##
    do_4 <- do_detectFeatures_centWaveWithPredIsoROIs(mz = mzVals, int = intVals,
                                                      scantime = xr@scantime,
                                                      valsPerSpect = valsPerSpect,
                                                      noise = 500,
                                                      maxIso = 7)
    fp_4 <- findPeaks.centWaveWithPredictedIsotopeROIs(xr, noise = 500,
                                                       maxiso = 7)
    xs_4 <- xcms:::.centWaveWithPredictedIsotopeROIs(xr, noise = 500, maxiso = 7)
    checkEquals(do_4, fp_4@.Data)
    checkEquals(do_4, xs_4@.Data)
}

## Test during implementation of do_define_isotopes to check whether results
## fit with the ones from the original code.
dontrun_test_do_define_isotopes <- function() {
    ## Subset to columns we want, transform to data.frame.
    req_cols <- c("mz", "mzmin", "mzmax", "scmin", "scmax", "intb", "scale")
    roiL <- split(res1, f = 1:nrow(res1))  ## is splitting a matrix faster?
    cns <- colnames(res1)
    roiL <- lapply(roiL, function(z) {
        names(z) <- cns
        return(z)
    })
    ## Original, returns a list which is appended to the original one.
    orig_ <- xcms:::do_define_isotopes_orig(roiList = roiL)
    new_ <- xcms:::do_define_isotopes(res1)
    checkIdentical(new_, orig_)
}

## Test during implementation of do_define_adducts to check whether results
## fit with the ones from the original code.
dontrun_test_do_define_adducts <- function() {
    ## Subset to columns we want, transform to data.frame.
    req_cols <- c("mz", "mzmin", "mzmax", "scmin", "scmax", "intb", "scale")
    roiL <- split(res1, f = 1:nrow(res1))  ## is splitting a matrix faster?
    cns <- colnames(res1)
    roiL <- lapply(roiL, function(z) {
        names(z) <- cns
        return(z)
    })
    ## Original, returns a list which is appended to the original one.
    orig_ <- xcms:::do_define_adducts_orig(roiList = roiL)
    new_ <- xcms:::do_define_adducts(res1)
    checkIdentical(new_, orig_)
}

