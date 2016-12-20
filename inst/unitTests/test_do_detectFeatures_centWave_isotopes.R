## Tests related to the findPeaks.centWaveWithAddIsotopeROIs.
library(faahKO)
fs <- c(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
        system.file('cdf/KO/ko16.CDF', package = "faahKO"),
        system.file('cdf/KO/ko18.CDF', package = "faahKO"),
        system.file('cdf/KO/ko19.CDF', package = "faahKO"))
xr <- xcmsRaw(fs[1], profstep = 0)
mzVals <- xr@env$mz
intVals <- xr@env$intensity

test_detectFeatures_centWaveWithPredIsoROIs <- function() {
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

dontrun_test_impl <- function() {
    ## Default parameter:
    ppm <- 25
    peakwidth <- c(20, 50)
    prefilter <- c(3, 100)
    mzCenterFun <- "wMean"
    integrate <- 1
    mzdiff <- -0.001
    fitgauss <- FALSE
    noise <- 0
    snthresh <- 10  ## That's for the first centWave.
    snthreshIsoROIs <- 6.25  ## That"s for the second centWave.
    verboseColumns <- FALSE
    firstBaselineCheck <- TRUE
    ## New ones
    maxCharge <- 3
    maxIso <- 5
    mzIntervalExtension <- TRUE

    ## First file
    xr <- xcmsRaw(fs[1], profstep = 0)
    mzVals <- xr@env$mz
    intVals <- xr@env$intensity
    ## Define the values per spectrum:
    valsPerSpect <- diff(c(xr@scanindex, length(mzVals)))

    orig <- xcms:::.centWaveWithPredictedIsotopeROIs(xr, ppm = ppm,
                                                       peakwidth = peakwidth,
                                                       snthresh = snthresh,
                                                       prefilter = prefilter,
                                                       mzdiff = mzdiff,
                                                       fitgauss = fitgauss,
                                                       noise = noise,
                                                       snthreshIsoROIs = snthreshIsoROIs,
                                                       maxcharge = maxCharge,
                                                       maxiso = maxIso,
                                                       mzIntervalExtension = TRUE)
    d <- do_detectFeatures_centWaveWithPredIsoROIs(mz = mzVals,
                                                   int = intVals,
                                                   scantime = xr@scantime,
                                                   valsPerSpect = valsPerSpect,
                                                   ppm = ppm,
                                                   peakwidth = peakwidth,
                                                   snthresh = snthresh,
                                                   prefilter = prefilter,
                                                   mzdiff = mzdiff,
                                                   fitgauss = fitgauss,
                                                   noise = noise,
                                                   snthreshIsoROIs = snthreshIsoROIs,
                                                   maxCharge = maxCharge,
                                                   maxIso = maxIso,
                                                   mzIntervalExtension = mzIntervalExtension)
    checkEquals(orig@.Data, d)
    ## Modify settings
    snthresh <- 50
    snthreshIsoROIs <- 12
    orig <- xcms:::.centWaveWithPredictedIsotopeROIs(xr, ppm = ppm,
                                                     peakwidth = peakwidth,
                                                     snthresh = snthresh,
                                                     prefilter = prefilter,
                                                     mzdiff = mzdiff,
                                                     fitgauss = fitgauss,
                                                     noise = noise,
                                                     snthreshIsoROIs = snthreshIsoROIs,
                                                     maxcharge = maxCharge,
                                                     maxiso = maxIso,
                                                     mzIntervalExtension = TRUE)
    d <- do_detectFeatures_centWaveWithPredIsoROIs(mz = mzVals,
                                                   int = intVals,
                                                   scantime = xr@scantime,
                                                   valsPerSpect = valsPerSpect,
                                                   ppm = ppm,
                                                   peakwidth = peakwidth,
                                                   snthresh = snthresh,
                                                   prefilter = prefilter,
                                                   mzdiff = mzdiff,
                                                   fitgauss = fitgauss,
                                                   noise = noise,
                                                   snthreshIsoROIs = snthreshIsoROIs,
                                                   maxCharge = maxCharge,
                                                   maxIso = maxIso,
                                                   mzIntervalExtension = mzIntervalExtension)
    checkEquals(orig@.Data, d)
    ##
    maxCharge <- 5
    maxIso <- 8
    orig <- xcms:::.centWaveWithPredictedIsotopeROIs(xr, ppm = ppm,
                                                     peakwidth = peakwidth,
                                                     snthresh = snthresh,
                                                     prefilter = prefilter,
                                                     mzdiff = mzdiff,
                                                     fitgauss = fitgauss,
                                                     noise = noise,
                                                     snthreshIsoROIs = snthreshIsoROIs,
                                                     maxcharge = maxCharge,
                                                     maxiso = maxIso,
                                                     mzIntervalExtension = TRUE)
    d <- do_detectFeatures_centWaveWithPredIsoROIs(mz = mzVals,
                                                   int = intVals,
                                                   scantime = xr@scantime,
                                                   valsPerSpect = valsPerSpect,
                                                   ppm = ppm,
                                                   peakwidth = peakwidth,
                                                   snthresh = snthresh,
                                                   prefilter = prefilter,
                                                   mzdiff = mzdiff,
                                                   fitgauss = fitgauss,
                                                   noise = noise,
                                                   snthreshIsoROIs = snthreshIsoROIs,
                                                   maxCharge = maxCharge,
                                                   maxIso = maxIso,
                                                   mzIntervalExtension = mzIntervalExtension)
    checkEquals(orig@.Data, d)
    ##
    mzIntervalExtension <- FALSE
    orig <- xcms:::.centWaveWithPredictedIsotopeROIs(xr, ppm = ppm,
                                                     peakwidth = peakwidth,
                                                     snthresh = snthresh,
                                                     prefilter = prefilter,
                                                     mzdiff = mzdiff,
                                                     fitgauss = fitgauss,
                                                     noise = noise,
                                                     snthreshIsoROIs = snthreshIsoROIs,
                                                     maxcharge = maxCharge,
                                                     maxiso = maxIso,
                                                     mzIntervalExtension = mzIntervalExtension)
    d <- do_detectFeatures_centWaveWithPredIsoROIs(mz = mzVals,
                                                   int = intVals,
                                                   scantime = xr@scantime,
                                                   valsPerSpect = valsPerSpect,
                                                   ppm = ppm,
                                                   peakwidth = peakwidth,
                                                   snthresh = snthresh,
                                                   prefilter = prefilter,
                                                   mzdiff = mzdiff,
                                                   fitgauss = fitgauss,
                                                   noise = noise,
                                                   snthreshIsoROIs = snthreshIsoROIs,
                                                   maxCharge = maxCharge,
                                                   maxIso = maxIso,
                                                   mzIntervalExtension = mzIntervalExtension)
    checkEquals(orig@.Data, d)
}

