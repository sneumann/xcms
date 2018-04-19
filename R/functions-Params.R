## Functions related to the Param class and sub-classes.
#' @include DataClasses.R

##
#' @description Extract all slot values and put them into a list, names being
#'     the slot names. If a slot \code{addParams} exist its content will be
#'     appended to the returned list.
#'
#' @param x A Param class.
#' 
#' @author Johannes Rainer
#' 
#' @noRd
.param2list <- function(x) {
    ## Get all slot names, skip those matching the provided pattern.
    sNames <- slotNames(x)
    skipSome <- grep(sNames, pattern = "^\\.")
    if(length(skipSome) > 0)
        sNames <- sNames[-skipSome]
    ## handle a slot called "addParams" differently: this is thougth to contain
    ## ... arguments thus we have to skip this one too.
    if (any(sNames == "addParams")) {
        sNames <- sNames[sNames != "addParams"]
        addP <- x@addParams
    } else {
        addP <- list()
    }
    if(length(sNames) > 0){
        resL <- vector("list", length(sNames))
        for(i in 1:length(sNames))
            resL[[i]] <- slot(x, name = sNames[i])
        names(resL) <- sNames
        resL <- c(resL, addP)
        return(resL)
    }else{
          return(list())
    }
}

## Just get the name of the algorithm for each Parameter class.
.param2string <- function(x) {
    if (is(x, "CentWaveParam"))
        return("centWave")
    if (is(x, "MatchedFilterParam"))
        return("matchedFilter")
    if (is(x, "MassifquantParam"))
        return("massifquant")
    if (is(x, "MSWParam"))
        return("MSW")
    if (is(x, "CentWavePredIsoParam"))
        return("centWave with predicted isotope ROIs")
    if (is(x, "PeakDensityParam"))
        return("chromatographic peak density")
    if (is(x, "MzClustParam"))
        return("mzClust")
    if (is(x, "NearestPeaksParam"))
        return("nearest peaks")
    if (is(x, "PeakGroupsParam"))
        return("peak groups")
    if (is(x, "ObiwarpParam"))
        return("obiwarp")
    return("unknown")
}

############################################################
## GenericParam
#' @return The \code{GenericParam} function returns a \code{GenericParam}
#'     object.
#' 
#' @param fun \code{character} representing the name of the function.
#' 
#' @param args \code{list} (ideally named) with the arguments to the function.
#' 
#' @rdname GenericParam
GenericParam <- function(fun = character(), args = list()) {
    return(new("GenericParam", fun = fun, args = args))
}

#' @return The \code{CentWaveParam} function returns a \code{CentWaveParam}
#'     class instance with all of the settings specified for chromatographic
#'     peak detection by the centWave method.
#'
#' @rdname findChromPeaks-centWave
CentWaveParam <- function(ppm = 25, peakwidth = c(20, 50), snthresh = 10,
                          prefilter = c(3, 100), mzCenterFun = "wMean",
                          integrate = 1L, mzdiff = -0.001, fitgauss = FALSE,
                          noise = 0, verboseColumns = FALSE, roiList = list(),
                          firstBaselineCheck = TRUE, roiScales = numeric()) {
    return(new("CentWaveParam", ppm = ppm, peakwidth = peakwidth,
               snthresh = snthresh, prefilter = prefilter,
               mzCenterFun = mzCenterFun, integrate = as.integer(integrate),
               mzdiff = mzdiff, fitgauss = fitgauss, noise = noise,
               verboseColumns = verboseColumns, roiList = roiList,
               firstBaselineCheck = firstBaselineCheck, roiScales = roiScales))
}

#' @return The \code{MatchedFilterParam} function returns a
#'     \code{MatchedFilterParam} class instance with all of the settings
#'     specified for chromatographic detection by the \emph{matchedFilter}
#'     method.
#'
#' @rdname findChromPeaks-matchedFilter
MatchedFilterParam <- function(binSize = 0.1, impute = "none",
                               baseValue = numeric(), distance = numeric(),
                               fwhm = 30, sigma = fwhm / 2.3548,
                               max = 5, snthresh = 10, steps = 2,
                               mzdiff = 0.8 - binSize * steps, index = FALSE) {
    return(new("MatchedFilterParam", binSize = binSize, impute = impute,
               baseValue = baseValue, distance = distance, fwhm = fwhm,
               sigma = sigma, max = max, snthresh = snthresh, steps = steps,
               mzdiff = mzdiff, index = index))
}
#' Convert the impute method to the old-style method name (e.g. for profMat
#'     calls)
#' 
#' @noRd
.impute2method <- function(x) {
    if (impute(x) == "none")
        return("bin")
    if (impute(x) == "lin")
        return("binlin")
    if (impute(x) == "linbase")
        return("binlinbase")
    return("intlin")
}

#' @return The \code{MassifquantParam} function returns a
#'     \code{MassifquantParam} class instance with all of the settings
#'     specified for chromatographic peak detection by the \emph{massifquant}
#'     method.
#'
#' @rdname findChromPeaks-massifquant
MassifquantParam <- function(ppm = 25, peakwidth = c(20, 50), snthresh = 10,
                             prefilter = c(3, 100), mzCenterFun = "wMean",
                             integrate = 1L, mzdiff = -0.001, fitgauss = FALSE,
                             noise = 0, verboseColumns = FALSE,
                             criticalValue = 1.125, consecMissedLimit = 2,
                             unions = 1, checkBack = 0, withWave = FALSE) {
    return(new("MassifquantParam", ppm = ppm, peakwidth = peakwidth,
               snthresh = snthresh, prefilter = prefilter,
               mzCenterFun = mzCenterFun, integrate = as.integer(integrate),
               mzdiff = mzdiff, fitgauss = fitgauss, noise = noise,
               verboseColumns = verboseColumns, criticalValue = criticalValue,
               consecMissedLimit = as.integer(consecMissedLimit),
               unions = as.integer(unions), checkBack = as.integer(checkBack),
               withWave = withWave))
}

#' @inheritParams findChromPeaks-centWave
#'
#' @param scales Numeric defining the scales of the continuous wavelet
#'     transform (CWT).
#'
#' @param nearbyPeak logical(1) whether to include nearby peaks of
#'     major peaks.
#'
#' @param peakScaleRange numeric(1) defining the scale range of the
#'     peak (larger than 5 by default).
#'
#' @param ampTh numeric(1) defining the minimum required relative
#'     amplitude of the peak (ratio of the maximum of CWT coefficients).
#'
#' @param minNoiseLevel numeric(1) defining the minimum noise level
#'     used in computing the SNR.
#'
#' @param ridgeLength numeric(1) defining the minimum highest scale
#'     of the peak in 2-D CWT coefficient matrix.
#'
#' @param peakThr numeric(1) with the minimum absolute intensity
#'     (above baseline) of peaks to be picked. If provided, the smoothing
#'     function \code{\link{sav.gol}} function (in the \code{MassSpecWavelet})
#'     package is called to estimate the local intensity.
#'
#' @param tuneIn logical(1) whther to tune in the parameter
#'     estimation of the detected peaks.
#'
#' @param ... Additional parameters to be passed to the
#'     \code{\link{identifyMajorPeaks}} and
#'     \code{\link{sav.gol}} functions from the
#'     \code{MassSpecWavelet} package.
#'
#' @return The \code{MSWParam} function returns a \code{MSWParam}
#'     class instance with all of the settings specified for peak detection by
#'     the \emph{MSW} method.
#'
#' @rdname findPeaks-MSW
MSWParam <- function(snthresh = 3, verboseColumns = FALSE,
                     scales = c(1, seq(2, 30, 2), seq(32, 64, 4)),
                     nearbyPeak = TRUE, peakScaleRange = 5,
                     ampTh = 0.01, minNoiseLevel = ampTh / snthresh,
                     ridgeLength = 24, peakThr = NULL, tuneIn = FALSE,
                     ... ) {
    addParams <- list(...)
    if (is.null(peakThr))
        peakThr <- numeric()
    return(new("MSWParam", snthresh = snthresh, verboseColumns = verboseColumns,
               scales = scales, nearbyPeak = nearbyPeak,
               peakScaleRange = peakScaleRange, ampTh = ampTh,
               minNoiseLevel = minNoiseLevel, ridgeLength = ridgeLength,
               peakThr = peakThr, tuneIn = tuneIn, addParams = addParams))
}

#' @return The \code{CentWavePredIsoParam} function returns a
#'     \code{CentWavePredIsoParam} class instance with all of the settings
#'     specified for the two-step centWave-based peak detection considering also
#'     isotopes.
#'
#' @rdname findChromPeaks-centWaveWithPredIsoROIs
CentWavePredIsoParam <- function(ppm = 25, peakwidth = c(20, 50), snthresh = 10,
                          prefilter = c(3, 100), mzCenterFun = "wMean",
                          integrate = 1L, mzdiff = -0.001, fitgauss = FALSE,
                          noise = 0, verboseColumns = FALSE, roiList = list(),
                          firstBaselineCheck = TRUE, roiScales = numeric(),
                          snthreshIsoROIs = 6.25, maxCharge = 3, maxIso = 5,
                          mzIntervalExtension = TRUE, polarity = "unknown") {
    return(new("CentWavePredIsoParam", ppm = ppm, peakwidth = peakwidth,
               snthresh = snthresh, prefilter = prefilter,
               mzCenterFun = mzCenterFun, integrate = as.integer(integrate),
               mzdiff = mzdiff, fitgauss = fitgauss, noise = noise,
               verboseColumns = verboseColumns, roiList = roiList,
               firstBaselineCheck = firstBaselineCheck, roiScales = roiScales,
               snthreshIsoROIs = snthreshIsoROIs, maxIso = as.integer(maxIso),
               maxCharge = as.integer(maxCharge),
               mzIntervalExtension = mzIntervalExtension, polarity = polarity))
}

#' @return The \code{PeakDensityParam} function returns a
#'     \code{PeakDensityParam} class instance with all of the settings
#'     specified for chromatographic peak alignment based on peak densities.
#'     Note that argument `sampleGroups` is mandatory and should represent
#'     either the sample grouping in the experiment. It's length has to match
#'     the number of sample in the experiments.
#' 
#' @rdname groupChromPeaks-density
PeakDensityParam <- function(sampleGroups = numeric(), bw = 30,
                                minFraction = 0.5, minSamples = 1,
                             binSize = 0.25, maxFeatures = 50) {
    if (length(sampleGroups) == 0 | any(is.na(sampleGroups)))
        stop("Argument 'sampleGroups' has to be defined. It should not ",
             "contain 'NA's")
    new("PeakDensityParam", sampleGroups = sampleGroups, bw = bw,
        minFraction = minFraction, minSamples = minSamples,
        binSize = binSize, maxFeatures = maxFeatures)
}

#' @return The \code{MzClustParam} function returns a
#'     \code{MzClustParam} class instance with all of the settings
#'     specified for high resolution single spectra peak alignment.
#' 
#' @rdname groupChromPeaks-mzClust
MzClustParam <- function(sampleGroups = numeric(), ppm = 20, absMz = 0,
                                minFraction = 0.5, minSamples = 1) {
    return(new("MzClustParam", sampleGroups = sampleGroups, ppm = ppm,
               absMz = absMz, minFraction = minFraction,
               minSamples = minSamples))
}

#' @return The \code{NearestPeaksParam} function returns a
#'     \code{NearestPeaksParam} class instance with all of the settings
#'     specified for peak alignment based on peak proximity.
#' 
#' @rdname groupChromPeaks-nearest
NearestPeaksParam <- function(sampleGroups = numeric(), mzVsRtBalance = 10,
                              absMz = 0.2, absRt = 15, kNN = 10) {
    return(new("NearestPeaksParam", sampleGroups = sampleGroups,
               mzVsRtBalance = mzVsRtBalance, absMz = absMz, absRt = absRt,
               kNN = kNN))
}

#' @return The \code{PeakGroupsParam} function returns a
#'     \code{PeakGroupsParam} class instance with all of the settings
#'     specified for retention time adjustment based on \emph{house keeping}
#'     features/peak groups.
#' 
#' @rdname adjustRtime-peakGroups
PeakGroupsParam <- function(minFraction = 0.9, extraPeaks = 1,
                               smooth = "loess", span = 0.2,
                            family = "gaussian",
                            peakGroupsMatrix = matrix(nrow = 0, ncol = 0)) {
    return(new("PeakGroupsParam", minFraction = minFraction,
               extraPeaks = extraPeaks, smooth = smooth, span = span,
               family = family, peakGroupsMatrix = peakGroupsMatrix))
}

#' @return The \code{ObiwarpParam} function returns a
#'     \code{ObiwarpParam} class instance with all of the settings
#'     specified for obiwarp retention time adjustment and alignment.
#' 
#' @rdname adjustRtime-obiwarp
ObiwarpParam <- function(binSize = 1, centerSample = integer(), response = 1L,
                         distFun = "cor_opt", gapInit = numeric(),
                         gapExtend = numeric(), factorDiag = 2, factorGap = 1,
                         localAlignment = FALSE, initPenalty = 0) {
    return(new("ObiwarpParam", binSize = binSize,
               centerSample = as.integer(centerSample),
               response = as.integer(response), distFun = distFun,
               gapInit = gapInit, gapExtend = gapExtend, factorDiag = factorDiag,
               factorGap = factorGap, localAlignment = localAlignment,
               initPenalty = initPenalty))
}

#' @return The \code{FillChromPeaksParam} function returns a
#'     \code{FillChromPeaksParam} object.
#' 
#' @rdname fillChromPeaks
FillChromPeaksParam <- function(expandMz = 0, expandRt = 0, ppm = 0) {
    return(new("FillChromPeaksParam", expandMz = expandMz, expandRt = expandRt,
               ppm = ppm))
}


#' @return The `CalibrantMassParam` function returns an instance of
#'     the `CalibrantMassParam` class with all settings and properties set.
#'
#' @md
#' 
#' @rdname calibrate-calibrant-mass
CalibrantMassParam <- function(mz = list(), mzabs = 0.0001, mzppm = 5,
                               neighbors = 3, method = "linear") {
    if (!is.list(mz))
        mz <- list(mz)
    mz <- lapply(mz, sort)
    new("CalibrantMassParam", mz = mz, mzabs = mzabs, mzppm = mzppm,
        neighbors = as.integer(neighbors), method = method)
}

.mzabs <- function(x)
    x@mzabs

.mzppm <- function(x)
    x@mzppm

.neighbors <- function(x)
    x@neighbors

.method <- function(x)
    x@method

.mz <- function(x)
    x@mz
