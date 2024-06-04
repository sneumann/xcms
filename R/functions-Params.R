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
    if (is(x, "LamaParama"))
        return("lama")
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
                          firstBaselineCheck = TRUE, roiScales = numeric(),
                          extendLengthMSW = FALSE, verboseBetaColumns = FALSE) {
    return(new("CentWaveParam", ppm = ppm, peakwidth = peakwidth,
               snthresh = snthresh, prefilter = prefilter,
               mzCenterFun = mzCenterFun, integrate = as.integer(integrate),
               mzdiff = mzdiff, fitgauss = fitgauss, noise = noise,
               verboseColumns = verboseColumns, roiList = roiList,
               firstBaselineCheck = firstBaselineCheck, roiScales = roiScales,
               extendLengthMSW = extendLengthMSW,
               verboseBetaColumns=verboseBetaColumns))
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
#'     Savitzky-Golay filter is used (in the \code{MassSpecWavelet})
#'     package to estimate the local intensity.
#'
#' @param tuneIn logical(1) whther to tune in the parameter
#'     estimation of the detected peaks.
#'
#' @param ... Additional parameters to be passed to the
#'     \code{\link{peakDetectionCWT}} and
#'     \code{\link{identifyMajorPeaks}} functions from the
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
                          extendLengthMSW = FALSE, verboseBetaColumns = FALSE,
                          snthreshIsoROIs = 6.25, maxCharge = 3, maxIso = 5,
                          mzIntervalExtension = TRUE, polarity = "unknown") {
    return(new("CentWavePredIsoParam", ppm = ppm, peakwidth = peakwidth,
               snthresh = snthresh, prefilter = prefilter,
               mzCenterFun = mzCenterFun, integrate = as.integer(integrate),
               mzdiff = mzdiff, fitgauss = fitgauss, noise = noise,
               verboseColumns = verboseColumns, roiList = roiList,
               firstBaselineCheck = firstBaselineCheck, roiScales = roiScales,
               extendLengthMSW = extendLengthMSW,
               verboseBetaColumns = verboseBetaColumns,
               snthreshIsoROIs = snthreshIsoROIs, maxIso = as.integer(maxIso),
               maxCharge = as.integer(maxCharge),
               mzIntervalExtension = mzIntervalExtension, polarity = polarity))
}

#' @rdname groupChromPeaks
PeakDensityParam <- function(sampleGroups = numeric(), bw = 30,
                             minFraction = 0.5, minSamples = 1,
                             binSize = 0.25, ppm = 0, maxFeatures = 50) {
    if (length(sampleGroups) == 0)
        stop("Argument 'sampleGroups' has to be defined.")
    new("PeakDensityParam", sampleGroups = sampleGroups, bw = bw,
        minFraction = minFraction, minSamples = minSamples,
        binSize = binSize, ppm = ppm, maxFeatures = maxFeatures)
}

#' @rdname groupChromPeaks
MzClustParam <- function(sampleGroups = numeric(), ppm = 20, absMz = 0,
                                minFraction = 0.5, minSamples = 1) {
    return(new("MzClustParam", sampleGroups = sampleGroups, ppm = ppm,
               absMz = absMz, minFraction = minFraction,
               minSamples = minSamples))
}

#' @rdname groupChromPeaks
NearestPeaksParam <- function(sampleGroups = numeric(), mzVsRtBalance = 10,
                              absMz = 0.2, absRt = 15, kNN = 10) {
    return(new("NearestPeaksParam", sampleGroups = sampleGroups,
               mzVsRtBalance = mzVsRtBalance, absMz = absMz, absRt = absRt,
               kNN = kNN))
}

#' @rdname adjustRtime
PeakGroupsParam <- function(minFraction = 0.9, extraPeaks = 1,
                               smooth = "loess", span = 0.2,
                            family = "gaussian",
                            peakGroupsMatrix = matrix(nrow = 0, ncol = 0),
                            subset = integer(),
                            subsetAdjust = c("average", "previous")) {
    subsetAdjust <- match.arg(subsetAdjust)
    new("PeakGroupsParam", minFraction = minFraction,
        extraPeaks = extraPeaks, smooth = smooth, span = span,
        family = family, peakGroupsMatrix = peakGroupsMatrix,
        subset = as.integer(subset), subsetAdjust = subsetAdjust)
}

#' @rdname LamaParama
LamaParama <- function(lamas = matrix(ncol = 2, nrow = 0,
                                      dimnames = list(NULL, c("mz", "rt"))),
                       method = c("loess", "gam"),
                       span = 0.5,
                       outlierTolerance = 3,
                       zeroWeight = 10,
                       ppm = 20,
                       tolerance = 0,
                       toleranceRt = 5,
                       bs = "tp") {
    method <- match.arg(method)
    if (method == "gam")
        .check_gam_library()
    if (is.data.frame(lamas))
        lamas <- as.matrix(lamas)
    if (ncol(lamas) != 2)
        stop("the 'lamas' matrix needs to have two columns, composed of m/z, ",
        "and retention time of the peaks from the reference dataset, in this ",
        "order")
    new("LamaParama", lamas = lamas,
        method = method,
        span = span,
        outlierTolerance = outlierTolerance,
        zeroWeight = zeroWeight,
        ppm = ppm,
        tolerance = tolerance,
        toleranceRt = toleranceRt,
        bs = bs)
}

#' @rdname adjustRtime
ObiwarpParam <- function(binSize = 1, centerSample = integer(), response = 1L,
                         distFun = "cor_opt", gapInit = numeric(),
                         gapExtend = numeric(), factorDiag = 2, factorGap = 1,
                         localAlignment = FALSE, initPenalty = 0,
                         subset = integer(),
                         subsetAdjust = c("average", "previous"),
                         rtimeDifferenceThreshold = 5) {
    subsetAdjust <- match.arg(subsetAdjust)
    new("ObiwarpParam", binSize = binSize,
        centerSample = as.integer(centerSample),
        response = as.integer(response), distFun = distFun,
        gapInit = gapInit, gapExtend = gapExtend, factorDiag = factorDiag,
        factorGap = factorGap, localAlignment = localAlignment,
        initPenalty = initPenalty, subset = as.integer(subset),
        subsetAdjust = subsetAdjust,
        rtimeDifferenceThreshold = rtimeDifferenceThreshold[1L])
}

#' @return The \code{FillChromPeaksParam} function returns a
#'     \code{FillChromPeaksParam} object.
#'
#' @rdname fillChromPeaks
FillChromPeaksParam <- function(expandMz = 0, expandRt = 0, ppm = 0,
                                fixedMz = 0, fixedRt = 0) {
    new("FillChromPeaksParam", expandMz = expandMz, expandRt = expandRt,
        ppm = ppm, fixedMz = fixedMz, fixedRt = fixedRt)
}

#' @rdname fillChromPeaks
fixedRt <- function(object) object@fixedRt

#' @rdname fillChromPeaks
fixedMz <- function(object) object@fixedMz

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

#' @rdname refineChromPeaks
CleanPeaksParam <- function(maxPeakwidth = 10) {
    new("CleanPeaksParam", maxPeakwidth = as.numeric(maxPeakwidth))
}

#' @rdname refineChromPeaks
MergeNeighboringPeaksParam <- function(expandRt = 2, expandMz = 0, ppm = 10,
                                       minProp = 0.75) {
    new("MergeNeighboringPeaksParam", expandRt = as.numeric(expandRt),
        expandMz = as.numeric(expandMz), ppm = as.numeric(ppm),
        minProp = as.numeric(minProp))
}

#' @rdname fillChromPeaks
ChromPeakAreaParam <-
    function(mzmin = function(z) quantile(z, probs = 0.25, names = FALSE),
             mzmax = function(z) quantile(z, probs = 0.75, names = FALSE),
             rtmin = function(z) quantile(z, probs = 0.25, names = FALSE),
             rtmax = function(z) quantile(z, probs = 0.75, names = FALSE)) {
        new("ChromPeakAreaParam", mzmin = mzmin, mzmax = mzmax, rtmin = rtmin,
            rtmax = rtmax)
}

#' @rdname refineChromPeaks
FilterIntensityParam <- function(threshold = 0, nValues = 1L, value = "maxo") {
    new("FilterIntensityParam", threshold = as.numeric(threshold),
        nValues = as.integer(nValues), value = value)
}
