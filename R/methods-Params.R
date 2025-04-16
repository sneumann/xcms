## Methods for the Param class and sub-classes
#' @include functions-Params.R

#' @rdname hidden_aliases
setMethod("ppm", "CentWaveParam", function(object){ return(object@ppm)})
#' @rdname hidden_aliases
setReplaceMethod("ppm", "CentWaveParam", function(object, value) {
    object@ppm <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("peakwidth", "CentWaveParam", function(object)
    return(object@peakwidth))
#' @rdname hidden_aliases
setReplaceMethod("peakwidth", "CentWaveParam", function(object, value) {
    object@peakwidth <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("snthresh", "CentWaveParam", function(object)
    return(object@snthresh))
#' @rdname hidden_aliases
setReplaceMethod("snthresh", "CentWaveParam", function(object, value) {
    object@snthresh <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("prefilter", "CentWaveParam", function(object)
    return(object@prefilter))
#' @rdname hidden_aliases
setReplaceMethod("prefilter", "CentWaveParam", function(object, value) {
    object@prefilter <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("mzCenterFun", "CentWaveParam", function(object)
    return(object@mzCenterFun))
#' @rdname hidden_aliases
setReplaceMethod("mzCenterFun", "CentWaveParam", function(object, value) {
    object@mzCenterFun <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("integrate", signature(f = "CentWaveParam"), function(f)
    return(f@integrate))
#' @rdname hidden_aliases
setReplaceMethod("integrate", "CentWaveParam", function(object, value) {
    object@integrate <- as.integer(value)
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("mzdiff", "CentWaveParam", function(object)
    return(object@mzdiff))
#' @rdname hidden_aliases
setReplaceMethod("mzdiff", "CentWaveParam", function(object, value) {
    object@mzdiff <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("fitgauss", "CentWaveParam", function(object)
    return(object@fitgauss))
#' @rdname hidden_aliases
setReplaceMethod("fitgauss", "CentWaveParam", function(object, value) {
    object@fitgauss <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("noise", "CentWaveParam", function(object)
    return(object@noise))
#' @rdname hidden_aliases
setReplaceMethod("noise", "CentWaveParam", function(object, value) {
    object@noise <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("verboseColumns", "CentWaveParam", function(object)
    return(object@verboseColumns))
#' @rdname hidden_aliases
setReplaceMethod("verboseColumns", "CentWaveParam", function(object, value) {
    object@verboseColumns <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("roiList", "CentWaveParam", function(object)
    return(object@roiList))
#' @rdname hidden_aliases
setReplaceMethod("roiList", "CentWaveParam", function(object, value) {
    object@roiList <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("firstBaselineCheck", "CentWaveParam", function(object)
    return(object@firstBaselineCheck))
#' @rdname hidden_aliases
setReplaceMethod("firstBaselineCheck", "CentWaveParam", function(object, value) {
    object@firstBaselineCheck <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("roiScales", "CentWaveParam", function(object)
    return(object@roiScales))
#' @rdname hidden_aliases
setReplaceMethod("roiScales", "CentWaveParam", function(object, value) {
    object@roiScales <- value
    if (validObject(object))
        return(object)
})


############################################################
## MatchedFilterParam

#' @rdname hidden_aliases
setMethod("binSize", "MatchedFilterParam", function(object)
    return(object@binSize))
#' @rdname hidden_aliases
setReplaceMethod("binSize", "MatchedFilterParam", function(object, value) {
    object@binSize <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("impute", "MatchedFilterParam", function(object)
    return(object@impute))
#' @rdname hidden_aliases
setReplaceMethod("impute", "MatchedFilterParam", function(object, value) {
    object@impute <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("baseValue", "MatchedFilterParam", function(object)
    return(object@baseValue))
#' @rdname hidden_aliases
setReplaceMethod("baseValue", "MatchedFilterParam", function(object, value) {
    object@baseValue <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("distance", "MatchedFilterParam", function(object)
    return(object@distance))
#' @rdname hidden_aliases
setReplaceMethod("distance", "MatchedFilterParam", function(object, value) {
    object@distance <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("fwhm", "MatchedFilterParam", function(object)
    return(object@fwhm))
#' @rdname hidden_aliases
setReplaceMethod("fwhm", "MatchedFilterParam", function(object, value) {
    object@fwhm <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("sigma", "MatchedFilterParam", function(object)
    return(object@sigma))
#' @rdname hidden_aliases
setReplaceMethod("sigma", "MatchedFilterParam", function(object, value) {
    object@sigma <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("max", signature(x="MatchedFilterParam"),
          function(x) return(x@max))
#' @rdname hidden_aliases
setReplaceMethod("max", "MatchedFilterParam", function(object, value) {
    object@max <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("snthresh", "MatchedFilterParam", function(object)
    return(object@snthresh))
#' @rdname hidden_aliases
setReplaceMethod("snthresh", "MatchedFilterParam", function(object, value) {
    object@snthresh <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("steps", "MatchedFilterParam", function(object)
    return(object@steps))
#' @rdname hidden_aliases
setReplaceMethod("steps", "MatchedFilterParam", function(object, value) {
    object@steps <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("mzdiff", "MatchedFilterParam", function(object)
    return(object@mzdiff))
#' @rdname hidden_aliases
setReplaceMethod("mzdiff", "MatchedFilterParam", function(object, value) {
    object@mzdiff <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("index", "MatchedFilterParam", function(object)
    return(object@index))
#' @rdname hidden_aliases
setReplaceMethod("index", "MatchedFilterParam", function(object, value) {
    object@index <- value
    if (validObject(object))
        return(object)
})

############################################################
## MassifquantParam
###

#' @rdname hidden_aliases
setMethod("ppm", "MassifquantParam", function(object){ return(object@ppm)})
#' @rdname hidden_aliases
setReplaceMethod("ppm", "MassifquantParam", function(object, value) {
    object@ppm <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("peakwidth", "MassifquantParam", function(object)
    return(object@peakwidth))
#' @rdname hidden_aliases
setReplaceMethod("peakwidth", "MassifquantParam", function(object, value) {
    object@peakwidth <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("snthresh", "MassifquantParam", function(object)
    return(object@snthresh))
#' @rdname hidden_aliases
setReplaceMethod("snthresh", "MassifquantParam", function(object, value) {
    object@snthresh <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("prefilter", "MassifquantParam", function(object)
    return(object@prefilter))
#' @rdname hidden_aliases
setReplaceMethod("prefilter", "MassifquantParam", function(object, value) {
    object@prefilter <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("mzCenterFun", "MassifquantParam", function(object)
    return(object@mzCenterFun))
#' @rdname hidden_aliases
setReplaceMethod("mzCenterFun", "MassifquantParam", function(object, value) {
    object@mzCenterFun <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("integrate", signature(f = "MassifquantParam"), function(f)
    return(f@integrate))
#' @rdname hidden_aliases
setReplaceMethod("integrate", "MassifquantParam", function(object, value) {
    object@integrate <- as.integer(value)
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("mzdiff", "MassifquantParam", function(object)
    return(object@mzdiff))
#' @rdname hidden_aliases
setReplaceMethod("mzdiff", "MassifquantParam", function(object, value) {
    object@mzdiff <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("fitgauss", "MassifquantParam", function(object)
    return(object@fitgauss))
#' @rdname hidden_aliases
setReplaceMethod("fitgauss", "MassifquantParam", function(object, value) {
    object@fitgauss <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("noise", "MassifquantParam", function(object)
    return(object@noise))
#' @rdname hidden_aliases
setReplaceMethod("noise", "MassifquantParam", function(object, value) {
    object@noise <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("verboseColumns", "MassifquantParam", function(object)
    return(object@verboseColumns))
#' @rdname hidden_aliases
setReplaceMethod("verboseColumns", "MassifquantParam", function(object, value) {
    object@verboseColumns <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("criticalValue", "MassifquantParam", function(object)
    return(object@criticalValue))
#' @rdname hidden_aliases
setReplaceMethod("criticalValue", "MassifquantParam", function(object, value) {
    object@criticalValue <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("consecMissedLimit", "MassifquantParam", function(object)
    return(object@consecMissedLimit))
#' @rdname hidden_aliases
setReplaceMethod("consecMissedLimit", "MassifquantParam",
                 function(object, value) {
                     object@consecMissedLimit <- as.integer(value)
                     if (validObject(object))
                         return(object)
                 })

#' @rdname hidden_aliases
setMethod("unions", "MassifquantParam", function(object)
    return(object@unions))
#' @rdname hidden_aliases
setReplaceMethod("unions", "MassifquantParam", function(object, value) {
    object@unions <- as.integer(value)
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("checkBack", "MassifquantParam", function(object)
    return(object@checkBack))
#' @rdname hidden_aliases
setReplaceMethod("checkBack", "MassifquantParam", function(object, value) {
    object@checkBack <- as.integer(value)
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("withWave", "MassifquantParam", function(object)
    return(object@withWave))
#' @rdname hidden_aliases
setReplaceMethod("withWave", "MassifquantParam", function(object, value) {
    object@withWave <- value
    if (validObject(object))
        return(object)
})


############################################################
## MSWParam
###

#' @rdname hidden_aliases
setMethod("snthresh", "MSWParam", function(object){ return(object@snthresh)})
#' @rdname hidden_aliases
setReplaceMethod("snthresh", "MSWParam", function(object, value) {
    object@snthresh <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("verboseColumns", "MSWParam", function(object){
    return(object@verboseColumns)})
#' @rdname hidden_aliases
setReplaceMethod("verboseColumns", "MSWParam", function(object, value) {
    object@verboseColumns <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("scales", "MSWParam", function(object){ return(object@scales)})
#' @rdname hidden_aliases
setReplaceMethod("scales", "MSWParam", function(object, value) {
    object@scales <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("nearbyPeak", "MSWParam", function(object){ return(object@nearbyPeak)})
#' @rdname hidden_aliases
setReplaceMethod("nearbyPeak", "MSWParam", function(object, value) {
    object@nearbyPeak <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("peakScaleRange", "MSWParam", function(object){
    return(object@peakScaleRange)})
#' @rdname hidden_aliases
setReplaceMethod("peakScaleRange", "MSWParam", function(object, value) {
    object@peakScaleRange <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("ampTh", "MSWParam", function(object){ return(object@ampTh)})
#' @rdname hidden_aliases
setReplaceMethod("ampTh", "MSWParam", function(object, value) {
    object@ampTh <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("minNoiseLevel", "MSWParam", function(object){
    return(object@minNoiseLevel)})
#' @rdname hidden_aliases
setReplaceMethod("minNoiseLevel", "MSWParam", function(object, value) {
    object@minNoiseLevel <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("ridgeLength", "MSWParam", function(object){
    return(object@ridgeLength)})
#' @rdname hidden_aliases
setReplaceMethod("ridgeLength", "MSWParam", function(object, value) {
    object@ridgeLength <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("peakThr", "MSWParam", function(object){ return(object@peakThr)})
#' @rdname hidden_aliases
setReplaceMethod("peakThr", "MSWParam", function(object, value) {
    object@peakThr <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("tuneIn", "MSWParam", function(object){ return(object@tuneIn)})
#' @rdname hidden_aliases
setReplaceMethod("tuneIn", "MSWParam", function(object, value) {
    object@tuneIn <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("addParams", "MSWParam", function(object){ return(object@addParams)})
#' @rdname hidden_aliases
setReplaceMethod("addParams", "MSWParam", function(object, value) {
    object@addParams <- value
    if (validObject(object))
        return(object)
})
## The 'setAs' method.
setAs("MSWParam" ,"list", function(from){
    L <- .param2list(from)
    ## Rename ampTh to amp.Th
    names(L) <- sub(names(L), pattern = "ampTh", replacement = "amp.Th",
                    fixed = TRUE)
    L
})
setMethod("as.list", signature(x = "MSWParam"), function(x, ...) {
    L <- .param2list(x)
    ## Rename ampTh to amp.Th
    names(L) <- sub(names(L), pattern = "ampTh", replacement = "amp.Th",
                    fixed = TRUE)
    L
})

############################################################
## CentWavePredIsoParam
###

#' @rdname hidden_aliases
setMethod("snthreshIsoROIs", "CentWavePredIsoParam", function(object){
    return(object@snthreshIsoROIs)})
#' @rdname hidden_aliases
setReplaceMethod("snthreshIsoROIs", "CentWavePredIsoParam", function(object, value) {
    object@snthreshIsoROIs <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("maxCharge", "CentWavePredIsoParam", function(object){
    return(object@maxCharge)})
#' @rdname hidden_aliases
setReplaceMethod("maxCharge", "CentWavePredIsoParam", function(object, value) {
    object@maxCharge <- as.integer(value)
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("maxIso", "CentWavePredIsoParam", function(object){
    return(object@maxIso)})
#' @rdname hidden_aliases
setReplaceMethod("maxIso", "CentWavePredIsoParam", function(object, value) {
    object@maxIso <- as.integer(value)
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("mzIntervalExtension", "CentWavePredIsoParam", function(object){
    return(object@mzIntervalExtension)})
#' @rdname hidden_aliases
setReplaceMethod("mzIntervalExtension", "CentWavePredIsoParam",
                 function(object, value) {
                     object@mzIntervalExtension <- value
                     if (validObject(object))
                         return(object)
                 })

#' @rdname hidden_aliases
setMethod("polarity", "CentWavePredIsoParam", function(object){
    return(object@polarity)})
#' @rdname hidden_aliases
setReplaceMethod("polarity", "CentWavePredIsoParam", function(object, value) {
    object@polarity <- value
    if (validObject(object))
        return(object)
})


############################################################
## PeakDensityParam

#' @rdname hidden_aliases
setMethod("sampleGroups", "PeakDensityParam", function(object){
    return(object@sampleGroups)})
#' @rdname hidden_aliases
setReplaceMethod("sampleGroups", "PeakDensityParam", function(object, value) {
    if (length(value) == 0 | any(is.na(value)))
        stop("'sampleGroups' should not contain any NAs and its length has ",
             "to be >= 1")
    object@sampleGroups <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("bw", "PeakDensityParam", function(object){
    return(object@bw)})
#' @rdname hidden_aliases
setReplaceMethod("bw", "PeakDensityParam", function(object, value) {
    object@bw <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("minFraction", "PeakDensityParam", function(object){
    return(object@minFraction)})
#' @rdname hidden_aliases
setReplaceMethod("minFraction", "PeakDensityParam", function(object, value) {
    object@minFraction <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("minSamples", "PeakDensityParam", function(object){
    return(object@minSamples)})
#' @rdname hidden_aliases
setReplaceMethod("minSamples", "PeakDensityParam", function(object, value) {
    object@minSamples <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("binSize", "PeakDensityParam", function(object){
    return(object@binSize)})
#' @rdname hidden_aliases
setReplaceMethod("binSize", "PeakDensityParam", function(object, value) {
    object@binSize <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("maxFeatures", "PeakDensityParam", function(object){
    return(object@maxFeatures)})
#' @rdname hidden_aliases
setReplaceMethod("maxFeatures", "PeakDensityParam", function(object, value) {
    object@maxFeatures <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("ppm", "PeakDensityParam", function(object) {
    if (.hasSlot(object, "ppm"))
        object@ppm
    else 0.0
})


############################################################
## MzClustParam

#' @rdname hidden_aliases
setMethod("sampleGroups", "MzClustParam", function(object){
    return(object@sampleGroups)})
#' @rdname hidden_aliases
setReplaceMethod("sampleGroups", "MzClustParam", function(object, value) {
    object@sampleGroups <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("ppm", "MzClustParam", function(object){
    return(object@ppm)})
#' @rdname hidden_aliases
setReplaceMethod("ppm", "MzClustParam", function(object, value) {
    object@ppm <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("absMz", "MzClustParam", function(object){
    return(object@absMz)})
#' @rdname hidden_aliases
setReplaceMethod("absMz", "MzClustParam", function(object, value) {
    object@absMz <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("minFraction", "MzClustParam", function(object){
    return(object@minFraction)})
#' @rdname hidden_aliases
setReplaceMethod("minFraction", "MzClustParam", function(object, value) {
    object@minFraction <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("minSamples", "MzClustParam", function(object){
    return(object@minSamples)})
#' @rdname hidden_aliases
setReplaceMethod("minSamples", "MzClustParam", function(object, value) {
    object@minSamples <- value
    if (validObject(object))
        return(object)
})


############################################################
## NearestPeaksParam

#' @rdname hidden_aliases
setMethod("sampleGroups", "NearestPeaksParam", function(object){
    return(object@sampleGroups)})
#' @rdname hidden_aliases
setReplaceMethod("sampleGroups", "NearestPeaksParam", function(object, value) {
    object@sampleGroups <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("mzVsRtBalance", "NearestPeaksParam", function(object){
    return(object@mzVsRtBalance)})
#' @rdname hidden_aliases
setReplaceMethod("mzVsRtBalance", "NearestPeaksParam", function(object, value) {
    object@mzVsRtBalance <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("absMz", "NearestPeaksParam", function(object){
    return(object@absMz)})
#' @rdname hidden_aliases
setReplaceMethod("absMz", "NearestPeaksParam", function(object, value) {
    object@absMz <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("absRt", "NearestPeaksParam", function(object){
    return(object@absRt)})
#' @rdname hidden_aliases
setReplaceMethod("absRt", "NearestPeaksParam", function(object, value) {
    object@absRt <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("kNN", "NearestPeaksParam", function(object){
    return(object@kNN)})
#' @rdname hidden_aliases
setReplaceMethod("kNN", "NearestPeaksParam", function(object, value) {
    object@kNN <- value
    if (validObject(object))
        return(object)
})


############################################################
## PeakGroupsParam

#' @rdname hidden_aliases
setMethod("minFraction", "PeakGroupsParam", function(object){
    return(object@minFraction)})
#' @rdname hidden_aliases
setReplaceMethod("minFraction", "PeakGroupsParam", function(object, value) {
    object@minFraction <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("extraPeaks", "PeakGroupsParam", function(object){
    return(object@extraPeaks)})
#' @rdname hidden_aliases
setReplaceMethod("extraPeaks", "PeakGroupsParam", function(object, value) {
    object@extraPeaks <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("smooth", "PeakGroupsParam", function(x){
    return(x@smooth)})
#' @rdname hidden_aliases
setReplaceMethod("smooth", "PeakGroupsParam", function(object, value) {
    object@smooth <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("span", "PeakGroupsParam", function(object){
    return(object@span)})
#' @rdname hidden_aliases
setReplaceMethod("span", "PeakGroupsParam", function(object, value) {
    object@span <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("family", "PeakGroupsParam", function(object){
    return(object@family)})
#' @rdname hidden_aliases
setReplaceMethod("family", "PeakGroupsParam", function(object, value) {
    object@family <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("peakGroupsMatrix", "PeakGroupsParam", function(object){
    return(object@peakGroupsMatrix)})
#' @rdname hidden_aliases
setReplaceMethod("peakGroupsMatrix", "PeakGroupsParam", function(object, value) {
    object@peakGroupsMatrix <- value
    if (validObject(object))
        return(object)
})
#' @rdname hidden_aliases<
setMethod("subset", "PeakGroupsParam", function(x){
    return(x@subset)})
#' @rdname hidden_aliases
setReplaceMethod("subset", "PeakGroupsParam", function(object, value) {
    object@subset <- value
    if (validObject(object))
        return(object)
})
#' @rdname hidden_aliases
setMethod("subsetAdjust", "PeakGroupsParam", function(object){
    return(object@subsetAdjust)})
#' @rdname hidden_aliases
setReplaceMethod("subsetAdjust", "PeakGroupsParam", function(object, value) {
    object@subsetAdjust <- value
    if (validObject(object))
        return(object)
})

############################################################
## LamaParama

#' @rdname LamaParama
setMethod("plot", signature(x = "LamaParama"),
          function(x, index = 1L,
                   colPoints = "#00000060",
                   colFit = "#00000080",
                   xlab = "Matched Chromatographic peaks",
                   ylab = "Lamas",...){
    model <- .rt_model(method = x@method,
                              rt_map= x@rtMap[[index]], span = x@span,
                              resid_ratio = x@outlierTolerance,
                              zero_weight = x@zeroWeight,
                              bs = x@bs)
    datap <- x@rtMap[[index]]
    plot(datap[, 2L], datap[, 1L], type = "p", xlab = xlab, ylab = ylab,
         col = colPoints, ...)
    points(model, type = "l", col = colFit)
})


############################################################
## ObiwarpParam

#' @rdname hidden_aliases
setMethod("binSize", "ObiwarpParam", function(object){
    return(object@binSize)})
#' @rdname adjustRtime
setReplaceMethod("binSize", "ObiwarpParam", function(object, value) {
    object@binSize <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("centerSample", "ObiwarpParam", function(object){
    return(object@centerSample)})
#' @rdname hidden_aliases
setReplaceMethod("centerSample", "ObiwarpParam", function(object, value) {
    object@centerSample <- as.integer(value)
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("response", "ObiwarpParam", function(object){
    return(object@response)})
#' @rdname hidden_aliases
setReplaceMethod("response", "ObiwarpParam", function(object, value) {
    object@response <- as.integer(value)
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("distFun", "ObiwarpParam", function(object){
    return(object@distFun)})
#' @rdname hidden_aliases
setReplaceMethod("distFun", "ObiwarpParam", function(object, value) {
    object@distFun <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("gapInit", "ObiwarpParam", function(object){
    if (length(object@gapInit) == 0) {
        if (object@distFun == "cor" | object@distFun == "cor_opt")
            return(0.3)
        if (object@distFun == "cov" | object@distFun == "prd")
            return(0)
        if (object@distFun == "euc")
            return(0.9)
    }
    return(object@gapInit)})
#' @rdname hidden_aliases
setReplaceMethod("gapInit", "ObiwarpParam", function(object, value) {
    object@gapInit <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("gapExtend", "ObiwarpParam", function(object){
    if (length(object@gapExtend) == 0) {
        if (object@distFun == "cor" | object@distFun == "cor_opt")
            return(2.4)
        if (object@distFun == "cov")
            return(11.7)
        if (object@distFun == "euc")
            return(1.8)
        if (object@distFun == "prd")
            return(7.8)
    }
    return(object@gapExtend)})
#' @rdname hidden_aliases
setReplaceMethod("gapExtend", "ObiwarpParam", function(object, value) {
    object@gapExtend <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("factorDiag", "ObiwarpParam", function(object){
    return(object@factorDiag)})
#' @rdname hidden_aliases
setReplaceMethod("factorDiag", "ObiwarpParam", function(object, value) {
    object@factorDiag <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("factorGap", "ObiwarpParam", function(object){
    return(object@factorGap)})
#' @rdname hidden_aliases
setReplaceMethod("factorGap", "ObiwarpParam", function(object, value) {
    object@factorGap <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("localAlignment", "ObiwarpParam", function(object){
    return(object@localAlignment)})
#' @rdname hidden_aliases
setReplaceMethod("localAlignment", "ObiwarpParam", function(object, value) {
    object@localAlignment <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("initPenalty", "ObiwarpParam", function(object){
    return(object@initPenalty)})
#' @rdname hidden_aliases
setReplaceMethod("initPenalty", "ObiwarpParam", function(object, value) {
    object@initPenalty <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("subset", "ObiwarpParam", function(x){
    return(x@subset)})
#' @rdname hidden_aliases
setReplaceMethod("subset", "ObiwarpParam", function(object, value) {
    object@subset <- value
    if (validObject(object))
        return(object)
})
#' @rdname hidden_aliases
setMethod("subsetAdjust", "ObiwarpParam", function(object){
    return(object@subsetAdjust)})
#' @rdname hidden_aliases
setReplaceMethod("subsetAdjust", "ObiwarpParam", function(object, value) {
    object@subsetAdjust <- value
    if (validObject(object))
        return(object)
})

############################################################
## FillChromPeaksParam
###

#' @rdname hidden_aliases
setMethod("expandMz", "FillChromPeaksParam", function(object){
    return(object@expandMz)})
#' @rdname hidden_aliases
setReplaceMethod("expandMz", "FillChromPeaksParam", function(object, value) {
    object@expandMz <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("expandRt", "FillChromPeaksParam", function(object){
    return(object@expandRt)})
#' @rdname hidden_aliases
setReplaceMethod("expandRt", "FillChromPeaksParam", function(object, value) {
    object@expandRt <- value
    if (validObject(object))
        return(object)
})

#' @rdname hidden_aliases
setMethod("ppm", "FillChromPeaksParam", function(object){
    return(object@ppm)})
#' @rdname hidden_aliases
setReplaceMethod("ppm", "FillChromPeaksParam", function(object, value) {
    object@ppm <- value
    if (validObject(object))
        return(object)
})

#' @rdname findChromPeaks-centWave
setMethod("as.list", "CentWaveParam", function(x, ...) {
    x <- updateObject(x)
    callNextMethod(x)
})

#' @rdname groupChromPeaks
setMethod("as.list", "PeakDensityParam", function(x, ...) {
    x <- updateObject(x)
    callNextMethod(x)
})
