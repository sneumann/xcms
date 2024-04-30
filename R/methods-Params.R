## Methods for the Param class and sub-classes
#' @include functions-Params.R

#' @aliases ppm
#'
#' @description \code{ppm},\code{ppm<-}: getter and setter for the \code{ppm}
#'     slot of the object.
#'
#' @rdname findChromPeaks-centWave
setMethod("ppm", "CentWaveParam", function(object){ return(object@ppm)})
#' @aliases ppm<-
#'
#' @param value The value for the slot.
#'
#' @rdname findChromPeaks-centWave
setReplaceMethod("ppm", "CentWaveParam", function(object, value) {
    object@ppm <- value
    if (validObject(object))
        return(object)
})

#' @aliases peakwidth
#'
#' @description \code{peakwidth},\code{peakwidth<-}: getter and setter for the
#'     \code{peakwidth} slot of the object.
#'
#' @rdname findChromPeaks-centWave
setMethod("peakwidth", "CentWaveParam", function(object)
    return(object@peakwidth))
#' @aliases peakwidth<-
#'
#' @rdname findChromPeaks-centWave
setReplaceMethod("peakwidth", "CentWaveParam", function(object, value) {
    object@peakwidth <- value
    if (validObject(object))
        return(object)
})

#' @aliases snthresh
#'
#' @description \code{snthresh},\code{snthresh<-}: getter and setter for the
#'     \code{snthresh} slot of the object.
#'
#' @rdname findChromPeaks-centWave
setMethod("snthresh", "CentWaveParam", function(object)
    return(object@snthresh))
#' @aliases snthresh<-
#'
#' @rdname findChromPeaks-centWave
setReplaceMethod("snthresh", "CentWaveParam", function(object, value) {
    object@snthresh <- value
    if (validObject(object))
        return(object)
})

#' @aliases prefilter
#'
#' @description \code{prefilter},\code{prefilter<-}: getter and setter for the
#'     \code{prefilter} slot of the object.
#'
#' @rdname findChromPeaks-centWave
setMethod("prefilter", "CentWaveParam", function(object)
    return(object@prefilter))
#' @aliases prefilter<-
#'
#' @rdname findChromPeaks-centWave
setReplaceMethod("prefilter", "CentWaveParam", function(object, value) {
    object@prefilter <- value
    if (validObject(object))
        return(object)
})

#' @aliases mzCenterFun
#'
#' @description \code{mzCenterFun},\code{mzCenterFun<-}: getter and setter for the
#'     \code{mzCenterFun} slot of the object.
#'
#' @rdname findChromPeaks-centWave
setMethod("mzCenterFun", "CentWaveParam", function(object)
    return(object@mzCenterFun))
#' @aliases mzCenterFun<-
#'
#' @rdname findChromPeaks-centWave
setReplaceMethod("mzCenterFun", "CentWaveParam", function(object, value) {
    object@mzCenterFun <- value
    if (validObject(object))
        return(object)
})

#' @description \code{integrate},\code{integrate<-}: getter and setter for the
#'     \code{integrate} slot of the object.
#'
#' @param f For \code{integrate}: a \code{CentWaveParam} object.
#'
#' @rdname findChromPeaks-centWave
setMethod("integrate", signature(f = "CentWaveParam"), function(f)
    return(f@integrate))
#' @aliases integrate<-
#'
#' @rdname findChromPeaks-centWave
setReplaceMethod("integrate", "CentWaveParam", function(object, value) {
    object@integrate <- as.integer(value)
    if (validObject(object))
        return(object)
})

#' @aliases mzdiff
#'
#' @description \code{mzdiff},\code{mzdiff<-}: getter and setter for the
#'     \code{mzdiff} slot of the object.
#'
#' @rdname findChromPeaks-centWave
setMethod("mzdiff", "CentWaveParam", function(object)
    return(object@mzdiff))
#' @aliases mzdiff<-
#'
#' @rdname findChromPeaks-centWave
setReplaceMethod("mzdiff", "CentWaveParam", function(object, value) {
    object@mzdiff <- value
    if (validObject(object))
        return(object)
})

#' @aliases fitgauss
#'
#' @description \code{fitgauss},\code{fitgauss<-}: getter and setter for the
#'     \code{fitgauss} slot of the object.
#'
#' @rdname findChromPeaks-centWave
setMethod("fitgauss", "CentWaveParam", function(object)
    return(object@fitgauss))
#' @aliases fitgauss<-
#'
#' @rdname findChromPeaks-centWave
setReplaceMethod("fitgauss", "CentWaveParam", function(object, value) {
    object@fitgauss <- value
    if (validObject(object))
        return(object)
})

#' @aliases noise
#'
#' @description \code{noise},\code{noise<-}: getter and setter for the
#'     \code{noise} slot of the object.
#'
#' @rdname findChromPeaks-centWave
setMethod("noise", "CentWaveParam", function(object)
    return(object@noise))
#' @aliases noise<-
#'
#' @rdname findChromPeaks-centWave
setReplaceMethod("noise", "CentWaveParam", function(object, value) {
    object@noise <- value
    if (validObject(object))
        return(object)
})

#' @aliases verboseColumns
#'
#' @description \code{verboseColumns},\code{verboseColumns<-}: getter and
#'     setter for the \code{verboseColumns} slot of the object.
#'
#' @rdname findChromPeaks-centWave
setMethod("verboseColumns", "CentWaveParam", function(object)
    return(object@verboseColumns))
#' @aliases verboseColumns<-
#'
#' @rdname findChromPeaks-centWave
setReplaceMethod("verboseColumns", "CentWaveParam", function(object, value) {
    object@verboseColumns <- value
    if (validObject(object))
        return(object)
})

#' @aliases roiList
#'
#' @description \code{roiList},\code{roiList<-}: getter and setter for the
#'     \code{roiList} slot of the object.
#'
#' @rdname findChromPeaks-centWave
setMethod("roiList", "CentWaveParam", function(object)
    return(object@roiList))
#' @aliases roiList<-
#'
#' @rdname findChromPeaks-centWave
setReplaceMethod("roiList", "CentWaveParam", function(object, value) {
    object@roiList <- value
    if (validObject(object))
        return(object)
})

#' @aliases firstBaselineCheck
#'
#' @description \code{fistBaselineCheck},\code{firstBaselineCheck<-}: getter
#'     and setter for the \code{firstBaselineCheck} slot of the object.
#'
#' @rdname findChromPeaks-centWave
setMethod("firstBaselineCheck", "CentWaveParam", function(object)
    return(object@firstBaselineCheck))
#' @aliases firstBaselineCheck<-
#'
#' @rdname findChromPeaks-centWave
setReplaceMethod("firstBaselineCheck", "CentWaveParam", function(object, value) {
    object@firstBaselineCheck <- value
    if (validObject(object))
        return(object)
})

#' @aliases roiScales
#'
#' @description \code{roiScales},\code{roiScales<-}: getter and setter for the
#'     \code{roiScales} slot of the object.
#'
#' @rdname findChromPeaks-centWave
setMethod("roiScales", "CentWaveParam", function(object)
    return(object@roiScales))
#' @aliases roiScales<-
#'
#' @rdname findChromPeaks-centWave
setReplaceMethod("roiScales", "CentWaveParam", function(object, value) {
    object@roiScales <- value
    if (validObject(object))
        return(object)
})


############################################################
## MatchedFilterParam

#' @aliases binSize
#'
#' @description \code{binSize},\code{binSize<-}: getter and setter for the
#'     \code{binSize} slot of the object.
#'
#' @rdname findChromPeaks-matchedFilter
setMethod("binSize", "MatchedFilterParam", function(object)
    return(object@binSize))
#' @aliases binSize<-
#'
#' @param value The value for the slot.
#'
#' @rdname findChromPeaks-matchedFilter
setReplaceMethod("binSize", "MatchedFilterParam", function(object, value) {
    object@binSize <- value
    if (validObject(object))
        return(object)
})

#' @description \code{impute},\code{impute<-}: getter and setter for the
#'     \code{impute} slot of the object.
#'
#' @rdname findChromPeaks-matchedFilter
setMethod("impute", "MatchedFilterParam", function(object)
    return(object@impute))
#' @aliases impute<-
#'
#' @rdname findChromPeaks-matchedFilter
setReplaceMethod("impute", "MatchedFilterParam", function(object, value) {
    object@impute <- value
    if (validObject(object))
        return(object)
})

#' @aliases baseValue
#'
#' @description \code{baseValue},\code{baseValue<-}: getter and setter for the
#'     \code{baseValue} slot of the object.
#'
#' @rdname findChromPeaks-matchedFilter
setMethod("baseValue", "MatchedFilterParam", function(object)
    return(object@baseValue))
#' @aliases baseValue<-
#'
#' @rdname findChromPeaks-matchedFilter
setReplaceMethod("baseValue", "MatchedFilterParam", function(object, value) {
    object@baseValue <- value
    if (validObject(object))
        return(object)
})

#' @aliases distance
#'
#' @description \code{distance},\code{distance<-}: getter and setter for the
#'     \code{distance} slot of the object.
#'
#' @rdname findChromPeaks-matchedFilter
setMethod("distance", "MatchedFilterParam", function(object)
    return(object@distance))
#' @aliases distance<-
#'
#' @rdname findChromPeaks-matchedFilter
setReplaceMethod("distance", "MatchedFilterParam", function(object, value) {
    object@distance <- value
    if (validObject(object))
        return(object)
})

#' @aliases fwhm
#'
#' @description \code{fwhm},\code{fwhm<-}: getter and setter for the
#'     \code{fwhm} slot of the object.
#'
#' @rdname findChromPeaks-matchedFilter
setMethod("fwhm", "MatchedFilterParam", function(object)
    return(object@fwhm))
#' @aliases fwhm<-
#'
#' @rdname findChromPeaks-matchedFilter
setReplaceMethod("fwhm", "MatchedFilterParam", function(object, value) {
    object@fwhm <- value
    if (validObject(object))
        return(object)
})

#' @aliases sigma
#'
#' @description \code{sigma},\code{sigma<-}: getter and setter for the
#'     \code{sigma} slot of the object.
#'
#' @rdname findChromPeaks-matchedFilter
setMethod("sigma", "MatchedFilterParam", function(object)
    return(object@sigma))
#' @aliases sigma<-
#'
#' @rdname findChromPeaks-matchedFilter
setReplaceMethod("sigma", "MatchedFilterParam", function(object, value) {
    object@sigma <- value
    if (validObject(object))
        return(object)
})

#' @description \code{max},\code{max<-}: getter and setter for the
#'      \code{max} slot of the object.
#'
#' @param x For \code{max}: a \code{MatchedFilterParam} object.
#'
#' @rdname findChromPeaks-matchedFilter
setMethod("max", signature(x="MatchedFilterParam"),
          function(x) return(x@max))
#' @aliases max<-
#'
#' @rdname findChromPeaks-matchedFilter
setReplaceMethod("max", "MatchedFilterParam", function(object, value) {
    object@max <- value
    if (validObject(object))
        return(object)
})

#' @description \code{snthresh},\code{snthresh<-}: getter and setter for the
#'     \code{snthresh} slot of the object.
#'
#' @rdname findChromPeaks-matchedFilter
setMethod("snthresh", "MatchedFilterParam", function(object)
    return(object@snthresh))
#' @rdname findChromPeaks-matchedFilter
setReplaceMethod("snthresh", "MatchedFilterParam", function(object, value) {
    object@snthresh <- value
    if (validObject(object))
        return(object)
})

#' @aliases steps
#'
#' @description \code{steps},\code{steps<-}: getter and setter for the
#'     \code{steps} slot of the object.
#'
#' @rdname findChromPeaks-matchedFilter
setMethod("steps", "MatchedFilterParam", function(object)
    return(object@steps))
#' @aliases steps<-
#'
#' @rdname findChromPeaks-matchedFilter
setReplaceMethod("steps", "MatchedFilterParam", function(object, value) {
    object@steps <- value
    if (validObject(object))
        return(object)
})

#' @description \code{mzdiff},\code{mzdiff<-}: getter and setter for the
#'      \code{mzdiff} slot of the object.
#'
#' @rdname findChromPeaks-matchedFilter
setMethod("mzdiff", "MatchedFilterParam", function(object)
    return(object@mzdiff))
#' @rdname findChromPeaks-matchedFilter
setReplaceMethod("mzdiff", "MatchedFilterParam", function(object, value) {
    object@mzdiff <- value
    if (validObject(object))
        return(object)
})

#' @aliases index
#'
#' @description \code{index},\code{index<-}: getter and setter for the
#'     \code{index} slot of the object.
#'
#' @rdname findChromPeaks-matchedFilter
setMethod("index", "MatchedFilterParam", function(object)
    return(object@index))
#' @aliases index<-
#'
#' @rdname findChromPeaks-matchedFilter
setReplaceMethod("index", "MatchedFilterParam", function(object, value) {
    object@index <- value
    if (validObject(object))
        return(object)
})

############################################################
## MassifquantParam
###

#' @description \code{ppm},\code{ppm<-}: getter and setter for the \code{ppm}
#'     slot of the object.
#'
#' @rdname findChromPeaks-massifquant
setMethod("ppm", "MassifquantParam", function(object){ return(object@ppm)})
#' @param value The value for the slot.
#'
#' @rdname findChromPeaks-massifquant
setReplaceMethod("ppm", "MassifquantParam", function(object, value) {
    object@ppm <- value
    if (validObject(object))
        return(object)
})

#' @description \code{peakwidth},\code{peakwidth<-}: getter and setter for the
#'     \code{peakwidth} slot of the object.
#'
#' @rdname findChromPeaks-massifquant
setMethod("peakwidth", "MassifquantParam", function(object)
    return(object@peakwidth))
#' @rdname findChromPeaks-massifquant
setReplaceMethod("peakwidth", "MassifquantParam", function(object, value) {
    object@peakwidth <- value
    if (validObject(object))
        return(object)
})

#' @description \code{snthresh},\code{snthresh<-}: getter and setter for the
#'     \code{snthresh} slot of the object.
#'
#' @rdname findChromPeaks-massifquant
setMethod("snthresh", "MassifquantParam", function(object)
    return(object@snthresh))
#' @rdname findChromPeaks-massifquant
setReplaceMethod("snthresh", "MassifquantParam", function(object, value) {
    object@snthresh <- value
    if (validObject(object))
        return(object)
})

#' @description \code{prefilter},\code{prefilter<-}: getter and setter for the
#'     \code{prefilter} slot of the object.
#'
#' @rdname findChromPeaks-massifquant
setMethod("prefilter", "MassifquantParam", function(object)
    return(object@prefilter))
#' @rdname findChromPeaks-massifquant
setReplaceMethod("prefilter", "MassifquantParam", function(object, value) {
    object@prefilter <- value
    if (validObject(object))
        return(object)
})

#' @description \code{mzCenterFun},\code{mzCenterFun<-}: getter and setter for the
#'     \code{mzCenterFun} slot of the object.
#'
#' @rdname findChromPeaks-massifquant
setMethod("mzCenterFun", "MassifquantParam", function(object)
    return(object@mzCenterFun))
#' @rdname findChromPeaks-massifquant
setReplaceMethod("mzCenterFun", "MassifquantParam", function(object, value) {
    object@mzCenterFun <- value
    if (validObject(object))
        return(object)
})

#' @description \code{integrate},\code{integrate<-}: getter and setter for the
#'     \code{integrate} slot of the object.
#'
#' @param f For \code{integrate}: a \code{MassifquantParam} object.
#'
#' @rdname findChromPeaks-massifquant
setMethod("integrate", signature(f = "MassifquantParam"), function(f)
    return(f@integrate))
#' @rdname findChromPeaks-massifquant
setReplaceMethod("integrate", "MassifquantParam", function(object, value) {
    object@integrate <- as.integer(value)
    if (validObject(object))
        return(object)
})

#' @description \code{mzdiff},\code{mzdiff<-}: getter and setter for the
#'     \code{mzdiff} slot of the object.
#'
#' @rdname findChromPeaks-massifquant
setMethod("mzdiff", "MassifquantParam", function(object)
    return(object@mzdiff))
#' @rdname findChromPeaks-massifquant
setReplaceMethod("mzdiff", "MassifquantParam", function(object, value) {
    object@mzdiff <- value
    if (validObject(object))
        return(object)
})

#' @description \code{fitgauss},\code{fitgauss<-}: getter and setter for the
#'     \code{fitgauss} slot of the object.
#'
#' @rdname findChromPeaks-massifquant
setMethod("fitgauss", "MassifquantParam", function(object)
    return(object@fitgauss))
#' @rdname findChromPeaks-massifquant
setReplaceMethod("fitgauss", "MassifquantParam", function(object, value) {
    object@fitgauss <- value
    if (validObject(object))
        return(object)
})

#' @description \code{noise},\code{noise<-}: getter and setter for the
#'     \code{noise} slot of the object.
#'
#' @rdname findChromPeaks-massifquant
setMethod("noise", "MassifquantParam", function(object)
    return(object@noise))
#' @rdname findChromPeaks-massifquant
setReplaceMethod("noise", "MassifquantParam", function(object, value) {
    object@noise <- value
    if (validObject(object))
        return(object)
})

#' @description \code{verboseColumns},\code{verboseColumns<-}: getter and
#'     setter for the \code{verboseColumns} slot of the object.
#'
#' @rdname findChromPeaks-massifquant
setMethod("verboseColumns", "MassifquantParam", function(object)
    return(object@verboseColumns))
#' @rdname findChromPeaks-massifquant
setReplaceMethod("verboseColumns", "MassifquantParam", function(object, value) {
    object@verboseColumns <- value
    if (validObject(object))
        return(object)
})

#' @aliases criticalValue
#'
#' @description \code{criticalValue},\code{criticalValue<-}: getter and
#'     setter for the \code{criticalValue} slot of the object.
#'
#' @rdname findChromPeaks-massifquant
setMethod("criticalValue", "MassifquantParam", function(object)
    return(object@criticalValue))
#' @aliases criticalValue<-
#'
#' @rdname findChromPeaks-massifquant
setReplaceMethod("criticalValue", "MassifquantParam", function(object, value) {
    object@criticalValue <- value
    if (validObject(object))
        return(object)
})

#' @aliases consecMissedLimit
#'
#' @description \code{consecMissedLimit},\code{consecMissedLimit<-}: getter and
#'     setter for the \code{consecMissedLimit} slot of the object.
#'
#' @rdname findChromPeaks-massifquant
setMethod("consecMissedLimit", "MassifquantParam", function(object)
    return(object@consecMissedLimit))
#' @aliases consecMissedLimit<-
#'
#' @rdname findChromPeaks-massifquant
setReplaceMethod("consecMissedLimit", "MassifquantParam",
                 function(object, value) {
                     object@consecMissedLimit <- as.integer(value)
                     if (validObject(object))
                         return(object)
                 })

#' @aliases unions
#'
#' @description \code{unions},\code{unions<-}: getter and
#'     setter for the \code{unions} slot of the object.
#'
#' @rdname findChromPeaks-massifquant
setMethod("unions", "MassifquantParam", function(object)
    return(object@unions))
#' @aliases unions<-
#'
#' @rdname findChromPeaks-massifquant
setReplaceMethod("unions", "MassifquantParam", function(object, value) {
    object@unions <- as.integer(value)
    if (validObject(object))
        return(object)
})

#' @aliases checkBack
#'
#' @description \code{checkBack},\code{checkBack<-}: getter and
#'     setter for the \code{checkBack} slot of the object.
#'
#' @rdname findChromPeaks-massifquant
setMethod("checkBack", "MassifquantParam", function(object)
    return(object@checkBack))
#' @aliases checkBack<-
#'
#' @rdname findChromPeaks-massifquant
setReplaceMethod("checkBack", "MassifquantParam", function(object, value) {
    object@checkBack <- as.integer(value)
    if (validObject(object))
        return(object)
})

#' @aliases withWave
#'
#' @description \code{withWave},\code{withWave<-}: getter and
#'     setter for the \code{withWave} slot of the object.
#'
#' @rdname findChromPeaks-massifquant
setMethod("withWave", "MassifquantParam", function(object)
    return(object@withWave))
#' @aliases withWave<-
#'
#' @rdname findChromPeaks-massifquant
setReplaceMethod("withWave", "MassifquantParam", function(object, value) {
    object@withWave <- value
    if (validObject(object))
        return(object)
})


############################################################
## MSWParam
###

#' @description \code{snthresh},\code{snthresh<-}: getter and setter for the
#'     \code{snthresh} slot of the object.
#'
#' @rdname findPeaks-MSW
setMethod("snthresh", "MSWParam", function(object){ return(object@snthresh)})
#' @param value The value for the slot.
#'
#' @rdname findPeaks-MSW
setReplaceMethod("snthresh", "MSWParam", function(object, value) {
    object@snthresh <- value
    if (validObject(object))
        return(object)
})

#' @description \code{verboseColumns},\code{verboseColumns<-}: getter and setter
#'     for the \code{verboseColumns} slot of the object.
#'
#' @rdname findPeaks-MSW
setMethod("verboseColumns", "MSWParam", function(object){
    return(object@verboseColumns)})
#' @rdname findPeaks-MSW
setReplaceMethod("verboseColumns", "MSWParam", function(object, value) {
    object@verboseColumns <- value
    if (validObject(object))
        return(object)
})

#' @aliases scales
#'
#' @description \code{scales},\code{scales<-}: getter and setter for the
#'     \code{scales} slot of the object.
#'
#' @rdname findPeaks-MSW
setMethod("scales", "MSWParam", function(object){ return(object@scales)})
#' @aliases scales<-
#'
#' @rdname findPeaks-MSW
setReplaceMethod("scales", "MSWParam", function(object, value) {
    object@scales <- value
    if (validObject(object))
        return(object)
})

#' @aliases nearbyPeak
#'
#' @description \code{nearbyPeak},\code{nearbyPeak<-}: getter and setter for the
#'     \code{nearbyPeak} slot of the object.
#'
#' @rdname findPeaks-MSW
setMethod("nearbyPeak", "MSWParam", function(object){ return(object@nearbyPeak)})
#' @aliases nearbyPeak<-
#'
#' @rdname findPeaks-MSW
setReplaceMethod("nearbyPeak", "MSWParam", function(object, value) {
    object@nearbyPeak <- value
    if (validObject(object))
        return(object)
})

#' @aliases peakScaleRange
#'
#' @description \code{peakScaleRange},\code{peakScaleRange<-}: getter and setter
#'     for the \code{peakScaleRange} slot of the object.
#'
#' @rdname findPeaks-MSW
setMethod("peakScaleRange", "MSWParam", function(object){
    return(object@peakScaleRange)})
#' @aliases peakScaleRange<-
#'
#' @rdname findPeaks-MSW
setReplaceMethod("peakScaleRange", "MSWParam", function(object, value) {
    object@peakScaleRange <- value
    if (validObject(object))
        return(object)
})

#' @aliases ampTh
#'
#' @description \code{ampTh},\code{ampTh<-}: getter and setter for the
#'     \code{ampTh} slot of the object.
#'
#' @rdname findPeaks-MSW
setMethod("ampTh", "MSWParam", function(object){ return(object@ampTh)})
#' @aliases ampTh<-
#'
#' @rdname findPeaks-MSW
setReplaceMethod("ampTh", "MSWParam", function(object, value) {
    object@ampTh <- value
    if (validObject(object))
        return(object)
})

#' @aliases minNoiseLevel
#'
#' @description \code{minNoiseLevel},\code{minNoiseLevel<-}: getter and setter
#'     for the \code{minNoiseLevel} slot of the object.
#'
#' @rdname findPeaks-MSW
setMethod("minNoiseLevel", "MSWParam", function(object){
    return(object@minNoiseLevel)})
#' @aliases minNoiseLevel<-
#'
#' @rdname findPeaks-MSW
setReplaceMethod("minNoiseLevel", "MSWParam", function(object, value) {
    object@minNoiseLevel <- value
    if (validObject(object))
        return(object)
})

#' @aliases ridgeLength
#'
#' @description \code{ridgeLength},\code{ridgeLength<-}: getter and setter for
#'     the \code{ridgeLength} slot of the object.
#'
#' @rdname findPeaks-MSW
setMethod("ridgeLength", "MSWParam", function(object){
    return(object@ridgeLength)})
#' @aliases ridgeLength<-
#'
#' @rdname findPeaks-MSW
setReplaceMethod("ridgeLength", "MSWParam", function(object, value) {
    object@ridgeLength <- value
    if (validObject(object))
        return(object)
})

#' @aliases peakThr
#'
#' @description \code{peakThr},\code{peakThr<-}: getter and setter for the
#'     \code{peakThr} slot of the object.
#'
#' @rdname findPeaks-MSW
setMethod("peakThr", "MSWParam", function(object){ return(object@peakThr)})
#' @aliases peakThr<-
#'
#' @rdname findPeaks-MSW
setReplaceMethod("peakThr", "MSWParam", function(object, value) {
    object@peakThr <- value
    if (validObject(object))
        return(object)
})

#' @aliases tuneIn
#'
#' @description \code{tuneIn},\code{tuneIn<-}: getter and setter for the
#'     \code{tuneIn} slot of the object.
#'
#' @rdname findPeaks-MSW
setMethod("tuneIn", "MSWParam", function(object){ return(object@tuneIn)})
#' @aliases tuneIn<-
#'
#' @rdname findPeaks-MSW
setReplaceMethod("tuneIn", "MSWParam", function(object, value) {
    object@tuneIn <- value
    if (validObject(object))
        return(object)
})

#' @aliases addParams
#'
#' @description \code{addParams},\code{addParams<-}: getter and setter for the
#'     \code{addParams} slot of the object. This slot stores optional additional
#'     parameters to be passed to the
#'     \code{\link{identifyMajorPeaks}} and
#'     \code{\link{peakDetectionCWT}} functions from the
#'     \code{MassSpecWavelet} package.
#'
#' @rdname findPeaks-MSW
setMethod("addParams", "MSWParam", function(object){ return(object@addParams)})
#' @aliases addParams<-
#'
#' @rdname findPeaks-MSW
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

#' @aliases snthreshIsoROIs
#'
#' @description \code{snthreshIsoROIs},\code{snthreshIsoROIs<-}: getter and
#'     setter for the \code{snthreshIsoROIs} slot of the object.
#'
#' @rdname findChromPeaks-centWaveWithPredIsoROIs
setMethod("snthreshIsoROIs", "CentWavePredIsoParam", function(object){
    return(object@snthreshIsoROIs)})
#' @aliases snthreshIsoROIs<-
#'
#' @rdname findChromPeaks-centWaveWithPredIsoROIs
setReplaceMethod("snthreshIsoROIs", "CentWavePredIsoParam", function(object, value) {
    object@snthreshIsoROIs <- value
    if (validObject(object))
        return(object)
})

#' @aliases maxCharge
#'
#' @description \code{maxCharge},\code{maxCharge<-}: getter and
#'     setter for the \code{maxCharge} slot of the object.
#'
#' @rdname findChromPeaks-centWaveWithPredIsoROIs
setMethod("maxCharge", "CentWavePredIsoParam", function(object){
    return(object@maxCharge)})
#' @aliases maxCharge<-
#'
#' @rdname findChromPeaks-centWaveWithPredIsoROIs
setReplaceMethod("maxCharge", "CentWavePredIsoParam", function(object, value) {
    object@maxCharge <- as.integer(value)
    if (validObject(object))
        return(object)
})

#' @aliases maxIso
#'
#' @description \code{maxIso},\code{maxIso<-}: getter and
#'     setter for the \code{maxIso} slot of the object.
#'
#' @rdname findChromPeaks-centWaveWithPredIsoROIs
setMethod("maxIso", "CentWavePredIsoParam", function(object){
    return(object@maxIso)})
#' @aliases maxIso<-
#'
#' @rdname findChromPeaks-centWaveWithPredIsoROIs
setReplaceMethod("maxIso", "CentWavePredIsoParam", function(object, value) {
    object@maxIso <- as.integer(value)
    if (validObject(object))
        return(object)
})

#' @aliases mzIntervalExtension
#'
#' @description \code{mzIntervalExtension},\code{mzIntervalExtension<-}: getter
#'     and setter for the \code{mzIntervalExtension} slot of the object.
#'
#' @rdname findChromPeaks-centWaveWithPredIsoROIs
setMethod("mzIntervalExtension", "CentWavePredIsoParam", function(object){
    return(object@mzIntervalExtension)})
#' @aliases mzIntervalExtension<-
#'
#' @rdname findChromPeaks-centWaveWithPredIsoROIs
setReplaceMethod("mzIntervalExtension", "CentWavePredIsoParam",
                 function(object, value) {
                     object@mzIntervalExtension <- value
                     if (validObject(object))
                         return(object)
                 })

#' @description \code{polarity},\code{polarity<-}: getter and
#'     setter for the \code{polarity} slot of the object.
#'
#' @rdname findChromPeaks-centWaveWithPredIsoROIs
setMethod("polarity", "CentWavePredIsoParam", function(object){
    return(object@polarity)})
#' @aliases polarity<-
#'
#' @rdname findChromPeaks-centWaveWithPredIsoROIs
setReplaceMethod("polarity", "CentWavePredIsoParam", function(object, value) {
    object@polarity <- value
    if (validObject(object))
        return(object)
})


############################################################
## PeakDensityParam

#' @aliases sampleGroups
#'
#' @rdname groupChromPeaks
setMethod("sampleGroups", "PeakDensityParam", function(object){
    return(object@sampleGroups)})
#' @aliases sampleGroups<-
#'
#' @rdname groupChromPeaks
setReplaceMethod("sampleGroups", "PeakDensityParam", function(object, value) {
    if (length(value) == 0 | any(is.na(value)))
        stop("'sampleGroups' should not contain any NAs and its length has ",
             "to be >= 1")
    object@sampleGroups <- value
    if (validObject(object))
        return(object)
})

#' @aliases bw
#'
#' @rdname groupChromPeaks
setMethod("bw", "PeakDensityParam", function(object){
    return(object@bw)})
#' @aliases bw<-
#'
#' @rdname groupChromPeaks
setReplaceMethod("bw", "PeakDensityParam", function(object, value) {
    object@bw <- value
    if (validObject(object))
        return(object)
})

#' @aliases minFraction
#'
#' @rdname groupChromPeaks
setMethod("minFraction", "PeakDensityParam", function(object){
    return(object@minFraction)})
#' @aliases minFraction<-
#'
#' @rdname groupChromPeaks
setReplaceMethod("minFraction", "PeakDensityParam", function(object, value) {
    object@minFraction <- value
    if (validObject(object))
        return(object)
})

#' @aliases minSamples
#'
#' @rdname groupChromPeaks
setMethod("minSamples", "PeakDensityParam", function(object){
    return(object@minSamples)})
#' @aliases minSamples<-
#'
#' @rdname groupChromPeaks
setReplaceMethod("minSamples", "PeakDensityParam", function(object, value) {
    object@minSamples <- value
    if (validObject(object))
        return(object)
})

#' @rdname groupChromPeaks
setMethod("binSize", "PeakDensityParam", function(object){
    return(object@binSize)})
#' @rdname groupChromPeaks
setReplaceMethod("binSize", "PeakDensityParam", function(object, value) {
    object@binSize <- value
    if (validObject(object))
        return(object)
})

#' @aliases maxFeatures
#'
#' @rdname groupChromPeaks
setMethod("maxFeatures", "PeakDensityParam", function(object){
    return(object@maxFeatures)})
#' @aliases maxFeatures<-
#'
#' @rdname groupChromPeaks
setReplaceMethod("maxFeatures", "PeakDensityParam", function(object, value) {
    object@maxFeatures <- value
    if (validObject(object))
        return(object)
})

#' @rdname groupChromPeaks
setMethod("ppm", "PeakDensityParam", function(object) {
    if (.hasSlot(object, "ppm"))
        object@ppm
    else 0.0
})


############################################################
## MzClustParam

#' @rdname groupChromPeaks
setMethod("sampleGroups", "MzClustParam", function(object){
    return(object@sampleGroups)})
#' @param value The value for the slot.
#'
#' @rdname groupChromPeaks
setReplaceMethod("sampleGroups", "MzClustParam", function(object, value) {
    object@sampleGroups <- value
    if (validObject(object))
        return(object)
})

#' @rdname groupChromPeaks
setMethod("ppm", "MzClustParam", function(object){
    return(object@ppm)})
#' @rdname groupChromPeaks
setReplaceMethod("ppm", "MzClustParam", function(object, value) {
    object@ppm <- value
    if (validObject(object))
        return(object)
})

#' @aliases absMz
#'
#' @rdname groupChromPeaks
setMethod("absMz", "MzClustParam", function(object){
    return(object@absMz)})
#' @aliases absMz<-
#'
#' @rdname groupChromPeaks
setReplaceMethod("absMz", "MzClustParam", function(object, value) {
    object@absMz <- value
    if (validObject(object))
        return(object)
})

#' @rdname groupChromPeaks
setMethod("minFraction", "MzClustParam", function(object){
    return(object@minFraction)})
#' @rdname groupChromPeaks
setReplaceMethod("minFraction", "MzClustParam", function(object, value) {
    object@minFraction <- value
    if (validObject(object))
        return(object)
})

#' @rdname groupChromPeaks
setMethod("minSamples", "MzClustParam", function(object){
    return(object@minSamples)})
#' @rdname groupChromPeaks
setReplaceMethod("minSamples", "MzClustParam", function(object, value) {
    object@minSamples <- value
    if (validObject(object))
        return(object)
})


############################################################
## NearestPeaksParam

#' @rdname groupChromPeaks
setMethod("sampleGroups", "NearestPeaksParam", function(object){
    return(object@sampleGroups)})
#' @param value The value for the slot.
#'
#' @rdname groupChromPeaks
setReplaceMethod("sampleGroups", "NearestPeaksParam", function(object, value) {
    object@sampleGroups <- value
    if (validObject(object))
        return(object)
})

#' @aliases mzVsRtBalance
#'
#' @rdname groupChromPeaks
setMethod("mzVsRtBalance", "NearestPeaksParam", function(object){
    return(object@mzVsRtBalance)})
#' @aliases mzVsRtBalance<-
#'
#' @rdname groupChromPeaks
setReplaceMethod("mzVsRtBalance", "NearestPeaksParam", function(object, value) {
    object@mzVsRtBalance <- value
    if (validObject(object))
        return(object)
})

#' @rdname groupChromPeaks
setMethod("absMz", "NearestPeaksParam", function(object){
    return(object@absMz)})
#' @rdname groupChromPeaks
setReplaceMethod("absMz", "NearestPeaksParam", function(object, value) {
    object@absMz <- value
    if (validObject(object))
        return(object)
})

#' @aliases absRt
#'
#' @rdname groupChromPeaks
setMethod("absRt", "NearestPeaksParam", function(object){
    return(object@absRt)})
#' @aliases absRt<-
#'
#' @rdname groupChromPeaks
setReplaceMethod("absRt", "NearestPeaksParam", function(object, value) {
    object@absRt <- value
    if (validObject(object))
        return(object)
})

#' @aliases kNN
#'
#' @rdname groupChromPeaks
setMethod("kNN", "NearestPeaksParam", function(object){
    return(object@kNN)})
#' @aliases kNN<-
#'
#' @rdname groupChromPeaks
setReplaceMethod("kNN", "NearestPeaksParam", function(object, value) {
    object@kNN <- value
    if (validObject(object))
        return(object)
})


############################################################
## PeakGroupsParam

#' @rdname adjustRtime
setMethod("minFraction", "PeakGroupsParam", function(object){
    return(object@minFraction)})
#' @param value The value for the slot.
#'
#' @rdname adjustRtime
setReplaceMethod("minFraction", "PeakGroupsParam", function(object, value) {
    object@minFraction <- value
    if (validObject(object))
        return(object)
})

#' @aliases extraPeaks
#'
#' @rdname adjustRtime
setMethod("extraPeaks", "PeakGroupsParam", function(object){
    return(object@extraPeaks)})
#' @aliases extraPeaks<-
#'
#' @rdname adjustRtime
setReplaceMethod("extraPeaks", "PeakGroupsParam", function(object, value) {
    object@extraPeaks <- value
    if (validObject(object))
        return(object)
})

#' @aliases smooth
#'
#' @rdname adjustRtime
setMethod("smooth", "PeakGroupsParam", function(x){
    return(x@smooth)})
#' @aliases smooth<-
#'
#' @rdname adjustRtime
setReplaceMethod("smooth", "PeakGroupsParam", function(object, value) {
    object@smooth <- value
    if (validObject(object))
        return(object)
})

#' @aliases span
#'
#' @rdname adjustRtime
setMethod("span", "PeakGroupsParam", function(object){
    return(object@span)})
#' @aliases span<-
#'
#' @rdname adjustRtime
setReplaceMethod("span", "PeakGroupsParam", function(object, value) {
    object@span <- value
    if (validObject(object))
        return(object)
})

#' @aliases family
#'
#' @rdname adjustRtime
setMethod("family", "PeakGroupsParam", function(object){
    return(object@family)})
#' @aliases family<-
#'
#' @rdname adjustRtime
setReplaceMethod("family", "PeakGroupsParam", function(object, value) {
    object@family <- value
    if (validObject(object))
        return(object)
})

#' @aliases peakGroupsMatrix
#'
#' @rdname adjustRtime
setMethod("peakGroupsMatrix", "PeakGroupsParam", function(object){
    return(object@peakGroupsMatrix)})
#' @aliases peakGroupsMatrix<-
#'
#' @rdname adjustRtime
setReplaceMethod("peakGroupsMatrix", "PeakGroupsParam", function(object, value) {
    object@peakGroupsMatrix <- value
    if (validObject(object))
        return(object)
})
#' @aliases subset
#'
#' @rdname adjustRtime
setMethod("subset", "PeakGroupsParam", function(x){
    return(x@subset)})
#' @aliases subset<-
#'
#' @rdname adjustRtime
setReplaceMethod("subset", "PeakGroupsParam", function(object, value) {
    object@subset <- value
    if (validObject(object))
        return(object)
})
#' @aliases subsetAdjust
#'
#' @rdname adjustRtime
setMethod("subsetAdjust", "PeakGroupsParam", function(object){
    return(object@subsetAdjust)})
#' @aliases subsetAdjust<-
#'
#' @rdname adjustRtime
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

#' @rdname adjustRtime
setMethod("binSize", "ObiwarpParam", function(object){
    return(object@binSize)})
#' @rdname adjustRtime
setReplaceMethod("binSize", "ObiwarpParam", function(object, value) {
    object@binSize <- value
    if (validObject(object))
        return(object)
})

#' @aliases centerSample
#'
#' @rdname adjustRtime
setMethod("centerSample", "ObiwarpParam", function(object){
    return(object@centerSample)})
#' @aliases centerSample<-
#'
#' @rdname adjustRtime
setReplaceMethod("centerSample", "ObiwarpParam", function(object, value) {
    object@centerSample <- as.integer(value)
    if (validObject(object))
        return(object)
})

#' @aliases response
#'
#' @rdname adjustRtime
setMethod("response", "ObiwarpParam", function(object){
    return(object@response)})
#' @aliases response<-
#'
#' @rdname adjustRtime
setReplaceMethod("response", "ObiwarpParam", function(object, value) {
    object@response <- as.integer(value)
    if (validObject(object))
        return(object)
})

#' @aliases distFun
#'
#' @rdname adjustRtime
setMethod("distFun", "ObiwarpParam", function(object){
    return(object@distFun)})
#' @aliases distFun<-
#'
#' @rdname adjustRtime
setReplaceMethod("distFun", "ObiwarpParam", function(object, value) {
    object@distFun <- value
    if (validObject(object))
        return(object)
})

#' @aliases gapInit
#'
#' @rdname adjustRtime
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
#' @aliases gapInit<-
#'
#' @rdname adjustRtime
setReplaceMethod("gapInit", "ObiwarpParam", function(object, value) {
    object@gapInit <- value
    if (validObject(object))
        return(object)
})

#' @aliases gapExtend
#'
#' @rdname adjustRtime
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
#' @aliases gapExtend<-
#'
#' @rdname adjustRtime
setReplaceMethod("gapExtend", "ObiwarpParam", function(object, value) {
    object@gapExtend <- value
    if (validObject(object))
        return(object)
})

#' @aliases factorDiag
#'
#' @rdname adjustRtime
setMethod("factorDiag", "ObiwarpParam", function(object){
    return(object@factorDiag)})
#' @aliases factorDiag<-
#'
#' @rdname adjustRtime
setReplaceMethod("factorDiag", "ObiwarpParam", function(object, value) {
    object@factorDiag <- value
    if (validObject(object))
        return(object)
})

#' @aliases factorGap
#'
#' @rdname adjustRtime
setMethod("factorGap", "ObiwarpParam", function(object){
    return(object@factorGap)})
#' @aliases factorGap<-
#'
#' @rdname adjustRtime
setReplaceMethod("factorGap", "ObiwarpParam", function(object, value) {
    object@factorGap <- value
    if (validObject(object))
        return(object)
})

#' @aliases localAlignment
#'
#' @rdname adjustRtime
setMethod("localAlignment", "ObiwarpParam", function(object){
    return(object@localAlignment)})
#' @aliases localAlignment<-
#'
#' @rdname adjustRtime
setReplaceMethod("localAlignment", "ObiwarpParam", function(object, value) {
    object@localAlignment <- value
    if (validObject(object))
        return(object)
})

#' @aliases initPenalty
#'
#' @rdname adjustRtime
setMethod("initPenalty", "ObiwarpParam", function(object){
    return(object@initPenalty)})
#' @aliases initPenalty<-
#'
#' @rdname adjustRtime
setReplaceMethod("initPenalty", "ObiwarpParam", function(object, value) {
    object@initPenalty <- value
    if (validObject(object))
        return(object)
})

#' @rdname adjustRtime
setMethod("subset", "ObiwarpParam", function(x){
    return(x@subset)})
#' @rdname adjustRtime
setReplaceMethod("subset", "ObiwarpParam", function(object, value) {
    object@subset <- value
    if (validObject(object))
        return(object)
})
#' @rdname adjustRtime
setMethod("subsetAdjust", "ObiwarpParam", function(object){
    return(object@subsetAdjust)})
#' @rdname adjustRtime
setReplaceMethod("subsetAdjust", "ObiwarpParam", function(object, value) {
    object@subsetAdjust <- value
    if (validObject(object))
        return(object)
})

############################################################
## FillChromPeaksParam
###

#' @aliases expandMz
#'
#' @description \code{expandMz},\code{expandMz<-}: getter and setter
#'     for the \code{expandMz} slot of the object.
#'
#' @param value The value for the slot.
#'
#' @rdname fillChromPeaks
setMethod("expandMz", "FillChromPeaksParam", function(object){
    return(object@expandMz)})
#' @aliases expandMz<-
#'
#' @rdname fillChromPeaks
setReplaceMethod("expandMz", "FillChromPeaksParam", function(object, value) {
    object@expandMz <- value
    if (validObject(object))
        return(object)
})

#' @aliases expandRt
#'
#' @description \code{expandRt},\code{expandRt<-}: getter and setter
#'     for the \code{expandRt} slot of the object.
#'
#' @rdname fillChromPeaks
setMethod("expandRt", "FillChromPeaksParam", function(object){
    return(object@expandRt)})
#' @aliases expandRt<-
#'
#' @rdname fillChromPeaks
setReplaceMethod("expandRt", "FillChromPeaksParam", function(object, value) {
    object@expandRt <- value
    if (validObject(object))
        return(object)
})

#' @description \code{ppm},\code{ppm<-}: getter and setter
#'     for the \code{ppm} slot of the object.
#'
#' @rdname fillChromPeaks
setMethod("ppm", "FillChromPeaksParam", function(object){
    return(object@ppm)})
#' @rdname fillChromPeaks
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
