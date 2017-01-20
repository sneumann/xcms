## Methods for the Param class and sub-classes
#' @include functions-Params.R

############################################################
## Param
###
setMethod("as.list", signature(x = "Param"), function(x, ...) {
    return(.param2list(x))
})
## The 'setAs' method.
setAs("Param" ,"list", function(from){
    return(.param2list(from))
})
setMethod("initialize", "Param", function(.Object, ...) {
    classVersion(.Object)["Param"] <- "0.0.1"
    callNextMethod(.Object, ...)
})


############################################################
## CentWaveParam
###
setMethod("initialize", "CentWaveParam", function(.Object, ...) {
    classVersion(.Object)["CentWaveParam"] <- "0.0.1"
    callNextMethod(.Object, ...)
})

## ##' @rdname featureDetection-centWave
## setMethod("print", "CentWaveParam", function(x, ...) show(x))
##' @rdname featureDetection-centWave
setMethod("show", "CentWaveParam", function(object) {
    cat("Object of class: ", class(object), "\n")
    cat("Parameters:\n")
    cat(" ppm:", ppm(object), "\n")
    cat(" peakwidth:", paste(peakwidth(object), collapse = ", "), "\n")
    cat(" snthresh:", snthresh(object), "\n")
    cat(" prefilter:", paste(prefilter(object), collapse = ", "), "\n")
    cat(" mzCenterFun:", mzCenterFun(object), "\n")
    cat(" integrate:", integrate(object), "\n")
    cat(" mzdiff:", mzdiff(object), "\n")
    cat(" fitgauss:", fitgauss(object), "\n")
    cat(" noise:", noise(object), "\n")
    cat(" verboseColumns:", verboseColumns(object), "\n")
    cat(" roiList length:", length(roiList(object)), "\n")
    cat(" firstBaselineCheck", firstBaselineCheck(object), "\n")
    cat(" roiScales length:", length(roiScales(object)), "\n")
})

##' @aliases ppm
##' @description \code{ppm},\code{ppm<-}: getter and setter for the \code{ppm}
##' slot of the object.
##' @rdname featureDetection-centWave
setMethod("ppm", "CentWaveParam", function(object){ return(object@ppm)})
##' @aliases ppm<-
##' @param value The value for the slot.
##' @rdname featureDetection-centWave
setReplaceMethod("ppm", "CentWaveParam", function(object, value) {
    object@ppm <- value
    if (validObject(object))
        return(object)
})

##' @aliases peakwidth
##' @description \code{peakwidth},\code{peakwidth<-}: getter and setter for the
##' \code{peakwidth} slot of the object.
##' @rdname featureDetection-centWave
setMethod("peakwidth", "CentWaveParam", function(object)
    return(object@peakwidth))
##' @aliases peakwidth<-
##' @rdname featureDetection-centWave
setReplaceMethod("peakwidth", "CentWaveParam", function(object, value) {
    object@peakwidth <- value
    if (validObject(object))
        return(object)
})

##' @aliases snthresh
##' @description \code{snthresh},\code{snthresh<-}: getter and setter for the
##' \code{snthresh} slot of the object.
##' @rdname featureDetection-centWave
setMethod("snthresh", "CentWaveParam", function(object)
    return(object@snthresh))
##' @aliases snthresh<-
##' @rdname featureDetection-centWave
setReplaceMethod("snthresh", "CentWaveParam", function(object, value) {
    object@snthresh <- value
    if (validObject(object))
        return(object)
})

##' @aliases prefilter
##' @description \code{prefilter},\code{prefilter<-}: getter and setter for the
##' \code{prefilter} slot of the object.
##' @rdname featureDetection-centWave
setMethod("prefilter", "CentWaveParam", function(object)
    return(object@prefilter))
##' @aliases prefilter<-
##' @rdname featureDetection-centWave
setReplaceMethod("prefilter", "CentWaveParam", function(object, value) {
    object@prefilter <- value
    if (validObject(object))
        return(object)
})

##' @aliases mzCenterFun
##' @description \code{mzCenterFun},\code{mzCenterFun<-}: getter and setter for the
##' \code{mzCenterFun} slot of the object.
##' @rdname featureDetection-centWave
setMethod("mzCenterFun", "CentWaveParam", function(object)
    return(object@mzCenterFun))
##' @aliases mzCenterFun<-
##' @rdname featureDetection-centWave
setReplaceMethod("mzCenterFun", "CentWaveParam", function(object, value) {
    object@mzCenterFun <- value
    if (validObject(object))
        return(object)
})

##' @description \code{integrate},\code{integrate<-}: getter and setter for the
##' \code{integrate} slot of the object.
##' @param f For \code{integrate}: a \code{CentWaveParam} object.
##'
##' @rdname featureDetection-centWave
setMethod("integrate", signature(f = "CentWaveParam"), function(f)
    return(f@integrate))
##' @aliases integrate<-
##' @rdname featureDetection-centWave
setReplaceMethod("integrate", "CentWaveParam", function(object, value) {
    object@integrate <- as.integer(value)
    if (validObject(object))
        return(object)
})

##' @aliases mzdiff
##' @description \code{mzdiff},\code{mzdiff<-}: getter and setter for the
##' \code{mzdiff} slot of the object.
##' @rdname featureDetection-centWave
setMethod("mzdiff", "CentWaveParam", function(object)
    return(object@mzdiff))
##' @aliases mzdiff<-
##' @rdname featureDetection-centWave
setReplaceMethod("mzdiff", "CentWaveParam", function(object, value) {
    object@mzdiff <- value
    if (validObject(object))
        return(object)
})

##' @aliases fitgauss
##' @description \code{fitgauss},\code{fitgauss<-}: getter and setter for the
##' \code{fitgauss} slot of the object.
##' @rdname featureDetection-centWave
setMethod("fitgauss", "CentWaveParam", function(object)
    return(object@fitgauss))
##' @aliases fitgauss<-
##' @rdname featureDetection-centWave
setReplaceMethod("fitgauss", "CentWaveParam", function(object, value) {
    object@fitgauss <- value
    if (validObject(object))
        return(object)
})

##' @aliases noise
##' @description \code{noise},\code{noise<-}: getter and setter for the
##' \code{noise} slot of the object.
##' @rdname featureDetection-centWave
setMethod("noise", "CentWaveParam", function(object)
    return(object@noise))
##' @aliases noise<-
##' @rdname featureDetection-centWave
setReplaceMethod("noise", "CentWaveParam", function(object, value) {
    object@noise <- value
    if (validObject(object))
        return(object)
})

##' @aliases verboseColumns
##' @description \code{verboseColumns},\code{verboseColumns<-}: getter and
##' setter for the \code{verboseColumns} slot of the object.
##' @rdname featureDetection-centWave
setMethod("verboseColumns", "CentWaveParam", function(object)
    return(object@verboseColumns))
##' @aliases verboseColumns<-
##' @rdname featureDetection-centWave
setReplaceMethod("verboseColumns", "CentWaveParam", function(object, value) {
    object@verboseColumns <- value
    if (validObject(object))
        return(object)
})

##' @aliases roiList
##' @description \code{roiList},\code{roiList<-}: getter and setter for the
##' \code{roiList} slot of the object.
##' @rdname featureDetection-centWave
setMethod("roiList", "CentWaveParam", function(object)
    return(object@roiList))
##' @aliases roiList<-
##' @rdname featureDetection-centWave
setReplaceMethod("roiList", "CentWaveParam", function(object, value) {
    object@roiList <- value
    if (validObject(object))
        return(object)
})

##' @aliases firstBaselineCheck
##' @description \code{fistBaselineCheck},\code{firstBaselineCheck<-}: getter
##' and setter for the \code{firstBaselineCheck} slot of the object.
##' @rdname featureDetection-centWave
setMethod("firstBaselineCheck", "CentWaveParam", function(object)
    return(object@firstBaselineCheck))
##' @aliases firstBaselineCheck<-
##' @rdname featureDetection-centWave
setReplaceMethod("firstBaselineCheck", "CentWaveParam", function(object, value) {
    object@firstBaselineCheck <- value
    if (validObject(object))
        return(object)
})

##' @aliases roiScales
##' @description \code{roiScales},\code{roiScales<-}: getter and setter for the
##' \code{roiScales} slot of the object.
##' @rdname featureDetection-centWave
setMethod("roiScales", "CentWaveParam", function(object)
    return(object@roiScales))
##' @aliases roiScales<-
##' @rdname featureDetection-centWave
setReplaceMethod("roiScales", "CentWaveParam", function(object, value) {
    object@roiScales <- value
    if (validObject(object))
        return(object)
})


############################################################
## MatchedFilterParam
setMethod("initialize", "MatchedFilterParam", function(.Object, ...) {
    classVersion(.Object)["MatchedFilterParam"] <- "0.0.1"
    callNextMethod(.Object, ...)
})
##' @rdname featureDetection-matchedFilter
setMethod("show", "MatchedFilterParam", function(object) {
    cat("Object of class: ", class(object), "\n")
    cat("Parameters:\n")
    cat(" binSize:", binSize(object), "\n")
    cat(" impute:", impute(object), "\n")
    cat(" baseValue:", baseValue(object), "\n")
    cat(" distance:", distance(object), "\n")
    cat(" fwhm:", fwhm(object), "\n")
    cat(" sigma:", sigma(object), "\n")
    cat(" max:", max(object), "\n")
    cat(" snthresh:", snthresh(object), "\n")
    cat(" steps:", steps(object), "\n")
    cat(" mzdiff:", mzdiff(object), "\n")
    cat(" index:", index(object), "\n")
})

##' @aliases binSize
##' @description \code{binSize},\code{binSize<-}: getter and setter for the
##' \code{binSize} slot of the object.
##' @rdname featureDetection-matchedFilter
setMethod("binSize", "MatchedFilterParam", function(object)
    return(object@binSize))
##' @aliases binSize<-
##' @param value The value for the slot.
##' @rdname featureDetection-matchedFilter
setReplaceMethod("binSize", "MatchedFilterParam", function(object, value) {
    object@binSize <- value
    if (validObject(object))
        return(object)
})

##' @description \code{impute},\code{impute<-}: getter and setter for the
##' \code{impute} slot of the object.
##' @rdname featureDetection-matchedFilter
setMethod("impute", "MatchedFilterParam", function(object)
    return(object@impute))
##' @aliases impute<-
##' @rdname featureDetection-matchedFilter
setReplaceMethod("impute", "MatchedFilterParam", function(object, value) {
    object@impute <- value
    if (validObject(object))
        return(object)
})

##' @aliases baseValue
##' @description \code{baseValue},\code{baseValue<-}: getter and setter for the
##' \code{baseValue} slot of the object.
##' @rdname featureDetection-matchedFilter
setMethod("baseValue", "MatchedFilterParam", function(object)
    return(object@baseValue))
##' @aliases baseValue<-
##' @rdname featureDetection-matchedFilter
setReplaceMethod("baseValue", "MatchedFilterParam", function(object, value) {
    object@baseValue <- value
    if (validObject(object))
        return(object)
})

##' @aliases distance
##' @description \code{distance},\code{distance<-}: getter and setter for the
##' \code{distance} slot of the object.
##' @rdname featureDetection-matchedFilter
setMethod("distance", "MatchedFilterParam", function(object)
    return(object@distance))
##' @aliases distance<-
##' @rdname featureDetection-matchedFilter
setReplaceMethod("distance", "MatchedFilterParam", function(object, value) {
    object@distance <- value
    if (validObject(object))
        return(object)
})

##' @aliases fwhm
##' @description \code{fwhm},\code{fwhm<-}: getter and setter for the
##' \code{fwhm} slot of the object.
##' @rdname featureDetection-matchedFilter
setMethod("fwhm", "MatchedFilterParam", function(object)
    return(object@fwhm))
##' @aliases fwhm<-
##' @rdname featureDetection-matchedFilter
setReplaceMethod("fwhm", "MatchedFilterParam", function(object, value) {
    object@fwhm <- value
    if (validObject(object))
        return(object)
})

##' @aliases sigma
##' @description \code{sigma},\code{sigma<-}: getter and setter for the
##' \code{sigma} slot of the object.
##' @rdname featureDetection-matchedFilter
setMethod("sigma", "MatchedFilterParam", function(object)
    return(object@sigma))
##' @aliases sigma<-
##' @rdname featureDetection-matchedFilter
setReplaceMethod("sigma", "MatchedFilterParam", function(object, value) {
    object@sigma <- value
    if (validObject(object))
        return(object)
})

##' @description \code{max},\code{max<-}: getter and setter for the
##' \code{max} slot of the object.
##' @param x For \code{max}: a \code{MatchedFilterParam} object.
##' @rdname featureDetection-matchedFilter
setMethod("max", signature(x="MatchedFilterParam"),
          function(x) return(x@max))
##' @aliases max<-
##' @rdname featureDetection-matchedFilter
setReplaceMethod("max", "MatchedFilterParam", function(object, value) {
    object@max <- value
    if (validObject(object))
        return(object)
})

##' @description \code{snthresh},\code{snthresh<-}: getter and setter for the
##' \code{snthresh} slot of the object.
##' @rdname featureDetection-matchedFilter
setMethod("snthresh", "MatchedFilterParam", function(object)
    return(object@snthresh))
##' @rdname featureDetection-matchedFilter
setReplaceMethod("snthresh", "MatchedFilterParam", function(object, value) {
    object@snthresh <- value
    if (validObject(object))
        return(object)
})

##' @aliases steps
##' @description \code{steps},\code{steps<-}: getter and setter for the
##' \code{steps} slot of the object.
##' @rdname featureDetection-matchedFilter
setMethod("steps", "MatchedFilterParam", function(object)
    return(object@steps))
##' @aliases steps<-
##' @rdname featureDetection-matchedFilter
setReplaceMethod("steps", "MatchedFilterParam", function(object, value) {
    object@steps <- value
    if (validObject(object))
        return(object)
})

##' @description \code{mzdiff},\code{mzdiff<-}: getter and setter for the
##' \code{mzdiff} slot of the object.
##' @rdname featureDetection-matchedFilter
setMethod("mzdiff", "MatchedFilterParam", function(object)
    return(object@mzdiff))
##' @rdname featureDetection-matchedFilter
setReplaceMethod("mzdiff", "MatchedFilterParam", function(object, value) {
    object@mzdiff <- value
    if (validObject(object))
        return(object)
})

##' @aliases index
##' @description \code{index},\code{index<-}: getter and setter for the
##' \code{index} slot of the object.
##' @rdname featureDetection-matchedFilter
setMethod("index", "MatchedFilterParam", function(object)
    return(object@index))
##' @aliases index<-
##' @rdname featureDetection-matchedFilter
setReplaceMethod("index", "MatchedFilterParam", function(object, value) {
    object@index <- value
    if (validObject(object))
        return(object)
})

############################################################
## MassifquantParam
###
setMethod("initialize", "MassifquantParam", function(.Object, ...) {
    classVersion(.Object)["MassifquantParam"] <- "0.0.1"
    callNextMethod(.Object, ...)
})

##' @rdname featureDetection-massifquant
setMethod("show", "MassifquantParam", function(object) {
    cat("Object of class: ", class(object), "\n")
    cat("Parameters:\n")
    cat(" ppm:", ppm(object), "\n")
    cat(" peakwidth:", paste(peakwidth(object), collapse = ", "), "\n")
    cat(" snthresh:", snthresh(object), "\n")
    cat(" prefilter:", paste(prefilter(object), collapse = ", "), "\n")
    cat(" mzCenterFun:", mzCenterFun(object), "\n")
    cat(" integrate:", integrate(object), "\n")
    cat(" mzdiff:", mzdiff(object), "\n")
    cat(" fitgauss:", fitgauss(object), "\n")
    cat(" noise:", noise(object), "\n")
    cat(" verboseColumns:", verboseColumns(object), "\n")
    cat(" criticalValue:", criticalValue(object), "\n")
    cat(" consecMissedLimit:", consecMissedLimit(object), "\n")
    cat(" unions:", unions(object), "\n")
    cat(" checkBack:", checkBack(object), "\n")
    cat(" withWave:", withWave(object), "\n")
})

##' @description \code{ppm},\code{ppm<-}: getter and setter for the \code{ppm}
##' slot of the object.
##' @rdname featureDetection-massifquant
setMethod("ppm", "MassifquantParam", function(object){ return(object@ppm)})
##' @param value The value for the slot.
##' @rdname featureDetection-massifquant
setReplaceMethod("ppm", "MassifquantParam", function(object, value) {
    object@ppm <- value
    if (validObject(object))
        return(object)
})

##' @description \code{peakwidth},\code{peakwidth<-}: getter and setter for the
##' \code{peakwidth} slot of the object.
##' @rdname featureDetection-massifquant
setMethod("peakwidth", "MassifquantParam", function(object)
    return(object@peakwidth))
##' @rdname featureDetection-massifquant
setReplaceMethod("peakwidth", "MassifquantParam", function(object, value) {
    object@peakwidth <- value
    if (validObject(object))
        return(object)
})

##' @description \code{snthresh},\code{snthresh<-}: getter and setter for the
##' \code{snthresh} slot of the object.
##' @rdname featureDetection-massifquant
setMethod("snthresh", "MassifquantParam", function(object)
    return(object@snthresh))
##' @rdname featureDetection-massifquant
setReplaceMethod("snthresh", "MassifquantParam", function(object, value) {
    object@snthresh <- value
    if (validObject(object))
        return(object)
})

##' @description \code{prefilter},\code{prefilter<-}: getter and setter for the
##' \code{prefilter} slot of the object.
##' @rdname featureDetection-massifquant
setMethod("prefilter", "MassifquantParam", function(object)
    return(object@prefilter))
##' @rdname featureDetection-massifquant
setReplaceMethod("prefilter", "MassifquantParam", function(object, value) {
    object@prefilter <- value
    if (validObject(object))
        return(object)
})

##' @description \code{mzCenterFun},\code{mzCenterFun<-}: getter and setter for the
##' \code{mzCenterFun} slot of the object.
##' @rdname featureDetection-massifquant
setMethod("mzCenterFun", "MassifquantParam", function(object)
    return(object@mzCenterFun))
##' @rdname featureDetection-massifquant
setReplaceMethod("mzCenterFun", "MassifquantParam", function(object, value) {
    object@mzCenterFun <- value
    if (validObject(object))
        return(object)
})

##' @description \code{integrate},\code{integrate<-}: getter and setter for the
##' \code{integrate} slot of the object.
##' @param f For \code{integrate}: a \code{MassifquantParam} object.
##'
##' @rdname featureDetection-massifquant
setMethod("integrate", signature(f = "MassifquantParam"), function(f)
    return(f@integrate))
##' @rdname featureDetection-massifquant
setReplaceMethod("integrate", "MassifquantParam", function(object, value) {
    object@integrate <- as.integer(value)
    if (validObject(object))
        return(object)
})

##' @description \code{mzdiff},\code{mzdiff<-}: getter and setter for the
##' \code{mzdiff} slot of the object.
##' @rdname featureDetection-massifquant
setMethod("mzdiff", "MassifquantParam", function(object)
    return(object@mzdiff))
##' @rdname featureDetection-massifquant
setReplaceMethod("mzdiff", "MassifquantParam", function(object, value) {
    object@mzdiff <- value
    if (validObject(object))
        return(object)
})

##' @description \code{fitgauss},\code{fitgauss<-}: getter and setter for the
##' \code{fitgauss} slot of the object.
##' @rdname featureDetection-massifquant
setMethod("fitgauss", "MassifquantParam", function(object)
    return(object@fitgauss))
##' @rdname featureDetection-massifquant
setReplaceMethod("fitgauss", "MassifquantParam", function(object, value) {
    object@fitgauss <- value
    if (validObject(object))
        return(object)
})

##' @description \code{noise},\code{noise<-}: getter and setter for the
##' \code{noise} slot of the object.
##' @rdname featureDetection-massifquant
setMethod("noise", "MassifquantParam", function(object)
    return(object@noise))
##' @rdname featureDetection-massifquant
setReplaceMethod("noise", "MassifquantParam", function(object, value) {
    object@noise <- value
    if (validObject(object))
        return(object)
})

##' @description \code{verboseColumns},\code{verboseColumns<-}: getter and
##' setter for the \code{verboseColumns} slot of the object.
##' @rdname featureDetection-massifquant
setMethod("verboseColumns", "MassifquantParam", function(object)
    return(object@verboseColumns))
##' @rdname featureDetection-massifquant
setReplaceMethod("verboseColumns", "MassifquantParam", function(object, value) {
    object@verboseColumns <- value
    if (validObject(object))
        return(object)
})

##' @aliases criticalValue
##' @description \code{criticalValue},\code{criticalValue<-}: getter and
##' setter for the \code{criticalValue} slot of the object.
##' @rdname featureDetection-massifquant
setMethod("criticalValue", "MassifquantParam", function(object)
    return(object@criticalValue))
##' @aliases criticalValue<-
##' @rdname featureDetection-massifquant
setReplaceMethod("criticalValue", "MassifquantParam", function(object, value) {
    object@criticalValue <- value
    if (validObject(object))
        return(object)
})

##' @aliases consecMissedLimit
##' @description \code{consecMissedLimit},\code{consecMissedLimit<-}: getter and
##' setter for the \code{consecMissedLimit} slot of the object.
##' @rdname featureDetection-massifquant
setMethod("consecMissedLimit", "MassifquantParam", function(object)
    return(object@consecMissedLimit))
##' @aliases consecMissedLimit<-
##' @rdname featureDetection-massifquant
setReplaceMethod("consecMissedLimit", "MassifquantParam",
                 function(object, value) {
                     object@consecMissedLimit <- as.integer(value)
                     if (validObject(object))
                         return(object)
                 })

##' @aliases unions
##' @description \code{unions},\code{unions<-}: getter and
##' setter for the \code{unions} slot of the object.
##' @rdname featureDetection-massifquant
setMethod("unions", "MassifquantParam", function(object)
    return(object@unions))
##' @aliases unions<-
##' @rdname featureDetection-massifquant
setReplaceMethod("unions", "MassifquantParam", function(object, value) {
    object@unions <- as.integer(value)
    if (validObject(object))
        return(object)
})

##' @aliases checkBack
##' @description \code{checkBack},\code{checkBack<-}: getter and
##' setter for the \code{checkBack} slot of the object.
##' @rdname featureDetection-massifquant
setMethod("checkBack", "MassifquantParam", function(object)
    return(object@checkBack))
##' @aliases checkBack<-
##' @rdname featureDetection-massifquant
setReplaceMethod("checkBack", "MassifquantParam", function(object, value) {
    object@checkBack <- as.integer(value)
    if (validObject(object))
        return(object)
})

##' @aliases withWave
##' @description \code{withWave},\code{withWave<-}: getter and
##' setter for the \code{withWave} slot of the object.
##' @rdname featureDetection-massifquant
setMethod("withWave", "MassifquantParam", function(object)
    return(object@withWave))
##' @aliases withWave<-
##' @rdname featureDetection-massifquant
setReplaceMethod("withWave", "MassifquantParam", function(object, value) {
    object@withWave <- value
    if (validObject(object))
        return(object)
})


############################################################
## MSWParam
###
setMethod("initialize", "MSWParam", function(.Object, ...) {
    classVersion(.Object)["MSWParam"] <- "0.0.1"
    callNextMethod(.Object, ...)
})

##' @rdname featureDetection-MSW
setMethod("show", "MSWParam", function(object) {
    cat("Object of class: ", class(object), "\n")
    cat("Parameters:\n")
    cat(" snthresh:", snthresh(object), "\n")
    cat(" verboseColumns:", verboseColumns(object), "\n")
    cat(" scales:", paste(scales(object), collapse = ","), "\n")
    cat(" nearbyPeak:", nearbyPeak(object), "\n")
    cat(" peakScaleRange:", peakScaleRange(object), "\n")
    cat(" ampTh:", ampTh(object), "\n")
    cat(" minNoiseLevel:", minNoiseLevel(object), "\n")
    cat(" ridgeLength:", ridgeLength(object), "\n")
    cat(" peakThr:", peakThr(object), "\n")
    cat(" tuneIn:", tuneIn(object), "\n")
    parms <- addParams(object)
    if (length(parms) > 0) {
        cat(" additional parameters:\n")
        for (i in 1:length(parms)) {
            cat("  ", names(parms)[i], ": ", parms[[i]], "\n", sep = "")
        }
    }
})

##' @description \code{snthresh},\code{snthresh<-}: getter and setter for the
##' \code{snthresh} slot of the object.
##' @rdname featureDetection-MSW
setMethod("snthresh", "MSWParam", function(object){ return(object@snthresh)})
##' @param value The value for the slot.
##' @rdname featureDetection-MSW
setReplaceMethod("snthresh", "MSWParam", function(object, value) {
    object@snthresh <- value
    if (validObject(object))
        return(object)
})

##' @description \code{verboseColumns},\code{verboseColumns<-}: getter and setter
##' for the \code{verboseColumns} slot of the object.
##' @rdname featureDetection-MSW
setMethod("verboseColumns", "MSWParam", function(object){
    return(object@verboseColumns)})
##' @rdname featureDetection-MSW
setReplaceMethod("verboseColumns", "MSWParam", function(object, value) {
    object@verboseColumns <- value
    if (validObject(object))
        return(object)
})

##' @aliases scales
##' @description \code{scales},\code{scales<-}: getter and setter for the
##' \code{scales} slot of the object.
##' @rdname featureDetection-MSW
setMethod("scales", "MSWParam", function(object){ return(object@scales)})
##' @aliases scales<-
##' @rdname featureDetection-MSW
setReplaceMethod("scales", "MSWParam", function(object, value) {
    object@scales <- value
    if (validObject(object))
        return(object)
})

##' @aliases nearbyPeak
##' @description \code{nearbyPeak},\code{nearbyPeak<-}: getter and setter for the
##' \code{nearbyPeak} slot of the object.
##' @rdname featureDetection-MSW
setMethod("nearbyPeak", "MSWParam", function(object){ return(object@nearbyPeak)})
##' @aliases nearbyPeak<-
##' @rdname featureDetection-MSW
setReplaceMethod("nearbyPeak", "MSWParam", function(object, value) {
    object@nearbyPeak <- value
    if (validObject(object))
        return(object)
})

##' @aliases peakScaleRange
##' @description \code{peakScaleRange},\code{peakScaleRange<-}: getter and setter
##' for the \code{peakScaleRange} slot of the object.
##' @rdname featureDetection-MSW
setMethod("peakScaleRange", "MSWParam", function(object){
    return(object@peakScaleRange)})
##' @aliases peakScaleRange<-
##' @rdname featureDetection-MSW
setReplaceMethod("peakScaleRange", "MSWParam", function(object, value) {
    object@peakScaleRange <- value
    if (validObject(object))
        return(object)
})

##' @aliases ampTh
##' @description \code{ampTh},\code{ampTh<-}: getter and setter for the
##' \code{ampTh} slot of the object.
##' @rdname featureDetection-MSW
setMethod("ampTh", "MSWParam", function(object){ return(object@ampTh)})
##' @aliases ampTh<-
##' @rdname featureDetection-MSW
setReplaceMethod("ampTh", "MSWParam", function(object, value) {
    object@ampTh <- value
    if (validObject(object))
        return(object)
})

##' @aliases minNoiseLevel
##' @description \code{minNoiseLevel},\code{minNoiseLevel<-}: getter and setter
##' for the \code{minNoiseLevel} slot of the object.
##' @rdname featureDetection-MSW
setMethod("minNoiseLevel", "MSWParam", function(object){
    return(object@minNoiseLevel)})
##' @aliases minNoiseLevel<-
##' @rdname featureDetection-MSW
setReplaceMethod("minNoiseLevel", "MSWParam", function(object, value) {
    object@minNoiseLevel <- value
    if (validObject(object))
        return(object)
})

##' @aliases ridgeLength
##' @description \code{ridgeLength},\code{ridgeLength<-}: getter and setter for
##' the \code{ridgeLength} slot of the object.
##' @rdname featureDetection-MSW
setMethod("ridgeLength", "MSWParam", function(object){
    return(object@ridgeLength)})
##' @aliases ridgeLength<-
##' @rdname featureDetection-MSW
setReplaceMethod("ridgeLength", "MSWParam", function(object, value) {
    object@ridgeLength <- value
    if (validObject(object))
        return(object)
})

##' @aliases peakThr
##' @description \code{peakThr},\code{peakThr<-}: getter and setter for the
##' \code{peakThr} slot of the object.
##' @rdname featureDetection-MSW
setMethod("peakThr", "MSWParam", function(object){ return(object@peakThr)})
##' @aliases peakThr<-
##' @rdname featureDetection-MSW
setReplaceMethod("peakThr", "MSWParam", function(object, value) {
    object@peakThr <- value
    if (validObject(object))
        return(object)
})

##' @aliases tuneIn
##' @description \code{tuneIn},\code{tuneIn<-}: getter and setter for the
##' \code{tuneIn} slot of the object.
##' @rdname featureDetection-MSW
setMethod("tuneIn", "MSWParam", function(object){ return(object@tuneIn)})
##' @aliases tuneIn<-
##' @rdname featureDetection-MSW
setReplaceMethod("tuneIn", "MSWParam", function(object, value) {
    object@tuneIn <- value
    if (validObject(object))
        return(object)
})

##' @aliases addParams
##' @description \code{addParams},\code{addParams<-}: getter and setter for the
##' \code{addParams} slot of the object. This slot stores optional additional
##' parameters to be passed to the
##' \code{\link[MassSpecWavelet]{identifyMajorPeaks}} and
##' \code{\link[MassSpecWavelet]{sav.gol}} functions from the
##' \code{MassSpecWavelet} package.
##'
##' @rdname featureDetection-MSW
setMethod("addParams", "MSWParam", function(object){ return(object@addParams)})
##' @aliases addParams<-
##' @rdname featureDetection-MSW
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
    return(L)
})
setMethod("as.list", signature(x = "MSWParam"), function(x, ...) {
    L <- .param2list(x)
    ## Rename ampTh to amp.Th
    names(L) <- sub(names(L), pattern = "ampTh", replacement = "amp.Th",
                    fixed = TRUE)
    return(L)
})

############################################################
## CentWavePredIsoParam
###
setMethod("initialize", "CentWavePredIsoParam", function(.Object, ...) {
    classVersion(.Object)["CentWavePredIsoParam"] <- "0.0.1"
    callNextMethod(.Object, ...)
})

##' @rdname featureDetection-centWaveWithPredIsoROIs
setMethod("show", "CentWavePredIsoParam", function(object) {
    cat("Object of class: ", class(object), "\n")
    cat("Parameters:\n")
    cat(" ppm:", ppm(object), "\n")
    cat(" peakwidth:", paste(peakwidth(object), collapse = ", "), "\n")
    cat(" snthresh:", snthresh(object), "\n")
    cat(" prefilter:", paste(prefilter(object), collapse = ", "), "\n")
    cat(" mzCenterFun:", mzCenterFun(object), "\n")
    cat(" integrate:", integrate(object), "\n")
    cat(" mzdiff:", mzdiff(object), "\n")
    cat(" fitgauss:", fitgauss(object), "\n")
    cat(" noise:", noise(object), "\n")
    cat(" verboseColumns:", verboseColumns(object), "\n")
    cat(" roiList length:", length(roiList(object)), "\n")
    cat(" firstBaselineCheck", firstBaselineCheck(object), "\n")
    cat(" roiScales length:", length(roiScales(object)), "\n")
    cat(" snthreshIsoROIs:", snthreshIsoROIs(object), "\n")
    cat(" maxCharge:", maxCharge(object), "\n")
    cat(" maxIso:", maxIso(object), "\n")
    cat(" mzIntervalExtension:", mzIntervalExtension(object), "\n")
    cat(" polarity:", polarity(object), "\n")
})

##' @aliases snthreshIsoROIs
##' @description \code{snthreshIsoROIs},\code{snthreshIsoROIs<-}: getter and
##' setter for the \code{snthreshIsoROIs} slot of the object.
##' @rdname featureDetection-centWaveWithPredIsoROIs
setMethod("snthreshIsoROIs", "CentWavePredIsoParam", function(object){
    return(object@snthreshIsoROIs)})
##' @aliases snthreshIsoROIs<-
##' @rdname featureDetection-centWaveWithPredIsoROIs
setReplaceMethod("snthreshIsoROIs", "CentWavePredIsoParam", function(object, value) {
    object@snthreshIsoROIs <- value
    if (validObject(object))
        return(object)
})

##' @aliases maxCharge
##' @description \code{maxCharge},\code{maxCharge<-}: getter and
##' setter for the \code{maxCharge} slot of the object.
##' @rdname featureDetection-centWaveWithPredIsoROIs
setMethod("maxCharge", "CentWavePredIsoParam", function(object){
    return(object@maxCharge)})
##' @aliases maxCharge<-
##' @rdname featureDetection-centWaveWithPredIsoROIs
setReplaceMethod("maxCharge", "CentWavePredIsoParam", function(object, value) {
    object@maxCharge <- as.integer(value)
    if (validObject(object))
        return(object)
})

##' @aliases maxIso
##' @description \code{maxIso},\code{maxIso<-}: getter and
##' setter for the \code{maxIso} slot of the object.
##' @rdname featureDetection-centWaveWithPredIsoROIs
setMethod("maxIso", "CentWavePredIsoParam", function(object){
    return(object@maxIso)})
##' @aliases maxIso<-
##' @rdname featureDetection-centWaveWithPredIsoROIs
setReplaceMethod("maxIso", "CentWavePredIsoParam", function(object, value) {
    object@maxIso <- as.integer(value)
    if (validObject(object))
        return(object)
})

##' @aliases mzIntervalExtension
##' @description \code{mzIntervalExtension},\code{mzIntervalExtension<-}: getter
##' and setter for the \code{mzIntervalExtension} slot of the object.
##' @rdname featureDetection-centWaveWithPredIsoROIs
setMethod("mzIntervalExtension", "CentWavePredIsoParam", function(object){
    return(object@mzIntervalExtension)})
##' @aliases mzIntervalExtension<-
##' @rdname featureDetection-centWaveWithPredIsoROIs
setReplaceMethod("mzIntervalExtension", "CentWavePredIsoParam",
                 function(object, value) {
                     object@mzIntervalExtension <- value
                     if (validObject(object))
                         return(object)
                 })

##' @description \code{polarity},\code{polarity<-}: getter and
##' setter for the \code{polarity} slot of the object.
##' @rdname featureDetection-centWaveWithPredIsoROIs
setMethod("polarity", "CentWavePredIsoParam", function(object){
    return(object@polarity)})
##' @aliases polarity<-
##' @rdname featureDetection-centWaveWithPredIsoROIs
setReplaceMethod("polarity", "CentWavePredIsoParam", function(object, value) {
    object@polarity <- value
    if (validObject(object))
        return(object)
})


############################################################
## FeatureDensityParam
setMethod("initialize", "FeatureDensityParam", function(.Object, ...) {
    classVersion(.Object)["FeatureDensityParam"] <- "0.0.1"
    callNextMethod(.Object, ...)
})

## ##' @rdname groupFeatures-density
## setMethod("print", "CentWaveParam", function(x, ...) show(x))
##' @rdname groupFeatures-density
setMethod("show", "FeatureDensityParam", function(object) {
    cat("Object of class: ", class(object), "\n")
    cat("Parameters:\n")
    cat(" sampleGroups:", class(object@sampleGroups), "of length",
        length(object@sampleGroups), "\n")
    cat(" bw:", object@bw, "\n")
    cat(" minFraction:", minFraction(object), "\n")
    cat(" minSamples:", minSamples(object), "\n")
    cat(" binSize:", binSize(object), "\n")
    cat(" maxFeatures:", maxFeatures(object), "\n")
})

##' @aliases sampleGroups
##' @description \code{sampleGroups},\code{sampleGroups<-}: getter and setter
##' for the \code{sampleGroups} slot of the object.
##' @rdname groupFeatures-density
setMethod("sampleGroups", "FeatureDensityParam", function(object){
    return(object@sampleGroups)})
##' @aliases sampleGroups<-
##' @param value The value for the slot.
##' @rdname groupFeatures-density
setReplaceMethod("sampleGroups", "FeatureDensityParam", function(object, value) {
    object@sampleGroups <- value
    if (validObject(object))
        return(object)
})

##' @aliases bw
##' @description \code{bw},\code{bw<-}: getter and setter for the \code{bw} slot
##' of the object.
##' @rdname groupFeatures-density
setMethod("bw", "FeatureDensityParam", function(object){
    return(object@bw)})
##' @aliases bw<-
##' @rdname groupFeatures-density
setReplaceMethod("bw", "FeatureDensityParam", function(object, value) {
    object@bw <- value
    if (validObject(object))
        return(object)
})

##' @aliases minFraction
##' @description \code{minFraction},\code{minFraction<-}: getter and setter for
##' the \code{minFraction} slot of the object.
##' @rdname groupFeatures-density
setMethod("minFraction", "FeatureDensityParam", function(object){
    return(object@minFraction)})
##' @aliases minFraction<-
##' @rdname groupFeatures-density
setReplaceMethod("minFraction", "FeatureDensityParam", function(object, value) {
    object@minFraction <- value
    if (validObject(object))
        return(object)
})

##' @aliases minSamples
##' @description \code{minSamples},\code{minSamples<-}: getter and setter for the
##' \code{minSamples} slot of the object.
##' @rdname groupFeatures-density
setMethod("minSamples", "FeatureDensityParam", function(object){
    return(object@minSamples)})
##' @aliases minSamples<-
##' @rdname groupFeatures-density
setReplaceMethod("minSamples", "FeatureDensityParam", function(object, value) {
    object@minSamples <- value
    if (validObject(object))
        return(object)
})

##' @description \code{binSize},\code{binSize<-}: getter and setter for the
##' \code{binSize} slot of the object.
##' @rdname groupFeatures-density
setMethod("binSize", "FeatureDensityParam", function(object){
    return(object@binSize)})
##' @rdname groupFeatures-density
setReplaceMethod("binSize", "FeatureDensityParam", function(object, value) {
    object@binSize <- value
    if (validObject(object))
        return(object)
})

##' @aliases maxFeatures
##' @description \code{maxFeatures},\code{maxFeatures<-}: getter and setter for
##' the \code{maxFeatures} slot of the object.
##' @rdname groupFeatures-density
setMethod("maxFeatures", "FeatureDensityParam", function(object){
    return(object@maxFeatures)})
##' @aliases maxFeatures<-
##' @rdname groupFeatures-density
setReplaceMethod("maxFeatures", "FeatureDensityParam", function(object, value) {
    object@maxFeatures <- value
    if (validObject(object))
        return(object)
})


############################################################
## MzClustParam
setMethod("initialize", "MzClustParam", function(.Object, ...) {
    classVersion(.Object)["MzClustParam"] <- "0.0.1"
    callNextMethod(.Object, ...)
})

## ##' @rdname groupFeatures-mzClust
## setMethod("print", "CentWaveParam", function(x, ...) show(x))
##' @rdname groupFeatures-mzClust
setMethod("show", "MzClustParam", function(object) {
    cat("Object of class: ", class(object), "\n")
    cat("Parameters:\n")
    cat(" sampleGroups:", class(object@sampleGroups), "of length",
        length(object@sampleGroups), "\n")
    cat(" ppm:", object@ppm, "\n")
    cat(" absMz:", object@absMz, "\n")
    cat(" minFraction:", minFraction(object), "\n")
    cat(" minSamples:", minSamples(object), "\n")
})

##' @description \code{sampleGroups},\code{sampleGroups<-}: getter and setter
##' for the \code{sampleGroups} slot of the object.
##' @rdname groupFeatures-mzClust
setMethod("sampleGroups", "MzClustParam", function(object){
    return(object@sampleGroups)})
##' @param value The value for the slot.
##' @rdname groupFeatures-mzClust
setReplaceMethod("sampleGroups", "MzClustParam", function(object, value) {
    object@sampleGroups <- value
    if (validObject(object))
        return(object)
})

##' @description \code{ppm},\code{ppm<-}: getter and setter for the \code{ppm}
##' slot of the object.
##' @rdname groupFeatures-mzClust
setMethod("ppm", "MzClustParam", function(object){
    return(object@ppm)})
##' @rdname groupFeatures-mzClust
setReplaceMethod("ppm", "MzClustParam", function(object, value) {
    object@ppm <- value
    if (validObject(object))
        return(object)
})

##' @aliases absMz
##' @description \code{absMz},\code{absMz<-}: getter and setter for the
##' \code{absMz} slot of the object.
##' @rdname groupFeatures-mzClust
setMethod("absMz", "MzClustParam", function(object){
    return(object@absMz)})
##' @aliases absMz<-
##' @rdname groupFeatures-mzClust
setReplaceMethod("absMz", "MzClustParam", function(object, value) {
    object@absMz <- value
    if (validObject(object))
        return(object)
})

##' @description \code{minFraction},\code{minFraction<-}: getter and setter for
##' the \code{minFraction} slot of the object.
##' @rdname groupFeatures-mzClust
setMethod("minFraction", "MzClustParam", function(object){
    return(object@minFraction)})
##' @rdname groupFeatures-mzClust
setReplaceMethod("minFraction", "MzClustParam", function(object, value) {
    object@minFraction <- value
    if (validObject(object))
        return(object)
})

##' @description \code{minSamples},\code{minSamples<-}: getter and setter for the
##' \code{minSamples} slot of the object.
##' @rdname groupFeatures-mzClust
setMethod("minSamples", "MzClustParam", function(object){
    return(object@minSamples)})
##' @rdname groupFeatures-mzClust
setReplaceMethod("minSamples", "MzClustParam", function(object, value) {
    object@minSamples <- value
    if (validObject(object))
        return(object)
})


############################################################
## NearestFeaturesParam
setMethod("initialize", "NearestFeaturesParam", function(.Object, ...) {
    classVersion(.Object)["NearestFeaturesParam"] <- "0.0.1"
    callNextMethod(.Object, ...)
})

##' @rdname groupFeatures-nearest
setMethod("show", "NearestFeaturesParam", function(object) {
    cat("Object of class: ", class(object), "\n")
    cat("Parameters:\n")
    cat(" sampleGroups:", class(object@sampleGroups), "of length",
        length(object@sampleGroups), "\n")
    cat(" mzVsRtBalance:", object@mzVsRtBalance, "\n")
    cat(" absMz:", object@absMz, "\n")
    cat(" absRt:", object@absRt, "\n")
    cat(" kNN:", object@kNN, "\n")
})

##' @description \code{sampleGroups},\code{sampleGroups<-}: getter and setter
##' for the \code{sampleGroups} slot of the object.
##' @rdname groupFeatures-nearest
setMethod("sampleGroups", "NearestFeaturesParam", function(object){
    return(object@sampleGroups)})
##' @param value The value for the slot.
##' @rdname groupFeatures-nearest
setReplaceMethod("sampleGroups", "NearestFeaturesParam", function(object, value) {
    object@sampleGroups <- value
    if (validObject(object))
        return(object)
})

##' @aliases mzVsRtBalance
##' @description \code{mzVsRtBalance},\code{mzVsRtBalance<-}: getter and setter
##' for the \code{mzVsRtBalance} slot of the object.
##' @rdname groupFeatures-nearest
setMethod("mzVsRtBalance", "NearestFeaturesParam", function(object){
    return(object@mzVsRtBalance)})
##' @aliases mzVsRtBalance<-
##' @rdname groupFeatures-nearest
setReplaceMethod("mzVsRtBalance", "NearestFeaturesParam", function(object, value) {
    object@mzVsRtBalance <- value
    if (validObject(object))
        return(object)
})

##' @description \code{absMz},\code{absMz<-}: getter and setter for the
##' \code{absMz} slot of the object.
##' @rdname groupFeatures-nearest
setMethod("absMz", "NearestFeaturesParam", function(object){
    return(object@absMz)})
##' @rdname groupFeatures-nearest
setReplaceMethod("absMz", "NearestFeaturesParam", function(object, value) {
    object@absMz <- value
    if (validObject(object))
        return(object)
})

##' @aliases absRt
##' @description \code{absRt},\code{absRt<-}: getter and setter for the
##' \code{absRt} slot of the object.
##' @rdname groupFeatures-nearest
setMethod("absRt", "NearestFeaturesParam", function(object){
    return(object@absRt)})
##' @aliases absRt<-
##' @rdname groupFeatures-nearest
setReplaceMethod("absRt", "NearestFeaturesParam", function(object, value) {
    object@absRt <- value
    if (validObject(object))
        return(object)
})

##' @aliases kNN
##' @description \code{kNN},\code{kNN<-}: getter and setter for the
##' \code{kNN} slot of the object.
##' @rdname groupFeatures-nearest
setMethod("kNN", "NearestFeaturesParam", function(object){
    return(object@kNN)})
##' @aliases kNN<-
##' @rdname groupFeatures-nearest
setReplaceMethod("kNN", "NearestFeaturesParam", function(object, value) {
    object@kNN <- value
    if (validObject(object))
        return(object)
})

