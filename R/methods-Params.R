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
## GenericParam
###
setMethod("initialize", "GenericParam", function(.Object, ...) {
    classVersion(.Object)["GenericParam"] <- "0.0.1"
    callNextMethod(.Object, ...)
})
#' @param object \code{GenericParam} object.
#' 
#' @rdname GenericParam
setMethod("show", "GenericParam", function(object) {
    cat("Object of class: ", class(object), "\n")
    cat(" fun:", object@fun, "\n")
    cat(" arguments:\n")
    if (length(object@args) > 0) {
        for (i in 1:length(object@args)) {
            if (!is.null(names(object@args)))
                cat(" ", names(object@args)[i], "= ")
            cat(object@args[[i]], "\n")
        }
    }
})

############################################################
## CentWaveParam
###
setMethod("initialize", "CentWaveParam", function(.Object, ...) {
    classVersion(.Object)["CentWaveParam"] <- "0.0.1"
    callNextMethod(.Object, ...)
})

#' @rdname findChromPeaks-centWave
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
setMethod("initialize", "MatchedFilterParam", function(.Object, ...) {
    classVersion(.Object)["MatchedFilterParam"] <- "0.0.1"
    callNextMethod(.Object, ...)
})
#' @rdname findChromPeaks-matchedFilter
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
setMethod("initialize", "MassifquantParam", function(.Object, ...) {
    classVersion(.Object)["MassifquantParam"] <- "0.0.1"
    callNextMethod(.Object, ...)
})

#' @rdname findChromPeaks-massifquant
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
setMethod("initialize", "MSWParam", function(.Object, ...) {
    classVersion(.Object)["MSWParam"] <- "0.0.1"
    callNextMethod(.Object, ...)
})

#' @rdname findPeaks-MSW
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
#'     \code{\link{sav.gol}} functions from the
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

#' @rdname findChromPeaks-centWaveWithPredIsoROIs
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
setMethod("initialize", "PeakDensityParam", function(.Object, ...) {
    classVersion(.Object)["PeakDensityParam"] <- "0.0.1"
    callNextMethod(.Object, ...)
})

#' @rdname groupChromPeaks-density
setMethod("show", "PeakDensityParam", function(object) {
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

#' @aliases sampleGroups
#' 
#' @description \code{sampleGroups},\code{sampleGroups<-}: getter and setter
#'     for the \code{sampleGroups} slot of the object. Its length should match
#'     the number of samples in the experiment and it should not contain
#'     \code{NA}s.
#' 
#' @rdname groupChromPeaks-density
setMethod("sampleGroups", "PeakDensityParam", function(object){
    return(object@sampleGroups)})
#' @aliases sampleGroups<-
#' 
#' @param value The value for the slot.
#' 
#' @rdname groupChromPeaks-density
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
#' @description \code{bw},\code{bw<-}: getter and setter for the \code{bw} slot
#'     of the object.
#' 
#' @rdname groupChromPeaks-density
setMethod("bw", "PeakDensityParam", function(object){
    return(object@bw)})
#' @aliases bw<-
#' 
#' @rdname groupChromPeaks-density
setReplaceMethod("bw", "PeakDensityParam", function(object, value) {
    object@bw <- value
    if (validObject(object))
        return(object)
})

#' @aliases minFraction
#' 
#' @description \code{minFraction},\code{minFraction<-}: getter and setter for
#'     the \code{minFraction} slot of the object.
#' 
#' @rdname groupChromPeaks-density
setMethod("minFraction", "PeakDensityParam", function(object){
    return(object@minFraction)})
#' @aliases minFraction<-
#' 
#' @rdname groupChromPeaks-density
setReplaceMethod("minFraction", "PeakDensityParam", function(object, value) {
    object@minFraction <- value
    if (validObject(object))
        return(object)
})

#' @aliases minSamples
#' 
#' @description \code{minSamples},\code{minSamples<-}: getter and setter for the
#'     \code{minSamples} slot of the object.
#' 
#' @rdname groupChromPeaks-density
setMethod("minSamples", "PeakDensityParam", function(object){
    return(object@minSamples)})
#' @aliases minSamples<-
#' 
#' @rdname groupChromPeaks-density
setReplaceMethod("minSamples", "PeakDensityParam", function(object, value) {
    object@minSamples <- value
    if (validObject(object))
        return(object)
})

#' @description \code{binSize},\code{binSize<-}: getter and setter for the
#'     \code{binSize} slot of the object.
#' 
#' @rdname groupChromPeaks-density
setMethod("binSize", "PeakDensityParam", function(object){
    return(object@binSize)})
#' @rdname groupChromPeaks-density
setReplaceMethod("binSize", "PeakDensityParam", function(object, value) {
    object@binSize <- value
    if (validObject(object))
        return(object)
})

#' @aliases maxFeatures
#' 
#' @description \code{maxFeatures},\code{maxFeatures<-}: getter and setter for
#'     the \code{maxFeatures} slot of the object.
#' 
#' @rdname groupChromPeaks-density
setMethod("maxFeatures", "PeakDensityParam", function(object){
    return(object@maxFeatures)})
#' @aliases maxFeatures<-
#' 
#' @rdname groupChromPeaks-density
setReplaceMethod("maxFeatures", "PeakDensityParam", function(object, value) {
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

#' @rdname groupChromPeaks-mzClust
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

#' @description \code{sampleGroups},\code{sampleGroups<-}: getter and setter
#'     for the \code{sampleGroups} slot of the object.
#' 
#' @rdname groupChromPeaks-mzClust
setMethod("sampleGroups", "MzClustParam", function(object){
    return(object@sampleGroups)})
#' @param value The value for the slot.
#' 
#' @rdname groupChromPeaks-mzClust
setReplaceMethod("sampleGroups", "MzClustParam", function(object, value) {
    object@sampleGroups <- value
    if (validObject(object))
        return(object)
})

#' @description \code{ppm},\code{ppm<-}: getter and setter for the \code{ppm}
#'     slot of the object.
#' 
#' @rdname groupChromPeaks-mzClust
setMethod("ppm", "MzClustParam", function(object){
    return(object@ppm)})
#' @rdname groupChromPeaks-mzClust
setReplaceMethod("ppm", "MzClustParam", function(object, value) {
    object@ppm <- value
    if (validObject(object))
        return(object)
})

#' @aliases absMz
#' 
#' @description \code{absMz},\code{absMz<-}: getter and setter for the
#'     \code{absMz} slot of the object.
#' 
#' @rdname groupChromPeaks-mzClust
setMethod("absMz", "MzClustParam", function(object){
    return(object@absMz)})
#' @aliases absMz<-
#' 
#' @rdname groupChromPeaks-mzClust
setReplaceMethod("absMz", "MzClustParam", function(object, value) {
    object@absMz <- value
    if (validObject(object))
        return(object)
})

#' @description \code{minFraction},\code{minFraction<-}: getter and setter for
#'     the \code{minFraction} slot of the object.
#' 
#' @rdname groupChromPeaks-mzClust
setMethod("minFraction", "MzClustParam", function(object){
    return(object@minFraction)})
#' @rdname groupChromPeaks-mzClust
setReplaceMethod("minFraction", "MzClustParam", function(object, value) {
    object@minFraction <- value
    if (validObject(object))
        return(object)
})

#' @description \code{minSamples},\code{minSamples<-}: getter and setter for the
#'     \code{minSamples} slot of the object.
#' 
#' @rdname groupChromPeaks-mzClust
setMethod("minSamples", "MzClustParam", function(object){
    return(object@minSamples)})
#' @rdname groupChromPeaks-mzClust
setReplaceMethod("minSamples", "MzClustParam", function(object, value) {
    object@minSamples <- value
    if (validObject(object))
        return(object)
})


############################################################
## NearestPeaksParam
setMethod("initialize", "NearestPeaksParam", function(.Object, ...) {
    classVersion(.Object)["NearestPeaksParam"] <- "0.0.1"
    callNextMethod(.Object, ...)
})

#' @rdname groupChromPeaks-nearest
setMethod("show", "NearestPeaksParam", function(object) {
    cat("Object of class: ", class(object), "\n")
    cat("Parameters:\n")
    cat(" sampleGroups:", class(object@sampleGroups), "of length",
        length(object@sampleGroups), "\n")
    cat(" mzVsRtBalance:", object@mzVsRtBalance, "\n")
    cat(" absMz:", object@absMz, "\n")
    cat(" absRt:", object@absRt, "\n")
    cat(" kNN:", object@kNN, "\n")
})

#' @description \code{sampleGroups},\code{sampleGroups<-}: getter and setter
#'     for the \code{sampleGroups} slot of the object.
#' 
#' @rdname groupChromPeaks-nearest
setMethod("sampleGroups", "NearestPeaksParam", function(object){
    return(object@sampleGroups)})
#' @param value The value for the slot.
#' 
#' @rdname groupChromPeaks-nearest
setReplaceMethod("sampleGroups", "NearestPeaksParam", function(object, value) {
    object@sampleGroups <- value
    if (validObject(object))
        return(object)
})

#' @aliases mzVsRtBalance
#' 
#' @description \code{mzVsRtBalance},\code{mzVsRtBalance<-}: getter and setter
#'     for the \code{mzVsRtBalance} slot of the object.
#' 
#' @rdname groupChromPeaks-nearest
setMethod("mzVsRtBalance", "NearestPeaksParam", function(object){
    return(object@mzVsRtBalance)})
#' @aliases mzVsRtBalance<-
#' 
#' @rdname groupChromPeaks-nearest
setReplaceMethod("mzVsRtBalance", "NearestPeaksParam", function(object, value) {
    object@mzVsRtBalance <- value
    if (validObject(object))
        return(object)
})

#' @description \code{absMz},\code{absMz<-}: getter and setter for the
#'     \code{absMz} slot of the object.
#' 
#' @rdname groupChromPeaks-nearest
setMethod("absMz", "NearestPeaksParam", function(object){
    return(object@absMz)})
#' @rdname groupChromPeaks-nearest
setReplaceMethod("absMz", "NearestPeaksParam", function(object, value) {
    object@absMz <- value
    if (validObject(object))
        return(object)
})

#' @aliases absRt
#' 
#' @description \code{absRt},\code{absRt<-}: getter and setter for the
#'     \code{absRt} slot of the object.
#' 
#' @rdname groupChromPeaks-nearest
setMethod("absRt", "NearestPeaksParam", function(object){
    return(object@absRt)})
#' @aliases absRt<-
#' 
#' @rdname groupChromPeaks-nearest
setReplaceMethod("absRt", "NearestPeaksParam", function(object, value) {
    object@absRt <- value
    if (validObject(object))
        return(object)
})

#' @aliases kNN
#' 
#' @description \code{kNN},\code{kNN<-}: getter and setter for the
#'     \code{kNN} slot of the object.
#' 
#' @rdname groupChromPeaks-nearest
setMethod("kNN", "NearestPeaksParam", function(object){
    return(object@kNN)})
#' @aliases kNN<-
#' 
#' @rdname groupChromPeaks-nearest
setReplaceMethod("kNN", "NearestPeaksParam", function(object, value) {
    object@kNN <- value
    if (validObject(object))
        return(object)
})


############################################################
## PeakGroupsParam
setMethod("initialize", "PeakGroupsParam", function(.Object, ...) {
    classVersion(.Object)["PeakGroupsParam"] <- "0.0.2"
    callNextMethod(.Object, ...)
})

#' @rdname adjustRtime-peakGroups
setMethod("show", "PeakGroupsParam", function(object) {
    cat("Object of class: ", class(object), "\n")
    cat("Parameters:\n")
    cat(" minFraction:", object@minFraction, "\n")
    cat(" extraPeaks:", object@extraPeaks, "\n")
    cat(" smooth:", object@smooth, "\n")
    cat(" span:", object@span, "\n")
    cat(" family:", object@family, "\n")
    pgm <- peakGroupsMatrix(object)
    if (nrow(pgm))
        cat(" number of peak groups:", nrow(pgm), "\n")
})

#' @description \code{minFraction},\code{minFraction<-}: getter and setter
#'     for the \code{minFraction} slot of the object.
#' 
#' @rdname adjustRtime-peakGroups
setMethod("minFraction", "PeakGroupsParam", function(object){
    return(object@minFraction)})
#' @param value The value for the slot.
#' 
#' @rdname adjustRtime-peakGroups
setReplaceMethod("minFraction", "PeakGroupsParam", function(object, value) {
    object@minFraction <- value
    if (validObject(object))
        return(object)
})

#' @aliases extraPeaks
#' 
#' @description \code{extraPeaks},\code{extraPeaks<-}: getter and setter
#'     for the \code{extraPeaks} slot of the object.
#' 
#' @rdname adjustRtime-peakGroups
setMethod("extraPeaks", "PeakGroupsParam", function(object){
    return(object@extraPeaks)})
#' @aliases extraPeaks<-
#' 
#' @rdname adjustRtime-peakGroups
setReplaceMethod("extraPeaks", "PeakGroupsParam", function(object, value) {
    object@extraPeaks <- value
    if (validObject(object))
        return(object)
})

#' @aliases smooth
#' 
#' @description \code{smooth},\code{smooth<-}: getter and setter
#'     for the \code{smooth} slot of the object.
#' 
#' @param x a \code{PeakGroupsParam} object.
#' 
#' @rdname adjustRtime-peakGroups
setMethod("smooth", "PeakGroupsParam", function(x){
    return(x@smooth)})
#' @aliases smooth<-
#' 
#' @rdname adjustRtime-peakGroups
setReplaceMethod("smooth", "PeakGroupsParam", function(object, value) {
    object@smooth <- value
    if (validObject(object))
        return(object)
})

#' @aliases span
#' 
#' @description \code{span},\code{span<-}: getter and setter
#'     for the \code{span} slot of the object.
#' 
#' @rdname adjustRtime-peakGroups
setMethod("span", "PeakGroupsParam", function(object){
    return(object@span)})
#' @aliases span<-
#' 
#' @rdname adjustRtime-peakGroups
setReplaceMethod("span", "PeakGroupsParam", function(object, value) {
    object@span <- value
    if (validObject(object))
        return(object)
})

#' @aliases family
#' 
#' @description \code{family},\code{family<-}: getter and setter
#'     for the \code{family} slot of the object.
#' 
#' @rdname adjustRtime-peakGroups
setMethod("family", "PeakGroupsParam", function(object){
    return(object@family)})
#' @aliases family<-
#' 
#' @rdname adjustRtime-peakGroups
setReplaceMethod("family", "PeakGroupsParam", function(object, value) {
    object@family <- value
    if (validObject(object))
        return(object)
})

#' @aliases peakGroupsMatrix
#' 
#' @description \code{peakGroupsMatrix},\code{peakGroupsMatrix<-}: getter and
#'     setter for the \code{peakGroupsMatrix} slot of the object.
#' 
#' @rdname adjustRtime-peakGroups
setMethod("peakGroupsMatrix", "PeakGroupsParam", function(object){
    return(object@peakGroupsMatrix)})
#' @aliases peakGroupsMatrix<-
#' 
#' @rdname adjustRtime-peakGroups
setReplaceMethod("peakGroupsMatrix", "PeakGroupsParam", function(object, value) {
    object@peakGroupsMatrix <- value
    if (validObject(object))
        return(object)
})


############################################################
## ObiwarpParam
setMethod("initialize", "ObiwarpParam", function(.Object, ...) {
    classVersion(.Object)["ObiwarpParam"] <- "0.0.1"
    callNextMethod(.Object, ...)
})

#' @rdname adjustRtime-obiwarp
setMethod("show", "ObiwarpParam", function(object) {
    cat("Object of class: ", class(object), "\n")
    cat("Parameters:\n")
    cat(" binSize:", binSize(object), "\n")
    cat(" centerSample:", centerSample(object), "\n")
    cat(" response:", response(object), "\n")
    cat(" distFun:", distFun(object), "\n")
    cat(" gapInit:", gapInit(object), "\n")
    cat(" gapExtend:", gapExtend(object), "\n")
    cat(" factorDiag:", factorDiag(object), "\n")
    cat(" factorGap:", factorGap(object), "\n")
    cat(" localAlignment:", localAlignment(object), "\n")
    cat(" initPenalty:", initPenalty(object), "\n")
})

#' @description \code{binSize},\code{binSize<-}: getter and setter
#'     for the \code{binSize} slot of the object.
#' 
#' @rdname adjustRtime-obiwarp
setMethod("binSize", "ObiwarpParam", function(object){
    return(object@binSize)})
#' @param value The value for the slot.
#' 
#' @rdname adjustRtime-obiwarp
setReplaceMethod("binSize", "ObiwarpParam", function(object, value) {
    object@binSize <- value
    if (validObject(object))
        return(object)
})

#' @aliases centerSample
#' 
#' @description \code{centerSample},\code{centerSample<-}: getter and setter
#'     for the \code{centerSample} slot of the object.
#' 
#' @rdname adjustRtime-obiwarp
setMethod("centerSample", "ObiwarpParam", function(object){
    return(object@centerSample)})
#' @aliases centerSample<-
#' 
#' @rdname adjustRtime-obiwarp
setReplaceMethod("centerSample", "ObiwarpParam", function(object, value) {
    object@centerSample <- as.integer(value)
    if (validObject(object))
        return(object)
})

#' @aliases response
#' 
#' @description \code{response},\code{response<-}: getter and setter
#'     for the \code{response} slot of the object.
#' 
#' @rdname adjustRtime-obiwarp
setMethod("response", "ObiwarpParam", function(object){
    return(object@response)})
#' @aliases response<-
#' 
#' @rdname adjustRtime-obiwarp
setReplaceMethod("response", "ObiwarpParam", function(object, value) {
    object@response <- as.integer(value)
    if (validObject(object))
        return(object)
})

#' @aliases distFun
#' 
#' @description \code{distFun},\code{distFun<-}: getter and setter
#'     for the \code{distFun} slot of the object.
#' 
#' @rdname adjustRtime-obiwarp
setMethod("distFun", "ObiwarpParam", function(object){
    return(object@distFun)})
#' @aliases distFun<-
#' 
#' @rdname adjustRtime-obiwarp
setReplaceMethod("distFun", "ObiwarpParam", function(object, value) {
    object@distFun <- value
    if (validObject(object))
        return(object)
})

#' @aliases gapInit
#' 
#' @description \code{gapInit},\code{gapInit<-}: getter and setter
#'     for the \code{gapInit} slot of the object.
#' 
#' @rdname adjustRtime-obiwarp
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
#' @rdname adjustRtime-obiwarp
setReplaceMethod("gapInit", "ObiwarpParam", function(object, value) {
    object@gapInit <- value
    if (validObject(object))
        return(object)
})

#' @aliases gapExtend
#' 
#' @description \code{gapExtend},\code{gapExtend<-}: getter and setter
#'     for the \code{gapExtend} slot of the object.
#' 
#' @rdname adjustRtime-obiwarp
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
#' @rdname adjustRtime-obiwarp
setReplaceMethod("gapExtend", "ObiwarpParam", function(object, value) {
    object@gapExtend <- value
    if (validObject(object))
        return(object)
})

#' @aliases factorDiag
#' 
#' @description \code{factorDiag},\code{factorDiag<-}: getter and setter
#'      for the \code{factorDiag} slot of the object.
#' 
#' @rdname adjustRtime-obiwarp
setMethod("factorDiag", "ObiwarpParam", function(object){
    return(object@factorDiag)})
#' @aliases factorDiag<-
#' 
#' @rdname adjustRtime-obiwarp
setReplaceMethod("factorDiag", "ObiwarpParam", function(object, value) {
    object@factorDiag <- value
    if (validObject(object))
        return(object)
})

#' @aliases factorGap
#' 
#' @description \code{factorGap},\code{factorGap<-}: getter and setter
#'     for the \code{factorGap} slot of the object.
#' 
#' @rdname adjustRtime-obiwarp
setMethod("factorGap", "ObiwarpParam", function(object){
    return(object@factorGap)})
#' @aliases factorGap<-
#' 
#' @rdname adjustRtime-obiwarp
setReplaceMethod("factorGap", "ObiwarpParam", function(object, value) {
    object@factorGap <- value
    if (validObject(object))
        return(object)
})

#' @aliases localAlignment
#' 
#' @description \code{localAlignment},\code{localAlignment<-}: getter and setter
#'     for the \code{localAlignment} slot of the object.
#' 
#' @rdname adjustRtime-obiwarp
setMethod("localAlignment", "ObiwarpParam", function(object){
    return(object@localAlignment)})
#' @aliases localAlignment<-
#' 
#' @rdname adjustRtime-obiwarp
setReplaceMethod("localAlignment", "ObiwarpParam", function(object, value) {
    object@localAlignment <- value
    if (validObject(object))
        return(object)
})

#' @aliases initPenalty
#' 
#' @description \code{initPenalty},\code{initPenalty<-}: getter and setter
#'     for the \code{initPenalty} slot of the object.
#' 
#' @rdname adjustRtime-obiwarp
setMethod("initPenalty", "ObiwarpParam", function(object){
    return(object@initPenalty)})
#' @aliases initPenalty<-
#' 
#' @rdname adjustRtime-obiwarp
setReplaceMethod("initPenalty", "ObiwarpParam", function(object, value) {
    object@initPenalty <- value
    if (validObject(object))
        return(object)
})

############################################################
## FillChromPeaksParam
###
setMethod("initialize", "FillChromPeaksParam", function(.Object, ...) {
    classVersion(.Object)["FillChromPeaksParam"] <- "0.0.1"
    callNextMethod(.Object, ...)
})
#' @rdname fillChromPeaks
setMethod("show", "FillChromPeaksParam", function(object) {
    cat("Object of class: ", class(object), "\n")
    cat("Parameters:\n")
    cat(" expandMz:", object@expandMz, "\n")
    cat(" expandRt:", object@expandRt, "\n")
    cat(" ppm:", object@ppm, "\n")
})

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
