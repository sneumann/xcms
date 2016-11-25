## Methods for the Param class and sub-classes
#' @include functions-Params.R

############################################################
## Param
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
## All the getter and setter methods.
## o ppm
##' @description \code{ppm},\code{ppm<-}: getter and setter for the \code{ppm}
##' slot of the object.
##'
##' @rdname featureDetection-centWave
setMethod("ppm", "CentWaveParam", function(object){ return(object@ppm)})
##' @param value The value for the slot.
##'
##' @rdname featureDetection-centWave
setReplaceMethod("ppm", "CentWaveParam", function(object, value) {
    object@ppm <- value
    OK <- validObject(object)
    if (is.character(OK))
        stop(OK)
    return(object)
})
## o peakwidth
##' @description \code{peakwidth},\code{peakwidth<-}: getter and setter for the
##' \code{peakwidth} slot of the object.
##'
##' @rdname featureDetection-centWave
setMethod("peakwidth", "CentWaveParam", function(object)
    return(object@peakwidth))
##' @rdname featureDetection-centWave
setReplaceMethod("peakwidth", "CentWaveParam", function(object, value) {
    object@peakwidth <- value
    OK <- validObject(object)
    if (is.character(OK))
        stop(OK)
    return(object)
})
## o snthresh
##' @description \code{snthresh},\code{snthresh<-}: getter and setter for the
##' \code{snthresh} slot of the object.
##'
##' @rdname featureDetection-centWave
setMethod("snthresh", "CentWaveParam", function(object)
    return(object@snthresh))
##' @rdname featureDetection-centWave
setReplaceMethod("snthresh", "CentWaveParam", function(object, value) {
    object@snthresh <- value
    OK <- validObject(object)
    if (is.character(OK))
        stop(OK)
    return(object)
})
## o prefilter
##' @description \code{prefilter},\code{prefilter<-}: getter and setter for the
##' \code{prefilter} slot of the object.
##'
##' @rdname featureDetection-centWave
setMethod("prefilter", "CentWaveParam", function(object)
    return(object@prefilter))
##' @rdname featureDetection-centWave
setReplaceMethod("prefilter", "CentWaveParam", function(object, value) {
    object@prefilter <- value
    OK <- validObject(object)
    if (is.character(OK))
        stop(OK)
    return(object)
})
## o mzCenterFun
##' @description \code{mzCenterFun},\code{mzCenterFun<-}: getter and setter for the
##' \code{mzCenterFun} slot of the object.
##'
##' @rdname featureDetection-centWave
setMethod("mzCenterFun", "CentWaveParam", function(object)
    return(object@mzCenterFun))
##' @rdname featureDetection-centWave
setReplaceMethod("mzCenterFun", "CentWaveParam", function(object, value) {
    object@mzCenterFun <- value
    OK <- validObject(object)
    if (is.character(OK))
        stop(OK)
    return(object)
})
## o integrate
##' @description \code{integrate},\code{integrate<-}: getter and setter for the
##' \code{integrate} slot of the object.
##'
##' @param f For \code{integrate}: a \code{CentWaveParam} object.
##'
##' @rdname featureDetection-centWave
setMethod("integrate", signature(f = "CentWaveParam"), function(f)
    return(f@integrate))
##' @rdname featureDetection-centWave
setReplaceMethod("integrate", "CentWaveParam", function(object, value) {
    if (!is.integer(value))
        value <- as.integer(value)
    object@integrate <- value
    OK <- validObject(object)
    if (is.character(OK))
        stop(OK)
    return(object)
})
## o mzdiff
##' @description \code{mzdiff},\code{mzdiff<-}: getter and setter for the
##' \code{mzdiff} slot of the object.
##'
##' @rdname featureDetection-centWave
setMethod("mzdiff", "CentWaveParam", function(object)
    return(object@mzdiff))
##' @rdname featureDetection-centWave
setReplaceMethod("mzdiff", "CentWaveParam", function(object, value) {
    object@mzdiff <- value
    OK <- validObject(object)
    if (is.character(OK))
        stop(OK)
    return(object)
})
## o fitgauss
##' @description \code{fitgauss},\code{fitgauss<-}: getter and setter for the
##' \code{fitgauss} slot of the object.
##'
##' @rdname featureDetection-centWave
setMethod("fitgauss", "CentWaveParam", function(object)
    return(object@fitgauss))
##' @rdname featureDetection-centWave
setReplaceMethod("fitgauss", "CentWaveParam", function(object, value) {
    object@fitgauss <- value
    OK <- validObject(object)
    if (is.character(OK))
        stop(OK)
    return(object)
})
## o noise
##' @description \code{noise},\code{noise<-}: getter and setter for the
##' \code{noise} slot of the object.
##'
##' @rdname featureDetection-centWave
setMethod("noise", "CentWaveParam", function(object)
    return(object@noise))
##' @rdname featureDetection-centWave
setReplaceMethod("noise", "CentWaveParam", function(object, value) {
    object@noise <- value
    OK <- validObject(object)
    if (is.character(OK))
        stop(OK)
    return(object)
})
## o verboseColumns
##' @description \code{verboseColumns},\code{verboseColumns<-}: getter and
##' setter for the \code{verboseColumns} slot of the object.
##'
##' @rdname featureDetection-centWave
setMethod("verboseColumns", "CentWaveParam", function(object)
    return(object@verboseColumns))
##' @rdname featureDetection-centWave
setReplaceMethod("verboseColumns", "CentWaveParam", function(object, value) {
    object@verboseColumns <- value
    OK <- validObject(object)
    if (is.character(OK))
        stop(OK)
    return(object)
})
## o roiList
##' @description \code{roiList},\code{roiList<-}: getter and setter for the
##' \code{roiList} slot of the object.
##'
##' @rdname featureDetection-centWave
setMethod("roiList", "CentWaveParam", function(object)
    return(object@roiList))
##' @rdname featureDetection-centWave
setReplaceMethod("roiList", "CentWaveParam", function(object, value) {
    object@roiList <- value
    OK <- validObject(object)
    if (is.character(OK))
        stop(OK)
    return(object)
})
## o firstBaselineCheck
##' @description \code{fistBaselineCheck},\code{firstBaselineCheck<-}: getter
##' and setter for the \code{firstBaselineCheck} slot of the object.
##'
##' @rdname featureDetection-centWave
setMethod("firstBaselineCheck", "CentWaveParam", function(object)
    return(object@firstBaselineCheck))
##' @rdname featureDetection-centWave
setReplaceMethod("firstBaselineCheck", "CentWaveParam", function(object, value) {
    object@firstBaselineCheck <- value
    OK <- validObject(object)
    if (is.character(OK))
        stop(OK)
    return(object)
})
## o roiScales
##' @description \code{roiScales},\code{roiScales<-}: getter and setter for the
##' \code{roiScales} slot of the object.
##'
##' @rdname featureDetection-centWave
setMethod("roiScales", "CentWaveParam", function(object)
    return(object@roiScales))
##' @rdname featureDetection-centWave
setReplaceMethod("roiScales", "CentWaveParam", function(object, value) {
    object@roiScales <- value
    OK <- validObject(object)
    if (is.character(OK))
        stop(OK)
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
##' @description \code{binSize},\code{binSize<-}: getter and setter for the
##' \code{binSize} slot of the object.
##'
##' @rdname featureDetection-matchedFilter
setMethod("binSize", "MatchedFilterParam", function(object)
    return(object@binSize))
##' @param value The value for the slot.
##'
##' @rdname featureDetection-matchedFilter
setReplaceMethod("binSize", "MatchedFilterParam", function(object, value) {
    object@binSize <- value
    OK <- validObject(object)
    if (is.character(OK))
        stop(OK)
    return(object)
})
##' @description \code{impute},\code{impute<-}: getter and setter for the
##' \code{impute} slot of the object.
##'
##' @rdname featureDetection-matchedFilter
setMethod("impute", "MatchedFilterParam", function(object)
    return(object@impute))
##' @rdname featureDetection-matchedFilter
setReplaceMethod("impute", "MatchedFilterParam", function(object, value) {
    object@impute <- value
    OK <- validObject(object)
    if (is.character(OK))
        stop(OK)
    return(object)
})
##' @description \code{baseValue},\code{baseValue<-}: getter and setter for the
##' \code{baseValue} slot of the object.
##'
##' @rdname featureDetection-matchedFilter
setMethod("baseValue", "MatchedFilterParam", function(object)
    return(object@baseValue))
##' @rdname featureDetection-matchedFilter
setReplaceMethod("baseValue", "MatchedFilterParam", function(object, value) {
    object@baseValue <- value
    OK <- validObject(object)
    if (is.character(OK))
        stop(OK)
    return(object)
})
##' @description \code{distance},\code{distance<-}: getter and setter for the
##' \code{distance} slot of the object.
##'
##' @rdname featureDetection-matchedFilter
setMethod("distance", "MatchedFilterParam", function(object)
    return(object@distance))
##' @rdname featureDetection-matchedFilter
setReplaceMethod("distance", "MatchedFilterParam", function(object, value) {
    object@distance <- value
    OK <- validObject(object)
    if (is.character(OK))
        stop(OK)
    return(object)
})
##' @description \code{fwhm},\code{fwhm<-}: getter and setter for the
##' \code{fwhm} slot of the object.
##'
##' @rdname featureDetection-matchedFilter
setMethod("fwhm", "MatchedFilterParam", function(object)
    return(object@fwhm))
##' @rdname featureDetection-matchedFilter
setReplaceMethod("fwhm", "MatchedFilterParam", function(object, value) {
    object@fwhm <- value
    OK <- validObject(object)
    if (is.character(OK))
        stop(OK)
    return(object)
})
##' @description \code{sigma},\code{sigma<-}: getter and setter for the
##' \code{sigma} slot of the object.
##'
##' @rdname featureDetection-matchedFilter
setMethod("sigma", "MatchedFilterParam", function(object)
    return(object@sigma))
##' @rdname featureDetection-matchedFilter
setReplaceMethod("sigma", "MatchedFilterParam", function(object, value) {
    object@sigma <- value
    OK <- validObject(object)
    if (is.character(OK))
        stop(OK)
    return(object)
})
##' @description \code{max},\code{max<-}: getter and setter for the
##' \code{max} slot of the object.
##'
##' @param x For \code{max}: a \code{MatchedFilterParam} object.
##'
##' @rdname featureDetection-matchedFilter
setMethod("max", signature(x="MatchedFilterParam"),
          function(x) return(x@max))
##' @rdname featureDetection-matchedFilter
setReplaceMethod("max", "MatchedFilterParam", function(object, value) {
    object@max <- value
    OK <- validObject(object)
    if (is.character(OK))
        stop(OK)
    return(object)
})
##' @description \code{snthresh},\code{snthresh<-}: getter and setter for the
##' \code{snthresh} slot of the object.
##'
##' @rdname featureDetection-matchedFilter
setMethod("snthresh", "MatchedFilterParam", function(object)
    return(object@snthresh))
##' @rdname featureDetection-matchedFilter
setReplaceMethod("snthresh", "MatchedFilterParam", function(object, value) {
    object@snthresh <- value
    OK <- validObject(object)
    if (is.character(OK))
        stop(OK)
    return(object)
})
##' @description \code{steps},\code{steps<-}: getter and setter for the
##' \code{steps} slot of the object.
##'
##' @rdname featureDetection-matchedFilter
setMethod("steps", "MatchedFilterParam", function(object)
    return(object@steps))
##' @rdname featureDetection-matchedFilter
setReplaceMethod("steps", "MatchedFilterParam", function(object, value) {
    object@steps <- value
    OK <- validObject(object)
    if (is.character(OK))
        stop(OK)
    return(object)
})
##' @description \code{mzdiff},\code{mzdiff<-}: getter and setter for the
##' \code{mzdiff} slot of the object.
##'
##' @rdname featureDetection-matchedFilter
setMethod("mzdiff", "MatchedFilterParam", function(object)
    return(object@mzdiff))
##' @rdname featureDetection-matchedFilter
setReplaceMethod("mzdiff", "MatchedFilterParam", function(object, value) {
    object@mzdiff <- value
    OK <- validObject(object)
    if (is.character(OK))
        stop(OK)
    return(object)
})
##' @description \code{index},\code{index<-}: getter and setter for the
##' \code{index} slot of the object.
##'
##' @rdname featureDetection-matchedFilter
setMethod("index", "MatchedFilterParam", function(object)
    return(object@index))
##' @rdname featureDetection-matchedFilter
setReplaceMethod("index", "MatchedFilterParam", function(object, value) {
    object@index <- value
    OK <- validObject(object)
    if (is.character(OK))
        stop(OK)
    return(object)
})
