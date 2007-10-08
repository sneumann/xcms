# The pipeline is a series of protocols

# Get a protocol instance for the given type and method
# eg xcmsProtocol("findPeaks", "matchedFilter") yields an instance of
# "xcmsProtoFindPeaksMatchedFilter"
xcmsProtocol <- function(type, method = "", ...)
{
  class <- xcmsProtocolClass(type, method)
  new(class, ...)
}
xcmsProtocolClass <- function(type, method = "")
{
  paste("xcmsProto", capitalize(type), capitalize(method), sep = "")
}

setProtocolClass <- function(Class, representation, prototype, ...)
{ 
  # Transform representation to allow language objects (delayed evaluation)
  repr <- lapply(representation, function(cl) {
    union <- paste(cl, "language", sep="OR")
    if (!isClassUnion(union))
      setClassUnion(union, c(cl, "language"))
    union
  })
  # create prototype without forcing argument evaluation
  if (!inherits(prototype, "classPrototypeDef"))
    prototype <- do.call("prototype", prototype, TRUE)
  setClass(Class, repr, prototype, ...)
}

# Base protocol class

setClass("xcmsProtocol", 
  representation(name = "character", desc = "character"), 
  contains = "VIRTUAL")

setGeneric("perform", function(object, data, ...) standardGeneric("perform"))

setMethod("perform", "xcmsProtocol", function(object, data, ...) {
  # by default, protocols should dispatch on their type/method
  # this may not be the best idea in general, since there is a redundancy
  # between the protocol and the underlying method
  cl <- class(object)
  cldef <- getClass(cl)
  supers <- names(cldef@contains)
  type_cl <- supers[which(supers == "xcmsProtocol")-1]
  if (!length(type_cl))
    stop("Cannot find a perform method for class '", cl, "'")
  type <- uncapitalize(sub("xcmsProto", "", type_cl))
  method <- uncapitalize(sub(type_cl, "", cl))
  slots <- lapply(slotNames(cl), function(slot_name) slot(object, slot_name))
  names(slots) <- slotNames(cl)
  # leave out base slots
  slots <- slots[!(names(slots) %in% slotNames("xcmsProtocol"))]
  # pass data and any extra args (eg subset specifications) to function
  args <- c(list(object = data), slots, list(...))
  # ensure all slots evaluated
  slots <- lapply(slots, eval, args)
  do.call(paste(type, method, sep="."), c(list(data), slots))
})

# returns a widget for controlling and viewing this object
setGeneric("widget", function(object, ...) standardGeneric("widget"))

# Base profile generation protocol

setClass("xcmsProtoGenProfile", 
    representation(profstep = "numeric", profmethod = "character", 
        naok = "logical", baselevel = "numeric", basespace = "numeric"),
    prototype = list(profstep = 1, profmethod = "intlin", naok = TRUE),
    contains = "xcmsProtocol")

.profFunctions <- list(intlin = "profIntLinM", binlin = "profBinLinM", 
                       binlinbase = "profBinLinBaseM", bin = "profBinM",
                       maxidx = "profMaxIdxM")

setGeneric("profMethod", function(object) standardGeneric("profMethod"))

setMethod("profMethod", "xcmsProtoGenProfile", function(object) {
    object@profmethod
})

setGeneric("profMethod<-", function(object, value) standardGeneric("profMethod<-"))

setReplaceMethod("profMethod", "xcmsProtoGenProfile", function(object, value) {

    if (! (value %in% names(.profFunctions)))
        stop("Invalid profile method")
    
    object@profmethod <- value
    
    object
})

setGeneric("profStep", function(object) standardGeneric("profStep"))

setMethod("profStep", "xcmsProtoGenProfile", function(object) {

    object@profstep
})

setGeneric("profStep<-", function(object, value) standardGeneric("profStep<-"))

setReplaceMethod("profStep", "xcmsProtoGenProfile", function(object, value) {

    object@profstep <- value
    object
})

setGeneric("baseLevel", function(object) standardGeneric("baseLevel"))

setMethod("baseLevel", "xcmsProtoGenProfile", function(object) object@baselevel)

setGeneric("baseLevel<-", function(object, value) standardGeneric("baseLevel<-"))

setReplaceMethod("baseLevel", "xcmsProtoGenProfile", function(object, value) {

    stopifnot(object@profmethod == "binlinbase")
    object@baselevel <- value
    object
})

setGeneric("baseSpace", function(object) standardGeneric("baseSpace"))

setMethod("baseSpace", "xcmsProtoGenProfile", function(object) object@basespace)

setGeneric("baseSpace<-", function(object, value) standardGeneric("baseSpace<-"))

setReplaceMethod("baseSpace", "xcmsProtoGenProfile", function(object, value) {

    stopifnot(object@profmethod == "binlinbase")
    object@basespace <- value
    object
})

setMethod("perform", c("xcmsProtoGenProfile", "xcmsRaw"), function(object, data, ...)
{
  performRaw(object, data@env$mz, data@env$intensity, data@scanindex, ...)
})

setGeneric("performRaw", function(object, mz, intensity, scanindex, ...)
    standardGeneric("performRaw"))
# FIXME: only mzrange is currently supported
setMethod("performRaw", "xcmsProtoGenProfile", 
    function(object, mz, intensity, scanindex, mzindexrange = numeric(), 
        scanrange = numeric(), mzrange = numeric(), rtrange = numeric())
{
    step <- profStep(object)
    if (!step)
        return(NULL)
    
    if (length(mzrange) == 2) {
        minmass <- mzrange[1]
        maxmass <- mzrange[2]
    } else {
        minmass <- min(mz)
        maxmass <- max(mz)
    }
    
    minmass <- round(minmass/step)*step
    maxmass <- round(maxmass/step)*step
    num <- round((maxmass - minmass)/step) + 1
    
    profFun <- match.fun(.profFunctions[[profMethod(object)]])
    
    # FIXME: Could factor this out to make this more extensible
    params <- list()
    if (length(object@baselevel))
      params <- c(params, baselevel = object@baselevel)
    if (length(object@basespace))
      params <- c(params, basespace = object@basespace)
    
    prof <- profFun(mz, intensity, scanindex, num, minmass, maxmass, 
        object@naok, params)
    
    # FIXME: copies the profile matrix; need to modify this in-place via C
    attr(prof, "profmz") <- minmass+profStep(object)*(0:(dim(prof)[1]-1))
    prof
})

setMethod("show", "xcmsProtoGenProfile", function(object) {
    
    cat("Profile generation protocol of class '", class(object), "'\n", sep="")
    if (length(object@profmethod)) {
        cat("Method:", object@profmethod, "\n")
        if (object@profmethod == "binlinbase") {
            cat("Baselevel:", object@baselevel, "\n")
            cat("Basespace:", object@basespace, "\n")
        }
    }
    if (length(object@profstep))
      cat("Step:", object@profstep, "m/z\n")
})

# Base baseline removal class

setClass("xcmsProtoRemoveBaseline", contains = "xcmsProtocol")

# Provide necessary margins (as list) in profile matrix for given ranges
# This is to avoid edge effects when processing subsets of the matrix
setGeneric("profMargins", function(object, ...) standardGeneric("profMargins"))

# Base peak finding class

setClass("xcmsProtoFindPeaks", contains = c("xcmsProtocol", "VIRTUAL"))
    
# The raw (per-sample) pipeline

setClass("xcmsRawPipeline", 
  representation(rawprotos = "list", genprofproto = "xcmsProtoGenProfile", 
    profprotos = "list"))

setGeneric("genProfProto", function(object) standardGeneric("genProfProto"))

setMethod("genProfProto", "xcmsRawPipeline", function(object) object@genprofproto)

setGeneric("genProfProto<-", function(object, value) standardGeneric("genProfProto<-"))

setReplaceMethod("genProfProto", "xcmsRawPipeline", function(object, value) {
  object@genprofproto <- value
  object
})

setGeneric("profProtos", function(object) standardGeneric("profProtos"))

setMethod("profProtos", "xcmsRawPipeline", function(object) object@profprotos)

setGeneric("profProtos<-", function(object, value) standardGeneric("profProtos<-"))

setReplaceMethod("profProtos", "xcmsRawPipeline", function(object, value) {
  object@profprotos <- value
  object
})

setGeneric("addProfProtos", function(object, value) standardGeneric("addProfProtos"))

setMethod("addProfProtos", "xcmsRawPipeline", function(object, value) {
  profProtos(object) <- c(profProtos(object), value)
  object
})

setMethod("show", "xcmsRawPipeline", function(object) {
  cat("Pipeline with", 
    length(object@genprofproto) + length(object@profprotos), "protocol(s)\n\n")
  if (length(object@genprofproto))
    show(object@genprofproto)
  if (length(object@profprotos)) {
    cat("\nProfile Processing Protocol(s):\n\n")
    for (proto in object@profprotos)
      show(proto)
  }
})

# The top-level pipelines

setClass("xcmsPipeline", 
  representation(rawpipeline = "xcmsRawPipeline", 
    findpeaksproto = "xcmsProtoFindPeaks", featureprotos = "list"))

setGeneric("rawPipeline", function(object) standardGeneric("rawPipeline"))

setMethod("rawPipeline", "xcmsPipeline", function(object) object@rawpipeline)

setGeneric("rawPipeline<-", function(object, value) standardGeneric("rawPipeline<-"))
setReplaceMethod("rawPipeline", "xcmsPipeline", function(object, value) {
  object@rawpipeline <- value
  object
})

setGeneric("findPeaksProto", function(object) standardGeneric("findPeaksProto"))

setMethod("findPeaksProto", "xcmsPipeline", function(object) object@findpeaksproto)

setGeneric("findPeaksProto<-", function(object, value) standardGeneric("findPeaksProto<-"))

setReplaceMethod("findPeaksProto", "xcmsPipeline", function(object, value) {
  object@findpeaksproto <- value
  object
})

setGeneric("featureProtos", function(object) standardGeneric("featureProtos"))

setMethod("featureProtos", "xcmsPipeline", function(object) object@featureprotos)

setGeneric("featureProtos<-", function(object, value) standardGeneric("featureProtos<-"))

setReplaceMethod("featureProtos", "xcmsPipeline", function(object, value) {
  object@featureprotos <- value
  object
})

setGeneric("addFeatureProtos", function(object, value) standardGeneric("addFeatureProtos"))

setMethod("addFeatureProtos", "xcmsPipeline", function(object, value) {
  featureProtos(object) <- c(featureProtos(object), value)
  object
})

# some utilities
capitalize <- function(str) {
  substring(str, 1, 1) <- toupper(substring(str, 1, 1))
  str
}
uncapitalize <- function(str) {
  substring(str, 1, 1) <- tolower(substring(str, 1, 1))
  str
}
