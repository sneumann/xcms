# The pipeline is a series of protocols

# Get a protocol instance for the given type and method
# eg xcmsProtocol("findPeaks", "matchedFilter") yields an instance of
# "xcmsProtoFindPeaksMatchedFilter"
xcmsProtocol <- function(type, method = xcmsMethodDefault(type), ...)
{
  class <- xcmsProtocolClass(type, method)
  new(class, ...)
}
xcmsProtocolClass <- function(type, method = xcmsMethodDefault(type))
{
  paste("xcmsProto", capitalize(type), capitalize(method), sep = "")
}
xcmsMethodDefault <- function(type)
{
  def <- getOption("BioC")$xcms[[paste(uncapitalize(type), "method", sep=".")]]
  if (is.null(def))
    def <- ""
  def
}
xcmsProtocolDefault <- function(type) xcmsProtocol(type, xcmsMethodDefault(type))

setProtocolClass <- function(Class, representation, prototype, ...)
{ 
  # Transform representation to allow language objects (delayed evaluation)
  mc <- tail(as.list(match.call()),-1)
  if (!missing(representation)) {
    mc$representation <- lapply(representation, function(cl) {
      union <- paste(cl, "language", sep="OR")
      if (!isClassUnion(union))
        setClassUnion(union, c(cl, "language"))
      union
    })
  }
  # create prototype without forcing argument evaluation
  if (!missing(prototype) && !inherits(prototype, "classPrototypeDef"))
    mc$prototype <- do.call("prototype", prototype, TRUE)
  do.call("setClass", mc)
}

# Base protocol class

setClass("xcmsProtocol", 
  representation(disptype = "character", dispname = "character", dispdesc = "character"), 
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
  slots <- parameters(object)
  # pass data and any extra args (eg subset specifications) to function
  args <- c(list(object = data), slots, list(...))
  fun_name <- paste(type, method, sep=".")

  # quote literal strings 
  quotedSlots <- ifelse(sapply(slots, is.character), shQuote(slots),slots)

  wrapper <- paste("function(", paste(names(slots), slots, sep="=", collapse=","), 
    ") { ", fun_name, "(", paste(c("obj", names(slots)), collapse=","), ") }", sep = "")
  eval(parse(text=wrapper), list(obj = data))()
  #slots <- lapply(slots, eval, args)
  #do.call(fun_name, c(list(data), slots))
})

# what information is provided by this protocol that may be required by others?
setGeneric("provides", function(object, ...) standardGeneric("provides"))

# what information does this protocol require that may be provided by others?
setGeneric("requires", function(object, ...) standardGeneric("requires"))

# returns a widget for controlling and viewing this object
setGeneric("widget", function(object, ...) standardGeneric("widget"))

# returns a widget containing an interactive visualization of the
# specified input and output in the context of this protocol
setGeneric("explore", function(object, ...) standardGeneric("explore"))

setGeneric("parameters", function(object) standardGeneric("parameters"))
setMethod("parameters", "xcmsProtocol", function(object) {
  slots <- lapply(slotNames(object), function(slot_name) slot(object, slot_name))
  names(slots) <- slotNames(object)
  # leave out base slots
  slots[!(names(slots) %in% slotNames("xcmsProtocol"))]
})

setGeneric("dispType", function(object) standardGeneric("dispType"))
setMethod("dispType", "xcmsProtocol", function(object) object@disptype)

setGeneric("dispName", function(object) standardGeneric("dispName"))
setMethod("dispName", "xcmsProtocol", function(object) object@dispname)

setGeneric("dispDesc", function(object) standardGeneric("dispDesc"))
setMethod("dispDesc", "xcmsProtocol", function(object) object@dispdesc)

setMethod("show", "xcmsProtocol", function(object)
{
  cat("Protocol [", dispType(object), ": ", dispName(object), "]\n\n", sep="")
  if (length(dispDesc(object)))
    cat(dispDesc(object), "\n")
  params <- parameters(object)
  if (length(params)) {
    cat("Parameters:\n")
    show(params)
  }
})

# Base profile generation protocol

setProtocolClass("xcmsProtoGenProfile", 
    representation(profstep = "numeric", profmethod = "character", 
        naok = "logical", baselevel = "numeric", basespace = "numeric"),
    list(profstep = 1, profmethod = "bin", naok = TRUE,
      baselevel = quote(0.5*min(object@env$intensity)), basespace = 0.075,
      disptype = "Generate Profile Matrix"),
    "xcmsProtocol")

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

setGeneric("performProfile", function(object, mz, intensity, scanindex, scantime, ...)
    standardGeneric("performProfile"))
# FIXME: only mzrange is currently supported
setMethod("performProfile", "xcmsProtoGenProfile", 
    function(object, mz, intensity, scanindex, scantime, 
        mzindexrange = numeric(), scanrange = numeric(), 
        mzrange = numeric(), rtrange = numeric())
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
    attr(prof, "proftime") <- scantime
    
    prof
})

setMethod("show", "xcmsProtoGenProfile", function(object) {
    
    cat("Profile matrix generation protocol\n")
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

# Base profile filter class

setClass("xcmsProtoFilterProfile", , prototype(disptype = "Profile Matrix Filter"),
  c("xcmsProtocol", "VIRTUAL"))

# Provide necessary margins (as list) in profile matrix for given ranges
# This is to avoid edge effects when processing subsets of the matrix
setGeneric("profMargins", function(object, ...) standardGeneric("profMargins"))

setClass("xcmsProtoFilterProfileSubtract", 
  representation(filter = "xcmsProtoFilterProfile"),
  contains = "xcmsProtoFilterProfile")

setMethod("perform", "xcmsProtoFilterProfileSubtract", 
  function(object, data, ...)
{
  data - perform(object@filter, data, ...)
})

setMethod("dispType", "xcmsProtoFilterProfileSubtract", function(object)
  paste("Subtract", dispType(object@filter)))

setMethod("dispName", "xcmsProtoFilterProfileSubtract", function(object)
  dispName(object@filter))

setMethod("dispDesc", "xcmsProtoFilterProfileSubtract", function(object)
  dispDesc(object@filter))

setMethod("show", "xcmsProtoFilterProfileSubtract", function(object)
{
  cat("Subtraction of ")
  show(object@filter)
})

# Base peak finding class

setClass("xcmsProtoFindPeaks", , prototype(disptype = "Find Peaks"),
  c("xcmsProtocol", "VIRTUAL"))

# Feature-level protocols

setClass("xcmsProtoFindComps", , prototype(disptype = "Find Components"), 
  c("xcmsProtocol", "VIRTUAL"))

setClass("xcmsProtoGroup", , prototype(disptype = "Group Components"), 
  c("xcmsProtocol", "VIRTUAL"))

setClass("xcmsProtoRetcor", , prototype(disptype = "Correct Retention Time"), 
  c("xcmsProtocol", "VIRTUAL"))

setClass("xcmsProtoFillPeaks", , prototype(disptype = "Impute Missing Peaks"),
  c("xcmsProtocol", "VIRTUAL"))
  
setClass("xcmsProtoSummarize", , prototype(disptype = "Summarize Quantities"),
  c("xcmsProtocol", "VIRTUAL"))

setClass("xcmsProtoNorm", , prototype(disptype = "Normalize Quantities"),
  c("xcmsProtocol", "VIRTUAL"))

setClass("xcmsProtoIdent", , prototype(disptype = "Identify Compounds"),
  c("xcmsProtocol", "VIRTUAL"))

# The raw (per-sample) pipeline

setClass("xcmsRawPipeline", 
  representation(rawprotos = "list", genprofproto = "xcmsProtoGenProfile", 
    filtprofprotos = "list"))

setGeneric("genProfProto", function(object) standardGeneric("genProfProto"))

setMethod("genProfProto", "xcmsRawPipeline", function(object) object@genprofproto)

setGeneric("genProfProto<-", function(object, value) standardGeneric("genProfProto<-"))

setReplaceMethod("genProfProto", "xcmsRawPipeline", function(object, value) {
  object@genprofproto <- value
  object
})

setGeneric("filtProfProtos", function(object) standardGeneric("filtProfProtos"))

setMethod("filtProfProtos", "xcmsRawPipeline", function(object) object@filtprofprotos)

setGeneric("filtProfProtos<-", function(object, value) standardGeneric("filtProfProtos<-"))

setReplaceMethod("filtProfProtos", "xcmsRawPipeline", function(object, value) {
  object@filtprofprotos <- value
  object
})

setGeneric("addFiltProfProtos", function(object, value) standardGeneric("addFiltProfProtos"))

setMethod("addFiltProfProtos", "xcmsRawPipeline", function(object, value) {
  filtProfProtos(object) <- c(filtProfProtos(object), value)
  object
})

setMethod("show", "xcmsRawPipeline", function(object) {
  cat("A raw data pipeline containing:\n\n")
  if (length(object@genprofproto))
    show(object@genprofproto)
  if (length(object@filtprofprotos)) {
    cat("\nProfile Processing Protocol(s):\n\n")
    for (proto in object@filtprofprotos)
      show(proto)
  }
})

# The top-level pipelines

setClass("xcmsPipeline", 
  representation(rawpipeline = "xcmsRawPipeline", 
    findpeaksproto = "xcmsProtoFindPeaks", featureprotos = "list"),
  prototype(findpeaksproto = NULL))

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

setMethod("show", "xcmsPipeline", function(object) {
  cat("A pipeline containing:\n\n----\n")
  if (length(object@rawpipeline))
    show(object@rawpipeline)
  cat("----\n")
  if (length(object@findpeaksproto))
    show(object@findpeaksproto)
  cat("----\n")
  if (length(object@featureprotos)) {
    cat("\nFeature Processing Protocol(s):\n\n")
    for (proto in object@featureprotos)
      show(proto)
  }
})
