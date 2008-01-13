# The pipeline is a series of protocols that implement stages.

# Common methods

# get the computational identifier of an object
setGeneric("name", function(object, ...) standardGeneric("name"))
# by default, the class name of the instance
setMethod("name", "ANY", function(object) class(object))

# how an object is identified in a user interface
setGeneric("dispName", function(object, ...) standardGeneric("dispName"))
# by default, the name of the object
setMethod("dispName", "ANY", function(object) name(object))

# name of specific method performed by a protocol
setGeneric("methodName", function(object, ...) standardGeneric("methodName"))

# name of type to accept as input
# this is a single class, but it could be a class union
setGeneric("inType", function(object, ...) standardGeneric("inType"))

# name of type to produce
setGeneric("outType", function(object, ...) standardGeneric("outType"))

# parameters controlling protocol behavior
setGeneric("parameters", function(object) standardGeneric("parameters"))

# returns a widget for controlling and viewing this object
setGeneric("widget", function(object, ...) standardGeneric("widget"))

# returns a widget containing an interactive visualization of the
# specified data in the context of the specified protocol
setGeneric("explore", function(object, protocol, ...)
           standardGeneric("explore"))

# perform an operation on a data structure and return the result
setGeneric("perform", function(object, data, ...) standardGeneric("perform"))

# protocol accessors
setGeneric("protocol", function(object, ...) standardGeneric("protocol"))
setGeneric("protocol<-", function(object, ..., value)
           standardGeneric("protocol<-"))

# extract a pipeline from an object
setGeneric("pipeline", function(object, ...) standardGeneric("pipeline"))

# Each protocol performs a role defined by a stage.
# An xcmsStage object is a factory for its implementing protocols.
setClass("xcmsStage", contains = "VIRTUAL")

# A protocol performs a stage in a particular way.
setClass("xcmsProtocol", contains = "VIRTUAL")

# The pipeline is a stage, protocol and a list of protocols.
setClass("xcmsPipeline",
         representation(dispName = "character"),
         contains = "list")

### PIPELINE METHODS

# The display name accessor
setMethod("dispName", "xcmsPipeline", function(object) {
  if (!length(object@dispName))
    object@name
  else object@dispName
})

# Perform all component protocols
setMethod("perform", "xcmsPipeline", function(object, data, ...)
{
  for (proto in object)
    data <- perform(proto, data)
  data
})

setMethod("inType", "xcmsPipeline", function(object) {
  first <- head(object, 1)
  if (length(first))
    inType(first[[1]])
  else NULL
})
setMethod("outType", "xcmsPipeline", function(object) {
  last <- tail(object, 1)
  if (length(last))
    outType(last[[1]])
  else NULL
})

setMethod("parameters", "xcmsPipeline", function(object) {
  lapply(object, parameters)
})

# return each contiguous range of protocols with same in and out types
setGeneric("pipeline", function(object, ...) standardGeneric("pipeline"))
setMethod("pipeline", "xcmsPipeline",
  function(object, intype = "ANY", outtype = "ANY")
{
  inmatch <- sapply(sapply(object, inType), extends, intype)
  outmatch <- sapply(sapply(object, outType), extends, outtype)
  ranges <- diff(c(FALSE, inmatch & outmatch, FALSE))
  rangemat <- cbind(which(ranges > 0), which(ranges < 0) - 1)
  apply(rangemat, 1, function(range) object[range[1]:range[2]])
})

setGeneric("findProtocols", function(object, ...)
  standardGeneric("findProtocols"))
setMethod("findProtocols", "xcmsPipeline",
  function(object, stage, method = character())
{
  which(sapply(object, is, protocolName(stage, method)))
})

setMethod("protocol", "xcmsPipeline",
          function(object, stage, method = character(), ...)
{
  protos <- findProtocols(object, stage, method)
  if (!length(protos))
    NULL
  else object[[protos[1]]]
})

setReplaceMethod("protocol", "xcmsPipeline",
                 function(object, stage, value)
{
  protos <- findProtocols(object, stage)
  if (length(protos))
    object@.Data[[protos[1]]] <- value
  else object@.Data <- c(object, value)
  object
})

setMethod("show", "xcmsPipeline", function(object) {
  cat("A pipeline with", length(object@.Data), "protocol(s).\n\n")
  if (length(object@.Data))
    show(object@.Data)
})

# Stage methods

setMethod("name", "xcmsStage",
          function(object) dequalifyStageName(class(object)))

# the factory method - creates a protocol given a method name
setMethod("protocol", "xcmsStage",
          function(object, method = xcmsMethodDefault(name(object)), ...)
{
  protocol <- NULL
  me <- name(object)
  class <- protocolName(me, method)
  if (extends(class, qualifyProtocolName(me)))
    protocol <- new(class, ...)
  protocol
})

# get an xcmsStage instance
xcmsStage <- function(name) {
  new(qualifyStageName(name))
}
# private: do not export these name manipulation functions
dequalifyStageName <- function(name) {
  decapitalize(sub("^xcmsStage", "", name))
}
qualifyStageName <- function(name) {
  paste("xcmsStage", capitalize(name), sep="")
}

# to circumvent spurious warnings resulting from not using a character
# literal in the call to standardGeneric().
.dyngeneric <- function(name, args = alist(object=, ...=))
{
  as.function(c(args, substitute(standardGeneric(name), list(name=name))))
}

# stage registration
setStage <- function(name, dispname = name,
                     intype = "xcmsSet", outtype = intype,
                     where = topenv(parent.frame()))
{
  name <- decapitalize(name)
  # register this stage as a class in the 'where' environment
  class <- setClass(qualifyStageName(name), contains = "xcmsStage",
                    where = where)
  # create accessors for 'dispname' and 'inType'
  setMethod("dispName", class, function(object) dispname, where = where)
  setMethod("inType", class, function(object) intype, where = where)
  setMethod("outType", class, function(object) outtype, where = where)
  # create the API for performing a method of this stage
  performFunc <- function(object, method = xcmsMethodDefault(name), ...)
    {
      # need to resolve the arguments against the stage.method() function
      generic <- paste(name, decapitalize(method), sep=".")
      call <- as.call(list(as.name(generic), object, ...))
      args <- as.list(match.call(getMethod(generic, intype), call))
      args <- tail(args, -2)
      slots <- names(args) %in% slotNames(protocolName(name, method))
      proto <- do.call("xcmsProtocol", c(list(name, method), args[slots]))
      do.call("perform", c(list(proto, object), args[!slots]))
    }
  setGeneric(name, .dyngeneric(name, formals(performFunc)), where = where)
  setMethod(name, intype, performFunc, where = where)
  # create a base protocol class for this stage
  protoclass <- setClass(qualifyProtocolName(name), contains = "xcmsProtocol",
                 where = where)
  setMethod("inType", protoclass, function(object) intype, where = where)
  setMethod("outType", protoclass, function(object) outtype, where = where)
  # create methods for getting and setting pipeline protocols
  # not sure if this is necessary
  accessor <- paste(name, "Proto", sep="")
  setGeneric(accessor, .dyngeneric(accessor), where = where)
  setMethod(accessor, "xcmsPipeline", function(object, method = character())
            protocol(object, name, method),
            where = where)
  setMethod(accessor, outtype, function(object, method = character())
            protocol(object@pipeline, name, method),
            where = where)
  replacer <- paste(accessor, "<-", sep="")
  setGeneric(replacer, .dyngeneric(replacer, alist(object=,value=)),
             where = where)
  setReplaceMethod(accessor, "xcmsPipeline",
                   function(object, value)
                   {
                     protocol(object, name) <- value
                     object
                   }, where = where)
  name
}

# Protocol methods

# this is a high-level wrapper that may be overriden for customizing
# protocol performance.
setMethod("perform", "xcmsProtocol", function(object,data,...) {
  performDelegate(object, data, ...)
})

# actually invokes the protocol function delegate
# private -- should NOT be exported from namespace
setGeneric("performDelegate",
           function(object, data, ...) standardGeneric("performDelegate"))

# get an instance of the stage to which this protocol belongs
setGeneric("stage", function(object, ...) standardGeneric("stage"))
# for the base classes
setMethod("stage", "xcmsProtocol",
          function(object) xcmsStageForProtocol(class(object)))

setMethod("parameters", "xcmsProtocol", function(object) {
  # simply return slots as a list
  slots <- lapply(slotNames(object),function(slot_name) slot(object, slot_name))
  names(slots) <- slotNames(object)
  slots
})

setMethod("name", "xcmsProtocol", function(object) 
  dequalifyProtocolName(class(object)))

setMethod("show", "xcmsProtocol", function(object)
{
  stage <- stage(object)
  cat("Stage:", dispName(stage), "\n")
  cat("Protocol:", dispName(object), "\n")
  cat("In/out: ", inType(object), "/", outType(object), "\n\n", sep="")
  params <- parameters(object)
  if (length(params)) {
    cat("Parameters:\n")
    show(params)
  }
})

setMethod("pipeline", "xcmsProtocol", function(object)
{
  if ("pipeline" %in% slotNames(object))
    object@pipeline
  else NULL
})

# Get a protocol instance for the given stage and method
# eg xcmsProtocol("findPeaks", "matchedFilter") yields an instance of
# "xcmsProtoFindPeaksMatchedFilter"

xcmsProtocol <- function(stage, method = xcmsMethodDefault(stage), ...)
{
  new(protocolName(stage, method), ...)
}

protocolName <- function(stage, method, qualify = TRUE)
{
  name <- paste(decapitalize(stage), capitalize(method), sep="")
  if (qualify)
    name <- qualifyProtocolName(name)
  name
}

qualifyProtocolName <- function(name)
{
  paste("xcmsProto", capitalize(name), sep = "")
}
dequalifyProtocolName <- function(name)
{
  decapitalize(sub("^xcmsProto", "", name))
}

# FIXME: need to support setting default methods for stages
xcmsMethodDefault <- function(stage)
{
  def <- getOption("BioC")$xcms[[paste(decapitalize(stage), "method", sep=".")]]
  if (is.null(def))
    def <- ""
  def
}

xcmsStageForProtocol <- function(name) {
  if (!extends(name, "xcmsProtocol"))
    stop("Class '", name, "' is not a protocol class")
  ancestors <- names(getClass(name)@contains)
  protos <- sapply(c(name, ancestors), dequalifyProtocolName)
  stages <- sapply(names(getClass("xcmsStage")@subclasses), dequalifyStageName)
  stages <- stages[stages %in% protos]
  if (length(stages) == 1)
    xcmsStage(stages[[1]])
  else if (length(stages) > 1)
    stop("Protocol '", name, "' inherits from multiple stages: ",
         paste("'", stages, "'", sep="", collapse=", "))
  else NULL
}

# Registration of protocols
setProtocol <- function(method, dispname = method, representation = list(),
                        fun, parent, prototype = list(), validity = NULL,
                        where = topenv(parent.frame()))
{ 
  method <- decapitalize(method)
  # resolve ancestors and find stage
  parent <- qualifyProtocolName(parent)
  stage <- xcmsStageForProtocol(parent)
  if (is.null(stage))
    stop("Failed to derive a stage from parent class: '", parent, '"')
  stagename <- name(stage)
  # class name directly computed from 'stage' and 'method'
  class <- protocolName(stagename, method)
  if (dequalifyProtocolName(class) == stagename)
    stop("Protocol name conflicts with existing stage name '", stagename, "'")
  contains <- parent
  if (missing(fun)) # no function, not pipeline, protocol is abstract
    contains <- c(contains, "VIRTUAL")
  
  # Transform representation to allow language objects (delayed evaluation)
  representation <- lapply(representation, function(cl) {
    union <- paste(cl, "language", sep="OR")
    if (!isClassUnion(union))
      setClassUnion(union, c(cl, "language"))
    union
  })
  
  # add function formals to prototype
  if (!missing(fun)) {
    slots <- c(slotNames(parent), names(representation))
    params <- names(formals(fun)) %in% slots
    nonmissing <- params & nchar(sapply(formals(fun), deparse)) > 0
    prototype[names(formals(fun))[nonmissing]] <- formals(fun)[nonmissing]
  }
  # create prototype without forcing argument evaluation
  prototype <- do.call("prototype", prototype, TRUE)
  setClass(class, representation, prototype, contains, validity, where = where)
  if (!missing(dispname))
    setMethod("dispName", class, function(object) dispname, where = where)
  # remember the 'stage' of the protocol
  # creates a new instance since stages can be redefined
  setMethod("stage", class, function(object) xcmsStage(stagename), where=where)
  # remember the method name
  setMethod("methodName", class, function(object) method, where = where)
  if (!missing(fun)) {
    .fun <- fun
    formal <- slotNames(class) %in% names(formals(fun))
    # set a 'performDelegate' method that calls 'fun'
    setMethod("performDelegate", class, function(object, data, ...)
    {
      .data <- data
      slots <- parameters(object)
      result <- as.function(c(slots[formal], quote({
        do.call(.fun, c(list(.data), 
          sapply(names(formals(sys.function())), as.name), list(...)))
      })))()
      if (!is.null(result)) {           # FIXME: need c.xcmsPipeline()
        result@pipeline@.Data <- c(data@pipeline, object)
        #names(result@pipeline)[length(names(result@pipeline))] <- name(object)
      }
      result
    }, where = where)
    # set a method with the same formals
    generic <- paste(stagename, method, sep=".")
    args <- formals(fun)
    parslots <- slots[!(slots %in% names(args))]
    args[parslots] <- attributes(getClass(parent)@prototype)[parslots]
    names(args)[1] <- "object"
    genargs <- c(args, alist(...=))
    setGeneric(generic, .dyngeneric(generic, genargs), where = where)
    .proto <- list(stagename, method)
    .slots <- slots
    setMethod(generic, inType(stage), 
      as.function(c(args, quote(
    {
      mc <- as.list(match.call())[[-1]]
      slots <- names(mc) %in% .slots
      protocol <- do.call("xcmsProtocol", c(.proto, mc[slots]))
      do.call("perform", c(protocol, mc[!slots]))
    }))), where = where)
  }
  dequalifyProtocolName(class)
}
