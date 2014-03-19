setClass("xcmsSource", representation("VIRTUAL"))
setGeneric("loadRaw", function(object, ...) standardGeneric("loadRaw"))
setGeneric("xcmsSource", function(object, ...) standardGeneric("xcmsSource"))

## If given an xcmsSource object, simply return it unchanged
setMethod("xcmsSource", "xcmsSource", function(object) object)
