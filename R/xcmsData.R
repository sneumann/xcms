setClass("xcmsData", representation(pipeline = "xcmsPipeline"))

setMethod("pipeline", "xcmsData", function(object, local = FALSE)
          {
            pipeline <- object@pipeline
            if (local) {
              me <- sapply(sapply(pipeline, outType), extends, class(object))
              if (any(!me)) {
                pipeline <- object@pipeline
                pipeline@.Data <- tail(pipeline, -tail(which(!me),1))
              }
            }
            pipeline
          })

# Perform all component protocols
setMethod("perform", c("xcmsPipeline", "xcmsData"), function(object, data, ...)
{
  for (proto in object)
    data <- perform(proto, data)
  data
})

# this is a high-level wrapper that may be overriden for customizing
# protocol performance.
setMethod("perform", c("xcmsProtocol", "xcmsData"), function(object,data,...) {
  performDelegate(object, data, ...)
})

# explore the data in the context of the last applied protocol
setMethod("explore", c("xcmsData", "missing"), function(object, protocol, ...)
          {
            proto <- NULL
            if (length(object@pipeline@.Data))
              proto <- tail(object@pipeline,1)[[1]]
            explore(object, proto)
          })
