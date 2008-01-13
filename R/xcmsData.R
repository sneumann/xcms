setClass("xcmsData", representation(pipeline = "xcmsPipeline"))

setMethod("pipeline", "xcmsData", function(object, local = FALSE)
          {
            if (local)
              tail(pipeline(object@pipeline, outtype = class(object)), 1)[[1]]
            else object@pipeline
          })

# explore the data in the context of the last applied protocol
setMethod("explore", c("xcmsData", "missing"), function(object, protocol, ...)
          {
            proto <- NULL
            if (length(object@pipeline@.Data))
              proto <- tail(object@pipeline,1)[[1]]
            explore(proto, object)
          })
