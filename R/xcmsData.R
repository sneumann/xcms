setClass("xcmsData", representation(pipeline = "xcmsPipeline"))

setMethod("pipeline", "xcmsData", function(object, local = FALSE)
          {
            pipeline <- object@pipeline
            if (local) {
              me <- sapply(sapply(pipeline, outType), is, class(object))
              if (any(!me)) {
                pipeline <- object@pipeline
                pipeline@.Data <- tail(pipeline, -tail(which(!me),1))
              }
            }
            pipeline
          })

# explore the data in the context of the last applied protocol
setMethod("explore", c("xcmsData", "missing"), function(object, protocol, ...)
          {
            proto <- NULL
            if (length(object@pipeline@.Data))
              proto <- tail(object@pipeline,1)[[1]]
            explore(proto, object)
          })
