require(methods) || stop("Couldnt load package methods")

setClass("xcmsFragments", representation(peaks = "matrix",
                                         MS2spec = "list",
                                         specinfo = "matrix"
                                         ##, pipeline = "xcmsRawPipeline"
                                         ),
         prototype(peaks = matrix(nrow = 0, ncol = 6),
                   MS2spec=NULL,
                   specinfo=NULL
                   ##, pipeline = new("xcmsRawPipeline")
                   ))
