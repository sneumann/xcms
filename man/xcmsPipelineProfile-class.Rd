\name{xcmsPipelineProfile-class}
\docType{class}
\alias{xcmsPipelineProfile-class}

\title{Pipeline for profile matrix generation}
\description{A pipeline with specific logic for the creation and
  filtering of profile matrices.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("xcmsPipelineProfile", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{The list of protocols}
    \item{\code{dispName}:}{The human-readable name}
  }
}
\section{Extends}{
Class \code{"\linkS4class{xcmsPipeline}"}, directly.
Class \code{"\linkS4class{list}"}, by class "xcmsPipeline", distance 2.
Class \code{"\linkS4class{vector}"}, by class "xcmsPipeline", distance 3.
}
\section{Methods}{
  \describe{
    \item{\link{perform}}{
      \code{signature(object = "xcmsPipelineProfile", data = "xcmsRaw")}:
      Applies the pipeline to a dataset.
    }
  }
}
\details{
  As profile matrices tend to be quite large, it is often more efficient
  to generate profile subsets. This class supports the specification of
  subsets (see \code{\link{perform}}) and (in the future) will
  try to minimize the "edge effects" when filtering a subset.
}
\author{Michael Lawrence }
\seealso{
  Calling \code{\link{pipeline}} on an \code{\linkS4class{xcmsProfile}}
  is a common way to obtain an instance of this class.
}
\keyword{classes}
