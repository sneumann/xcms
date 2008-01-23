\name{perform-methods}
\docType{methods}
\alias{perform}
\alias{performDelegate}
\alias{perform-methods}
\alias{perform,xcmsPipeline,xcmsData-method}
\alias{perform,xcmsProtocol,xcmsData-method}
\alias{perform,xcmsPipelineProfile,xcmsRaw-method}
\alias{perform,xcmsProtoFilterProfileBase,xcmsProfile-method}
\title{ Perform an analysis strategy  }
\description{
  Perform an analysis strategy on a dataset.
}
\section{Methods}{
\describe{

  \item{object = "xcmsPipeline", data = "xcmsData"}{
    Perform the protocols of a pipeline in order.
  }

  \item{object = "xcmsProtocol", data = "xcmsData"}{
    Perform the protocol on the dataset. This delegates to
    \code{performDelegate}, which has a method for every concrete protocol.
  }

  \item{object = "xcmsPipelineProfile", data = "xcmsRaw"}{
    Special logic for generating and filtering a profile matrix through
    a pipeline.
  }

  \item{object = "xcmsProtoFilterProfileBase", data = "xcmsProfile"}{
    Subtracts the filter result from the profile matrix if the
    \code{subtract} parameter in the protocol is \code{TRUE}.
  }
}}
\arguments{
  \item{object}{The object representing the analytical strategy
    performed on the data.}
  \item{data}{The input data to the analysis}
}
\author{Michael Lawrence}
\keyword{methods}