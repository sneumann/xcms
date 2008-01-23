\name{protocol-methods}
\docType{methods}
\alias{protocol}
\alias{protocol-methods}
\alias{protocol,xcmsPipeline-method}
\alias{protocol,xcmsStage-method}
\title{ Get a protocol object. }
\description{
  Retrieves a protocol from the given object.
}
\section{Methods}{
\describe{

  \item{object = "xcmsPipeline"}{
    \code{protocol(object, role, method = character())} - Returns the
    first protocol in the pipeline with the given \code{role} and,
    optionally, \code{method}.
  }

  \item{object = "xcmsStage"}{
    \code{protocol(object, method = defaultMethod(object))} -
    Creates a protocol for the stage that implements the specified
    \code{method}.
  }
}}
\arguments
{
  \item{role}{The role represented by a stage, like "findPeaks".}
  \item{method}{The method implemented by a protocol, like "matchedFilter".}
}
\value
{
  An instance of \code{\linkS4class{xcmsProtocol}} or \code{NULL} if
  none matches the specified criteria.
}
\author{Michael Lawrence}
\seealso{\code{\link{protocol<-}} for setting (instead of getting) a protocol.}
\keyword{methods}
\keyword{internal}