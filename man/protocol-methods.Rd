\name{protocol-methods}
\docType{methods}
\alias{protocol}
\alias{protocol-methods}
\alias{protocol,xcmsPipeline-method}
\alias{protocol,xcmsStage-method}
\alias{protocol<-}
\alias{protocol<--methods}
\alias{protocol<-,xcmsPipeline-method}
\title{ Get/set a protocol object. }
\description{
  Set / retrieves a protocol from the given object.
}
\section{Methods}{
\describe{
  \item{object = "xcmsPipeline"}{
    \code{protocol(object, role, method = character())} - Returns the
    first protocol in the pipeline with the given \code{role} and,
    optionally, \code{method}.

    \code{protocol(object, role) <- value} - Set the protocol that performs the
    specified \code{role} in the pipeline. The first protocol found to
    implement the stage identified by \code{role} is replaced by \code{value}.

  }

  \item{object = "xcmsStage"}{
    \item{object}{The protocol container}
    \code{protocol(object, method = defaultMethod(object))} -
    Creates a protocol for the stage that implements the specified
    \code{method}.
      \item{role}{A character vector identifying the role/stage performed by
    the protocol}

  }
}}
\arguments{
  \item{role}{The role represented by a stage, like "findPeaks".}
  \item{method}{The method implemented by a protocol, like "matchedFilter".}
}
\value{
  An instance of \code{\linkS4class{xcmsProtocol}} or \code{NULL} if
  none matches the specified criteria.
}
\author{Michael Lawrence}
\seealso{\code{\link{protocol<-}} for setting (instead of getting) a protocol.}
\keyword{methods}
\keyword{internal}