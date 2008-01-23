\name{protocol<--methods}
\docType{methods}
\alias{protocol<-}
\alias{protocol<--methods}
\alias{protocol<-,xcmsPipeline-method}
\title{ Set a protocol. }
\description{
 Set a protocol into an object.
}
\section{Methods}{
\describe{
  \item{object = "xcmsPipeline"}{
    \code{protocol(object, role) <- value} - Set the protocol that performs the
    specified \code{role} in the pipeline. The first protocol found to
    implement the stage identified by \code{role} is replaced by \code{value}.
  }
}}
\arguments{
  \item{object}{The protocol container}
  \item{role}{A character vector identifying the role/stage performed by
    the protocol}
}
\author{Michael Lawrence}
\seealso{\code{\link{protocol}} for getting (instead of setting) a protocol.}
\keyword{methods}
\keyword{internal}