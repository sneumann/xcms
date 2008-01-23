\name{xcmsProtocol}
\alias{xcmsProtocol}
\alias{protocolClass}
\title{ Create a protocol }
\description{
  The \code{xcmsProtocol} function is a utility for creating instances
  of \code{\linkS4class{xcmsProtocol}} given a \code{role} and \code{method}.
}
\usage{
xcmsProtocol(role, method = xcmsMethodDefault(role), ...)
protocolClass(role, method)
}

\arguments{
  \item{role}{Identifies the stage with which the protocol is associated.}
  \item{method}{Identifies the analytical method of the protocol,
    distinguishing the protocol from other protocols of the same stage.}
  \item{\dots}{Arguments to pass to \code{\link{new}} when creating the
    instance; e.g., slot (parameter) values.}
}
\value{
  For \code{xcmsProtocol}, an instance of
  \code{\linkS4class{xcmsProtocol}}.
  For \code{protocolClass}, the name of the protocol class. This is
  useful for setting extra methods on a particular protocol class.
}

\author{ Michael Lawrence }
\seealso{ \code{\link{protocol}} for getting a protocol from a container. }
\examples{
xcmsProtocol("findPeaks", "matchedFilter", step = 1)
}
\keyword{classes}
