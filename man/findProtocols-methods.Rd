\name{findProtocols-methods}
\docType{methods}
\alias{findProtocols}
\alias{findProtocols-methods}
\alias{findProtocols,xcmsPipeline-method}
\title{ Find protocols in a container }
\description{
  Returns the indices of the protocols matching the
  specified criteria.
}
\section{Methods}{
\describe{

  \item{object = "xcmsPipeline"}{
    \code{findProtocols(object, role, method = character())}:
    Returns the indices of the protocols in the pipeline that play
    \code{role} and, optionally, implement \code{method}.
  }
}}
\arguments{
  \item{object}{The (ordered) container of protocols}
  \item{role}{The string identifier of the \code{\link{role}},
    corresponding to a stage.}
  \item{method}{The string identifier of the analytical
    \code{\link{method}}, corresponding to a specific protocol type.}
}
\value{
  The indices of the hits, as in \code{\link{grep}} and \code{\link{match}}.
}
\author{Michael Lawrence}
\keyword{methods}
\keyword{internal}
