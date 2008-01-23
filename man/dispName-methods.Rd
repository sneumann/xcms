\name{dispName-methods}
\docType{methods}
\alias{dispName}
\alias{dispName-methods}
\alias{dispName,ANY-method}
\alias{dispName,xcmsPipeline-method}
\title{Get human-readable identifier}
\description{
  Returns the human-readable name for an object. Typically, this is used
  for display in a user interface.
}
\section{Methods}{
\describe{
  \item{object = "ANY"}{ Returns \code{\link{class}(object)}. }
  \item{object = "xcmsPipeline"}{ Returns value of \code{dispName}
    slot. }

  Every protocol class and stage class has a method, as well.
}}
\author{Michael Lawrence}
\keyword{methods}
\keyword{internal}