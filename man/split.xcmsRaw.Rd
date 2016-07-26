\name{split.xcmsRaw}
\docType{methods}
\alias{split.xcmsRaw}
\title{Divide an xcmsRaw object}
\description{
  Divides the scans from a \code{xcmsRaw} object into
  a list of multiple objects. MS$^n$ data is discarded.
}
\section{Methods}{
\describe{
\item{xr = "xcmsRaw"}{
  \code{
    split(x, f, drop = TRUE, ...)
  }
}
}}

\arguments{
  \item{x}{\code{xcmsRaw} object}
  \item{f}{factor such that \code{factor(f)} defines the scans which go
    into the new \code{xcmsRaw} objects}
  \item{drop}{
    logical indicating if levels that do not occur should be dropped
    (if 'f' is a 'factor' or a list).
 }
 \item{\dots}{further potential arguments passed to methods.}
}
\value{
  A list of \code{xcmsRaw} objects.
}
\author{Steffen Neumann, \email{sneumann(at)ipb-halle.de}}
\seealso{
  \code{\link{xcmsRaw-class}}
}
\keyword{manip}
