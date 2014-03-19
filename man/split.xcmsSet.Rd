\name{split.xcmsSet}
\docType{methods}
\alias{split.xcmsSet}
\alias{split, split-methods}
\title{Divide an xcmsSet object}
\description{
  Divides the samples and peaks from a \code{xcmsSet} object into
  a list of multiple objects. Group data is discarded.
}
\section{Methods}{
\describe{
\item{xs = "xcmsSet"}{
  \code{
    split(x, f, drop = TRUE, ...)
  }
}
}}

\usage{

}
\arguments{
  \item{xs}{\code{xcmsSet} object}
  \item{f}{factor such that \code{factor(f)} defines the grouping}
  \item{drop}{
    logical indicating if levels that do not occur should be dropped
    (if 'f' is a 'factor' or a list).
 }
 \item{\dots}{further potential arguments passed to methods.}
}
\value{
  A list of \code{xcmsSet} objects.
}
\author{Colin A. Smith, \email{csmith@scripps.edu}}
\seealso{
  \code{\link{xcmsSet-class}}
}
\keyword{manip}
