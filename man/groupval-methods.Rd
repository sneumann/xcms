\name{groupval-methods}
\docType{methods}
%\alias{groupval-methods}
\alias{groupval}
\alias{groupval,xcmsSet-method}
\title{Extract a matrix of peak values for each group}
\description{
  Generate a matrix of peak values with rows for every group and
  columns for every sample. The value included in the matrix can
  be any of the columns from the \code{xcmsSet} \code{peaks} slot
  matrix. Collisions where more than one peak from a single sample
  are in the same group get resolved with one of several user-selectable
  methods.
}
\section{Methods}{
\describe{
\item{object = "xcmsSet"}{
  \code{groupval(object, method = c("medret", "maxint"),
                 value = "index", intensity = "into")}
}
}}
\arguments{
  \item{object}{the \code{xcmsSet} object}
  \item{method}{
    conflict resolution method, \code{"medret"} to use the peak
    closest to the median retention time or \code{"maxint"} to use
    the peak with the highest intensity
  }
  \item{value}{
    name of peak column to enter into returned matrix, or \code{"index"}
    for index to the corresponding row in the \code{peaks} slot
    matrix
  }
  \item{intensity}{
    if \code{method == "maxint"}, name of peak column to use for
    intensity
  }
}
\value{
  A matrix with with rows for every group and columns for every
  sample. Missing peaks have \code{NA} values.
}
\seealso{
  \code{\link{xcmsSet-class}}
}
\keyword{methods}
\keyword{manip}
