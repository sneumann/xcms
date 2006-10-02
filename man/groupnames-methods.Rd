\name{groupnames-methods}
\docType{methods}
%\alias{groupnames-methods}
\alias{groupnames}
\alias{groupnames,xcmsEIC-method}
\alias{groupnames,xcmsSet-method}
\title{Generate unque names for peak groups}
\description{
  Allow linking of peak group data between classes using unique
  group names that remain the same as long as no re-grouping occurs.
}
\section{Methods}{
\describe{
\item{object = "xcmsSet"}{
  \code{(object, mzdec = 0, rtdec = 0, template = NULL)}
}
\item{object = "xcmsEIC"}{
  \code{(object)}
}
}}
\arguments{
  \item{object}{the \code{xcmsSet} or \code{xcmsEIC} object}
  \item{mzdec}{number of decimal places to use for m/z}
  \item{rtdec}{number of decimal places to use for retention time}
  \item{template}{
    a character vector with existing group names whose format should
    be emulated
  }
}
\value{
  A character vector with unique names for each peak group in the
  object. The format is \code{M[m/z]T[time in seconds]}.
}
\seealso{
  \code{\link{xcmsSet-class}},
  \code{\link{xcmsEIC-class}}
}
\keyword{methods}
\keyword{file}
