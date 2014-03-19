\name{absent-methods}
\docType{methods}
\alias{absent}
\alias{absent,xcmsSet-method}
\alias{present}
\alias{present,xcmsSet-method}
\title{Determine which peaks are absent / present in a sample class}
\description{
  Determine which peaks are absent / present in a sample class
}
\section{Methods}{
\describe{
\item{object = "xcmsSet"}{
  \code{
    absent(object, ...)
    present(object, ...)
  }
}
}}
\arguments{
  \item{object}{\code{\link{xcmsSet-class}} object}
  \item{class}{Name of a sample class from \code{\link{sampclass}}}
  \item{minfrac}{
    minimum fraction of samples necessary in the class to be absent/present
  }
}
\details{
  Determine which peaks are absent / present in a sample class
  The functions treat peaks that are only present because
  of \code{\link{fillPeaks}} correctly, i.e. does not count them as present.
}

\value{
  An logical vector with the same length as \code{nrow(groups(object))}.
}
\seealso{
  \code{\link{group}}
  \code{\link{diffreport}}
}
\keyword{methods}
