\name{clustunique}
\alias{clustunique}
\title{Create unique matrix of m/z and retention time windows}
\description{
  Find indicies of nonoverlapping m/z and retention time windows.
  THIS FUNCTION HAS BEEN DEPRICATED IN FAVOR OF rectUnique AND WILL
  BE REMOVED IN A FUTURE RELEASE.
}
\usage{
clustunique(clust, priority = seq(length = nrow(clust)), mzdiff = 0, max = 5)
}
\arguments{
  \item{clust}{
    matrix with \code{mzmin}, \code{mzmax}, \code{retmin}, and
    \code{retmax} columns
  }
  \item{priority}{integer vector with indicies ordered by priority}
  \item{mzdiff}{minimal difference in m/z}
  \item{max}{maximum number of sorted rows to check for uniqueness}
}
\value{
  An integer vector with indicies of rows to keep
}
\author{Colin A. Smith, \email{csmith@scripps.edu}}
\keyword{array}
\keyword{internal}
