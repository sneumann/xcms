\name{xcmsFragments}
\alias{xcmsFragments}
\title{Constructor for xcmsFragments objects which holds Tandem MS
  peaks}
\description{
  EXPERIMANTAL FEATURE

  xcmsFragments is an object similar to xcmsSet, which holds peaks
  picked (or collected) from one or several xcmsRaw objects.

  There are still discussions going on about the exact API for MS$^n$
  data, so this is likely to change in the future. The code is not yet
  pipeline-ified.

}
\usage{
xcmsFragments(xs, ...)
}
\arguments{
  \item{xs}{A \code{\link{xcmsSet-class}} object which contains
    picked ms1-peaks from one or several experiments}
  \item{\dots}{
    further arguments to the \code{collect} method
  }
}
\details{
  After running collect(xFragments,xSet) The peaktable of the
  xcmsFragments includes the ms1Peaks from all experinemts stored in a
  xcmsSet-object. Further it contains the relevant MSn-peaks from the
  xcmsRaw-objects, which were created temporarily with the paths in
  xcmsSet.
}
\value{
  An \code{xcmsFragments} object.
}
\author{Joachim Kutzera, Steffen Neumann, \email{sneumann@ipb-halle.de}}
\seealso{
  \code{\link{xcmsFragments-class}},
  \code{\link{collect}}
}
\keyword{file}
