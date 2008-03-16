\name{xcmsFragments}
\alias{xcmsFragments}
\title{Constructor for xcmsFragments objects which collects MS$^n$
  spectra from mzXML/mzData files}
\description{
  This function handles the construction of xcmsFragments objects. 
}
\usage{
xcmsFragments(xs, ...)
}
\arguments{
  \item{xs}{A \code{\link{xcmsSet-class}} object which contains
    picked ms1-peaks from several experiments}  
  \item{\dots}{
    further arguments to the \code{collect} method
  }
}
\details{
  After running collect(xFragments,xSet) The peaktable of the
  xcmsFragments includes the ms1Peaks from all experinemts stored in a
  xcmsSet-object. Further it contains the relevant msN-peaks from the
  xcmsRaw-objects, which were created temporarily with the paths in
  xcmsSet.  
}
\value{
  An \code{xcmsFragments} object.
}
\author{S. Neumann, J. Kutzera}
\seealso{
  \code{\link{xcmsFragments-class}},
  \code{\link{collect}}
}
\keyword{file}
