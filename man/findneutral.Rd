\name{findneutral}
\alias{findneutral}
\alias{findneutral,xcmsFragments-method}
\title{Find neutral losses in xcmsFragment objects}
\description{
  This is a method to find a neutral loss with a ppm window in a xcmsFragment object
}
\usage{
findneutral(object, find, ppmE=25, print=TRUE)
}
\arguments{
  \item{object}{xcmsFragment object type}
  \item{find}{The neutral loss to be found}
  \item{ppmE}{the ppm error window for searching}
  \item{print}{If we should print a nice little report}
}
\details{
  The method searches for a given neutral loss in an xcmsFragment object type given a certain ppm error window. The neutral losses are generated between neighbouring ions. The resulting data frame shows the whole scan in which the neutral loss was found.
}
\value{
  A data frame with the following columns:

\item{PrecursorMz}{The precursor m/z of the neutral losses}
\item{MSnParentPeakID}{An index ID of the location of the precursor peak in the xcmsFragment object}
\item{msLevel}{The level of the found fragment ion}
\item{rt}{the Retention time of the found ion}
\item{mz}{the actual m/z of the found fragment ion}
\item{intensity}{The intensity of the fragment ion}
\item{sample}{Which sample the fragment ion came from}
\item{GroupPeakMSn}{an ID if the peaks were grouped by an xcmsSet grouping}
\item{CollisionEnergy}{The collision energy of the precursor scan}

}
\references{
  H. Paul Benton, D.M. Wong, S.A.Strauger, G. Siuzdak "XC\eqn{MS^2}"
  Analytical Chemistry 2008
}
\seealso{
\code{\link{findMZ}},
}
\examples{
\dontrun{
    library(msdata)
    mzMLpath <- system.file("iontrap", package = "msdata")
    mzMLfiles<-list.files(mzMLpath, pattern = "extracted.mzML",
                            recursive = TRUE, full.names = TRUE)
    xs <- xcmsSet(mzMLfiles, method = "MS1")
    ##takes only one file from the file set
    xfrag <- xcmsFragments(xs)
    found<-findneutral(xfrag, 58.1455, 50)
}

}
\author{H. Paul Benton, \email{hpbenton@scripps.edu}}
\keyword{methods}
