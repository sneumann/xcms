\name{msn2xcmsRaw}
\alias{msn2xcmsRaw}
\title{Copy MSn data in an xcmsRaw to the MS slots}
\description{
  The MS2 and MSn data is stored in separate slots,
  and can not directly be used by e.g. findPeaks().
  \code{msn2xcmsRaw()} will copy the MSn spectra
  into the "normal" \code{xcmsRaw} slots.  
}
\usage{
msn2xcmsRaw(xmsn)
}
\arguments{
  \item{xmsn}{an object of class \code{xcmsRaw} that contains
    spectra read with includeMSn=TRUE}
}
\details{

  The default gap value is determined from the 90th percentile of
  the pair-wise differences between adjacent mass values.
}
\value{
  An xcmsRaw object
}
\author{Steffen Neumann \email{sneumann@ipb-halle.de}}
\examples{
 msnfile <- system.file("microtofq/MSMSpos20_6.mzML", package = "msdata")
 xrmsn <- xcmsRaw(msnfile, includeMSn=TRUE)
 xr <- msn2xcmsRaw(xrmsn)
 p <- findPeaks(xr, method="centWave")
}

\seealso{
  \code{\link{xcmsRaw}},
}
\keyword{manip}
