\name{AutoLockMass-methods}
\docType{methods}
\alias{AutoLockMass}
\alias{AutoLockMass-methods}
\alias{AutoLockMass,xcmsRaw-method}
\title{ Automatic parameter for Lock mass fixing \code{AutoLockMass}  ~~}
\description{
 \code{AutoLockMass}  - This function decides where the lock mass scans are
	in the xcmsRaw object. This is done by using the scan time differences.
}
\section{Methods}{
\describe{
	\item{object = "xcmsRaw"}{
		\code{signature(object = "xcmsRaw")}}
	}
}
\arguments{
  \item{object}{An \code{\link{xcmsRaw-class}} object}
}

\value{
	\code{AutoLockMass} A numeric vector of scan locations corresponding to lock Mass scans
}
\author{Paul Benton, \email{hpaul.benton08@imperial.ac.uk}}
\examples{
\dontrun{library(xcms)
    library(faahKO)
    ## These files do not have this problem
    ## to correct for but just for an example
    cdfpath <- system.file("cdf", package = "faahKO")
    cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)
    xr<-xcmsRaw(cdffiles[1])
    xr
    ##Lets assume that the lockmass starts at 1 and is every 100 scans
    lockMass<-xcms:::makeacqNum(xr, freq=100, start=1)
    ## these are equalvent
    lockmass2<-AutoLockMass(xr)
    all((lockmass == lockmass2) == TRUE)

    ob<-stitch(xr, lockMass)
}
}

\keyword{methods}
\keyword{manip}
\keyword{lockmass}
