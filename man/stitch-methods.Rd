\name{stitch-methods}
\docType{methods}
\alias{stitch}
\alias{stitch-methods}
\alias{stitch,xcmsRaw-method}
\alias{stitch.netCDF}
\alias{stitch.xml}
\alias{makeacqNum}
\alias{makeacqNum, xcmsRaw-method}

\title{Correct gaps in data}
\description{
  Fixes gaps in data due to calibration scans or lock mass. Automatically detects file type
  and calls the relevant method. The mzXML file keeps the data the same length in time but
  overwrites the lock mass scans. The netCDF version adds the scans back into the data thereby
  increasing the length of the data and correcting for the unseen gap.
}
\section{Methods}{
	\describe{
	\item{object = "xcmsRaw"}{
		\code{
		stitch(object, lockMass=numeric())
		}
		}
	}
	\describe{
	\item{object = "xcmsRaw"}{
	\code{
	makeacqNum(object, freq=numeric(), start=1)
		}
		}
	}

}

\arguments{
  \item{object}{An \code{\link{xcmsRaw-class}} object}
  \item{lockMass}{A dataframe of locations of the gaps}
  \item{freq}{The intervals of the lock mass scans}
  \item{start}{The starting lock mass scan location, default is 1}
}
\details{
	\code{makeacqNum} takes locates the gap using the starting lock mass scan and it's intervals. This data frame is then used in
	\code{stitch} to correct for the gap caused by the lock mass. Correction works by using scans from either side of the gap to fill it in.
}
\value{
	\code{stitch} A corrected \code{xcmsRaw-class} object
	\code{makeacqNum} A numeric vector of scan locations corresponding to lock Mass scans

}
\author{Paul Benton, \email{hpaul.benton08@imperial.ac.uk}}
\examples{
\dontrun{library(xcms)
    library(faahKO)
    ## These files do not have this problem to correct for but just
    ## for an example
    cdfpath <- system.file("cdf", package = "faahKO")
    cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)
    xr<-xcmsRaw(cdffiles[1])
    xr
    ##Lets assume that the lockmass starts at 1 and is every 100 scans
    lockMass<-xcms:::makeacqNum(xr, freq=100, start=1)
    ## these are equcal
    lockmass<-AutoLockMass(xr)
    ob<-stitch(xr, lockMass)
    ob

    ## plot the old data before correction
    foo<-rawEIC(xr, m=c(200,210), scan=c(80,140))
    plot(foo$scan, foo$intensity, type="h")

    ## plot the new corrected data to see what changed
    foo<-rawEIC(ob, m=c(200,210), scan=c(80,140))
    plot(foo$scan, foo$intensity, type="h")
}
}
\keyword{manip}
\keyword{methods}
