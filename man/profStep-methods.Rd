\name{profStep-methods}
\docType{methods}
%\alias{profStep-methods}
%\alias{profStep<--methods}
\alias{profStep}
\alias{profStep<-}
\alias{profStep,xcmsRaw-method}
\alias{profStep<-,xcmsRaw-method}
\title{Get and set m/z step for generating profile data}
\description{
  These methods get and set the m/z step for generating profile
  (matrix) data from raw mass spectral data. Smaller steps yield
  more precision at the cost of greater memory usage.
}
\section{Methods}{
\describe{
\item{object = "xcmsRaw"}{
  \code{profStep(object)}
}
}}
\seealso{
  \code{\link{xcmsRaw-class}},
  \code{\link{profMethod}}
}
\examples{
	\dontrun{
		library(faahKO)
		cdfpath <- system.file("cdf", package = "faahKO")
		cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)
		xset <- xcmsRaw(cdffiles[1])

		xset
		plotSurf(xset, mass=c(200,500))

		profStep(xset)<-0.1 ## decrease the bin size to get better resolution
		plotSurf(xset, mass=c(200, 500))
		##works nicer on high resolution data.
	}
}
\keyword{methods}
\keyword{manip}
