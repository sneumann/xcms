\name{peakPlots}
\alias{peakPlots}
\alias{peakPlots, xcmsSet-methods}
\title{Feature list plotter}
\description{
  Plots a picked peak feature plot
}
\section{Methods}{
\describe{
\item{object = "xcmsSet"}{
		\code{peakPlots(object, mzrange, rtrange, density=c(TRUE, FALSE), ...)}
}
}}

\arguments{
  \item{object}{Grouped \code{xcmsSet} object}
  \item{mzrange}{m/z range to be used}
  \item{rtrange}{rt range to be used}
  \item{density}{Should a density plot be plotted as well}
  \item{...}{Any other parameters for plotting}
}
\value{
  Plots each peak within a defined window in a different colour
}

\examples{
	\dontrun{
		library(xcms)
		library(faahKO)
		
		cdfpath <- system.file("cdf", package = "faahKO")
		cdffiles <- list.files(cdfpath, recursive = TRUE,full=T)

		xset <- xcmsSet(cdffiles)
		xsg <- group(xset)
		peakPlots(xsg, mz=c(200,210), rt=c(100,200), density=TRUE)
	}
	
}
\author{H. Paul Benton, \email{hpaul.beonton08@imperial.ac.uk}}
\keyword{methods}
