\name{plotRaw-methods}
\docType{methods}
%\alias{plotRaw-methods}
\alias{plotRaw}
\alias{plotRaw,xcmsRaw-method}
\title{Scatterplot of raw data points}
\description{
  Produce a scatterplot showing raw data point location in retention
  time and m/z. This plot is more useful for centroided data than
  continuum data.
}
\section{Methods}{
\describe{
\item{object = "xcmsRaw"}{
  \code{plotRaw(object, massrange = numeric(), timerange = numeric(),
                scanrange = numeric())}
}
}}
\arguments{
  \item{object}{the \code{xcmsRaw} object}
  \item{massrange}{
    numeric vector of length >= 2 whose range will be used to select
    the masses to plot
  }
  \item{timerange}{
    numeric vector of length >= 2 whose range will be used to select
    the retention times to plot
  }
  \item{scanrange}{
    numeric vector of length >= 2 whose range will be used to select
    scans to plot
  }
}
\value{
  A matrix with the points plotted.
}
\seealso{
  \code{\link{xcmsRaw-class}}}
}
\keyword{methods}
\keyword{hplot}
