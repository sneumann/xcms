\name{rawMat-methods}
\docType{methods}
\alias{rawMat}
\alias{rawMat-methods}
\alias{rawMat,xcmsRaw-method}
\title{Get a raw data matrix}
\description{
  Returns a matrix with columns for time, m/z, and intensity that
  represents the raw data from a chromatography mass spectrometry experiment.
}
\section{Methods}{
\describe{

  \item{object = "xcmsRaw"}{
    \code{rawMat(object, mzrange = numeric(), rtrange = numeric(),
      scanrange = numeric(), log=FALSE)}
  }
}}
\arguments{
  \item{object}{The container of the raw data}
  \item{mzrange}{Subset by m/z range}
  \item{rtrange}{Subset by retention time range}
  \item{scanrange}{Subset by scan index range}
  \item{log}{Whether to log transform the intensities}
}
\value{
  A numeric matrix with three columns: time, mz and intensity.
}
\seealso{\code{\link{plotRaw}} for plotting the raw intensities}
\author{Michael Lawrence}
\keyword{methods}
