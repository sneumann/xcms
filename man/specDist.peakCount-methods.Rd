\name{specDist.peakCount-methods}
\docType{methods}
\alias{specDist.peakCount-methods}
\alias{specDist.peakCount}
\alias{specDist.peakCount,matrix,matrix-method}

\title{a Distance function based on matching peaks}
\description{This method calculates the distance of two sets of peaks by just returning the number of matching peaks (m/z-values).}
\section{Methods}{
\describe{
  \item{peakTable1 = "matrix", peakTable2 = "matrix"}{
  \code{
    specDist.peakCount(peakTable1, peakTable2, mzppm=10,symmetric=FALSE )
    }}
}}
\arguments{
	\item{peakTable1}{a Matrix containing at least m/z-values,
	row must be called "mz"}
	\item{peakTable2}{the matrix for the other mz-values}
	\item{mzabs}{maximum absolute deviation for two matching peaks}
	\item{mzppm}{relative deviations in ppm for two matching peaks}
	\item{symmetric}{use symmetric pairwise m/z-matches only, or each match}
}
\usage{
specDist.peakCount(peakTable1, peakTable2, mzabs=0.001, mzppm=10, symmetric=FALSE)
}
\author{Joachim Kutzera, \email{jkutzer@ipb-halle.de}}
\keyword{methods}
