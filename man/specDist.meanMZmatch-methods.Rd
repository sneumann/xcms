\name{specDist.meanMZmatch}
\docType{methods}
\alias{specDist.meanMZmatch}
\alias{specDist.meanMZmatch,matrix,matrix-method}
\title{a Distance function based on matching peaks}
\description{This method calculates the distance of two sets of peaks.}
\section{Methods}{
\describe{
\item{peakTable1 = "matrix", peakTable2 = "matrix"}{
  \code{
    specDist.meanMZmatch(peakTable1, peakTable2,
    matchdist=1, matchrate=1,
    mzabs=0.001, mzppm=10, symmetric=TRUE)
    }
}
}}
\details{
The result of the calculation is a weighted sum of two values.
Value one is the mean absolute difference of the matching peaks,
value two is the relation of matching peaks and non matching peaks.
if no distance is calculated (for example because no matching peaks
were found) the return-value is NA.
}
\arguments{
	\item{peakTable1}{a Matrix containing at least m/z-values,
	row must be called "mz"}
	\item{peakTable2}{the matrix for the other mz-values}
	\item{mzabs}{maximum absolute deviation for two matching peaks}
	\item{mzppm}{relative deviations in ppm for two matching peaks}
	\item{symmetric}{use symmetric pairwise m/z-matches only, or each match}
	\item{matchdist}{the weight for value one (see details) }
	\item{matchrate}{the weight for value two}
}
\usage{
specDist.meanMZmatch(peakTable1, peakTable2, matchdist=1, matchrate=1,
                     mzabs=0.001, mzppm=10, symmetric=TRUE)
}
\author{Joachim Kutzera, \email{jkutzer@ipb-halle.de}}
\keyword{methods}
