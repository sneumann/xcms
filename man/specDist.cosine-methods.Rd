\name{specDist.cosine}
\docType{methods}
\alias{specDist.cosine}
\alias{specDist.cosine,matrix,matrix-method}
\title{a Distance function based on matching peaks}
\description{This method calculates the distance of two sets of peaks
  using the cosine-distance.}
\section{Methods}{
\describe{
  \item{peakTable1 = "matrix", peakTable2 = "matrix"}{
    \code{
      specDist.cosine(peakTable1, peakTable2, mzabs = 0.001, mzppm = 10,
      mzExp = 0.6, intExp = 3, nPdiff = 2, nPmin = 8,
      symmetric = FALSE)
    }
  }
}}
\details{
The result is the cosine-distance of the product from weighted factors
of mz and intensity from matching peaks in the two peaktables.
The factors are calculated as wFact = mz^mzExp * int^intExp.
if no distance is calculated (for example because no matching
peaks were found) the return-value is NA.
}

\arguments{
	\item{peakTable1}{a Matrix containing at least m/z-values,
	row must be called "mz"}
	\item{peakTable2}{the matrix for the other mz-values}
	\item{mzabs}{maximum absolute deviation for two matching peaks}
	\item{mzppm}{relative deviations in ppm for two matching peaks}
	\item{symmetric}{use symmetric pairwise m/z-matches only, or each match}
	\item{mzExp}{the exponent used for mz}
	\item{intExp}{the exponent used for intensity}
	\item{nPdiff}{the maximum nrow-difference of the two peaktables}
	\item{nPmin}{the minimum absolute sum of peaks from both praktables}
}
\usage{
specDist.cosine(peakTable1, peakTable2, mzabs=0.001, mzppm=10, mzExp=0.6,
                intExp=3, nPdiff=2, nPmin=8, symmetric=FALSE)
}
\author{Joachim Kutzera, \email{jkutzer@ipb-halle.de}}
\keyword{methods}
