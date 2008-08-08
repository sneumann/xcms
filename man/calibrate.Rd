\name{calibrate-methods}
\docType{methods}
\alias{calibrate}
\alias{calibrate,xcmsSet-method}
\title{Calibrate peaks for correcting unprecise m/z values}
\description{
 Calibrate peaks of a xcmsSet via a set of known masses
}
\section{Methods}{
\describe{
\item{object = "xcmsSet"}{
  \code{
    calibrate(object, wishlist,method="linear",
    mzabs=0.0001, mzppm=5,
    neighbours=3, plotres=FALSE)
  }
}
}}
\arguments{
  \item{object}{a \code{xcmsSet} object with uncalibrated mz}
    \item{wishlist}{a vector or a list of vectors with reference m/z-values}
    \item{method}{the used calibrating-method, see below}
    \item{mzppm}{the relative error used for matching peaks in ppm (parts per million)}
    \item{mzabs}{the absolute error used for matching peaks in Da}
    \item{neighbours}{the number of neighbours from wich the one with the highest intensity is used (instead of the nearest)}
    \item{plotres}{can be set to TRUE if wanted a result-plot showing the found m/z with the distances and the regression}
}
\value{
    \item{object}{a \code{xcmsSet} with one ore more samples}
\item{wishlist}{for each sample can be used a different wishlist if there is given a list of m/z-vectors. the length of the list must be the same as the number of samples, alternatively a single vector of masses can be given which is used for all samples.}
   \item{method}{ "shift" for shifting each m/z, "linear" does a linear regression and adds a linear term to each m/z. "edgeshift" does a linear regression within the range of the mz-wishlist and a shift outside.}
}
\seealso{
  \code{\link{xcmsSet-class}},
}
\keyword{methods}
\keyword{file}
