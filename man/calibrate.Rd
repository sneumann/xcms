\name{calibrate-methods}
\docType{methods}
\alias{calibrate}
\alias{calibrate,xcmsSet-method}
\title{Calibrate peaks via a set of known masses}
\description{
 Calibrate peaks 
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
  \item{object}{the \code{xcmsSet} object}
    \item{wishlist}{}
    \item{method}{}
    \item{mzppm}{the relative error used for matching peaks in ppm (parts per million)}
    \item{mzabs}{the absolute error used for matching peaks in Da}
    \item{neighbours}{}
    \item{plotres}{}
}
\value{
  An \code{xcmsSet} object with calibrated peaks.
}
\seealso{
  \code{\link{xcmsSet-class}},
}
\keyword{methods}
\keyword{file}
