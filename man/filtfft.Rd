\name{filtfft}
\alias{filtfft}
\title{Apply an convolution filter using an FFT}
\description{
  Expands a vector to the length of the filter and then convolutes
  it using two successive FFTs.
}
\usage{
filtfft(y, filt)
}
\arguments{
  \item{y}{numeric vector of data to be filtered}
  \item{filt}{filter with length \code{nextn(length(y))}}
}
\value{
  A numeric vector the same length as y.
}
\author{Colin A. Smith, \email{csmith@scripps.edu}}
\keyword{math}
\keyword{internal}
