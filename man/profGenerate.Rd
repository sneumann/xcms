\name{profGenerate}
\alias{profBin}
\alias{profBinM}
\alias{profBinLin}
\alias{profBinLinM}
\alias{profBinLinBase}
\alias{profBinLinBaseM}
\alias{profIntLin}
\alias{profIntLinM}
\alias{profMaxIdx}
\alias{profMaxIdxM}
\title{Generation of profile data}
\description{
  Generates profile (binned) data in a given range from an indexed pair
  of vectors.
}
\usage{
profBin(x, y, num, xstart = min(x), xend = max(x), param = list())
profBinM(x, y, zidx, num, xstart = min(x), xend = max(x), NAOK = FALSE,
         param = list())
profBinLin(x, y, num, xstart = min(x), xend = max(x), param = list())
profBinLinM(x, y, zidx, num, xstart = min(x), xend = max(x), NAOK = FALSE,
            param = list())
profBinLinBase(x, y, num, xstart = min(x), xend = max(x), param = list())
profBinLinBaseM(x, y, zidx, num, xstart = min(x), xend = max(x), NAOK = FALSE,
                param = list())
profIntLin(x, y, num, xstart = min(x), xend = max(x), param = list())
profIntLinM(x, y, zidx, num, xstart = min(x), xend = max(x), NAOK = FALSE,
            param = list())
profMaxIdx(x, y, num, xstart = min(x), xend = max(x), param = list())
profMaxIdxM(x, y, zidx, num, xstart = min(x), xend = max(x), NAOK = FALSE,
            param = list())
}
\arguments{
  \item{x}{numeric vector of value positions}
  \item{y}{numeric vector of values to bin}
  \item{zidx}{starting position of each new segment}
  \item{num}{number of equally spaced x bins}
  \item{xstart}{starting x value}
  \item{xend}{ending x value}
  \item{NAOK}{allow NA values (faster)}
  \item{param}{parameters for profile generation}
}
\details{
  These functions take a vector of unequally spaced \code{y} values
  and transform them into either a vector or matrix, depending on
  whether there is an index or not. Each point in the vector or
  matrix represents the data for the point centered at its corresponding
  \code{x} value, plus or minus half the \code{x} step size
  (\code{xend-xstart/(num-1)}).

  The \code{Bin} functions set each matrix or vector value to the
  maximal point that gets binned into it.

  The \code{BinLin} functions do the same except that they linearly
  interpolate values into which nothing was binned.

  The \code{BinLinBase} functions do the same except that they populate
  empty parts of spectra with a base value. They take to two parameters:
  1) \code{baselevel}, the intensity level to fill in for empty parts
  of the spectra. It defaluts to half of the minimum intensity. 2)
  \code{basespace}, the m/z length after which the signal will drop to
  the base level. Linear interpolation will be used between consecuitive
  data points falling within \code{2*basespace} of eachother. It defaluts
  to 0.075.

  The \code{IntLin} functions set each matrix or vector value to
  the integral of the linearly interpolated data from plus to minus
  half the step size.

  The \code{MaxIdx} functions work similarly to the \code{Bin}
  functions execpt that the return the integer index of which x,y
  pair would be placed in a particular cell.
}
\note{
  There are some issues with the \code{profBinLin} method, see
  \url{https://github.com/sneumann/xcms/issues/46} and
  \url{https://github.com/sneumann/xcms/issues/49}. Thus it is suggested
  to use the functions \code{\link{binYonX}} in combination with
  \code{\link{imputeLinInterpol}} instead.
}
\value{
  For \code{prof*}, a numeric vector of length \code{num}.

  For \code{prof*M}, a matrix with dimensions \code{num} by
  \code{length(zidx)}.

  For \code{MaxIdx}, the data type is integer, for all others it
  is double.
}
\examples{
\dontrun{
    library(faahKO)
    cdfpath <- system.file("cdf", package = "faahKO")
    cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)
    xraw <- xcmsRaw(cdffiles[1])

    image(xraw) ## not how with intLin the intensity's blur
    profMethod(xraw) <- "bin"
    image(xraw) ## now with 'bin' there is no blurring good for centroid data
    ##try binlinbase for profile data
}
}
\author{Colin A. Smith, \email{csmith@scripps.edu}}
\keyword{manip}
\keyword{internal}
