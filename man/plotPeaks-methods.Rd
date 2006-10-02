\name{plotPeaks-methods}
\docType{methods}
%\alias{plotPeaks-methods}
\alias{plotPeaks}
\alias{plotPeaks,xcmsRaw-method}
\title{Plot a grid of a large number of peaks}
\description{
  Plot extracted ion chromatograms for many peaks simultaneously, 
  indicating peak integration start and end points with vertical
  grey lines.
}
\section{Methods}{
\describe{
\item{object = "xcmsRaw"}{
  \code{plotPeaks(object, peaks, figs, width = 200)}
}
}}
\arguments{
  \item{object}{the \code{xcmsRaw} object}
  \item{peaks}{
    matrix with peak information as produced by \code{\link{findPeaks}}
  }
  \item{figs}{
    two-element vector describing the number of rows and the number
    of columns of peaks to plot, if missing then an approximately
    square grid that will fit the number of peaks supplied
  }
  \item{width}{
    width of chromatogram retention time to plot for each peak
  }
}
\details{
  This function is intended to help graphically analyze the results
  of peak picking. It can help estimate the number of false positives
  and improper integration start and end points. Its output is very
  compact and tries to waste as little space as possible. Each plot
  is labeled with rounded m/z and retention time separated by a
  space.
}
\seealso{
  \code{\link{xcmsRaw-class}},
  \code{\link{findPeaks}},
  \code{\link{split.screen}}
}
\keyword{methods}
\keyword{hplot}
