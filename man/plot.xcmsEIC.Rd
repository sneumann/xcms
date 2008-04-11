\name{plot.xcmsEIC}
\docType{methods}
\alias{plot.xcmsEIC}
\alias{plot, plot-methods}
\title{Plot extracted ion chromatograms from multiple files}
\description{
  Batch plot a list of extracted ion chromatograms to the current
  graphics device.
}
\section{Methods}{
\describe{
\item{x = "xcmsEIC"}{
  \code{
plot.xcmsEIC(x, y, groupidx = groupnames(x), sampleidx = sampnames(x), rtrange = x@rtrange, 
             col = rep(1, length(sampleidx)), legtext = NULL, peakint = TRUE, sleep = 0, ...)}
  }

}}

\value{
  A \code{xcmsSet} object.
}
\arguments{
  \item{x}{the \code{xcmsEIC} object}
  \item{y}{optional \code{xcmsSet} object with peak integration data}
  \item{groupidx}{
    either character vector with names or integer vector with
    indicies of peak groups for which to plot EICs
  }
  \item{sampleidx}{
    either character vector with names or integer vector with
    indicies of samples for which to plot EICs
  }
  \item{rtrange}{
    a two column matrix with minimum and maximum retention times
    between which to return EIC data points

    if it has the same number of rows as the number groups in the
    \code{xcmsEIC} object, then \code{sampleidx} is used to subset
    it. otherwise, it is repeated over the length of \code{sampleidx}

    it may also be a single number specifying the time window around
    the peak for which to plot EIC data
  }
  \item{col}{
    color to use for plotting extracted ion chromatograms. if missing
    and \code{y} is specified, colors are taken from
    \code{unclass(sampclass(y))} and the default palette

    if it is the same length as the number groups in the \code{xcmsEIC}
    object, then \code{sampleidx} is used to subset it. otherwise,
    it is repeated over the length of \code{sampleidx}
  }
  \item{legtext}{
    text to use for legend. if \code{NULL} and \code{y} is specified,
    legend text is taken from the sample class information found
    in the \code{xcmsSet}
  }
  \item{peakint}{
    logical, plot integrated peak area with darkened lines (requires
    that \code{y} also be specified)
  }
  \item{sleep}{seconds to pause between plotting EICs}
  \item{...}{other graphical parameters}
}
\author{Colin A. Smith, \email{csmith@scripps.edu}}
\seealso{
  \code{\link{xcmsEIC-class}},
  \code{\link{png}},
  \code{\link{pdf}},
  \code{\link{postscript}},
}
\keyword{methods}
\keyword{hplot}
