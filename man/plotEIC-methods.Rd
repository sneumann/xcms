\name{plotEIC-methods}
\docType{methods}
%\alias{plotEIC-methods}
\alias{plotEIC}
%\alias{plotEIC,xcmsSet-method}
\title{Plot a list of EICs for each sample}
\description{
  Batch plot a list of extracted ion chromatograms to the current
  graphics device or to files.
}
\section{Methods}{
\describe{
\item{object = "xcmsSet"}{
  \code{plotEIC(object, eics, peakrange,
                nums = seq(length = nrow(peakrange)), 
                classlabel = sampclass(object), peakindex = NULL,
                filebase = character(), wh = c(640,480), 
                sleep = 0)}
}
}}
\arguments{
  \item{object}{the \code{xcmsSet} object}
  \item{eics}{
    list of EIC matricies as produced by \code{link{getEIC}}
  }
  \item{peakrange}{
    matrix with named columns: \code{mzmin}, \code{mzmax}, \code{rtmin},
    and \code{rtmax}
  }
  \item{nums}{
    integer vector of indecies of EICs to plot
  }
  \item{classlabel}{
    factor with class labels for each EIC matrix in the list
  }
  \item{peakindex}{
    integer vector with indecies to the rows of the \code{xcmsSet}
    \code{group} matrix corresponding to the rows of the \code{peakrange}
    matrix
  }
  \item{filebase}{
    base file name under which to save the numbered EICs, no image
    files are saved if missing
  }
  \item{wh}{
    integer vector of length two with pixel width and height of the
    output PNG images
  }
  \item{sleep}{
    seconds to pause between plotting EICs, only useful if
    \code{filebase} is missing and the plots are interactively
    output
  }
}
\seealso{
  \code{\link{xcmsSet-class}}
}
\keyword{methods}
\keyword{hplot}
\keyword{file}
