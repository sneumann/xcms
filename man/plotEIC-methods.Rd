\name{plotEIC-methods}
\docType{methods}
\alias{plotEIC}
\alias{plotEIC,xcmsRaw-method}
\title{Plot extracted ion chromatograms for specified m/z range}
\description{
    Plot extracted ion chromatogram for m/z values of interest. The raw data is used in contrast to  \code{\link{plotChrom}} which uses data from the profile matrix.
}
\section{Methods}{
\describe{
\item{object = "xcmsRaw"}{
  \code{
    plotEIC(object, mzrange = numeric(), rtrange = numeric(),
    scanrange = numeric(), mzdec=2, type="l", add=FALSE, ...)
  }
}
}}
\arguments{
  \item{object}{\code{xcmsRaw} object}
  \item{mzrange}{
    m/z range for EIC. Uses the full m/z range by default.
  }
  \item{rtrange}{
    retention time range for EIC. Uses the full retention time range by default.
  }
  \item{scanrange}{
    scan range for EIC
  }
  \item{mzdec}{
    Number of decimal places of title m/z values in the eic plot.
  }
  \item{type}{
    Speficies how the data should be plotted (by default as a line).
  }
  \item{add}{
    If the EIC should be added to an existing plot.
  }
  \item{...}{
    Additional parameters passed to the plotting function
    (e.g. \code{col} etc).
  }
}
\value{
 A two-column matrix with the plotted points.
}

\author{Ralf Tautenhahn}

\seealso{
  \code{\link{rawEIC},\link{xcmsRaw-class}}
}
\keyword{methods}
