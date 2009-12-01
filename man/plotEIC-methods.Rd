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
    plotEIC(object, mzrange = numeric(), rtrange = numeric(), scanrange = numeric())
  }
}
}}
\arguments{
  \item{object}{\code{xcmsRaw} object}
  \item{mzrange}{
    m/z range for EIC
  }
  \item{rtrange}{
    retention time range for EIC
  }  
  \item{scanrange}{
    scan range for EIC
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
