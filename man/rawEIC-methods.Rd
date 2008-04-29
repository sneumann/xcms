\name{rawEIC-methods}
\docType{methods}
\alias{rawEIC}
\alias{rawEIC,xcmsRaw-method}
\title{Get extracted ion chromatograms for specified m/z range}
\description{
    Generate extracted ion chromatogram for m/z values of interest. The raw data is used in contrast to  \code{\link{getEIC}} which uses data from the profile matrix.
}
\section{Methods}{
\describe{
\item{object = "xcmsRaw"}{
  \code{
    rawEIC(object, massrange, scanrange=c(1,length(object@scantime)))
  }
}
}}
\arguments{
  \item{object}{\code{xcmsRaw} object}
  \item{massrange}{
    m/z range for EIC
  }
  \item{scanrange}{
    scan range for EIC
  }
}
\value{
  A list of :
  
  \item{scan}{
    scan number
  }
  \item{intensity}{
    added intensity values
  }

}

\author{Ralf Tautenhahn, \email{rtautenh@ipb-halle.de}}

\seealso{
  \code{\link{xcmsRaw-class}}
}
\keyword{methods}
