\name{rawEIC-methods}
\docType{methods}
\alias{rawEIC}
\alias{rawEIC,xcmsRaw-method}
\title{Get extracted ion chromatograms for specified m/z range}
\description{
    Generate extracted ion chromatogram for m/z values of interest. The
    raw data is used in contrast to  \code{\link{getEIC}} which uses
    data from the profile matrix (i.e. values binned along the M/Z
    dimension).
}
\section{Methods}{
\describe{
\item{object = "xcmsRaw"}{
  \code{
    rawEIC(object, mzrange = numeric(), rtrange = numeric(), scanrange = numeric())
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
  A list of :

  \item{scan}{
    scan number
  }
  \item{intensity}{
    added intensity values
  }

}

\author{Ralf Tautenhahn}

\seealso{
  \code{\link{xcmsRaw-class}}
}
\keyword{methods}
