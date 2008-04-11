\name{findMZBoxes-methods}
\docType{methods}
\alias{findMZBoxes}
\alias{findMZBoxes,xcmsRaw-method}
\title{Find areas of high peak density in  LC/MS Data}
\description{
  Find regions of interest (ROI) in centroided LC/MS data
  for further consideration by peak-pickers.
}
\section{Methods}{
\describe{
\item{object = "xcmsRaw"}{
  \code{
    findMZBoxes(object, massrange=c(0.0,0.0),
      scanrange=c(1,length(object@scantime)), dev, minEntries)
  }
}
}}
\arguments{
  \item{object}{\code{\link{xcmsRaw-class}} object}
  \item{massrange}{Not implemented.}
  \item{scanrange}{scan range to check for peaks}
  \item{dev}{maximum width of m/z box in ppm}
  \item{minEntries}{minimum number of peaks(centroids) within the range of \code{dev} ppm}
}
\details{
  Areas of high peak density (characterised by having at least \code{minEntries} peaks
  within \code{dev} ppm in consecutive scans) in the LC/MS map are
  located. The regions need not be "real" peaks, this is only a
  filtering step applied before peak picking.  
}
\value{
  A list with entries :

  \item{scan}{
    scan number
  }
  \item{mz}{
    m/z values
  }
  \item{intensity}{
    intensity values
  }
}

\author{Ralf Tautenhahn, \email{rtautenh@ipb-halle.de}}

\seealso{
  \code{\link{findPeaks.centWave}}
  \code{\link{xcmsRaw-class}}
}
\keyword{methods}
