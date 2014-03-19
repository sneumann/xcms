\name{write.mzdata-methods}
\docType{methods}
\alias{write.mzdata}
\alias{write.mzdata,xcmsRaw-method}
\title{Save an xcmsRaw object to a file}
\description{
  Write the raw data to a (simple) mzData file.
}
\section{Methods}{
\describe{
\item{object = "xcmsRaw"}{
  \code{write.mzdata(object, filename)}
}
}}
\arguments{
  \item{object}{the \code{xcmsRaw} object}
  \item{filename}{filename (may include full path)
    for the mzData file. Pipes or URLs are not allowed.
  }
}
\details{
  This function will export a given xcmsRaw object to an mzData file.
  The mzData file will contain a <spectrumList> containing the
  <spectrum> with mass and intensity values in 32 bit precision. Other
  formats are currently not supported. Any header information
  (e.g. additional <software> information or <cvParams>) will be lost.
  Currently, also any MSn information will not be stored.
}
\value{
  None.
}
\seealso{
  \code{\link{xcmsRaw-class}},
  \code{\link{xcmsRaw}},
}
\keyword{methods}
\keyword{file}
