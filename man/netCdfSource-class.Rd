\name{netCdfSource-class}
\Rdversion{1.1}
\docType{class}
\alias{netCdfSource-class}
\alias{initialize,netCdfSource-method}
\alias{loadRaw,netCdfSource-method}

\title{Class \code{"netCdfSource"}, for reading raw data from
  \code{netCDF} files}
\description{
  An implementation of \code{\linkS4class{xcmsSource}} which can load
  from \code{netCDF} files.
}
\section{Objects from the Class}{
  Objects can be created by calling \code{new("netCdfSource", path)}.
  Upon creation, the presence and format of the file at \code{path} is
  checked.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"character"}. File path
      of a file from which to read raw data as the object's data part.}
    \item{\code{cdf}:}{Object of class \code{"integer"}. Internal file
      descriptor of \code{netCDF} file. There is no need to access this
      slot.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{xcmsFileSource}"}, directly.
Class \code{"\linkS4class{character}"}, by class "xcmsFileSource", distance 2.
Class \code{"\linkS4class{xcmsSource}"}, by class "xcmsFileSource", distance 2.
}
\section{Methods}{
  \describe{
    \item{initialize}{
      \code{signature(.Object = "netCdfSource", path = "character")}:
      Create a new \code{netCdfSource} and try to open the file at \code{path}}
    \item{loadRaw}{\code{signature(object = "netCdfSource")}:
      Read and return the data from the file. Value is that of
      \code{\link{loadRaw}}}
  }
}

\author{Daniel Hackney, \email{dan@haxney.org}}
\seealso{
  \code{\linkS4class{xcmsSource}}, \code{\link{loadRaw-methods}}
}
\examples{
\dontrun{
library(xcms)
library(faahKO)
cdfpath <- system.file("cdf", package = "faahKO")
cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)
## nc is a reference to the netCDF file
nc <- new("netCdfSource", cdffiles[1])

rawdata <- loadRaw(nc)
## The raw data has over 1,000 scans, so don't try to print them all out!
summary(rawdata)
}
}
\keyword{classes}
\keyword{file}
\keyword{internal}
