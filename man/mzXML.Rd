\name{mzXML}
\alias{mzXMLIsFile}
\alias{mzXMLOpen}
\alias{mzXMLClose}
\alias{mzXMLNumScans}
\alias{mzXMLRawData}
\title{High-performance, low-level access to mzXML data files}
\description{
  Functions for low-level access to data in mzXML files.  They have
  been optimized for maximum speed and memory efficiency.
}
\usage{
mzXMLIsFile(filename)
mzXMLOpen(filename)
mzXMLClose(mxid)
mzXMLNumScans(mxid)
mzXMLRawData(mxid)
}
\arguments{
  \item{filename}{path to mxXML file}
  \item{mxid}{mzXML file ID}
}
\details{
  The mechanisms and structure of mzXML files are not discussed
  here, please se the mzXML documentation for more information.
}
\value{
  \code{mzXMLIsFile} returns \code{TRUE} if the referenced file
  is an mzXML file and \code{FALSE} otherwise.

  \code{mzXMLOpen} returns the \code{mxid} of the opened file.

  \code{mzXMLClose} returns the closing status.

  \code{mzXMLNumScans} returns the number of scans with an msLevel
  of 1.

  \code{mzXMLRawData} returns a named list with components \code{rt},
  \code{tic}, \code{scanindex}, \code{mz}, and \code{intensity}.

  An error in any of the functions will return a (negative) integer
  error code.
}
\references{
  mzXML file format:
  \url{http://sashimi.sourceforge.net/software_glossolalia.html}
}
\author{Colin A. Smith, \email{csmith@scripps.edu}}
\keyword{file}
\keyword{internal}
