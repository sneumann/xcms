\name{ramp}
\alias{rampInit}
\alias{rampPrintFiles}
\alias{rampIsFile}
\alias{rampOpen}
\alias{rampClose}
\alias{rampCloseAll}
\alias{rampNumScans}
\alias{rampScanHeaders}
\alias{rampSIPeaks}
\alias{rampRawData}
\title{High-performance, low-level access to mzXML/mzData files}
\description{
  Functions for low-level access to data in mzXML/mzData files.
  They have been optimized for maximum speed and memory efficiency.
}
\usage{
rampInit()
rampPrintFiles()
rampIsFile(filename)
rampOpen(filename)
rampClose(rampid)
rampCloseAll()
rampNumScans(rampid)
rampScanHeaders(rampid)
rampSIPeaks(rampid, seqNum, peaksCount)
rampRawData(rampid)
}
\arguments{
  \item{filename}{path to mxXML/mzData file}
  \item{rampid}{RAMP file ID}
  \item{seqNum}{integer vector with indecies of scans to read}
  \item{peaksCount}{
    integer vector with number of peaks in each scan to read
  }
}
\details{
  The mechanisms and structure of mzXML/mzData files are not discussed
  here, please see the mzXML/mzData documentation for more information.
  These functions make use of the RAMP (Random Access Minimal Parser)
  code from the SASHIMI open-source project.

  \code{rampInit} initializes the data structures used for holding
  RAMP file pointers and indices.

  \code{rampPrintFiles} prints out a list of all files currently
  open.

  \code{rampCloseAll} closes all mzXML/mzData files and frees memory
  associated with their indices.
}
\value{
  \code{rampIsFile} returns \code{TRUE} if the referenced file is
  an mzXML or mzData file and \code{FALSE} otherwise.

  \code{rampOpen} returns the \code{rampid} of the opened file.

  \code{rampClose} returns the closing status.

  \code{rampNumScans} returns the number of scans with an msLevel
  of 1.

  \code{rampSIPeaks} returns a named list with components
  \code{scanindex}, \code{mz}, and \code{intensity}.

  \code{rampScanHeaders} returns a data frame with header information
  for each scan.

  \code{rampRawData} returns a named list with components \code{rt},
  \code{tic}, \code{scanindex}, \code{mz}, and \code{intensity}.

  An error in any of the functions will return a (negative) integer
  error code.
}
\references{
  mzXML file format:
  \url{http://sashimi.sourceforge.net/software_glossolalia.html}
  
  mzData file format:
  \url{http://psidev.sourceforge.net/ms/index.html}
}
\author{Colin A. Smith, \email{csmith@scripps.edu}}
\keyword{file}
\keyword{internal}
