\name{netCDF}
\alias{netCDFStrError}
\alias{netCDFIsFile}
\alias{netCDFOpen}
\alias{netCDFClose}
\alias{netCDFVarID}
\alias{netCDFVarLen}
\alias{netCDFVarDouble}
\alias{netCDFVarInt}
\alias{netCDFMSPoints}
\alias{netCDFRawData}
\title{High-performance, low-level access to NetCDF data files}
\description{
  Functions for low-level access to numerical data in NetCDF files.
  They have been optimized for maximum speed and memory efficiency.
}
\usage{
netCDFStrError(ncerr)
netCDFIsFile(filename)
netCDFOpen(filename)
netCDFClose(ncid)
netCDFVarID(ncid, var)
netCDFVarLen(ncid, var)
netCDFVarDouble(ncid, var)
netCDFVarInt(ncid, var)
netCDFMSPoints(ncid, scanIndex)
netCDFRawData(ncid)
}
\arguments{
  \item{filename}{path to NetCDF file}
  \item{ncid}{NetCDF file ID}
  \item{var}{NetCDF variable ID or variable name}
  \item{scanIndex}{double vector with \code{scan_index} from CDF file}
}
\details{
  The mechanisms and structure of NetCDF files are not discussed
  here, please se the NetCDF documentation for more information.
}
\value{
  \code{netCDFStrError} returns the textual description of an error code.

  \code{netCDFIsFile} returns \code{TRUE} if the referenced file
  is a NetCDF file and \code{FALSE} otherwise.

  \code{netCDFOpen} returns the \code{ncid} of the opened file.

  \code{netCDFClose} returns the closing status.

  \code{netCDFVarID} returns the variable ID.

  \code{netCDFVarLen} returns the variable length.

  \code{netCDFVarDouble} and \code{netCDFVarInt} return either
  double or integer numeric vectors.

  \code{netCDFMSPoints} returns a list with elements \code{massValues}
  and \code{intensityValues}, with the masses for each scan in
  ascending order.

  \code{netCDFRawData} returns a named list with components \code{rt},
  \code{tic}, \code{scanindex}, \code{mz}, and \code{intensity}.

  An error in any of the functions will return an integer error
  code with attribute \code{errortext} containing the error description
}
\references{
  NetCDF file format:
  \url{http://my.unidata.ucar.edu/content/software/netcdf/}
}
\author{Colin A. Smith, \email{csmith@scripps.edu}}
\keyword{file}
\keyword{internal}
