\name{xcmsRaw}
\alias{xcmsRaw}
\title{Constructor for xcmsRaw objects which reads NetCDF/mzXML files}
\description{
  This function handles the task of reading a NetCDF/mzXML file containing
  LC/MS or GC/MS data into a new \code{xcmsRaw} object. It also
  transforms the data into profile (maxrix) mode for efficient
  plotting and data exploration.
}
\usage{
xcmsRaw(filename, profstep = 1, profmethod = "intlin", profparam =
list(), includeMSn = FALSE)
}
\arguments{
  \item{filename}{path name of the NetCDF or mzXML file to read}
  \item{profstep}{step size (in m/z) to use for profile generation}
  \item{profmethod}{method to use for profile generation}
  \item{profparam}{extra parameters to use for profile generation}
  \item{includeMSn}{whether to read additional MS$^n$ spectra (mzData
    and mzXML only)}
}
\details{
  If \code{profstep} is set to 0, no profile matrix is generated.
  Only first level MS data is read, not MS/MS, etc.)
}
\value{
  A \code{xcmsRaw} object.
}
\references{
  NetCDF file format:
  \url{http://my.unidata.ucar.edu/content/software/netcdf/}

  mzXML file format:
  \url{http://sashimi.sourceforge.net/software_glossolalia.html}
}
\author{Colin A. Smith, \email{csmith@scripps.edu}}
\note{
  
}
\seealso{
  \code{\link{xcmsRaw-class}},
  \code{\link{profStep}},
  \code{\link{profMethod}}
}
\keyword{file}
