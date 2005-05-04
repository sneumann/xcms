\name{xcmsRaw}
\alias{xcmsRaw}
\title{Constructor for xcmsRaw objects which reads NetCDF files}
\description{
  This function handles the task of reading a NetCDF file containing
  LC/MS or GC/MS data into a new \code{xcmsRaw} object. It also
  transforms the data into profile (maxrix) mode for efficient
  plotting and data exploration.
}
\usage{
xcmsRaw(cdfname, profstep = 1, profmethod = "intlin", profparam = list(), msmsfile = NA)
}
\arguments{
  \item{cdfname}{path name of the NetCDF file to read}
  \item{profstep}{step size (in m/z) to use for profile generation}
  \item{profmethod}{method to use for profile generation}
  \item{profparam}{extra parameters to use for profile generation}
  \item{msmsfile}{path name of text file with additional MS/MS information}
}
\details{
  If \code{profstep} is set to 0, no profile matrix is generated.
}
\value{
  A \code{xcmsRaw} object.
}
\references{
  NetCDF file format:
  \url{http://my.unidata.ucar.edu/content/software/netcdf/}
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
