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
list(), includeMSn=FALSE)
}
\arguments{
  \item{filename}{path name of the NetCDF or mzXML file to read}
  \item{profstep}{step size (in m/z) to use for profile generation}
  \item{profmethod}{method to use for profile generation}
  \item{profparam}{extra parameters to use for profile generation}
  \item{includeMSn}{only for XML file formats: also read MS$^n$
    (Tandem-MS of Ion-/Orbi- Trap spectra)}
}
\details{
  If \code{profstep} is set to 0, no profile matrix is generated.
  Unless includeMSn=TRUE only first level MS data is read, not MS/MS, etc.)
}
\value{
  A \code{xcmsRaw} object.
}
\references{
  NetCDF file format:
  \url{http://my.unidata.ucar.edu/content/software/netcdf/}
  \url{http://www.astm.org/Standards/E2077.htm}
  \url{http://www.astm.org/Standards/E2078.htm}

  mzXML file format:
  \url{http://sashimi.sourceforge.net/software_glossolalia.html}

  PSI-MS working group who developed mzData and mzML file formats:
  \url{http://www.psidev.info/index.php?q=node/80}

  Parser used for XML file formats:
  \url{http://tools.proteomecenter.org/wiki/index.php?title=Software:RAMP}

}
\author{Colin A. Smith, \email{csmith@scripps.edu}}
\note{
  
}
\seealso{
  \code{\link{xcmsRaw-class}},
  \code{\link{profStep}},
  \code{\link{profMethod}}
  \code{\link{xcmsFragments}}
}
\keyword{file}
