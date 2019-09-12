\name{xcmsRaw}
\alias{xcmsRaw}
\alias{deepCopy}
\alias{deepCopy,xcmsRaw-method}
\title{Constructor for xcmsRaw objects which reads NetCDF/mzXML files}
\description{
  This function handles the task of reading a NetCDF/mzXML file containing
  LC/MS or GC/MS data into a new \code{xcmsRaw} object. It also
  transforms the data into profile (maxrix) mode for efficient
  plotting and data exploration.
}
\usage{
xcmsRaw(filename, profstep = 1, profmethod = "bin", profparam =
list(), includeMSn=FALSE, mslevel=NULL, scanrange=NULL)

deepCopy(object)
}
\arguments{
  \item{filename}{path name of the NetCDF or mzXML file to read}
  \item{profstep}{step size (in m/z) to use for profile generation}
  \item{profmethod}{method to use for profile generation. See
    \code{\link{profile-matrix}} for details and supported values.}
  \item{profparam}{extra parameters to use for profile generation}
  \item{includeMSn}{only for XML file formats: also read MS$^n$
    (Tandem-MS of Ion-/Orbi- Trap spectra)}
  \item{mslevel}{move data from mslevel into normal MS1 slots, e.g. for
    peak picking and visualisation}
  \item{scanrange}{scan range to read}
  \item{object}{An xcmsRaw object}
}
\details{
  See \code{\link{profile-matrix}} for details on profile matrix
  generation methods and settings.

  The scanrange to import can be restricted, otherwise all MS1 data
  is read. If \code{profstep} is set to 0, no profile matrix is generated.
  Unless \code{includeMSn = TRUE} only first level MS data is read, not MS/MS,
  etc.

  deepCopy(xraw) will create a copy of the xcmsRaw object with its own
  copy of mz and intensity data in \code{xraw@env}.
}
\value{
  A \code{xcmsRaw} object.
}
\references{
  NetCDF file format:
  \url{https://www.unidata.ucar.edu/software/netcdf/}
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
\examples{
\dontrun{
    library(xcms)
    library(faahKO)
    cdfpath <- system.file("cdf", package = "faahKO")
    cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)
    xr<-xcmsRaw(cdffiles[1])
    xr
    ##This gives some information about the file
    names(attributes(xr))
    ## Lets have a look at the structure of the object

    str(xr)
    ##same but with a preview of each slot in the object
    ##SO... lets have a look at how this works
    head(xr@scanindex)
    ##[1]    0  429  860 1291 1718 2140
    xr@env$mz[425:430]
    ##[1] 596.3 597.0 597.3 598.1 599.3 200.1
    ##We can see that the 429 index is the last mz of scan 1 therefore...

    mz.scan1<-xr@env$mz[(1+xr@scanindex[1]):xr@scanindex[2]]
    intensity.scan1<-xr@env$intensity[(1+xr@scanindex[1]):xr@scanindex[2]]
    plot(mz.scan1, intensity.scan1, type="h",
         main=paste("Scan 1 of file", basename(cdffiles[1]), sep=""))
    ##the easier way :p
    scan1<-getScan(xr, 1)
    head(scan1)
    plotScan(xr, 1)
	}

}
\seealso{
  \code{\link{xcmsRaw-class}},
  \code{\link{profStep}},
  \code{\link{profMethod}}
  \code{\link{xcmsFragments}}
}
\keyword{file}
