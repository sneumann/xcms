\name{loadRaw-methods}
\alias{loadRaw}
\title{Read binary data from a source}
\description{
  This function extracts the raw data which will be used an
  \code{\link{xcmsRaw-class}} object. Further processing of data is done
  in the \code{\link{xcmsRaw}} constructor.
}
\section{Methods}{
\describe{
\item{object = "character"}{
  \code{
    loadRaw(object, includeMSn = FALSE, ...)
  }
}
}}
\arguments{
  \item{object}{Specification of a data source such as a file name or
    database query.}
}
\details{
  The implementing methods can decide how to gather the data.
}
\value{
  TODO
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
\author{Daniel Hackney, \email{dan@haxney.org}}
\examples{
        TODO
}
\seealso{
  \code{\link{xcmsRaw-class}},
  \code{\link{profStep}},
  \code{\link{profMethod}}
  \code{\link{xcmsFragments}}
}
\keyword{file}
