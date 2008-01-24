\name{xcmsRaw}
\alias{xcmsRaw}
\title{Constructor for xcmsRaw objects which reads NetCDF/mzXML/mzData files}
\description{
  This function handles the task of reading a NetCDF/mzXML/mzData file containing
  LC/MS or GC/MS data into a new \code{xcmsRaw} object. It also
  transforms the data into profile (maxrix) mode for efficient
  plotting and data exploration.
}
\usage{
xcmsRaw(filename, profstep = 1, profmethod = "intlin", profparam =
list(), genprof = TRUE, includeMSn = FALSE, pipeline = NULL)
}
\arguments{
  \item{filename}{path name of the NetCDF or mzXML file to read}
  \item{profstep}{step size (in m/z) to use for profile generation}
  \item{profmethod}{method to use for profile generation}
  \item{profparam}{extra parameters to use for profile generation}
  \item{genprof}{whether to generate profile matrix using the above parameters}
  \item{includeMSn}{whether to read additional MS$^n$ spectra (mzData
    and mzXML only)}
  \item{pipeline}{\code{\linkS4class{xcmsPipeline}} to apply to the
    nascent \code{\linkS4class{xcmsRaw}}. The pipeline should have
    \code{xcmsRaw} as its input and output type. This parameter
    overrides the profile generation parameters above (\code{genprof},
    \code{profparam}, \code{profmethod} and \code{profstep}).}
}
\details{
  In the canonical workflow, the analyst begins by importing a raw
  sample as an \code{\linkS4class{xcmsRaw}} using this function. There
  are various utilities for performing an exploratory analysis on the
  raw data. Several of these utilities, as well as some peak detection
  methods, depend on the generation of a profile matrix. For convenience, it is
  possible to specify profile generation parameters in the call to this
  function.

  Depending on the \code{step} parameter, profile matrices can
  be quite large. Instead of generating the entire matrix and adding it
  to the \code{xcmsRaw} object, the user may desire to generate subsets
  of the profile matrix as a separate object of class
  \code{\linkS4class{xcmsProfile}}, which offers various visualization
  utilities. This is achived by setting \code{genprof} to \code{FALSE}
  in this function and then calling \code{\link{profileMatrix}}. Working
  with smaller pieces allows the analyst to compare
  multiple versions of the matrix, such as before and after background
  subtraction. After the analyst is satisfied with a matrix, the
  \code{\link{pipeline}} may be extracted from the \code{xcmsProfile}
  and used to generate the full matrix in the \code{xcmsRaw} using
  \code{\link{genProfile}}. The pipeline could also be passed to
  algorithms like \code{\link{findPeaks.matchedFilter}} which generate
  the profile matrix incrementally. 

  If \code{profstep} is set to 0, a pipeline generation protocol will be
  applied, but it will not generate a profile matrix. This differs from
  setting the \code{genprof} parameter to \code{FALSE}, which does not apply a
  profile generation protocol and thus does not save the profile
  generation parameters. Setting \code{profstep} to 0 avoids generating
  the entire (space consuming) matrix but still provides parameters to
  algorithms that generate and process the profile matrix incrementally.

  To repeat and check the same analysis with a different raw sample, the
  analyst may extract the pipeline from the \code{xcmsRaw} object and pass it
  to this function when importing the new sample.
}
\value{
  A \code{xcmsRaw} object.
}
\note{
  Only first level MS data is read, not MS/MS, etc.
}
\references{
  NetCDF file format:
  \url{http://my.unidata.ucar.edu/content/software/netcdf/}

  mzXML file format:
  \url{http://sashimi.sourceforge.net/software_glossolalia.html}
}
\author{Colin A. Smith, \email{csmith@scripps.edu},
  Michael Lawrence \email{mflawren@fhcrc.org}}
\seealso{
  \code{\link{xcmsRaw-class}},
  \code{\link{profStep}},
  \code{\link{profMethod}},
  \code{\link{genProfile}}
}
\keyword{file}
