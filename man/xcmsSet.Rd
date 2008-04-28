\name{xcmsSet}
\alias{xcmsSet}
\title{Constructor for xcmsSet objects which finds peaks in NetCDF/mzXML files}
\description{
  This function handles the construction of xcmsSet objects. It
  finds peaks in batch mode and pre-sorts files from subdirectories
  into different classes suitable for grouping.
}
\usage{
xcmsSet(files = NULL, snames = NULL, sclass = NULL,
        profmethod = "bin", profparam = list(), ..., rawpipeline = NULL,
        pipeline = NULL)
}
\arguments{
  \item{files}{path names of the NetCDF/mzXML files to read}
  \item{snames}{sample names}
  \item{sclass}{sample classes}
  \item{profmethod}{method to use for profile generation}
  \item{profparam}{parameters to use for profile generation}
  \item{\dots}{
    further arguments to the \code{findPeaks} method of the
    \code{xcmsRaw} class
  }
  \item{rawpipeline}{
    A pipeline that transforms an \code{\linkS4class{xcmsRaw}} to
    an \code{\linkS4class{xcmsPeaks}} (as in
    \code{\link{processRaws}}). Overrides any profile parameters
    (\code{profmethod} and \code{profparam}) and peak parameters in \code{...}. 
  }
  \item{pipeline}{
    A pipeline that transforms the nascent \code{xcmsSet} (i.e. one
    without any detected peaks) to another \code{xcmsSet} (i.e. one with
    peaks and possibly additional results). Overrides the
    \code{rawpipeline} parameter, as well as the profile and peak parameters.
  }
}
\details{
  The default values of the \code{files}, \code{snames}, and
  \code{sclass} arguments cause the function to recursively search
  for readable files. The filename without extention is used for the
  sample name. The subdirectory path is used for the sample class.

  In the canonical workflow, one will call this function after having
  developed, presumably through direct use of the
  \code{\linkS4class{xcmsRaw}} class, a pipeline that generates an
  \code{\linkS4class{xcmsPeaks}} from an \code{xcmsRaw}. This pipeline
  may be passed as the the \code{rawpipeline} parameter.

  The \code{pipeline} parameter may be a convenient shortcut when a
  compatible pipeline has been developed during prior analysis of this
  or a very similar dataset.
}
\value{
  A \code{xcmsSet} object.
}
\author{Colin A. Smith, \email{csmith@scripps.edu}}
\note{
  
}
\seealso{
  \code{\link{xcmsSet-class}},
  \code{\link{findPeaks}},
  \code{\link{profStep}},
  \code{\link{profMethod}}
}
\keyword{file}
