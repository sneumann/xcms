\name{xcmsSet}
\alias{xcmsSet}
\title{Constructor for xcmsSet objects which finds peaks in NetCDF/mzXML files}
\description{
  This function handles the construction of xcmsSet objects. It
  finds peaks in batch mode and pre-sorts files from subdirectories
  into different classes suitable for grouping.
}
\usage{
xcmsSet(files = list.files(pattern = "\.[Cc][Dd][Ff]$|\.[Nn][Cc]$|\.([Mm][Zz])?[Xx][Mm][Ll]$", 
                           recursive = TRUE), 
        snames = gsub("\.[^.]*$", "", basename(files)), 
        sclass = gsub("^\.$", "sample", dirname(files)),
        profmethod = "bin", profparam = list(), ...)
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
}
\details{
  The default values of the \code{files}, \code{snames}, and
  \code{sclass} arguments cause the function to recursively search
  for readable files. The filename without extention is used for the
  sample name. The subdirectory path is used for the sample class.
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
