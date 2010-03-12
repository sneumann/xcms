\name{xcmsSet}
\alias{xcmsSet}
\title{Constructor for xcmsSet objects which finds peaks in NetCDF/mzXML files}
\description{
  This function handles the construction of xcmsSet objects. It
  finds peaks in batch mode and pre-sorts files from subdirectories
  into different classes suitable for grouping.
}
\usage{
xcmsSet(files = NULL, snames = NULL, sclass = NULL, phenoData = NULL,
        profmethod = "bin", profparam = list(),
        polarity = NULL, mslevel=NULL,
        nSlaves=0, progressCallback=NULL,...)
}
\arguments{
  \item{files}{path names of the NetCDF/mzXML files to read}
  \item{snames}{sample names}
  \item{sclass}{sample classes}
    \item{phenoData}{sample names and classes}
  \item{profmethod}{method to use for profile generation}
  \item{profparam}{parameters to use for profile generation}
  \item{polarity}{filter raw data for positive/negative scans}
    \item{mslevel}{perform peak picking on data of given mslevel}
    \item{nSlaves}{number of MPI-slaves to use for parallel peak detection,
                 works only if \code{Rmpi} is installed properly, see \code{\link{xcmsPapply}}.}
  \item{progressCallback}{function to be called, when progressInfo changes (useful for GUIs)}
  \item{\dots}{
    further arguments to the \code{findPeaks} method of the
    \code{xcmsRaw} class
  }
}
\details{
  The default values of the \code{files}, \code{snames}, \code{sclass}, and
  \code{phenoData} arguments cause the function to recursively search
  for readable files. The filename without extention is used for the
  sample name. The subdirectory path is used for the sample class.
  If the files contain both positive and negative spectra, the polarity
  can be selected explicitely. The default (NULL) is to read all scans.
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
