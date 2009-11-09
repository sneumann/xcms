\name{peakTable-methods}
\docType{methods}
\alias{peakTable}
\alias{peakTable,xcmsSet-method}
\title{Create report of aligned peak intensities}
\description{
  Create a report showing all aligned peaks.
}
\section{Methods}{
\describe{
\item{object = "xcmsSet"}{
  \code{peakTable(object, 
                   filebase = character(), 
                   value=c("into","maxo","intb"))}
}
}}
\arguments{
  \item{object}{the \code{xcmsSet} object}
  \item{filebase}{
    base file name to save report, \code{.tsv} file and \code{_eic}
    will be appended to this name for the tabular report and EIC
    directory, respectively. if blank nothing will be saved
  }
  \item{value}{
    intensity values to be used for the peakTable. \cr
    If \code{value="into"}, integrated peak intensities are used. \cr
    If \code{value="maxo"}, maximum peak intensities are used. \cr
    If \code{value="intb"}, baseline corrected integrated peak intensities are used (only available if peak detection was done by \code{\link{findPeaks.centWave}}).
  }
}
\details{
  
  This method handles creation of summary reports similar to
  \code{\link{diffreport}}. It returns a summary report that can
  optionally be written out to a tab-separated file.

  If a base file name is provided, the report (see Value section)
  will be saved to a tab separated file. 
}

\value{
  A data frame with the following columns:
  
  \item{mzmed}{median m/z of peaks in the group}
  \item{mzmin}{minimum m/z of peaks in the group}
  \item{mzmax}{maximum m/z of peaks in the group}
  \item{rtmed}{median retention time of peaks in the group}
  \item{rtmin}{minimum retention time of peaks in the group}
  \item{rtmax}{maximum retention time of peaks in the group}
  \item{npeaks}{number of peaks assigned to the group}
  \item{Sample Classes}{
    number samples from each sample class represented in the group
  }
  \item{...}{one column for every sample class}
  \item{Sample Names}{integrated intensity value for every sample}
  \item{...}{one column for every sample}
}
\seealso{
  \code{\link{xcmsSet-class}},
}
\keyword{methods}
\keyword{file}
