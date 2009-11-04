\name{group.nearest}
\docType{methods}
\alias{group.nearest}
\alias{group.nearest,xcmsSet-method}
\title{Group peaks from different samples together}
\description{
 Group peaks together across samples by creating a master peak list and
 assigning corresponding peaks from all samples. It is inspired by the
 alignment algorithm of mzMine. For further details check
 \url{http://mzmine.sourceforge.net/} and

 Katajamaa M, Miettinen J, Oresic M: MZmine: Toolbox for processing and
 visualization of mass spectrometry based molecular profile
 data. Bioinformatics (Oxford, England) 2006, 22:634?636.

 Currently, there is no equivalent to \code{minfrac} or \code{minsamp}.
}
\section{Methods}{
\describe{
\item{object = "xcmsSet"}{
  \code{    
    group(object, mzVsRTbalance=10, mzCheck=0.2, rtCheck=15, kNN=10)
  }
}
}}
\arguments{
  \item{object}{the \code{xcmsSet} object}
  \item{mzVsRTbalance}{Multiplicator for mz value before calculating the
    (euclidean) distance between two peaks.}  
  \item{mzCheck}{Maximum tolerated distance for mz.}
  \item{rtCheck}{Maximum tolerated distance for RT.}
  \item{kNN}{Number of nearest Neighbours to check}  
}
\value{
  An \code{xcmsSet} object with peak group assignments and statistics.
}
\seealso{
  \code{\link{xcmsSet-class}},
  \code{\link{group.density}} and 
  \code{\link{group.mzClust}}
}
\keyword{methods}
\keyword{file}
