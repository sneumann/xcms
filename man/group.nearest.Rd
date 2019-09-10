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

\examples{
\dontrun{library(xcms)
    library(faahKO)
    ## These files do not have this problem to correct for
    ## but just for an example
    cdfpath <- system.file("cdf", package = "faahKO")
    cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)

    xset<-xcmsSet(cdffiles)

    gxset<-group(xset, method="nearest")
    nrow(gxset@groups) == 1096 ## the number of features before minFrac

    post.minFrac<-function(object, minFrac=0.5){
        ix.minFrac<-sapply(1:length(unique(sampclass(object))),
                           function(x, object, mf){
                               meta<-groups(object)
                               minFrac.idx<-numeric(length=nrow(meta))
                               idx<-which(
                                   meta[,levels(sampclass(object))[x]] >=
                                   mf*length(which(levels(sampclass(object))[x]
                                                   == sampclass(object)) ))
                               minFrac.idx[idx]<-1
                               return(minFrac.idx)
                           }, object, minFrac)
        ix.minFrac<-as.logical(apply(ix.minFrac, 1, sum))
        ix<-which(ix.minFrac == TRUE)
        return(ix)
    }

    ## using the above function we can get a post processing minFrac
    idx<-post.minFrac(gxset)

    gxset.post<-gxset ## copy the xcmsSet object
    gxset.post@groupidx<-gxset@groupidx[idx]
    gxset.post@groups<-gxset@groups[idx,]

    nrow(gxset.post@groups) == 465 ## number of features after minFrac
}
}


\keyword{methods}
\keyword{file}
