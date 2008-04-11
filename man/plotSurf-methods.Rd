\name{plotSurf-methods}
\docType{methods}
%\alias{plotSurf-methods}
\alias{plotSurf}
\alias{plotSurf,xcmsRaw-method}
\title{Plot profile matrix 3D surface using OpenGL}
\description{
  This method uses the rgl package to create interactive three
  dimensonal representations of the profile matrix. It uses the terrain
  color scheme.
}
\section{Methods}{
\describe{
\item{object = "xcmsRaw"}{
  \code{plotSurf(object, log = FALSE, aspect = c(1, 1, .5), ...)}
}
}}
\arguments{
  \item{object}{the \code{xcmsRaw} object}
  \item{log}{logical, log transform intensity}
  \item{aspect}{
    numeric vector with aspect ratio of the m/z, retention time and
    intensity components of the plot
  }
  \item{...}{arguments passed to \code{\link{profRange}}}
}
\details{
  The rgl package is still in development and imposes some limitations
  on the output format. A bug in the axis label code means that the
  axis labels only go from 0 to the aspect ratio constant of that
  axis. Additionally the axes are not labeled with what they are.
  
  It is important to only plot a small portion of the profile matrix.
  Large portions can quickly overwhelm your CPU and memory.
}
\seealso{
  \code{\link{xcmsRaw-class}}
}
\keyword{methods}
\keyword{hplot}
\keyword{iplot}
