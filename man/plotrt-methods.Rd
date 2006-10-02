\name{plotrt-methods}
\docType{methods}
%\alias{plotrt-methods}
\alias{plotrt}
\alias{plotrt,xcmsSet-method}
\title{Plot retention time deviation profiles}
\description{
  Use corrected retention times for each sample to calculate retention
  time deviation profiles and plot each on the same graph.
}
\section{Methods}{
\describe{
\item{object = "xcmsSet"}{
  \code{plotrt(object, col = NULL, ty = NULL, leg = TRUE, 
               densplit = FALSE)}
}
}}
\arguments{
  \item{object}{the \code{xcmsSet} object}
  \item{col}{vector of colors for plotting each sample}
  \item{ty}{vector of line and point types for plotting each sample}
  \item{leg}{logical plot legend with sample labels}
  \item{densplit}{logical, also plot peak overall peak density}
}
\seealso{
  \code{\link{xcmsSet-class}},
  \code{\link{retcor}}
}
\keyword{methods}
\keyword{hplot}
