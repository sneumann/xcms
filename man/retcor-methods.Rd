\name{retcor-methods}
\docType{methods}
%\alias{retcor-methods}
\alias{retcor}
\alias{retcor,xcmsSet-method}
\title{Align retention times across samples}
\description{
  Use ``well behaved'' peak groups to calculate retention time
  deviations for every time point of each sample. Use smoothed
  deviations to align retention times.
}
\section{Methods}{
\describe{
\item{object = "xcmsSet"}{
  \code{retcor(object, missing = 1, extra = 1,
               method = c("loess", "linear"), span = .2,
               family = c("gaussian", "symmetric"),
               plottype = c("none", "deviation", "mdevden"),
               col = NULL, ty = NULL)}
}
}}
\arguments{
  \item{object}{the \code{xcmsSet} object}
  \item{missing}{
    number of missing samples to allow in retention time correction
    groups
  }
  \item{extra}{
    number of extra peaks to allow in retention time correction
    correction groups
  }
  \item{method}{
    either \code{"loess"} for non-linear alignment or \code{"linear"}
    for linear alignment
  }
  \item{span}{
    degree of smoothing for local polynomial regression fitting
  }
  \item{family}{
    if \code{gaussian} fitting is by least-squares with no outlier
    removal, and if \code{symmetric} a re-descending M estimator
    is used with Tukey's biweight function, allowing outlier removal
  }
  \item{plottype}{
    if \code{deviation} plot retention time deviation points and
    regression fit, and if \code{mdevden} also plot peak overall
    peak density and retention time correction peak density
  }
  \item{col}{vector of colors for plotting each sample}
  \item{ty}{vector of line and point types for plotting each sample}
}
\value{
  An \code{xcmsSet} object
}
\seealso{
  \code{\link{xcmsSet-class}},
  \code{\link{loess}}
}
\keyword{methods}
\keyword{file}
