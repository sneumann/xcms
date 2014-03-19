\name{group-methods}
\docType{methods}
\alias{group-methods}
\alias{group}
\alias{group,xcmsSet-method}
\title{Group peaks from different samples together}
\description{
  A number of grouping (or alignment) methods exist in XCMS. \code{group}
  is the generic method.
}
\section{Methods}{
\describe{
\item{object = "xcmsSet"}{
  \code{
    group(object, ...)
  }
}
}}
\arguments{
  \item{object}{\code{\link{xcmsSet-class}} object}
  \item{method}{Method to use for grouping. See details.}
  \item{...}{Optional arguments to be passed along}
}
\details{
  Different algorithms can be used by specifying them with the
  \code{method} argument. For example to use the density-based
  approach described by Smith et al (2006) one would use:
  \code{group(object, method="density")}. This is also the default.

  Further arguments given by \code{...} are
  passed through to the function implementing
  the \code{method}.

  A character vector of \emph{nicknames} for the
  algorithms available is returned by
  \code{getOption("BioC")$xcms$group.methods}.
  If the nickname of a method is called "mzClust",
  the help page for that specific method can
  be accessed with \code{?group.mzClust}.
}

\value{
  An \code{xcmsSet} object with peak group assignments and statistics.
}
\seealso{
  \code{\link{group.density}}
  \code{\link{group.mzClust}}
  \code{\link{group.nearest}}
  \code{\link{xcmsSet-class}},
}
\keyword{methods}
