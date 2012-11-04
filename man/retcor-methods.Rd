\name{retcor-methods}
\docType{methods}
\alias{retcor-methods}
\alias{retcor}
\alias{retcor,xcmsSet-method}
\title{Correct retention time from different samples}
\description{
  To correct differences between retention times between different
  samples, a number of of methods exist in XCMS. \code{retcor}
  is the generic method.
}
\section{Methods}{
\describe{
\item{object = "xcmsSet"}{
  \code{
    retcor(object, ...)
  }
}
}}
\arguments{
  \item{object}{\code{\link{xcmsSet-class}} object}
  \item{method}{Method to use for retention time correction. See details.}
  \item{...}{Optional arguments to be passed along}
}
\details{
    Different algorithms can be used by specifying them with the
  \code{method} argument. For example to use the approach described by
  Smith et al (2006) one would use: \code{retcor(object,
  method="loess")}. This is also the default.

  Further arguments given by \code{...} are
  passed through to the function implementing
  the \code{method}.

  A character vector of \emph{nicknames} for the
  algorithms available is returned by
  \code{getOption("BioC")$xcms$retcor.methods}.
  If the nickname of a method is called "loess",
  the help page for that specific method can
  be accessed with \code{?retcor.loess}.
}

\value{
  An \code{xcmsSet} object with corrected retntion times.
}
\seealso{
  \code{\link{retcor.loess}}
  \code{\link{retcor.obiwarp}}
  \code{\link{xcmsSet-class}},
}
\keyword{methods}
