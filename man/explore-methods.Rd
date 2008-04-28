\name{explore-methods}
\docType{methods}
\alias{explore}
\alias{explore-methods}
\alias{explore,xcmsData,missing-method}
\alias{explore,xcmsRaw,NULL-method}
\alias{explore,xcmsSet,NULL-method}
\title{ Graphically exploring objects }
\description{
  Launch a GUI for exploring the object using interactive graphics and controls.
}
\section{Methods}{
\describe{

\item{object = "xcmsData", protocol = "missing"}{ Invokes the
  \code{explore} method on the last-performed protocol of this data
  object. If no protocol has been performed, invokes \code{explore} on
  \code{NULL}. Subclasses can hook into the \code{NULL} signature, as
  shown below, to provide an exploration of the raw data. }

\item{object = "xcmsRaw", protocol = "NULL"}{ Provides exploration of
  the raw mass spectral data. }

\item{object = "xcmsSet", protocol = "NULL"}{ Provides exploration of
  each raw dataset in the experiment. }
}}
\arguments{
  \item{object}{The object to be explored.}
  \item{protocol}{The protocol providing the analysis context for the
    exploration.}
}
\author{Michael Lawrence}
\keyword{methods}
