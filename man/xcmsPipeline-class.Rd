\name{xcmsPipeline-class}
\docType{class}
\alias{xcmsPipeline-class}
\alias{show,xcmsPipeline-method}

\title{ A sequence of protocols }
\description{ A pipeline is a list of protocols that is labeled with a
  human-readable name. Convenience methods exist for accessing protocols
  and performing them in sequence. }
\section{Objects from the Class}{
  Objects can be created by calls of the form
  \code{new("xcmsPipeline",  ...)}. As a list, the various utilities for
  accessing and manipulating lists in R also apply to a pipeline
  object. For example, one could retrieve the first protocol with
  \code{pipeline[[1]]}.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{ Actual list in which the protocols are stored }
    \item{\code{dispName}:}{ Human-readable name of the pipeline }
  }
}
\section{Extends}{
Class \code{"\linkS4class{list}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
}
\section{Methods}{
  \describe{
    \item{dispName}{\code{signature(object = "xcmsPipeline")}: Returns
      human-readable name, see \code{\link{dispName}}.}
    \item{findProtocols}{\code{signature(object = "xcmsPipeline")}: Find
      the indices of protocols that match the specified criteria, see
      \code{\link{findProtocols}}.}
    \item{inType}{\code{signature(object = "xcmsPipeline")}: Returns
      input type of the first protocol, see \code{\link{inType}}.}
    \item{outType}{\code{signature(object = "xcmsPipeline")}: Returns
      output type of the last protocol, see \code{\link{outType}}.}
    \item{perform}{
      \code{signature(object = "xcmsPipeline", data = "xcmsData")}:
      performs each protocol in order, see \code{\link{perform}}.}
    \item{pipeline}{\code{signature(object = "xcmsPipeline")}: Extract a
      subpipeline. See \code{\link{pipeline}}.}
    \item{protocol}{\code{signature(object = "xcmsPipeline")}: Extract a
      single protocol, see \code{\link{protocol}}.}
    \item{protocol<-}{\code{signature(object = "xcmsPipeline")}: Replace
      a single protocol, see \code{\link{protocol<-}}.}
    \item{show}{\code{signature(object = "xcmsPipeline")}: Display a
      textual description of the pipeline and its protocols.}
  }

  There are also methods for retrieving a protocol for a given
  stage. See the documentation for the stage class of interest.
}
\details{
  A pipeline is a high-level data structure representing a series of
  steps in a data analysis. It may serve multiple roles, including:
  \itemize{
    \item Modular representation of an analysis strategy that may be
    passed conveniently between analytical routines.
    \item Container for the history of data analysis, as in
    \code{\linkS4class{xcmsData}}.
    \item Data model in an model-view-controller (MVC)
    framework. Of specific interest are visualizations and GUIs, see
    \code{\link{explore}} and \code{\link{widget}}.
  }
}
\author{ Michael Lawrence }
\seealso{
  \code{\link{pipeline}} for obtaining a pipeline instance from a container.  
}
\keyword{classes}
\keyword{internal}