\name{xcmsData-class}
\docType{class}
\alias{xcmsData-class}

\title{Base class for datasets}
\description{A dataset analyzed through
  \code{\linkS4class{xcmsProtocol}} objects. The history of the analysis
  is stored as an \code{\linkS4class{xcmsPipeline}} in the
  \code{pipeline} slot. The class also provides a convenient means to
  \code{\link{explore}} analysis results.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("xcmsData", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{pipeline}:}{The pipeline that represents the history of
    the analysis of this dataset. This pipeline can be extracted and
    applied to other datasets, providing a convenient means for
    consistent data analysis.}
  }
}
\section{Methods}{
  \describe{
    \item{\link{explore}}{
      \code{signature(object = "xcmsData", protocol = "missing")}:
      Explore the data in the context of the analysis that produced
      it, if any.
    }
    \item{\link{perform}}{
      \code{signature(object = "xcmsPipeline", data = "xcmsData")}:
      Apply a pipeline on this dataset, producing a new one.
    }
    \item{\link{perform}}{
      \code{signature(object = "xcmsProtocol", data = "xcmsData")}:
      Apply a single protocol on the dataset, producing a new one.
    }
    \item{\link{pipeline}}{
      \code{signature(object = "xcmsData")}:
      Retrieve the \code{pipeline} slot (the analysis history).
    }
  }
}
\details{
  This class supports a step-wise data analysis strategy where the
  results of each step are visually inspected before moving on to the next
  step. Calling the \code{\link{perform}} method on an \code{xcmsData}
  performs an analysis as specified by the provided protocol or
  pipeline. The result is another \code{xcmsData}. Calling the
  \code{\link{explore}} method on the result instance displays an
  interactive exploratory visualization, the nature of which depends on
  the last applied protocol. The analyst may use the visualization to
  ensure that the results at least appear valid before moving on to the
  next step in the analysis. If the \code{explore} method is called
  prior to the application of any protocols, a visualization of the raw
  data is displayed. The type of that visualization depends on the type
  of the data (the subclass of \code{xcmsData}). This is useful for
  initial data checking.
}
\author{ Michael Lawrence }
\seealso{
  Subclasses include: \code{\linkS4class{xcmsRaw}} and
  \code{\linkS4class{xcmsSet}}.
}
\keyword{classes}
