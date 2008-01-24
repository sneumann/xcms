\name{xcmsPeaks-class}
\docType{class}
\alias{xcmsPeaks-class}

\title{A matrix of peaks}
\description{This is essentially a matrix of peak information with an
  associated pipeline that describes how the peaks were detected.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("xcmsPeaks", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{The matrix holding the peak information}
    \item{\code{pipeline}:}{The pipeline holding the analysis history}
  }
}
\section{Extends}{
Class \code{"\linkS4class{xcmsData}"}, directly.
Class \code{"\linkS4class{matrix}"}, from data part.
Class \code{"\linkS4class{array}"}, by class "matrix", distance 2.
Class \code{"\linkS4class{structure}"}, by class "matrix", distance 3.
Class \code{"\linkS4class{vector}"}, by class "matrix", distance 4, with explicit coerce.
}
\section{Methods}{
  None yet. Some utilities for working with peak data would be nice.
}
\author{Michael Lawrence}
\seealso{
  \code{\link{findPeaks}} for detecting peaks in an
  \code{\linkS4class{xcmsRaw}}.
}
\keyword{classes}
