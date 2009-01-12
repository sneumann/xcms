\name{xcmsPeaks-class}
\docType{class}
\alias{xcmsPeaks-class}
\alias{show,xcmsPeaks-method}

\title{A matrix of peaks}
\description{A matrix of peak information. The actual columns depend on
  how it is generated (i.e. the \code{\link{findPeaks}} method).}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("xcmsPeaks", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{The matrix holding the peak information}
  }
}
\section{Extends}{
%Class \code{"\linkS4class{xcmsData}"}, directly.
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
