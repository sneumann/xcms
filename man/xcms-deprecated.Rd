\name{xcms-deprecated}
\alias{xcms-deprecated}
\title{Deprecated functions in package \sQuote{xcms}}

\description{
  These functions are provided for compatibility with older versions of
  \sQuote{xcms} only, and will be defunct at the next release.
}

\details{

  The following functions/methods are deprecated.

  \itemize{
    \item{\code{profBin}, \code{profBinM}, \code{profBinLin},
      \code{profBinLinM}, \code{profBinLinBase}, \code{profBinLinBaseM}
      have been deprecated and \code{\link{binYonX}}} in combination
    with \code{\link{imputeLinInterpol}} should be used instead.

    \item{\code{extractMsData}}: replaced by \code{as(x, "data.frame")}.

    \item{\code{plotMsData}}: replaced by \code{plot(x, type = "XIC")}.
  }
}
