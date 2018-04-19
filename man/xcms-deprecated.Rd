\name{xcms-deprecated}
\alias{xcms-deprecated}
\alias{extractChromatograms}
\alias{extractChromatograms,OnDiskMSnExp-method}
\alias{extractChromatograms,XCMSnExp-method}
\alias{plotChromatogram}
\title{Deprecated functions in package \sQuote{xcms}}

\description{
  These functions are provided for compatibility with older versions of
  \sQuote{xcms} only, and will be defunct at the next release.
}

\details{

  The following functions/methods are deprecated.

  \itemize{
    \item{\code{xcmsPapply}: this function is no longer available and the use
      of \code{\link{bplapply}} is suggested.}

    \item{\code{profBin}, \code{profBinM}, \code{profBinLin},
      \code{profBinLinM}, \code{profBinLinBase}, \code{profBinLinBaseM}
      have been deprecated and \code{\link{binYonX}}} in combination
    with \code{\link{imputeLinInterpol}} should be used instead.

    \item{\code{extractChromatograms}}: replaced by \code{\link{chromatogram}}.

    \item{\code{plotChromatogram}}: replaced by \code{plot} method for
    \code{\link{Chromatogram}} or
    \code{\link{Chromatograms}} objects.

    \item{\code{extractMsData}}: replaced by \code{as(x, "data.frame")}.

    \item{\code{plotMsData}}: replaced by \code{plot(x, type = "XIC")}.
  }
}


