\name{write.cdf-methods}
\docType{methods}
\alias{write.cdf}
\alias{write.cdf,xcmsRaw-method}
\title{Save an xcmsRaw object to file}
\description{
  Write the raw data to a (simple) CDF file.
}
\section{Methods}{
\describe{
\item{object = "xcmsRaw"}{
  \code{write.cdf(object, filename)}
}
}}
\arguments{
  \item{object}{the \code{xcmsRaw} object}
  \item{filename}{filename (may include full path)
    for the CDF file. Pipes or URLs are not allowed.
  }
}
\details{
  Currently the only application known to read the resulting file is
  XCMS. Others, especially those which build on the AndiMS library,
  will refuse to load the output.
}
\value{
  None.
}
\seealso{
  \code{\link{xcmsRaw-class}},
  \code{\link{xcmsRaw}},
}
\keyword{methods}
\keyword{file}
