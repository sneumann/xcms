\name{verify.mzQuantM}
\alias{verify.mzQuantML}
\title{Verify an mzQuantML file}
\description{
  Export in XML data formats: verify the written data
}

\usage{
  verify.mzQuantML(filename, xsdfilename)
}

\arguments{
  \item{filename}{filename (may include full path)
    for the output file. Pipes or URLs are not allowed.
  }
  \item{xsdfilename}{Filename of the XSD to verify against (may include full path)
  }
}


\details{
  The verify.mzQuantML() function will verify an PSI standard format mzQuantML
  document against the XSD schemda, see \url{http://www.psidev.info/mzquantml}
  
}
\value{
  None.
}
\seealso{
  \code{\link{write.mzQuantML}}
}
\keyword{file}
