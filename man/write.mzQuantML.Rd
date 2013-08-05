\name{write.mzQuantML-methods}
\docType{methods}
\alias{write.mzQuantML}
\alias{write.mzQuantML,xcmsSet-method}
\title{Save an xcmsSet object to an PSI mzQuantML file}
\description{
  Export in XML data formats: 
  Write the processed data in an xcmsSet to mzQuantML.
}

\section{Methods}{
\describe{
\item{object = "xcmsSet"}{
  \code{write.mzQuantML(object, filename)}
}
}
}
\arguments{
  \item{object}{the \code{xcmsRaw} or \code{xcmsSet} object}
  \item{filename}{filename (may include full path)
    for the output file. Pipes or URLs are not allowed.
  }
}


\details{
  The write.mzQuantML() function will write a (grouped) xcmsSet
  into the PSI standard format mzQuantML, see \url{http://www.psidev.info/mzquantml}
  
}
\value{
  None.
}
\seealso{
  \code{\link{xcmsSet-class}},
  \code{\link{xcmsSet}},
  \code{\link{verify.mzQuantML}},
}
\keyword{methods}
\keyword{file}
