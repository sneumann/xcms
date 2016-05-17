\name{writeMzTab}
\alias{writeMzTab}
\title{Save a grouped xcmsSet object in mzTab-1.1 format file}
\description{
  Write the grouped xcmsSet to an mzTab file.
}
\section{Usage}{
\describe{
\item{object = "xcmsSet"}{
  \code{writeMzTab(object, filename)}
}
}}
\arguments{
  \item{object}{the \code{xcmsSet} object}
  \item{filename}{filename (may include full path)
    for the mzTab file. Pipes or URLs are not allowed.
  }
}
\details{
   The mzTab file format for MS-based metabolomics (and proteomics)
   is a lightweight
   supplement to the existing standard XML-based file formats (mzML,
   mzIdentML, mzQuantML), providing a comprehensive summary, similar in
   concept to the supplemental material of a scientific
   publication. mzTab files from xcms contain small
   molecule sections together with experimental metadata and
   basic quantitative information. The format is intended to store
   a simple summary of the final results.
}
\value{
  None.
}
\examples{
    library(faahKO)
    xs <- group(faahko)

    mzt <- data.frame(character(0))
    mzt <- xcms:::mzTabHeader(mzt,
                       version="1.1.0", mode="Complete", type="Quantification",
                       description="faahKO",
                       xset=xs)
    mzt <- xcms:::mzTabAddSME(mzt, xs)
    
    xcms:::writeMzTab(mzt, "faahKO.mzTab")
}

\seealso{
  \code{\link{xcmsSet-class}},
  \code{\link{xcmsSet}},
}
\keyword{file}
