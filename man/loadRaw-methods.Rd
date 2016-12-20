\name{loadRaw-methods}
\alias{loadRaw}
\alias{loadRaw-methods}
\alias{loadRaw,xcmsSource-method}
\alias{loadRaw,xcmsFileSource-method}
\title{Read binary data from a source}
\description{
  This function extracts the raw data which will be used an
  \code{\linkS4class{xcmsRaw}} object. Further processing of data is
  done in the \code{\link{xcmsRaw}} constructor.
}
\section{Methods}{
  \describe{
    \item{\code{signature(object = "xcmsSource")}}{
      Uses \code{\link{loadRaw,xcmsSource-method}} to extract raw data.
      Subclasses of \code{\linkS4class{xcmsSource}} can provide different
      ways of fetching data.
    }
  }
}

\arguments{
  \item{object}{Specification of a data source (such as a file name or
    database query)}
}
\details{
  The implementing methods decide how to gather the data.
}
\value{

  A list containing elements describing the data source. The \code{rt},
  \code{scanindex}, \code{tic}, and \code{acquisitionNum} components
  each have one entry per scan. They are \emph{parallel} in the sense that
  \code{rt[1]}, \code{scanindex[1]}, and \code{acquisitionNum[1]} all
  refer to the same scan. The list containst the following components:

  \item{\code{rt}}{
    Numeric vector with acquisition time (in seconds) for each scan
  }
  \item{\code{tic}}{Numeric vector with Total Ion Count for each scan}
  \item{\code{scanindex}}{
    Integer vector with starting positions of each scan in the \code{mz}
    and \code{intensity} components. It is an exclusive offset, so
    \code{scanindex[i]} is the offset in \code{mz} and \code{intensity}
    \emph{before} the beginning of scan \code{i}. This means that the
    \code{mz} (respectively \code{intensity}) values for scan \code{i}
    would be from \code{scanindex[i] + 1} to \code{scanindex[i + 1]}
  }
  \item{\code{mz}}{Concatenated vector of m/z values for all scans}
  \item{\code{intensity}}{Concatenated vector of intensity values for all scans}
}

\author{Daniel Hackney, \email{dan@haxney.org}}
\seealso{
  \code{\link{xcmsRaw-class}},
  \code{\linkS4class{xcmsSource}}
}
\keyword{methods}
