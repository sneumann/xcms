\name{group-methods}
\docType{methods}
\alias{group-methods}
\alias{group}
\alias{group,xcmsSet-method}
\title{Group peaks from different samples together}
\description{
 Group peaks together across samples using overlapping m/z bins and
 calculation of smoothed peak distributions in chromatographic time.
}
\section{Methods}{
\describe{
\item{object = "xcmsSet"}{
  \code{group(object, mzwid = 0.25, bw = 30, max = 5, sleep = 0)}
}
}}
\arguments{
  \item{object}{the \code{xcmsSet} object}
  \item{mzwid}{
    width of overlapping m/z slices to use for creating peak density
    chromatograms and grouping peaks across samples
  }
  \item{bw}{
    bandwidth of gaussian smoothing kernel to apply to the peak
    density chromatogram
  }
  \item{max}{
    maximum number of groups to identify in a single m/z slice
  }
  \item{sleep}{
    seconds to pause between plotting successive steps of the peak
    grouping algorithm. peaks are plotted as points showing relative
    intensity. identified groups are flanked by dotted vertical
    lines.
  }
}
\value{
  An \code{xcmsSet} object with peak group assignments and statistics.
}
\seealso{
  \code{\link{xcmsSet-class}},
  \code{\link{density}}
}
\keyword{methods}
\keyword{file}
