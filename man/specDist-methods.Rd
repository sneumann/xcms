\name{specDist-methods}
\docType{methods}
\alias{specDist-methods}
\alias{specDist}
\alias{specDist,xcmsSet-method}
\title{Distance methods for xcmsSet, xcmsRaw and xsAnnotate}
\description{
  There are several methods for calculating a distance between two sets of peaks in xcms. \code{specDist}
  is the generic method.
}
\section{Methods}{
	\describe{
	\item{object = "xcmsSet"}{
		\code{
		specDist(object, peakIDs1, peakIDs2,...)
		}
		}
	}
	\describe{
	\item{object = "xsAnnotate"}{
	\code{
	specDist(object, PSpec1, PSpec2,...)
		}
		}
	}

}


\arguments{
  \item{object}{a xcmsSet or xcmsRaw.}
  \item{method}{Method to use for distance calculation. See details.}
  \item{...}{mzabs, mzppm and parameters for the distance function.}
}
\details{
  Different algorithms can be used by specifying them with the
  \code{method} argument. For example to use the "meanMZmatch"
  approach with xcmsSet one would use:
  \code{specDist(object, peakIDs1, peakIDs2, method="meanMZmatch")}. This is also
  the default.

  Further arguments given by \code{...} are
  passed through to the function implementing
  the \code{method}.

  A character vector of \emph{nicknames} for the
  algorithms available is returned by
  \code{getOption("BioC")$xcms$specDist.methods}.
  If the nickname of a method is called "meanMZmatch",
  the help page for that specific method can
  be accessed with \code{?specDist.meanMZmatch}.
}
\value{
\item{mzabs}{
    maximum absolute deviation for two matching peaks
  }
  \item{mzppm}{
    relative deviations in ppm for two matching peaks
  }
  \item{symmetric}{
    use symmetric pairwise m/z-matches only, or each match
    }
}

\author{Joachim Kutzera, \email{jkutzer@ipb-halle.de}}
\keyword{methods}
