\name{xcmsFileSource-class}
\Rdversion{1.1}
\docType{class}
\alias{xcmsFileSource-class}
\alias{xcmsSource,character-method}
\title{Base class for loading raw data from a file}
\description{
  Data sources which read data from a file should inherit from this
  class. The \code{xcms} package provides classes to read from
  \code{netCDF}, \code{mzData}, \code{mzXML}, and \code{mzML} files
  using \code{xcmsFileSource}.

  This class should be considered virtual and will not work if passed to
  \code{\link{loadRaw-methods}}. The reason it is not explicitly
  virtual is that there does not appear to be a way for a class to be
  both virtual and have a data part (which lets functions treat objects
  as if they were character strings).

  This class validates that a file exists at the path given.
}
\section{Objects from the Class}{
  \code{xcmsFileSource} objects should not be instantiated directly.
  Instead, create subclasses and instantiate those.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"character"}. File path
      of a file from which to read raw data as the object's data part}
  }
}
\section{Extends}{
  Class \code{"\linkS4class{character}"}, from data part.
  Class \code{"\linkS4class{xcmsSource}"}, directly.
}
\section{Methods}{
  \describe{
    \item{\code{xcmsSource}}{
      \code{signature(object = "character")}: Create an
      \code{xcmsFileSource} object referencing the given file name.
    }
  }
}

\author{Daniel Hackney \email{dan@haxney.org}}
\seealso{\code{\linkS4class{xcmsSource}}}
\keyword{classes}
\keyword{file}
