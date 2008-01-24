\name{xcmsProfile-class}
\docType{class}
\alias{xcmsProfile-class}
\alias{mzBreaks}
\alias{mzBreaks,xcmsProfile-method}
\alias{selectRange}
\alias{selectRange,xcmsProfile-method}
\alias{show,xcmsProfile-method}

\title{A profile matrix}
\description{Associates a profile matrix with an
  \code{\linkS4class{xcmsPipeline}} and other metadata.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("xcmsProfile", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{The matrix holding the profile data}
    \item{\code{mzrange}:}{The range in m/z for the matrix}
    \item{\code{step}:}{The m/z step between each matrix row}
    \item{\code{scantime}:}{The times for each column in the matrix}
    \item{\code{pipeline}:}{The pipeline holding the analysis history}
  }
}
\section{Extends}{
Class \code{"\linkS4class{xcmsData}"}, directly.
Class \code{"\linkS4class{matrix}"}, from data part.
Class \code{"\linkS4class{array}"}, by class "matrix", distance 2.
Class \code{"\linkS4class{structure}"}, by class "matrix", distance 3.
Class \code{"\linkS4class{vector}"}, by class "matrix", distance 4, with explicit coerce.
}
\section{Methods}{
  \describe{
    \item{mzBreaks}{\code{signature(object = "xcmsProfile")}: The m/z
      values for each row. See also \code{\link{profMz}}.}
    \item{\link{perform}}{
      \code{signature(object = "xcmsProtoFilterProfileBase",
        data = "xcmsProfile")}:
      Filter the profile matrix.
    }
    \item{\link{pipeline}}{
      \code{signature(object = "xcmsProfile")}:
      Get the pipeline slot as an
      \code{\linkS4class{xcmsPipelineProfile}}.
    }
    \item{\link{plotChrom}}{
      \code{signature(object = "xcmsProfile")}:
      Plot a chromatogram.
    }
    \item{\link{plotSpec}}{
      \code{signature(object = "xcmsProfile")}:
      Plot mass spectra.
    }
    \item{\link{plotSurf}}{
      \code{signature(object = "xcmsProfile")}:
      Plot 3D surface of the matrix.
    }
    \item{selectRange}{
      \code{signature(object = "xcmsProfile")}:
      Utility for obtaining a subset specification for a profile matrix.
      Same behavior and interface as \code{\link{profRange}}.
    }
    \item{show}{
      \code{signature(object = "xcmsProfile")}:
      Display a textual description of an instance.
    }
  }
}
\author{ Michael Lawrence }
\seealso{
  \code{\link{profileMatrix}} for generating a profile matrix.
}
\keyword{classes}
