\name{xcmsProtocol-class}
\docType{class}
\alias{method}
\alias{method,xcmsProtocol-method}
\alias{stage}
\alias{stage,xcmsProtocol-method}
\alias{xcmsProtocol-class}
\alias{show,xcmsProtocol-method}

\title{A step in a pipeline}
\description{Encapsulates an analysis strategy that is performed on a
  dataset. The primary component of an \code{\linkS4class{xcmsPipeline}}.}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Methods}{
  \describe{
    \item{\link{perform}}{
      \code{signature(object = "xcmsProtocol", data = "xcmsData")}:
      Performs the analysis strategy on the given \code{data}. The
      protocol is added to the pipeline (analysis history) of the
      returned dataset.
    }
    \item{stage}{
      \code{signature(object = "xcmsProtocol")}:
      Returns an instance of the stage associated with the protocol.
    }
    \item{method}{
      \code{signature(object = "xcmsProtocol")}:
      Returns the identifier of the analytical method implemented by the
      protocol. This distinguishes the protocol from the other protocols
      for the same stage.
    }
    \item{\link{pipeline}}{
      \code{signature(object = "xcmsProtocol")}:
      If the protocol subclass has a \code{pipeline} field, return
      it. Otherwise, return \code{NULL}.
    }
    \item{show}{
      \code{signature(object = "xcmsProtocol")}:
      Display the protocol description and parameters.
    }
  }
}
\details{
  This class is virtual and does not itself implement an analysis
  stategy. By convention, each class that directly extends this class is
  also virtual and corresponds to a stage (see
  \code{\linkS4class{xcmsStage}}). The concrete protocols, those that
  actually implement an analysis method, extend, directly or indirectly,
  one of those stage base classes. This associates the protocol with its
  stage (role). Note that only single inheritance is assumed.

  An instance of each concrete (non-virtual) subclass provides:
  \itemize{
    \item Analysis strategy that may be \code{\link{perform}}ed on a dataset.
    \item Set of parameters (each corresponding to a slot in the class)
    controlling the analysis strategy.
    \item Character vector identifying the analysis strategy, see
    \code{\link{method}}.
    \item Human-readable name of the strategy, see
    \code{\link{dispName}}.
    \item Stage instance representing the role played by the protocol,
    see \code{\link{stage}}.
  }
}
\author{ Michael Lawrence }
\seealso{
  \code{\link{protocol}} and \code{\link{xcmsProtocol}} for obtaining a
  protocol instance, \code{\link{setProtocol}} for registering a new
  protocol type.
}
\keyword{classes}
\keyword{internal}
