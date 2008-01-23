\name{xcmsStage-class}
\docType{class}
\alias{xcmsStage-class}
\alias{defaultMethod}
\alias{defaultMethod,xcmsStage-method}
\alias{role,xcmsStage-method}
\alias{role}

\title{A stage in a pipeline}
\description{Represents a role, such as finding peaks or normalizing
  quanitities, played by one or more \code{\linkS4class{xcmsProtocol}}
  objects in a pipeline. Stores metadata, such as the input and output
  type, for the role, and also serves as a factory for protocols that
  play the role.}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Methods}{
  \describe{
    \item{defaultMethod}{
      \code{signature(object = "xcmsStage")}: Returns
      the name of the analytical method (corresponding to a protocol
      type) that is registered as the default for the stage in the global
      \code{\link{options}}.
    }
    \item{defaultMethod<-}{
      \code{signature(object = "xcmsStage")}:
      Replaces the name of the analytical method registered as the default
      for the stage in the global options.
    }
    \item{protocol}{
      \code{signature(object = "xcmsStage")}: Creates a protocol that
      implements this stage given a method name. See \code{\link{protocol}}.
    }
    \item{role}{
      \code{signature(object = "xcmsStage")}: Gets the string
      identifying the role of this stage.
    }
  }
}
\details{
  Each stage is represented by a direct subclass of this virtual base class.
  An instance of a subclass provides the following:
  \itemize{
    \item String identifying its role, see above.
    \item Human-readable name of its role, see \code{\link{dispName}}.
    \item Default method name (protocol), see above.
    \item Names of the input and output types, see \code{\link{inType}}
    and \code{\link{outType}}.
    \item Protocol instances that implement the stage using a
    specified analytical method, see \code{\link{protocol}}.
  }
  A stage should be thought of as a singleton, in that all of its
  state is static and thus all stage instances of the same class are essentially
  identical. This tries to ensure that the parameters are the same for a
  given role across the R session.
}
\author{ Michael Lawrence }
\note{
  One could think of stages as extension points in the Eclipse
  framework or more generally as interfaces in OOP languages.
}

\seealso{
  \code{\link{stage}} and \code{\link{xcmsStage}} for obtaining stage
  instances, \code{\link{setClass}} for registering a new stage.
}
\keyword{classes}
\keyword{internal}