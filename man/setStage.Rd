\name{setStage}
\alias{setStage}
\title{ Define a new stage }
\description{
  Defines a stage - a role for protocols to perform in a pipeline - in
  the specified environment. The new stage is represented by a subclass
  of \code{\linkS4class{xcmsStage}}.
}
\usage{
setStage(name, dispname = name, intype = "xcmsSet", outtype = intype, where = topenv(parent.frame()))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{name}{ The role of the stage. }
  \item{dispname}{ The human-readable name, see \code{\link{dispName}}. }
  \item{intype}{ The input data type for all protocols implementing the
    stage, see \code{\link{inType}}. }
  \item{outtype}{ The output data type for all protocols implementing
    the stage, see \code{\link{outType}}. }
  \item{where}{ The environment in which the stage is defined. Not
    normally specified, see \code{\link{setClass}}.}
}
\details{
  This function has the following side-effects:
  \itemize{
    \item Registers a stage in the pipeline framework by defining
      a subclass of \code{\linkS4class{xcmsStage}} named according to
      \code{xcmsStageName} where \code{Name} is specified by the \code{name}
      parameter. The subclass represents an analytical role that may be
      played by an \code{\linkS4class{xcmsProtocol}}.
    \item Registers a subclass of \code{\linkS4class{xcmsProtocol}} that
      serves as the base class for all protocols performing the role
      represented by the stage. It is named according to
      \code{xcmsProtoName} where \code{Name} is from the \code{name}
      parameter.
    \item Defines methods for \code{\link{dispName}}, on the stage
      class, and \code{\link{inType}} and \code{\link{outType}} on the
      stage class and protocol base class.
    \item Defines a generic method for performing the stage via the command
      line interface with the usage:
      \code{name(object, method = \link{defaultMethod}(name))},
      where \code{name} is name of the stage, \code{object} is the data
      (of type \code{intype}) on which the protocol specified by
      \code{method} is performed.
    \item Defines a generic, named according to \code{nameProto}, where
      \code{name} is the name of the stage. The purpose is to retrieve
      a protocol that performs the stage. Methods are defined for
      \code{\linkS4class{xcmsPipeline}} and \code{outtype}. Also defines
      a replace method of the same name on \code{xcmsPipeline}.
  }
}
\value{
  The \code{name} parameter.
}

\author{ Michael Lawrence }
\seealso{ \code{\link{setProtocol}} for defining protocol types,
  \code{\link{setClass}} to which this function is similar in spirit. }
\keyword{methods}
\keyword{internal}