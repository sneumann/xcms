\name{setProtocol}
\alias{setProtocol}
\title{ Define a new protocol type }
\description{
  Utility for defining a new type of protocol - an analysis strategy -
  in the specified environment. Represented by a class derived from
  \code{\linkS4class{xcmsProtocol}}.
}
\usage{
setProtocol(method, dispname = method, representation = list(), fun, parent, prototype = list(), validity = NULL, where = topenv(parent.frame()))
}

\arguments{
  \item{method}{ Character vector identifying the analytical
    method. This identifier distinguishes the protocol from the other
    protocols associated with the same stage. See \code{\link{method}}. }
  \item{dispname}{Human-readable name for the protocol. See
    \code{\link{dispName}}.}
  \item{representation}{A list with an element for each parameter of the
    protocol. The value of an element is a character vector specifying the
    class for the parameter that matches the name of the element. The
    parameters are represented as slots in the protocol class. See
    \code{\link{representation}}.}
  \item{fun}{A function that implements the analysis strategy. It should
    take the dataset as its first argument. The other arguments should be
    named according to the names of the parameters in the representation,
    though it is not necessary for the set of function arguments and set of
    protocol parameters to be identical. The defaults (i.e. the class prototype)
    of the parameters are taken from the function formals. If \code{fun}
    is omitted, the protocol is virtual and is meant to serve as a base for
    other protocols.}
  \item{parent}{A string representing the single parent protocol
    class. This is analogous to the \code{contains} parameter in
    \code{\link{setClass}} except multiple inheritance is not allowed
    and it is not necessary to supply the full parent class name.  In
    particular, it is possible to specify the role identifier of the stage.}
  \item{prototype}{An optional list, see \code{\link{prototype}}, specifying
    default values for the parameters (slots). The formals of \code{fun}
    override this argument, which is mostly meant for defining virtual base
    protocols.}
  \item{validity}{An optional function for checking validity of
    instances, see \code{\link{setClass}}.}
  \item{where}{Environment in which the protocol type is
    registered. Defaults to the top-level environment. See
    \code{\link{setClass}}.}
}
\details{
  Defining a protocol type has the following side effects:
  \itemize{
    \item Defines a representative class derived from the class
    indicated by \code{parent}. The slots in the class may be of the
    type given in \code{representation} or a language type (call, name,
    etc). This is to permit lazy evaluation, as is commonly used in
    function formals. The prototype is based on the \code{prototype}
    argument, which is overriden by the formals of \code{fun}.
    \item Sets a method for \code{\link{dispName}} on the class to
    provide a human-readable name.
    \item Sets a method for \code{\link{performDelegate}}, which invokes
    the \code{fun} argument and is normally called by
    \code{\link{perform}}.
    \item Sets a generic named according to \code{role.method}, where
    \code{role} is the \code{\link{role}} of the associated stage and
    \code{method} is the same as the argument. This is the most clearly
    specified interface to the protocol and should be the focal point in
    its documentation.
  }
}
\value{
  The name of the protocol class.
}
\author{Michael Lawrence}
\seealso{\code{\link{xcmsProtocol}} and \code{\link{protocol}} for
  obtaining instances of protocols. \code{\link{setClass}} for
  registering stages.}
\examples{
\dontrun{
  setProtocol("median", "Median",
              representation(mzrad = "numeric", scanrad = "numeric"),
              function(object, mzrad = 0, scanrad = 0) {
                object@.Data <- medianFilter(object, mzrad, scanrad)
                object
              }, "filterProfile")
}
}
\keyword{internal}
