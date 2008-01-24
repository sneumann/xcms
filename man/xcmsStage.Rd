\name{xcmsStage}
\alias{xcmsStage}
\title{ Create a stage instance }
\description{
  Returns an instance of \code{\linkS4class{xcmsStage}} given a role
  identifier.
}
\usage{
xcmsStage(role)
}
\arguments{
  \item{role}{ The character vector identifying the role of the stage. }
}
\value{
  An instance of \code{\linkS4class{xcmsStage}}.
}
\author{ Michael Lawrence }
\seealso{ \code{\link{stage}} for getting a stage from a container. }
\examples{
\dontrun{
  xcmsStage("findPeaks")
}
}
\keyword{internal}
