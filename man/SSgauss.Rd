\name{SSgauss}
\alias{SSgauss}
\title{Gaussian Model}
\description{
  This \code{selfStart} model evalueates the Gaussian model and its
  gradient. It has an \code{initial} attribute that will evalueate
  the inital estimates of the parameters \code{mu}, \code{sigma},
  and \code{h}.
}
\usage{
SSgauss(x, mu, sigma, h)
}
\arguments{
  \item{x}{a numeric vector of values at which to evaluate the model}
  \item{mu}{mean  of the distribution function}
  \item{sigma}{standard deviation of the distribution fuction}
  \item{h}{height of the distribution function}
}
\details{
  Initial values for \code{mu} and \code{h} are chosen from the
  maximal value of \code{x}. The initial value for \code{sigma} is
  determined from the area under \code{x} divided by \code{h*sqrt(2*pi)}.
}
\value{
  A numeric vector of the same length as \code{x}. It is the value
  of the expression \code{h*exp(-(x-mu)^2/(2*sigma^2)}, which is a
  modified gaussian function where the maximum height is treated
  as a separate parameter not dependent on \code{sigma}. If arguments
  \code{mu}, \code{sigma}, and \code{h} are names of objects, the
  gradient matrix with respect to these names is attached as an
  attribute named \code{gradient}.
}
\author{Colin A. Smith, \email{csmith@scripps.edu}}
\seealso{
  \code{\link{nls}},
  \code{\link{selfStart}}
}
\keyword{nonlinear}
