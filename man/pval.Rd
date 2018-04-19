\name{pval}
\alias{pval}
\title{Generate p-values for a vector of t-statistics}
\description{
  Generate p-values for a vector of Welch's two-sample t-statistics
  based on the t distribution.
}
\usage{
pval(X, classlabel, teststat)
}
\arguments{
  \item{X}{original data matrix}
  \item{classlabel}{integer vector with classlabel}
  \item{teststat}{numeric vector with Welch's two-sample t-statistics}
}
\value{
  A numeric vector of p-values.
}
\author{Colin A. Smith, \email{csmith@scripps.edu}}
\seealso{
  \code{\link{mt.teststat}}
}
\keyword{univar}
\keyword{internal}
