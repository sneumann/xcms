\name{etg}
\alias{etg}
\title{Empirically Transformed Gaussian function}
\description{
  A general function for asymmetric chromatographic peaks.
}
\usage{
etg(x, H, t1, tt, k1, kt, lambda1, lambdat, alpha, beta)
}
\arguments{
  \item{x}{times to evaluate function at}
  \item{H}{peak height}
  \item{t1}{time of leading edge inflection point}
  \item{tt}{time of trailing edge inflection point}
  \item{k1}{leading edge parameter}
  \item{kt}{trailing edge parameter}
  \item{lambda1}{leading edge parameter}
  \item{lambdat}{trailing edge parameter}
  \item{alpha}{leading edge parameter}
  \item{beta}{trailing edge parameter}
}
\value{
  The function evaluated at times \code{x}.
}
\references{
  Jianwei Li. Development and Evaluation of Flexible Empirical Peak
  Functions for Processing Chromatographic Peaks. Anal. Chem., 69
  (21), 4452-4462, 1997.
  \url{http://dx.doi.org/10.1021/ac970481d}
}
\author{Colin A. Smith, \email{csmith@scripps.edu}}
\keyword{models}
