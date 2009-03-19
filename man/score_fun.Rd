\name{score_fun}
\alias{score_fun}
\title{Scoring for MS/MS spectra}
\description{
  A similarity scoring function for MS/MS spectra against a reference.
}
\usage{
score_fun(ref, exp, method="distMatrix", ...)
}
\arguments{
  \item{ref}{An array of numbers for the refference}
  \item{exp}{An array of numbers for the test}
  \item{method}{A string of either 'distMatrix' or 'cor'}
  \item{...}{Any other arguments to be passed to other functions}
}
\details{
  This is a starting method to score MS/MS data. Which reports the parameters to the various one of the scoreing functions.
}
\value{
  \item{score}{Percentage score between the two arrays}
}
\seealso{
  \code{\link{score_fun.distMatrix}},
  \code{\link{score_fun.cor}}
}

\references{
  H. Paul Benton, D.M. Wong, S.A.Strauger, G. Siuzdak "XC\eqn{MS^2}:"
  Analytical Chemistry 2008 DOI:\url{http://pubs.acs.org/doi/abs/10.1021/ac800795f/}
}
\author{H. Paul Benton, \email{hpbenton@scripps.edu}}
