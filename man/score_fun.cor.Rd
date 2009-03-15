\name{score_fun.cor}
\alias{score_fun.cor}
\title{Scoring for MS/MS spectra Via correlation}
\description{
  A similarity scoring function for MS/MS spectra against a reference via correlation.
}
\usage{
  score_fun.cor(ref, exp)
}
\arguments{
  \item{ref}{An array of numbers for the refference}
  \item{exp}{An array of numbers for the test}
}
\details{
  A score system using correlation analysis to correlate two arrays. If the arrays are a different length then the system will buffer the shorter one with 0s.
}
\value{
  \item{score}{Correlation between the two arrays}
}
\seealso{
  \code{\link{score_fun.distMatrix}}
}

\references{
  H. Paul Benton, D.M. Wong, S.A.Strauger, G. Siuzdak "XC\eqn{MS^2}:"
  Analytical Chemistry 2008 DOI:\url{http://pubs.acs.org/doi/abs/10.1021/ac800795f/}{10.1021/ac800795f}
}
\author{H. Paul Benton, \email{hpbenton@scripps.edu}}

\examples{
\dontrun{
score<-xcms:::score_fun.cor(rnorm(10), rnorm(10))
score

a<-abs(rnorm(5))
a[2]<-xcms:::ppmDev(a[2], 30)
score<-xmcs:::score_fun.cor()
score
}
}