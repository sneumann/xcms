\name{score_fun.distMatrix}
\alias{score_fun.distMatrix}
\title{Scoring for MS/MS spectra}
\description{
  A similarity scoring function for MS/MS spectra against a reference via a distance matrix.
}
\usage{
  score_fun.distMatrix(ref, exp, ppmfrag)
}
\arguments{
  \item{ref}{An array of numbers for the refference}
  \item{exp}{An array of numbers for the test}
  \item{ppmfrag}{A numerical string for the amount of error in }
}
\details{
   A simple scoring function to score two arrays of numbers and give a percentage match between the two. Uses a a distance and similarity matrix score system. When the two scores are calculated the percentage score is calculated from the theoritical maximum score and the theoritical minimum score.
}
\value{
  \item{score}{Percentage score between the two arrays}
}
\seealso{
  \code{\link{score_fun.cor}}
}

\references{
  H. Paul Benton, D.M. Wong, S.A.Strauger, G. Siuzdak "XC\eqn{MS^2}:"
  Analytical Chemistry 2008 DOI:\url{http://pubs.acs.org/doi/abs/10.1021/ac800795f/}{10.1021/ac800795f}
}
\author{H. Paul Benton, \email{hpbenton@scripps.edu}}

\examples{
\dontrun{
score<-xcms:::score_fun.distMatrix(rnorm(10), rnorm(10), 20)
score

a<-abs(rnorm(5))
a[2]<-xcms:::ppmDev(a[2], 30)
score<-xmcs:::score_fun.distMatrix()
score
}
}