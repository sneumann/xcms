\name{score_fun}
\alias{score_fun}
\title{Scoring for MS/MS spectra}
\description{
  A similarity scoring function for MS/MS spectra against a reference.
}
\usage{
score_fun(ref, exp, ppmfrag)
}
\arguments{
  \item{ref}{An array of numbers for the refference}
  \item{exp}{An array of numbers for the test}
  \item{ppmfrag}{Error to match to in ppm}
}
\details{
  A simple scoring function to score two arrays of numbers and give a percentage match between the two. Uses a a distance and similarity matrix score system. When the two scores are calculated the percentage score is calculated from the theoritical maximum score and the theoritical minimum score.
}
\value{
  \item{score}{Percentage score between the two arrays}
}
\references{
  H. Paul Benton, D.M. Wong, S.A.Strauger, G. Siuzdak "XC\eqn{MS^2}:"
  Analytical Chemistry 2008
}
\author{H. Paul Benton, \email{hpbenton@scripps.edu}}
