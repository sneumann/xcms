\name{doubleMatrix}
\alias{doubleMatrix}
\alias{logicalMatrix}
\title{Allocate double or logical matricies}
\description{
  Allocate double or logical matricies in one step without copying
  memory around.
}
\usage{
doubleMatrix(nrow = 0, ncol = 0)
logicalMatrix(nrow = 0, ncol = 0)
}
\arguments{
  \item{nrow}{number of matrix rows}
  \item{ncol}{number of matrix columns}
}
\value{
  Matrix of double or logical values.
}
\author{Colin A. Smith, \email{csmith@scripps.edu}}
\keyword{programming}
\keyword{internal}
