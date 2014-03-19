\name{rectUnique}
\alias{rectUnique}
\title{Determine a subset of rectangles with unique, non-overlapping areas}
\description{
  Given a matrix of rectangular areas, this function determines a subset
  of those rectangles that do not overlap. Rectangles are preserved on
  a first come, first served basis, with user control over the order in
  which the rectangles are processed.
}
\usage{
rectUnique(m, order = seq(length = nrow(m)), xdiff = 0, ydiff = 0)
}
\arguments{
  \item{m}{four column matrix defining rectangular areas}
  \item{order}{order in which matrix columns should be scanned}
  \item{xdiff}{maximum space between overlapping rectangles in x dimension}
  \item{ydiff}{maximum space between overlapping rectangles in y dimension}
}
\details{
  The \code{m} matrix must contain four colums defining the position of
  rectangle sides in the folloing order: left, right, bottom, top. This
  function is currently implemented in \code{C} using a an algorithm
  with quadratic running time.
}
\value{
  A logical vector indicating which rows should be kept.
}
\author{Colin A. Smith, \email{csmith@scripps.edu}}
\examples{
m <- rbind(c(0,4,0,3), c(1,3,2,6), c(3,6,4,6))
plot(0, 0, type = "n", xlim=range(m[,1:2]), ylim=range(m[,3:4]))
rect(m[,1], m[,3], m[,2], m[,4])
xcms:::rectUnique(m)
# Changing order of processing
xcms:::rectUnique(m, c(2,1,3))
# Requiring border spacing
xcms:::rectUnique(m, ydiff = 1)
# Allowing adjacent boxes
xcms:::rectUnique(m, c(2,1,3), xdiff = -0.00001)
# Allowing interpenetration
xcms:::rectUnique(m, xdiff = -1.00001, ydiff = -1.00001)
}
\keyword{array}
\keyword{internal}
