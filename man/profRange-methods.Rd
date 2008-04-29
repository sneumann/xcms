\name{profRange-methods}
\docType{methods}
%\alias{profRange-methods}
\alias{profRange}
\alias{profRange,xcmsRaw-method}
\title{Specify a subset of profile mode data}
\description{
  Specify a subset of the profile mode matrix given a mass, time,
  or scan range. Allow flexible user entry for other functions.
}
\section{Methods}{
\describe{
\item{object = "xcmsRaw"}{
  \code{profRange(object, massrange = numeric(),
                  timerange = numeric(), scanrange = numeric(),
                  ...)}
}
}}
\arguments{
  \item{object}{the \code{xcmsRaw} object}
  \item{massrange}{single numeric mass or vector of masses}
  \item{timerange}{
    single numeric time (in seconds) or vector of times
  }
  \item{scanrange}{single integer scan index or vector of indecies}
  \item{...}{arguments to other functions}
}
\details{
  This function handles selection of mass/time subsets of the profile
  matrix for other functions. It allows the user to specify such
  subsets in a variety of flexible ways with minimal typing.

  Because \R does partial argument matching, \code{massrange},
  \code{scanrange}, and \code{timerange} can be specified in short
  form using \code{m=}, \code{s=}, and \code{t=}, respectively. If
  both a \code{scanrange} and \code{timerange} are specified, then
  the \code{timerange} specification takes precedence.

  When specifying ranges, you may either enter a single number or
  a numeric vector. If a single number is entered, then the closest
  single scan or mass value is selected. If a vector is entered,
  then the range is set to the \code{range()} of the values entered.
  That allows specification of ranges using shortened, slightly
  non-standard syntax. For example, one could specify 400 to 500
  seconds using any of the following: \code{t=c(400,500)},
  \code{t=c(500,400)}, or \code{t=400:500}. Use of the sequence
  operator (\code{:}) can save several keystrokes when specifying
  ranges. However, while the sequence operator works well for
  specifying integer ranges, fractional ranges do not always work
  as well.
}
\value{
  A list with the folloing items:
  
  \item{massrange}{numeric vector with start and end mass}
  \item{masslab}{textual label of mass range}
  \item{massidx}{integer vector of mass indecies}
  \item{scanrange}{integer vector with stat ane end scans}
  \item{scanlab}{textual label of scan range}
  \item{scanidx}{integer vector of scan range}
  \item{timerange}{numeric vector of start and end times}
  \item{timelab}{textual label of time range}
}
\seealso{
  \code{\link{xcmsRaw-class}}
}
\keyword{methods}
\keyword{manip}
