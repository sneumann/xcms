\name{xcmsPapply}
\alias{xcmsPapply}
\title{\emph{Deprecated}: xcmsPapply }
\description{
  This function is deprecated, use \code{\link{bplapply}}
  instead.

  An apply-like function which uses Rmpi to distribute
  the processing evenly across a cluster.  Will use a non-MPI
  version if distributed processing is not available.
}
\usage{
    xcmsPapply(arg_sets, papply_action, papply_commondata = list(),
           show_errors = TRUE, do_trace = FALSE, also_trace = c())
}
\arguments{
  \item{arg_sets}{ a list, where each item will be given as an argument
                   to papply\_action }
  \item{papply_action}{ A function which takes one argument.  It will be
                        called on each element of arg\_sets }
  \item{papply_commondata}{ A list containing the names and values of
                        variables to be accessible to the papply\_action.
                        'attach' is used locally to import this list. }
  \item{show_errors}{ If set to TRUE, overrides Rmpi's default, and
                      messages for errors which occur in R slaves are
                      produced. }
  \item{do_trace}{ If set to TRUE, causes the papply\_action function to
                   be traced. i.e. Each statement is output before it is
                   executed by the slaves. }
  \item{also_trace}{ If supplied an array of function names, as strings,
                     tracing will also occur for the specified functions. }

}
\details{
  Similar to apply and lapply, applies a function to all items
   of a list, and returns a list with the corresponding results.

   Uses Rmpi to implement a pull idiom in order to distribute
   the processing evenly across a cluster.
   If Rmpi is not available, or there are no slaves,
   implements this as a non-parallel algorithm.

  \code{xcmsPapply} is a modified version of the papply function from package papply 0.2 (Duane Currie).
  Parts of the slave function were wrapped in try() to make it failsafe and progress output was added.

  Make sure \code{Rmpi} was installed properly by executing the example below.
  \code{Rmpi} was tested with
  \itemize{
  \item OpenMPI : Unix, \url{http://www.open-mpi.org/}, don't forget to export \code{MPI_ROOT} before installing Rmpi e.g. \code{export MPI_ROOT=/usr/lib/openmpi}
  \item DeinoMPI : Windows, \url{http://mpi.deino.net/}, also see \url{http://www.stats.uwo.ca/faculty/yu/Rmpi/}
  }
}

\value{
  A list of return values from papply\_action.  Each value corresponds to
  the element of arg\_sets used as a parameter to papply\_action
}
\references{ \url{http://ace.acadiau.ca/math/ACMMaC/software/papply/} }
\author{ Duane Currie <duane.currie@acadiau.ca>, modified by Ralf Tautenhahn <rtautenh@ipb-halle.de>. }
\note{ Does not support distributing recursive calls in parallel.  If
  papply is used inside papply\_action, it will call a non-parallel
  version }

\examples{
\dontrun{
library(Rmpi)
library(xcms)

number_lists <- list(1:10,4:40,2:27)

mpi.spawn.Rslaves(nslaves=2)

results <- xcmsPapply(number_lists,sum)
results

mpi.close.Rslaves()

}
}
