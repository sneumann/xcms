\name{pipeline-methods}
\docType{methods}
\alias{pipeline}
\alias{pipeline-methods}
\alias{pipeline,xcmsPipeline-method}
\alias{pipeline,xcmsProtocol-method}
\alias{pipeline,xcmsData-method}
\alias{pipeline,xcmsProfile-method}
\title{ Get a pipeline }
\description{
  Returns an \code{\linkS4class{xcmsPipeline}} from a container.
}
\section{Methods}{
\describe{

  \item{object = "xcmsPipeline"}{
    \code{pipeline(object, intype = "ANY", outtype = "ANY")} - extracts
    as a sub-pipeline the longest sequence of protocols in
    \code{object} that begins with a protocol accepting \code{intype} as
    input and ends with a protocol with \code{outtype} as its output.
  }

  \item{object = "xcmsProtocol"}{
    If \code{object} contains a \code{pipeline} slot return it;
    otherwise, return \code{NULL}. This supports "compound"
    protocols that contain (and delegate) to pipelines.
  }

  \item{object = "xcmsData"}{
    \code{pipeline(object, local = FALSE)} - If
    \code{local} is \code{FALSE}, return the \code{pipeline} slot (the
    entire pipeline that produced the object). If \code{local} is
    \code{TRUE}, return only the protocols that were performed during or
    after the conversion of the dataset to the current data type. That is,
    return the last contiguous sequence of protocols where the output type
    of each protocol matches the type of \code{object}.
  } 

  \item{object = "xcmsProfile"}{
    Returns its \code{pipeline} slot as an instance of
    \code{\linkS4class{xcmsPipelineProfile}}.
  }
  
}}
\author{Michael Lawrence}
\keyword{methods}
