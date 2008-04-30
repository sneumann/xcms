\name{xcmsRaw-class}
\docType{class}
\alias{xcmsRaw-class}
\alias{profMz,xcmsRaw-method}
\alias{profMz}
\alias{revMz,xcmsRaw-method}
\alias{revMz}
\alias{show,xcmsRaw-method}
\alias{sortMz,xcmsRaw-method}
\alias{sortMz}
\alias{write.cdf,xcmsRaw-method}
\alias{write.cdf}

\title{Class xcmsRaw, a class for handling raw data}
\description{
  This class handles processing and visualization of the raw data
  from a single LC/MS or GS/MS run. It includes methods for producing
  a standard suite of plots including individual spectra, multi-scan
  average spectra, TIC, and EIC. It will also produce a feature
  list of significant peaks using matched filtration.
}
\section{Objects from the Class}{
  Objects can be created with the \code{\link{xcmsRaw}} constructor
  which reads data from a NetCDF file into a new object.
}
\section{Slots}{
  \describe{
    \item{\code{env}:}{
      environment with three variables: \code{mz} - concatenated
      m/z values for all scans, \code{intensity} - corresponding
      signal intensity for each m/z value, and \code{profile} -
      matrix represention of the intensity values with columns
      representing scans and rows representing equally spaced m/z
      values.

      For objects with MS$^n$ spectra also two matrices \code{msninfo}
      and \code{MSNpeaks} containing the level (e.g.\ 2 for Tandem MS)
      and other precursor information. The \code{MSNpeaks} matrix
      contains the raw spectra with their m/z and intensity.
    }
    \item{\code{tic}:}{
      numeric vector with total ion count (intensity) for each scan
    }
    \item{\code{scantime}:}{
      numeric vector with acquisition time (in seconds) for each
      scan
    }
    \item{\code{scanindex}:}{
      integer vector with starting positions of each scan in the
      \code{mz} and \code{intensity} variables (note that index
      values are based off a 0 initial position instead of 1)
    }
    \item{\code{profmethod}:}{
      characer value with name of method used for generating the
      profile matrix
    }
    \item{\code{mzrange}:}{
      numeric vector of length 2 with minimum and maximum m/z values
      represented in the profile matrix
    }
    \item{\code{gradient}:}{
      matrix with first row, \code{time}, containing the time point
      for interpolation and successive columns representing solvent
      fractions at each point
    }
    \item{\code{Rt}:}{
      a numeric vector of retention times of the MS$^n$ scan
    }
    \item{\code{msnLevel}:}{
      a numeric vector giving the fragmentation level eg MS$^2$ msLevel = 2
    }
    \item{\code{msnIntensity}:}{
      a numeric vector giving the intensity of the precursor ion.
    }
    \item{\code{msnPrecursorMZ}:}{
      a numeric vector of mass/charge of the precursor ion.
    }
    \item{\code{msnCollisionEnergy}:}{
      a numeric vector giving all collision energys for the MSn scans
    }
    \item{\code{filepath}:}{
      path to the file from which this object was created
    }
  }
}
\section{Methods}{
  \describe{
    \item{\link{findPeaks}}{
      \code{signature(object = "xcmsRaw")}: feature detection using
      matched filtration in the chromatographic time domain
    }
    \item{\link{getEIC}}{
      \code{signature(object = "xcmsRaw")}: get extracted ion
      chromatograms in specified m/z ranges
    }
    \item{\link{getPeaks}}{
      \code{signature(object = "xcmsRaw")}: get data for peaks in
      specified m/z and time ranges
    }
    \item{\link{getScan}}{
      \code{signature(object = "xcmsRaw")}: get m/z and intensity
      values for a single mass scan
    }
    \item{\link{getSpec}}{
      \code{signature(object = "xcmsRaw")}: get average m/z and
      intensity values for multiple mass scans
    }
    \item{\link[xcms:image.xcmsRaw]{image}}{
      \code{signature(x = "xcmsRaw")}: get data for peaks in
      specified m/z and time ranges
    }
    \item{\link{plotChrom}}{
      \code{signature(object = "xcmsRaw")}: plot a chromatogram
      from profile data
    }
    \item{\link{plotRaw}}{
      \code{signature(object = "xcmsRaw")}: plot locations of raw
      intensity data points
    }
    \item{\link{plotScan}}{
      \code{signature(object = "xcmsRaw")}: plot a mass spectrum
      of an individual scan from the raw data
    }
    \item{\link{plotSpec}}{
      \code{signature(object = "xcmsRaw")}: plot a mass spectrum
      from profile data
    }
    \item{\link{plotSurf}}{
      \code{signature(object = "xcmsRaw")}: experimental method for
      plotting 3D surface of profile data with \code{rgl}.
    }
    \item{\link{plotTIC}}{
      \code{signature(object = "xcmsRaw")}: plot total ion count
      chromatogram
    }
    \item{\link{profMedFilt}}{
      \code{signature(object = "xcmsRaw")}: median filter profile
      data in time and m/z dimensions
    }
    \item{\link{profMethod<-}}{
      \code{signature(object = "xcmsRaw")}: change the method of
      generating the \code{profile} matrix
    }
    \item{\link{profMethod}}{
      \code{signature(object = "xcmsRaw")}: get the method of
      generating the \code{profile} matrix
    }
    \item{profMz}{
      \code{signature(object = "xcmsRaw")}: get vector of m/z values
      for each row of the \code{profile} matrix
    }
    \item{\link{profRange}}{
      \code{signature(object = "xcmsRaw")}: interpret flexible ways
      of specifying subsets of the \code{profile} matrix
    }
    \item{\link{profStep<-}}{
      \code{signature(object = "xcmsRaw")}: change the m/z step
      used for generating the \code{profile} matrix
    }
    \item{\link{profStep}}{
      \code{signature(object = "xcmsRaw")}: get the m/z step used
      for generating the \code{profile} matrix
    }
    \item{\link{rawMat}}{
      \code{signature(object = "xcmsRaw")}: get a matrix of the raw data
      (columns time, mz and intensity)
    }
    \item{revMz}{
      \code{signature(object = "xcmsRaw")}: reverse the order of the
      data points for each scan
    }
    \item{sortMz}{
      \code{signature(object = "xcmsRaw")}: sort the data points
      by increasing m/z for each scan
    }
    \item{write.cdf}{
      \code{signature(object = "xcmsRaw")}: write the xcmsRaw object to
      a NetCDF file. This is useful e.g. for converting the XML formats
      mzData/mzXML/mzML to CDF.
    }
  }
}
\references{
  A parallel effort in metabolite profiling data sharing:
  \url{http://metlin.scripps.edu/}
}
\author{Colin A. Smith, \email{csmith@scripps.edu}}
\note{
  No notes yet.
}
\seealso{
  \code{\link{xcmsRaw}}
}
\keyword{classes}
