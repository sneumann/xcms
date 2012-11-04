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
    \item{\code{acquisitionNum}:}{acquisitionNum}
    \item{\code{env}:}{
      environment with three variables: \code{mz} - concatenated
      m/z values for all scans, \code{intensity} - corresponding
      signal intensity for each m/z value, and \code{profile} -
      matrix represention of the intensity values with columns
      representing scans and rows representing equally spaced m/z
      values
    }
    \item{\code{filepath}:}{Path to the raw data file}
    \item{\code{gradient}:}{
      matrix with first row, \code{time}, containing the time point
      for interpolation and successive columns representing solvent
      fractions at each point
    }
    \item{\code{msnAcquisitionNum}:}{for each scan a unique acquisition
  number as reported via "spectrum id" (mzData) or "<scan num=...>" and
  "<scanOrigin num=...>" (mzXML)}
    \item{\code{msnCollisionEnergy}:}{"CollisionEnergy" (mzData) or
  "collisionEnergy" (mzXML)}
    \item{\code{msnLevel}:}{for each scan the "msLevel" (both mzData and mzXML)}
%    \item{\code{msnPeakCount}:}{msnPeakCount is equal to scan length} OBSOLETED
    \item{\code{msnPrecursorCharge}:}{"ChargeState" (mzData) and
  "precursorCharge" (mzXML)}
    \item{\code{msnPrecursorIntensity}:}{"Intensity" (mzData) or "precursorIntensity" (mzXML)}
    \item{\code{msnPrecursorMz}:}{"MassToChargeRatio" (mzData) or
  "precursorMz" (mzXML)}
    \item{\code{msnPrecursorScan}:}{"spectrumRef" (both mzData and mzXML)}
    \item{\code{msnRt}:}{Retention time of the scan}
    \item{\code{msnScanindex}:}{msnScanindex}
    \item{\code{mzrange}:}{
      numeric vector of length 2 with minimum and maximum m/z values
      represented in the profile matrix
    }
    \item{\code{polarity}:}{polarity}
    \item{\code{profmethod}:}{
      characer value with name of method used for generating the
      profile matrix
    }
    \item{\code{profparam}:}{profparam}
    \item{\code{scanindex}:}{
      integer vector with starting positions of each scan in the
      \code{mz} and \code{intensity} variables (note that index
      values are based off a 0 initial position instead of 1)
    }
    \item{\code{scantime}:}{
      numeric vector with acquisition time (in seconds) for each
      scan
    }
    \item{\code{tic}:}{
      numeric vector with total ion count (intensity) for each scan
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
    \item{revMz}{
      \code{signature(object = "xcmsRaw")}: reverse the order of the
      data points for each scan
    }
    \item{sortMz}{
      \code{signature(object = "xcmsRaw")}: sort the data points
      by increasing m/z for each scan
	}
	\item{stitch}{
	  \code{signature(object = "xcmsRaw")}: Raw data correction for
	  lock mass calibration gaps.
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
