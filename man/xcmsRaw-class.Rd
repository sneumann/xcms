\name{xcmsRaw-class}
\docType{class}
\alias{xcmsRaw-class}
\alias{mslevel,xcmsRaw-method}
\alias{profinfo,xcmsRaw-method}
\alias{profMz,xcmsRaw-method}
\alias{profMz}
\alias{revMz,xcmsRaw-method}
\alias{revMz}
\alias{scanrange,xcmsRaw-method}
\alias{show,xcmsRaw-method}
\alias{sortMz,xcmsRaw-method}
\alias{sortMz}
\alias{levelplot}
\alias{findmzROI}
\alias{findmzROI,xcmsRaw-method}

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
    \item{\code{acquisitionNum}:}{Numeric representing the acquisition
      number of the individual scans/spectra. Length of
      \code{acquisitionNum} is equal to the number of spectra/scans in the
      object and hence equal to the \code{scantime} slot. Note however that
      this information is only available in mzML files.}
    \item{\code{env}:}{
      environment with three variables: \code{mz} - concatenated
      m/z values for all scans, \code{intensity} - corresponding
      signal intensity for each m/z value, and \code{profile} -
      matrix represention of the intensity values with columns
      representing scans and rows representing equally spaced m/z
      values. The profile matrix should be extracted with the
      \code{\link{profMat}} method.
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
    \item{\code{msnPrecursorIntensity}:}{"Intensity" (mzData) or
      "precursorIntensity" (mzXML)}
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
      profile matrix.
    }
    \item{\code{profparam}:}{list to store additional profile matrix
      generation settings. Use the \code{\link{profinfo}} method to
      extract all profile matrix creation relevant information.}
    \item{\code{scanindex}:}{
      integer vector with starting positions of each scan in the
      \code{mz} and \code{intensity} variables (note that index
      values are based off a 0 initial position instead of 1).
    }
    \item{\code{scantime}:}{
      numeric vector with acquisition time (in seconds) for each
      scan.
    }
    \item{\code{tic}:}{
      numeric vector with total ion count (intensity) for each scan
    }
    \item{\code{mslevel}:}{
      Numeric representing the MS level that is present in MS1
      slot. This slot should be accessed through its getter method
      \code{mslevel}.
    }
    \item{\code{scanrange}:}{
      Numeric of length 2 specifying the scan range (or \code{NULL} for
      the full range). This slot should be accessed through its getter
      method \code{scanrange}. Note that the \code{scanrange} will
      always be 1 to the number of scans within the \code{xcmsRaw}
      object, which does not necessarily have to match to the scan index in
      the original mzML file (e.g. if the original data was sub-setted). The
      \code{acquisitionNum} information can be used to track the
      original \emph{position} of each scan in the mzML file.
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
      chromatograms in specified m/z ranges. This will return the total
  ion chromatogram (TIC) if the m/z range corresponds to the full m/z
  range (i.e. sum of all signals per retention time across all m/z).
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
    \item{\link[xcms:levelplot.xcmsRaw]{levelplot}}{
      Create an image of the raw (profile) data m/z against retention
  time, with the intensity color coded.
    }
    \item{mslevel}{
      Getter method for the \code{mslevel} slot.
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
    \item{\link{profinfo}}{
      \code{signature(object = "xcmsRaw")}: returns a list containing
  the profile generation method and step (profile m/z step size) and
  eventual additional parameters to the profile function.
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
    \item{scanrange}{
      Getter method for the \code{scanrange} slot. See slot description
      above for more information.
    }
    \item{sortMz}{
      \code{signature(object = "xcmsRaw")}: sort the data points
      by increasing m/z for each scan
    }
    \item{stitch}{
      \code{signature(object = "xcmsRaw")}: Raw data correction for
      lock mass calibration gaps.
    }
    \item{findmzROI}{
      \code{signature(object = "xcmsRaw")}:
	internal function to identify regions of interest in the raw
	data as part of the first step of centWave-based peak detection.
      }
    }
}
\author{Colin A. Smith, \email{csmith@scripps.edu},
  Johannes Rainer \email{johannes.rainer@eurac.edu}}
\seealso{
  \code{\link{xcmsRaw}}, \code{\link{subset-xcmsRaw}} for subsetting by spectra.
}
\keyword{classes}
