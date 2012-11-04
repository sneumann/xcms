\name{findPeaks.MS1-methods}
\docType{methods}
\alias{findPeaks.MS1}
\alias{findPeaks.MS1,xcmsRaw-method}
\title{Collecting MS1 precursor peaks}
\description{
  Collecting Tandem MS or MS$^n$ Mass Spectrometry
  precursor peaks as annotated in XML raw file}
\section{Methods}{
\describe{
\item{object = "xcmsRaw"}{
  \code{
    findPeaks.MS1(object)
    }
}
}}
\details{
  Some mass spectrometers can acquire MS1 and MS2 (or MS$^n$ scans)
  quasi simultaneously, e.g. in data dependent tandem MS or DDIT mode.

  Since xcmsFragments attaches \emph{all} MS$^n$ peaks to MS1 peaks in xcmsSet,
  it is important that findPeaks and xcmsSet do not miss any MS1 precursor peak.

  To be sure that all MS1 precursor peaks are in an xcmsSet,
  findPeaks.MS1 does not do an actual peak picking, but simply uses the
  annotation stored in mzXML, mzData or mzML raw files.

  This relies on the following XML tags:

  mzData:
  \code{
  <spectrum id="463">
        <spectrumInstrument msLevel="2">
	  <cvParam cvLabel="psi" accession="PSI:1000039" name="TimeInSeconds" value="92.7743"/>
        </spectrumInstrument>
        <precursor msLevel="1" spectrumRef="461">
            <cvParam cvLabel="psi" accession="PSI:1000040" name="MassToChargeRatio" value="462.091"/>
            <cvParam cvLabel="psi" accession="PSI:1000042" name="Intensity" value="366.674"/>
        </precursor>
  </spectrum>
}

  mzXML:
  \code{
     <scan num="17" msLevel="2" retentionTime="PT1.5224S">
     <precursorMz precursorIntensity="125245">220.1828003</precursorMz>
     </scan>
   }

   Several mzXML and mzData converters are known to create incomplete
  files, either without intensities (they will be set to 0) or without
  the precursor retention time (then a reasonably close rt will be
  chosen. NYI).

}
\arguments{
  \item{object}{\code{xcmsRaw} object}
}
\value{
  A matrix with columns:

  \item{mz, mzmin, mzmax}{
    annotated MS1 precursor selection mass
  }

  \item{rt, rtmin, rtmax}{
    annotated MS1 precursor retention time
  }

  \item{into, maxo, sn}{
    annotated MS1 precursor intensity
  }

}
\author{Steffen Neumann, \email{sneumann@ipb-halle.de}}
\seealso{
  \code{\link{findPeaks-methods}}
  \code{\link{xcmsRaw-class}}
}
\keyword{methods}
