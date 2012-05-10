setGeneric("write.mzdata", function(object, ...) standardGeneric("write.mzdata"))

setMethod("write.mzdata", "xcmsRaw", function(object, filename) {

  require(XML) || stop("We need library(XML) to write mzData")
  require(caTools) || stop("We need library(caTools) to encode base64 values in mzData")

  mzdata <- buildMzdata(object)

  ## the sink() workaround seems to be needed for proper indenting.
  sink(file=filename)
  cat(saveXML(mzdata))
  sink(NULL)
})

setGeneric("getMsnScan", function(object, ...) standardGeneric("getMsnScan"))
setMethod("getMsnScan", "xcmsRaw", function(object, scan, mzrange = numeric()) {

    if (scan < 0)
        scan <- length(object@msnRt) + 1 + scan

    idx <- seq(object@msnScanindex[scan]+1, min(object@msnScanindex[scan+1],
                                             length(object@env$msnMz), na.rm=TRUE))

    if (length(mzrange) >= 2) {
        mzrange <- range(mzrange)
        idx <- idx[object@env$msnMz[idx] >= mzrange[1] & object@env$msnMz[idx] <= mzrange[2]]
    }

    points <- cbind(mz = object@env$msnMz[idx], intensity = object@env$msnIntensity[idx])

    invisible(points)
})

buildMzdata <- function(xr) {

  mzdata <- xmlTree(namespaces = c(xsi="http://www.w3.org/2001/XMLSchema-instance"))
  
  mzdata$addNode("mzData",
                 attrs=c(
                   version="1.05",
                   "xsi:noNamespaceSchemaLocation"="http://psidev.sourceforge.net/ms/xml/mzdata/mzdata.xsd"),
                 close=FALSE)
  
  mzdata$addNode("spectrumList",
                 attrs=c(count=length(xr@scanindex)),
                 close = FALSE)

  mslevel <- 1; # Write MS1 first
  offset <- 0;
  
  if (length(xr@scanindex) > 0) {
  for (id in 1:length(xr@scanindex)) {

    mzdata$addNode("spectrum",
                   attrs=c(id=id),
                   close = FALSE)
    mzdata$addNode("spectrumDesc",
                   attrs=c(),
                   close = FALSE)
    mzdata$addNode("spectrumSettings",
                   attrs=c(),
                   close = FALSE)
    mzdata$addNode("spectrumInstrument",
                   attrs=c(msLevel="1"),
                   close = FALSE)

    mzdata$addNode("cvParam",
                   attrs=c(
                     cvLabel="PSI", accession="PSI:1000039",
                     name="TimeInSeconds", value=xr@scantime[id]))
    mzdata$closeTag() ## </spectrumInstrument>
    mzdata$closeTag() ## </spectrumSettings>

    mzdata$closeTag() ## </spectrumDesc>

    target <- new("raw")
    peaks <- getScan(xr, id)
    
    if (is.unsorted(peaks[,"mz"])) { ## fix "bad" scans
        peaks <- peaks[order(peaks[,"mz"]),]
    }    
    
    mzdata$addNode("mzArrayBinary",
                   close = FALSE)
    mzdata$addNode("data",
                   attrs=c(
                     precision="32", endian="big", length=nrow(peaks)),
                   value=base64encode(writeBin(peaks[,"mz"], con=target, size=4, endian="big")))

    mzdata$closeTag() ## </mzArrayBinary>

    mzdata$addNode("intenArrayBinary",
                   close = FALSE)
    mzdata$addNode("data",
                   attrs=c(
                     precision="32", endian="big", length=nrow(peaks)),
                   value=base64encode(writeBin(peaks[,"intensity"], con=target, size=4, endian="big")))
    mzdata$closeTag() ## </intenArrayBinary>

    mzdata$closeTag() ## </spectrum>
    offset <- id  
    }  
}
  
  mslevel <- xr@msnLevel
  if (length(xr@msnScanindex) >0 ) {
  for (id in seq(1, length(xr@msnScanindex))) {

    mzdata$addNode("spectrum",
                   attrs=c(id=id+offset),
                   close = FALSE)
    mzdata$addNode("spectrumDesc",
                   attrs=c(),
                   close = FALSE)
    mzdata$addNode("spectrumSettings",
                   attrs=c(),
                   close = FALSE)
    mzdata$addNode("spectrumInstrument",
                   attrs=c(msLevel=mslevel[id]),
                   close = FALSE)

    mzdata$addNode("cvParam",
                   attrs=c(
                     cvLabel="PSI", accession="PSI:1000039",
                     name="TimeInSeconds", value=xr@msnRt[id]))
    mzdata$closeTag() ## </spectrumInstrument>
    mzdata$closeTag() ## </spectrumSettings>

    if (mslevel[id] > 1) { ## only for mslevel >1, here always true
      mzdata$addNode("precursorList",
                     attrs=c(count=1),
                     close = FALSE)
      
      mzdata$addNode("precursor",
                     attrs=c(msLevel=mslevel[id]-1, spectrumRef=max(0, xr@msnPrecursorScan[id], na.rm=T)),
                     close = FALSE)
      
      mzdata$addNode("ionSelection",
                     close = FALSE)
      mzdata$addNode("cvParam",
                     attrs=c(
                       cvLabel="PSI", accession="PSI:1000040",
                       name="m/z", value=xr@msnPrecursorMz[id] ))
      mzdata$closeTag() ## </ionSelection>
      
      mzdata$addNode("activation",
                     close = FALSE)
      mzdata$addNode("cvParam",
                     attrs=c(
                       cvLabel="PSI", accession="PSI:1000045",
                       name="Collision Energy", value=xr@msnCollisionEnergy[id] ))
      mzdata$closeTag() ## </activation>
      
      mzdata$closeTag() ## </precursor>
      mzdata$closeTag() ## </precursorList>         
    }
    mzdata$closeTag() ## </spectrumDesc>

    target <- new("raw")
    peaks <- getMsnScan(xr, id)
    
    if (is.unsorted(peaks[,"mz"])) { ## fix "bad" scans
        peaks <- peaks[order(peaks[,"mz"]),]
    }    
    
    mzdata$addNode("mzArrayBinary",
                   close = FALSE)
    mzdata$addNode("data",
                   attrs=c(
                     precision="32", endian="big", length=nrow(peaks)),
                   value=base64encode(writeBin(peaks[,"mz"], con=target, size=4, endian="big")))

    mzdata$closeTag() ## </mzArrayBinary>

    mzdata$addNode("intenArrayBinary",
                   close = FALSE)
    mzdata$addNode("data",
                   attrs=c(
                     precision="32", endian="big", length=nrow(peaks)),
                   value=base64encode(writeBin(peaks[,"intensity"], con=target, size=4, endian="big")))
    mzdata$closeTag() ## </intenArrayBinary>

    mzdata$closeTag() ## </spectrum>
  }  
}
  
  mzdata$closeTag() ## </spectrumList>
  mzdata$closeTag() ## </mzData>

  mzdata
}
