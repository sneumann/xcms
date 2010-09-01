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

getMsnScan <- function(object, scan, mzrange = numeric()) {

    if (scan < 0)
        scan <- length(object@msnRt) + 1 + scan ## shouldn't that be -scan ? 

    idx <- seq(object@msnScanindex[scan]+1, min(object@msnScanindex[scan+1],
                                             length(object@env$msnMz), na.rm=TRUE))

    if (length(mzrange) >= 2) {
        mzrange <- range(mzrange)
        idx <- idx[object@env$msnMz[idx] >= mzrange[1] & object@env$msnMz[idx] <= mzrange[2]]
    }

    points <- cbind(mz = object@env$msnMz[idx], intensity = object@env$msnIntensity[idx])

    invisible(points)
}

buildMzdata <- function(xr) {

  mslevel <- 1; # This can be enhanced to MSn one day
  
  mzdata = xmlTree(tag="mzData",
    attrs=c(
      version="1.05",
      "xsi:noNamespaceSchemaLocation"="http://psidev.sourceforge.net/ms/xml/mzdata/mzdata.xsd"),
    namespaces = c(
      xsi="http://www.w3.org/2001/XMLSchema-instance")
    )

  mzdata$addNode("spectrumList",
                 attrs=c(count=length(xr@scanindex-1)),
                 close = FALSE)

  for (id in 1:length(xr@scanindex-1)) {

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

    if (mslevel>1) { ## only for mslevel >1
      mzdata$addNode("precursorList",
                     attrs=c(count=1),
                     close = FALSE)
      
      mzdata$addNode("precursor",
                     attrs=c(msLevel=1, spectrumRef=0),
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
    peaks <- getScan(xr, id)
    
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
  mzdata$closeTag() ## </spectrumList>
  mzdata$closeTag() ## </mzData>

  mzdata
}
