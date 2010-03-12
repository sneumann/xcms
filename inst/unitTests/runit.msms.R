test.xcmsRawms2 <- function() {


  filename <- system.file('iontrap/extracted.mzData', package = "msdata")

  x<-xcmsRaw(filename, includeMSn=TRUE)

  x2<-xcmsRaw(filename, includeMSn=TRUE, mslevel=2)

  ## STN: add tests
}


##library(xcms);library(faahKO);library(msdata);filename <- "/vol/R/R-devel/lib64/R/library/msdata/iontrap/extracted.mzData";xs <- xcmsSet(files=filename, mslevel=2)

