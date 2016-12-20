testFilledFlag <- function() {
    xsg <- group(faahko)
    xsgf <- fillPeaks(xsg, method="chrom")

    checkEqualsNumeric(nrow(peaks(xsg)) + length(xsgf@filled), nrow(peaks(xsgf)))
}

testFillPeaksPar <- function() {
    xsg <- group(faahko)

    xsgfSerial <- fillPeaks(xsg, method="chrom")
    xsgfParallel <- fillPeaks(xsg, method="chrom") # parallel disabled: , nSlaves=2)

    checkEqualsNumeric(nrow(peaks(xsgfSerial)),
                       nrow(peaks(xsgfParallel)))
}



test.fillPeaksColumns <- function() {
  xsg <- group(faahko)
  peaks(xsg) <- cbind(peaks(xsg), anotherColumn=4711)

  oldCnames <- colnames(peaks(xsg))
  xsgf <- fillPeaks(xsg) # parallel disabled: , nSlaves=2)

  newCnames <- colnames(peaks(xsgf))
  checkEquals(oldCnames, newCnames)

  ## Check dims if nothing to do
  oldDims <- dim(peaks(xsgf))
  xsgf2 <- fillPeaks(xsgf) # parallel disabled: , nSlaves=2)
  newDims <- dim(peaks(xsgf2))
  checkEquals(oldDims, newDims)

  ## Case where only some samples have NA values
  xsg <- group(faahko, minfrac=1)
  xsgf <- fillPeaks(xsg) # parallel disabled: , nSlaves=2)
  sampclass(xsgf) <- c(rep("KO", 5), rep("WT", 7))
  xsgf <- group(xsgf, minfrac=1)
  xsgf <- fillPeaks(xsgf) # parallel disabled: , nSlaves=2)

}


## testFilledFlagMSW <- function() {

##   xsg <- group(ham)
##   xsgf <- fillPeaks(xsg)

##   checkEqualsNumeric(nrow(peaks(xsg)) + length(xsgf@filled, nrow(peaks(xsgf))) )
## }
