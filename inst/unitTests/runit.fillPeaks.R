testFilledFlag <- function() {
  xsg <- group(faahko)
  xsgf <- fillPeaks(xsg, method="chrom")

  checkEqualsNumeric(nrow(peaks(xsg)) + length(xsgf@filled), nrow(peaks(xsgf)))
}

## testFilledFlagMSW <- function() {
  
##   xsg <- group(ham)
##   xsgf <- fillPeaks(xsg)

##   checkEqualsNumeric(nrow(peaks(xsg)) + length(xsgf@filled, nrow(peaks(xsgf))) )
## }


