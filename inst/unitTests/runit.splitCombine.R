testSplit <- function() {
    xsl <- split(faahko,sampclass(faahko))
    checkEqualsNumeric(length(xsl), 2)
}

testSplitAll <- function() {
    xsl <- split(faahko,sampnames(faahko))
    checkEqualsNumeric(length(xsl), length(sampnames(faahko)))
}

testSplitFirst <- function() {
    xsl <- split(faahko,c(1,2,2,2,2,2,2,2,2,2,2,2))
    checkEqualsNumeric(length(xsl), 2)
}
testSplitLast <- function() {
    xsl <- split(faahko,c(2,2,2,2,2,2,2,2,2,2,2,1))
    checkEqualsNumeric(length(xsl), 2)
}
testSplitMiddle <- function() {
    xsl <- split(faahko,c(2,2,2,2,2,1,2,2,2,2,2,2))
    checkEqualsNumeric(length(xsl), 2)
}
testSplitNone <- function() {
    xsl <- split(faahko,c(2,2,2,2,2,2,2,2,2,2,2,2))
    checkEqualsNumeric(length(xsl), 1)
}


testCombine <- function() {
    xsl <- split(faahko,sampclass(faahko))
    checkEqualsNumeric(length(sampnames(c(xsl[[1]], xsl[[2]]))), 12)
}


testSplitFactorShort = function() {
  f = c(1,2,2)
  xsl <- split(faahko, f)
  for (x in unique(f)) {
    num.samps = sum(x==f)
    checkEqualsNumeric(length(xsl[[x]]@filepaths), num.samps)
    checkEqualsNumeric(nrow(xsl[[x]]@phenoData), num.samps)
    checkEqualsNumeric(length(xsl[[x]]@rt$raw), num.samps)
  }
}

testSplitFactorLongDrop = function() {
  f = c(1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,3)
  actual.f = f[1:nrow(faahko@phenoData)]
  xsl <- split(faahko, f)
  for (x in unique(actual.f)) {
    num.samps = sum(x==actual.f)
    checkEqualsNumeric(length(xsl[[x]]@filepaths), num.samps)
    checkEqualsNumeric(nrow(xsl[[x]]@phenoData), num.samps)
    checkEqualsNumeric(length(xsl[[x]]@rt$raw), num.samps)
  }
}