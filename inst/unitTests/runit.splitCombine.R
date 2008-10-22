testSplit <- function() {
    xsl <- split(faahko,sampclass(faahko))
    checkEqualsNumeric(length(xsl), 2)
}

testCombine <- function() {
    xsl <- split(faahko,sampclass(faahko))
    checkEqualsNumeric(length(sampnames(c(xsl[[1]], xsl[[2]]))), 12)
}
