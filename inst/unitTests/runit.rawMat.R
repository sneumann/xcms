test.rawMat <- function() {
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xraw <- xcmsRaw(file)
    rawmat <- rawMat(xraw, c(200,300), c(2500,3000))
    checkEqualsNumeric(nrow(rawmat),31770)
}
