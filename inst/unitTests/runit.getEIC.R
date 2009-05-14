test.getEICxraw <- function() {
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xraw <- xcmsRaw(file)
    e <- getEIC(xraw, rtrange=cbind(3000,3500), mzrange=cbind(200,201))
}

test.getEICxset <- function() {
    xset <- fillPeaks(group(faahko))
    e <- getEIC(xset, sampleidx=c(1,2), groupidx=c(1,2), rtrange=200)
    plot(e)
}

test.plotEIC <- function() {
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xraw <- xcmsRaw(file)
    e <- getEIC(xraw, rtrange=cbind(3000,3500), mzrange=cbind(200,201))
    plot(e)
}
