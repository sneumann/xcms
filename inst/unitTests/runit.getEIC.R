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

test.getEICretcor <- function() {
    xset <- fillPeaks(group(retcor(group(faahko))))
    opt.warn <- options("warn")$warn
    options("warn" = 2) ## turns warning into errors
    e <- getEIC(xset, sampleidx=c(1,2), groupidx=c(1,2),
                rt="corrected", rtrange=200)
    options("warn" = opt.warn)
    plot(e)
}

test.plotEIC <- function() {
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xraw <- xcmsRaw(file)
    e <- getEIC(xraw, rtrange=cbind(3000,3500), mzrange=cbind(200,201))
    plot(e)
}
