## just plain function that reads the raw data...
test.getXcmsRaw <- function(){
    xset <- fillPeaks(group(retcor(group(faahko))))
    ## get the first as raw data file.
    xr <- getXcmsRaw(xset, sampleidx=1)
    ## apply the rt correction
    xr <- getXcmsRaw(xset, sampleidx=1, rt="raw")
    ## get the first 4
    xr <- getXcmsRaw(xset, sampleidx=1:4)
    ## check if the settings are translated correctly
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xs <- xcmsSet(file, step=2)
    xr <- getXcmsRaw(xs)
    ## check the prof step:
    checkEqualsNumeric(profStep(xs), 2)
    ## check all *new* methods for xcmsSet
    checkEquals(mslevel(xs), mslevel(xr))
    checkEquals(profMethod(xs), profMethod(xr))
    checkEquals(profStep(xs), profStep(xr))
    checkEquals(scanrange(xs), scanrange(xr))
    profinfo(xs)
    profinfo(xr)
    ## testing alternative scan range.
    xr2 <- getXcmsRaw(xs, scanrange=c(5, 100))
    scanrange(xr2)
    checkEquals(scanrange(xr2), c(5, 100))
}



