## just plain function that reads the raw data...
test.getXcmsRaw <- function(){
    xset <- fillPeaks(group(retcor(group(faahko))))
    ## get the first as raw data file.
    xr <- getXcmsRaw(xset, sampleidx=1)
    ## apply the rt correction
    xr <- getXcmsRaw(xset, sampleidx=1, rt="raw")
    ## get the first 4
    xr <- getXcmsRaw(xset, sampleidx=1:4)
}
