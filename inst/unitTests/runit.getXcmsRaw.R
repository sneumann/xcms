## just plain function that reads the raw data...
test.getXcmsRaw <- function(){

    xsetRaw <- updateObject(faahko)
    xset <- fillPeaks(group(retcor(group(xsetRaw))))

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


test.getXcmsRawIssue44 <- function() {


    ## Subset to two files.
    xsetRaw <- updateObject(faahko)
    xsetRaw <- xsetRaw[, 1:2]

    ## First sample is reference, i.e. no rt adjustment performed
    xs <- retcor(group(xsetRaw), method = "obiwarp", center = 1)
    ## Second is corrected, first is center:
    checkIdentical(xs@rt$raw[[1]], xs@rt$corrected[[1]])
    checkTrue(!all(xs@rt$raw[[2]] == xs@rt$corrected[[2]]))

    ## Now, if we get the second raw file we expect to get the raw times.
    all(xs@rt$corrected[[1]] == xs@rt$raw[[1]]) ## TRUE, wouldn't do correction.
    all(xs@rt$corrected[[2]] == xs@rt$raw[[2]]) ## FALSE, would do correction.

    ## We get the raw data for the second file; this one was corrected and
    ## thus it's retention time is expected to be different from raw.
    xr2 <- getXcmsRaw(xs, sampleidx = 2, rt = "corrected")
    checkIdentical(xr2@scantime, xs@rt$corrected[[2]])

    all(xr2@scantime == xs@rt$raw[[2]])  ## That should be FALSE!

}
