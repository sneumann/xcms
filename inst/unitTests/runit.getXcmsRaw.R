## just plain function that reads the raw data...
test.getXcmsRaw <- function(){

    ## xsetRaw <- updateObject(faahko)
    ## xset <- fillPeaks(group(retcor(group(xsetRaw))))
    xset <- faahko_grouped_retcor_filled
    
    ## get the first as raw data file.
    xr <- getXcmsRaw(xset, sampleidx = 1)
    ## apply the rt correction
    xr <- getXcmsRaw(xset, sampleidx = 1, rt="raw")
    ## get the first 4
    xr <- getXcmsRaw(xset, sampleidx = 1:4)
    ## check if the settings are translated correctly
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xs <- xcmsSet(file, step = 2)
    xr <- getXcmsRaw(xs)
    xr_orig <- xcmsRaw(file, profstep = 2)
    checkEquals(xr, xr_orig)
    ## check the prof step:
    checkEqualsNumeric(profStep(xs), 2)
    ## check all *new* methods for xcmsSet
    checkEquals(mslevel(xs), mslevel(xr))
    checkEquals(profMethod(xs), profMethod(xr))
    checkEquals(profStep(xs), profStep(xr))
    ## scanrange for the xcmsSet is NULL which means we're reading all data from
    ## the raw data files, while the one of the xcmsRaw is always
    ## (1, lenght(scans)).
    ## checkEquals(scanrange(xs), scanrange(xr))
    profinfo(xs)
    profinfo(xr)
    ## testing alternative scan range.
    xr2 <- getXcmsRaw(xs, scanrange = c(5, 100))
    scanrange(xr2)
    xr2_orig <- xcmsRaw(file, scanrange = c(5, 100), profstep = 2)
    checkEquals(xr2, xr2_orig)
    ## This scanrange is expected to be from 1 to length(xr@scantime)
    checkEquals(scanrange(xr2), c(1, length(xr2@scantime)))
    ## Test xcmsSet with scanrange:
    xs <- xcmsSet(file, step = 2, scanrange = c(5, 100))
    checkEquals(scanrange(xs), c(5, 100))
    ## BUT: if we extract the xcmsRaw from this xcmsSet object we will get
    ## (1, 96) instead of (5, 100).
    xr <- getXcmsRaw(xs)
    checkEquals(scanrange(xr), c(1, length(xr@scantime)))
    xr_2 <- xcmsRaw(file, scanrange = c(5, 100), profstep = 2)
    checkEquals(xr, xr_2)
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
