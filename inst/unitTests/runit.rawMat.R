############################################################
## Test the rawMat method.
test.rawMat <- function() {
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xraw <- xcmsRaw(file, profstep = 0)
    ## That's the old test; pretty arbitrary
    ##rawmat <- rawMat(xraw, mzrange = c(200,300), rtrange = c(2500,3000))
    ##checkEqualsNumeric(nrow(rawmat),31770)

    ## New evaluations:
    ## Check with manually generated content:
    scnt <- xraw@scantime
    scnrange <- range(which((scnt >= 2500) & (scnt <= 3000)))
    ## Just to ensure: are the scantimes within our range?
    checkTrue(all((scnt[scnrange] >= 2500 & scnt[scnrange] <= 3000)))
    ## Now get the values falling into this scan range:
    ## scanindex records the start index for each scan/spectrum -1.
    scanStart <- xraw@scanindex[scnrange[1]] + 1
    ## To get the end index we add also the total number of scans and
    ## pick the start index of the next scan outside our scnrange
    scanEnd <- c(xraw@scanindex, length(xraw@env$mz))[scnrange[2] + 1]
    mzVals <- xraw@env$mz[scanStart:scanEnd]
    ## Check if that's what we get from the function.
    res_scnr <- rawMat(xraw, scanrange = scnrange)
    checkEquals(res_scnr[, "mz"], mzVals)
    ## Is that the same if we used the rtrange instead?
    res_rtr <- rawMat(xraw, rtrange = c(2500, 3000))
    checkIdentical(res_rtr[, "mz"], mzVals)
    ## Hm. Actually we have a discrepancy here:
    checkIdentical(res_scnr, res_rtr)
    ## Compare the retention time for the subset: get the retention time
    ## for each scan and "rep" it the number of values per scan.
    valsPerSpect <- diff(c(xraw@scanindex, length(xraw@env$mz)))
    rtimes <- rep(scnt[scnrange[1]:scnrange[2]],
                  valsPerSpect[scnrange[1]:scnrange[2]])
    checkIdentical(res_scnr[, "time"], rtimes)
    ## This failing test was reported as issue #58 on github.

    ## Add also an mzrange:
    mzr <- c(300, 330)
    mzInRange <- xraw@env$mz >= mzr[1] & xraw@env$mz <= mzr[2]
    res_mzr <- rawMat(xraw, mzrange = mzr)
    checkEquals(res_mzr[, "mz"], xraw@env$mz[mzInRange])
    ## Combine with scnrange:
    scnInRange <- rep(FALSE, length(mzInRange))
    scnInRange[scanStart:scanEnd] <- TRUE
    res_mzr_scnr <- rawMat(xraw, scanrange = scnrange, mzrange = mzr)
    checkIdentical(res_mzr_scnr[, "mz"], xraw@env$mz[mzInRange & scnInRange])
    ## And the retention time.
    allRtimes <- rep(xraw@scantime, diff(c(xraw@scanindex, length(xraw@env$mz))))
    checkIdentical(res_mzr_scnr[, "time"], allRtimes[mzInRange & scnInRange])

    ## And a second scan range:
    scnrange <- c(13, 54)
    scanStart <- xraw@scanindex[scnrange[1]] + 1
    scanEnd <- c(xraw@scanindex, length(xraw@env$mz))[scnrange[2] + 1]
    mzVals <- xraw@env$mz[scanStart:scanEnd]
    res_scnr <- rawMat(xraw, scanrange = scnrange)
    checkIdentical(res_scnr[, "mz"], mzVals)
    valsPerSpect <- diff(c(xraw@scanindex, length(xraw@env$mz)))
    rtimes <- rep(scnt[scnrange[1]:scnrange[2]],
                  valsPerSpect[scnrange[1]:scnrange[2]])
    checkIdentical(res_scnr[, "time"], rtimes)
}

