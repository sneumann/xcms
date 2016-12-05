## Splitting of xcmsRaw objects.

file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
xraw <- xcmsRaw(file, profstep = 0)

testSplitRawEven <- function() {
    xrl <- split(xraw, f = xraw@scanindex%%2)
    checkEqualsNumeric(length(xrl), 2)

    checkTrue(length(xraw@scanindex) == length(xrl[[1]]@scanindex)
              + length(xrl[[2]]@scanindex))
    checkTrue(length(xraw@env$mz) == length(xrl[[1]]@env$mz)
              + length(xrl[[2]]@env$mz))
}

testSplitRawOdd <- function() {
    xrl <- split(xraw, f = (xraw@scanindex+1)%%2)
    checkEqualsNumeric(length(xrl), 2)

    checkTrue(length(xraw@scanindex) == length(xrl[[1]]@scanindex) +
              length(xrl[[2]]@scanindex))
    checkTrue(length(xraw@env$mz) == length(xrl[[1]]@env$mz) +
              length(xrl[[2]]@env$mz))
}

testSplitRawFirst <- function() {
    xrl <- split(xraw, f=c(1,rep(2,length(xraw@scanindex)-1)))

    checkEqualsNumeric(length(xrl), 2)

}

testSplitRawLast <- function() {
    xrl <- split(xraw, f=c(rep(1,length(xraw@scanindex)-1), 2))

    checkEqualsNumeric(length(xrl), 2)
}

testSplitRawNone <- function() {
    xrl <- split(xraw, f=rep(1,length(xraw@scanindex)))

    checkEqualsNumeric(length(xrl), 1)

    xraw2 <- xrl[[1]]
    checkTrue(length(xraw@scanindex) == length(xrl[[1]]@scanindex))
    checkTrue(all(xraw@scanindex == xraw2@scanindex))

}

############################################################
## Subset an xcmsRaw object by scan index.
test_bracket_subset_xcmsRaw <- function() {
    ## Get scans 1:10
    xsub <- xraw[1:10]
    checkIdentical(xraw@scantime[1:10], xsub@scantime)
    checkIdentical(xraw@scanindex[1:10], xsub@scanindex[1:10])
    checkIdentical(xraw@env$mz[1:xraw@scanindex[11]], xsub@env$mz)
    checkIdentical(xraw@env$intensity[1:xraw@scanindex[11]], xsub@env$intensity)

    ## Check if mz is sorted:

    ## Check using logical
    bm <- rep(FALSE, length(xraw@scanindex))
    bm[1:10] <- TRUE
    xsub_2 <- xraw[bm]
    checkEquals(xsub, xsub_2)

    ## Get none.
    xempty <- xraw[rep(FALSE, length(xraw@scanindex))]
    checkTrue(length(xempty@scanindex) == 0)

    ## Get the full one:
    xsub <- xraw[1:length(xraw@scanindex)]
    checkEquals(xsub@env$mz, xraw@env$mz)
    checkEquals(xsub@env$intensity, xraw@env$intensity)

    ## Get some scans in the middle somewhere
    i <- c(5, 99, 317)
    xsub <- xraw[i]
    checkIdentical(xsub@scantime, xraw@scantime[i])
    ## scanindex:
    vps <- diff(c(xraw@scanindex, length(xraw@env$mz)))
    scnidx <- xcms:::valueCount2ScanIndex(vps[i])
    checkEquals(scnidx, xsub@scanindex)
    whichIdx <- c(((xraw@scanindex[5] + 1):xraw@scanindex[6]),
    ((xraw@scanindex[99] + 1):xraw@scanindex[100]),
    ((xraw@scanindex[317] + 1):xraw@scanindex[318]))
    checkIdentical(xsub@env$mz, xraw@env$mz[whichIdx])
    checkIdentical(xsub@env$intensity, xraw@env$intensity[whichIdx])

    ## Finally check that we get the same object by calling getXcmsRaw from an
    ## xcmsSet and by subsetting or using xcmsRaw with scanrange:
    xraw <- xcmsRaw(file, profstep = 20)
    xraw_sub <- xraw[5:100]
    suppressWarnings(
        xset <- xcmsSet(file, scanrange = c(5, 100), step = 20)
    )
    xraw_xset <- getXcmsRaw(xset)
    checkEquals(scanrange(xraw_sub), scanrange(xraw_xset))
    ## Load the object using xcmsRaw and scanrange
    xraw_2 <- xcmsRaw(file, scanrange = c(5, 100), profstep = 20)
    checkEquals(scanrange(xraw_sub), scanrange(xraw_2))
    ## Compare objects
    checkEquals(xraw_sub, xraw_xset)
    checkEquals(xraw_sub, xraw_2)
    checkEquals(xraw_xset, xraw_2)
}
