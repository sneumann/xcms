test.xcmsRawms1 <- function() {

    filename <- system.file('microtofq/MM14.mzdata', package = "msdata")

    ## This file has no MS/MS data at all, but should not fail
    x1 <- xcmsRaw(filename, includeMSn=TRUE, profstep = 0)

}

test.xcmsRawms123 <- function() {

    filename <- system.file('iontrap/extracted.mzData', package = "msdata")

    x1 <- xcmsRaw(filename, includeMSn=TRUE, profstep = 0)
    x2 <- xcmsRaw(filename, includeMSn=TRUE, mslevel=2, profstep = 0)
    x3 <- xcmsRaw(filename, includeMSn=TRUE, mslevel=3, profstep = 0)

    checkTrue(length(x1@env$msnMz) == length(x2@env$mz) + length(x3@env$mz))

    checkTrue(all(x1@msnLevel[1:6]==2))

    checkTrue(all(x1@msnScanindex[1:6] == x2@scanindex[1:6]))

    checkEqualsNumeric(nrow(getMsnScan(x1, scan=1)), 278)

    ## This would fail, since mslevel=2 above seems to use split(),
    ## which does drop MSn information ?
    ## checkEqualsNumeric(nrow(xcms:::getMsnScan(x2, scan=1)), 278)

}

test.xcmsSetms2mf <- function() {

    filename <- system.file('iontrap/extracted.mzData', package = "msdata")
    xs2 <- xcmsSet(filename, snthresh = 4, mslevel = 2)
}

test.xcmsSetms2cw <- function() {

    filename <- system.file('iontrap/extracted.mzData', package = "msdata")
    xs2 <- xcmsSet(filename, method="centWave", mslevel = 2)

}

test.msn2xcmsRaw <- function() {
 msnfile <- system.file("microtofq/MSMSpos20_6.mzML", package = "msdata")
 xrmsn <- xcmsRaw(msnfile, includeMSn=TRUE)
 xr <- msn2xcmsRaw(xrmsn)

 checkEqualsNumeric(length(xr@env$mz), 3132)
 checkEqualsNumeric(length(xr@env$intensity), 3132)
 ## In reality, it seems there are 1612 MS2 spectra in the file, just that
 ## 1121 have a peaksCount > 0
 checkEqualsNumeric(length(xr@scantime), 1121)
}
