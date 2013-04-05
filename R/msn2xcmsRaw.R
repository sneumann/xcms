msn2xcmsRaw <-
function(xmsn) {
    x <- xmsn

    x@tic <- x@msnAcquisitionNum
                                        
    x@scantime <- x@msnRt          # Fake time in secs
    x@acquisitionNum <- x@msnAcquisitionNum
    x@scanindex <- x@msnScanindex

    x@env$mz <- x@env$msnMz
    x@env$intensity <- x@env$msnIntensity
    invisible(x)
}

