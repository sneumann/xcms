############################################################
## Tests related to data subset extraction from the full (raw) data.
## Most of these functions and methods are internal methods. Some
## of the tests here are "notrun" as they are used to validate the
## "transition phase" from using the original code and some unified
## new functions.
## library(xcms)
## library(RUnit)

## library(msdata)
## mzf <- c(system.file("microtofq/MM14.mzML", package = "msdata"),
##          system.file("microtofq/MM8.mzML", package = "msdata"))

## xraw <- xcmsRaw(mzf[1], profstep = 0)
xraw <- deepCopy(faahko_xr_1)

dontrun_test_rawMat <- function() {

    mz <- xraw@env$mz
    int <- xraw@env$intensity
    vps <- diff(c(xraw@scanindex, length(mz)))
    scantime <- xraw@scantime

    ## Full range
    mzr <- numeric()
    rtr <- numeric()
    scr <- numeric()
    rm <- rawMat(xraw, mzrange = mzr, rtrange = rtr, scanrange = scr)
    rm_2 <- .rawMat(mz = mz, int = int, scantime = scantime,
                    valsPerSpect = vps, mzrange = mzr,
                    rtrange = rtr, scanrange = scr)
    checkEquals(rm, rm_2)
    ## mzr
    mzr <- c(130, 150)
    rm <- rawMat(xraw, mzrange = mzr, rtrange = rtr, scanrange = scr)
    rm_2 <- .rawMat(mz = mz, int = int, scantime = scantime,
                    valsPerSpect = vps, mzrange = mzr,
                    rtrange = rtr, scanrange = scr)
    checkEquals(rm, rm_2)
    ## scr
    scr <- c(3, 9)
    rm <- rawMat(xraw, mzrange = mzr, rtrange = rtr, scanrange = scr)
    rm_2 <- .rawMat(mz = mz, int = int, scantime = scantime,
                    valsPerSpect = vps, mzrange = mzr,
                    rtrange = rtr, scanrange = scr)
    checkEquals(rm, rm_2)
    ## rtr
    rtr <- c(280, 290)
    rm <- rawMat(xraw, mzrange = mzr, rtrange = rtr, scanrange = scr)
    rm_2 <- .rawMat(mz = mz, int = int, scantime = scantime,
                    valsPerSpect = vps, mzrange = mzr,
                    rtrange = rtr, scanrange = scr)
    checkEquals(rm, rm_2)
}

############################################################
## benchmark: performance is essentially the same.
dontrun_benchmark_rawMat <- function() {
    library(microbenchmark)
    mz <- xraw@env$mz
    int <- xraw@env$intensity
    vps <- diff(c(xraw@scanindex, length(mz)))
    scantime <- xraw@scantime
    mzr <- numeric()
    rtr <- numeric()
    scr <- numeric()
    microbenchmark(rawMat(xraw, mzrange = mzr, rtrange = rtr, scanrange = scr),
                   xcms:::.rawMat(mz = mz, int = int, scantime = scantime,
                                  valsPerSpect = vps, mzrange = mzr,
                                  rtrange = rtr, scanrange = scr))
    ## mzr
    mzr <- c(130, 150)
    microbenchmark(rawMat(xraw, mzrange = mzr, rtrange = rtr, scanrange = scr),
                   xcms:::.rawMat(mz = mz, int = int, scantime = scantime,
                                  valsPerSpect = vps, mzrange = mzr,
                                  rtrange = rtr, scanrange = scr))
    ## scnr
    scr <- c(13, 54)
    microbenchmark(rawMat(xraw, mzrange = mzr, rtrange = rtr, scanrange = scr),
                   xcms:::.rawMat(mz = mz, int = int, scantime = scantime,
                                  valsPerSpect = vps, mzrange = mzr,
                                  rtrange = rtr, scanrange = scr))
    ## rtr
    rtr <- c(280, 290)
    scr <- numeric()
    microbenchmark(rawMat(xraw, mzrange = mzr, rtrange = rtr, scanrange = scr),
                   xcms:::.rawMat(mz = mz, int = int, scantime = scantime,
                                  valsPerSpect = vps, mzrange = mzr,
                                  rtrange = rtr, scanrange = scr))
    ## Compare performance with getEIC C call.
    microbenchmark(.Call("getEIC", mz, int, xraw@scanindex, as.double(mzr),
                   as.integer(scr), as.integer(length(scantime)),
                   package = "xcms"),
                   rawMat(xraw, mzrange = mzr, scanrange = scr))
    ## Although getEIC does more stuff than just subsetting, it is considerably faster.

}

