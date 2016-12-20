test.rawMat <- function() {
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xraw <- xcmsRaw(file, profstep = 0)
    xrawCopy <- deepCopy(xraw)

    checkTrue(all(xraw@env$mz == xrawCopy@env$mz))

    xrawCopy@env$mz <-     xrawCopy@env$mz+1
    checkTrue(all(xraw@env$mz != xrawCopy@env$mz))

}
