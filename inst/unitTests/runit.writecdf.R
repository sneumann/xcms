test.write.cdf <- function() {
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xraw <- xcmsRaw(file, profstep = 0)

    cdffile <- paste(tempdir(), "ko15.cdf", sep="/")

    write.cdf(xraw, cdffile)

    xrawCopy <- xcmsRaw(cdffile)

    checkTrue(all(xraw@env$mz == xrawCopy@env$mz))
    checkTrue(all(xraw@env$intensity == xrawCopy@env$intensity))

}

if (FALSE) {
    library(xcms)
    library(RUnit)
}
