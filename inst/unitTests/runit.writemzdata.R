test.write.mzdata <- function() {
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xraw <- xcmsRaw(file)

    mzdataFile <- paste(tempdir(), "ko15.mzData", sep="/")
    
    write.mzdata(xraw, mzdataFile)
    
    xrawCopy <- xcmsRaw(mzdataFile)

    checkTrue(all(xraw@env$mz == xrawCopy@env$mz))

}
