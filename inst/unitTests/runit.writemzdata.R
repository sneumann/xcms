test.write.mzdata <- function() {
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xraw <- xcmsRaw(file, profstep = 0)

    mzdataFile <- paste(tempdir(), "ko15.mzData", sep="/")

    write.mzdata(xraw, mzdataFile)

    xrawCopy <- xcmsRaw(mzdataFile)

    checkTrue(all(xraw@env$mz == xrawCopy@env$mz))
    checkTrue(all(xraw@env$intensity == xrawCopy@env$intensity))

}

test.writeMS2.mzdata <- function() {
    file <- system.file('microtofq/MSMSpos20_6.mzML', package = "msdata")
    xraw <- xcmsRaw(file, includeMSn=TRUE, profstep = 0)

    mzdataFile <- paste(tempdir(), "MSMSpos20_6.mzData", sep="/")

    write.mzdata(xraw, mzdataFile)

    xrawCopy <- xcmsRaw(mzdataFile)

    checkTrue(all(xraw@env$intensity == xrawCopy@env$intensity))
    checkTrue(all(xraw@env$msnIntensity == xrawCopy@env$msnIntensity))
}

test.writeMSn.mzdata <- function() {
    file <- system.file('threonine/threonine_i2_e35_pH_tree.mzXML', package = "msdata")
    xraw <- xcmsRaw(file, includeMSn=TRUE, profstep = 0)

    mzdataFile <- paste(tempdir(), "threonine_i2_e35_pH_tree.mzData", sep="/")

    write.mzdata(xraw, mzdataFile)

    xrawCopy <- xcmsRaw(mzdataFile)

    checkTrue(all(xraw@env$intensity == xrawCopy@env$intensity))
    checkTrue(all(xraw@env$msnIntensity == xrawCopy@env$msnIntensity))
}

test.writePolarity.mzdata <- function() {
    file <- system.file('microtofq/MM14.mzdata', package = "msdata")
    xraw <- xcmsRaw(file, profstep = 0)

    oldpolarity <- xraw@polarity

    mzdataFile <- paste(tempdir(), "MM14.mzdata", sep="/")

    write.mzdata(xraw, mzdataFile)

    xrawCopy <- xcmsRaw(mzdataFile)

    checkTrue(all(xraw@polarity == xrawCopy@polarity))
}
