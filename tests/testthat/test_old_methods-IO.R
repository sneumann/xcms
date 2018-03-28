test_that("write.cdf,xcmsRaw works", {
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xraw <- xcmsRaw(file, profstep = 0)
    cdffile <- paste(tempdir(), "ko15.cdf", sep="/")
    write.cdf(xraw, cdffile)
    xrawCopy <- xcmsRaw(cdffile)
    expect_true(all(xraw@env$mz == xrawCopy@env$mz))
    expect_true(all(xraw@env$intensity == xrawCopy@env$intensity))
})

test_that("write.mzdata,xcmsRaw works", {
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xraw <- xcmsRaw(file, profstep = 0)
    mzdataFile <- paste(tempdir(), "ko15.mzData", sep="/")
    write.mzdata(xraw, mzdataFile)
    xrawCopy <- xcmsRaw(mzdataFile)
    expect_true(all(xraw@env$mz == xrawCopy@env$mz))
    expect_true(all(xraw@env$intensity == xrawCopy@env$intensity))
})

test_that("write.mzdata,xcmsRaw works with MS2 data", {
    file <- system.file('microtofq/MSMSpos20_6.mzML', package = "msdata")
    xraw <- xcmsRaw(file, includeMSn=TRUE, profstep = 0)
    mzdataFile <- paste(tempdir(), "MSMSpos20_6.mzData", sep="/")
    write.mzdata(xraw, mzdataFile)
    xrawCopy <- xcmsRaw(mzdataFile)
    expect_true(all(xraw@env$intensity == xrawCopy@env$intensity))
    expect_true(all(xraw@env$msnIntensity == xrawCopy@env$msnIntensity))
})

test_that("write.mzdata,xcmsRaw works with MSn data", {
    file <- system.file('threonine/threonine_i2_e35_pH_tree.mzXML', package = "msdata")
    xraw <- xcmsRaw(file, includeMSn=TRUE, profstep = 0)
    mzdataFile <- paste(tempdir(), "threonine_i2_e35_pH_tree.mzData", sep="/")
    write.mzdata(xraw, mzdataFile)
    xrawCopy <- xcmsRaw(mzdataFile)
    expect_true(all(xraw@env$intensity == xrawCopy@env$intensity))
    expect_true(all(xraw@env$msnIntensity == xrawCopy@env$msnIntensity))
})

test_that("write.mzdata,xcmsRaw writes polarity", {
    file <- system.file('microtofq/MM14.mzdata', package = "msdata")
    xraw <- xcmsRaw(file, profstep = 0)
    oldpolarity <- xraw@polarity
    mzdataFile <- paste(tempdir(), "MM14.mzdata", sep="/")
    write.mzdata(xraw, mzdataFile)
    xrawCopy <- xcmsRaw(mzdataFile)
    expect_true(all(xraw@polarity == xrawCopy@polarity))
})

test_that("write.mzQuantML,xcmsSet works", {
    xsg <- group(faahko)
    mzqFile <- paste(tempdir(), "faahKO.mzq.xml", sep="/")
    expect_warning(write.mzQuantML(xsg, mzqFile))
    v <- verify.mzQuantML(filename=mzqFile)
    expect_true(v$status == "0")    
})
