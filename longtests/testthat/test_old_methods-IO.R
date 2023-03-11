test_that("write.cdf,xcmsRaw works", {
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xraw <- xcmsRaw(file, profstep = 0)
    cdffile <- paste(tempdir(), "ko15.cdf", sep="/")
    write.cdf(xraw, cdffile)
    xrawCopy <- xcmsRaw(cdffile)
    expect_true(all(xraw@env$mz == xrawCopy@env$mz))
    expect_true(all(xraw@env$intensity == xrawCopy@env$intensity))
})

test_that("write.mzQuantML,xcmsSet works", {
    xsg <- group(faahko)
    mzqFile <- paste(tempdir(), "faahKO.mzq.xml", sep="/")
    write.mzQuantML(xsg, mzqFile)
    v <- verify.mzQuantML(filename=mzqFile)
    expect_true(v$status == "0")    
})

