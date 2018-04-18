test_that("xcmsSource works", {
    mz_file <- system.file("microtofq/MM8.mzML", package = "msdata")
    src <- xcms:::xcmsSource(mz_file)
    expect_true(is(src, "xcmsFileSource"))
    tmp <- loadRaw(src)
    expect_equal(names(tmp), c("rt", "acquisitionNum", "tic", "scanindex",
                              "mz", "intensity", "polarity"))

    cdf_file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    src <- xcms:::xcmsSource(cdf_file)
    expect_true(is(src, "xcmsFileSource"))
    tmp <- loadRaw(src)
    expect_equal(names(tmp), c("rt", "acquisitionNum", "tic", "scanindex",
                              "mz", "intensity", "polarity"))

    ## MSn:
    mzdatapath <- system.file("iontrap", package = "msdata")
    mzdatafiles <- list.files(mzdatapath, pattern="extracted.mzData",
                              recursive = TRUE, full.names = TRUE)
    src <- xcms:::xcmsSource(mzdatafiles[1])
    tmp <- loadRaw(src, includeMSn = TRUE)

    ## OLD code:
    rid <- mzR:::rampOpen(mzdatafiles[1])
    rawdata <- mzR:::rampRawData(rid)
    rawdata$MSn <- mzR:::rampRawDataMSn(rid)
    mzR:::rampClose(rid)
    rm(rid)
    ## Ramp does not read polarity!
    tmp$polarity <- rawdata$polarity
    expect_equal(rawdata, tmp)

    ## Next example:
    msnfile <- system.file("microtofq/MSMSpos20_6.mzML", package = "msdata")
    src <- xcms:::xcmsSource(msnfile)
    tmp <- loadRaw(src, includeMSn = TRUE)
    ## expect_true(all(tmp$polarity == 1))
    ## OLD code:
    rid <- mzR:::rampOpen(msnfile)
    rawdata <- mzR:::rampRawData(rid)
    rawdata$MSn <- mzR:::rampRawDataMSn(rid)
    mzR:::rampClose(rid)
    rm(rid)
    rawdata$polarity <- tmp$polarity
    expect_equal(rawdata, tmp)
})
