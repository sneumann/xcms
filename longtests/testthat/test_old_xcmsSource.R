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
    mzmlpath <- system.file("iontrap", package = "msdata")
    mzmlfiles <- list.files(mzmlpath, pattern="extracted.mzML",
                              recursive = TRUE, full.names = TRUE)
    src <- xcms:::xcmsSource(mzmlfiles[1])
    tmp <- loadRaw(src, includeMSn = TRUE)

})
