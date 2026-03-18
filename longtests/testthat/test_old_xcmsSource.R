test_that("xcmsSource works", {
    mz_file <- faahko_3_files[1L]
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
})
