test_that("readRawData works", {
    cdf_file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")

    ## loadRaw
    lr_res <- loadRaw(xcmsSource(cdf_file))
    rr_res <- readRawData(cdf_file)
    expect_equal(lr_res, rr_res[names(lr_res)])

})
