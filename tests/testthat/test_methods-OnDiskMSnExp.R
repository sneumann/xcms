test_that("profMat,OnDiskMSnExp works", {
    ## Get it from all 3 files in one go.
    res <- profMat(faahko_od, step = 2)
    res_2 <- profMat(xcmsRaw(faahko_3_files[2], profstep = 0), step = 2)
    expect_equal(res_2, res[[2]])
    res_2 <- profMat(xcmsRaw(faahko_3_files[3], profstep = 0), step = 2)
    expect_equal(res_2, res[[3]])
    res_2 <- profMat(faahko_xod, step = 2)
    expect_equal(res, res_2)
    res <- profMat(faahko_od, step = 2, method = "binlin", fileIndex = 2)
    res_2 <- profMat(xcmsRaw(faahko_3_files[2], profstep = 0), step = 2,
                     method = "binlin")
    expect_equal(res_2, res[[1]])
})

