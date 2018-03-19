test_that("profMat,xcmsRaw works", {
    fs <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    ## Create a new profile matrix:
    xr <- deepCopy(faahko_xr_1)
    expect_error(pm <- profMat(xr))
    xr_2 <- xcmsRaw(fs, profstep = 2)
    expect_equal(profMat(xr_2), profMat(xr, step = 2))
    xr_3 <- xcmsRaw(fs, profstep = 2, profmethod = "binlinbase",
                    profparam = list(baselevel = 666666))
    expect_equal(xr_3@env$profile, profMat(xr, step = 2, method = "binlinbase",
                                           baselevel = 666666))
    expect_equal(xr_3@env$profile, profMat(xr_3))
})
 
test_that("profStep<-,xcmsRaw works", {
    ## Profile matrix will be generated/replaced if the step parameter is > 0
    ## and differs from the one within the object.
    ## xr <- xcmsRaw(fs, profstep = 0)
    xr <- deepCopy(faahko_xr_1)
    xr_2 <- xr
    expect_true(length(xr_2@env$profile) == 0)
    profStep(xr_2) <- 2
    expect_true(length(xr_2@env$profile) > 0)
    expect_equal(profMat(xr_2), profMat(xr, step = 2))
    expect_equal(profMat(xr_2, step = 2), xr_2@env$profile)
})

test_that("profMethod<-,xcmsRaw works", {
    ## Profile matrix will be generated/replaced if profMethod is changed.
    xr <- deepCopy(faahko_xr_1)
    xr_2 <- xr
    expect_true(length(xr_2@env$profile) == 0)
    ## Just setting profMethod doesn't help here
    profMethod(xr_2) <- "binlin"
    expect_equal(profMethod(xr_2), "binlin")
    expect_true(length(xr_2@env$profile) == 0)
    profStep(xr_2) <- 2
    xr_3 <- xcmsRaw(fs, profstep = 2, profmethod = "binlin")
    expect_equal(profMat(xr_3), profMat(xr_2))
    ## binlinbase
    xr_3@profparam <- list(baselevel = 666666)
    profMethod(xr_3) <- "binlinbase"
    expect_equal(xr_3@env$profile, profMat(xr_3, method = "binlinbase",
                                          baselevel = 666666))
    xr_4 <- xcmsRaw(fs, profstep = 2, profmethod = "binlinbase",
                    profparam = list(baselevel = 666666))
    expect_equal(xr_4@env$profile, xr_3@env$profile)
})
