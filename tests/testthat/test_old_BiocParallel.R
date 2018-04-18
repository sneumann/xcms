test_that("BiocParallel works with xcmsSet", {
    fs <- c(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
            system.file('cdf/KO/ko19.CDF', package = "faahKO"))

    ## Default.
    bpp <- bpparam()
    ## bpprogressbar(bpp) <- TRUE
    xs2 <- xcmsSet(fs, BPPARAM=bpp, method = "centWave", noise = 10000,
                   snthresh = 50, profparam = list(profstep = 0))
    expect_equal(filepaths(xs2), fs)
    bpp <- SerialParam()
    xs3 <- xcmsSet(fs, BPPARAM=bpp, method = "centWave", noise = 10000,
                   snthresh = 50, profparam = list(profstep = 0))
    expect_equal(filepaths(xs3), fs)
    ## Drop the processHistory
    xs2@.processHistory <- list()
    xs3@.processHistory <- list()
    expect_equal(xs2, xs3)
})
