test_that(".createProfileMatrix works", {
    xr <- deepCopy(faahko_xr_1)
    mz <- xr@env$mz
    int <- xr@env$intensity
    numPerSc <- diff(c(xr@scanindex, length(xr@env$mz)))
    ## Testing all properties.
    ## o bin
    pm <- .createProfileMatrix(mz = mz, int = int,
                               valsPerSpect = numPerSc,
                               method = "bin", step = 2)
    ## o binlin
    pm_2 <- .createProfileMatrix(mz = mz, int = int,
                                 valsPerSpect = numPerSc,
                                 method = "binlin", step = 2)
    expect_equal(dim(pm), dim(pm_2))
    ## o binlinbase
    pm_3 <- .createProfileMatrix(mz = mz, int = int,
                                 valsPerSpect = numPerSc,
                                 method = "binlinbase", step = 2)
    expect_equal(dim(pm), dim(pm_3))
    expect_equal(sum(pm == 0), sum(pm_3 == 35))
    ##   setting parameter: baselevel
    pm_3_2 <- .createProfileMatrix(mz = mz, int = int,
                                   valsPerSpect = numPerSc,
                                   method = "binlinbase", step = 2,
                                   baselevel = 666666)
    expect_equal(sum(pm_3_2 == 666666), sum(pm == 0))
    ##   setting parameter: basespace
    pm_3_3 <- .createProfileMatrix(mz = mz, int = int,
                                   valsPerSpect = numPerSc,
                                   method = "binlinbase", step = 2,
                                   basespace = 0.5)
    expect_equal(pm_3_3, pm_3)
    pm_3_3 <- .createProfileMatrix(mz = mz, int = int,
                                   valsPerSpect = numPerSc,
                                   method = "binlinbase", step = 2,
                                   basespace = 300)
    expect_true(!all(pm_3_3 == pm_3))
    ## o intlin
    pm_4 <- .createProfileMatrix(mz = mz, int = int,
                                 valsPerSpect = numPerSc,
                                 method = "intlin", step = 2)
    expect_equal(dim(pm), dim(pm_4))
})

test_that("plotMsData works", {
    msd <- extractMsData(faahko_od, mz = c(334.9, 335.1), rt = c(2700, 2900))
    plotMsData(msd[[1]])
})
