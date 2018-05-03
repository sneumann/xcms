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

test_that("rla, rowRla work", {
    x <- c(3, 4, 5, 1, 2, 3, 7, 8, 9)
    grp <- c(1, 1, 1, 2, 2, 2, 3, 3, 3)
    res <- rla(x, grp)
    expect_equal(unname(res[1]), log2(3) - median(log2(c(3, 4, 5))))
    expect_equal(unname(res[2]), log2(4) - median(log2(c(3, 4, 5))))
    expect_equal(unname(res[3]), log2(5) - median(log2(c(3, 4, 5))))
    expect_equal(unname(res[4]), log2(1) - median(log2(c(1, 2, 3))))
    expect_equal(unname(res[5]), log2(2) - median(log2(c(1, 2, 3))))
    expect_equal(unname(res[6]), log2(3) - median(log2(c(1, 2, 3))))

    idx <- c(1, 5, 3, 8, 9, 2, 6, 4, 7)
    res_2 <- rla(x[idx], grp[idx])
    expect_identical(res_2, res[idx])

    mat <- rbind(x, x, x, x)

    res_mat <- rowRla(mat, grp)
    expect_equal(res_mat[2, ], res)
})
