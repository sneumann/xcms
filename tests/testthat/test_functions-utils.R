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

test_that(".rect_overlap works", {
    xl <- c(1, 3, 1.5, 4, 4, 5.5, 7, 6)
    xr <- c(2, 4, 3.5, 5, 5, 6.5, 8, 7.5)
    yb <- c(1, 2, 3.5, 4.5, 7, 8, 9.5, 10.5)
    yt <- c(3, 4, 5, 6, 9, 10, 11, 12)
    names(xl) <- c("a", "b", "c", "d", "e", "f", "g", "h")

    ## plot(3, 3, pch = NA, xlim = range(c(xl, xr)), ylim = range(c(yb, yt)),
    ##      xlab = "x", ylab = "y")
    ## rect(xleft = xl, xright = xr, ybottom = yb, ytop = yt)
    ## text(x = rowMeans(cbind(xl, xr)), y = rowMeans(cbind(yb, yt)),
    ##      labels = names(xl))
    res <- .rect_overlap(xl, xr, yb, yt)
    ## Expecting that b and c as well as g and h are overlapping
    expect_equal(res, list(c(2, 3), c(7, 8)))

    ## Expand them
    xl_2 <- xl - 1
    xr_2 <- xr + 1
    yb_2 <- yb
    yt_2 <- yt
    plot(3, 3, pch = NA, xlim = range(c(xl_2, xr_2)),
         ylim = range(c(yb_2, yt_2)), xlab = "x", ylab = "y")
    rect(xleft = xl_2, xright = xr_2, ybottom = yb_2, ytop = yt_2)
    text(x = rowMeans(cbind(xl_2, xr_2)), y = rowMeans(cbind(yb_2, yt_2)),
         labels = names(xl_2))
    res <- .rect_overlap(xl_2, xr_2, yb_2, yt_2)
    expect_equal(res, list(c(1:4), 5:8))

    idx <- sample(1:length(xl_2), length(xl_2))
    xl_2 <- xl[idx]
    xr_2 <- xr[idx]
    yb_2 <- yb[idx]
    yt_2 <- yt[idx]
    plot(3, 3, pch = NA, xlim = range(c(xl_2, xr_2)),
         ylim = range(c(yb_2, yt_2)), xlab = "x", ylab = "y")
    rect(xleft = xl_2, xright = xr_2, ybottom = yb_2, ytop = yt_2)
    text(x = rowMeans(cbind(xl_2, xr_2)), y = rowMeans(cbind(yb_2, yt_2)),
         labels = names(xl_2))
    res <- .rect_overlap(xl_2, xr_2, yb_2, yt_2)
    expect_equal(length(res), 2)
    expect_true(all(lengths(res) %in% c(2, 2)))
    expect_true((all(names(xl_2)[res[[1]]] %in% c("c", "b")) |
                 (all(names(xl_2)[res[[2]]] %in% c("c", "b")))))
    expect_true((all(names(xl_2)[res[[1]]] %in% c("g", "h")) |
                 (all(names(xl_2)[res[[2]]] %in% c("g", "h")))))


    ## Expand them
    xl_2 <- xl - 0.4
    xr_2 <- xr + 0.4
    yb_2 <- yb
    yt_2 <- yt
    plot(3, 3, pch = NA, xlim = range(c(xl_2, xr_2)),
         ylim = range(c(yb_2, yt_2)), xlab = "x", ylab = "y")
    rect(xleft = xl_2, xright = xr_2, ybottom = yb_2, ytop = yt_2)
    text(x = rowMeans(cbind(xl_2, xr_2)), y = rowMeans(cbind(yb_2, yt_2)),
         labels = names(xl_2))
    res <- .rect_overlap(xl_2, xr_2, yb_2, yt_2)
    expect_equal(res, list(c(2:4), 5:8))

    ## Expand them
    xl_2 <- xl - 0.1
    xr_2 <- xr + 0.1
    yb_2 <- yb - 0.3
    yt_2 <- yt + 0.3
    plot(3, 3, pch = NA, xlim = range(c(xl_2, xr_2)),
         ylim = range(c(yb_2, yt_2)), xlab = "x", ylab = "y")
    rect(xleft = xl_2, xright = xr_2, ybottom = yb_2, ytop = yt_2)
    text(x = rowMeans(cbind(xl_2, xr_2)), y = rowMeans(cbind(yb_2, yt_2)),
         labels = names(xl_2))
    res <- .rect_overlap(xl_2, xr_2, yb_2, yt_2)
    expect_equal(res, list(c(1:4), 6:8))
    
    idx <- sample(1:length(xl_2), length(xl_2))
    xl_2 <- xl_2[idx]
    xr_2 <- xr_2[idx]
    yb_2 <- yb_2[idx]
    yt_2 <- yt_2[idx]
    plot(3, 3, pch = NA, xlim = range(c(xl_2, xr_2)),
         ylim = range(c(yb_2, yt_2)), xlab = "x", ylab = "y")
    rect(xleft = xl_2, xright = xr_2, ybottom = yb_2, ytop = yt_2)
    text(x = rowMeans(cbind(xl_2, xr_2)), y = rowMeans(cbind(yb_2, yt_2)),
         labels = names(xl_2))
    res <- .rect_overlap(xl_2, xr_2, yb_2, yt_2)
    expect_equal(length(res), 2)
    expect_true(all(lengths(res) %in% c(3, 4)))
    expect_true((all(names(xl_2)[res[[1]]] %in% c("a", "b", "c", "d")) |
                 (all(names(xl_2)[res[[2]]] %in% c("a", "b", "c", "d")))))
    expect_true((all(names(xl_2)[res[[1]]] %in% c("f", "g", "h")) |
                 (all(names(xl_2)[res[[2]]] %in% c("f", "g", "h")))))

})
