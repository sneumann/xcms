test_that(".createProfileMatrix works", {
    skip_on_os(os = "windows", arch = "i386")

    xr <- filterFile(faahko_od, 1)
    mz <- mz(xr)
    int <- unlist(intensity(xr), use.names = FALSE)
    numPerSc <- lengths(mz)
    mz <- unlist(mz, use.names = FALSE)
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
    skip_on_os(os = "windows", arch = "i386")

    msd <- extractMsData(faahko_od, mz = c(334.9, 335.1), rt = c(2700, 2900))
    plotMsData(msd[[1]])
})

test_that(".featureIDs works", {
    skip_on_os(os = "windows", arch = "i386")

    res <- .featureIDs(200)
    expect_equal(length(res), 200)
    expect_true(length(unique(res)) == 200)
    res <- .featureIDs(221495, prefix = "CP")
    expect_true(length(unique(res)) == 221495)
})

test_that("rla, rowRla work", {
    skip_on_os(os = "windows", arch = "i386")

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
    skip_on_os(os = "windows", arch = "i386")

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

test_that(".insertColumn works", {
    skip_on_os(os = "windows", arch = "i386")

    mat <- matrix(1:100, ncol = 5)
    expect_equal(.insertColumn(mat), mat)

    expect_error(.insertColumn(mat, 3))
    expect_error(.insertColumn(mat, 3, 3:4))

    res <- .insertColumn(mat, 3, 5)
    expect_true(all(res[, 3] == 5))
    expect_equal(res[, -3], mat)

    res <- .insertColumn(mat, c(2, 4), 6)
    expect_true(ncol(res) == ncol(mat) + 2)
    expect_equal(mat, res[, -c(2, 4)])
    expect_true(all(res[, 2] == 6))
    expect_true(all(res[, 4] == 6))

    res <- .insertColumn(mat, c(2, 4), list(101:120))
    expect_true(ncol(res) == ncol(mat) + 2)
    expect_equal(res[, 2], 101:120)
    expect_equal(res[, 4], 101:120)
})

test_that(".ppm_range works", {
    skip_on_os(os = "windows", arch = "i386")

    res <- .ppm_range(100)
    expect_equal(res[1], 100)
    expect_equal(res[2], 100)
    res <- .ppm_range(100, 100)
    expect_equal(res[1], 100 - 5000 / 1e6)
    expect_equal(res[2], 100 + 5000 / 1e6)
})

test_that(".update_feature_definitions works", {
    skip_on_os(os = "windows", arch = "i386")

    cps <- matrix(nrow = 22, ncol = 3)
    rownames(cps) <- 1:22
    fts <- DataFrame(a = letters[1:6])
    pidx <- list(
        c(1, 2, 3, 6, 9, 12),
        c(5, 10, 22),
        c(4, 9, 13, 14, 15, 16, 17),
        c(11, 15, 18, 19),
        c(17, 20, 21, 22),
        c(5, 13, 17)
    )
    fts$peakidx <- pidx
    cps_sub <- cps[c(4, 6, 17, 19), ]
    res <- .update_feature_definitions(fts, rownames(cps), rownames(cps_sub))
    expect_equal(res$a, c("a", "c", "d", "e", "f"))
    expect_equal(res$peakidx[[1]], c(2))
    expect_equal(res$peakidx[[2]], c(1, 3))
    expect_equal(res$peakidx[[3]], c(4))
    expect_equal(res$peakidx[[4]], c(3))
    cps_sub <- cps[1:10, ]
    res <- .update_feature_definitions(fts, rownames(cps), rownames(cps_sub))
    expect_equal(res$a, c("a", "b", "c", "f"))
    expect_equal(res$peakidx[[1]], c(1, 2, 3, 6, 9))
    expect_equal(res$peakidx[[2]], c(5, 10))
    expect_equal(res$peakidx[[3]], c(4, 9))
    expect_equal(res$peakidx[[4]], 5)

    ## Real data set:
    orig_names <- rownames(chromPeaks(xod_xgrg))
    sub_names <- sample(orig_names, (length(orig_names) / 2))
    fts <- featureDefinitions(xod_xgrg)
    res <- xcms:::.update_feature_definitions(fts, orig_names, sub_names)
    expect_s4_class(res, "DataFrame")
    expect_true(all(lengths(res$peakidx) > 0))
    tmp <- lapply(res$peakidx, function(z) sub_names[z])
    tmp <- unlist(tmp, use.names = FALSE)

    onames <- intersect(orig_names[unlist(fts$peakidx, use.names = FALSE)],
                        sub_names)
    expect_true(all(onames %in% tmp))
    expect_true(all(tmp %in% sub_names))
})

## test_that(".chrom_peak_id works", {
    ## skip_on_os(os = "windows", arch = "i386")

##     res <- .chrom_peak_id(matrix(nrow = 0, ncol = 5))
##     expect_equal(res, character())
##     cpks <- rbind(c(3, 2, 4, 12, 13),
##                   c(4, 2, 4, 123, 43),
##                   c(3, 2, 4, 12, 13),
##                   c(5, 4, 6, 123, 45))
##     colnames(cpks) <- c("rt", "rtmin", "rtmax", "into", "maxo")
##     expect_error(.chrom_peak_id(cpks))
##     res <- .chrom_peak_id(cpks[-3, ])
##     expect_equal(res, c("3-2-4-12-13", "4-2-4-123-43", "5-4-6-123-45"))
##     cpks <- chromPeaks(xod_x)
##     res <- .chrom_peak_id(cpks)
## })

test_that(".rbind_fill works", {
    skip_on_os(os = "windows", arch = "i386")

    ## matrix
    a <- matrix(1:9, nrow = 3, ncol = 3)
    colnames(a) <- c("a", "b", "c")
    b <- matrix(1:12, nrow = 3, ncol = 4)
    colnames(b) <- c("b", "a", "d", "e")
    res <- .rbind_fill(a, b)
    expect_equal(colnames(res), c("a", "b", "c", "d", "e"))
    expect_equal(class(res), class(a))
    expect_equal(res[, "a"], c(a[, "a"], b[, "a"]))
    expect_equal(res[, "b"], c(a[, "b"], b[, "b"]))
    expect_equal(res[, "d"], c(NA, NA, NA, b[, "d"]))

    res <- .rbind_fill(a, b[, c("b", "a")])
    expect_equal(colnames(res), c("a", "b", "c"))
    expect_equal(res[, "a"], c(a[, "a"], b[, "a"]))

    ## DataFrame
    a <- DataFrame(a = 1:4, b = FALSE, c = letters[1:4])
    b <- DataFrame(d = 1:4, b = TRUE)
    res <- .rbind_fill(a, b)
    expect_equal(colnames(res), c("a", "b", "c", "d"))
    expect_equal(res$a, c(1:4, NA, NA, NA, NA))
    expect_equal(res$b, rep(c(FALSE, TRUE), each = 4))
})

test_that(".reduce works", {
    skip_on_os(os = "windows", arch = "i386")

    a <- c(1.23, 1.431, 2.43, 5.44, 6)
    b <- c(1.33, 2.43, 5, 6, 7)
    res <- .reduce(a, b)
    expect_true(nrow(res) == 3)
    expect_equal(res[, 1], c(1.23, 1.431, 5.44))
    expect_equal(res[, 2], c(1.33, 5, 7))

    idx <- sample(1:length(a))
    res_2 <- .reduce(a[idx], b[idx])
    expect_identical(res, res_2)

    res <- .reduce(a[1], b[1])
    expect_equal(res, cbind(start = a[1], end = b[1]))

    res <- .reduce(numeric(), numeric())
    expect_equal(nrow(res), 0)

    res <- .reduce(a - 0.1, b + 0.1)
    expect_equal(res[, 1], c(1.13, 5.34))
    expect_equal(res[, 2], c(5.1, 7.1))

    a <- c(4, 4)
    b <- c(5, 5)
    res <- .reduce(a, b)
    expect_true(nrow(res) == 1)
    expect_equal(res[1, 1], c(start = 4))
    expect_equal(res[1, 2], c(end = 5))

    a <- c(3, 4, 8)
    b <- c(7, 5, 10)
    res <- .reduce(a, b)
    expect_equal(res[, 1], c(3, 8))
    expect_equal(res[, 2], c(7, 10))

    a <- c(3, 4, 6)
    b <- c(7, 5, 10)
    res <- .reduce(a, b)
    expect_equal(unname(res[, 1]), 3)
    expect_equal(unname(res[, 2]), 10)
})

test_that("groupOverlaps works", {
    skip_on_os(os = "windows", arch = "i386")

    x <- c(12.2, 13, 5)
    y <- c(16, 15, 6)
    res <- groupOverlaps(x, y)
    expect_true(is.list(res))
    expect_equal(length(res), 2)
    expect_equal(res, list(3, 1:2))

    expect_error(groupOverlaps(x, 1:2), "lengths differ")
})

test_that(".require_spectra works", {
    skip_on_os(os = "windows", arch = "i386")

    if (requireNamespace("Spectra", quietly = TRUE))
        expect_true(.require_spectra())
    else expect_error("installed.")
})

test_that(".i2index works", {
    skip_on_os(os = "windows", arch = "i386")

    ids <- c("a", "b", "c", "d")
    res <- .i2index(c("c", "d"), ids)
    expect_equal(res, c(3L, 4L))
    expect_error(.i2index(12, ids), "out of bounds")
})
