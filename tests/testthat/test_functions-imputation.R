test_that("imputeRowMin works", {
    skip_on_os(os = "windows", arch = "i386")

    mat <- cbind(c(4, 2, 4, NA, NA), c(NA, NA, NA, NA, NA), c(4, NA, 6, 3, 9),
                 c(6, 3, NA, 6, NA))
    mat_imp <- imputeRowMin(mat)
    expect_equal(mat_imp[, 2], c(2, 1, 2, 1.5, 4.5))

    mat_imp <- imputeRowMin(mat[, -3])
    expect_equal(mat_imp[, 2], c(2, 1, 2, 3, NA))
})

test_that("imputeRowMinRand works", {
    skip_on_os(os = "windows", arch = "i386")

    set.seed(123)
    mat <- cbind(c(4, 2, 4, NA, NA), c(NA, NA, NA, NA, NA), c(4, NA, 6, 3, 9),
                 c(6, 3, NA, 6, NA))
    mat_imp <- imputeRowMinRand(mat)
    rmin <- apply(mat, 1, min, na.rm = TRUE)
    rmin_imp <- apply(mat_imp, 1, min, na.rm = TRUE)
    expect_true(all(rmin_imp < rmin))

    mat_2 <- imputeRowMinRand(mat, method = "from_to")
    expect_true(all(!is.na(mat_2)))
    rmin_imp2 <- apply(mat_2, 1, min, na.rm = TRUE)
    expect_true(all(rmin_imp2 < rmin))
    expect_true(all(rmin_imp2 < rmin / 2))
    expect_true(all(rmin_imp2 > rmin / 1000))
    expect_true(all(rmin_imp != rmin_imp2))

    mat_imp <- imputeRowMinRand(mat[, -3])
    expect_true(all(is.na(mat_imp[5, ])))
})
