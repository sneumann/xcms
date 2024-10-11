library(rhdf5)
xmse_h5 <- .xcms_experiment_to_hdf5(loadXcmsData("faahko_sub2"), tempfile())

test_that(".xcms_experiment_to_hdf5 works", {
    expect_error(.xcms_experiment_to_hdf5(4), "Parameter 'h5_file'")

    h5f <- tempfile()
    ref <- XcmsExperiment()
    res <- .xcms_experiment_to_hdf5(ref, h5f)
    expect_true(validObject(res))
    expect_true(is(res, "XcmsExperimentHdf5"))
    expect_true(nrow(res@chromPeaks) == 0)
    expect_true(nrow(res@chromPeakData) == 0)
    expect_true(length(res@sample_id) == 0)
    expect_equal(res@chromPeakData, ref@chromPeakData)
    expect_equal(res@chromPeaks, ref@chromPeaks)
    expect_equal(res@hdf5_mod_count, 0L)
    h5 <- rhdf5::H5Fopen(h5f)
    hdr <- rhdf5::h5read(h5, "header")
    expect_equal(as.vector(hdr$modcount), 0)
    expect_equal(as.vector(hdr$package), "package:xcms")
    expect_false(hasChromPeaks(res))

    rhdf5::H5Fclose(h5)
    file.remove(h5f)

    ref <- loadXcmsData()
    ref <- dropFeatureDefinitions(ref)
    res <- .xcms_experiment_to_hdf5(ref, h5f)
    expect_true(validObject(res))
    expect_true(is(res, "XcmsExperimentHdf5"))
    expect_true(nrow(res@chromPeaks) == 0)
    expect_true(nrow(res@chromPeakData) == 0)
    expect_true(length(res@sample_id) == length(ref))
    expect_true(res@hdf5_mod_count > 0)
    h5 <- rhdf5::H5Fopen(h5f)
    ## General checks
    expect_equal(.h5_ms_levels(h5, "S1"), 1L)
    expect_true(hasChromPeaks(res))
    expect_true(hasChromPeaks(res, 1L))
    expect_false(hasChromPeaks(res, 2L))
    expect_false(hasChromPeaks(res, c(1L, 2L)))
    ##
    hdr <- rhdf5::h5read(h5, "header")
    expect_equal(as.vector(hdr$modcount), res@hdf5_mod_count)
    expect_equal(as.vector(hdr$package), "package:xcms")
    pks <- rhdf5::h5read(h5, "/S8/ms_1/chrom_peaks")
    pks_rn <- as.vector(rhdf5::h5read(h5, "/S8/ms_1/chrom_peaks_rownames"))
    expect_true(!anyDuplicated(pks_rn))
    pks_cn <- as.vector(rhdf5::h5read(h5, "/S8/ms_1/chrom_peaks_colnames"))
    colnames(pks) <- pks_cn
    pks_ref <- chromPeaks(ref)[chromPeaks(ref)[, "sample"] == 8, ]
    rownames(pks_ref) <- NULL
    expect_equal(colnames(pks_ref), c(pks_cn, "sample"))
    expect_equal(pks_ref[, colnames(pks_ref) != "sample"], pks)
    pkd <- as.data.frame(rhdf5::h5read(h5, "/S8/ms_1/chrom_peak_data"))
    pkd_ref <- chromPeakData(
        ref, return.type = "data.frame")[chromPeaks(ref)[, "sample"] == 8, ]
    expect_equal(colnames(pkd_ref), c("ms_level", colnames(pkd)))
    expect_equal(pkd_ref$is_filled, pkd$is_filled)
    expect_equal(pkd_ref$is_merged, pkd$is_merged)
    rhdf5::H5Fclose(h5)
    file.remove(h5f)
    expect_error(validObject(res), "Data storage file")
})

test_that(".h5_dataset_names works", {
    h5 <- rhdf5::H5Fopen(xmse_h5@hdf5_file)
    res <- .h5_dataset_names("/", h5)
    expect_true(length(res) > 1)
    expect_true("header" %in% res)
    res <- .h5_dataset_names("/header", h5)
    expect_equal(res, c("modcount", "package"))
    rhdf5::H5Fclose(h5)
})

test_that(".h5_ms_levels works", {
    h5 <- rhdf5::H5Fopen(xmse_h5@hdf5_file)
    res <- .h5_ms_levels(h5, "S1")
    expect_equal(res, 1L)
    rhdf5::H5Fclose(h5)
})

test_that(".h5_chrom_peak_ms_levels works", {
    res <- .h5_chrom_peak_ms_levels(xmse_h5@hdf5_file, "S1")
    expect_equal(res, 1L)
})

test_that(".h5_subset_xcms_experiment works", {
    a <- new("XcmsExperimentHdf5")
    res <- .h5_subset_xcms_experiment(a)
    expect_equal(a, res)
    a <- xmse_h5
    res <- .h5_subset_xcms_experiment(a, c(1, 3))
    expect_equal(length(res), 2)
    expect_equal(res@sample_id, a@sample_id[c(1L, 3L)])
    expect_equal(sampleData(res), sampleData(a)[c(1, 3), ])
    expect_true(hasChromPeaks(res))
    res <- .h5_subset_xcms_experiment(a, c(1, 3), keepChromPeaks = FALSE)
    expect_equal(length(res), 2)
    expect_equal(res@sample_id, a@sample_id[c(1L, 3L)])
    expect_equal(sampleData(res), sampleData(a)[c(1, 3), ])
    expect_false(hasChromPeaks(res))
    expect_true(length(res@processHistory) == 0)
    res <- .h5_subset_xcms_experiment(a, c(1, 3), keepChromPeaks = FALSE,
                                      ignoreHistory = TRUE)
    expect_equal(length(res), 2)
    expect_equal(res@sample_id, a@sample_id[c(1L, 3L)])
    expect_equal(sampleData(res), sampleData(a)[c(1, 3), ])
    expect_false(hasChromPeaks(res))
    expect_equal(res@processHistory, a@processHistory)
})

test_that(".h5_xmse_merge_neighboring_peaks works", {
    h5f <- tempfile()
    ref <- loadXcmsData("faahko_sub2")
    x <- .xcms_experiment_to_hdf5(ref, h5f)
    ref <- .h5_read_data(x@hdf5_file, index = x@sample_id,
                         ms_level = rep(1L, length(x)),
                         read_colnames = TRUE, read_rownames = TRUE)
    .h5_xmse_merge_neighboring_peaks(x)
    mod_count <- as.vector(rhdf5::h5read(h5f, "/header/modcount"))
    expect_true(mod_count > x@hdf5_mod_count)
    ## Check that content was changed.
    res <- .h5_read_data(x@hdf5_file, index = x@sample_id,
                         ms_level = rep(1L, length(x)),
                         read_colnames = TRUE, read_rownames = TRUE)
    expect_true(nrow(ref[[1L]]) > nrow(res[[1L]]))
    expect_true(nrow(ref[[2L]]) > nrow(res[[2L]]))
    expect_true(nrow(ref[[3]]) > nrow(res[[3L]]))
    expect_true(!anyDuplicated(rownames(res[[1L]])))
    expect_true(!anyDuplicated(rownames(res[[2L]])))
    expect_true(!anyDuplicated(rownames(res[[3L]])))
    same <- intersect(rownames(ref[[1L]]), rownames(res[[1L]]))
    expect_equal(ref[[1L]][same, ], res[[1L]][same, ])
    same <- intersect(rownames(ref[[2L]]), rownames(res[[2L]]))
    expect_equal(ref[[2L]][same, ], res[[2L]][same, ])
    same <- intersect(rownames(ref[[3L]]), rownames(res[[3L]]))
    expect_equal(ref[[3L]][same, ], res[[3L]][same, ])
    file.remove(h5f)

    ## Compare with reference results.
    h5f <- tempfile()
    ref <- loadXcmsData("faahko_sub2")
    res <- .xcms_experiment_to_hdf5(ref, h5f)

    .h5_xmse_merge_neighboring_peaks(res)
    res <- .h5_read_data(res@hdf5_file, index = res@sample_id,
                         ms_level = rep(1L, length(res)),
                         read_colnames = TRUE, read_rownames = TRUE)
    ref <- .xmse_merge_neighboring_peaks(ref)
    ref <- split.data.frame(ref[[1L]], ref[[1L]][, "sample"])
    expect_equal(unname(ref[[1L]][, colnames(ref[[1L]]) != "sample"]),
                 unname(res[[1L]]))
    expect_equal(unname(ref[[2L]][, colnames(ref[[2L]]) != "sample"]),
                 unname(res[[2L]]))
    expect_equal(unname(ref[[3L]][, colnames(ref[[3L]]) != "sample"]),
                 unname(res[[3L]]))
})

test_that(".h5_read_chrom_peaks works", {
    h5f <- tempfile()
    .h5_initialize_file(h5f)

    a <- cbind(a = c(1.2, 1.4), b = c(3.5, 3.6), c = c(5.3, 5.1))
    rownames(a) <- c("CP1", "CP2")
    b <- cbind(a = c(12.4, 13.1, 3.2), b = c(1.3, 1.3, 1.4), c(4.2, 5.1, 4.6))
    rownames(b) <- c("CP3", "CP4", "CP5")
    l <- list(a, b)
    names(l) <- c(1, 2)
    expect_equal(.h5_write_data(h5f, l, name = "chrom_peaks",
                                ms_level = c(1L, 1L)), 1L)

    h5 <- H5Fopen(h5f)
    res <- .h5_read_chrom_peaks("/1/ms_1/chrom_peaks", h5)
    expect_equal(res, unname(a))
    res <- .h5_read_chrom_peaks("/1/ms_1/chrom_peaks", h5,
                                read_colnames = TRUE)
    expect_equal(unname(res), unname(a))
    expect_equal(colnames(res), colnames(a))
    expect_equal(rownames(res), NULL)
    res <- .h5_read_chrom_peaks("/1/ms_1/chrom_peaks", h5,
                                read_colnames = TRUE, read_rownames = TRUE)
    expect_equal(res, a)
    res <- .h5_read_chrom_peaks("/1/ms_1/chrom_peaks", h5, index = 3,
                                read_colnames = TRUE, read_rownames = TRUE)
    expect_equal(res, a[, 3, drop = FALSE])
    res <- .h5_read_chrom_peaks("/1/ms_1/chrom_peaks", h5, index = c(1, 3),
                                read_colnames = FALSE, read_rownames = FALSE)
    expect_equal(res, unname(a[, c(1, 3)]))
    H5Fclose(h5)
    file.remove(h5f)
})

test_that(".h5_read_chrom_peak_data works", {
    h5f <- tempfile()
    .h5_initialize_file(h5f)

    a <- data.frame(is_filled = c(TRUE, FALSE), other_col = "c")
    b <- data.frame(is_filled = c(FALSE, FALSE, TRUE), other_col = "d")
    l <- list(a, b)
    names(l) <- 1:2
    .h5_write_data(h5f, l, name = "chrom_peak_data", ms_level = c(2L, 2L))
    a <- cbind(a = c(1.2, 1.4), b = c(3.5, 3.6), c = c(5.3, 5.1))
    rownames(a) <- c("CP1", "CP2")
    b <- cbind(a = c(12.4, 13.1, 3.2), b = c(1.3, 1.3, 1.4), c(4.2, 5.1, 4.6))
    rownames(b) <- c("CP3", "CP4", "CP5")
    l <- list(a, b)
    names(l) <- c(1, 2)
    .h5_write_data(h5f, l, name = "chrom_peaks", ms_level = c(2L, 2L))

    h5 <- H5Fopen(h5f)
    res <- .h5_read_chrom_peak_data("/1/ms_2/chrom_peak_data", h5,
                                    read_rownames = TRUE)
    expect_equal(rownames(res), rownames(a))
    expect_equal(colnames(res), c("is_filled", "other_col"))
    res <- .h5_read_chrom_peak_data("/1/ms_2/chrom_peak_data", h5,
                                    read_rownames = FALSE)
    expect_equal(colnames(res), c("is_filled", "other_col"))
    expect_equal(rownames(res), c("1", "2"))

    ## Read single column
    res <- .h5_dataset_names("/1/ms_2/chrom_peak_data", h5)
    expect_equal(res, c("is_filled", "other_col"))
    res <- .h5_read_chrom_peak_data("/2/ms_2/chrom_peak_data/other_col", h5,
                                    read_rownames = FALSE)
    expect_equal(res[, 1L], c("d", "d", "d"))
    H5Fclose(h5)
    file.remove(h5f)
})

test_that(".h5_read_data works", {
    h5f <- tempfile()
    .h5_initialize_file(h5f)

    a <- data.frame(is_filled = c(TRUE, FALSE), other_col = "c")
    b <- data.frame(is_filled = c(FALSE, FALSE, TRUE), other_col = "d")
    l <- list(a, b)
    names(l) <- 1:2
    .h5_write_data(h5f, l, name = "chrom_peak_data", ms_level = c(2L, 2L))
    a2 <- cbind(a = c(1.2, 1.4), b = c(3.5, 3.6), c = c(5.3, 5.1))
    rownames(a2) <- c("CP1", "CP2")
    b2 <- cbind(a = c(12.4, 13.1, 3.2), b = c(1.3, 1.3, 1.4),
                c = c(4.2, 5.1, 4.6))
    rownames(b2) <- c("CP3", "CP4", "CP5")
    l <- list(a2, b2)
    names(l) <- c(1, 2)
    .h5_write_data(h5f, l, name = "chrom_peaks", ms_level = c(2L, 2L))

    ## chrom peaks
    res <- .h5_read_data(h5f)
    expect_equal(res, list())
    res <- .h5_read_data(h5f, index = 2, name = "chrom_peaks", ms_level = 2L)
    expect_equal(length(res), 1L)
    expect_equal(res[[1L]], unname(b2))
    res <- .h5_read_data(h5f, index = 1, name = "chrom_peaks", ms_level = 2L,
                         read_colnames = TRUE)
    expect_equal(unname(res[[1L]]), unname(a2))
    expect_equal(colnames(res[[1L]]), colnames(a2))
    expect_true(is.null(rownames(res[[1L]])))
    res <- .h5_read_data(h5f, index = 1, name = "chrom_peaks", ms_level = 2L,
                         read_rownames = TRUE)
    expect_equal(unname(res[[1L]]), unname(a2))
    expect_equal(rownames(res[[1L]]), rownames(a2))
    expect_true(is.null(colnames(res[[1L]])))
    ## single column
    res <- .h5_read_data(h5f, index = c(2, 1), name = "chrom_peaks",
                         ms_level = c(2L, 2L), column = 2)
    expect_equal(length(res), 2L)
    expect_true(ncol(res[[1L]]) == 1L)
    expect_equal(res[[1L]][, 1], unname(b2[, 2]))
    res <- .h5_read_data(h5f, index = c(2, 1), name = "chrom_peaks",
                         ms_level = c(2L, 2L), column = 2, read_colnames = TRUE,
                         read_rownames = TRUE)
    expect_equal(length(res), 2L)
    expect_true(ncol(res[[1L]]) == 1L)
    expect_equal(res[[1L]][, 1, drop = FALSE], b2[, 2, drop = FALSE])
    res <- .h5_read_data(h5f, index = c(1, 2, 1), name = "chrom_peaks",
                         ms_level = c(2, 2, 2), column = 1L,
                         read_colnames = TRUE,
                         read_rownames = TRUE)
    expect_equal(length(res), 3)
    expect_equal(res[[1]], res[[3]])
    expect_equal(res[[2]], b2[, 1, drop = FALSE])

    ## chrom peak data
    res <- .h5_read_data(h5f, index = c(2, 1), name = "chrom_peak_data",
                         ms_level = c(2L, 2L), read_colnames = TRUE,
                         read_rownames = TRUE)
    expect_equal(length(res), 2)
    rownames(b) <- c("CP3", "CP4", "CP5")
    expect_equal(unname(res[[1L]]), unname(b))
    res <- .h5_read_data(h5f, index = 1, name = "chrom_peak_data",
                         ms_level = 2L, column = "is_filled")
    expect_equal(length(res), 1L)
    expect_equal(res[[1L]][, 1], a$is_filled)

    file.remove(h5f)
})

test_that(".h5_compression_level works", {
    expect_equal(.h5_compression_level(), 0L)
})

test_that(".h5_initialize_file", {
    h5f <- tempfile()
    .h5_initialize_file(h5f, mod_count = 10L)
    expect_identical(h5read(h5f, "/header/modcount")[1L], 10L)
    expect_error(.h5_initialize_file(h5f), "already exists")
    file.remove(h5f)
})

test_that(".h5_increment_mod_count works", {
    h5f <- tempfile()
    .h5_initialize_file(h5f)
    h5 <- H5Fopen(h5f)
    res <- .h5_increment_mod_count(h5)
    expect_equal(res, 1L)
    res <- .h5_increment_mod_count(h5)
    expect_equal(res, 2L)
    H5Fclose(h5)
    file.remove(h5f)
})

test_that("HDF5 validity works", {
    h5f <- tempfile()
    h5 <- H5Fcreate(h5f)
    expect_error(.h5_valid_file(h5f), "not in correct format")
    h5createGroup(h5, "header")
    h5write("package:other", h5, "/header/package", level = 0L)
    expect_error(.h5_valid_file(h5f), "not in correct format")
    H5Fclose(h5)
    file.remove(h5f)

    .h5_initialize_file(h5f)
    expect_true(.h5_valid_file(h5f))

    h5 <- H5Fopen(h5f)
    expect_true(.h5_check_mod_count(h5, 0L))
    expect_error(.h5_check_mod_count(h5, 1L), "changed by a")
    H5Fclose(h5)
    file.remove(h5f)
})

test_that(".h5_write_chrom_peaks works", {
    h5f <- tempfile()
    .h5_initialize_file(h5f)
    h5 <- H5Fopen(h5f)

    mat <- cbind(mz = c(1.12, 1.34, 43.4), rt = c(23.2, 124.3, 123.5))
    rownames(mat) <- c("CP01", "CP02", "CP03")

    .h5_write_chrom_peaks(mat, h5, "test_1", 0L)
    expect_equal(h5read(h5, "test_1"), unname(mat))
    expect_equal(as.vector(h5read(h5, "test_1_rownames")), rownames(mat))
    expect_equal(as.vector(h5read(h5, "test_1_colnames")), colnames(mat))
    l <- h5ls(h5)
    expect_true(any(l$name == "test_1_colnames"))
    expect_true(any(l$name == "test_1_rownames"))

    .h5_write_chrom_peaks(mat, h5, "test_2", 0L, FALSE, FALSE)
    l <- h5ls(h5)
    expect_false(any(l$name == "test_2_colnames"))
    expect_false(any(l$name == "test_2_rownames"))
    H5Fclose(h5)
    file.remove(h5f)
})

test_that(".h5_write_chrom_peak_data works", {
    h5f <- tempfile()
    .h5_initialize_file(h5f)
    h5 <- H5Fopen(h5f)

    df <- data.frame(ms_level = c(1L, 2L), is_filled = FALSE)

    .h5_write_chrom_peak_data(df, h5, "test_1", 0L)
    res <- h5read(h5, "test_1")
    expect_true(is.list(res))
    expect_equal(names(res), c("is_filled", "ms_level"))
    expect_equal(as.vector(res$ms_level), unname(df$ms_level))
    expect_equal(as.vector(res$is_filled), unname(df$is_filled))
    res <- as.data.frame(res)
    expect_equal(res[, colnames(df)], df)
    H5Fclose(h5)
    file.remove(h5f)
})

test_that(".h5_write_data works", {
    h5f <- tempfile()
    .h5_initialize_file(h5f)

    a <- cbind(a = c(1.2, 1.4), b = c(3.5, 3.6), c = c(5.3, 5.1))
    rownames(a) <- c("CP1", "CP2")
    b <- cbind(a = c(12.4, 13.1, 3.2), b = c(1.3, 1.3, 1.4), c(4.2, 5.1, 4.6))
    rownames(b) <- c("CP3", "CP4", "CP5")
    l <- list(a, b)
    names(l) <- c(1, 2)

    ## chrom peaks
    expect_equal(.h5_write_data(h5f, l, name = "chrom_peaks",
                                ms_level = c(1L, 1L)), 1L)
    res <- h5read(h5f, "/1/ms_1/chrom_peaks")
    expect_equal(res, unname(a))
    res <- h5read(h5f, "/2/ms_1/chrom_peaks")
    expect_equal(res, unname(b))
    ## Update the first data set.
    a[1, 1] <- 10.4
    expect_equal(.h5_write_data(h5f, list(`1` = a), name = "chrom_peaks",
                                ms_level = 1L, replace = FALSE), 2L)
    res <- h5read(h5f, "/1/ms_1/chrom_peaks")
    expect_equal(res, unname(a))

    ## chrom peak data
    a <- data.frame(is_filled = c(TRUE, FALSE))
    b <- data.frame(is_filled = c(FALSE, FALSE, TRUE))
    l <- list(a, b)
    names(l) <- 1:2
    expect_equal(.h5_write_data(h5f, l, name = "chrom_peak_data",
                                ms_level = c(1L, 1L)), 3L)
    res <- h5read(h5f, "/1/ms_1/chrom_peak_data")
    expect_equal(a, as.data.frame(res))

    file.remove(h5f)
})
