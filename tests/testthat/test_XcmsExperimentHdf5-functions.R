library(rhdf5)
xmse_h5 <- xcms:::.xcms_experiment_to_hdf5(loadXcmsData("faahko_sub2"), tempfile())

## correspondence
h5f_full_g <- tempfile()
xmseg_full_h5 <- loadXcmsData("xmse") |>
    dropFeatureDefinitions() |>
    xcms:::.xcms_experiment_to_hdf5(h5f_full_g)
pdp <- PeakDensityParam(sampleGroups = sampleData(xmseg_full_h5)$sample_group,
                        minFraction = 0.4, bw = 30)
xmseg_full_h5 <- groupChromPeaks(xmseg_full_h5, pdp, msLevel = 1L)

test_that(".h5_chrom_peaks works", {
    res <- .h5_chrom_peaks(xmse_h5)
    expect_true(is.list(res))
    expect_true(length(res) == 0)

    res <- .h5_chrom_peaks(xmse_h5, msLevel = 1L)
    expect_true(is.list(res))
    expect_true(length(res) == 3)
    expect_equal(names(res), xmse_h5@sample_id)

    res <- .h5_chrom_peaks(xmse_h5, msLevel = 1L, by_sample = FALSE)
    expect_true(is.matrix(res))
    expect_true(any(colnames(res) == "sample"))
    expect_true(is.numeric(res))

    res_2 <- .h5_chrom_peaks(xmse_h5, msLevel = 1L, by_sample = FALSE,
                             columns = c("mz", "rt", "maxo"))
    expect_true(is.matrix(res_2))
    expect_true(is.numeric(res_2))
    expect_equal(colnames(res_2), c("mz", "rt", "maxo", "sample"))
    expect_equal(res[, colnames(res_2)], res_2)

    expect_error(
        .h5_chrom_peaks(xmse_h5, msLevel = 1L, by_sample = FALSE,
                        columns = c("mz", "rt", "maxo", "other")),
        "not found")

    ## with rt
    res <- .h5_chrom_peaks(xmse_h5, rt = c(2500, 2700), type = "apex_within",
                           msLevel = 1L, by_sample = FALSE)
    expect_true(all(res[, "rt"] > 2500))
    expect_true(all(res[, "rt"] < 2700))

    ## with mz and rt
    res <- .h5_chrom_peaks(xmse_h5, rt = c(2500, 2700), type = "apex_within",
                           msLevel = 1L, by_sample = FALSE, mz = c(400, 600))
    expect_true(all(res[, "rt"] > 2500))
    expect_true(all(res[, "rt"] < 2700))
    expect_true(all(res[, "mz"] > 400))
    expect_true(all(res[, "mz"] < 600))
})

test_that(".h5_chrom_peak_data works", {
    cp <- chromPeaks(xmse_h5)
    res <- .h5_chrom_peak_data(xmse_h5, msLevel = 1L, by_sample = FALSE)
    expect_true(is.data.frame(res))
    expect_equal(nrow(cp), nrow(res))
    expect_equal(colnames(res), c("is_filled", "ms_level"))

    res <- .h5_chrom_peak_data(xmse_h5, msLevel = 1L, by_sample = TRUE,
                               peaks = rownames(cp)[c(4, 10, 100, 200)])
    expect_true(is.list(res))
    expect_equal(names(res), xmse_h5@sample_id)
    expect_equal(rownames(res[[1L]]), rownames(cp)[c(4, 10)])
    expect_equal(rownames(res[[2L]]), rownames(cp)[100])
    expect_equal(rownames(res[[3L]]), rownames(cp)[200])
    res <- .h5_chrom_peak_data(xmse_h5, msLevel = 1L, by_sample = FALSE,
                               peaks = rownames(cp)[c(4, 10, 100, 200)])
    expect_equal(rownames(res), rownames(cp)[c(4, 10, 100, 200)])
})

test_that(".h5_chrom_peak_data_colnames works", {
    res <- .h5_chrom_peak_data_colnames(xmse_h5)
    expect_equal(res, c("is_filled", "ms_level"))
})

test_that(".h5_chrom_peaks_chunk works", {
    p <- CentWaveParam(noise = 10000, snthresh = 40, prefilter = c(3, 10000))
    sps <- spectra(xmse_h5)[xmse_h5@sampleDataLinks[["spectra"]][, 2]]
    sps$.SAMPLE_IDX <- xmse_h5@sampleDataLinks[["spectra"]][, 1L]
    sps <- sps[sps$.SAMPLE_IDX %in% 2:3]

    h5_file <- tempfile()
    .h5_initialize_file(h5_file)

    res <- .h5_find_chrom_peaks_chunk(
        sps, msLevel = 1L, param = p, h5_file = h5_file,
        sample_id = xmse_h5@sample_id)
    expect_true(is.integer(res))
    expect_equal(res, 2L)
    h5 <- rhdf5::H5Fopen(h5_file)
    expect_equal(.h5_dataset_names("/", h5), c("S2", "S3", "header"))
    expect_equal(.h5_dataset_names("/S2/ms_1/", h5),
                 c("chrom_peak_data", "chrom_peaks",
                   "chrom_peaks_colnames", "chrom_peaks_rownames"))
    H5Fclose(h5)
    a <- .h5_read_data(h5_file, id = "S2", name = "chrom_peaks",
                       ms_level = 1L, read_colnames = TRUE,
                       read_rownames = TRUE)[[1L]]
    ## add = TRUE
    res <- .h5_find_chrom_peaks_chunk(
        sps, msLevel = 1L, param = p, h5_file = h5_file, add = TRUE,
        sample_id = xmse_h5@sample_id)
    expect_equal(res, 4L)
    b <- .h5_read_data(h5_file, id = "S2", name = "chrom_peaks",
                       ms_level = 1L, read_colnames = TRUE,
                       read_rownames = TRUE)[[1L]]
    expect_equal(nrow(b), 2 * nrow(a))
    expect_equal(anyDuplicated(rownames(b)), 0L)
    same <- rownames(b) %in% rownames(a)
    expect_equal(a, b[same, ])
    expect_equal(unname(a), unname(b[!same, ]))
    unlink(tempfile())
})

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

    ## With feature definitions.
    ref <- loadXcmsData("xmse")
    h5f <- tempfile()

    res <- .xcms_experiment_to_hdf5(ref, h5f)
    expect_true(validObject(res))
    expect_s4_class(res, "XcmsExperimentHdf5")
    expect_true(hasFeatures(res))
    expect_equal(res@features_ms_level, 1L)
    a <- featureDefinitions(ref)
    b <- featureDefinitions(res)
    rownames(a) <- NULL
    rownames(b) <- NULL
    expect_equal(a[, colnames(b)], b)

    a <- featureValues(ref, method = "sum")
    b <- featureValues(res, method = "sum")
    rownames(a) <- NULL
    rownames(b) <- NULL
    expect_equal(a, b)

    expect_equal(hasFilledChromPeaks(ref), hasFilledChromPeaks(res))
    a <- featureValues(ref, method = "sum", filled = FALSE)
    b <- featureValues(res, method = "sum", filled = FALSE)
    rownames(a) <- NULL
    rownames(b) <- NULL
    expect_equal(a, b)

    file.remove(h5f)
})

test_that("toXcmsExperimentHdf5 works", {
    res <- toXcmsExperimentHdf5(XcmsExperiment())
    expect_true(validObject(res))
    expect_s4_class(res, "XcmsExperimentHdf5")

    expect_error(toXcmsExperimentHdf5(5), "'XcmsExperiment'")

    tf <- tempfile()
    ref <- loadXcmsData("faahko_sub2")
    res <- toXcmsExperimentHdf5(ref, tf)
    expect_s4_class(res, "XcmsExperimentHdf5")
    rm(tf)
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
    ref <- .h5_read_data(x@hdf5_file, id = x@sample_id,
                         ms_level = rep(1L, length(x)),
                         read_colnames = TRUE, read_rownames = TRUE)
    .h5_xmse_merge_neighboring_peaks(x)
    mod_count <- .h5_mod_count(h5f)
    expect_true(mod_count > x@hdf5_mod_count)
    ## Check that content was changed.
    res <- .h5_read_data(x@hdf5_file, id = x@sample_id,
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
    res <- .h5_read_data(res@hdf5_file, id = res@sample_id,
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

test_that(".h5_read_matrix works", {
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
    res <- .h5_read_matrix("/1/ms_1/chrom_peaks", h5)
    expect_equal(res, unname(a))
    res <- .h5_read_matrix("/1/ms_1/chrom_peaks", h5,
                           read_colnames = TRUE)
    expect_equal(unname(res), unname(a))
    expect_equal(colnames(res), colnames(a))
    expect_equal(rownames(res), NULL)
    res2 <- .h5_read_matrix2("/1/ms_1/chrom_peaks", h5,
                             read_colnames = TRUE)
    expect_equal(res, res2)
    res <- .h5_read_matrix("/1/ms_1/chrom_peaks", h5,
                           read_colnames = TRUE, read_rownames = TRUE)
    expect_equal(res, a)
    res2 <- .h5_read_matrix2("/1/ms_1/chrom_peaks", h5,
                             read_colnames = TRUE, read_rownames = TRUE)
    expect_equal(res, res2)
    res <- .h5_read_matrix("/1/ms_1/chrom_peaks", h5, index = list(NULL, 3),
                           read_colnames = TRUE, read_rownames = TRUE)
    expect_equal(res, a[, 3, drop = FALSE])
    res2 <- .h5_read_matrix2("/1/ms_1/chrom_peaks", h5, index = list(NULL, 3),
                             read_colnames = TRUE, read_rownames = TRUE)
    expect_equal(res, res2)
    res <- .h5_read_matrix("/1/ms_1/chrom_peaks", h5,
                           index = list(NULL, c(1, 3)),
                           read_colnames = FALSE, read_rownames = FALSE)
    expect_equal(res, unname(a[, c(1, 3)]))
    res2 <- .h5_read_matrix2("/1/ms_1/chrom_peaks", h5,
                             index = list(NULL, c(1, 3)),
                             read_colnames = FALSE, read_rownames = FALSE)
    expect_equal(res, res2)
    res <- .h5_read_matrix("/1/ms_1/chrom_peaks", h5, index = list(2, c(1, 3)),
                           read_colnames = TRUE, read_rownames = TRUE)
    expect_equal(res, a[2, c(1, 3), drop = FALSE])
    res2 <- .h5_read_matrix2("/1/ms_1/chrom_peaks", h5, index = list(2, c(1, 3)),
                             read_colnames = TRUE, read_rownames = TRUE)
    expect_equal(res, res2)
    res <- .h5_read_matrix("/1/ms_1/chrom_peaks", h5,
                           index = list(c(2, 1, 2), c(1, 3)),
                           read_colnames = TRUE, read_rownames = TRUE)
    expect_equal(res, a[c(2, 1, 2), c(1, 3), drop = FALSE])
    res2 <- .h5_read_matrix2("/1/ms_1/chrom_peaks", h5,
                             index = list(c(2, 1, 2), c(1, 3)),
                             read_colnames = TRUE, read_rownames = TRUE)
    expect_equal(res, res2)
    res <- .h5_read_matrix("/1/ms_1/chrom_peaks", h5,
                           index = list(c(2, 1, 2), c(1, 3)),
                           read_colnames = FALSE, read_rownames = FALSE)
    expect_equal(res, unname(a[c(2, 1, 2), c(1, 3), drop = FALSE]))
    res2 <- .h5_read_matrix2("/1/ms_1/chrom_peaks", h5,
                             index = list(c(2, 1, 2), c(1, 3)),
                             read_colnames = FALSE, read_rownames = FALSE)
    expect_equal(res, res2)

    H5Fclose(h5)
    file.remove(h5f)
})

test_that(".h5_read_chrom_peaks_matrix works", {
    res <- .h5_read_chrom_peaks_matrix(
        "/S2/ms_1/chrom_peaks", xmse_h5@hdf5_file,
        read_colnames = FALSE, read_rownames = FALSE)
    expect_true(is.matrix(res))
    expect_true(is.numeric(res))
    res <- .h5_read_chrom_peaks_matrix(
        "/S2/ms_1/chrom_peaks", xmse_h5@hdf5_file,
        read_colnames = TRUE, read_rownames = FALSE,
        mz = c(300, 350), rt = c(3000, 3500), type = "within")
    expect_true(all(res[, "mz"] > 300 & res[, "mz"] < 350))
    expect_true(all(res[, "rt"] > 3000 & res[, "rt"] < 3500))

    res <- .h5_read_chrom_peaks_matrix(
        "/S2/ms_1/chrom_peaks", xmse_h5@hdf5_file,
        read_colnames = TRUE, read_rownames = FALSE,
        mz = c(300, 350), type = "within")
    expect_true(all(res[, "mz"] > 300 & res[, "mz"] < 350))

    ## with sample_index
    sidx <- c(4L, 5L, 6L)
    names(sidx) <- c("/S1/ms_1/chrom_peaks", "/S2/ms_1/chrom_peaks",
                     "/S3/ms_1/chrom_peaks")
    res <- .h5_read_chrom_peaks_matrix(
        "/S2/ms_1/chrom_peaks", xmse_h5@hdf5_file, read_colnames = TRUE,
        read_rownames = FALSE, sample_index = sidx)
    expect_true(is.matrix(res))
    expect_true(is.numeric(res))
    expect_true(nrow(res) > 0)
    expect_true(any(colnames(res) %in% "sample"))
    expect_true(all(res[, "sample"] == 5))
})

test_that(".h5_read_data_frame works", {
    h5f <- tempfile()
    .h5_initialize_file(h5f)

    a <- cbind(a = c(1.2, 1.4), b = c(3.5, 3.6), c = c(5.3, 5.1))
    rownames(a) <- c("CP1", "CP2")
    b <- cbind(a = c(12.4, 13.1, 3.2), b = c(1.3, 1.3, 1.4), c(4.2, 5.1, 4.6))
    rownames(b) <- c("CP3", "CP4", "CP5")
    l <- list(a, b)
    names(l) <- c(1, 2)
    .h5_write_data(h5f, l, name = "chrom_peaks", ms_level = c(2L, 2L))
    a <- data.frame(is_filled = c(TRUE, FALSE), other_col = "c")
    b <- data.frame(is_filled = c(FALSE, FALSE, TRUE), other_col = "d")
    l <- list(a, b)
    names(l) <- 1:2
    .h5_write_data(h5f, l, name = "chrom_peak_data", ms_level = c(2L, 2L))

    h5 <- H5Fopen(h5f)
    res <- .h5_read_chrom_peak_data("/1/ms_2/chrom_peak_data", h5,
                                    read_rownames = TRUE)
    expect_equal(rownames(res), c("CP1", "CP2"))
    expect_equal(colnames(res), c("is_filled", "other_col"))
    rownames(res) <- NULL
    expect_equal(res, a)
    res <- .h5_read_data_frame("/1/ms_2/chrom_peak_data", h5,
                               read_rownames = FALSE)
    expect_equal(colnames(res), c("is_filled", "other_col"))
    expect_equal(rownames(res), c("1", "2"))
    expect_equal(res, a)

    ## Read selected rows
    res <- .h5_read_chrom_peak_data("/1/ms_2/chrom_peak_data", h5,
                                    read_rownames = FALSE,
                                    index = list(2, NULL))
    expect_equal(res, a[2, ])
    res <- .h5_read_chrom_peak_data("/1/ms_2/chrom_peak_data", h5,
                                    read_rownames = FALSE,
                                    index = list(c(2, 1, 2), NULL))
    expect_equal(res, a[c(2, 1, 2), ])
    res <- .h5_read_chrom_peak_data("/1/ms_2/chrom_peak_data", h5,
                                    read_rownames = TRUE,
                                    index = list(c(2, 1, 2), NULL))
    expect_equal(rownames(res), c("CP2", "CP1", "CP2.1"))

    ## Read single column
    res <- .h5_dataset_names("/1/ms_2/chrom_peak_data", h5)
    expect_equal(res, c("is_filled", "other_col"))
    res <- .h5_read_chrom_peak_data("/1/ms_2/chrom_peak_data/other_col", h5)
    expect_true(is.data.frame(res))
    expect_equal(res[, 1L], c("c", "c"))

    ## Read selected rows
    res <- .h5_read_chrom_peak_data("/1/ms_2/chrom_peak_data/is_filled", h5,
                                    index = list(2, NULL))
    expect_true(nrow(res) == 1L)
    expect_equal(res[, 1L], FALSE)

    ## Read selected rows using chromPeak IDs.
    cp <- chromPeaks(xmse_h5)
    res <- .h5_read_chrom_peak_data(
        "/S1/ms_1/chrom_peak_data", xmse_h5@hdf5_file,
        peaks = rownames(cp)[4:10])
    expect_equal(nrow(res), 7)
    expect_equal(rownames(res), as.character(1:7))

    res <- .h5_read_chrom_peak_data(
        "/S1/ms_1/chrom_peak_data", xmse_h5@hdf5_file,
        peaks = rownames(cp)[c(4, 20, 100, 200)], read_rownames = TRUE)
    expect_equal(nrow(res), 2L)
    expect_equal(rownames(res), rownames(cp)[c(4, 20)])

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
    res <- .h5_read_data(h5f, id = 2, name = "chrom_peaks", ms_level = 2L)
    expect_equal(length(res), 1L)
    expect_equal(res[[1L]], unname(b2))
    res <- .h5_read_data(h5f, id = 1, name = "chrom_peaks", ms_level = 2L,
                         read_colnames = TRUE)
    expect_equal(unname(res[[1L]]), unname(a2))
    expect_equal(colnames(res[[1L]]), colnames(a2))
    expect_true(is.null(rownames(res[[1L]])))
    res <- .h5_read_data(h5f, id = 1, name = "chrom_peaks", ms_level = 2L,
                         read_rownames = TRUE)
    expect_equal(unname(res[[1L]]), unname(a2))
    expect_equal(rownames(res[[1L]]), rownames(a2))
    expect_true(is.null(colnames(res[[1L]])))
    ## single column
    res <- .h5_read_data(h5f, id = c(2, 1), name = "chrom_peaks",
                         ms_level = c(2L, 2L), j = 2)
    expect_equal(length(res), 2L)
    expect_true(ncol(res[[1L]]) == 1L)
    expect_equal(res[[1L]][, 1], unname(b2[, 2]))
    res <- .h5_read_data(h5f, id = c(2, 1), name = "chrom_peaks",
                         ms_level = c(2L, 2L), j = 2, read_colnames = TRUE,
                         read_rownames = TRUE)
    expect_equal(length(res), 2L)
    expect_true(ncol(res[[1L]]) == 1L)
    expect_equal(res[[1L]][, 1, drop = FALSE], b2[, 2, drop = FALSE])
    res <- .h5_read_data(h5f, id = c(1, 2, 1), name = "chrom_peaks",
                         ms_level = c(2, 2, 2), j = 1L,
                         read_colnames = TRUE,
                         read_rownames = TRUE)
    expect_equal(length(res), 3)
    expect_equal(res[[1]], res[[3]])
    expect_equal(res[[2]], b2[, 1, drop = FALSE])

    ## selected rows.
    res <- .h5_read_data(h5f, id = c(1, 2, 1), name = "chrom_peaks",
                         ms_level = c(2, 2, 2), j = 1L, i = c(2, 1, 2),
                         read_colnames = TRUE,
                         read_rownames = TRUE)
    expect_equal(length(res), 3)
    expect_equal(res[[1]], res[[3]])
    expect_equal(res[[2]], b2[c(2, 1, 2), 1, drop = FALSE])
    expect_equal(res[[1]], a2[c(2, 1, 2), 1, drop = FALSE])

    ## chrom peak data
    res <- .h5_read_data(h5f, id = c(2, 1), name = "chrom_peak_data",
                         ms_level = c(2L, 2L), read_colnames = TRUE,
                         read_rownames = TRUE)
    expect_equal(length(res), 2)
    rownames(b) <- c("CP3", "CP4", "CP5")
    expect_equal(unname(res[[1L]]), unname(b))
    res <- .h5_read_data(h5f, id = 1, name = "chrom_peak_data",
                         ms_level = 2L, j = "is_filled")
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

test_that(".h5_mod_count works", {
    h5f <- tempfile()
    .h5_initialize_file(h5f, mod_count = 7L)
    expect_identical(.h5_mod_count(h5f), 7L)
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
    expect_error(.h5_valid_file(), "missing with no default")
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

test_that(".h5_write_matrix works", {
    h5f <- tempfile()
    xcms:::.h5_initialize_file(h5f)
    h5 <- H5Fopen(h5f)

    mat <- cbind(mz = c(1.12, 1.34, 43.4), rt = c(23.2, 124.3, 123.5))
    rownames(mat) <- c("CP01", "CP02", "CP03")

    xcms:::.h5_write_matrix(mat, h5, "test_1", 0L)
    expect_equal(h5read(h5, "test_1"), unname(mat))
    expect_equal(as.vector(h5read(h5, "test_1_rownames")), rownames(mat))
    expect_equal(as.vector(h5read(h5, "test_1_colnames")), colnames(mat))
    l <- h5ls(h5)
    expect_true(any(l$name == "test_1_colnames"))
    expect_true(any(l$name == "test_1_rownames"))

    xcms:::.h5_write_matrix(mat, h5, "test_2", 0L, FALSE, FALSE)
    l <- h5ls(h5)
    expect_false(any(l$name == "test_2_colnames"))
    expect_false(any(l$name == "test_2_rownames"))

    xcms:::.h5_write_matrix(mat[integer(), , drop = FALSE], h5, "test_3", 0L,
                            TRUE, TRUE)
    l <- h5ls(h5)
    expect_true(any(l$name == "test_3_colnames"))
    expect_true(any(l$name == "test_3_rownames"))
    res <- h5read(h5, "/test_3")
    expect_true(is.matrix(res))
    expect_equal(ncol(res), 2)
    expect_equal(nrow(res), 0)
    res <- h5read(h5, "/test_3_colnames")
    expect_equal(as.vector(res), c("mz", "rt"))
    res <- h5read(h5, "/test_3_rownames")
    expect_equal(as.vector(res), character())

    H5Fclose(h5)
    file.remove(h5f)
})

test_that(".h5_write_data_frame works", {
    h5f <- tempfile()
    xcms:::.h5_initialize_file(h5f)
    h5 <- H5Fopen(h5f)

    df <- data.frame(ms_level = c(1L, 2L), is_filled = FALSE)

    .h5_write_data_frame(df, h5, "test_1", 0L)
    res <- h5read(h5, "test_1")
    expect_true(is.list(res))
    expect_equal(names(res), c("is_filled", "ms_level"))
    expect_equal(as.vector(res$ms_level), unname(df$ms_level))
    expect_equal(as.vector(res$is_filled), unname(df$is_filled))
    res <- as.data.frame(res)
    expect_equal(res[, colnames(df)], df)

    ## Write an empty data.frame
    df <- df[integer(), , drop = FALSE]
    .h5_write_data_frame(df, h5, "test_1", 0L)
    res <- h5read(h5, "test_1")
    expect_true(is.list(res))
    expect_equal(names(res), c("is_filled", "ms_level"))
    expect_equal(as.vector(res$ms_level), unname(df$ms_level))
    expect_equal(as.vector(res$is_filled), unname(df$is_filled))
    res <- as.data.frame(res)
    expect_equal(res[, colnames(df)], df)

    ## write an empty data.frame
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

test_that(".h5_feature_values_sample", {
    ## unit test in test_XcmsExperimentHdf5.R
    expect_true(TRUE)
})

test_that(".h5_feature_values_ms_level", {
    ## unit test in test_XcmsExperimentHdf5.R
    expect_true(TRUE)
})

test_that(".h5_feature_values_ms_level", {
    ## unit test in test_XcmsExperimentHdf5.R
    expect_true(TRUE)
})

test_that(".h5_chrom_peaks_colnames works", {
    res <- .h5_chrom_peaks_colnames(xmse_h5, 1L)
    expect_true(is.character(res))
    expect_true(all(c("mz", "mzmin", "mzmax", "rt") %in% res))
})

test_that(".h5_chrom_peaks_rownames works", {
    res <- .h5_chrom_peaks_rownames(xmse_h5)
    expect_true(is.list(res))
    expect_equal(length(res), 3L)
    expect_true(is.character(res[[1L]]))
})

test_that(".h5_filter works", {
    expect_equal(.h5_filter(), "NONE")
})

test_that(".h5_update_rt_chrom_peaks_sample works", {
    ## Unit test is in text_XcmsExperimentHdf5.R @adjustRtime,XcmsExperimentHdf5
    expect_true(TRUE)
})

test_that(".h5_x_chromatogram works", {
    rt <- matrix(c(2600, 2700), nrow = 1)
    mz <- matrix(c(340, 400), nrow = 1)

    res <- .h5_x_chromatograms(xmse_h5, mz = mz, rt = rt, chromPeaks = "any")
    expect_s4_class(res, "XChromatograms")
    expect_true(nrow(res) == 1L)
    expect_equal(ncol(res), length(xmse_h5))
    expect_true(nrow(chromPeaks(res)) > 0)
    expect_true(all(chromPeaks(res)[, "mz"] >= 340 &
                    chromPeaks(res)[, "mz"] <= 400))
    expect_true(all(chromPeaks(res[1, 1])[, "sample"] == 1L))
    expect_true(all(chromPeaks(res[1, 2])[, "sample"] == 2L))
    expect_true(all(chromPeaks(res[1, 3])[, "sample"] == 3L))
    ref <- chromatogram(loadXcmsData("faahko_sub2"), mz = mz, rt = rt,
                        chromPeaks = "any")
    expect_equal(unname(chromPeaks(res)), unname(chromPeaks(ref)))

    ## Multiple rows.
    res <- .h5_x_chromatograms(
        xmse_h5, mz = chromPeaks(xmse_h5)[1:10, c("mzmin", "mzmax")],
        rt = chromPeaks(xmse_h5)[1:10, c("rtmin", "rtmax")],
        ms_level = 1L,  chromPeaks = "apex_within", BPPARAM = bpparam())
    expect_true(nrow(res) == 10L)

    ## With features.
    rt <- rbind(rt, c(3000, 3500))
    mz <- rbind(mz, c(450, 500))
    tmp <- loadXcmsData("xmse") |>
        dropFeatureDefinitions() |>
        groupChromPeaks(pdp, msLevel = 1L)

    ref <- chromatogram(tmp, mz = mz, rt = rt, chromPeaks = "any")
    res <- .h5_x_chromatograms(xmseg_full_h5, mz = mz, rt = rt,
                               chromPeaks = "any")
    a <- chromPeaks(ref)
    b <- chromPeaks(res)
    expect_equal(colnames(a), colnames(b))
    rownames(a) <- NULL
    rownames(b) <- NULL
    expect_equal(a, b)

    a <- featureDefinitions(ref)
    b <- featureDefinitions(res)
    expect_true(all(colnames(a) %in% colnames(b)))
    ## We don't have exactly the same number of features, because XcmsExperiment
    ## selects features based on mz and rt range, not based on included chrom
    ## peaks.
    expect_true(all(a$mzmed %in% b$mzmed))
    b <- b[b$mzmed %in% a$mzmed, ]
    expect_equal(a$peakidx, b$peakidx)

    ## With overlapping ranges -> duplicated features and chrom peaks.
    mz <- rbind(mz, mz[1, ])
    rt <- rbind(rt, rt[1, ])
    res <- .h5_x_chromatograms(xmseg_full_h5, mz = mz, rt = rt,
                               chromPeaks = "apex_within")
    cp <- chromPeaks(res)
    expect_equal(cp[cp[, "row"] == 1L, colnames(cp) != "row"],
                 cp[cp[, "row"] == 3L, colnames(cp) != "row"])
    fts <- featureDefinitions(res)
})

test_that(".h5_features_ms_region_values works", {
    tmpl <- rep(NA_real_, 20)
    v <- as.numeric(1:14)
    m <- c(2, 2, 3, 4, 5, 6, 7, 8, 8, 9, 9, 9, 10, 11)
    res <- .h5_features_ms_region_values(tmpl, v, m)
    expect_true(is.na(res[1L]))
    expect_true(all(is.na(res[12:20])))
    expect_equal(res[2:11], c(1, 3, 4, 5, 6, 7, 8, 10, 13, 14))
    res <- .h5_features_ms_region_values(tmpl, v, m, max)
    expect_true(is.na(res[1L]))
    expect_true(all(is.na(res[12:20])))
    expect_equal(res[2:11], c(2, 3, 4, 5, 6, 7, 9, 12, 13, 14))
})

test_that(".h5_features_ms_region works", {
    res <- .h5_features_ms_region(
        xmseg_full_h5, mzmin = min, mzmax = max, rtmin = min, rtmax = max,
        features = rownames(featureDefinitions(xmseg_full_h5)), ms_level = 1L)
    expect_true(is.matrix(res))
    expect_equal(colnames(res), c("mzmin", "mzmax", "rtmin", "rtmax"))
    expect_equal(nrow(res), nrow(featureDefinitions(xmseg_full_h5)))
    expect_equal(rownames(res), rownames(featureDefinitions(xmseg_full_h5)))

    res_sub <- .h5_features_ms_region(
        xmseg_full_h5, mzmin = min, mzmax = max, rtmin = min, rtmax = max,
        features = rownames(res)[c(4, 10, 12)], ms_level = 1L)
    expect_true(is.matrix(res_sub))
    expect_equal(colnames(res), c("mzmin", "mzmax", "rtmin", "rtmax"))
    expect_equal(nrow(res_sub), 3)
    expect_equal(res[c(4, 10, 12), ], res_sub)

    res_sub <- .h5_features_ms_region(
        xmseg_full_h5, mzmin = min, mzmax = max, rtmin = min, rtmax = max,
        features = rownames(res)[c(10, 12, 10, 4)], ms_level = 1L)
    expect_true(is.matrix(res_sub))
    expect_equal(colnames(res), c("mzmin", "mzmax", "rtmin", "rtmax"))
    expect_equal(nrow(res_sub), 4)
    expect_equal(res[c(10, 12, 10, 4), ], res_sub)

    ## errors
    expect_error(.h5_features_ms_region(
        xmseg_full_h5, mzmin = min, mzmax = max, rtmin = min, rtmax = max,
        features = c(rownames(res)[c(10, 12, 10, 4)], "a"), ms_level = 1L),
        "not found")
})

test_that(".h5_xmse_integrate_chrom_peaks works", {
    tf <- tempfile()
    file.copy(xmseg_full_h5@hdf5_file, tf)
    tmp <- xmseg_full_h5
    tmp@hdf5_file <- tf

    x <- .h5_subset_xcms_experiment(
        tmp, 1:2, keepChromPeaks = TRUE, keepAdjustedRtime = TRUE,
        keepFeatures = TRUE, ignoreHistory = TRUE)
    pal <- split.data.frame(chromPeaks(x), chromPeaks(x)[, "sample"])
    res <- .h5_xmse_integrate_chrom_peaks(x, pal, param = CentWaveParam())
    expect_equal(res, 22L)
    ## check file content.
    a <- .h5_read_data(x@hdf5_file, x@sample_id, name = "chrom_peaks",
                       ms_level = rep(1L, 2), read_rownames = TRUE,
                       read_colnames = TRUE)
    expect_equal(nrow(a[[1L]]), 2 * nrow(pal[[1L]]))
    expect_equal(nrow(a[[2L]]), 2 * nrow(pal[[2L]]))
    expect_equal(a[[1L]][seq_len(nrow(pal[[1L]])), ],
                 pal[[1L]][, colnames(pal[[1L]]) != "sample"])
    expect_equal(a[[2L]][seq_len(nrow(pal[[2L]])), ],
                 pal[[2L]][, colnames(pal[[2L]]) != "sample"])
    a <- .h5_read_data(x@hdf5_file, x@sample_id, name = "chrom_peak_data",
                       ms_level = rep(1L, 2))
    expect_equal(nrow(a[[1L]]), 2 * nrow(pal[[1L]]))
    expect_equal(nrow(a[[2L]]), 2 * nrow(pal[[2L]]))
    expect_true(all(!a[[1L]]$is_filled[seq_len(nrow(pal[[1L]]))]))
    expect_true(all(a[[1L]]$is_filled[seq(nrow(pal[[1L]]) + 1L,
                                          length.out = nrow(pal[[1L]]))]))

    ## Test with two other samples, get chrom peak area for some features
    ## present in a file. then integrate. check that featureValues are doubled
    ## with sum.
    tmp@hdf5_mod_count <- 22L
    x <- .h5_subset_xcms_experiment(
        tmp, 3:4, keepChromPeaks = TRUE, keepAdjustedRtime = TRUE,
        keepFeatures = TRUE, ignoreHistory = TRUE)
    fvals_a <- featureValues(x, method = "sum")
    idx <- c(2, 4, 5, 6)
    fvals <- fvals_a[idx, ]
    a <- .h5_read_data(x@hdf5_file, id = x@sample_id,
                       "feature_to_chrom_peaks", ms_level = c(1L, 1L))
    pks <- .h5_read_data(x@hdf5_file, id = x@sample_id, "chrom_peaks",
                         ms_level = c(1L, 1L), read_colnames = TRUE)
    pal <- mapply(pks, a, FUN = function(z, i) {
        z[i[i[, 1L] %in% idx, 2L], c("mzmin", "mzmax", "rtmin", "rtmax")]
    }, SIMPLIFY = FALSE, USE.NAMES = FALSE)
    pal <- lapply(pal, function(z) {
        rownames(z) <- rownames(fvals_a)[idx]
        z
    })
    names(pal) <- seq_along(pal)
    res <- .h5_xmse_integrate_chrom_peaks(
        x, pal, param = CentWaveParam(), BPPARAM = SerialParam(),
        update_features = TRUE)
    expect_equal(res, 28L)
    b <- .h5_read_data(x@hdf5_file, id = x@sample_id, "feature_to_chrom_peaks",
                       ms_level = c(1L, 1L))
    expect_equal(vapply(a, nrow, NA_integer_) + 4L,
                 vapply(b, nrow, NA_integer_))

    fvals_b <- featureValues(x, method = "sum")
    expect_equal(fvals_a[-idx, ], fvals_b[-idx, ])
    expect_true(all(fvals_b[idx, ] > fvals_a[idx, ]))
    rm(tf)
})

test_that(".h5_feature_definitions_rownames works", {
    res <- .h5_feature_definitions_rownames(xmseg_full_h5)
    expect_true(is.list(res))
    expect_equal(res[[1L]], rownames(featureDefinitions(xmseg_full_h5)))
})

test_that(".h5_filter_chrom_peaks,XcmsExperimentHdf5 works", {
    tf <- tempfile()
    file.copy(xmse_h5@hdf5_file, tf)
    x <- xmse_h5
    x@hdf5_file <- tf
    res <- .h5_filter_chrom_peaks(x, 1L, .which_in_range,
                                  range = c(-Inf, Inf), column = "rt")
    expect_equal(res, x@hdf5_mod_count)
    res <- .h5_filter_chrom_peaks(x, 1L, .which_in_range,
                                  range = c(3000, 3300), column = "rt")
    expect_equal(res, x@hdf5_mod_count + 6)
    x@hdf5_mod_count <- res
    pks <- chromPeaks(x)
    expect_true(all(pks[, "rt"] > 3000 & pks[, "rt"] < 3300))
    rm(tf)

    ## With features.
    tf <- tempfile()
    file.copy(xmseg_full_h5@hdf5_file, tf)
    x <- xmseg_full_h5
    x@hdf5_file <- tf
    res <- .h5_filter_chrom_peaks(
        x, 1L, .which_in_range, range = c(3000, 3300), column = "rt")
    expect_equal(res, x@hdf5_mod_count * 2)
    fts <- .h5_read_data(
        x@hdf5_file, "features", "feature_definitions", 1L)[[1L]]
    expect_true(all(fts[, "rtmax"] > 3000))
    expect_true(all(fts[, "rtmin"] < 3300))
    map <- .h5_read_data(x@hdf5_file, x@sample_id, "feature_to_chrom_peaks",
                         rep(1L, length(x)))
    lapply(map, function(z) {
        expect_true(!anyNA(z))
        expect_true(all(z[, 1L] %in% seq_len(nrow(fts))))
    })
})

test_that(".h5_to_xcms_experiment works", {
    res <- toXcmsExperiment(xmse_h5)
    expect_equal(chromPeaks(res), chromPeaks(xmse_h5))
    expect_equal(rtime(res), rtime(xmse_h5))
    expect_equal(res@processHistory, xmse_h5@processHistory)

    res <- xcms:::.h5_to_xcms_experiment(xmseg_full_h5)
    expect_equal(chromPeaks(res), chromPeaks(xmseg_full_h5))
    expect_equal(rtime(res), rtime(xmseg_full_h5))
    expect_equal(res@processHistory, xmseg_full_h5@processHistory)
    expect_equal(chromPeakData(res), chromPeakData(xmseg_full_h5))
    a <- featureDefinitions(res)
    b <- featureDefinitions(xmseg_full_h5)
    expect_equal(a[, colnames(b)], b)
    a <- featureValues(res, method = "sum")
    b <- featureValues(xmseg_full_h5, method = "sum")
    expect_equal(a, b)
})

test_that(".h5_chrom_peak_spectra_sample works", {
    s <- spectra(xmse_h5[2L])
    cp <- chromPeaks(xmse_h5)
    cp2 <- cp[cp[, "sample"] == 2L, ]

    all <- xcms:::.h5_chrom_peak_spectra_sample(xmse_h5@hdf5_file, "S2", s,
                                                "all", msLevel = 1L, ppm = 10,
                                                expandMz = 0.01, expandRt = 1,
                                                skipFilled = TRUE)
    expect_s4_class(all, "Spectra")
    expect_true(all(rownames(cp)[cp[, "sample"] == 2] %in% all$chrom_peak_id))

    ## peaks as character
    peaks <- rownames(cp)[1:4]
    sel <- xcms:::.h5_chrom_peak_spectra_sample(xmse_h5@hdf5_file, "S2", s,
                                                "closest_rt", msLevel = 1L,
                                                skipFilled = TRUE,
                                                peaks = peaks)
    expect_s4_class(sel, "Spectra")
    expect_true(length(sel) == 0L)

    ## peaks as integer
    peaks <- c(6000, 60043) # should throw an error
    expect_error(.h5_chrom_peak_spectra_sample(
        xmse_h5@hdf5_file, "S2", s, "closest_rt", msLevel = 1L,
        skipFilled = TRUE, peaks = peaks), "out of bounds")

    peaks <- c(4, 2, 5, 2, 7)
    sel <- xcms:::.h5_chrom_peak_spectra_sample(xmse_h5@hdf5_file, "S2", s,
                                                "closest_rt", msLevel = 1L,
                                                skipFilled = TRUE,
                                                peaks = peaks)
    expect_s4_class(sel, "Spectra")
    expect_true(length(sel) == length(peaks))
    expect_equal(sel$chrom_peak_id, rownames(cp2)[peaks])
})

rm(h5f_full_g)
