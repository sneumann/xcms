h5f <- tempfile()
xmse_h5 <- xcms:::.xcms_experiment_to_hdf5(loadXcmsData("faahko_sub2"), h5f)
h5f_full <- tempfile()
a <- dropFeatureDefinitions(loadXcmsData("xmse"))
xmse_full_h5 <- xcms:::.xcms_experiment_to_hdf5(a, h5f_full)

## correspondence
h5f_full_g <- tempfile()
xmseg_full_h5 <- xcms:::.xcms_experiment_to_hdf5(a, h5f_full_g)
param <- PeakDensityParam(sampleGroups = sampleData(xmseg_full_h5)$sample_group,
                          minFraction = 0.4, bw = 30)
xmseg_full_h5 <- groupChromPeaks(xmseg_full_h5, param, msLevel = 1L)

test_that("XcmsExperimentHdf5 validation works", {
    a <- new("XcmsExperimentHdf5")
    expect_true(validObject(a))
    expect_equal(a@sample_id, character())
    expect_equal(a@chrom_peaks_ms_level, integer())
    expect_equal(a@features_ms_level, integer())

    expect_true(validObject(xmse_h5))
    a <- xmse_h5
    a@sample_id <- a@sample_id[c(1L, 3L)]
    expect_error(validObject(a), "number of samples does not match")
})

test_that("[,XcmsExperimentHdf5 works", {
    a <- new("XcmsExperiment")
    expect_error(a[1, 2], "subsetting by j")
})

test_that("chromPeaks,XcmsExperiementHdf5 works", {
    a <- new("XcmsExperimentHdf5")
    res <- chromPeaks(a)
    expect_equal(res, a@chromPeaks)
    expect_error(chromPeaks(a, isFilledColumn = TRUE), "not supported")

    a <- xmse_h5
    res <- chromPeaks(a, msLevel = c(1, 3))
    expect_equal(res, a@chromPeaks)
    res <- chromPeaks(a)
    ref <- chromPeaks(loadXcmsData("faahko_sub2"))
    expect_equal(colnames(res), colnames(ref))
    expect_equal(unname(res), unname(ref))

    res <- chromPeaks(a, msLevel = 1, columns = c("mz", "mzmin", "mzmax"))
    expect_equal(colnames(res), c("mz", "mzmin", "mzmax", "sample"))
    expect_equal(unname(res), unname(ref[, c("mz","mzmin","mzmax","sample")]))

    ## providing mz and rt
    res <- chromPeaks(a, msLevel = 1, type = "apex_within", rt = c(2500, 2600))
    expect_true(all(res[, "rt"] > 2500))
    expect_true(all(res[, "rt"] < 2600))
})

test_that("findChromPeaks,XcmsExperimentHdf5 works", {
    a <- as(xmse_h5, "MsExperiment")
    a <- as(a, "XcmsExperimentHdf5")
    h5_file <- tempfile()
    xcms:::.h5_initialize_file(h5_file)
    a@hdf5_file <- h5_file
    a@sample_id <- c("S1", "S2", "S3")
    p <- xmse_h5@processHistory[[1L]]@param
    a <- findChromPeaks(a, param = p, msLevel = 1L)
    expect_true(hasChromPeaks(a))
    res <- chromPeaks(a)
    ref <- chromPeaks(xmse_h5)
    expect_equal(res, ref)

    ## Errors/warnings
    expect_error(findChromPeaks(a, param = p, msLevel = 1:2), "single MS level")
    expect_warning(findChromPeaks(a, param = p, msLevel = 2), "No spectra of")

    ## Add to existing.
    a <- findChromPeaks(a, param = p, msLevel = 1L, add = TRUE)
    res <- chromPeaks(a)
    expect_equal(nrow(res), 2 * nrow(ref))
    expect_equal(anyDuplicated(rownames(res)), 0L)
    expect_equal(a@chrom_peaks_ms_level, 1L)

    ## Remove previous detection.
    a <- findChromPeaks(a, param = p, msLevel = 1L, add = FALSE)
    expect_true(hasChromPeaks(a))
    expect_equal(a@chrom_peaks_ms_level, 1L)
    expect_equal(chromPeaks(a), ref)

    unlink(h5_file)

    ## Test that it works with object being a MsExperiment
    a <- as(xmse_h5, "MsExperiment")
    a <- findChromPeaks(a, param = p, chunkSize = 3L, hdf5File = h5_file)
    expect_s4_class(a, "XcmsExperimentHdf5")
    expect_true(hasChromPeaks(a))
    expect_equal(chromPeaks(a), ref)
    a <- as(a, "MsExperiment")
    expect_error(findChromPeaks(a, param = p, hdf5File = h5_file),
                 "already exists")
    unlink(h5_file)
})

test_that("dropChromPeaks,XcmsExperimentHdf5 works", {
    res <- dropChromPeaks(xmse_h5)
    expect_false(hasChromPeaks(res))
    expect_equal(res@chrom_peaks_ms_level, integer())
    expect_true(validObject(res))
    ## With adjusted retention times
    ## With features
})

test_that("refineChromPeaks,XcmsExperimentHdf5,MergeNeighboringPeaksParam", {
    a <- new("XcmsExperimentHdf5")

    expect_warning(res <- refineChromPeaks(a, MergeNeighboringPeaksParam()),
                   "No chromatographic")
    expect_equal(a, res)

    af <- tempfile()
    ref <- loadXcmsData("faahko_sub2")
    a <- .xcms_experiment_to_hdf5(ref, af)
    res <- refineChromPeaks(a, MergeNeighboringPeaksParam())
    expect_error(validObject(a))
    expect_true(validObject(res))
    ## Compare results from both. Need chromPeaks() function first.
    ref <- refineChromPeaks(ref, MergeNeighboringPeaksParam())
    ref_pks <- chromPeaks(ref)
    res_pks <- .h5_read_data(res@hdf5_file, index = res@sample_id,
                             ms_level = rep(1L, length(res)),
                             read_colnames = TRUE, read_rownames = TRUE)
    res_pks <- do.call(
        rbind, mapply(FUN = function(x, i) cbind(x, sample = rep(i, nrow(x))),
                      res_pks, seq_along(res_pks)))
    expect_equal(dim(res_pks), dim(ref_pks))
    expect_equal(unname(res_pks), unname(ref_pks))

    file.remove(af)
})

test_that("groupChromPeaks,featureDefinitions,XcmsExperimentHdf5 works", {
    x <- xmse_full_h5
    param <- PeakDensityParam(sampleGroups = sampleData(x)$sample_group,
                              minFraction = 0.4, bw = 30)
    expect_error(groupChromPeaks(x, param, msLevel = 1:2), "one MS level")
    x2 <- x
    x2@chrom_peaks_ms_level <- integer()
    expect_error(groupChromPeaks(x2, param),
                 "No chromatographic")
    expect_error(groupChromPeaks(x, param, msLevel = 2L),
                 "No chromatographic")
    x <- groupChromPeaks(x, param)
    expect_true(validObject(x))
    expect_true(hasFeatures(x))
    expect_true(hasFeatures(x, 1L))
    expect_false(hasFeatures(x, 2L))
    expect_false(hasFeatures(x, 1:2))
    a <- .h5_read_data_frame("/features/ms_1/feature_definitions",
                             x@hdf5_file, read_rownames = TRUE)
    ref <- featureDefinitions(loadXcmsData("xmse"))
    ref$peakidx <- NULL
    ref$ms_level <- NULL
    rownames(a) <- NULL
    rownames(ref) <- NULL
    expect_true(all(colnames(ref) %in% colnames(a)))
    expect_equal(ref, a[, colnames(ref)])
    pks <- .h5_chrom_peaks(x, msLevel = 1L)
    for (i in seq_along(pks)) {
        b <- .h5_read_matrix(paste0("/S", i, "/ms_1/feature_to_chrom_peaks"),
                             x@hdf5_file)
        expect_true(all(b[, 2L] <= nrow(pks[[i]])))
    }
    expect_error(groupChromPeaks(x, param, msLevel = 1L, add = TRUE),
                 "currently not supported")

    ## featureDefinitions
    res <- featureDefinitions(x, msLevel = 2L)
    expect_true(is.data.frame(res))
    expect_true(nrow(res) == 0)
    res <- featureDefinitions(x, msLevel = 1:2)
    expect_true(is.data.frame(res))
    expect_true(nrow(res) > 0)
    ref$ms_level <- 1L
    rownames(res) <- NULL
    expect_true(all(colnames(res) %in% colnames(ref)))
    expect_equal(ref, res[, colnames(ref)])
})

test_that("hasFeatures,XcmsExperimentHdf5 works", {
    expect_false(hasFeatures(xmse_h5))
    expect_false(hasFeatures(new("XcmsExperimentHdf5")))
    expect_false(hasFeatures(new("XcmsExperimentHdf5"), msLevel = 2))
})

test_that("featureDefinitions,XcmsExperimentHdf5 works", {
    expect_error(featureDefinitions(xmse_full_h5) <- 4, "Not implemented")
})

test_that("featureValues,XcmsExperimentHdf5 etc works", {
    ref <- loadXcmsData("xmse")
    a <- featureDefinitions(ref)
    a$peakidx <- NULL
    b <- featureDefinitions(xmseg_full_h5)
    rownames(a) <- NULL
    rownames(b) <- NULL
    all(colnames(a) %in% colnames(b))
    expect_equal(a, b[, colnames(a)])
    nf <- nrow(b)
    rtmed <- b$rtmed
    ## .h5_feature_values_sample
    a <- .h5_feature_values_sample(
        xmseg_full_h5@hdf5_file, sample_id = "S1", ms_level = 1L,
        n_features = nf, method = "sum", filled = FALSE, col_idx = 9L)
    b <- unname(featureValues(ref, method = "sum", value = "maxo",
                              filled = FALSE)[, 1L])
    expect_equal(a, b)
    a <- .h5_feature_values_sample(
        xmseg_full_h5@hdf5_file, sample_id = "S4", ms_level = 1L,
        n_features = nf, filled = FALSE, method = "maxint", col_idx = c(7L, 9L))
    b <- unname(featureValues(ref, method = "maxint", value = "into",
                              filled = FALSE, intensity = "maxo")[, 4L])
    expect_equal(a, b)
    a <- .h5_feature_values_sample(
        xmseg_full_h5@hdf5_file, sample_id = "S4", ms_level = 1L,
        n_features = nf, filled = FALSE, method = "medret", col_idx = c(8L, 4L),
        rtmed = rtmed)
    b <- unname(featureValues(ref, method = "medret", value = "intb",
                              filled = FALSE)[, 4L])
    expect_equal(a, b)

    ## .h5_feature_values_ms_level
    a <- .h5_feature_values_ms_level(1L, xmseg_full_h5, method = "medret",
                                     value = "into", filled = FALSE)
    b <- featureValues(ref, method = "medret", value = "into", filled = FALSE)
    expect_equal(unname(a), unname(b))
    a <- .h5_feature_values_ms_level(1L, xmseg_full_h5, method = "sum",
                                     value = "maxo", filled = FALSE)
    b <- featureValues(ref, method = "sum", value = "maxo", filled = FALSE)
    expect_equal(unname(a), unname(b))
    a <- .h5_feature_values_ms_level(1L, xmseg_full_h5, method = "maxint",
                                     value = "sn", intensity = "into",
                                     filled = FALSE)
    b <- featureValues(ref, method = "maxint", value = "sn", intensity = "into",
                       filled = FALSE)
    expect_equal(unname(a), unname(b))

    ## featureValues
    expect_error(featureValues(xmse_h5), "No feature definitions")
    expect_error(featureValues(xmseg_full_h5, value = "index"),
                 "does not support")
    expect_error(featureValues(xmseg_full_h5, missing = "other"),
                 "or a numeric")
    ## column that does not exist.
    expect_error(featureValues(xmseg_full_h5, msLevel = 1L, value = "other"),
                 "Not all requested columns available.")
    ## check column names, missing values.
    fv_ref <- featureValues(ref, value = "into", method = "maxint",
                            intensity = "maxo", filled = FALSE)
    res <- featureValues(xmseg_full_h5, value = "into", method = "maxint",
                         intensity = "maxo", filled = FALSE)
    rownames(fv_ref) <- NULL
    rownames(res) <- NULL
    expect_equal(res, fv_ref)
    fv_ref <- featureValues(ref, value = "into", method = "maxint",
                            intensity = "maxo", filled = FALSE,
                            missing = "rowmin_half")
    res <- featureValues(xmseg_full_h5, value = "into", method = "maxint",
                         intensity = "maxo", filled = FALSE,
                         missing = "rowmin_half")
    rownames(fv_ref) <- NULL
    rownames(res) <- NULL
    expect_equal(res, fv_ref)
})

unlink(h5f)
unlink(h5f_full)
unlink(h5f_full_g)