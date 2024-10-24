h5f <- tempfile()
xmse_h5 <- xcms:::.xcms_experiment_to_hdf5(loadXcmsData("faahko_sub2"), h5f)
h5f_full <- tempfile()
a <- dropFeatureDefinitions(loadXcmsData("xmse"))
xmse_full_h5 <- xcms:::.xcms_experiment_to_hdf5(a, h5f_full)

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

test_that("groupChromPeaks,XcmsExperimentHdf5 works", {
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
    expect_equal(ref, a)
    pks <- .h5_chrom_peaks(x, msLevel = 1L)
    for (i in seq_along(pks)) {
        b <- .h5_read_matrix(paste0("/S", i, "/ms_1/feature_to_chrom_peaks"),
                             x@hdf5_file)
        expect_true(all(b[, 2L] <= nrow(pks[[i]])))
    }
    expect_error(groupChromPeaks(x, param, msLevel = 1L, add = TRUE),
                 "currently not supported")
})

test_that("hasFeatures,XcmsExperimentHdf5 works", {
    expect_false(hasFeatures(xmse_h5))
    expect_false(hasFeatures(new("XcmsExperimentHdf5")))
    expect_false(hasFeatures(new("XcmsExperimentHdf5"), msLevel = 2))
})
