h5f <- tempfile()
xmse_h5 <- .xcms_experiment_to_hdf5(loadXcmsData("faahko_sub2"), h5f)

test_that("XcmsExperimentHdf5 validation works", {
    a <- new("XcmsExperimentHdf5")
    expect_true(validObject(a))
    expect_equal(a@sample_id, integer())
    expect_false(a@has_chrom_peaks)
    expect_false(a@has_features)

    expect_true(validObject(xmse_h5))
    a <- xmse_h5
    a@sample_id <- c(1L, 3L)
    expect_error(validObject(a), "number of samples does not match")
})
