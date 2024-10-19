h5f <- tempfile()
xmse_h5 <- .xcms_experiment_to_hdf5(loadXcmsData("faahko_sub2"), h5f)

test_that("XcmsExperimentHdf5 validation works", {
    a <- new("XcmsExperimentHdf5")
    expect_true(validObject(a))
    expect_equal(a@sample_id, character())
    expect_false(a@has_chrom_peaks)
    expect_false(a@has_features)

    expect_true(validObject(xmse_h5))
    a <- xmse_h5
    a@sample_id <- a@sample_id[c(1L, 3L)]
    expect_error(validObject(a), "number of samples does not match")
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
