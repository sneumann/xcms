h5f <- tempfile()
xmse_h5 <- xcms:::.xcms_experiment_to_hdf5(loadXcmsData("faahko_sub2"), h5f)
h5f_full <- tempfile()
a <- loadXcmsData("xmse") |>
    dropFeatureDefinitions()
xmse_full_h5 <- xcms:::.xcms_experiment_to_hdf5(a, h5f_full)

## correspondence
h5f_full_g <- tempfile()
xmseg_full_h5 <- xcms:::.xcms_experiment_to_hdf5(a, h5f_full_g)
pdp <- PeakDensityParam(sampleGroups = sampleData(xmseg_full_h5)$sample_group,
                        minFraction = 0.4, bw = 30)
xmseg_full_h5 <- groupChromPeaks(xmseg_full_h5, pdp, msLevel = 1L)
## reference
xmseg_full_ref <- dropFeatureDefinitions(loadXcmsData("xmse"))
xmseg_full_ref <- groupChromPeaks(xmseg_full_ref, pdp, msLevel = 1L)

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

test_that("chromPeakData,XcmsExperimentHdf5 works", {
    a <- new("XcmsExperimentHdf5")
    res <- chromPeakData(a)
    expect_s4_class(res, "DataFrame")
    expect_true(nrow(res) == 0)
    res <- chromPeakData(a, return.type = "data.frame")
    expect_true(is.data.frame(res))
    expect_true(nrow(res) == 0)

    a <- xmse_h5
    res <- chromPeakData(a, msLevel = 3L)
    expect_s4_class(res, "DataFrame")
    expect_true(nrow(res) == 0)

    cp <- chromPeaks(a)
    res <- chromPeakData(a, return.type = "data.frame")
    expect_true(is.data.frame(res))
    expect_equal(colnames(res), c("is_filled", "ms_level"))
    expect_equal(nrow(res), nrow(cp))
    expect_equal(rownames(res), rownames(cp))

    res <- chromPeakData(a, peaks = rownames(cp)[3:10])
    expect_s4_class(res, "DataFrame")
    expect_true(nrow(res) == 8)
    expect_equal(rownames(res), rownames(cp)[3:10])

    res <- chromPeakData(a, msLevel = c(1L, 3L))
    expect_true(nrow(res) == 0)
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
    a <- xcms:::.xcms_experiment_to_hdf5(ref, af)
    res <- refineChromPeaks(a, MergeNeighboringPeaksParam())
    expect_error(validObject(a))
    expect_true(validObject(res))
    ## Compare results from both. Need chromPeaks() function first.
    ref <- refineChromPeaks(ref, MergeNeighboringPeaksParam())
    ref_pks <- chromPeaks(ref)
    res_pks <- xcms:::.h5_read_data(res@hdf5_file, id = res@sample_id,
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
    a <- xcms:::.h5_read_data_frame("/features/ms_1/feature_definitions",
                             x@hdf5_file, read_rownames = TRUE)
    ref <- featureDefinitions(loadXcmsData("xmse"))
    ref$peakidx <- NULL
    ref$ms_level <- NULL
    rownames(a) <- NULL
    rownames(ref) <- NULL
    expect_true(all(colnames(ref) %in% colnames(a)))
    expect_equal(ref, a[, colnames(ref)])
    pks <- xcms:::.h5_chrom_peaks(x, msLevel = 1L)
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
    a <- xcms:::.h5_feature_values_sample(
        xmseg_full_h5@hdf5_file, sample_id = "S1", ms_level = 1L,
        n_features = nf, method = "sum", filled = FALSE, col_idx = 9L)
    b <- unname(featureValues(ref, method = "sum", value = "maxo",
                              filled = FALSE)[, 1L])
    expect_equal(a, b)
    a <- xcms:::.h5_feature_values_sample(
        xmseg_full_h5@hdf5_file, sample_id = "S4", ms_level = 1L,
        n_features = nf, filled = FALSE, method = "maxint", col_idx = c(7L, 9L))
    b <- unname(featureValues(ref, method = "maxint", value = "into",
                              filled = FALSE, intensity = "maxo")[, 4L])
    expect_equal(a, b)
    a <- xcms:::.h5_feature_values_sample(
        xmseg_full_h5@hdf5_file, sample_id = "S4", ms_level = 1L,
        n_features = nf, filled = FALSE, method = "medret", col_idx = c(8L, 4L),
        rtmed = rtmed)
    b <- unname(featureValues(ref, method = "medret", value = "intb",
                              filled = FALSE)[, 4L])
    expect_equal(a, b)

    ## .h5_feature_values_ms_level
    a <- xcms:::.h5_feature_values_ms_level(1L, xmseg_full_h5, method = "medret",
                                     value = "into", filled = FALSE)
    b <- featureValues(ref, method = "medret", value = "into", filled = FALSE)
    expect_equal(unname(a), unname(b))
    a <- xcms:::.h5_feature_values_ms_level(1L, xmseg_full_h5, method = "sum",
                                     value = "maxo", filled = FALSE)
    b <- featureValues(ref, method = "sum", value = "maxo", filled = FALSE)
    expect_equal(unname(a), unname(b))
    a <- xcms:::.h5_feature_values_ms_level(1L, xmseg_full_h5, method = "maxint",
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

test_that("adjustRtimePeakGroups works", {
    ref <- xmseg_full_ref

    a <- featureValues(
        ref, method = "maxint", intensity = "into", value = "rt")
    b <- featureValues(
        xmseg_full_h5, method = "maxint", intensity = "into", value = "rt")
    expect_equal(unname(a), unname(b))

    p <- PeakGroupsParam(minFraction = 0.7, extraPeaks = 100,
                         subset = c(1, 2, 20, 23))
    expect_error(adjustRtimePeakGroups(xmse_full_h5, PeakGroupsParam()),
                 "No features present")
    expect_error(adjustRtimePeakGroups(xmseg_full_h5, p, msLevel = 2L),
                 "No features present")
    expect_error(adjustRtimePeakGroups(xmseg_full_h5, p),
                 "out of bounds")
    p@subset <- c(1L, 3L, 4L, 7L, 8L)

    apeaks_ref <- adjustRtimePeakGroups(ref, p)
    apeaks <- adjustRtimePeakGroups(xmseg_full_h5, p)
    expect_equal(unname(apeaks_ref), unname(apeaks))
    expect_equal(colnames(apeaks_ref), colnames(apeaks))

    p <- PeakGroupsParam(minFraction = 0.3, extraPeaks = 100,
                         subset = c(1, 2, 3, 4, 7, 8))
    apeaks_ref <- adjustRtimePeakGroups(ref, p)
    apeaks <- adjustRtimePeakGroups(xmseg_full_h5, p)
    expect_equal(unname(apeaks_ref), unname(apeaks))
    expect_equal(colnames(apeaks_ref), colnames(apeaks))
})

test_that("adjustRtime,XcmsExperimentHdf5 and related function work", {
    ## Note: using `extraPeaks = 100` because that parameter is not supported
    ## for XcmsExperimentHdf5
    p <- PeakGroupsParam(span = 0.4, minFraction = 0.7, subset = c(1, 3, 5, 7),
                         extraPeaks = 100)
    ## Define the reference data
    ref <- loadXcmsData("xmse") |>
        dropFeatureDefinitions() |>
        applyAdjustedRtime()
    res_h5 <- tempfile()
    res <- xcms:::.xcms_experiment_to_hdf5(ref, res_h5)
    ## Create a single sample XcmsExperimentHdf5
    a <- ref[3L]
    a_h5 <- tempfile()
    a <- xcms:::.xcms_experiment_to_hdf5(a, a_h5)
    ## Perform retention time alignment on reference data
    ref <- ref |>
        groupChromPeaks(pdp, msLevel = 1L) |>
        adjustRtime(param = p)
    rt_raw <- rtime(ref, adjusted = FALSE)
    rt_raw <- split(rt_raw, fromFile(ref))[[3L]]
    rt_adj <- rtime(ref, adjusted = TRUE)
    rt_adj <- split(rt_adj, fromFile(ref))[[3L]]
    cp_raw <- chromPeaks(a)

    ############################################################################
    ## .h5_update_rt_chrom_peaks_sample: adjust rt of chrom peaks:
    cnt <- xcms:::.h5_update_rt_chrom_peaks_sample(
        a@sample_id[1L], rt_raw, rt_adj, 1L, a@hdf5_file)
    expect_equal(cnt, a@hdf5_mod_count + 1L)
    a@hdf5_mod_count <- cnt
    cp_adj <- chromPeaks(a)
    expect_true(all(cp_raw[, "rt"] != cp_adj[, "rt"]))
    expect_true(all(cp_raw[, "rtmin"] != cp_adj[, "rtmin"]))
    expect_true(all(cp_raw[, "rtmax"] != cp_adj[, "rtmax"]))

    ############################################################################
    ## adjustRtime: retention time adjustment
    ## errors
    expect_error(adjustRtime(xmseg_full_h5, p), "Alignment results already")
    expect_error(adjustRtime(res, p, msLevel = 2L), "supported for MS level 1")
    expect_error(adjustRtime(res, p), "No feature definitions present")

    ## Perform alignment
    res <- groupChromPeaks(res, pdp, msLevel = 1L)
    expect_false(hasAdjustedRtime(res))
    cp_ref_raw <- chromPeaks(res)
    res <- adjustRtime(res, param = p)
    expect_true(hasAdjustedRtime(res))
    expect_equal(rtime(res), rtime(ref))
    expect_equal(unname(chromPeaks(ref)), unname(chromPeaks(res)))
    expect_true(all(chromPeaks(res)[, "rt"] != cp_ref_raw[, "rt"]))
    expect_true(validObject(res))

    ############################################################################
    ## dropAdjustedRtime: revert retention times
    cp_ref_adj <- chromPeaks(res)
    cnt <- res@hdf5_mod_count
    phl <- length(res@processHistory)
    res <- dropAdjustedRtime(res)
    expect_true(res@hdf5_mod_count > cnt)
    expect_true(length(res@processHistory) < phl)
    expect_false(hasAdjustedRtime(res))
    expect_true(all(chromPeaks(res)[, "rt"] != cp_ref_adj[, "rt"]))
    ref <- dropAdjustedRtime(ref)
    expect_equal(rtime(ref), rtime(res))
    expect_equal(unname(chromPeaks(ref)), unname(chromPeaks(res)))
    res <- dropAdjustedRtime(res)
    expect_false(hasAdjustedRtime(res))

    unlink(a_h5)
    unlink(res_h5)
})

test_that(".hasFilledPeaks works with XcmsExperimentHdf5", {
    expect_false(.hasFilledPeaks(xmse_h5))
})

test_that("chromatogram,XcmsExperimentHdf5 works", {
    expect_error(chromatogram(xmse_h5, adjustedRtime = FALSE), "unused")
    expect_warning(res <- chromatogram(xmse_h5, include = "apex_within",
                                       return.type = "MChromatograms"),
                   "deprecated")
    expect_s4_class(res, "MChromatograms")
    expect_true(nrow(res) == 1L)
    ref <- chromatogram(faahko_od)
    expect_equal(intensity(res[1, 1]), unname(intensity(ref[1, 1])))

    rtr <- c(2600, 2700)
    mzr <- c(340, 400)
    res <- chromatogram(xmse_h5, mz = mzr, rt = rtr)
    expect_s4_class(res, "XChromatograms")
    expect_true(nrow(res) == 1L)
    expect_true(nrow(chromPeaks(res)) > 0)
    expect_true(all(chromPeaks(res)[, "mz"] >= 340 &
                    chromPeaks(res)[, "mz"] <= 400))
    expect_true(all(chromPeaks(res[1, 1])[, "sample"] == 1L))
    expect_true(all(chromPeaks(res[1, 2])[, "sample"] == 2L))
    expect_true(all(chromPeaks(res[1, 3])[, "sample"] == 3L))
    ref <- chromatogram(xod_x, mz = mzr, rt = rtr)
    expect_equal(unname(chromPeaks(res)), unname(chromPeaks(ref)))

    ## with features
    res <- chromatogram(
        xmseg_full_h5, mz = chromPeaks(xmseg_full_h5)[1:5, c("mzmin", "mzmax")],
        rt = chromPeaks(xmseg_full_h5)[1:5, c("rtmin", "rtmax")],
        chunkSize = 2L, BPPARAM = bpparam(), msLevel = 1L,
        aggregationFun = "sum", isolationWindow = NULL,
        chromPeaks = "apex_within", return.type = "XChromatograms")
    expect_true(nrow(featureDefinitions(res)) == 2)
    expect_true(all(unlist(featureDefinitions(res)$peakidx) %in%
                    seq_len(nrow(chromPeaks(res)))))
    ref <- chromatogram(
        xmseg_full_ref,mz = chromPeaks(xmseg_full_h5)[1:5, c("mzmin", "mzmax")],
        rt = chromPeaks(xmseg_full_h5)[1:5, c("rtmin", "rtmax")],
        chunkSize = 2L, BPPARAM = bpparam(), msLevel = 1L,
        aggregationFun = "sum", isolationWindow = NULL,
        chromPeaks = "apex_within", return.type = "XChromatograms")
    a <- featureDefinitions(res)
    b <- featureDefinitions(ref)
    rownames(a) <- NULL
    rownames(b) <- NULL
    expect_true(all(colnames(a) %in% colnames(b)))
    expect_equal(a[, colnames(b)], b[, colnames(b)])

    ## MS2 data.
    res <- chromatogram(
        xmseg_full_h5, msLevel = 2L,
        mz = chromPeaks(xmseg_full_h5)[1:5, c("mzmin", "mzmax")],
        rt = chromPeaks(xmseg_full_h5)[1:5, c("rtmin", "rtmax")])
    expect_true(validObject(res))
    expect_true(length(intensity(res[[1L]])) == 0)
    expect_true(length(intensity(res[[2L]])) == 0)
    expect_s4_class(res, "XChromatograms")
    expect_true(nrow(chromPeaks(res)) == 0)

    ## Defining only mz or rt.
    rtr <- c(2600, 2700)
    mzr <- c(340, 400)
    res <- chromatogram(xmse_h5, mz = mzr)
    expect_s4_class(res, "XChromatograms")
    expect_true(nrow(res) == 1L)
    expect_true(nrow(chromPeaks(res)) > 0)
    expect_true(all(chromPeaks(res)[, "mz"] >= 340 &
                    chromPeaks(res)[, "mz"] <= 400))
    expect_true(all(chromPeaks(res[1, 1])[, "sample"] == 1L))
    expect_true(all(chromPeaks(res[1, 2])[, "sample"] == 2L))
    expect_true(all(chromPeaks(res[1, 3])[, "sample"] == 3L))
    rrt <- range(lapply(res, rtime))
    expect_true(rrt[1] < 2600)
    expect_true(rrt[2] > 4400)

    res <- chromatogram(xmse_h5, rt = rtr)
    expect_s4_class(res, "XChromatograms")
    expect_true(nrow(res) == 1L)
    expect_true(nrow(chromPeaks(res)) > 0)
    expect_true(any(chromPeaks(res)[, "mz"] < 340 |
                    chromPeaks(res)[, "mz"] > 400))
    expect_true(all(chromPeaks(res[1, 1])[, "sample"] == 1L))
    expect_true(all(chromPeaks(res[1, 2])[, "sample"] == 2L))
    expect_true(all(chromPeaks(res[1, 3])[, "sample"] == 3L))
    rrt <- range(lapply(res, rtime))
    expect_true(rrt[1] >= 2600)
    expect_true(rrt[2] <= 2700)
})

test_that("hasFilledChromPeaks,XcmsExperimentHdf5 works", {
    expect_false(hasFilledChromPeaks(new("XcmsExperimentHdf5")))
})

test_that("fillChromPeaks,XcmsExperimentHdf5 works", {
    expect_error(fillChromPeaks(new("XcmsExperimentHdf5"),
                                ChromPeakAreaParam(),
                                msLevel = 1:2), "one MS level at a time")
    expect_error(fillChromPeaks(xmse_full_h5, ChromPeakAreaParam()),
                 "No feature definitions")
})

test_that("featureArea,XcmsExperimentHdf5 works", {
    expect_error(featureArea(xmse_h5), "No correspondence")
    expect_error(featureArea(xmseg_full_h5, msLevel = 2L), "No correspondence")
    expect_error(featureArea(xmseg_full_h5, features = c("a", "b")),
                 "Some of the provided")
    res <- featureArea(xmseg_full_h5)
    ref <- featureArea(xmseg_full_ref)
    expect_equal(unname(res), unname(ref))
    res <- featureArea(xmseg_full_h5, features = rownames(res)[c(5, 12, 20)])
    ref <- featureArea(xmseg_full_ref, features = rownames(ref)[c(5, 12, 20)])
    expect_equal(unname(res), unname(ref))
})

test_that("fillChromPeaks,XcmsExperimentHdf5,PeakAreaParam", {
    tf <- tempfile()
    file.copy(xmseg_full_h5@hdf5_file, tf)
    x <- xmseg_full_h5
    x@hdf5_file <- tf
    fvals <- featureValues(x, msLevel = 1L)
    cps <- chromPeaks(x, msLevel = 1L)

    expect_error(
        fillChromPeaks(x, param = ChromPeakAreaParam(), msLevel = 1:2),
        "Can only perform peak filling")
    expect_error(
        fillChromPeaks(x, param = ChromPeakAreaParam(), msLevel = 2),
        "No feature definitions for MS level")

    p <- ChromPeakAreaParam(mzmin = min, mzmax = max, rtmin = min, rtmax = max)
    res <- fillChromPeaks(x, param = p)
    expect_true(res@hdf5_mod_count > x@hdf5_mod_count)
    expect_equal(res@gap_peaks_ms_level, 1L)
    res_cpd <- chromPeakData(res)
    res_cps <- chromPeaks(res)
    expect_true(sum(res_cpd$is_filled) > 0)
    expect_true(length(res@processHistory) > length(x@processHistory))

    ## Compare results with "reference"
    ref <- fillChromPeaks(xmseg_full_ref, p)
    ref_cpd <- chromPeakData(ref)
    ref_cps <- chromPeaks(ref)
    idx <- order(ref_cps[, "sample"])
    ref_cpd <- ref_cpd[idx, ]
    ref_cps <- ref_cps[idx, ]
    expect_equal(res_cpd$is_filled, ref_cpd$is_filled)
    expect_equal(unname(res_cps), unname(ref_cps))

    ## Compare feature values.
    fvals_res <- featureValues(res, msLevel = 1L)
    expect_equal(dim(fvals_res), dim(fvals))
    expect_equal(dimnames(fvals_res), dimnames(fvals))
    expect_true(sum(is.na(fvals_res)) < sum(is.na(fvals)))
    rownames(fvals_res) <- NULL

    fvals_ref <- featureValues(ref, msLevel = 1L)
    rownames(fvals_ref) <- NULL
    expect_equal(fvals_ref, fvals_res)

    ## Test featureValues with filled = FALSE
    tmp <- featureValues(res, msLevel = 1L, filled = FALSE)
    expect_equal(tmp, fvals)

    ## dropFilledChromPeaks
    res <- dropFilledChromPeaks(res)
    expect_equal(res@hdf5_mod_count, 66L)
    expect_false(hasFilledChromPeaks(res))
    cps_res <- chromPeaks(res)
    expect_equal(cps_res, cps)
    cpd_res <- chromPeakData(res)
    expect_true(all(!cpd_res$is_filled))

    fvals_res <- featureValues(res, msLevel = 1L)
    expect_equal(fvals, fvals_res)
    expect_equal(res@processHistory, x@processHistory)

    rm(tf)
})

test_that("filterMsLevel,XcmsExperimentHdf5 works", {
    res <- filterMsLevel(xmseg_full_h5, msLevel. = integer())
    expect_equal(res@chrom_peaks_ms_level, xmseg_full_h5@chrom_peaks_ms_level)
    res <- filterMsLevel(xmseg_full_h5)
    expect_equal(res@chrom_peaks_ms_level, xmseg_full_h5@chrom_peaks_ms_level)
    expect_equal(res@chrom_peaks_ms_level, 1L)
    expect_equal(msLevel(res@spectra), msLevel(xmseg_full_h5@spectra))
    expect_s4_class(res, "XcmsExperimentHdf5")

    res <- filterMsLevel(xmseg_full_h5, 2L)
    expect_equal(res@chrom_peaks_ms_level, integer())
    expect_equal(res@gap_peaks_ms_level, integer())
    expect_equal(res@features_ms_level, integer())
    expect_equal(length(res@spectra), 0L)
    expect_equal(sampleData(res), sampleData(xmseg_full_h5))
    expect_s4_class(res, "XcmsExperimentHdf5")
})

test_that("filterRt,XcmsExperimentHdf5 works", {
    tf <- tempfile()
    file.copy(xmse_h5@hdf5_file, tf)
    x <- xmse_h5
    x@hdf5_file <- tf
    x <- filterRt(x, rt = c(3300, 3500))
    expect_true(validObject(x))
    ref <- loadXcmsData("faahko_sub2")
    ref <- filterRt(ref, rt = c(3300, 3500))
    expect_equal(unname(chromPeaks(x)), unname(chromPeaks(ref)))
    rm(tf)

    ## with features
    tf <- tempfile()
    file.copy(xmseg_full_h5@hdf5_file, tf)
    x <- xmseg_full_h5
    x@hdf5_file <- tf
    x <- filterRt(x, rt = c(3300, 3500))
    expect_true(validObject(x))
    ref <- filterRt(xmseg_full_ref, rt = c(3300, 3500))
    expect_equal(unname(chromPeaks(x)), unname(chromPeaks(ref)))
    a <- chromPeakData(x)
    b <- chromPeakData(ref)
    rownames(a) <- NULL
    rownames(b) <- NULL
    expect_equal(a, b[, colnames(a)])
    a <- featureDefinitions(x)
    b <- featureDefinitions(ref)
    rownames(a) <- NULL
    rownames(b) <- NULL
    expect_equal(a, b[, colnames(a)])
    a <- featureValues(x, method = "sum")
    b <- featureValues(ref, method = "sum")
    expect_equal(unname(a), unname(b))

    rm(ft)
})

## test_that(".h5_feature_chrom_peaks_sample works", {
##     cn <- .h5_chrom_peaks_colnames(xmseg_full_h5, 1L)
##     res <- .h5_feature_chrom_peaks_sample("S3", xmseg_full_h5@hdf5_file,
##                                           1L, j = match("into", cn))
##     ref <- featureValues(xmseg_full_h5, method = "sum", value = "into")
##     vals <- split(res[, 2L], factor(res[, 1L], levels = seq_len(nrow(ref))))
##     vals <- vapply(vals, function(z) {
##         if (length(z))
##             sum(z)
##         else NA_real_
##     }, 2.2)
##     expect_equal(unname(vals), unname(ref[, 3L]))
##     ## With index in arbitrary order and with duplicates
##     i <- c(1, 4, 2, 3, 2)
##     res <- .h5_feature_chrom_peaks_sample("S3", xmseg_full_h5@hdf5_file,
##                                           1L, j = match("into", cn), i = i)
##     expect_equal(res[, 1L], c(4, 2, 2))
##     expect_equal(res[, 2L], unname(ref[c(4, 2, 2), 3L]))
## })

unlink(h5f)
unlink(h5f_full)
unlink(h5f_full_g)
