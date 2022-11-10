test_that(".mse_add_chrom_peaks works", {
    a <- new("XcmsExperiment")
    pks <- chromPeaks(faahko_xod)
    pkd <- as.data.frame(chromPeakData(faahko_xod))

    res <- .mse_add_chrom_peaks(a, pks, pkd)
    expect_equal(res@chromPeakData, pkd)
    expect_equal(res@chromPeaks, pks)

    pkd$add_col <- "a"
    res <- .mse_add_chrom_peaks(res, pks, pkd)
    expect_equal(rownames(res@chromPeaks), rownames(res@chromPeakData))
    expect_equal(rownames(res@chromPeaks),
                 .featureIDs(nrow(res@chromPeaks), prefix = "CP"))
    expect_equal(colnames(res@chromPeakData),
                 c("ms_level", "is_filled", "add_col"))
})

test_that(".mse_valid_chrom_peaks works", {
    a <- cbind(mz = 1:3, mzmin = 1:3)
    expect_match(.mse_valid_chrom_peaks(a), "missing")
    expect_equal(.mse_valid_chrom_peaks(chromPeaks(faahko_xod)), NULL)
})

test_that(".mse_valid_chrom_peak_data works", {
    expect_equal(.mse_valid_chrom_peak_data(chromPeakData(faahko_xod)), NULL)
    a <- data.frame(msLevel = 1:3, other = "a")
    expect_match(.mse_valid_chrom_peak_data(a), "missing")
})

test_that(".mse_valid_feature_def", {
    expect_equal(.mse_valid_feature_def(.empty_feature_definitions()), NULL)
    expect_match(.mse_valid_feature_def(3), "required columns")
    expect_match(.mse_valid_feature_def(data.frame(a = 1:4, b = 1:4)),
                 "required columns")
})

test_that(".mse_same_rownames works", {
    expect_equal(.mse_same_rownames(chromPeaks(faahko_xod),
                                    chromPeakData(faahko_xod)), NULL)
    a <- cbind(a = 1:3, b = 1:3)
    rownames(a) <- c("a", "b", "c")
    expect_match(.mse_same_rownames(a, data.frame(ms_level = 1:3)),
                 "don't match")
})

test_that(".aggregate_intensities works", {
    pks <- cbind(mz = c(12.3, 12.4, 12.5, 13, 13.2, 13.3, 13.4, 13.5),
                 intensity = 1:8)
    res <- .aggregate_intensities(pks)
    expect_equal(res, sum(1:8))

    res <- .aggregate_intensities(pks, mzr = c(12.4, 13.2))
    expect_equal(res, sum(2:5))

    res <- .aggregate_intensities(pks, mzr = c(12.4, 13.2), INTFUN = max)
    expect_equal(res, 5)
})

test_that("sumi works", {
    expect_equal(sumi(1:4), sum(1:4))
    expect_equal(sumi(c(1:4, NA)), sum(c(1:4, NA), na.rm = TRUE))
    expect_equal(sumi(c(NA, NA)), NA_real_)
})
