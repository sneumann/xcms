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

test_that(".mse_same_rownames works", {
    expect_equal(.mse_same_rownames(chromPeaks(faahko_xod),
                                    chromPeakData(faahko_xod)), NULL)
    a <- cbind(a = 1:3, b = 1:3)
    rownames(a) <- c("a", "b", "c")
    expect_match(.mse_same_rownames(a, data.frame(ms_level = 1:3)),
                 "don't match")
})
