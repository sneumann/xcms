library(MsExperiment)
mse <- MsExperiment()
fls <- normalizePath(faahko_3_files)
df <- data.frame(mzML_file = basename(fls),
                 dataOrigin = fls,
                 sample = c("ko15", "ko16", "ko18"))
sampleData(mse) <- DataFrame(df)
spectra(mse) <- Spectra::Spectra(fls)

## Link samples to spectra.
mse <- linkSampleData(mse, with = "sampleData.dataOrigin = spectra.dataOrigin")
p <- CentWaveParam(noise = 10000, snthresh = 40, prefilter = c(3, 10000))
xmse <- findChromPeaks(mse, param = p)

test_that(".empty_chrom_peaks works", {
    res <- .empty_chrom_peaks()
    expect_true(nrow(res) == 0)
    expect_equal(colnames(res), .REQ_PEAKS_COLS)

    res <- .empty_chrom_peaks(sample = FALSE)
    expect_true(nrow(res) == 0)
    expect_true(!any(colnames(res) == "sample"))
})

test_that("XcmsExperiment validation works", {
    a <- new("XcmsExperiment")
    expect_true(validObject(a))

    a@chromPeakData <- data.frame(a = 1:3, b = 1:3)
    expect_error(validObject(a), "chromPeakData")

    a <- new("XcmsExperiment")
    a@chromPeaks <- cbind(a = 1, b = 2)
    expect_error(validObject(a), "chromPeaks")
})

## All chrom peak related functions:
## - findChromPeaks
## - hasChromPeaks
## - dropChromPeaks
## - chromPeaks
## - chromPeakData
test_that("findChromPeaks,MsExperiment et al works", {
    expect_error(findChromPeaks(MsExperiment(), param = p), "No spectra")

    a <- MsExperiment()
    spectra(a) <- spectra(mse)
    expect_error(findChromPeaks(a, param = p), "No link")

    res <- xmse
    expect_equal(res@chromPeaks, chromPeaks(faahko_xod))
    expect_equal(res@chromPeakData, as.data.frame(chromPeakData(faahko_xod)))
    expect_true(hasChromPeaks(res))

    ## chromPeaks
    expect_equal(chromPeaks(res), res@chromPeaks)
    cp <- chromPeaks(res, isFilledColumn = TRUE)
    expect_true(any(colnames(cp) == "is_filled"))
    cp <- chromPeaks(res, msLevel = 2)
    expect_true(is.matrix(cp))
    expect_true(nrow(cp) == 0)
    cp <- chromPeaks(res, rt = c(3000, 3500), type = "within")
    expect_true(all(cp[, "rt"] >= 3000 & cp[, "rt"] <= 3500))
    cp <- chromPeaks(res, mz = c(300, 310), type = "within")
    expect_true(all(cp[, "mz"] >= 300 & cp[, "mz"] <= 310))

    ## chromPeakData
    expect_equal(chromPeakData(res), DataFrame(res@chromPeakData))

    ## dropChromPeaks
    rres <- dropChromPeaks(res)
    expect_true(length(rres@processHistory) == 0)
    expect_true(nrow(rres@chromPeaks) == 0)
    expect_false(hasChromPeaks(rres))

    res2 <- findChromPeaks(mse, param = p, msLevel = 2L)
    expect_true(nrow(res2@chromPeaks) == 0)
    expect_false(hasChromPeaks(res2))

    res2 <- findChromPeaks(res, param = p, msLevel = 2L, add = TRUE)
    expect_equal(res@chromPeaks, res2@chromPeaks)
    expect_equal(res@chromPeakData, res2@chromPeakData)
    expect_true(length(res2@processHistory) == 2)

    res2 <- findChromPeaks(res, param = p, msLevel = 2L, add = FALSE)
    expect_equal(nrow(res2@chromPeaks), 0)
    expect_true(length(res2@processHistory) == 1)

    res2 <- findChromPeaks(mse, param = p, chunkSize = -1)
    expect_equal(res@chromPeaks, res2@chromPeaks)

    expect_true(hasChromPeaks(res))
    expect_true(hasChromPeaks(res, msLevel = 1L))
    expect_true(hasChromPeaks(res, msLevel = 1:4))
    expect_false(hasChromPeaks(res, msLevel = 2))
})

test_that("subsetting,XcmsExperiment works", {
    expect_error(.subset_xcms_experiment(xmse, i = 1:4), "out of bounds")
    expect_error(.subset_xcms_experiment(xmse, i = c(1, 1, 2)), "Duplicated")

    res <- .subset_xcms_experiment(xmse, i = 2)
    expect_true(hasChromPeaks(res))
    expect_true(all(chromPeaks(res)[, "sample"] == 1L))
    expect_true(length(res) == 1)
    expect_equal(spectra(res), spectra(mse[2L]))
    cp <- chromPeaks(xmse)
    expect_equal(chromPeaks(res)[, colnames(chromPeaks(res)) != "sample"],
                 cp[cp[, "sample"] == 2L, colnames(cp) != "sample"])

    res <- .subset_xcms_experiment(xmse, i = c(3, 1))
    expect_true(hasChromPeaks(res))
    cp3 <- chromPeaks(res)[chromPeaks(res)[, "sample"] == 1, ]
    cp1 <- chromPeaks(res)[chromPeaks(res)[, "sample"] == 2, ]
    expect_equal(cp3[, colnames(cp3) != "sample"],
                 cp[cp[, "sample"] == 3, colnames(cp) != "sample"])
    expect_equal(cp1[, colnames(cp1) != "sample"],
                 cp[cp[, "sample"] == 1, colnames(cp) != "sample"])

    res <- .subset_xcms_experiment(xmse, i = 3, keepChromPeaks = FALSE)
    expect_false(hasChromPeaks(res))
    expect_true(length(res@processHistory) == 0)

    res <- .subset_xcms_experiment(xmse, i = 3, keepChromPeaks = FALSE,
                                   ignoreHistory = TRUE)
    expect_false(hasChromPeaks(res))
    expect_true(length(res@processHistory) == 1)

    expect_error(xmse[3, 4], "not supported")
    res <- xmse[c(3, 1)]
    expect_true(hasChromPeaks(res))
    cp3 <- chromPeaks(res)[chromPeaks(res)[, "sample"] == 1, ]
    cp1 <- chromPeaks(res)[chromPeaks(res)[, "sample"] == 2, ]
    expect_equal(cp3[, colnames(cp3) != "sample"],
                 cp[cp[, "sample"] == 3, colnames(cp) != "sample"])
    expect_equal(cp1[, colnames(cp1) != "sample"],
                 cp[cp[, "sample"] == 1, colnames(cp) != "sample"])
})

test_that("filterRt,XcmsExperiment works", {
    res <- filterRt(xmse)
    expect_equal(res, xmse)

    res <- filterRt(xmse, rt = c(3000, 3500))
    expect_true(all(rtime(spectra(res)) >= 3000 &
                    rtime(spectra(res)) <= 3500))
    expect_true(all(chromPeaks(res)[, "rt"] >= 3000 &
                    chromPeaks(res)[, "rt"] <= 3500))
    ## Check error: define msLevel
    expect_warning(res_2 <- filterRt(xmse, rt = c(3000, 3500), msLevel = 2L),
                   "ignored")
    expect_equal(chromPeaks(res), chromPeaks(res_2))
    expect_equal(rtime(spectra(res)), rtime(spectra(res_2)))
})

test_that("filterFile,XcmsExperiment works", {
    res <- filterFile(xmse)
    expect_s4_class(res, "XcmsExperiment")
    expect_true(length(res) == 0)
    expect_false(hasChromPeaks(res))

    res <- filterFile(xmse, 2)
    expect_equal(res, xmse[2])

    res <- filterFile(xmse, c(3, 1))
    expect_equal(res, xmse[c(1, 3)])
})

test_that("adjustRtime,MsExperiment,XcmsExperiment,ObiwarpParam works", {
    op <- ObiwarpParam(binSize = 35.5)
    ref <- adjustRtime(faahko_xod, param = op)

    res <- adjustRtime(mse, param = op)
    expect_equal(unname(rtime(ref)), spectra(res)$rtime_adjusted)
    expect_s4_class(res, "XcmsExperiment")
    expect_true(length(res@processHistory) == 1L)
    expect_true(hasAdjustedRtime(res))
    expect_equal(rtime(res), unname(rtime(ref)))
    expect_equal(rtime(res, adjusted = FALSE),
                 unname(rtime(ref, adjusted = FALSE)))

    ## applyAdjustedRtime
    res2 <- applyAdjustedRtime(res)
    expect_false(hasAdjustedRtime(res2))
    expect_equal(rtime(res2), rtime(res, adjusted = TRUE))
    res2 <- dropAdjustedRtime(res)
    expect_false(hasAdjustedRtime(res2))
    expect_true(length(res2@processHistory) == 0L)

    ## xcms object.
    res2 <- adjustRtime(xmse, param = op)
    expect_true(hasAdjustedRtime(res2))
    expect_equal(rtime(res, adjusted = FALSE), rtime(res2, adjusted = FALSE))
    expect_equal(rtime(res, adjusted = TRUE), rtime(res2, adjusted = TRUE))
    expect_equal(chromPeaks(ref), chromPeaks(res2))
    ## chrom peaks got adjusted too
    a <- chromPeaks(xmse)
    b <- chromPeaks(res2)
    expect_true(all(a[a[, "sample"] == 1L, "rt"] !=
                    b[b[, "sample"] == 1L, "rt"]))
    ## those of center sample are not changed
    expect_true(all(a[a[, "sample"] == 2L, "rt"] ==
                    b[b[, "sample"] == 2L, "rt"]))
    expect_true(length(res2@processHistory) == 2L)

    ## Order: peak detection, alignment.
    ## dropAdjustedRtime:
    res3 <- dropAdjustedRtime(res2)
    expect_false(hasAdjustedRtime(res3))
    expect_true(hasChromPeaks(res3))
    ## chrom peak rt gets reverted
    expect_equal(chromPeaks(res3), chromPeaks(xmse))
    expect_true(length(res3@processHistory) == 1L)
    ref2 <- dropAdjustedRtime(ref)
    expect_equal(chromPeaks(res3), chromPeaks(ref2))
    ## dropChromPeaks
    res3 <- dropChromPeaks(res2)
    expect_false(hasChromPeaks(res3))
    expect_false(hasAdjustedRtime(res3))
    expect_true(length(res3@processHistory) == 0L)
    res3 <- dropChromPeaks(res2, keepAdjustedRtime = TRUE)
    expect_false(hasChromPeaks(res3))
    expect_true(hasAdjustedRtime(res3))
    expect_true(length(res3@processHistory) == 1L)
    expect_equal(rtime(res3, adjusted = TRUE), rtime(res2, adjusted = TRUE))

    ## Order: alignment, peak detection.
    res3 <- findChromPeaks(res, param = p)
    expect_true(hasChromPeaks(res3))
    expect_true(hasAdjustedRtime(res3))
    expect_true(length(res3@processHistory) == 2L)
    ## dropAdjustedRtime
    res4 <- dropAdjustedRtime(res3)
    expect_true(hasChromPeaks(res4))
    expect_false(hasAdjustedRtime(res4))
    expect_true(length(res4@processHistory) == 1L)
    ## chrom peak rt should be "reverted" (as with raw data) - but they
    ## are not identical because of the interpolation
    ## expect_equal(chromPeaks(res4), chromPeaks(xmse))
    ## a <- chromPeaks(res4)
    ## b <- chromPeaks(xmse)
    ## expect_true(all(a[a[, "sample"] == 1L, "rt"] ==
    ##                 b[b[, "sample"] == 1L, "rt"]))
    ## expect_true(all(a[a[, "sample"] == 2L, "rt"] ==
    ##                 b[b[, "sample"] == 2L, "rt"]))
    ## expect_true(all(a[a[, "sample"] == 3L, "rt"] ==
    ##                 b[b[, "sample"] == 3L, "rt"]))
    ## dropChromPeaks
    res4 <- dropChromPeaks(res3)
    expect_false(hasChromPeaks(res4))
    expect_true(hasAdjustedRtime(res4))
    expect_true(length(res4@processHistory) == 1L)

    ## With spectra that are NOT all associated to a sample.
    mse2 <- MsExperiment()
    sampleData(mse2) <- DataFrame(df)

    sps <- spectra(mse)
    tmp <- sps[1:10]
    tmp$dataOrigin <- "a"
    spectra(mse2) <- c(tmp, sps)
    mse2 <- linkSampleData(
        mse2, with = "sampleData.dataOrigin = spectra.dataOrigin")
    expect_error(adjustRtime(mse2, param = op), "to a sample")
})
