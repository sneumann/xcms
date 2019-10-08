test_that("do_groupChromPeaks_density works", {
    fts <- peaks(faahko)
    grps <- sampclass(faahko)
    res <- do_groupChromPeaks_density(fts, sampleGroups = grps)
    res_2 <- do_groupChromPeaks_density(fts, sampleGroups = grps,
                                      minFraction = 0.9)
    expect_true(nrow(res) > nrow(res_2))
})

test_that("do_groupPeaks_mzClust works", {
    fts <- peaks(fticr_xs)
    res <- do_groupPeaks_mzClust(peaks = fts,
                                 sampleGroups = sampclass(fticr_xs))
    res_2 <- do_groupPeaks_mzClust(peaks = fts,
                                   sampleGroups = sampclass(fticr_xs),
                                   minFraction = 0, absMz = 2)
    expect_true(nrow(res$featureDefinitions) > nrow(res_2$featureDefinitions))

    res_x <- group(fticr_xs, method = "mzClust")
    expect_equal(res_x@groups, res$featureDefinitions)
    expect_equal(res_x@groupidx, res$peakIndex)

    ## Issue 416
    nas <- sample(1:nrow(fts), 10)
    fts[nas, "mz"] <- NA
    expect_warning(res <- .fix_mz_clust_peaks(fts), "Replaced them with the")
    expect_false(any(is.na(res[, "mz"])))
    fts <- fts[, !(colnames(fts) %in% c("mzmin", "mzmax"))]
    expect_error(res <- .fix_mz_clust_peaks(fts), "peaks with missing")
})

test_that("do_groupChromPeaks_nearest works", {
    xs <- faahko_xs
    features <- peaks(xs)
    sampleGroups <- sampclass(xs)
    mzVsRtBalance <- 10
    mzCheck <- 0.2
    rtCheck <- 15
    kNN <- 10

    res <- do_groupChromPeaks_nearest(features, sampleGroups)
    res_2 <- do_groupChromPeaks_nearest(features, sampleGroups, absRt = 3)
    expect_true(nrow(res$featureDefinitions) < nrow(res_2$featureDefinitions))
    res_x <- group(xs, method = "nearest")
    expect_equal(res_x@groups, res$featureDefinitions)
})

test_that(".group_peaks_density works", {
    x <- rbind(c(rt = 3.1, mz = 3, index = 1, sample = 1, into = 120),
               c(rt = 3.2, mz = 3, index = 2, sample = 2, into = 130),
               c(rt = 3.15, mz = 3, index = 3, sample = 3, into = 29),
               c(rt = 5, mz = 3, index = 4, sample = 4, into = 32),
               c(rt = 3, mz = 3, index = 5, sample = 5, into = 43),
               c(rt = 6, mz = 3, index = 6, sample = 6, into = 35))
    densFrom <- 1
    densTo <- 100
    densN <- 100
    sampleGroups <- c("b", "a", "b", "a", "b", "a")
    sampleGroupTable <- table(sampleGroups)
    bw <- 2
    maxFeatures <- 20
    minFraction <- 0.8
    minSamples <- 1
    res <- .group_peaks_density(x, bw = bw, densFrom = densFrom,
                                densTo = densTo, densN = densN,
                                sampleGroups = sampleGroups,
                                sampleGroupTable = sampleGroupTable,
                                minFraction = minFraction,
                                minSamples = minSamples,
                                maxFeatures = maxFeatures, sleep = 0)
    expect_true(nrow(res) == 1)
    expect_equal(res$peakidx, list(1:6))
    res <- .group_peaks_density(x, bw = bw, densFrom = densFrom,
                                densTo = densTo, densN = densN,
                                sampleGroups = sampleGroups,
                                sampleGroupTable = sampleGroupTable,
                                minFraction = minFraction,
                                minSamples = 7,
                                maxFeatures = maxFeatures, sleep = 0)
    expect_true(nrow(res) == 0)
    expect_true(is(res, "data.frame"))
})
