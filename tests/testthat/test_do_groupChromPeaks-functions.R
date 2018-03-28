test_that("do_groupChromPeaks_density works", {
    fts <- peaks(faahko)
    res <- do_groupChromPeaks_density(fts, sampleGroups = sampclass(faahko))
    res_2 <- do_groupChromPeaks_density(fts, sampleGroups = sampclass(faahko),
                                      minFraction = 0.9)
    expect_true(nrow(res$featureDefinitions) > nrow(res_2$featureDefinitions))
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
