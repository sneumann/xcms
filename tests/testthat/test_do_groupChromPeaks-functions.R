test_that("do_groupChromPeaks_density works", {
    skip_on_os(os = "windows", arch = "i386")

    fts <- chromPeaks(xod_xg)
    grps <- rep(1, 3)
    res <- do_groupChromPeaks_density(fts, sampleGroups = grps)
    res_2 <- do_groupChromPeaks_density(fts, sampleGroups = grps,
                                        minFraction = 0.9)
    expect_true(nrow(res) > nrow(res_2))

    res_3 <- do_groupChromPeaks_density(fts, sampleGroups = grps, ppm = 20)
    expect_equal(nrow(res), nrow(res_3))

    ## Different ways to calculate rtmed:
    expect_error(do_groupChromPeaks_density(fts, sampleGroups = grps,
                                            rtCenterFun = "sum"),
                 "should be one of")
    res_2 <- do_groupChromPeaks_density(fts, sampleGroups = grps,
                                        rtCenterFun = "mean")
    expect_true(sum(res_2$rtmed != res$rtmed) > 10)
    res_3 <- do_groupChromPeaks_density(fts, sampleGroups = grps,
                                        rtCenterFun = "wMean")
    expect_true(sum(res_2$rtmed != res_3$rtmed) > 10)
    expect_true(sum(res$rtmed != res_3$rtmed) > 10)
})

test_that("do_groupPeaks_mzClust works", {
    skip_on_os(os = "windows", arch = "i386")

    fts <- chromPeaks(fticr_xod)
    res <- do_groupPeaks_mzClust(peaks = fts,
                                 sampleGroups = c(1, 1))
    res_2 <- do_groupPeaks_mzClust(peaks = fts,
                                   sampleGroups = c(1, 1),
                                   minFraction = 0, absMz = 2)
    expect_true(nrow(res$featureDefinitions) > nrow(res_2$featureDefinitions))

    ## Issue 416
    nas <- sample(1:nrow(fts), 10)
    fts[nas, "mz"] <- NA
    expect_warning(res <- .fix_mz_clust_peaks(fts), "Replaced them with the")
    expect_false(any(is.na(res[, "mz"])))
    fts <- fts[, !(colnames(fts) %in% c("mzmin", "mzmax"))]
    expect_error(res <- .fix_mz_clust_peaks(fts), "peaks with missing")
})

test_that("do_groupChromPeaks_nearest works", {
    skip_on_os(os = "windows", arch = "i386")

    tmp <- filterFile(xod_xg)
    features <- chromPeaks(tmp)
    sampleGroups <- rep(1, length(fileNames(tmp)))
    mzVsRtBalance <- 10
    mzCheck <- 0.2
    rtCheck <- 15
    kNN <- 10

    res <- do_groupChromPeaks_nearest(features, sampleGroups)
    res_2 <- do_groupChromPeaks_nearest(features, sampleGroups, absRt = 3)
    expect_true(nrow(res$featureDefinitions) < nrow(res_2$featureDefinitions))
})

test_that(".group_peaks_density works", {
    skip_on_os(os = "windows", arch = "i386")

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
    res_2 <- .group_peaks_density(x, bw = bw, densFrom = densFrom,
                                  densTo = densTo, densN = densN,
                                  sampleGroups = sampleGroups,
                                  sampleGroupTable = sampleGroupTable,
                                  minFraction = minFraction,
                                  minSamples = minSamples,
                                  maxFeatures = maxFeatures, sleep = 0,
                                  rtFun = .feature_rt_wmean)
    expect_equal(res$mzmed, res_2$mzmed)
    expect_equal(res$mzmin, res_2$mzmin)
    expect_equal(res$mzmax, res_2$mzmax)
    expect_equal(res$rtmin, res_2$rtmin)
    expect_equal(res$rtmax, res_2$rtmax)
    expect_equal(res$peakidx, res_2$peakidx)
    expect_true(res$rtmed < res_2$rtmed)
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

test_that("do_groupChromPeaks_density works with skipping samples", {
    x <- loadXcmsData("xmse")
    pks <- chromPeaks(x)
    ## Errors
    expect_error(do_groupChromPeaks_density(pks), "sampleGroups")
    expect_error(do_groupChromPeaks_density(3, sampleGroups = 3), "matrix")
    expect_error(do_groupChromPeaks_density(pks[, 1:3], sampleGroups = 1),
                 "not found")
    expect_error(do_groupChromPeaks_density(pks, sampleGroups = 1:3),
                 "Sample indices")

    ## groups for all samples.
    grps <- sampleData(x)$sample_group
    res <- do_groupChromPeaks_density(pks, sampleGroups = grps,
                                      minFraction = 1, bw = 30)
    expect_true(all(res$WT == 4 | res$KO == 4))
    expect_true(all(res$WT <= 4))
    expect_true(all(res$KO <= 4))

    res_2 <- do_groupChromPeaks_density(
        pks, sampleGroups = rep(1, length(grps)), minFraction = 1)
    expect_true(nrow(res_2) < nrow(res))
    expect_true(all(res_2$`1` == 8))

    ## using only one sample group
    grps[grps == "KO"] <- NA
    res_3 <- do_groupChromPeaks_density(pks, sampleGroups = grps,
                                        minFraction = 1)
    expect_true(nrow(res_3) < nrow(res))
    expect_true(all(res_3$WT == 4))
    expect_equal(nrow(res_3), sum(res$WT == 4))
    tmp <- res[res$WT == 4, ]
})

test_that("rtCenterFun functions work", {
    x <- cbind(rt = c(200, 231.1, 323.1, 232.1),
               maxo = c(200.1, 240.1, 11322, 245.2))

    res <- .feature_rt_mean(x)
    expect_equal(res, mean(x[, "rt"]))
    resw <- .feature_rt_wmean(x, w = "maxo")
    expect_true(resw > res)
    expect_equal(resw, weighted.mean(x[, "rt"], x[, "maxo"]))
    res <- .feature_rt_median(x)
    expect_equal(res, median(x[, "rt"]))
})
