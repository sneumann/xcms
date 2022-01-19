xodg <- xod_xgrg
library(MsFeatures)
xodgg <- groupFeatures(xodg, param = SimilarRtimeParam(4))
xodgg <- groupFeatures(xodgg, param = AbundanceSimilarityParam(threshold = 0.3))


test_that("featureGroups,featureGroups<-,XCMSnExp works", {
    skip_on_os(os = "windows", arch = "i386")

    expect_error(featureGroups(xod_x), "Please run")
    res <- featureGroups(xodg)
    expect_true(all(is.na(res)))
    tmp <- xodg
    featureGroups(tmp) <- "a"
    expect_true(all(featureGroups(tmp) == "a"))

    expect_error(featureGroups(xod_x) <- "a", "Please run")
    expect_error(featureGroups(xodg) <- 1:2, "length")
})

test_that("SimilarRtimeParam works", {
    skip_on_os(os = "windows", arch = "i386")

    prm <- SimilarRtimeParam(3)

    expect_error(groupFeatures(xod_x, prm), "No feature definitions")
    expect_error(groupFeatures(xodg, prm, msLevel = 1:2), "single MS")
    res <- groupFeatures(xodg, prm)
    expect_true(any(colnames(featureDefinitions(res)) == "feature_group"))
    expect_false(any(is.na(featureGroups(res))))
    expect_true(is.character(featureGroups(res)))

    res2 <- groupFeatures(xodg,
                          SimilarRtimeParam(3, groupFun = MsCoreUtils::group))
    expect_true(length(table(featureGroups(res2))) <
                length(table(featureGroups(res))))

    ## Different MS levels
    tmp <- xodg
    idx <- c(1:3, 5, 45, 47)
    featureDefinitions(tmp)$ms_level[idx] <- 2
    res <- groupFeatures(tmp, prm)
    expect_true(all(is.na(featureGroups(res))[idx]))
    expect_false(any(is.na(featureGroups(res))[-idx]))
    res <- groupFeatures(tmp, prm, msLevel = 2L)
    expect_false(any(is.na(featureGroups(res))[idx]))
    expect_true(all(is.na(featureGroups(res))[-idx]))

    ## Pre-defined groups
    fgs <- rep("AB", nrow(featureDefinitions(xodg)))
    fgs[idx] <- NA
    tmp <- xodg
    featureGroups(tmp) <- fgs
    res <- groupFeatures(tmp, prm)
    expect_true(all(is.na(featureGroups(res))[idx]))
    expect_false(any(is.na(featureGroups(res))[-idx]))
})

test_that("AbundanceSimilarityParam works", {
    skip_on_os(os = "windows", arch = "i386")

    prm <- AbundanceSimilarityParam(threshold = 0.5, value = "maxo")
    expect_equal(prm@threshold, 0.5)
    expect_equal(prm@dots, list(value = "maxo"))

    expect_error(AbundanceSimilarityParam(subset = "4"), "integer")

    expect_error(groupFeatures(xod_x, AbundanceSimilarityParam()), "feature")
    expect_error(
        groupFeatures(xodg, AbundanceSimilarityParam(subset = c(1, 4, 5))),
        "should be between")

    res <- groupFeatures(xodg, AbundanceSimilarityParam())
    expect_true(any(colnames(featureDefinitions(res)) == "feature_group"))
    expect_true(length(unique(featureDefinitions(res)$feature_group)) <
                nrow(featureDefinitions(res)))
    res_2 <- groupFeatures(xodg, AbundanceSimilarityParam(subset = c(2, 3)))

    plotFeatureGroups(res_2)
    expect_error(plotFeatureGroups(res_2, featureGroups = "a"), "None of the")
    expect_error(plotFeatureGroups(xodg), "None of the")

    ## With pre-defined grps.
    tmp <- xodg
    featureDefinitions(tmp)$feature_group <- "FG.2"
    idx <- c(4, 12, 23, 46)
    featureDefinitions(tmp)$ms_level[idx] <- 2

    res <- groupFeatures(tmp, AbundanceSimilarityParam(), msLevel = 1)
    expect_true(all(featureGroups(res)[idx] == "FG.2"))
    expect_true(all(featureGroups(res)[-idx] != "FG.2"))
    res_2 <- groupFeatures(tmp, AbundanceSimilarityParam(), msLevel = 2)
    expect_true(all(featureGroups(res_2)[-idx] == "FG.2"))
    expect_true(all(featureGroups(res_2)[idx] != "FG.2"))

    tmp <- quantify(xodg, filled = TRUE, method = "sum", value = "maxo")
    res <- groupFeatures(xodg, AbundanceSimilarityParam(), filled = TRUE,
                         method = "sum", value = "maxo")
    res_2 <- groupFeatures(tmp, AbundanceSimilarityParam())
    expect_equal(featureGroups(res), featureGroups(res_2))
})

## test_that("featureGroupPseudoSpectrum works", {
    ## skip_on_os(os = "windows", arch = "i386")

##     fvals <- featureValues(xodgg, method = "maxint", value = "maxo")
##     ## 3 feature
##     ft_idx <- which(featureGroups(xodgg) == "FG.010.001")
##     res <- featureGroupPseudoSpectrum("FG.010.001", xodgg, fvals = fvals,
##                                       intensityFun = median)
##     expect_true(is(res, "Spectrum1"))
##     expect_true(peaksCount(res) == 3)
##     expect_true(validObject(res))
##     expect_equal(intensity(res), apply(fvals[ft_idx, ], MARGIN = 1,
##                                        median, na.rm = TRUE))
##     expect_equal(mz(res), featureDefinitions(xodgg)$mzmed[ft_idx])
##     expect_equal(rtime(res), median(featureDefinitions(xodgg)$rtmed[ft_idx]))

##     ## 1 feature
##     res <- featureGroupPseudoSpectrum("FG.010.002", xodgg, fvals = fvals,
##                                       intensityFun = median)
##     ft_idx <- which(featureGroups(xodgg) == "FG.010.002")
##     expect_true(is(res, "Spectrum1"))
##     expect_true(peaksCount(res) == 1)
##     expect_true(validObject(res))
##     expect_equal(unname(intensity(res)),
##                  unname(median(fvals[ft_idx, ], na.rm = TRUE)))
##     expect_equal(mz(res), featureDefinitions(xodgg)$mzmed[ft_idx])
##     expect_equal(rtime(res), median(featureDefinitions(xodgg)$rtmed[ft_idx]))

##     expect_error(
##         featureGroupPseudoSpectrum("FG.009.1", xodgg, fvals = fvals, n = 12),
##         "has to be an integer")
## })

## test_that("featureGroupFullScan works", {
    ## skip_on_os(os = "windows", arch = "i386")

##     fvals <- featureValues(xodgg, method = "maxint", value = "maxo")
##     ## 3 feature
##     res <- featureGroupFullScan("FG.010.001", xodgg, fvals = fvals)
##     ft_idx <- which(featureGroups(xodgg) == "FG.010.001")
##     expect_true(is(res, "Spectrum1"))
##     expect_true(
##         abs(rtime(res) -
##             median(featureDefinitions(xodgg)[ft_idx, "rtmed"])) < 0.1)
##     expect_true(all(featureDefinitions(xodgg)[ft_idx, "mzmed"] %in% mz(res)))

##     ## 1 feature
##     res <- featureGroupFullScan("FG.010.002", xodgg, fvals = fvals)
##     ft_idx <- which(featureGroups(xodgg) == "FG.010.002")
##     expect_true(is(res, "Spectrum1"))
##     expect_true(
##         abs(rtime(res) -
##             median(featureDefinitions(xodgg)[ft_idx, "rtmed"])) < 0.8)
##     expect_true(all(featureDefinitions(xodgg)[ft_idx, "mzmed"] %in% mz(res)))
## })

## test_that("featureGroupSpectra works", {
    ## skip_on_os(os = "windows", arch = "i386")

##     ## Errors
##     expect_error(featureGroupSpectra(xodgg, subset = 1:5), "an integer")
##     expect_error(featureGroupSpectra(xod), "feature definitions")
##     expect_error(featureGroupSpectra(xodgg, featureGroup = c("a")), "all feature")

##     ## Get all of them
##     res_all <- featureGroupSpectra(xodgg)
##     expect_true(is(res_all, "MSpectra"))
##     expect_equal(mcols(res_all)$feature_group, unique(featureGroups(xodgg)))
##     expect_equal(unname(peaksCount(res_all)),
##                  unname(lengths(mcols(res_all)$feature_id)))

##     ## Get them in a subset
##     res_sub <- featureGroupSpectra(xodgg, subset = c(1, 3))
##     expect_true(sum(is.na(rtime(res_sub))) == 59)

##     ## Get only selected ones
##     res <- featureGroupSpectra(
##         xodgg, featureGroup = c("FG.010.001", "FG.010.002"))
##     expect_true(length(res) == 2)
##     expect_equal(mcols(res)$feature_group, c("FG.010.001", "FG.010.002"))
##     idx <- which(mcols(res_all)$feature_group %in% c("FG.010.001", "FG.010.002"))
##     expect_equal(res[[1]], res_all[[idx[1]]])
##     expect_equal(res[[2]], res_all[[idx[2]]])
## })

test_that(".group_eic_similarity works", {
    skip_on_os(os = "windows", arch = "i386")

    set.seed(123)
    chr1 <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
                         intensity = c(5, 29, 50, NA, 100, 12, 3, 4, 1, 3))
    chr2 <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
                         intensity = c(80, 50, 20, 10, 9, 4, 3, 4, 1, 3))
    chr3 <- Chromatogram(rtime = 3:9 + rnorm(7, sd = 0.3),
                         intensity = c(53, 80, 130, 15, 5, 3, 2))
    chrs <- MChromatograms(list(chr1, chr2, chr3))

    res <- .group_eic_similarity(chrs, ALIGNFUNARGS = list(method = "closest"))
    expect_true(is.factor(res))
    expect_equal(res, factor(c(1L, 2L, 1L)))
    res <- .group_eic_similarity(
        chrs, ALIGNFUNARGS = list(method = "closest", tolerance = 0))
    expect_equal(res, factor(c(1L, 2L, 3L)))

    chrs <- MChromatograms(list(chr1, chr2, chr3, chr1, chr2, chr3), ncol = 2)
    res <- .group_eic_similarity(chrs, aggregationFun = mean,
                                 ALIGNFUNARGS = list(method = "closest"))
    expect_equal(res, factor(c(1L, 2L, 1L)))
    res <- .group_eic_similarity(chrs, aggregationFun = max,
                                 ALIGNFUNARGS = list(method = "closest"))
    expect_equal(res, factor(c(1L, 2L, 1L)))
    res <- .group_eic_similarity(chrs, aggregationFun = min,
                                 ALIGNFUNARGS = list(method = "closest"))
    expect_equal(res, factor(c(1L, 2L, 1L)))

    chrs <- MChromatograms(list(chr1, chr2, chr3, chr1, chr2, chr3,
                                chr2, chr3, chr1), ncol = 3)
    res <- .group_eic_similarity(chrs, ALIGNFUNARGS = list(method = "closest"))
    expect_true(is.factor(res))
    expect_equal(res, factor(c(1L, 2L, 3L)))

    res <- .group_eic_similarity(chrs, aggregationFun = max,
                                 threshold = 0.1,
                                 ALIGNFUNARGS = list(method = "closest"))
    expect_true(is.factor(res))
    expect_equal(res, factor(c(1L, 1L, 1L)))

    res <- .group_eic_similarity(chrs, aggregationFun = median,
                                 ALIGNFUNARGS = list(method = "closest"))
    expect_true(is.factor(res))
    expect_equal(res, factor(c(1L, 2L, 1L)))
})

test_that("EicSimilarityParam works", {
    skip_on_os(os = "windows", arch = "i386")

    res <- EicSimilarityParam()
    expect_equal(res@threshold, 0.9)
    expect_equal(res@ALIGNFUNARGS, list(tolerance = 0, method = "closest"))
    expect_equal(res@ALIGNFUN, alignRt)
    expect_equal(res@FUN, stats::cor)
    expect_equal(res@FUNARGS, list(use = "pairwise.complete.obs"))
    expect_equal(res@n, 1L)
    expect_equal(res@onlyPeak, TRUE)
    expect_equal(res@dots, list())

    res <- EicSimilarityParam(FUN = dist)
    expect_equal(res@FUN, dist)
    res <- EicSimilarityParam(ALIGNFUN = sum)
    expect_equal(res@ALIGNFUN, sum)
    res <- EicSimilarityParam(groupFun = max)
    expect_equal(res@groupFun, max)
    res <- EicSimilarityParam(threshold = 0, n = 10, onlyPeak = FALSE)
    expect_equal(res@threshold, 0)
    expect_equal(res@n, 10)
    expect_equal(res@onlyPeak, FALSE)
    res <- EicSimilarityParam(ALIGNFUNARGS = list(a = 4))
    expect_equal(res@ALIGNFUNARGS, list(a = 4))
    res <- EicSimilarityParam(FUNARGS = list(b = 5))
    expect_equal(res@FUNARGS, list(b = 5))
    res <- EicSimilarityParam(someother = 5)
    expect_equal(res@dots, list(someother = 5))

    expect_error(EicSimilarityParam(threshold = c(4, 3)), "positive numeric")
    expect_error(EicSimilarityParam(n = 1:2), "positive numeric")
    expect_error(EicSimilarityParam(onlyPeak = c(TRUE, FALSE)), "length 1")
    expect_error(EicSimilarityParam(value = "other"))

})

test_that("groupFeatures,EicSimilarityParam works", {
    skip_on_os(os = "windows", arch = "i386")

    ## n bigger than 3
    expect_error(groupFeatures(xodg, param = EicSimilarityParam(n = 5)),
                 "smaller than or")
    ## no feature definitions
    expect_error(groupFeatures(xod_x, param = EicSimilarityParam()), "No")
    ## MS level length > 1
    expect_error(
        groupFeatures(xodg, param = EicSimilarityParam(), msLevel = 1:2),
        "single MS level")

    tmp <- xodg
    res_all <- groupFeatures(tmp, param = EicSimilarityParam())
    expect_true(is.character(featureGroups(res_all)))

    idx <- c(3, 12, 13, 34, 39, 40)
    tmp <- xodg
    featureDefinitions(tmp)$feature_group <- NA
    featureDefinitions(tmp)$feature_group[idx] <- "FG"
    res <- groupFeatures(tmp, param = EicSimilarityParam())
    expect_true(all(is.na(featureGroups(res)[-idx])))
    expect_true(length(unique(featureGroups(res))) < length(idx))
    a <- featureGroups(res)[idx]
    b <- featureGroups(res_all)[idx]
    expect_equal(as.integer(factor(a, levels = unique(a))),
                 as.integer(factor(b, levels = unique(b))))

    featureDefinitions(tmp)$feature_group <- NULL
    featureDefinitions(tmp)$ms_level[idx] <- 2

    res_2 <- groupFeatures(tmp, param = EicSimilarityParam(), msLevel = 2)
    expect_equal(featureDefinitions(res)$feature_group,
                 featureDefinitions(res_2)$feature_group)

    res <- groupFeatures(xodgg, param = EicSimilarityParam(threshold = 0.7))
    expect_true(length(table(featureGroups(xodgg))) <
                length(table(featureGroups(res))))
})
