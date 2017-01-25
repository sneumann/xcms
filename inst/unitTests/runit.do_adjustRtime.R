## retention time correction methods.

test_adjustRtime_FeatureGroups <- function() {
    xod <- faahko_xod
    xs <- faahko_xs

    ## Group these
    xsg <- group(xs)
    xodg <- groupFeatures(xod,
                          param = FeatureDensityParam(sampleGroups = xs$class))
    checkEquals(peaks(xsg), features(xodg))
    checkEquals(xsg@groupidx, featureGroups(xodg)$featureidx)
    checkTrue(length(processHistory(xodg,
                                    type = xcms:::.PROCSTEP.FEATURE.DETECTION)) == 1)
    checkTrue(length(processHistory(xodg,
                                    type = xcms:::.PROCSTEP.FEATURE.ALIGNMENT)) == 1)
    ## Now do the retention time correction
    xsr <- retcor(xsg, method = "peakgroups")
    minFr <- (length(fileNames(xod)) - 1) / length(fileNames(xod))
    p <- FeatureGroupsParam(minFraction = minFr)
    xodr <- adjustRtime(xodg, param = p)
    ## Check that we've got process histories.
    checkTrue(validObject(xodr))
    checkTrue(hasDetectedFeatures(xodr))
    checkTrue(!hasAlignedFeatures(xodr))
    ## But we would like to keep the related process history step:
    checkTrue(hasAdjustedRtime(xodr))
    checkTrue(hasAlignedFeatures(xodg))
    ## We want to keep the process history step of the feature alignment!
    checkTrue(length(processHistory(xodr,
                                    type = xcms:::.PROCSTEP.FEATURE.ALIGNMENT)) == 1)
    checkTrue(length(processHistory(xodr,
                                    type = xcms:::.PROCSTEP.RTIME.CORRECTION)) == 1)
    ## Different from original:
    checkTrue(sum(features(xod)[, "rt"] != features(xodr)[, "rt"]) > 200)
    checkTrue(sum(features(xod)[, "rtmin"] != features(xodr)[, "rtmin"]) > 200)
    checkTrue(sum(features(xod)[, "rtmax"] != features(xodr)[, "rtmax"]) > 200)
    ## between xcmsSet and XCMSnExp
    checkEquals(features(xodr), peaks(xsr))
    ## To compare the adjusted retention time we have to extract it by sample!
    ## Otherwise the ordering will not be the same, as rtime is ordered by
    ## retention time, but @rt$raw by sample.
    checkEquals(unlist(adjustedRtime(xodr, bySample = TRUE), use.names = FALSE),
                unlist(xsr@rt$corrected, use.names = FALSE))
    ## Just to ensure - are the raw rt the same?
    checkEquals(unlist(rtime(xod, bySample = TRUE), use.names = FALSE),
                unlist(xs@rt$raw, use.names = FALSE))
    ## Doing an additional grouping
    xodrg <- groupFeatures(xodr, param = FeatureDensityParam(sampleGroups =
                                                                 xs$class))
    checkTrue(length(processHistory(xodrg,
                                    type = xcms:::.PROCSTEP.FEATURE.ALIGNMENT)) == 2)
    checkTrue(hasAdjustedRtime(xodrg))
    checkTrue(hasAlignedFeatures(xodrg))
    xsrg <- group(xsr)
    checkEquals(xsrg@groupidx, featureGroups(xodrg)$featureidx)
    
    ## Mod settings:
    xsr <- retcor(xsg, method = "peakgroups", missing = 0, span = 1)
    xodr <- adjustRtime(xodg, param = FeatureGroupsParam(minFraction = 1,
                                                         span = 1))
    checkEquals(features(xodr), peaks(xsr))
    checkEquals(unlist(adjustedRtime(xodr, bySample = TRUE), use.names = FALSE),
                unlist(xsr@rt$corrected, use.names = FALSE))

    xsr <- retcor(xsg, method = "peakgroups", missing = 0, span = 1,
                  smooth = "linear")
    xodr <- adjustRtime(xodg, param = FeatureGroupsParam(minFraction = 1,
                                                         span = 1,
                                                         smooth = "linear"))
    checkEquals(features(xodr), peaks(xsr))
    checkEquals(unlist(adjustedRtime(xodr, bySample = TRUE), use.names = FALSE),
                unlist(xsr@rt$corrected, use.names = FALSE))

    xsr <- retcor(xsg, method = "peakgroups", missing = 0, span = 1,
                  family = "symmetric")
    xodr <- adjustRtime(xodg, param = FeatureGroupsParam(minFraction = 1,
                                                         span = 1,
                                                         family = "symmetric"))
    checkEquals(features(xodr), peaks(xsr))
    checkEquals(unlist(adjustedRtime(xodr, bySample = TRUE), use.names = FALSE),
                unlist(xsr@rt$corrected, use.names = FALSE))
}

## This is to ensure that the original code works with the new one using the
## do_ function
dontrun_test_retcor.peakgroups <- function() {
    xs <- faahko
    xsg <- group(xs)

    res_1 <- retcor.peakgroups(xsg)
    res_2 <- xcms:::.retcor.peakgroups_orig(xsg)
    checkEquals(unlist(res_1@rt$corrected, use.names = FALSE),
                unlist(res_2@rt$corrected, use.names = FALSE))
    checkEquals(res_1, res_2)
    
    res_1 <- retcor.peakgroups(xsg, missing = 2)
    res_2 <- xcms:::.retcor.peakgroups_orig(xsg, missing = 2)
    checkEquals(unlist(res_1@rt$corrected, use.names = FALSE),
                unlist(res_2@rt$corrected, use.names = FALSE))
    checkEquals(res_1, res_2)

    res_1 <- retcor.peakgroups(xsg, extra = 3)
    res_2 <- xcms:::.retcor.peakgroups_orig(xsg, extra = 3)
    checkEquals(unlist(res_1@rt$corrected, use.names = FALSE),
                unlist(res_2@rt$corrected, use.names = FALSE))
    checkEquals(res_1, res_2)

    res_1 <- retcor.peakgroups(xsg, smooth = "linear")
    res_2 <- xcms:::.retcor.peakgroups_orig(xsg, smooth = "linear")
    checkEquals(unlist(res_1@rt$corrected, use.names = FALSE),
                unlist(res_2@rt$corrected, use.names = FALSE))
    checkEquals(res_1, res_2)

    res_1 <- retcor.peakgroups(xsg, span = 1)
    res_2 <- xcms:::.retcor.peakgroups_orig(xsg, span = 1)
    checkEquals(unlist(res_1@rt$corrected, use.names = FALSE),
                unlist(res_2@rt$corrected, use.names = FALSE))
    checkEquals(res_1, res_2)

    res_1 <- retcor.peakgroups(xsg, family = "symmetric")
    res_2 <- xcms:::.retcor.peakgroups_orig(xsg, family = "symmetric")
    checkEquals(unlist(res_1@rt$corrected, use.names = FALSE),
                unlist(res_2@rt$corrected, use.names = FALSE))
    checkEquals(res_1, res_2)
    
    res_1 <- retcor.peakgroups(xsg, plottype = "deviation")
    res_2 <- xcms:::.retcor.peakgroups_orig(xsg, plottype = "deviation")
    checkEquals(unlist(res_1@rt$corrected, use.names = FALSE),
                unlist(res_2@rt$corrected, use.names = FALSE))
    checkEquals(res_1, res_2)

    res_1 <- retcor.peakgroups(xsg, plottype = "mdevden")
    res_2 <- xcms:::.retcor.peakgroups_orig(xsg, plottype = "mdevden")
    checkEquals(unlist(res_1@rt$corrected, use.names = FALSE),
                unlist(res_2@rt$corrected, use.names = FALSE))
    checkEquals(res_1, res_2)
}

## That's to evaluate the do_ function with the original code. Once the
## retcor.peakgroups calls the do_function we rename it to dontrun.
test_do_adjustRtime_featureGroups_implementation <- function() {
    xs <- faahko
    xsg <- group(xs)
    
    misSamp <- 1
    xsa <- retcor(xsg, method = "peakgroups", missing = misSamp)

    minFr <- (length(sampnames(xs)) - misSamp) / length(sampnames(xs))
    res <- do_adjustRtime_featureGroups(features = peaks(xs),
                                        featureIndex = xsg@groupidx,
                                        rtime = xsg@rt$raw,
                                        minFraction = minFr)
    checkEquals(xsa@rt$corrected, res)

    ## Change settings.
    misSamp <- 3
    xsa <- retcor(xsg, method = "peakgroups", missing = misSamp)

    minFr <- (length(sampnames(xs)) - misSamp) / length(sampnames(xs))
    res <- do_adjustRtime_featureGroups(features = peaks(xs),
                                        featureIndex = xsg@groupidx,
                                        rtime = xsg@rt$raw,
                                        minFraction = minFr)
    checkEquals(xsa@rt$corrected, res)

    misSamp <- 2
    xtr <- 2
    xsa <- retcor(xsg, method = "peakgroups", missing = misSamp, extra = xtr)

    minFr <- (length(sampnames(xs)) - misSamp) / length(sampnames(xs))
    res <- do_adjustRtime_featureGroups(features = peaks(xs),
                                        featureIndex = xsg@groupidx,
                                        rtime = xsg@rt$raw,
                                        minFraction = minFr, extraFeatures = xtr)
    checkEquals(xsa@rt$corrected, res)

    xsa <- retcor(xsg, method = "peakgroups", missing = misSamp, extra = xtr,
                  smooth = "linear")
    minFr <- (length(sampnames(xs)) - misSamp) / length(sampnames(xs))
    res <- do_adjustRtime_featureGroups(features = peaks(xs),
                                        featureIndex = xsg@groupidx,
                                        rtime = xsg@rt$raw,
                                        minFraction = minFr, extraFeatures = xtr,
                                        smooth = "linear")
    checkEquals(xsa@rt$corrected, res)

    xsa <- retcor(xsg, method = "peakgroups", missing = misSamp, extra = xtr,
                  family = "symmetric")
    minFr <- (length(sampnames(xs)) - misSamp) / length(sampnames(xs))
    res <- do_adjustRtime_featureGroups(features = peaks(xs),
                                        featureIndex = xsg@groupidx,
                                        rtime = xsg@rt$raw,
                                        minFraction = minFr, extraFeatures = xtr,
                                        family = "symmetric")
    checkEquals(xsa@rt$corrected, res)

    xsa <- retcor(xsg, method = "peakgroups", missing = misSamp, extra = xtr,
                  span = 1)
    minFr <- (length(sampnames(xs)) - misSamp) / length(sampnames(xs))
    res <- do_adjustRtime_featureGroups(features = peaks(xs),
                                        featureIndex = xsg@groupidx,
                                        rtime = xsg@rt$raw,
                                        minFraction = minFr, extraFeatures = xtr,
                                        span = 1)
    checkEquals(xsa@rt$corrected, res)
}

dontrun_do_adjustRtime_peakgroups_implementation <- function() {
    library(xcms)
    library(RUnit)
    faahko_files <- c(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
                      system.file('cdf/KO/ko16.CDF', package = "faahKO"),
                      system.file('cdf/KO/ko18.CDF', package = "faahKO"),
                      system.file('cdf/WT/wt15.CDF', package = "faahKO"),
                      system.file('cdf/WT/wt16.CDF', package = "faahKO"),
                      system.file('cdf/WT/wt18.CDF', package = "faahKO"))

    od <- readMSData2(faahko_files)

    xod <- detectFeatures(od, param = CentWaveParam(noise = 100,
                                                    snthresh = 20))
    xs <- xcmsSet(faahko_files, profparam = list(step = 0),
                  method = "centWave", noise = 100, snthresh = 20)

    checkEquals(features(xod), peaks(xs))
    ## feature grouping
    p <- FeatureDensityParam(sampleGroups = rep(c("KO", "WT"), each = 3))
    xod <- groupFeatures(xod, param = p)
    xs <- group(xs, method = "density")
    
    ## Feature alignment on those:
    xs <- retcor(xs, method = "peakgroups")

    ##
    minFr <- 5/6
    res <- xcms:::do_adjustRtime_featureGroups(features = features(xod),
                                               featureGroups(xod)$featureidx,
                                               rtime = rtime(xod, bySample = TRUE),
                                               minFraction = minFr)
    a <- unname(unlist(res, use.names = FALSE))
    b <- unlist(xs@rt$corrected, use.names = FALSE)
    ## Now, they are slightly different now, because we order by median rt of the
    ## actual peaks, and they by the median rt of the whole peak group.
    checkEquals(a, b)
    checkEquals(res, unname(xs@rt$corrected))

    ## Manually correcting the guys:
    rtr <- rtime(xod, bySample = TRUE)[[1]]
    rtc <- res[[1]]
    ## rtdevsmo <- rtr - rtc

    ## That's strange!
    ## cfun <- stepfun(rtr[-1] - diff(rtr) / 2, rtr - rtdevsmo)
    cfun <- stepfun(rtr[-1] - diff(rtr) / 2, rtc)
    corFeat <- features(xod)
    whichSamp <- which(corFeat[, "sample"] == 1)
    corFeat[whichSamp, c("rt", "rtmin", "rtmax")] <-
        cfun(corFeat[whichSamp, c("rt", "rtmin", "rtmax")])

    checkEquals(corFeat[whichSamp, ], peaks(xs)[whichSamp, ])
    checkTrue(any(peaks(xs) != features(xod)))
    
    ## Do the backwards correction.
    adjFun <- stepfun(rtc[-1] - diff(rtc) / 2, rtr)
    origFeats <- corFeat
    origFeats[whichSamp, c("rt", "rtmin", "rtmax")] <-
        adjFun(corFeat[whichSamp, c("rt", "rtmin", "rtmax")])
    checkEquals(features(xod)[whichSamp, ], origFeats[whichSamp, ])
    ## OK.
}

## Testing the internal .applyRtAdjustment function.
test_applyRtAdjustment <- function() {
    xs <- faahko
    ## group em.
    xsg <- group(xs)
    ## align em.
    xsa <- retcor(xsg, method = "peakgroups")

    pksAdj <- xcms:::.applyRtAdjToFeatures(peaks(xsg),
                                           rtraw = xsa@rt$raw,
                                           rtadj = xsa@rt$corrected)
    checkEquals(pksAdj, peaks(xsa))
    ## Reset em.
    pksRaw <- xcms:::.applyRtAdjToFeatures(pksAdj,
                                           rtraw = xsa@rt$corrected,
                                           rtadj = xsa@rt$raw)
    checkEquals(pksRaw, peaks(xsg))
}
