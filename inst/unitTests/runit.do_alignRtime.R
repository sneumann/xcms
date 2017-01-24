## retention time correction methods.

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

    pksAdj <- xcms:::.applyRtAdjustmentToFeatures(peaks(xsg),
                                                  rtraw = xsa@rt$raw,
                                                  rtadj = xsa@rt$corrected)
    checkEquals(pksAdj, peaks(xsa))
    ## Reset em.
    pksRaw <- xcms:::.applyRtAdjustmentToFeatures(pksAdj,
                                                  rtraw = xsa@rt$corrected,
                                                  rtadj = xsa@rt$raw)
    checkEquals(pksRaw, peaks(xsg))
}
