## retention time correction methods and functionality related to adjusted
## retention times.

test_adjustRtime_PeakGroups <- function() {
    xod <- faahko_xod
    xs <- faahko_xs

    ## Group these
    xsg <- group(xs)
    xodg <- groupChromPeaks(xod,
                            param = PeakDensityParam(sampleGroups = xs$class))
    checkEquals(peaks(xsg), chromPeaks(xodg)[, colnames(peaks(xsg))])
    checkEquals(xsg@groupidx, featureDefinitions(xodg)$peakidx)
    checkTrue(length(processHistory(xodg,
                                    type = xcms:::.PROCSTEP.PEAK.DETECTION)) == 1)
    checkTrue(length(processHistory(xodg,
                                    type = xcms:::.PROCSTEP.PEAK.GROUPING)) == 1)
    ## Now do the retention time correction
    xsr <- retcor(xsg, method = "peakgroups", missing = 0, span = 0.3)
    ## minFr <- (length(fileNames(xod)) - 1) / length(fileNames(xod))
    p <- PeakGroupsParam(minFraction = 1, span = 0.3)
    xodr <- adjustRtime(xodg, param = p)
    ## Check that we've got process histories.
    checkTrue(validObject(xodr))
    checkTrue(hasChromPeaks(xodr))
    checkTrue(!hasFeatures(xodr))
    ## But we would like to keep the related process history step:
    checkTrue(hasAdjustedRtime(xodr))
    checkTrue(hasFeatures(xodg))
    ## We want to keep the process history step of the feature alignment!
    checkTrue(length(processHistory(xodr,
                                    type = xcms:::.PROCSTEP.PEAK.GROUPING)) == 1)
    checkTrue(length(processHistory(xodr,
                                    type = xcms:::.PROCSTEP.RTIME.CORRECTION)) == 1)
    ## Different from original:
    checkTrue(sum(chromPeaks(xod)[, "rt"] != chromPeaks(xodr)[, "rt"]) > 200)
    checkTrue(sum(chromPeaks(xod)[, "rtmin"] != chromPeaks(xodr)[, "rtmin"]) > 200)
    checkTrue(sum(chromPeaks(xod)[, "rtmax"] != chromPeaks(xodr)[, "rtmax"]) > 200)
    ## between xcmsSet and XCMSnExp
    checkEquals(chromPeaks(xodr)[, colnames(peaks(xsr))], peaks(xsr))
    ## To compare the adjusted retention time we have to extract it by sample!
    ## Otherwise the ordering will not be the same, as rtime is ordered by
    ## retention time, but @rt$raw by sample.
    checkEquals(unlist(adjustedRtime(xodr, bySample = TRUE), use.names = FALSE),
                unlist(xsr@rt$corrected, use.names = FALSE))
    ## Just to ensure - are the raw rt the same?
    checkEquals(unlist(rtime(xod, bySample = TRUE), use.names = FALSE),
                unlist(xs@rt$raw, use.names = FALSE))
    ## Check that we get the same by supplying the peakGroupsMatrix.
    pgm <- adjustRtimePeakGroups(xodg, param = p)
    p_2 <- p
    minFraction(p_2) <- 0.5
    extraPeaks(p_2) <- 20
    peakGroupsMatrix(p_2) <- pgm
    xodr_2 <- adjustRtime(xodg, param = p_2)
    checkEquals(adjustedRtime(xodr), adjustedRtime(xodr_2))
    checkEquals(chromPeaks(xodr), chromPeaks(xodr_2))
    p_got <- processParam(
        processHistory(xodr, type = xcms:::.PROCSTEP.RTIME.CORRECTION)[[1]])
    peakGroupsMatrix(p_got) <- matrix(ncol = 0, nrow = 0)
    checkEquals(p_got, p)
    checkEquals(processParam(
        processHistory(xodr_2, type = xcms:::.PROCSTEP.RTIME.CORRECTION)[[1]]),
        p_2)
    ## Doing an additional grouping
    xodrg <- groupChromPeaks(xodr, param = PeakDensityParam(sampleGroups =
                                                                xs$class))
    checkTrue(length(processHistory(xodrg,
                                    type = xcms:::.PROCSTEP.PEAK.GROUPING)) == 2)
    checkTrue(hasAdjustedRtime(xodrg))
    checkTrue(hasFeatures(xodrg))
    xsrg <- group(xsr)
    checkEquals(xsrg@groupidx, featureDefinitions(xodrg)$peakidx)
    
    ## Mod settings:
    xsr <- retcor(xsg, method = "peakgroups", missing = 0, span = 1)
    xodr <- adjustRtime(xodg, param = PeakGroupsParam(minFraction = 1,
                                                         span = 1))
    checkEquals(chromPeaks(xodr)[, colnames(peaks(xsr))], peaks(xsr))
    checkEquals(unlist(adjustedRtime(xodr, bySample = TRUE), use.names = FALSE),
                unlist(xsr@rt$corrected, use.names = FALSE))

    xsr <- retcor(xsg, method = "peakgroups", missing = 0, span = 1,
                  smooth = "linear")
    xodr <- adjustRtime(xodg, param = PeakGroupsParam(minFraction = 1,
                                                         span = 1,
                                                         smooth = "linear"))
    checkEquals(chromPeaks(xodr)[, colnames(peaks(xsr))], peaks(xsr))
    checkEquals(unlist(adjustedRtime(xodr, bySample = TRUE), use.names = FALSE),
                unlist(xsr@rt$corrected, use.names = FALSE))

    xsr <- retcor(xsg, method = "peakgroups", missing = 0, span = 1,
                  family = "symmetric")
    xodr <- adjustRtime(xodg, param = PeakGroupsParam(minFraction = 1,
                                                         span = 1,
                                                         family = "symmetric"))
    checkEquals(chromPeaks(xodr)[, colnames(peaks(xsr))], peaks(xsr))
    checkEquals(unlist(adjustedRtime(xodr, bySample = TRUE), use.names = FALSE),
                unlist(xsr@rt$corrected, use.names = FALSE))
    ## Dropping results.
    tmp <- dropAdjustedRtime(xodr)
    checkEquals(tmp, xod)
}

test_getPeakGroupsRtMatrix <- function() {
    param <- PeakGroupsParam()
    nSamples <- length(fileNames(xod_xg))
    pkGrp <- xcms:::.getPeakGroupsRtMatrix(
        peaks = chromPeaks(xod_xg),
        peakIndex = xcms:::.peakIndex(xod_xg),
        nSamples = nSamples,
        missingSample = nSamples - (nSamples * minFraction(param)),
        extraPeaks = extraPeaks(param)
        )
    ## checkEquals(colnames(pkGrp), colnames(chromPeaks(xod_xg)))
    fts <- featureDefinitions(xod_xg)[rownames(pkGrp), ]
    checkTrue(all(pkGrp[, 1] >= fts$rtmin & pkGrp[, 1] <= fts$rtmax))
    checkTrue(all(pkGrp[, 2] >= fts$rtmin & pkGrp[, 2] <= fts$rtmax))
    checkTrue(all(pkGrp[, 3] >= fts$rtmin & pkGrp[, 3] <= fts$rtmax))
}

test_plotAdjustedRtime <- function() {
    plotAdjustedRtime(xod_xgr)
    plotAdjustedRtime(xod_xgrg)
    plotAdjustedRtime(xod_x)
    plotAdjustedRtime(xod_xg)
}

dontrun_issue146 <- function() {
    ## For some files it can happen that the adjusted retention times are no
    ## longer ordered increasingly.

    ## Using my data that caused the problems
    library(xcms)
    library(RUnit)
    load("/Users/jo/R-workspaces/2017/2017-03-Mitra-untargeted/data/RData/mitra-extraction/mitra.RData")
    mzWid <- 0.02
    bw_1 <- 1.5
    bw_2 <- 2
    
    ## Retention time adjustment using "PeakGroups"
    ## First grouping of samples. Setting minFraction
    pdp <- PeakDensityParam(sampleGroups = pData(mitra)$extraction_name,
                            bw = bw_1, binSize = mzWid, minFraction = 0.5,
                            maxFeatures = 200)
    mitra_pg <- groupChromPeaks(mitra, param = pdp)
    
    ## These are if we want to jump into the do_adjustRtime_peakGroups function.
    peaks <- chromPeaks(mitra_pg)
    peakIndex <- featureDefinitions(mitra_pg)$peakidx
    rtime <- rtime(mitra_pg, adjusted = FALSE, bySample = TRUE)
    minFraction <- 0.85
    extraPeaks <- 1
    span <- 0.2
    family <- "gaussian"

    ## Running the original code.
    res_o <- xcms:::do_adjustRtime_peakGroups_orig(peaks, peakIndex, rtime = rtime,
                                       minFraction = minFraction,
                                       extraPeaks = extraPeaks)
    sum(unlist(lapply(res_o, is.unsorted)))
    ## Alternative 1 - uh, does not finish???
    res_2 <- do_adjustRtime_peakGroups(peaks, peakIndex, rtime = rtime,
                                       minFraction = minFraction,
                                       extraPeaks = extraPeaks)
    sum(unlist(lapply(res_2, is.unsorted)))

    res <- adjustRtime(mitra_pg,
                       param = PeakGroupsParam(minFraction = minFraction,
                                               span = 1))
    tmp <- dropAdjustedRtime(res)
    checkEquals(chromPeaks(tmp), chromPeaks(mitra))
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
test_do_adjustRtime_peakGroups_implementation <- function() {
    xs <- faahko
    xsg <- group(xs)
    
    misSamp <- 1
    xsa <- retcor(xsg, method = "peakgroups", missing = misSamp)

    minFr <- (length(sampnames(xs)) - misSamp) / length(sampnames(xs))
    res <- do_adjustRtime_peakGroups(peaks = peaks(xs),
                                        peakIndex = xsg@groupidx,
                                        rtime = xsg@rt$raw,
                                        minFraction = minFr)
    checkEquals(xsa@rt$corrected, res)

    ## Change settings.
    misSamp <- 3
    xsa <- retcor(xsg, method = "peakgroups", missing = misSamp)

    minFr <- (length(sampnames(xs)) - misSamp) / length(sampnames(xs))
    res <- do_adjustRtime_peakGroups(peaks = peaks(xs),
                                        peakIndex = xsg@groupidx,
                                        rtime = xsg@rt$raw,
                                        minFraction = minFr)
    checkEquals(xsa@rt$corrected, res)

    misSamp <- 2
    xtr <- 2
    xsa <- retcor(xsg, method = "peakgroups", missing = misSamp, extra = xtr)

    minFr <- (length(sampnames(xs)) - misSamp) / length(sampnames(xs))
    res <- do_adjustRtime_peakGroups(peaks = peaks(xs),
                                        peakIndex = xsg@groupidx,
                                        rtime = xsg@rt$raw,
                                        minFraction = minFr, extraPeaks = xtr)
    checkEquals(xsa@rt$corrected, res)

    xsa <- retcor(xsg, method = "peakgroups", missing = misSamp, extra = xtr,
                  smooth = "linear")
    minFr <- (length(sampnames(xs)) - misSamp) / length(sampnames(xs))
    res <- do_adjustRtime_peakGroups(peaks = peaks(xs),
                                        peakIndex = xsg@groupidx,
                                        rtime = xsg@rt$raw,
                                        minFraction = minFr, extraPeaks = xtr,
                                        smooth = "linear")
    checkEquals(xsa@rt$corrected, res)

    xsa <- retcor(xsg, method = "peakgroups", missing = misSamp, extra = xtr,
                  family = "symmetric")
    minFr <- (length(sampnames(xs)) - misSamp) / length(sampnames(xs))
    res <- do_adjustRtime_peakGroups(peaks = peaks(xs),
                                        peakIndex = xsg@groupidx,
                                        rtime = xsg@rt$raw,
                                        minFraction = minFr, extraPeaks = xtr,
                                        family = "symmetric")
    checkEquals(xsa@rt$corrected, res)

    xsa <- retcor(xsg, method = "peakgroups", missing = misSamp, extra = xtr,
                  span = 1)
    minFr <- (length(sampnames(xs)) - misSamp) / length(sampnames(xs))
    res <- do_adjustRtime_peakGroups(peaks = peaks(xs),
                                        peakIndex = xsg@groupidx,
                                        rtime = xsg@rt$raw,
                                        minFraction = minFr, extraPeaks = xtr,
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

    od <- readMSData(faahko_files, mode = "onDisk")

    xod <- findChromPeaks(od, param = CentWaveParam(noise = 100,
                                                    snthresh = 20))
    xs <- xcmsSet(faahko_files, profparam = list(step = 0),
                  method = "centWave", noise = 100, snthresh = 20)

    checkEquals(chromPeaks(xod), peaks(xs))
    ## feature grouping
    p <- PeakDensityParam(sampleGroups = rep(c("KO", "WT"), each = 3))
    xod <- groupChromPeaks(xod, param = p)
    xs <- group(xs, method = "density")
    
    ## Feature alignment on those:
    xs <- retcor(xs, method = "peakgroups")

    ##
    minFr <- 5/6
    res <- xcms:::do_adjustRtime_peakGroups(peaks = chromPeaks(xod),
                                               featureDefinitions(xod)$peakidx,
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
    corFeat <- chromPeaks(xod)
    whichSamp <- which(corFeat[, "sample"] == 1)
    corFeat[whichSamp, c("rt", "rtmin", "rtmax")] <-
        cfun(corFeat[whichSamp, c("rt", "rtmin", "rtmax")])

    checkEquals(corFeat[whichSamp, ], peaks(xs)[whichSamp, ])
    checkTrue(any(peaks(xs) != chromPeaks(xod)))
    
    ## Do the backwards correction.
    adjFun <- stepfun(rtc[-1] - diff(rtc) / 2, rtr)
    origFeats <- corFeat
    origFeats[whichSamp, c("rt", "rtmin", "rtmax")] <-
        adjFun(corFeat[whichSamp, c("rt", "rtmin", "rtmax")])
    checkEquals(chromPeaks(xod)[whichSamp, ], origFeats[whichSamp, ])
    ## OK.
}

## Testing the internal .applyRtAdjustment function.
test_applyRtAdjustment <- function() {
    xs <- faahko
    ## group em.
    xsg <- group(xs)
    ## align em.
    xsa <- retcor(xsg, method = "peakgroups")

    pksAdj <- xcms:::.applyRtAdjToChromPeaks(peaks(xsg),
                                             rtraw = xsa@rt$raw,
                                             rtadj = xsa@rt$corrected)
    checkEquals(pksAdj, peaks(xsa))
    ## Reset em.
    pksRaw <- xcms:::.applyRtAdjToChromPeaks(pksAdj,
                                             rtraw = xsa@rt$corrected,
                                             rtadj = xsa@rt$raw)
    checkEquals(pksRaw, peaks(xsg))
}

## Obiwarp:
test_obiwarp <- function() {

    xs <- faahko_xs
    od <- faahko_od
    xod <- faahko_xod
    ## Feature alignment on those:
    ## object <- findChromPeaks(faahko_od, param = CentWaveParam(noise = 10000,
    ##                                                           snthresh = 40))
    prm <- ObiwarpParam(binSize = 1)
    xs_2 <- retcor.obiwarp(xs, profStep = binSize(prm))
    checkEquals(xs_2@rt$raw[[2]], xs_2@rt$corrected[[2]])
    checkTrue(sum(xs_2@rt$raw[[1]] != xs_2@rt$corrected[[1]]) > 500)
    checkTrue(sum(xs_2@rt$raw[[3]] != xs_2@rt$corrected[[3]]) > 500)
    
    ## And the OnDiskMSnExp implementation:
    res <- xcms:::.obiwarp(od, param = prm)
    checkEquals(xs_2@rt$corrected, res)
    res_2 <- adjustRtime(od, param = prm)
    res_3 <- adjustRtime(xod, param = prm)
    checkEquals(adjustedRtime(res_3), res_2)
    checkEquals(lapply(adjustedRtime(res_3, bySample = TRUE), unname), res)
    checkEquals(adjustedRtime(res_3), res_2)
    ## Check if peaks were corrected correctly
    checkTrue(sum(chromPeaks(res_3)[, "rt"] == chromPeaks(xod)) <
              nrow(chromPeaks(res_3)) / 2)
    ## Dropping the adjusted rtime on these
    hasAdjustedRtime(res_3)
    tmp <- dropAdjustedRtime(res_3)
    checkEquals(chromPeaks(tmp), chromPeaks(xod))
    
    ## File issue on that! retcor.obiwarp does use round for the adjustment of
    ## the peak!
    ## -> issue #122
    ## checkEquals(chromPeaks(res_3), peaks(xs_2))
    
    ## Manually specify center Sample
    centerSample(prm) <- 3
    xs_2 <- retcor.obiwarp(xs, profStep = binSize(prm), center = centerSample(prm))
    checkEquals(xs_2@rt$raw[[centerSample(prm)]],
                xs_2@rt$corrected[[centerSample(prm)]])
    res <- xcms:::.obiwarp(od, param = prm)
    checkEquals(xs_2@rt$corrected, res)
    ## change some settings
    gapInit(prm) <- 3.1
    gapExtend(prm) <- 0.9
    xs_2 <- retcor.obiwarp(xs, profStep = binSize(prm), gapInit = gapInit(prm),
                           center = centerSample(prm), gapExtend = gapExtend(prm))
    checkEquals(xs_2@rt$raw[[centerSample(prm)]],
                xs_2@rt$corrected[[centerSample(prm)]])
    res <- xcms:::.obiwarp(od, param = prm)
    checkEquals(xs_2@rt$corrected, res)
}

## Run this test manually to perform an exhaustive test to validate obiwarp
## Results.
exhaustive_test <- function() {
    ## Load test files...
    faahko_files <- c(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
                      system.file('cdf/KO/ko16.CDF', package = "faahKO"),
                      system.file('cdf/KO/ko18.CDF', package = "faahKO"),
                      system.file('cdf/WT/wt15.CDF', package = "faahKO"),
                      system.file('cdf/WT/wt16.CDF', package = "faahKO"),
                      system.file('cdf/WT/wt18.CDF', package = "faahKO"))
    library(RUnit)
    library(xcms)
    ob <- readMSData(faahko_files, mode = "onDisk")
    xs <- xcmsSet(faahko_files, profparam = list(step = 0), method = "centWave",
                  noise = 10000, snthresh = 40)
    prm <- ObiwarpParam(binSize = 1, centerSample = 2)
    xs_r <- retcor.obiwarp(xs, profStep = binSize(prm))
    checkEquals(xs_r@rt$raw[[2]], xs_r@rt$corrected[[2]])
    res <- xcms:::.obiwarp(ob, param = prm)
    checkEquals(res, xs_r@rt$corrected)
    ## binSize
    binSize(prm) <- 0.2
    xs_r <- retcor.obiwarp(xs, profStep = binSize(prm))
    checkEquals(xs_r@rt$raw[[2]], xs_r@rt$corrected[[2]])
    res <- xcms:::.obiwarp(ob, param = prm)
    checkEquals(res, xs_r@rt$corrected)
    ## centersampe
    binSize(prm) <- 2
    centerSample(prm) <- 4
    xs_r <- retcor.obiwarp(xs, profStep = binSize(prm),
                           center = centerSample(prm))
    checkEquals(xs_r@rt$raw[[centerSample(prm)]],
                xs_r@rt$corrected[[centerSample(prm)]])
    res <- xcms:::.obiwarp(ob, param = prm)
    checkEquals(res, xs_r@rt$corrected)
    ## distFun
    distFun(prm) <- "euc"
    xs_r <- retcor.obiwarp(xs, profStep = binSize(prm),
                           center = centerSample(prm), distFunc = distFun(prm))
    checkEquals(xs_r@rt$raw[[centerSample(prm)]],
                xs_r@rt$corrected[[centerSample(prm)]])
    res <- xcms:::.obiwarp(ob, param = prm)
    checkEquals(res, xs_r@rt$corrected)
    ## localAlignment
    localAlignment(prm) <- TRUE
    distFun(prm) <- "cor"
    ## Uh huh! GET AN C ALLOCATION ERROR with local, stepsize 2 and euc
    xs_r <- retcor.obiwarp(xs, profStep = binSize(prm), localAlignment = 1,
                           center = centerSample(prm), distFunc = distFun(prm))
    checkEquals(xs_r@rt$raw[[centerSample(prm)]],
                xs_r@rt$corrected[[centerSample(prm)]])
    res <- xcms:::.obiwarp(ob, param = prm)
    checkEquals(res, xs_r@rt$corrected)
    ## factorDiag
    factorDiag(prm) <- 2.7
    localAlignment(prm) <- FALSE
    xs_r <- retcor.obiwarp(xs, profStep = binSize(prm), factorDiag = factorDiag(prm),
                           center = centerSample(prm), distFunc = distFun(prm))
    checkEquals(xs_r@rt$raw[[centerSample(prm)]],
                xs_r@rt$corrected[[centerSample(prm)]])
    res <- xcms:::.obiwarp(ob, param = prm)
    checkEquals(res, xs_r@rt$corrected)
    
    ## And all again using some of my own files.
    fls <- dir("/Users/jo/data/2016/2016-11/NoSN/", pattern = "mzML",
               full.names = TRUE)
    if (length(fls)) {
        fls <- fls[1:20]
        xs <- xcmsSet(fls, profparam = list(step = 0), method = "centWave",
                      noise = 10000, snthresh = 40)
        ## Compare also the timings!
        ## HM, why, with binSize 1.2 I get a "Dimension of profile matrices do
        ## not match"!
        prm <- ObiwarpParam(centerSample = 11, binSize = 1)
        system.time(
            xs_2 <- retcor.obiwarp(xs, profStep = binSize(prm), center = 11)
        )
        od <- readMSData(fls, mode = "onDisk")
        ## ???? dimension of profile matrix does not match???
        ## z <- filterFile(od, file = 1)
        ## cntr <- filterFile(od, file = centerSample(prm))
        ## cntrPr <- profMat(cntr, step = binSize(prm), returnBreaks = TRUE)[[1]]
        ## parms <- prm
        ##
        system.time(
            res <- xcms:::.obiwarp(od, param = prm)
        )
        checkEquals(res, xs_2@rt$corrected)
    }
}
