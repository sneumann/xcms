## Unit tests for all do_groupChromPeaks_* functions and methods/functions related
## to feature grouping.

## General functions/methods
test_featureValues_XCMSnExp <- function() {
    od_x <- faahko_xod    
    xs <- faahko_xs
    
    p <- PeakDensityParam(sampleGroups = xs$class)
    od_x <- groupChromPeaks(od_x, param = p)

    xs <- group(xs, method = "density")

    checkEquals(unname(groupval(xs, value = "into")),
                unname(featureValues(od_x, value = "into")))
    checkEquals(unname(groupval(xs, method = "maxint", value = "into")),
                unname(featureValues(od_x, method = "maxint", value = "into")))
    ## Checking errors
    checkException(featureValues(od_x, value = "bla"))
    
}

############################################################
## density
##

test_groupChromPeaks_PeakDensityParam <- function() {
    od_x <- faahko_xod
    xs <- faahko_xs

    fdp <- PeakDensityParam(sampleGroups = xs$class)
    od_x <- groupChromPeaks(od_x, param = fdp)
    xs <- group(xs, method = "density")
    checkEquals(xs@groupidx, featureDefinitions(od_x)$peakidx)
    fg <- featureDefinitions(od_x)
    fg <- S4Vectors::as.matrix(fg[, -ncol(fg)])
    rownames(fg) <- NULL
    checkEquals(xs@groups, fg)
    checkTrue(length(processHistory(od_x)) == 2)
    ph <- processHistory(od_x, type = xcms:::.PROCSTEP.PEAK.GROUPING)[[1]]
    checkEquals(processParam(ph), fdp)
    checkEquals(rownames(featureDefinitions(od_x)),
                xcms:::.featureIDs(nrow(featureDefinitions(od_x))))
    
    fdp2 <- PeakDensityParam(sampleGroups = xs$class, binSize = 2,
                                minFraction = 0.8)
    od_x <- groupChromPeaks(od_x, param = fdp2)
    xs <- group(xs, method = "density", minfrac = 0.8, mzwid = 2)
    checkEquals(xs@groupidx, featureDefinitions(od_x)$peakidx)
    fg <- featureDefinitions(od_x)
    fg <- S4Vectors::as.matrix(fg[, -ncol(fg)])
    rownames(fg) <- NULL
    checkEquals(xs@groups, fg)
    checkTrue(length(processHistory(od_x)) == 2)
    ph <- processHistory(od_x, type = xcms:::.PROCSTEP.PEAK.GROUPING)[[1]]
    checkEquals(processParam(ph), fdp2)
    checkEquals(rownames(featureDefinitions(od_x)),
                xcms:::.featureIDs(nrow(featureDefinitions(od_x))))
}

test_do_groupChromPeaks_density <- function() {
    fts <- peaks(faahko)
    res <- do_groupChromPeaks_density(fts, sampleGroups = sampclass(faahko))
    res_2 <- do_groupChromPeaks_density(fts, sampleGroups = sampclass(faahko),
                                      minFraction = 0.9)
    checkTrue(nrow(res$featureDefinitions) > nrow(res_2$featureDefinitions))
}

dontrun_do_groupChromPeaks_density_parallel <- function() {
    library(xcms)
    library(faahKO)
    library(RUnit)
    data(faahko)
    fts <- peaks(faahko)

    res <- xcms:::do_groupChromPeaks_density_par(fts, sampclass(faahko))
    res_2 <- do_groupChromPeaks_density(fts, sampclass(faahko))
    checkEquals(res$featureDefinitions, res_2$featureDefinitions)
    checkEquals(res$peakIndex, res_2$peakIndex)

    res <- xcms:::do_groupChromPeaks_density_par(fts, sampclass(faahko), bw = 10)
    res_2 <- do_groupChromPeaks_density(fts, sampclass(faahko), bw = 10)
    checkEquals(res$featureDefinitions, res_2$featureDefinitions)
    checkEquals(res$peakIndex, res_2$peakIndex)

    res <- xcms:::do_groupChromPeaks_density_par(fts, sampclass(faahko),
                                               minFraction = 0.9)
    res_2 <- do_groupChromPeaks_density(fts, sampclass(faahko), minFraction = 0.9)
    checkEquals(res$featureDefinitions, res_2$featureDefinitions)
    checkEquals(res$peakIndex, res_2$peakIndex)
}

## This is to ensure that the do_groupChromPeaks_density yields identical results
## than the group.density method. Once we're sure of that we rename this to
## "dontrun" and replace the code within the group.density method.
dontrun_do_groupChromPeaks_density_compare <- function() {
    xs <- faahko
    
    fts <- peaks(xs)
    ## o Default
    smpGrp <- sampclass(xs)
    theBw <- 30
    minFr <- 0.5
    minSmp <- 1
    mzW <- 0.25
    maxFts <- 50
    ftGrps <- xcms:::do_groupChromPeaks_density(peaks = fts,
                                              sampleGroups = smpGrp,
                                              bw = theBw,
                                              minFraction = minFr,
                                              minSamples = minSmp,
                                              binSize = mzW,
                                              maxFeatures = maxFts)
    xs$class <- smpGrp
    res_x <- group.density(xs, bw = theBw, minfrac = minFr,
                           minsamp = minSmp, mzwid = mzW, max = maxFts)
    checkEquals(res_x@groups, ftGrps$featureDefinitions)
    checkEquals(res_x@groupidx, ftGrps$peakIndex)
    ## o All one group, different bw
    smpGrp <- rep(1, length(filepaths(xs)))
    theBw <- 10
    minFr <- 0.5
    minSmp <- 1
    mzW <- 0.25
    maxFts <- 50
    ftGrps <- xcms:::do_groupChromPeaks_density(peaks = fts,
                                              sampleGroups = smpGrp,
                                              bw = theBw,
                                              minFraction = minFr,
                                              minSamples = minSmp,
                                              binSize = mzW,
                                              maxFeatures = maxFts)
    xs$class <- smpGrp
    res_x <- group.density(xs, bw = theBw, minfrac = minFr,
                           minsamp = minSmp, mzwid = mzW, max = maxFts)
    checkEquals(res_x@groups, ftGrps$featureDefinitions)
    checkEquals(res_x@groupidx, ftGrps$peakIndex)
    ## o Three groups, minfrac
    smpGrp <- c(1, 1, 2, 2, 3, 3, 3, 3, 2, 2, 3, 3)
    theBw <- 30
    minFr <- 0.3
    minSmp <- 1
    mzW <- 0.4
    maxFts <- 50
    ftGrps <- xcms:::do_groupChromPeaks_density(peaks = fts,
                                              sampleGroups = smpGrp,
                                              bw = theBw,
                                              minFraction = minFr,
                                              minSamples = minSmp,
                                              binSize = mzW,
                                              maxFeatures = maxFts)
    xs$class <- smpGrp
    res_x <- group.density(xs, bw = theBw, minfrac = minFr,
                           minsamp = minSmp, mzwid = mzW, max = maxFts)
    checkEquals(res_x@groups, ftGrps$featureDefinitions)
    checkEquals(res_x@groupidx, ftGrps$peakIndex)
    ## o change also minSmp and maxFts
    smpGrp <- c(1, 1, 2, 2, 3, 3, 1, 1, 2, 2, 3, 3)
    theBw <- 30
    minFr <- 0.3
    minSmp <- 2
    mzW <- 0.4
    maxFts <- 10
    ftGrps <- xcms:::do_groupChromPeaks_density(peaks = fts,
                                              sampleGroups = smpGrp,
                                              bw = theBw,
                                              minFraction = minFr,
                                              minSamples = minSmp,
                                              binSize = mzW,
                                              maxFeatures = maxFts)
    xs$class <- smpGrp
    res_x <- group.density(xs, bw = theBw, minfrac = minFr,
                           minsamp = minSmp, mzwid = mzW, max = maxFts)
    checkEquals(res_x@groups, ftGrps$featureDefinitions)
    checkEquals(res_x@groupidx, ftGrps$peakIndex)
}


dontrun_groupChromPeaks_density_implementation <- function() {
    library(faahKO)
    data("faahko")
    ## 1) check whether the bins are really overlapping.
    pks <- peaks(faahko)
    pks <- pks[order(pks[, "mz"]), ]
    mzWid <- 2
    mass <- seq(pks[1, "mz"], pks[nrow(pks), "mz"] + mzWid, by = mzWid / 2)
    ## masspos is the index in pks for the feature with an mz value >= mass
    masspos <- xcms:::findEqualGreaterM(pks[, "mz"], mass)
    ## My approach:
    idx <- findInterval(pks[, "mz"], mass)
    ## That's not working with overlapping though.
    pks_my <- split.data.frame(pks, f = idx)

    ## Test 2:
    didx <- diff(idx)
    nds <- c(which(didx > 0), nrow(pks))
    strt <- c(1, nds[-length(nds)] + 1)
    ## Fix the starts, if diff is equals to one include that too.
    pks_2 <- mapply(strt, nds, FUN = function(a, b) {
        return(pks[a:b, , drop = FALSE])
    })
    library(RUnit)
    checkEquals(unname(pks_2), unname(pks_my))

    ## Idea is to add also the data from the previous bin, if its only an index
    ## of one away
    tmp <- pks_my
    nameNum <- as.numeric(names(pks_my))
    for (i in 2:length(tmp)) {
        if ((nameNum[i] - nameNum[i - 1]) == 1)
            tmp[[i]] <- rbind(pks_my[[i-1]], pks_my[[i]])
    }

    ## The overlapping mz ranges is tricky - eventually just stick with original
    ## code.

    ## xcms:
    res_x <- vector("list", length(mass))
    idxs <- res_x
    for (i in seq(length = length(mass) - 2)) {
        startidx <- masspos[i]
        endidx <- masspos[i + 2] - 1
        if (endidx - startidx < 0)
            next
        idxs[[i]] <- c(startidx, endidx)
        res_x[[i]] <- pks[startidx:endidx, , drop = FALSE]
    }
    res_x <- res_x[lengths(res_x) > 0]
    idxs <- idxs[lengths(idxs) > 0]
}


############################################################
## mzClust
##
## library(msdata)
## fticrf <- list.files(system.file("fticr", package = "msdata"),
##                      recursive = TRUE, full.names = TRUE)

## ## old
## ## new
## fticr_od <- readMSData(fticrf[1:2], msLevel. = 1, mode = "onDisk")
## p <- MSWParam(scales = c(1, 7), peakThr = 80000, ampTh = 0.005,
##               SNR.method = "data.mean", winSize.noise = 500)
## fticr_xod <- findChromPeaks(fticr_od, param = p)

test_do_groupPeaks_mzClust <- function() {
    fts <- peaks(fticr_xs)
    res <- do_groupPeaks_mzClust(peaks = fts,
                                 sampleGroups = sampclass(fticr_xs))
    res_2 <- do_groupPeaks_mzClust(peaks = fts,
                                   sampleGroups = sampclass(fticr_xs),
                                   minFraction = 0, absMz = 2)
    checkTrue(nrow(res$featureDefinitions) > nrow(res_2$featureDefinitions))

    res_x <- group(fticr_xs, method = "mzClust")
    checkEquals(res_x@groups, res$featureDefinitions)
    checkEquals(res_x@groupidx, res$peakIndex)
}

## This is to compare the function to the group.mzClust method. Once all is fine
## rename it to "dontrun"
dontrun_test_groupPeaks_mzClust_compare <- function() {
    library(RUnit)
    library(xcms)
    library(msdata)
    mzdatapath <- system.file("fticr", package = "msdata")
    mzdatafiles <- list.files(mzdatapath, recursive = TRUE, full.names = TRUE)

    ## old
    xs <- xcmsSet(method="MSW", files=mzdatafiles, scales=c(1,7),
                  SNR.method='data.mean' , winSize.noise=500,
                  peakThr=80000,  amp.Th=0.005)
    xsg <- group(xs, method="mzClust")

    ## new
    od <- readMSData(mzdatafiles, msLevel. = 1, mode = "onDisk")
    p <- MSWParam(scales = c(1, 7), ampTh = 0.005, peakThr = 80000,
                  SNR.method = 'data.mean', winSize.noise = 500)
    xod <- findChromPeaks(od, param = p)
    res <- do_groupPeaks_mzClust(chromPeaks(xod), sampleGroups = sampclass(xs))
    
    checkEquals(peaks(xs), chromPeaks(xod))
    checkEquals(res$featureDefinitions, xsg@groups)
    checkEquals(res$peakIndex, xsg@groupidx)
    
    ## Check with different class ordering!
    sc_orig <- sampclass(xs)
    sc <- c(2, 2, 1, 1, 4, 4, 4, 3, 3, 3)
    sampclass(xs) <- factor(sc)
    xsg <- group(xs, method="mzClust", minfrac = 0.2)
    res <- do_groupPeaks_mzClust(chromPeaks(xod), sampleGroups = sc,
                                    minFraction = 0.2)
    checkEquals(res$featureDefinitions, xsg@groups)
    checkEquals(res$peakIndex, xsg@groupidx)

    sc <- c("z", "z", "b", "a", "a", "a", "z", "b", "e", "e")
    sampclass(xs) <- factor(sc)
    xsg <- group(xs, method="mzClust", minfrac = 0.2, mzppm = 40)
    res <- do_groupPeaks_mzClust(chromPeaks(xod), sampleGroups = sc,
                                    minFraction = 0.1, ppm = 40)
    checkEquals(res$featureDefinitions, xsg@groups)
    checkEquals(res$peakIndex, xsg@groupidx)

    xsg <- group(xs, method="mzClust", minfrac = 0.2, mzppm = 40, mzabs = 0.1)
    res <- do_groupPeaks_mzClust(chromPeaks(xod), sampleGroups = sc,
                                    minFraction = 0.1, ppm = 40, absMz = 0.1)
    checkEquals(res$featureDefinitions, xsg@groups)
    checkEquals(res$peakIndex, xsg@groupidx)
}

test_groupPeaks_MzClustParam <- function() {

    p <- MzClustParam(sampleGroups = sampclass(fticr_xs))
    
    fticr_xod2 <- groupChromPeaks(fticr_xod, param = p)
    fticr_xs2 <- group(fticr_xs, method = "mzClust")
    checkEquals(fticr_xs2@groupidx, featureDefinitions(fticr_xod2)$peakidx)
    fg <- featureDefinitions(fticr_xod2)
    fg <- S4Vectors::as.matrix(fg[, -ncol(fg)])
    rownames(fg) <- NULL
    checkEquals(fticr_xs2@groups, fg)
    checkTrue(length(processHistory(fticr_xod2)) == 2)
    ph <- processHistory(fticr_xod2,
                         type = xcms:::.PROCSTEP.PEAK.GROUPING)[[1]]
    checkEquals(processParam(ph), p)
    checkEquals(rownames(featureDefinitions(fticr_xod2)),
                xcms:::.featureIDs(nrow(featureDefinitions(fticr_xod2))))
    
    p2 <- MzClustParam(sampleGroups = fticr_xs$class, absMz = 1,
                       minFraction = 0.8)
    fticr_xod2 <- groupChromPeaks(fticr_xod, param = p2)
    fticr_xs2 <- group(fticr_xs, method = "mzClust", minfrac = 0.8, mzabs = 1)
    checkEquals(fticr_xs2@groupidx, featureDefinitions(fticr_xod2)$peakidx)
    fg <- featureDefinitions(fticr_xod2)
    fg <- S4Vectors::as.matrix(fg[, -ncol(fg)])
    rownames(fg) <- NULL
    checkEquals(fticr_xs2@groups, fg)
    checkTrue(length(processHistory(fticr_xod2)) == 2)
    ph <- processHistory(fticr_xod2,
                         type = xcms:::.PROCSTEP.PEAK.GROUPING)[[1]]
    checkEquals(processParam(ph), p2)    
    checkEquals(rownames(featureDefinitions(fticr_xod2)),
                xcms:::.featureIDs(nrow(featureDefinitions(fticr_xod2))))
}

############################################################
## nearest
##

test_do_groupChromPeaks_nearest <- function() {
    ## xs <- faahko
    xs <- faahko_xs
    features <- peaks(xs)
    sampleGroups <- sampclass(xs)
    mzVsRtBalance <- 10
    mzCheck <- 0.2
    rtCheck <- 15
    kNN <- 10

    res <- do_groupChromPeaks_nearest(features, sampleGroups)
    res_2 <- do_groupChromPeaks_nearest(features, sampleGroups, absRt = 3)
    checkTrue(nrow(res$featureDefinitions) < nrow(res_2$featureDefinitions))
    res_x <- group(xs, method = "nearest")
    checkEquals(res_x@groups, res$featureDefinitions)
}

test_groupChromPeaks_NearestPeaksParam <- function() {
    od_x <- faahko_xod
    xs <- faahko_xs

    p <- NearestPeaksParam(sampleGroups = xs$class)
    od_x <- groupChromPeaks(od_x, param = p)
    xs <- group(xs, method = "nearest")
    checkEquals(xs@groupidx, featureDefinitions(od_x)$peakidx)
    fg <- featureDefinitions(od_x)
    fg <- S4Vectors::as.matrix(fg[, -ncol(fg)])
    rownames(fg) <- NULL
    checkEquals(xs@groups, fg)
    checkTrue(length(processHistory(od_x)) == 2)
    ph <- processHistory(od_x, type = xcms:::.PROCSTEP.PEAK.GROUPING)[[1]]
    checkEquals(processParam(ph), p)
    checkEquals(rownames(featureDefinitions(od_x)),
                xcms:::.featureIDs(nrow(featureDefinitions(od_x))))
    
    fdp2 <- NearestPeaksParam(sampleGroups = xs$class, kNN = 3)
    od_x <- groupChromPeaks(od_x, param = fdp2)
    xs <- group(xs, method = "nearest", kNN = 3)
    checkEquals(xs@groupidx, featureDefinitions(od_x)$peakidx)
    fg <- featureDefinitions(od_x)
    fg <- S4Vectors::as.matrix(fg[, -ncol(fg)])
    rownames(fg) <- NULL
    checkEquals(xs@groups, fg)
    checkTrue(length(processHistory(od_x)) == 2)
    ph <- processHistory(od_x, type = xcms:::.PROCSTEP.PEAK.GROUPING)[[1]]
    checkEquals(processParam(ph), fdp2)    
    checkEquals(rownames(featureDefinitions(od_x)),
                xcms:::.featureIDs(nrow(featureDefinitions(od_x))))
}


## That's to ensure that the do_ function yields identical results than the
## group.nearest method. Once we've replaced the code in the latter we rename
## this function to "dontrun".
dontrun_test_nearest_impl <- function() {
    library(RUnit)
    library(xcms)
    library(faahKO)
    data(faahko)
    xs <- faahko
    features <- peaks(xs)
    sampleGroups <- sampclass(xs)
    mzVsRtBalance <- 10
    absMz <- 0.2
    absRt <- 15
    kNN <- 10

    res <- group(xs, method = "nearest") ## Oh, nasty warnings! These were
    ## already there in version 1.51.0!
    ## 1.48.0: Yup.    
    res_2 <- do_groupChromPeaks_nearest(features, sampleGroups)
    res_n <- xcms:::do_groupChromPeaks_nearest_mod(features, sampleGroups)
    
    checkEquals(res@groups, res_2$featureDefinitions)
    checkEquals(res@groupidx, res_2$peakIndex)
    checkEquals(res_n$peakIndex, res_2$peakIndex)
    checkEquals(res_n$featureDefinitions, res_2$featureDefinitions)

    ## change sample grouping.
    sc <- c("b", "b", "a", "a", "z", "z", "a", "b", "e", "e", "e", "e")
    sampclass(xs) <- sc
    ## levels are NOT ordered.
    res <- group(xs, method = "nearest")
    res_2 <- do_groupChromPeaks_nearest(features, sc)
    checkEquals(res@groups, res_2$featureDefinitions)
    checkEquals(res@groupidx, res_2$peakIndex)
    res_n <- xcms:::do_groupChromPeaks_nearest_mod(features, sc)
    checkEquals(res_n$peakIndex, res_2$peakIndex)
    checkEquals(res_n$featureDefinitions, res_2$featureDefinitions)

    ## Use sample assignment with ordered levels.
    sampclass(xs) <- factor(sc) ## this re-orders levels!
    res_3 <- group(xs, method = "nearest")
    res_4 <- do_groupChromPeaks_nearest(features, factor(sc))
    checkEquals(res_3@groups, res@groups)
    ## Now, the do functions does NOT re-order levels!
    checkEquals(res@groups[, 1:7], res_4$featureDefinitions[, 1:7])
    checkEquals(res@groups[, levels(factor(sc))],
                res_4$featureDefinitions[, levels(factor(sc))])
    checkEquals(res@groupidx, res_4$peakIndex)
    res_n <- xcms:::do_groupChromPeaks_nearest_mod(features, factor(sc))
    checkEquals(res_n$peakIndex, res_4$peakIndex)
    checkEquals(res_n$featureDefinitions, res_4$featureDefinitions)
    
    ## Now change settings.
    sampclass(xs) <- sc
    res <- group(xs, method = "nearest", mzVsRTbalance = 5)
    res_2 <- do_groupChromPeaks_nearest(features, sc, mzVsRtBalance = 5)
    checkEquals(res@groups, res_2$featureDefinitions)
    checkEquals(res@groupidx, res_2$peakIndex)
    res_n <- xcms:::do_groupChromPeaks_nearest_mod(features, sc, mzVsRtBalance = 5)
    checkEquals(res_n$peakIndex, res_2$peakIndex)
    checkEquals(res_n$featureDefinitions, res_2$featureDefinitions)

    res <- group(xs, method = "nearest", kNN = 3)
    res_2 <- do_groupChromPeaks_nearest(features, sc, kNN = 3)
    checkEquals(res@groups, res_2$featureDefinitions)
    checkEquals(res@groupidx, res_2$peakIndex)
    res_n <- xcms:::do_groupChromPeaks_nearest_mod(features, sc, kNN)
    checkEquals(res_n$peakIndex, res_2$peakIndex)
    checkEquals(res_n$featureDefinitions, res_2$featureDefinitions)

    res <- group(xs, method = "nearest", mzCheck = 0.5)
    res_2 <- do_groupChromPeaks_nearest(features, sc, absMz = 0.5)
    checkEquals(res@groups, res_2$featureDefinitions)
    checkEquals(res@groupidx, res_2$peakIndex)
    res_n <- xcms:::do_groupChromPeaks_nearest_mod(features, sc, absMz = 0.5)
    checkEquals(res_n$peakIndex, res_2$peakIndex)
    checkEquals(res_n$featureDefinitions, res_2$featureDefinitions)

    res <- group(xs, method = "nearest", rtCheck = 3)
    res_2 <- do_groupChromPeaks_nearest(features, sc, absRt = 3)
    checkEquals(res@groups, res_2$featureDefinitions)
    checkEquals(res@groupidx, res_2$peakIndex)
    res_n <- xcms:::do_groupChromPeaks_nearest_mod(features, sc, absRt = 3)
    checkEquals(res_n$peakIndex, res_2$peakIndex)
    checkEquals(res_n$featureDefinitions, res_2$featureDefinitions)

    ## library(profvis)
    ## profvis({
    ##     res_2 <- xcms:::do_groupChromPeaks_nearest(features, sampleGroups)
    ## })
}
