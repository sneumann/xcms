## Unit tests for all do_groupFeatures_* functions.

############################################################
## density
##

test_groupFeatures_FeatureDensityParam <- function() {
    od_x <- faahko_xod
    xs <- faahko_xs

    fdp <- FeatureDensityParam(sampleGroups = xs$class)
    od_x <- groupFeatures(od_x, param = fdp)
    xs <- group(xs, method = "density")
    checkEquals(xs@groupidx, featureGroups(od_x)$featureidx)
    fg <- featureGroups(od_x)
    fg <- as.matrix(fg[, -ncol(fg)])
    checkEquals(xs@groups, fg)
    checkTrue(length(processHistory(od_x)) == 2)
    ph <- processHistory(od_x, type = xcms:::.PROCSTEP.FEATURE.ALIGNMENT)[[1]]
    checkEquals(processParam(ph), fdp)

    fdp2 <- FeatureDensityParam(sampleGroups = xs$class, binSize = 2,
                                minFraction = 0.8)
    od_x <- groupFeatures(od_x, param = fdp2)
    xs <- group(xs, method = "density", minfrac = 0.8, mzwid = 2)
    checkEquals(xs@groupidx, featureGroups(od_x)$featureidx)
    fg <- featureGroups(od_x)
    fg <- as.matrix(fg[, -ncol(fg)])
    checkEquals(xs@groups, fg)
    checkTrue(length(processHistory(od_x)) == 2)
    ph <- processHistory(od_x, type = xcms:::.PROCSTEP.FEATURE.ALIGNMENT)[[1]]
    checkEquals(processParam(ph), fdp2)    
}


test_do_groupFeatures_density <- function() {
    fts <- peaks(faahko)
    res <- do_groupFeatures_density(fts, sampleGroups = sampclass(faahko))
    res_2 <- do_groupFeatures_density(fts, sampleGroups = sampclass(faahko),
                                      minFraction = 0.9)
    checkTrue(nrow(res$featureGroups) > nrow(res_2$featureGroups))
}

dontrun_do_groupFeatures_density_parallel <- function() {
    library(xcms)
    library(faahKO)
    library(RUnit)
    data(faahko)
    fts <- peaks(faahko)

    res <- xcms:::do_groupFeatures_density_par(fts, sampclass(faahko))
    res_2 <- do_groupFeatures_density(fts, sampclass(faahko))
    checkEquals(res$featureGroups, res_2$featureGroups)
    checkEquals(res$featureIndex, res_2$featureIndex)

    res <- xcms:::do_groupFeatures_density_par(fts, sampclass(faahko), bw = 10)
    res_2 <- do_groupFeatures_density(fts, sampclass(faahko), bw = 10)
    checkEquals(res$featureGroups, res_2$featureGroups)
    checkEquals(res$featureIndex, res_2$featureIndex)

    res <- xcms:::do_groupFeatures_density_par(fts, sampclass(faahko),
                                               minFraction = 0.9)
    res_2 <- do_groupFeatures_density(fts, sampclass(faahko), minFraction = 0.9)
    checkEquals(res$featureGroups, res_2$featureGroups)
    checkEquals(res$featureIndex, res_2$featureIndex)
}

## This is to ensure that the do_groupFeatures_density yields identical results
## than the group.density method. Once we're sure of that we rename this to
## "dontrun" and replace the code within the group.density method.
dontrun_do_groupFeatures_density_compare <- function() {
    xs <- faahko
    
    fts <- peaks(xs)
    ## o Default
    smpGrp <- sampclass(xs)
    theBw <- 30
    minFr <- 0.5
    minSmp <- 1
    mzW <- 0.25
    maxFts <- 50
    ftGrps <- xcms:::do_groupFeatures_density(features = fts,
                                              sampleGroups = smpGrp,
                                              bw = theBw,
                                              minFraction = minFr,
                                              minSamples = minSmp,
                                              binSize = mzW,
                                              maxFeatures = maxFts)
    xs$class <- smpGrp
    res_x <- group.density(xs, bw = theBw, minfrac = minFr,
                           minsamp = minSmp, mzwid = mzW, max = maxFts)
    checkEquals(res_x@groups, ftGrps$featureGroups)
    checkEquals(res_x@groupidx, ftGrps$featureIndex)
    ## o All one group, different bw
    smpGrp <- rep(1, length(filepaths(xs)))
    theBw <- 10
    minFr <- 0.5
    minSmp <- 1
    mzW <- 0.25
    maxFts <- 50
    ftGrps <- xcms:::do_groupFeatures_density(features = fts,
                                              sampleGroups = smpGrp,
                                              bw = theBw,
                                              minFraction = minFr,
                                              minSamples = minSmp,
                                              binSize = mzW,
                                              maxFeatures = maxFts)
    xs$class <- smpGrp
    res_x <- group.density(xs, bw = theBw, minfrac = minFr,
                           minsamp = minSmp, mzwid = mzW, max = maxFts)
    checkEquals(res_x@groups, ftGrps$featureGroups)
    checkEquals(res_x@groupidx, ftGrps$featureIndex)
    ## o Three groups, minfrac
    smpGrp <- c(1, 1, 2, 2, 3, 3, 3, 3, 2, 2, 3, 3)
    theBw <- 30
    minFr <- 0.3
    minSmp <- 1
    mzW <- 0.4
    maxFts <- 50
    ftGrps <- xcms:::do_groupFeatures_density(features = fts,
                                              sampleGroups = smpGrp,
                                              bw = theBw,
                                              minFraction = minFr,
                                              minSamples = minSmp,
                                              binSize = mzW,
                                              maxFeatures = maxFts)
    xs$class <- smpGrp
    res_x <- group.density(xs, bw = theBw, minfrac = minFr,
                           minsamp = minSmp, mzwid = mzW, max = maxFts)
    checkEquals(res_x@groups, ftGrps$featureGroups)
    checkEquals(res_x@groupidx, ftGrps$featureIndex)
    ## o change also minSmp and maxFts
    smpGrp <- c(1, 1, 2, 2, 3, 3, 1, 1, 2, 2, 3, 3)
    theBw <- 30
    minFr <- 0.3
    minSmp <- 2
    mzW <- 0.4
    maxFts <- 10
    ftGrps <- xcms:::do_groupFeatures_density(features = fts,
                                              sampleGroups = smpGrp,
                                              bw = theBw,
                                              minFraction = minFr,
                                              minSamples = minSmp,
                                              binSize = mzW,
                                              maxFeatures = maxFts)
    xs$class <- smpGrp
    res_x <- group.density(xs, bw = theBw, minfrac = minFr,
                           minsamp = minSmp, mzwid = mzW, max = maxFts)
    checkEquals(res_x@groups, ftGrps$featureGroups)
    checkEquals(res_x@groupidx, ftGrps$featureIndex)
}


dontrun_groupFeatures_density_implementation <- function() {
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

## This is to compare the function to the group.mzClust method. Once all is fine
## rename it to "dontrun"
test_groupFeatures_mzClust_compare <- function() {
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
    od <- readMSData2(mzdatafiles, msLevel. = 1)
    p <- MSWParam(scales = c(1, 7), ampTh = 0.005, peakThr = 80000,
                  SNR.method = 'data.mean', winSize.noise = 500)
    xod <- detectFeatures(od, param = p)
    res <- do_groupFeatures_mzClust(features(xod), sampleGroups = sampclass(xs))
    
    checkEquals(peaks(xs), features(xod))
    checkEquals(res$featureGroups, xsg@groups)
    checkEquals(res$featureIndex, xsg@groupidx)
    
    ## Check with different class ordering!
    sc_orig <- sampclass(xs)
    sc <- c(2, 2, 1, 1, 4, 4, 4, 3, 3, 3)
    sampclass(xs) <- factor(sc)
    xsg <- group(xs, method="mzClust", minfrac = 0.2)
    res <- do_groupFeatures_mzClust(features(xod), sampleGroups = sc,
                                    minFraction = 0.2)
    checkEquals(res$featureGroups, xsg@groups)
    checkEquals(res$featureIndex, xsg@groupidx)

    sc <- c("z", "z", "b", "a", "a", "a", "z", "b", "e", "e")
    sampclass(xs) <- factor(sc)
    xsg <- group(xs, method="mzClust", minfrac = 0.2, mzppm = 40)
    res <- do_groupFeatures_mzClust(features(xod), sampleGroups = sc,
                                    minFraction = 0.1, ppm = 40)
    checkEquals(res$featureGroups, xsg@groups)
    checkEquals(res$featureIndex, xsg@groupidx)

    xsg <- group(xs, method="mzClust", minfrac = 0.2, mzppm = 40, mzabs = 0.1)
    res <- do_groupFeatures_mzClust(features(xod), sampleGroups = sc,
                                    minFraction = 0.1, ppm = 40, absMz = 0.1)
    checkEquals(res$featureGroups, xsg@groups)
    checkEquals(res$featureIndex, xsg@groupidx)
}
