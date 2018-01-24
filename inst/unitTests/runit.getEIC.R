## Testing the profEIC as well as the getEIC method.
test_profEIC <- function() {
    ## file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    ## Compare the results with manual calculations on the profile matrix.
    step <- 1
    ## xraw <- xcmsRaw(file, profstep = 0)
    xraw <- deepCopy(faahko_xr_1)
    profmat <- profMat(xraw, step = step)
    mzr <- c(200, 201)
    rtr <- c(3000, 3500)
    res_1 <- xcms:::profEIC(xraw, mzrange = mzr, rtrange = rtr, step = step)
    res_x <- getEIC(xraw, mzrange = mzr, rtrange = rtr, step = step)
    ## manual calculation:
    scn_idx <- which(xraw@scantime >= rtr[1] & xraw@scantime <= rtr[2])
    mass <- seq(floor(xraw@mzrange[1] / step) * step,
                ceiling(xraw@mzrange[2] / step) * step, by = step)
    mass_idx <- which(mass >= mzr[1] & mass <= mzr[2])
    max_mass <- apply(profmat[mass_idx, scn_idx], MARGIN = 2, max)
    rts <- xraw@scantime[scn_idx]
    checkEquals(res_1@eic[[1]][[1]][, 1], rts)
    checkEquals(res_1@eic[[1]][[1]][, 2], max_mass)
    checkEquals(res_1, res_x)

    ## OK, now with binlinbase
    profMethod(xraw) <- "binlinbase"
    xraw@profparam <- list(baselevel = 666666)
    profmat <- profMat(xraw, step = 1)
    res_2 <- xcms:::profEIC(xraw, mzrange = mzr, rtrange = rtr, step = step)
    res_x <- getEIC(xraw, mzrange = mzr, rtrange = rtr, step = step)
    max_mass <- apply(profmat[mass_idx, scn_idx], MARGIN = 2, max)
    checkEquals(res_2@eic[[1]][[1]][, 1], rts)
    checkEquals(res_2@eic[[1]][[1]][, 2], max_mass)
    checkTrue(sum(res_2@eic[[1]][[1]][, 2] == 666666) > 0)
    checkEquals(res_2, res_x)
    ## Results should be different to previous ones.
    checkTrue(any(res_2@eic[[1]][[1]][, 2] != res_1@eic[[1]][[1]][, 2]))
    checkTrue(all(res_2@eic[[1]][[1]][, 1] == res_1@eic[[1]][[1]][, 1]))
    ## Without pre-calculated profile matrix:
    profmat <- profMat(xraw, step = 0.2)
    res_3 <- xcms:::profEIC(xraw, mzrange = mzr, rtrange = rtr, step = 0.2)
    res_x <- getEIC(xraw, mzrange = mzr, rtrange = rtr, step = 0.2)
    max_mass <- apply(profmat[mass_idx, scn_idx], MARGIN = 2, max)
    checkEquals(res_3@eic[[1]][[1]][, 1], rts)
    checkEquals(res_3@eic[[1]][[1]][, 2], max_mass)
    checkEquals(res_3, res_x)
    ## Results should be different to previous ones.
    checkTrue(any(res_3@eic[[1]][[1]][, 2] != res_2@eic[[1]][[1]][, 2]))
    checkTrue(all(res_3@eic[[1]][[1]][, 1] == res_2@eic[[1]][[1]][, 1]))

    ## Now with intlin
    profMethod(xraw) <- "intlin"
    profmat <- profMat(xraw, step = 1)
    res_3 <- xcms:::profEIC(xraw, mzrange = mzr, rtrange = rtr, step = step)
    res_x <- getEIC(xraw, mzrange = mzr, rtrange = rtr, step = step)
    max_mass <- apply(profmat[mass_idx, scn_idx], MARGIN = 2, max)
    checkEquals(res_3@eic[[1]][[1]][, 1], rts)
    checkEquals(res_3@eic[[1]][[1]][, 2], max_mass)
    checkEquals(res_3, res_x)

    ## Now test with multiple mz ranges
    mzrm <- rbind(c(300, 302), mzr, c(205, 207), c(500, 507))
    profMethod(xraw) <- "bin"
    profmat <- profMat(xraw, step = 1)
    max_mass <- apply(profmat[mass_idx, scn_idx], MARGIN = 2, max)
    res_1 <- xcms:::profEIC(xraw, mzrange = mzrm, rtrange = rtr, step = step)
    res_x <- getEIC(xraw, mzrange = mzrm, rtrange = rtr, step = step)
    checkEquals(res_1@eic[[1]][[2]][, 2], max_mass)
    ## 4th range:
    mass_idx_2 <- which(mass >= 500 & mass <= 507)
    max_mass <- apply(profmat[mass_idx_2, scn_idx], MARGIN = 2, max)
    checkEquals(res_1@eic[[1]][[4]][, 2], max_mass)
    checkEquals(res_1, res_x)

    ## Different step:
    step <- 0.2
    res_1 <- xcms:::profEIC(xraw, mzrange = mzrm, rtrange = rtr, step = step)
    res_x <- getEIC(xraw, mzrange = mzrm, rtrange = rtr, step = step)
    mass <- seq(floor(xraw@mzrange[1] / step) * step,
                ceiling(xraw@mzrange[2] / step) * step, by = step)
    mass_idx <- which(mass >= mzr[1] & mass <= mzr[2])
    scn_idx <- which(xraw@scantime >= rtr[1] & xraw@scantime <= rtr[2])
    profmat <- profMat(xraw, step = step)
    max_mass <- apply(profmat[mass_idx, scn_idx], MARGIN = 2, max)
    checkEquals(res_1@eic[[1]][[2]][, 2], max_mass)
    checkEquals(res_1, res_x)
    ## 4th range:
    mass_idx_2 <- which(mass >= 500 & mass <= 507)
    max_mass <- apply(profmat[mass_idx_2, scn_idx], MARGIN = 2, max)
    checkEquals(res_1@eic[[1]][[4]][, 2], max_mass)

    ## Test with the whole mzrange.
    step <- 0.5
    profMethod(xraw) <- "bin"
    profmat <- profMat(xraw, step = step)
    res_1 <- xcms:::profEIC(xraw, rtrange = rtr, step = step)
    res_x <- getEIC(xraw, rtrange = rtr, step = step)
    scn_idx <- which(xraw@scantime >= rtr[1] & xraw@scantime <= rtr[2])
    max_mass <- apply(profmat[, scn_idx], MARGIN = 2, max)
    checkEquals(res_1@eic[[1]][[1]][, 2], max_mass) ## Compare max per rt
    checkEquals(res_1@eic[[1]][[1]][, 1], xraw@scantime[scn_idx])
    checkEquals(res_1, res_x)

    ## Whole rtrange:
    res_1 <- xcms:::profEIC(xraw, mzrange = mzr, step = step)
    res_x <- getEIC(xraw, mzrange = mzr, step = step)
    mass <- seq(floor(xraw@mzrange[1] / step) * step,
                ceiling(xraw@mzrange[2] / step) * step, by = step)
    mass_idx <- which(mass >= mzr[1] & mass <= mzr[2])
    max_mass <- apply(profmat[mass_idx, ], MARGIN = 2, max)
    checkEquals(res_1@eic[[1]][[1]][, 2], max_mass)
    checkEquals(res_1@eic[[1]][[1]][, 1], xraw@scantime)
    checkEquals(res_1, res_x)

    ## Everything.
    res_1 <- xcms:::profEIC(xraw, step = step)
    res_x <- getEIC(xraw, step = step)
    max_mass <- apply(profmat, MARGIN = 2, max)
    checkEquals(res_1@eic[[1]][[1]][, 2], max_mass)
    checkEquals(res_1, res_x)

    ## Error checking:
    ## mz range outside
    mzr <- c(0, 14)
    checkException(xcms:::profEIC(xraw, mzrange = mzr))
    checkException(getEIC(xraw, mzrange = mzr))
    ## rt range outside
    checkException(xcms:::profEIC(xraw, rtrange = c(6000, 6060)))
    checkException(getEIC(xraw, rtrange = c(6000, 6060)))
    ## mzrange with more than 2 columns
    checkException(xcms:::profEIC(xraw, mzrange = c(1, 3, 5)))
    checkException(getEIC(xraw, mzrange = c(1, 3, 5)))
    ## rtrange with more than 2 columns
    checkException(xcms:::profEIC(xraw, rtrange = c(2, 4, 6)))
    checkException(getEIC(xraw, rtrange = c(2, 4, 6)))
}

## Comparing the new profEIC function with the result from the
## getEICold and new functions.
notrun_test_profEIC_implementation <- function() {
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    library(xcms)
    library(RUnit)
    step <- 0.1
    xraw <- xcmsRaw(file, profstep=step)
    mzr <- c(200, 201)
    rtr <- c(3000, 3500)
    res_1 <- xcms:::profEIC(xraw, mzrange = mzr, rtrange = rtr)

    res_2 <- getEIC(xraw, rtrange = matrix(rtr, nrow = 1),
                    mzrange = matrix(mzr, nrow = 1))
    checkEquals(res_1, res_2)

    ## Changing settings
    res_1 <- xcms:::profEIC(xraw, mzrange = mzr, rtrange = rtr, step = 0.5)
    res_2 <- getEIC(xraw, rtrange = matrix(rtr, nrow = 1),
                    mzrange = matrix(mzr, nrow = 1), step = 0.5)
    checkEquals(res_1, res_2)
    profMethod(xraw) <- "binlinbase"
    xraw@profparam <- list(baselevel = 666666)
    res_3 <- xcms:::profEIC(xraw, mzrange = mzr, rtrange = rtr, step = 0.5)
    res_4 <- getEIC(xraw, rtrange = matrix(rtr, nrow = 1),
                    mzrange = matrix(mzr, nrow = 1), step = 0.5)
    checkEquals(res_3, res_4)
    ## intlin
    profMethod(xraw) <- "intlin"
    res_3 <- xcms:::profEIC(xraw, mzrange = mzr, rtrange = rtr, step = 0.5)
    res_4 <- getEIC(xraw, rtrange = matrix(rtr, nrow = 1),
                    mzrange = matrix(mzr, nrow = 1), step = 0.5)
    checkEquals(res_3, res_4)

    ## With multiple ranges
    mzr_m <- rbind(c(240, 242), mzr, c(255, 257), c(409, 410))
    rtr_m <- matrix(rep(rtr, each = 4), nrow = 4)
    res_1 <- xcms:::profEIC(xraw, mzrange = mzr_m, rtrange = rtr_m, step = 0.4)
    res_2 <- getEIC(xraw, rtrange = rtr_m, mzrange = mzr_m, step = 0.4)
    checkEquals(res_1, res_2)

    ## Performance tests.
    M <- rbind(mzr_m, mzr_m)
    R <- rbind(rtr_m, rtr_m + 200)
    profMethod(xraw) <- "bin"
    library(microbenchmark)
    checkEquals(getEIC(xraw, mzrange = M, rtrange = R),
                xcms:::profEIC(xraw, mzrange = M, rtrange = R))
    microbenchmark(getEIC(xraw, mzrange = M, rtrange = R),
                   xcms:::profEIC(xraw, mzrange = M, rtrange = R), times = 20)
    ## Speed: profMat considerably faster if the profile matrix is re-used.
    microbenchmark(getEIC(xraw, mzrange = M, rtrange = R, step = 2),
                   xcms:::profEIC(xraw, mzrange = M, rtrange = R, step = 2),
                   times = 20)
    ## Speed: getEIC faster if prof matrix calculated.
    ## Compare speed with existing profile matrix and without.
    xr_1 <- xcmsRaw(file, profstep = 0.1)
    xr_2 <- xcmsRaw(file, profstep = 0)
    res_1 <- xcms:::profEIC(xr_1, mzrange = M, rtrange = R, step = 0.1)
    res_2 <- xcms:::profEIC(xr_2, mzrange = M, rtrange = R, step = 0.1)
    checkEquals(res_1, res_2)
    microbenchmark(xcms:::profEIC(xr_1, mzrange = M, rtrange = R, step = 0.1),
                   xcms:::profEIC(xr_2, mzrange = M, rtrange = R, step = 0.1))
    ## Well, if profile matrix is re-used we're way faster.

    res_3 <- getEIC(xr_2, mzrange = M, rtrange = R, step = 0.1)
    checkEquals(res_3, res_1)
    checkEquals(res_3, res_2)  ## OK, on-the-fly calculation is the same.
    ## manually:
    rtidx <- which(xr_1@scantime >= R[4, 1] & xr_1@scantime <= R[4, 2])
    mass <- seq(floor(xr_1@mzrange[1] / 0.1) * 0.1,
                ceiling(xr_1@mzrange[2] / 0.1) * 0.1, by = 0.1)
    mzIdx <- which(mass >= M[4, 1] & mass <= M[4, 2])
    ecalc <- apply(xr_1@env$profile[mzIdx, rtidx, drop=FALSE], MARGIN=2, max)
}

## This test case checks for the problem reported by Alan Smith in issue #7
## on github.
## Apparently (and actually), the length of the retention time corrected and
## original retention times differ, e.g. when scanrange is used. Using the
## getXcmsRaw method fixes that, since within that method we ensure that all
## vectors, times etc are aligned (i.e. match in length and ordering).
dontrun_test.issue7 <- function(){
    library(xcms)
    library(faahKO)

    cdfpath <- system.file("cdf", package = "faahKO")
    list.files(cdfpath, recursive = TRUE)
    cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)

    xset<-xcmsSet(cdffiles, profmethod = "bin", method = "centWave",
                  ppm = 30, peakwidth = c(5,50), snthresh = 10,
                  prefilter = c(7,1000),
                  integrate = 1, mzdiff = -0.001, fitgauss = FALSE,
                  scanrange = c(0,500), mzCenterFun = "wMean", noise = 300,
                  sleep = 0, verbose.columns = TRUE)
    xset2 <- retcor(xset, method = "obiwarp",distFunc = "cor_opt",
                    profStep = 1, plot = "deviation",gapInit = 0.8,
                    gapExtend = 3.3)
    gxset2 <- group(xset2, method = "density", bw = 8, minfrac = .5,
                    minsamp = 2, mzwid = .01, max = 50, sleep = 0)
    xset3 <- fillPeaks(gxset2)

    peaktab <- data.frame(ID = groupnames(xset3), peakTable(xset3))
    rownames(peaktab) <- peaktab$ID

    pkid <- c("M205T2790","M392T2927")
    a <- getEIC(xset3, groupidx=pkid)
    ## Get the retention times and eic for the first file.
    rts <- a@eic[[1]][[1]]
    ## rt has to be increasing!
    rtidx <- order(rts[, "rt"])
    checkEquals(rtidx, 1:length(rtidx))

    ## Somehow the error below can no longer be reproduced, most likely due to
    ## changes in xcmsRaw.
    ## ## Reproduce what's happening within the getEIC for xcmsSet:
    ## ## That's how the xcmsRaw is loaded.
    ## braw <- xcmsRaw(cdffiles[1], profmethod=profMethod(xset3), profstep=0)
    ## ## Note that the vector lengths are different.
    ## checkTrue(length(braw@scantime) != length(xset3@rt$corrected[[1]]))
    ## ## That's the problem! That's how the corrected scantime is applied.
    ## braw@scantime <- xset3@rt$corrected[[1]]
    ## ## Now the lengths are different and we get the error.
    ## rtrange <- as.matrix(peaktab[pkid, c("rtmin", "rtmax")])
    ## rtrange[,1] <- rtrange[,1] - 200
    ## rtrange[,2] <- rtrange[,2] + 200
    ## beic <- getEIC(braw, mzrange=as.matrix(peaktab[pkid, c("mzmin", "mzmax")]),
    ##                rtrange=rtrange)
    ## plot(beic)
    ## rts <- beic@eic[[1]][[1]]
    ## ## rt has to be increasing!
    ## rtidx <- order(rts[, "rt"])
    ## checkTrue(all(rtidx != 1:length(rtidx)))

    ## Applying also the scanrange (as done in getXcmsRaw) fix it.
    craw <- getXcmsRaw(xset3)
    rtrange <- as.matrix(peaktab[pkid, c("rtmin", "rtmax")])
    rtrange[,1] <- rtrange[,1] - 200
    rtrange[,2] <- rtrange[,2] + 200
    ceic <- getEIC(craw, mzrange=as.matrix(peaktab[pkid, c("mzmin", "mzmax")]),
                   rtrange = rtrange)
    plot(ceic)
    rts <- ceic@eic[[1]][[1]]
    ## rt has to be increasing!
    rtidx <- order(rts[, "rt"])
    checkEquals(rtidx, 1:length(rtidx))

    ## Interestingly, using getEICNew also helps:
    ## BioC <- getOption("BioC")
    ## BioCBefore <- getOption("BioC")
    ## BioC$xcms$getEIC.method <- "getEICNew"
    ## options(BioC=BioC)
    ## zic <- getEIC(xset3, groupidx=pkid[1])
    ## plot(zic)
    ## rts <- zic@eic[[1]][[1]]
    ## ## rt has to be increasing!
    ## rtidx <- order(rts[, "rt"])
    ## checkEquals(rtidx, 1:length(rtidx))

    ## ## Reset options.
    ## options(BioC=BioCBefore)
}

test_getEICxset <- function() {
    ## xset <- fillPeaks(group(faahko))
    xset <- faahko_grouped_filled
    ## xset <- faahko_grouped_filled
    e <- getEIC(xset, sampleidx = c(1,2), groupidx = c(1,2), rtrange=200)
    checkEquals(sampnames(e), c("ko15", "ko16"))
    ## plot(e)
    ## Reproduce issue #92
    e <- getEIC(xset, sampleidx = c(5, 9), groupidx = c(1, 2), rtrange = 200)
    checkEquals(sampnames(e), sampnames(xset)[c(5, 9)])
    ## Compare with raw data.
    rtr <- matrix(c(2876, 2932), nrow = 1)
    mzr <- matrix(c(200.1, 200.1), nrow = 1)
    e <- getEIC(xset, sampleidx = c(1), mzrange = mzr,
                rtrange = rtr, rt = "raw")
    ## Read the raw data of file 1:
    xr <- xcmsRaw(filepaths(xset)[1], profstep = profStep(xset))
    e_2 <- getEIC(xr, mzrange = mzr, rtrange = rtr, step = 0.1)
    checkEquals(e_2@eic[[1]][[1]], e@eic[[1]][[1]])
    ## Check what happens if we select another -> issue #92
    e <- getEIC(xset, sampleidx = c(5, 9), mzrange = mzr,
                rtrange = rtr, rt = "raw")
    ## sample 5
    xr <- xcmsRaw(filepaths(xset)[5], profstep = profStep(xset))
    e_2 <- getEIC(xr, mzrange = mzr, rtrange = rtr, step = 0.1)
    checkEquals(e_2@eic[[1]][[1]], e@eic[[1]][[1]])
    ## sample 9
    xr <- xcmsRaw(filepaths(xset)[9], profstep = profStep(xset))
    e_2 <- getEIC(xr, mzrange = mzr, rtrange = rtr, step = 0.1)
    checkEquals(e_2@eic[[1]][[1]], e@eic[[2]][[1]])
}

test.getEICretcor <- function() {
    ## xset <- fillPeaks(group(retcor(group(faahko))))
    xset <- faahko_grouped_retcor_filled
    ## xset <- faahko_processed
    opt.warn <- options("warn")$warn
    options("warn" = 2) ## turns warning into errors
    e <- getEIC(xset, sampleidx=c(1,2), groupidx=c(1,2),
                rt="corrected", rtrange=200)
    options("warn" = opt.warn)
    plot(e)
}

test.plotEIC <- function() {
    ## file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    ## xraw <- xcmsRaw(file)
    xraw <- deepCopy(faahko_xr_1)
    e <- getEIC(xraw, rtrange=cbind(3000,3500), mzrange=cbind(200,201))
    plot(e)
}
