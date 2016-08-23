## Test detectFeatures matchedFilter

library(xcms)
library(RUnit)

library(faahKO)
fs <- c(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
        system.file('cdf/KO/ko16.CDF', package = "faahKO"),
        system.file('cdf/KO/ko18.CDF', package = "faahKO"),
        system.file('cdf/KO/ko19.CDF', package = "faahKO"))

test_do_detectFeatures_matchedFilter <- function() {
    xr <- xcmsRaw(fs[1])
    ## We expect that changing a parameter has an influence on the result.
}


fails_do_matchedFilter_tackle_differences <- function() {
    ## Try to tackle what causes the problem: the iterative buffering or the
    ## binning function. The puzzling thing is that there is a problem for
    ## step = 0.2, 0.4 and similar, but not for other step sizes!
    ##

    library(xcms)
    library(RUnit)

    library(faahKO)
    fs <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xr <- xcmsRaw(fs)
    mz <- xr@env$mz
    int <- xr@env$intensity
    scantime <- xr@scantime
    scanindex <- xr@scanindex

    ## do and do* new HAVE to yield same results since for both the profFun is
    ## used.
    step <- 0.1
    resDo <- xcms:::do_detectFeatures_matchedFilter(mz, int, scantime,
                                                    diff(c(scanindex, length(mz))),
                                                    step = step)
    resDo_new <- xcms:::do_detectFeatures_matchedFilter_no_iter(mz, int, scantime,
                                                                diff(c(scanindex, length(mz))),
                                                                step = step)
    checkEquals(resDo, resDo_new)
    ## newbin: new binning function but also iterative buffer filling
    resDo_newbin <- xcms:::do_detectFeatures_matchedFilter_binYonX_iter(mz, int,
                                                                        scantime,
                                                                        diff(c(scanindex,
                                                                               length(mz))),
                                                                        step = step)
    checkEquals(resDo, resDo_newbin)


    ## THAT'S PUZZLING!!!
    step <- 0.2
    ## Original code.
    resDo <- xcms:::do_detectFeatures_matchedFilter(mz, int, scantime,
                                                    diff(c(scanindex, length(mz))),
                                                    step = step)
    resDo2 <- xcms:::do_detectFeatures_matchedFilter(mz, int, scantime,
                                                    diff(c(scanindex, length(mz))),
                                                    step = step)
    ## Original code, but no iterative buffer binning.
    resDo_new <- xcms:::do_detectFeatures_matchedFilter_no_iter(mz, int, scantime,
                                                                diff(c(scanindex, length(mz))),
                                                                step = step)
    checkEquals(resDo, resDo_new)  ## FAILS!
    ## Where is the difference?
    for (i in 1:ncol(resDo)) {
        cat("* ", colnames(resDo)[i], ": ")
        cat(all.equal(sort(resDo[, i]), sort(resDo_new[, i])), "\n")
        ##cat(paste(which(resDo[, i] != resDo_new[, i]), collapse = ", "), "\n")
    }

    resDo_newbin <- xcms:::do_detectFeatures_matchedFilter_binYonX_iter(mz, int,
                                                                        scantime,
                                                                        diff(c(scanindex,
                                                                               length(mz))),
                                                                        step = step)
    checkEquals(resDo, resDo_newbin)  ## FAILS!
    checkEquals(resDo_new, resDo_newbin)  ## FAILS!!

    ## See dontrun_tackle_binning_differences function in
    ## test_binning.R for reasons of the difference!
}

############################################################
## This is only relevant during development of the do_ function
## to evaluate that results are identical.
dontrun_test_do_detectFeatures_matchedFilter_impl <- function() {

    library(xcms)
    library(RUnit)

    i <- 1

    ## Run this an all files...
    for (i in 1:length(fs)) {

        xr <- xcmsRaw(fs[i])
        profFun <- "bin"
        profparam <- list()
        fwhm <- 30
        sigma <- fwhm / 2.3548
        max <- 5
        snthresh <- 10
        step <- 0.1
        steps <- 2
        mzdiff <- 0.8 - step * steps
        index <- FALSE
        verbose.columns <- FALSE

        #######
        ## bin
        .runAndCompare_matchedFilter(xr, profFun, profparam, fwhm, sigma, max,
                                     snthresh, step, steps, mzdiff, index,
                                     verbose.columns)
        steps <- 4
        mzdiff <- 0.8 - step * steps
        .runAndCompare_matchedFilter(xr, profFun, profparam, fwhm, sigma, max,
                                     snthresh, step, steps, mzdiff, index,
                                     verbose.columns)
        steps <- 5
        mzdiff <- 0.8 - step * steps
        .runAndCompare_matchedFilter(xr, profFun, profparam, fwhm, sigma, max,
                                     snthresh, step, steps, mzdiff, index,
                                     verbose.columns)
        step <- 0.2  ## LLLL FIX THIS! difference between xrDo_new and xrPeaks
        ## This is a puzzling one!, it's the same for 0.4, 0.6, 0.8, 1.2
        ## Not for any other number up to 1.4 (including 1.0, 1.4).
        steps <- 2
        mzdiff <- 0.8 - step * steps
        .runAndCompare_matchedFilter(xr, profFun, profparam, fwhm, sigma, max,
                                     snthresh, step, steps, mzdiff, index,
                                     verbose.columns)

        step <- 0.13
        steps <- 2
        mzdiff <- 0.8 - step * steps
        .runAndCompare_matchedFilter(xr, profFun, profparam, fwhm, sigma, max,
                                     snthresh, step, steps, mzdiff, index,
                                     verbose.columns)

        snthresh <- 4
        .runAndCompare_matchedFilter(xr, profFun, profparam, fwhm, sigma, max,
                                     snthresh, step, steps, mzdiff, index,
                                     verbose.columns)
        fwhm <- 15
        sigma <- fwhm / 2.3548
        .runAndCompare_matchedFilter(xr, profFun, profparam, fwhm, sigma, max,
                                     snthresh, step, steps, mzdiff, index,
                                     verbose.columns)

        ## #####
        ## using binlinbase.
        profFun <- "binlinbase"
        profparam <- list()
        steps <- 2
        step <- 0.1
        mzdiff <- 0.8 - step * steps
        snthresh <- 10
        fwhm <- 30
        sigma <- fwhm / 2.3548
        .runAndCompare_matchedFilter(xr, profFun, profparam, fwhm, sigma, max,
                                     snthresh, step, steps, mzdiff, index,
                                     verbose.columns)
        ## Increase the basespace to 0.1 so that we're actually doing some
        ## interpolation (the value is actually slightly larger to ensure that
        ## really interpolation is taking place).
        profparam$basespace <- 0.101
        .runAndCompare_matchedFilter(xr, profFun, profparam, fwhm, sigma, max,
                                     snthresh, step, steps, mzdiff, index,
                                     verbose.columns)
        steps <- 4
        mzdiff <- 0.8 - step * steps
        .runAndCompare_matchedFilter(xr, profFun, profparam, fwhm, sigma, max,
                                     snthresh, step, steps, mzdiff, index,
                                     verbose.columns)
        step <- 0.4
        mzdiff <- 0.8 - step * steps
        .runAndCompare_matchedFilter(xr, profFun, profparam, fwhm, sigma, max,
                                     snthresh, step, steps, mzdiff, index,
                                     verbose.columns)
        snthresh <- 4
        .runAndCompare_matchedFilter(xr, profFun, profparam, fwhm, sigma, max,
                                     snthresh, step, steps, mzdiff, index,
                                     verbose.columns)
        fwhm <- 15
        sigma <- fwhm / 2.3548
        .runAndCompare_matchedFilter(xr, profFun, profparam, fwhm, sigma, max,
                                     snthresh, step, steps, mzdiff, index,
                                     verbose.columns)


        ## #####
        ## using binlin
        ##
        ## We're only checking if the iterative approach matches the non-iterative.
        ## We're not evaluating the binYonX and imputeLinInterpol!
        profFun <- "binlin"
        profparam <- list()
        steps <- 2
        step <- 0.1
        snthresh <- 10
        fwhm <- 30
        sigma <- fwhm / 2.3548
        .runAndCompare_matchedFilter(xr, profFun, profparam, fwhm, sigma, max,
                                     snthresh, step, steps, mzdiff, index,
                                     verbose.columns, skipNew = TRUE)
        steps <- 4
        mzdiff <- 0.8 - step * steps
        .runAndCompare_matchedFilter(xr, profFun, profparam, fwhm, sigma, max,
                                     snthresh, step, steps, mzdiff, index,
                                     verbose.columns, skipNew = TRUE)
        step <- 0.4
        mzdiff <- 0.8 - step * steps
        .runAndCompare_matchedFilter(xr, profFun, profparam, fwhm, sigma, max,
                                     snthresh, step, steps, mzdiff, index,
                                     verbose.columns, skipNew = TRUE)
        snthresh <- 4
        .runAndCompare_matchedFilter(xr, profFun, profparam, fwhm, sigma, max,
                                     snthresh, step, steps, mzdiff, index,
                                     verbose.columns, skipNew = TRUE)
        fwhm <- 15
        sigma <- fwhm / 2.3548
        .runAndCompare_matchedFilter(xr, profFun, profparam, fwhm, sigma, max,
                                     snthresh, step, steps, mzdiff, index,
                                     verbose.columns, skipNew = TRUE)

    }
}


.runAndCompare_matchedFilter <- function(xr, profFun, profparam, fwhm, sigma,
                                         max, snthresh, step, steps, mzdiff,
                                         index, verbose.columns,
                                         skipNew = FALSE) {
    require(RUnit)
    mz <- xr@env$mz
    int <- xr@env$intensity
    scantime <- xr@scantime
    scanindex <- xr@scanindex
    profMethod(xr) <- profFun
    xr@profparam <- profparam
    if (any(names(profparam) == "basespace")) {
        baseSpace <- profparam$basespace
    } else {
        baseSpace <- 0.075
    }
    ## Run the findPeaks.matchedFilter on the xcmsRaw.
    a <- system.time(
        xrPeaks <- findPeaks.matchedFilter(xr,
                                           fwhm = fwhm,
                                           sigma = sigma,
                                           max = max,
                                           snthresh = snthresh,
                                           step = step,
                                           steps = steps,
                                           mzdiff = mzdiff,
                                           index = index,
                                           verbose.columns = verboseColumns)
    )  ##
    cat("XCMS matchedFilter: ", a, "\n")
    b <- system.time(
        xrDo <- xcms:::do_detectFeatures_matchedFilter(mz, int, scantime,
                                                       diff(c(scanindex, length(mz))),
                                                       profFun = profFun,
                                                       profparam = profparam,
                                                       fwhm = fwhm,
                                                       sigma = sigma,
                                                       max = max,
                                                       snthresh = snthresh,
                                                       step = step,
                                                       steps = steps,
                                                       mzdiff = mzdiff,
                                                       index = index,
                                                       verboseColumns = verboseColumns)
    ) ##
    cat("do_ original code: ", b, "\n")
    if (!checkEquals(new("xcmsPeaks", xrDo), xrPeaks))
        stop("do_ and xcms yield different results!")
    c <- system.time(
        xrDo_new <- xcms:::do_detectFeatures_matchedFilter_new(mz, int, scantime,
                                                               diff(c(scanindex, length(mz))),
                                                               profFun = profFun,
                                                               profparam = profparam,
                                                               fwhm = fwhm,
                                                               sigma = sigma,
                                                               max = max,
                                                               snthresh = snthresh,
                                                               step = step,
                                                               steps = steps,
                                                               mzdiff = mzdiff,
                                                               index = index,
                                                               verboseColumns = verboseColumns)
    ) ##
    if (!checkEquals(new("xcmsPeaks", xrDo_new), xrPeaks))
        stop("do_*_new and xcms yield different results!")
    cat("do_*_new full matrix binning: ", c, "\n")
    if (!skipNew) {
        d <- system.time(
            xrDo_newer <- xcms:::do_detectFeatures_matchedFilter_newer(mz, int, scantime,
                                                                       diff(c(scanindex, length(mz))),
                                                                       profFun = profFun,
                                                                       baseSpace = baseSpace,
                                                                       fwhm = fwhm,
                                                                       sigma = sigma,
                                                                       max = max,
                                                                       snthresh = snthresh,
                                                                       step = step,
                                                                       steps = steps,
                                                                       mzdiff = mzdiff,
                                                                       index = index,
                                                                       verboseColumns = verboseColumns)
        ) ##
        cat("do_*_newer binYonX: ", d, "\n")
        if (!checkEquals(new("xcmsPeaks", xrDo_newer), xrPeaks))
            stop("do_*_newer and xcms yield different results!")
    }
}


## Some speed tests.
.otherTest <- function() {
    Testv <- c(2, 4.2, 34.1, 34.5, 6.4, 6.3, 1.2)
    RforM <- matrix(nrow = 0, ncol = length(Testv))
    system.time(
        for(i in 1:5000){
            RforM <- rbind(RforM, Testv)
        }
    ) ## 1.27
    ## with append to list.
    RforL <- vector("list", 0)
    system.time(
        for(i in 1:5000){
            RforL <- c(RforL, Testv)
        }
    ) ## 1.12
    system.time(
        RapplyL <- lapply(1:5000, function(z) {return(Testv)})
    ) ## 0.003
    RM <- matrix(nrow=5000, ncol = length(Testv))
    system.time(
        for (i in 1:5000) {
            RM[i, ] <- Testv
        }
    ) ## 0.006

    ## Compare adding to list instead of adding to existing. [[]]
    RexL <- vector("list", 5000)
    system.time(
        for (i in 1:5000){
            RexL[[i]] <- Testv
        }
    ) ## 0.005
    ## Dynamically...
    RexL <- list()
    system.time(
        for (i in 1:5000){
            RexL[[i]] <- Testv
        }
    ) ## 0.005

}

############################################################
## Lengthy explanation:
## Point is that the division taking place inside the ProfBinLinBase C-function
## to determine how many neighboring bins should be included in the linear
## interpolation can become unstable if step and basespace are set to the
## same value; this does seem to take place specifically if iterative
## buffer creation is used.
notrun_odd_matchedFilter_behaviour <- function() {

    library(xcms)
    library(RUnit)
    library(faahKO)
    fs <- c(system.file('cdf/KO/ko15.CDF', package = "faahKO"))
    xr <- xcmsRaw(fs)
    profparam <- list()
    ## That's stable
    profparam$basespace <- 0.10000001
    step <- 0.1
    xr@profparam <- profparam
    profMethod(xr) <- "binlinbase"
    res_stable <- findPeaks.matchedFilter(xr, step = step)

    ## Do the same with 0.1
    profparam$basespace <- 0.1
    xr@profparam <- profparam
    res_unstable <- findPeaks.matchedFilter(xr, step = step)

    checkTrue(nrow(res_stable) != nrow(res_unstable))

    ## Checking the number of peaks we would get if we're switching on or off
    ## the interpolation.
    profparam$basespace <- 0.05
    xr@profparam <- profparam
    res_no_inter <- findPeaks.matchedFilter(xr, step = step)
    profparam$basespace <- 0.13
    xr@profparam <- profparam
    res_with_inter <- findPeaks.matchedFilter(xr, step = step)

    ## Check what we've got:
    checkEquals(res_with_inter, res_stable)
    checkTrue(nrow(res_no_inter) > nrow(res_with_inter))
}
