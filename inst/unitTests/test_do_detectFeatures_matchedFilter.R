## Test detectFeatures matchedFilter
## LLLL clean up.
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

############################################################
## Compare each individual function to the original one changing
## settings.
## Comparing each of the functions to the original one:
## A: do_detectFeatures_matchedFilter
## B: do_detectFeatures_matchedFilter_binYonX_iter
## C: do_detectFeatures_matchedFilter_no_iter
## D: do_detectFeatures_matchedFilter_binYonX_no_iter
dontrun_test_do_detectFeatures_matchedFilter_impl <- function {

    library(xcms)
    library(RUnit)

    library(faahKO)
    fs <- c(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
            system.file('cdf/KO/ko16.CDF', package = "faahKO"),
            system.file('cdf/KO/ko18.CDF', package = "faahKO"),
            system.file('cdf/KO/ko19.CDF', package = "faahKO"))
    i <- 1

    for (i in 1:length(fs)) {

        xr <- xcmsRaw(fs[i])
        mz <- xr@env$mz
        int <- xr@env$intensity
        scantime <- xr@scantime
        scanindex <- xr@scanindex
        valsPerSpect <- diff(c(scanindex, length(mz)))

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

        .compare_em(xr, mz, int, scantime, valsPerSpect, profFun, profparam,
                    fwhm, sigma, max, snthresh, step, steps, mzdiff, index,
                    verbose.columns)    # OK
        steps <- 4
        mzdiff <- 0.8 - step * steps
        .compare_em(xr, mz, int, scantime, valsPerSpect, profFun, profparam,
                    fwhm, sigma, max, snthresh, step, steps, mzdiff, index,
                    verbose.columns)    # OK
        steps <- 5
        mzdiff <- 0.8 - step * steps
        .compare_em(xr, mz, int, scantime, valsPerSpect, profFun, profparam,
                    fwhm, sigma, max, snthresh, step, steps, mzdiff, index,
                    verbose.columns)    # OK
        ## THIS WILL FAIL!!! REASONE: rounding error in bin definition.
        step <- 0.2
        steps <- 2
        mzdiff <- 0.8 - step * steps
        .compare_em(xr, mz, int, scantime, valsPerSpect, profFun, profparam,
                    fwhm, sigma, max, snthresh, step, steps, mzdiff, index,
                    verbose.columns, stopOnError = FALSE)    #
        step <- 0.1999
        steps <- 2
        mzdiff <- 0.8 - step * steps
        .compare_em(xr, mz, int, scantime, valsPerSpect, profFun, profparam,
                    fwhm, sigma, max, snthresh, step, steps, mzdiff, index,
                    verbose.columns, stopOnError = TRUE)    # OK
        snthresh <- 4
        .compare_em(xr, mz, int, scantime, valsPerSpect, profFun, profparam,
                    fwhm, sigma, max, snthresh, step, steps, mzdiff, index,
                    verbose.columns, stopOnError = TRUE)    # OK
        fwhm <- 15
        sigma <- fwhm / 2.3548
        .compare_em(xr, mz, int, scantime, valsPerSpect, profFun, profparam,
                    fwhm, sigma, max, snthresh, step, steps, mzdiff, index,
                    verbose.columns, stopOnError = TRUE)    # OK

        ## ######
        ## binLin

        ## ######
        ## binLinBase
    }

}


.compare_em <- function(xr, mz, int, scantime, valsPerSpect, profFun, profparam,
                        fwhm, sigma, max, snthresh, step, steps, mzdiff, index,
                        verbose.columns, stopOnError = TRUE) {
    profMethod(xr) <- profFun
    xr@profparam <- profparam
    if (any(names(profparam) == "basespace")) {
        baseSpace <- profparam$basespace
    } else {
        baseSpace <- 0.075
    }
    ## The reference is the old code.
    orig <- findPeaks.matchedFilter(xr,
                                    fwhm = fwhm,
                                    sigma = sigma,
                                    max = max,
                                    snthresh = snthresh,
                                    step = step,
                                    steps = steps,
                                    mzdiff = mzdiff,
                                    index = index,
                                    verbose.columns)@.Data

    ## A
    A <- xcms:::do_detectFeatures_matchedFilter(mz, int, scantime,
                                                valsPerSpect,
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
                                                verbose.columns)
    res <- all.equal(orig, A)
    if (is.character(res)) {
        msg <- paste0("A vs original FAILED! ", res, "\n")
        cat(msg)
        if (stopOnError)
            stop(msg)
    } else {
        cat("A vs original OK\n")
    }
    ## B
    B <- xcms:::do_detectFeatures_matchedFilter_binYonX_iter(mz, int, scantime,
                                                             valsPerSpect,
                                                             profFun = profFun,
                                                             fwhm = fwhm,
                                                             sigma = sigma,
                                                             max = max,
                                                             snthresh = snthresh,
                                                             step = step,
                                                             steps = steps,
                                                             mzdiff = mzdiff,
                                                             index = index,
                                                             verbose.columns)
    res <- all.equal(orig, B)
    if (is.character(res)) {
        msg <- paste0("B vs original FAILED! ", res, "\n")
        cat(msg)
        if (stopOnError)
            stop(msg)
    } else {
        cat("B vs original OK\n")
    }
    ## C
    C <- xcms:::do_detectFeatures_matchedFilter_no_iter(mz, int, scantime,
                                                     valsPerSpect,
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
                                                     verbose.columns)
    res <- all.equal(orig, C)
    if (is.character(res)) {
        msg <- paste0("C vs original FAILED! ", res, "\n")
        cat(msg)
        if (stopOnError)
            stop(msg)
    } else {
        cat("C vs original OK\n")
    }
    ## D
    D <- xcms:::do_detectFeatures_matchedFilter_binYonX_no_iter(mz, int, scantime,
                                                                valsPerSpect,
                                                                profFun = profFun,
                                                                fwhm = fwhm,
                                                                sigma = sigma,
                                                                max = max,
                                                                snthresh = snthresh,
                                                                step = step,
                                                                steps = steps,
                                                                mzdiff = mzdiff,
                                                                index = index,
                                                                verbose.columns)
    res <- all.equal(orig, D)
    if (is.character(res)) {
        msg <- paste0("D vs original FAILED! ", res, "\n")
        cat(msg)
        if (stopOnError)
            stop(msg)
    } else {
        cat("D vs original OK\n")
    }
    ## Now compare between:
    res <- all.equal(B, C)
    if (is.character(res)) {
        msg <- paste0("B vs C FAILED! ", res, "\n")
        cat(msg)
    } else {
        cat("B vs C OK\n")
    }
    res <- all.equal(B, D)
    if (is.character(res)) {
        msg <- paste0("B vs D FAILED! ", res, "\n")
        cat(msg)
    } else {
        cat("B vs D OK\n")
    }
    res <- all.equal(C, D)
    if (is.character(res)) {
        msg <- paste0("C vs D FAILED! ", res, "\n")
        cat(msg)
    } else {
        cat("C vs D OK\n")
    }
}


############################################################
## This is only relevant during development of the do_ function
## to evaluate that results are identical.
## Actually, we don't expect all results to be identical:
## o profBinLin is buggy.
## o A step = 0.2 leads to slightly different bin definitions (C rounding)
##   and thus to different results.
dontrun_test_do_detectFeatures_matchedFilter_impl_old <- function() {

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
        .matchedFilter_iter(xr, profFun, profparam, fwhm, sigma, max,
                            snthresh, step, steps, mzdiff, index,
                            verbose.columns)
        .matchedFilter_no_iter(xr, profFun, profparam, fwhm, sigma, max,
                               snthresh, step, steps, mzdiff, index,
                               verbose.columns)
        ## OK
        steps <- 4
        mzdiff <- 0.8 - step * steps
        .matchedFilter_iter(xr, profFun, profparam, fwhm, sigma, max,
                            snthresh, step, steps, mzdiff, index,
                            verbose.columns)
        .matchedFilter_no_iter(xr, profFun, profparam, fwhm, sigma, max,
                               snthresh, step, steps, mzdiff, index,
                               verbose.columns)
        ## OK
        steps <- 5
        mzdiff <- 0.8 - step * steps
        .matchedFilter_iter(xr, profFun, profparam, fwhm, sigma, max,
                            snthresh, step, steps, mzdiff, index,
                            verbose.columns)
        .matchedFilter_no_iter(xr, profFun, profparam, fwhm, sigma, max,
                               snthresh, step, steps, mzdiff, index,
                               verbose.columns)
        ## OK
        step <- 0.2  ## LLLL FIX THIS! difference between xrDo_new and xrPeaks
        ## This is a puzzling one!, it's the same for 0.4, 0.6, 0.8, 1.2
        ## Not for any other number up to 1.4 (including 1.0, 1.4).
        steps <- 2
        mzdiff <- 0.8 - step * steps
        .matchedFilter_iter(xr, profFun, profparam, fwhm, sigma, max,
                            snthresh, step, steps, mzdiff, index,
                            verbose.columns)
        ## FAIL: These differences are due to small binning differences.

        step <- 0.13
        steps <- 2
        mzdiff <- 0.8 - step * steps
        .matchedFilter_iter(xr, profFun, profparam, fwhm, sigma, max,
                            snthresh, step, steps, mzdiff, index,
                            verbose.columns)
        ## OK
        snthresh <- 4
        .matchedFilter_iter(xr, profFun, profparam, fwhm, sigma, max,
                            snthresh, step, steps, mzdiff, index,
                            verbose.columns)
        ## OK
        fwhm <- 15
        sigma <- fwhm / 2.3548
        .matchedFilter_iter(xr, profFun, profparam, fwhm, sigma, max,
                            snthresh, step, steps, mzdiff, index,
                            verbose.columns)
        ## OK

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


############################################################
## This compares all iterative approaches:
## o Original code: findPeaks.matchedFilter
## o Original code within a "do_" function: do_detectFeatures_matchedFilter
## o New implementation using binYonX: do_detectFeatures_matchedFilter_binYonX_iter
.matchedFilter_iter <- function(xr, profFun, profparam, fwhm,
                                sigma, max, snthresh, step, steps,
                                mzdiff, index, verbose.columns,
                                stopOnError = TRUE) {
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
    ## The reference is the old code.
    ref <- findPeaks.matchedFilter(xr,
                                   fwhm = fwhm,
                                   sigma = sigma,
                                   max = max,
                                   snthresh = snthresh,
                                   step = step,
                                   steps = steps,
                                   mzdiff = mzdiff,
                                   index = index,
                                   verbose.columns = verboseColumns)
    ## Old code as-is.
    old_iter <- xcms:::do_detectFeatures_matchedFilter(mz, int, scantime,
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
    res <- all.equal(ref@.Data, old_iter)
    ## if (!checkEquals(new("xcmsPeaks", old_iter), ref)) {
    if (is.character(res)) {
        msg <- paste0("do_*_matchedFilter and findPeaks.matchedFilter yield",
                      " different results: ", res, "\n")
        if (stopOnError) {
            stop(msg)
        } else {
            cat(msg)
            warning(msg)
        }
    } else {
        cat("Same results for do_*_matchedFilter and findPeaks.matchedFilter.\n")
    }
    ## new code but with iteration.
    new_iter <- xcms:::do_detectFeatures_matchedFilter_binYonX_iter(mz, int, scantime,
                                                                    diff(c(scanindex,
                                                                           length(mz))),

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
    res <- all.equal(old_iter, new_iter)
    ## if (!checkEquals(old_iter, new_iter)){
    if (is.character(res)) {
        msg <- paste0("do_*_matchedFilter and do_*_matchedFilter_binYonX_iter",
                      " yield different results: ", res, "\n")
        if (stopOnError) {
            stop(msg)
        } else {
            cat(msg)
            warning(msg)
        }
    } else {
        cat("Same results for do_*_matchedFilter and",
            " do_*_matchedFilter_binYonX_iter.\n")
    }
}


############################################################
## This compares all iterative agains non-iterative approaches:
## o Original code: findPeaks.matchedFilter
## o Single matrix generation using original code:
##   do_detectFeatures_matchedFilter_no_iter
## o New implementation using binYonX:
##   do_detectFeatures_matchedFilter_binYonX_no_iter
.matchedFilter_no_iter <- function(xr, profFun, profparam, fwhm,
                                   sigma, max, snthresh, step, steps,
                                   mzdiff, index, verbose.columns,
                                   stopOnError = TRUE) {
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
    ## The reference is the old code.
    ref <- findPeaks.matchedFilter(xr,
                                   fwhm = fwhm,
                                   sigma = sigma,
                                   max = max,
                                   snthresh = snthresh,
                                   step = step,
                                   steps = steps,
                                   mzdiff = mzdiff,
                                   index = index,
                                   verbose.columns = verboseColumns)
    ## Old code as-is.
    old <- xcms:::do_detectFeatures_matchedFilter_no_iter(mz, int, scantime,
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
    ## if (!checkEquals(new("xcmsPeaks", old), ref)) {
    res <- all.equal(ref@.Data, old)
    if (is.character(res)) {
        msg <- paste0("do_*_matchedFilter_no_iter and findPeaks.matchedFilter yield",
                      " different results: ", res, "\n")
        if (stopOnError) {
            stop(msg)
        } else {
            cat(msg)
            warning(msg)
        }
    } else {
        cat("Same results for do_*_matchedFilter_no_iter",
            " and findPeaks.matchedFilter.\n")
    }
    ## new code but with iteration.
    new <- xcms:::do_detectFeatures_matchedFilter_binYonX_no_iter(mz, int, scantime,
                                                                  diff(c(scanindex,
                                                                         length(mz))),

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
    ## if (!checkEquals(old, new)){
    res <- all.equal(old, new)
    if (is.character(res)) {
        msg <- paste0("do_*_matchedFilter_no_iter and do_*_matchedFilter_binYonX_no_iter",
                      " yield different results: ", res, "\n")
        if (stopOnError) {
            stop(msg)
        } else {
            cat(msg)
            warning(msg)
        }
    } else {
        cat("Same results for do_*_matchedFilter_no_iter and",
            " do_*_matchedFilter_binYonX_no_iter.\n")
    }
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
