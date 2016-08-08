## Test detectFeatures matchedFilter

## library(xcms)
## library(RUnit)

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
## This is only relevant during development of the do_ function
## to evaluate that results are identical.
dontrun_test_do_detectFeatures_matchedFilter_impl <- function() {

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
    i <- 1

    ## Run this an all files...
    for (i in 1:length(fs)) {
        xr <- xcmsRaw(fs[i])

        ## Default settings
        .runAndCompare_matchedFilter(xr, profFun, profparam, fwhm, sigma, max,
                                     snthresh, step, steps, mzdiff, index,
                                     verbose.columns)
        profFun <- "binlin"
        .runAndCompare_matchedFilter(xr, profFun, profparam, fwhm, sigma, max,
                                     snthresh, step, steps, mzdiff, index,
                                     verbose.columns)
        step <- 0.2
        .runAndCompare_matchedFilter(xr, profFun, profparam, fwhm, sigma, max,
                                     snthresh, step, steps, mzdiff, index,
                                     verbose.columns)
        snthresh <- 4
        .runAndCompare_matchedFilter(xr, profFun, profparam, fwhm, sigma, max,
                                     snthresh, step, steps, mzdiff, index,
                                     verbose.columns)
        steps <- 5
        .runAndCompare_matchedFilter(xr, profFun, profparam, fwhm, sigma, max,
                                     snthresh, step, steps, mzdiff, index,
                                     verbose.columns)
        profFun <- "binlinbase"
        .runAndCompare_matchedFilter(xr, profFun, profparam, fwhm, sigma, max,
                                     snthresh, step, steps, mzdiff, index,
                                     verbose.columns)
        profFun <- "intlin"
        .runAndCompare_matchedFilter(xr, profFun, profparam, fwhm, sigma, max,
                                     snthresh, step, steps, mzdiff, index,
                                     verbose.columns)
    }
}


.runAndCompare_matchedFilter <- function(xr, profFun, profparam, fwhm, sigma, max, snthresh,
                                         step, steps, mzdiff, index, verbose.columns) {
    require(RUnit)
    mz <- xr@env$mz
    int <- xr@env$intensity
    scantime <- xr@scantime
    scanindex <- xr@scanindex
    profMethod(xr) <- profFun
    xr@profparam <- profparam
    a <- system.time(
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
    ## Run the findPeaks.centWave on the xcmsRaw.
    b <- system.time(
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
    ## Compare.
    ## For now the problem seems to be that arguments are not passed correctly!
    cat("DO: ", a, "\n")
    cat("XCMS: ", b, "\n")
    if (!checkEquals(new("xcmsPeaks", xrDo), xrPeaks))
        stop("do_ and xcms yield different results!")
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
