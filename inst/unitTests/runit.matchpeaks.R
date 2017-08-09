## Check functions in matchpeaks.R
dontrun_matchpeaks <- function() {
    library(xcms)
    library(RUnit)
    library(msdata)
    mzdatapath <- system.file("fticr", package = "msdata")
    mzdatafiles <- list.files(mzdatapath, recursive = TRUE, full.names = TRUE)

    xs4 <- xcmsSet(
        method = "MSW",
        files = mzdatafiles[1],
        scales = c(1,4, 9),
        nearbyPeak = T,
        verbose.columns = FALSE,
        winSize.noise = 500,
        SNR.method = "data.mean",
        snthr = 10)

    ## Define the calibrants.
    masslist <- xs4@peaks[c(1, 4, 7), "mz"]
    xs4@peaks[,"mz"] <- xs4@peaks[,"mz"] +
        0.00001*runif(1,0,0.4)*xs4@peaks[,"mz"] + 0.0001
    
    xs4c <- calibrate(xs4,
                      calibrants=masslist,
                      method="edgeshift",
                      mzabs=0.0001,
                      mzppm=5,
                      neighbours=3,
                      plotres=FALSE
                      )
    ## Do the steps separately
    pkl <- peaks(xs4)
    mpks <- xcms:::matchpeaks(pkl, masslist)
    estm <- xcms:::estimate(mpks)
}

test_matchpeaks <- function() {
    library(xcms)
    faahko_file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    
    faahko_xs <- xcmsSet(faahko_file, profparam = list(step = 0),
                         method = "centWave", noise = 10000, snthresh = 40)
    pks <- peaks(faahko_xs)

    calibs <- pks[c(3, 5, 7, 13, 17, 29), "mz"]

    res <- xcms:::matchpeaks(pks, calibs)
    res_2 <- xcms:::matchpeaks(pks[order(pks[, "mz"]), ], calibs)
}

dontrun_implementation_matchpeaks <- function() {
    ## Check that xcms:::matchpeaks and xcms:::.matchpeaks2 return the same
    ## results.
    set.seed(123)
    pks <- cbind(mz = sort(abs(rnorm(300, mean = 200, sd = 3))),
                 into = abs(rnorm(300, mean = 2000, sd = 700)))
    masses_idx <- sort(sample(1:nrow(pks), size = 50))
    masses <- pks[masses_idx, "mz"]
    res_old <- xcms:::matchpeaks(pks, masses)
    res_new <- xcms:::.matchpeaks2(pks, masses)
    checkEquals(res_old, res_new)

    res_old <- xcms:::matchpeaks(pks, masses, mzabs = 0, mzppm = 0)
    checkEquals(res_old[, "pos"], masses_idx)
    res_new <- xcms:::.matchpeaks2(pks, masses, mzabs = 0, mzppm = 0)
    checkEquals(res_old, res_new)

    res_old <- xcms:::matchpeaks(pks, masses, mzppm = 0, mzabs = 0.01)
    res_new <- xcms:::.matchpeaks2(pks, masses, mzppm = 0, mzabs = 0.01)
    checkEquals(res_old, res_new)
    
    
    ## Real data... peaks have to be mz sorted! Report that as issue?
    pks <- chromPeaks(faahko_xod)
    pks_2 <- pks[pks[, "sample"] == 2, ]
    masses_idx <- c(4, 13, 32, 33, 37, 41, 45, 53, 58, 67, 74, 88, 90)

    masses <- pks_2[masses_idx, "mz"]
    masses_order <- order(masses)
    res_old <- xcms:::matchpeaks(pks_2, masses)
    res_new <- xcms:::.matchpeaks2(pks_2, masses)

    ## Why the heck does matchpeaks not work???
    
    masses <- masses + 3 * masses / 1e6

    library(microbenchmark)
    microbenchmark(xcms:::matchpeaks(pks, masses),
                   xcms:::.matchpeaks2(pks, masses))
    ## 10x faster.
}
