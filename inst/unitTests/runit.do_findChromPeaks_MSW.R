############################################################
## do_findPeaks_MSW tests

## xraw <- deepCopy(microtofq_xr)

test_do_findPeaks_MSW <- function() {
    first_file <- filterFile(fticr, file = 1)
    spctr <- spectra(first_file)
    checkTrue(length(spctr) == 1)
    mzs <- unname(mz(spctr[[1]]))
    ints <- unname(intensity(spctr[[1]]))
    feats1 <- do_findPeaks_MSW(mz = mzs[10000:20000],
                               int = ints[10000:20000],
                               snthresh = 100)
    feats2 <- do_findPeaks_MSW(mz = mzs[10000:20000],
                               int = ints[10000:20000],
                               snthresh = 50)
    checkTrue(nrow(feats2) > nrow(feats1))
}

test_findChromPeaks_MSW <- function() {
    ## library(MSnbase)
    ## od <- readMSData(mzf, mode = "onDisk")
    od <- microtofq_od
    ## Restrict to first spectrum
    od1 <- od[1]
    sp1 <- od[[1]]
    res_1 <- do_findPeaks_MSW(mz = mz(sp1), int = intensity(sp1))
    mp <- MSWParam()
    checkException(findChromPeaks(od1, param = mp, msLevel = 2))
    res_2 <- findChromPeaks(od1, param = mp)
    checkEquals(res_1, chromPeaks(res_2)[, colnames(res_1), drop = FALSE])
    ## Changing settings.
    snthresh(mp) <- 1
    nearbyPeak(mp) <- FALSE
    res_1 <- do_findPeaks_MSW(mz = mz(sp1), int = intensity(sp1),
                              snthresh = 1, nearbyPeak = FALSE)
    res_2 <- findChromPeaks(od1, param = mp, return.type = "list")
    checkEquals(res_1, res_2[[1]][, colnames(res_1)])
    peakThr(mp) <- 200
    res_1 <- do_findPeaks_MSW(mz = mz(sp1), int = intensity(sp1),
                              snthresh = 1, nearbyPeak = FALSE,
                              peakThr = 200)
    res_2 <- findChromPeaks(od1, param = mp, return.type = "list")
    checkEquals(res_1, res_2[[1]][, colnames(res_1)])
    addParams(mp) <- list(forder = 2)
    res_3 <- do_findPeaks_MSW(mz = mz(sp1), int = intensity(sp1),
                              snthresh = 1, nearbyPeak = FALSE,
                              peakThr = 200, forder = 2)
    res_4 <- findChromPeaks(od1, param = mp, return.type = "list")
    checkEquals(res_3, res_4[[1]][, colnames(res_3)])
    addParams(mp) <- list(forder = 2, dorder = 1)
    res_3 <- do_findPeaks_MSW(mz = mz(sp1), int = intensity(sp1),
                              snthresh = 1, nearbyPeak = FALSE,
                              peakThr = 200, forder = 2, dorder = 1)
    res_4 <- findChromPeaks(od1, param = mp, return.type = "list")
    checkEquals(res_3, res_4[[1]][, colnames(res_3)])

    ## Compare old vs new:
    checkEquals(chromPeaks(fticr_xod)[, -ncol(chromPeaks(fticr_xod))],
                peaks(fticr_xs))
}


############################################################
## Test the implementation of the "do" function
dontrun_test_do_findPeaks_MSW_impl <- function() {

    library(xcms)
    library(RUnit)

    library(faahKO)
    fs <- c(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
            system.file('cdf/KO/ko16.CDF', package = "faahKO"),
            system.file('cdf/KO/ko18.CDF', package = "faahKO"),
            system.file('cdf/KO/ko19.CDF', package = "faahKO"))
    i <- 1
    for (i in 1:length(fs)) {
        cat("============================================================\n")
        cat("|   file", i, ":", fs[i], "\n")
        cat("------------------------------------------------------------")

        xr <- xcmsRaw(fs[i])
        mz <- xr@env$mz
        int <- xr@env$intensity

        ## Default
        a <- findPeaks.MSW(xr)
        b <- xcms:::.MSW_orig(int, mz)
        checkEquals(a@.Data, b)

        ## Slightly modified one:
        c <- xcms:::.MSW(int, mz)
        checkEquals(a@.Data, c)
    }
}

dontrun_benchmark_MSW <- function() {
    ## Actually, we don't expect much of a speed improvement here.
    library(microbenchmark)
    microbenchmark(xcms:::.MSW(xraw@env$intensity, xraw@env$mz),
                   xcms:::.MSW_orig(xraw@env$intensity, xraw@env$mz),
                   times = 4)
}
