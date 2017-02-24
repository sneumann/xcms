############################################################
## Test useOriginalCode

test_useOriginalCode <- function() {

    ## orig <- useOriginalCode()
    ## checkTrue(useOriginalCode())
    ## checkTrue(useOriginalCode(TRUE))
    ## checkTrue(useOriginalCode())
    ## checkTrue(!useOriginalCode(FALSE))
    ## checkTrue(!useOriginalCode())
    ## useOriginalCode(orig)
}

############################################################
## Check if we switch between new and old code correctly:
## The results for binSize = 0.2 should differ between old
## and new (rounding errors).
dontrun_test_matchedFilter_orig_code <- function() {
    require(faahKO)
    fs <- system.file('cdf/KO/ko15.CDF', package = "faahKO")

    xr <- xcmsRaw(fs, profstep = 0)
    mz <- xr@env$mz
    int <- xr@env$intensity
    scantime <- xr@scantime
    scanindex <- xr@scanindex
    valsPerSpect <- diff(c(scanindex, length(mz)))
    step <- 0.2

    orig <- useOriginalCode()
    res_new <- xcms:::do_findChromPeaks_matchedFilter(mz, int,
                                                      scantime,
                                                      valsPerSpect,
                                                      binSize = step)
    useOriginalCode(TRUE)
    res_old <- xcms:::do_findChromPeaks_matchedFilter(mz, int,
                                                      scantime,
                                                      valsPerSpect,
                                                      binSize = step)
    useOriginalCode(FALSE)
    ## Compare
    checkTrue(is.character(all.equal(res_new, res_old)))
    useOriginalCode(orig)
}
