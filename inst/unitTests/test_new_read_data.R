############################################################
## Testing the new import functions. These might replace the
## xcmsSource, loadRaw etc methods. The new functions use the
## mzR::openMSFile, mzR::peaks etc functions from the mzR
## package, which are newer than the ones used by loadRaw.

############################################################
## Compare the readRawData results to the loadRaw.
test_compare_readRawData <- function() {

    library(faahKO)
    library(xcms)
    library(RUnit)
    cdf_file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")

    ## loadRaw
    lr_res <- loadRaw(xcmsSource(cdf_file))
    rr_res <- xcms:::readRawData(cdf_file)
    checkEquals(lr_res, rr_res[names(lr_res)])

    mz_file <- system.file("microtofq/MM8.mzML", package = "msdata")
    lr_res <- loadRaw(xcmsSource(mz_file))
    rr_res <- xcms:::readRawData(mz_file)
    checkEquals(lr_res, rr_res[names(lr_res)])

    ## ## Check also MSn level import. Note: this can not be run automatically
    ## ## because the mzML file is gzipped; xcmsSource does not support gz
    ## ## input! The test was performed by manually unzipping the file and
    ## ## running on the mzML file.
    ## msn <- msdata::proteomics(full.names = TRUE, pattern = "TMT_Erwinia")
    ## lr_res <- loadRaw(xcmsSource(msn))  ## Can not read .gz files!!!
    ## rr_res <- xcms:::readRawData(msn)
    ## checkEquals(lr_res, rr_res[names(lr_res)])
    ## ## Include MSn
    ## lr_res <- loadRaw(xcmsSource(msn), includeMSn = TRUE)
    ## rr_res <- xcms:::readRawData(msn, includeMSn = TRUE)
    ## checkEquals(lr_res, rr_res[names(lr_res)])
    ## ## Rocks!
}

dontrun_benchmarks <- function() {
    library(microbenchmark)

    cdf_file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")

    ## CDF
    microbenchmark(loadRaw(xcmsSource(cdf_file)),
                   xcms:::readRawData(cdf_file), times = 30)
    ## loadRaw is faster.

    ## mzML
    mz_file <- system.file("microtofq/MM8.mzML", package = "msdata")
    microbenchmark(loadRaw(xcmsSource(mz_file)),
                   xcms:::readRawData(mz_file), times = 30)
    ## readRawData is faster.
}
