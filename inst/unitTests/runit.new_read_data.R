############################################################
## Testing the new import functions. These might replace the
## xcmsSource, loadRaw etc methods. The new functions use the
## mzR::openMSFile, mzR::peaks etc functions from the mzR
## package, which are newer than the ones used by loadRaw.

############################################################
## Compare the readRawData results to the loadRaw.
test_compare_readRawData <- function() {

    cdf_file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")

    ## loadRaw
    lr_res <- loadRaw(xcmsSource(cdf_file))
    rr_res <- xcms:::readRawData(cdf_file)
    checkEquals(lr_res, rr_res[names(lr_res)])

    mz_file <- system.file("microtofq/MM8.mzML", package = "msdata")
    lr_res <- loadRaw(xcmsSource(mz_file))
    rr_res <- xcms:::readRawData(mz_file)
    checkEquals(lr_res, rr_res[names(lr_res)])

    ## Check readRawData with and without dropEmptyScans:
    msnfile <- system.file("microtofq/MSMSpos20_6.mzML", package = "msdata")
    res <- loadRaw(xcmsSource(msnfile), includeMSn = TRUE)
    res_2 <- xcms:::readRawData(msnfile, includeMSn = TRUE)
    checkEquals(res, res_2)
    res_2 <- xcms:::readRawData(msnfile, includeMSn = TRUE,
                                dropEmptyScans = FALSE)
    ## Now I expect to have more data:
    checkTrue(length(res_2$MSn$precursorIntensity) >
              length(res$MSn$precursorIntensity))
    checkTrue(length(res_2$MSn$precursorIntensity) == 1612)
    checkTrue(length(res$MSn$precursorIntensity) == 1121)
    ## Now, the difference is supposed to represent spectra without peaks:
    empties <- res_2$MSn$peaksCount == 0
    checkTrue(sum(empties) == (1612 - 1121))

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

## Testing also the modified xcmsSource, loadRaw etc.
test_evaluate_xcmsSource <- function() {
    library(msdata)
    mz_file <- system.file("microtofq/MM8.mzML", package = "msdata")
    src <- xcms:::xcmsSource(mz_file)
    checkTrue(is(src, "pwizSource"))
    tmp <- loadRaw(src)
    checkEquals(names(tmp), c("rt", "acquisitionNum", "tic", "scanindex",
                              "mz", "intensity", "polarity"))

    cdf_file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    src <- xcms:::xcmsSource(cdf_file)
    checkTrue(is(src, "netCdfSource"))
    tmp <- loadRaw(src)
    checkEquals(names(tmp), c("rt", "acquisitionNum", "tic", "scanindex",
                              "mz", "intensity", "polarity"))

    ## MSn:
    mzdatapath <- system.file("iontrap", package = "msdata")
    mzdatafiles <- list.files(mzdatapath, pattern="extracted.mzData",
                              recursive = TRUE, full.names = TRUE)
    src <- xcms:::xcmsSource(mzdatafiles[1])
    tmp <- loadRaw(src, includeMSn = TRUE)

    ## OLD code:
    rid <- mzR:::rampOpen(mzdatafiles[1])
    rawdata <- mzR:::rampRawData(rid)
    rawdata$MSn <- mzR:::rampRawDataMSn(rid)
    mzR:::rampClose(rid)
    rm(rid)
    ## Ramp does not read polarity!
    tmp$polarity <- rawdata$polarity
    checkEquals(rawdata, tmp)

    ## Next example:
    msnfile <- system.file("microtofq/MSMSpos20_6.mzML", package = "msdata")
    src <- xcms:::xcmsSource(msnfile)
    tmp <- loadRaw(src, includeMSn = TRUE)
    ## checkTrue(all(tmp$polarity == 1))
    ## OLD code:
    rid <- mzR:::rampOpen(msnfile)
    rawdata <- mzR:::rampRawData(rid)
    rawdata$MSn <- mzR:::rampRawDataMSn(rid)
    mzR:::rampClose(rid)
    rm(rid)
    rawdata$polarity <- tmp$polarity
    checkEquals(rawdata, tmp)
}

dontrun_use_MSnbase <- function() {
    ## Check how many spectra MSnbase reads:
    library(MSnbase)
    msnfile <- system.file("microtofq/MSMSpos20_6.mzML", package = "msdata")
    ms2 <- readMSData(msnfile, msLevel. = 2)
    ## Reads 1612 MS2 spectra.
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

dontrun_evaluate_memory <- function() {
    ## Increased memory demand was observed using the original loadRaw compared
    ## to the newer function. Problem is that can not be easily quantified.
    ## Running feature detection in parallel in the same file.
    tmpDir <- tempdir()
    for (i in 1:40)
        file.copy(system.file("cdf/KO/ko15.CDF", package = "faahKO"),
                  to = paste0(tmpDir, "/", i, ".CDF"))
    fls <- dir(tmpDir, pattern = "CDF", full.names = TRUE)
    ## use the "original" code
    useOriginalCode(TRUE)
    ## 6.7 GB mem free. 1.6 swap used
    system.time(
        xs <- xcmsSet(fls)
    )
    useOriginalCode(FALSE)
    system.time(
        xs <- xcmsSet(fls)
    )
    ## Nothing abnormal observed; but 15 seconds faster.
    tmpDir <- tempdir()
    for (i in 1:40)
        file.copy("../../local_data/mzML-files/130616_10004685_PC_POS.mzML",
                  to = paste0(tmpDir, "/", i, ".mzML"))
    ##
    fls <- dir(tmpDir, pattern = "mzML", full.names = TRUE)
    useOriginalCode(TRUE)
    system.time(
        xs <- xcmsSet(fls, method = "centWave")
    ) ## 370 seconds. Lots of R processes open even after calculation finished;
    ## is BiocParallel not closing the connections properly?
    useOriginalCode(FALSE)
    system.time(
        xs <- xcmsSet(fls, method = "centWave")
    ) ## 380 seconds. There seem to be still some (2) R's open.
}
