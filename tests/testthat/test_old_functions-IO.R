test_that("readRawData works", {
    cdf_file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")

    ## loadRaw
    lr_res <- loadRaw(xcmsSource(cdf_file))
    rr_res <- readRawData(cdf_file)
    expect_equal(lr_res, rr_res[names(lr_res)])

    mz_file <- system.file("microtofq/MM8.mzML", package = "msdata")
    lr_res <- loadRaw(xcmsSource(mz_file))
    rr_res <- readRawData(mz_file)
    expect_equal(lr_res, rr_res[names(lr_res)])

    ## Check readRawData with and without dropEmptyScans:
    msnfile <- system.file("microtofq/MSMSpos20_6.mzML", package = "msdata")
    res <- loadRaw(xcmsSource(msnfile), includeMSn = TRUE)
    res_2 <- readRawData(msnfile, includeMSn = TRUE)
    expect_equal(res, res_2)
    res_2 <- readRawData(msnfile, includeMSn = TRUE,
                         dropEmptyScans = FALSE)
    ## Now I expect to have more data:
    expect_true(length(res_2$MSn$precursorIntensity) >
                length(res$MSn$precursorIntensity))
    expect_true(length(res_2$MSn$precursorIntensity) == 1612)
    expect_true(length(res$MSn$precursorIntensity) == 1121)
    ## Now, the difference is supposed to represent spectra without peaks:
    empties <- res_2$MSn$peaksCount == 0
    expect_true(sum(empties) == (1612 - 1121))

    ## ## Check also MSn level import. Note: this can not be run automatically
    ## ## because the mzML file is gzipped; xcmsSource does not support gz
    ## ## input! The test was performed by manually unzipping the file and
    ## ## running on the mzML file.
    ## msn <- msdata::proteomics(full.names = TRUE, pattern = "TMT_Erwinia")
    ## lr_res <- loadRaw(xcmsSource(msn))  ## Can not read .gz files!!!
    ## rr_res <- xcms:::readRawData(msn)
    ## expect_equal(lr_res, rr_res[names(lr_res)])
    ## ## Include MSn
    ## lr_res <- loadRaw(xcmsSource(msn), includeMSn = TRUE)
    ## rr_res <- xcms:::readRawData(msn, includeMSn = TRUE)
    ## expect_equal(lr_res, rr_res[names(lr_res)])
    ## ## Rocks!
})
