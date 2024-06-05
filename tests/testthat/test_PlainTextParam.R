library(xcms)
library(testthat)
library(Spectra)
xmse_full <- loadXcmsData("xmse")

s <- spectra(xmse_full)
b <- s@backend

test_that("storeResults,loadResults,PlainTextParam,MsBackendMzR works", {
    pth <- file.path(tempdir(), "test")
    param <- PlainTextParam(path = pth)
    storeResults(b, param = param)
    expect_true(dir.exists(pth))
    expect_true(file.exists(file.path(param@path, "backend_data.txt")))
    ## Loading data again
    b2 <- loadResults(object = MsBackendMzR(), param)
    expect_true(inherits(b2, "MsBackendMzR"))
    b <- dropNaSpectraVariables(b) #the function does this to be robust, is it a problem ? i should mention it in the doc
    expect_equal(b@spectraData, b2@spectraData)
    expect_equal(peaksVariables(b), peaksVariables(b2))  # true even without forcing the slot
})

test_that("storeResults,loadResults,PlainTextParam,Spectra works", {
    pth <- file.path(tempdir(), "test1")
    param <- PlainTextParam(path = pth)
    #add processingQueueVariables to test export
    s <- filterMzRange(s, c(200,300))
    storeResults(s, param = param)
    expect_true(dir.exists(pth))
    expect_true(file.exists(file.path(param@path, "backend_data.txt")))
    expect_true(file.exists(file.path(param@path, "spectra_slots.txt")))
    expect_true(file.exists(file.path(param@path, "spectra_processing_queue.json")))
    ## Loading data again
    s2 <- loadResults(object = Spectra(), param)
    expect_true(inherits(s2, "Spectra"))
    s <- dropNaSpectraVariables(s)
    expect_equal(s@processingQueue[[1L]]@ARGS, s2@processingQueue[[1L]]@ARGS)
    expect_equal(s@processingQueueVariables, s2@processingQueueVariables)
    expect_equal(s@processing, s2@processing)
    expect_equal(processingChunkSize(s), processingChunkSize(s2))
    expect_equal(s@backend@spectraData, s2@backend@spectraData)
    expect_equal(rtime(s), rtime(s2))
    expect_no_error(filterRt(s2, c(3000, 3500)))
})

test_that("storeResults,loadResults,PlainTextParam,MsExperiment works", {
    pth <- file.path(tempdir(), "test")
    param <- PlainTextParam(path = pth)
    param2 <- PlainTextParam()
    expect_false(is.null(param2))
    expect_error(new("PlainTextParam", path = c(tempdir(), tempdir())))
    mse <- filterMzRange(mse, c(200, 500))
    storeResults(mse, param = param)
    expect_true(dir.exists(pth))
    expect_true(file.exists(file.path(param@path, "sample_data.txt")))
    expect_true(file.exists(file.path(param@path, "backend_data.txt")))
    expect_true(file.exists(file.path(param@path, "spectra_slots.txt")))
    expect_true(file.exists(file.path(param@path, "spectra_processing_queue.json")))
    ## Loading data again
    load_mse <- loadResults(object = MsExperiment(), param)
    expect_true(inherits(load_mse, "MsExperiment"))
    expect_equal(sampleData(mse), sampleData(load_mse))
    a <- spectra(mse)
    b <- spectra(load_mse)
    expect_equal(a@processingQueue[[1L]]@ARGS, b@processingQueue[[1L]]@ARGS)
    expect_equal(rtime(a), rtime(b))
    expect_no_error(filterRt(load_mse, c(3000, 3500)))
})

test_that("storeResults,loadResults,PlainTextParam,XcmsExperiment works", {
    pth = file.path(tempdir(), "test")
    param <- PlainTextParam(path = pth)
    param2 <- PlainTextParam()
    expect_false(is.null(param2))
    xmse_full <- filterMzRange(xmse_full, c(200, 500))
    storeResults(xmse_full, param = param)
    expect_true(dir.exists(pth))
    expect_true(file.exists(file.path(param@path, "sample_data.txt")))
    expect_true(file.exists(file.path(param@path, "backend_data.txt")))
    expect_true(file.exists(file.path(param@path, "spectra_slots.txt")))
    expect_true(file.exists(file.path(param@path, "spectra_processing_queue.json")))
    expect_true(file.exists(file.path(param@path, "process_history.json")))
    expect_true(file.exists(file.path(param@path, "chrom_peaks.txt")))
    expect_true(file.exists(file.path(param@path, "chrom_peak_data.txt")))
    expect_true(file.exists(file.path(param@path, "rtime_adjusted.txt")))
    expect_true(file.exists(file.path(param@path, "feature_definitions.txt")))
    expect_true(file.exists(file.path(param@path, "feature_peak_index.txt")))
    ##load data again
    load_xmse <- loadResults(object = XcmsExperiment(), param)
    expect_true(inherits(load_xmse, "XcmsExperiment"))
    expect_equal(xmse_full@featureDefinitions,
                 load_xmse@featureDefinitions)
    expect_equal(featureValues(xmse_full), featureValues(load_xmse))
    expect_equal(adjustedRtime(xmse_full), adjustedRtime(load_xmse))
    expect_no_error(filterRt(load_xmse, c(3000, 3500)))
    ## not sure how to check for the processHistory slot
    ## still not sure how to create UT for spectraPath
})


