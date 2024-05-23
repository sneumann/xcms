xmse_full <- loadXcmsData("xmse")

test_that("storeResults,PlainTextParam,MsExperiment works", {
    pth <- file.path(tempdir(), "test")
    param <- PlainTextParam(path = pth)
    param2 <- PlainTextParam()
    expect_false(is.null(param2))
    expect_error(new("PlainTextParam", path = c(tempdir(), tempdir())))
    mse <- filterMzRange(mse, c(200, 500))
    storeResults(mse, param = param)
    expect_true(dir.exists(pth))
    expect_true(file.exists(file.path(param@path, "sample_data.txt")))
    expect_true(file.exists(file.path(param@path, "spectra_files.txt")))
    expect_true(file.exists(file.path(param@path, "spectra_processing_queue.json")))
})

test_that("storeResults,PlainTextParam,XcmsExperiment works", {
    pth = file.path(tempdir(), "test")
    param <- PlainTextParam(path = pth)
    param2 <- PlainTextParam()
    expect_false(is.null(param2))
    xmse_full <- filterMzRange(xmse_full, c(200, 500))
    storeResults(xmse_full, param = param)
    expect_true(dir.exists(pth))
    expect_true(file.exists(file.path(param@path, "sample_data.txt")))
    expect_true(file.exists(file.path(param@path, "spectra_files.txt")))
    expect_true(file.exists(file.path(param@path, "spectra_processing_queue.json")))
    expect_true(file.exists(file.path(param@path, "process_history.json")))
    expect_true(file.exists(file.path(param@path, "chrom_peaks.txt")))
    expect_true(file.exists(file.path(param@path, "chrom_peak_data.txt")))
    expect_true(file.exists(file.path(param@path, "rtime_adjusted.txt")))
    expect_true(file.exists(file.path(param@path, "feature_definitions.txt")))
    expect_true(file.exists(file.path(param@path, "feature_peak_index.txt")))
    pth = file.path(tempdir(), "test2")
    param <- PlainTextParam(path = pth, spectraExport = FALSE)
    storeResults(xmse_full, param = param)
    expect_false(file.exists(file.path(param@path, "spectra_files.txt")))
})

test_that("loadResults, PlainTextParam works", {
    ## test for MsExperiment object only
    pth = file.path(tempdir(), "test3")
    param <- PlainTextParam(path = pth)
    storeResults(mse, param = param)
    load_mse <- loadResults(object = MsExperiment(), param)
    expect_true(inherits(load_mse, "MsExperiment"))
    expect_equal(mse, load_mse)

    ## test for XcmsExperiment object
    pth = file.path(tempdir(), "test4")
    param <- PlainTextParam(path = pth)
    storeResults(xmse_full, param = param)
    load_xmse <- loadResults(object = XcmsExperiment(), param)
    expect_true(inherits(load_xmse, "XcmsExperiment"))
    expect_equal(xmse_full, load_xmse)
    expect_equal(processHistory(xmse_full), processHistory(load_xmse)) #fail why ?
    expect_equal(xmse_full@featureDefinitions,
                 load_xmse@featureDefinitions)
    expect_equal(adjustedRtime(xmse_full), adjustedRtime(load_xmse))
    expect_equal(xmse_full, load_xmse)
    # not sure how to check  for `spectraFilePath`
    })


