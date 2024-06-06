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
    b <- dropNaSpectraVariables(b)
    expect_equal(b@spectraData, b2@spectraData)
    expect_equal(peaksVariables(b), peaksVariables(b2))
    expect_equal(mz(b[1:20]), mz(b2[1:20]))

    ## Check the spectraPath parameter.
    bp <- dataStorageBasePath(b)
    ## manually change dataStorage path of backend
    sd <- read.table(file.path(param@path, "backend_data.txt"), header = TRUE)
    sd$dataStorage <- sub("faahKO", "other", sd$dataStorage)
    write.table(sd, file = file.path(param@path, "backend_data.txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    A <- loadResults(MsBackendMzR(), param)
    expect_error(validObject(A), "invalid class")
    A <- loadResults(MsBackendMzR(), param, spectraPath = bp)
    expect_true(validObject(A))

    param <- PlainTextParam(tempdir())
    expect_error(loadResults(MsBackendMzR(), param), "No 'backend_data")
})

test_that("storeResults,loadResults,PlainTextParam,Spectra works", {
    pth <- file.path(tempdir(), "test1")
    param <- PlainTextParam(path = pth)
    ## add processingQueueVariables to test export
    s@processingQueueVariables <- c(s@processingQueueVariables, "rtime")
    s <- filterMzRange(s, c(200,300))
    s <- filterRt(s, c(3000, 3500)) # to ensure subsetted object would work
    storeResults(s, param = param)
    expect_true(dir.exists(pth))
    expect_true(file.exists(file.path(param@path, "backend_data.txt")))
    expect_true(file.exists(file.path(param@path, "spectra_slots.txt")))
    expect_true(file.exists(file.path(param@path,
                                      "spectra_processing_queue.json")))
    ## Loading data again
    s2 <- loadResults(object = Spectra(), param)
    expect_true(inherits(s2, "Spectra"))
    expect_true(inherits(s2@backend, "MsBackendMzR"))
    s <- dropNaSpectraVariables(s)
    expect_equal(length(s@processingQueue), length(s2@processingQueue))
    expect_equal(s@processingQueue[[1L]]@ARGS, s2@processingQueue[[1L]]@ARGS)
    expect_equal(s@processingQueueVariables, s2@processingQueueVariables)
    expect_equal(s@processing, s2@processing)
    expect_equal(processingChunkSize(s), processingChunkSize(s2))
    expect_equal(s@backend@spectraData, s2@backend@spectraData)
    expect_equal(rtime(s), rtime(s2))
    expect_equal(mz(s[1:10]), mz(s2[1:10])) # that does NOT work without calling functions directly with MsCoreUtils:: in Spectra.
    expect_no_error(filterRt(s2, c(3000, 3500)))

    ## Check the spectraPath parameter.
    ## Changing the path in the MsBackendMzR to simulate moving the exported data
    bp <- dataStorageBasePath(s)
    sd <- read.table(file.path(param@path, "backend_data.txt"), header = TRUE)
    sd$dataStorage <- sub("faahKO", "other", sd$dataStorage)
    write.table(sd, file = file.path(param@path, "backend_data.txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    A <- loadResults(Spectra(), param)
    expect_error(validObject(A@backend), "invalid class")
    A <- loadResults(Spectra(), param, spectraPath = bp)
    expect_true(validObject(A@backend))

    param <- PlainTextParam(file.path(tempdir()))
    expect_error(loadResults(Spectra(), param), "No 'spectra_slots")
})

test_that("storeResults,loadResults,PlainTextParam,MsExperiment works", {
    pth <- file.path(tempdir(), "test3")
    param <- PlainTextParam(path = pth)
    param2 <- PlainTextParam()
    expect_false(is.null(param2))
    expect_error(new("PlainTextParam", path = c(tempdir(), tempdir())))
    tmp <- filterMzRange(mse, c(200, 500))
    tmp <- filterRt(tmp, c(3000, 3500))
    storeResults(tmp, param = param)
    expect_true(dir.exists(pth))
    expect_true(file.exists(file.path(param@path, "sample_data.txt")))
    expect_true(file.exists(file.path(param@path, "backend_data.txt")))
    expect_true(file.exists(file.path(param@path, "spectra_slots.txt")))
    expect_true(file.exists(file.path(param@path, "spectra_processing_queue.json")))
    ## Loading data again
    load_mse <- loadResults(object = MsExperiment(), param)
    expect_true(inherits(load_mse, "MsExperiment"))
    expect_equal(sampleData(tmp), sampleData(load_mse))
    a <- spectra(tmp)
    b <- spectra(load_mse)
    expect_equal(length(a@processingQueue), length(b@processingQueue))
    expect_equal(a@processingQueue[[1L]]@ARGS, b@processingQueue[[1L]]@ARGS)
    expect_equal(rtime(a), rtime(b))
    expect_no_error(filterRt(load_mse, c(3000, 3500)))

    ## Check the spectraPath parameter.
    bp <- dataStorageBasePath(tmp@spectra)
    ## manually change dataStorage path of backend
    sd <- read.table(file.path(param@path, "backend_data.txt"), header = TRUE)
    sd$dataStorage <- sub("faahKO", "other", sd$dataStorage)
    write.table(sd, file = file.path(param@path, "backend_data.txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    A <- loadResults(MsExperiment(), param)
    expect_error(validObject(spectra(A)@backend), "invalid class")
    A <- loadResults(MsExperiment(), param, spectraPath = bp)
    expect_true(validObject(spectra(A)@backend))

    param <- PlainTextParam(tempdir())
    expect_error(loadResults(MsExperiment(), param), "No 'sample_data")
})

test_that("storeResults,loadResults,PlainTextParam,XcmsExperiment works", {
    pth = file.path(tempdir(), "test4")
    param <- PlainTextParam(path = pth)
    param2 <- PlainTextParam()
    expect_false(is.null(param2))
    tmp <- filterMzRange(xmse_full, c(200, 500))
    tmp <- filterRt(tmp, c(3000, 4000))
    storeResults(tmp, param = param)
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

    ## load data again
    load_xmse <- loadResults(object = XcmsExperiment(), param)
    expect_true(inherits(load_xmse, "XcmsExperiment"))
    expect_equal(tmp@featureDefinitions,
                 load_xmse@featureDefinitions)
    expect_equal(featureValues(tmp), featureValues(load_xmse))
    expect_equal(adjustedRtime(tmp), adjustedRtime(load_xmse))
    expect_no_error(filterRt(load_xmse, c(3000, 3500)))
    expect_equal(tmp@chromPeaks, load_xmse@chromPeaks)
    expect_equal(tmp@chromPeakData, load_xmse@chromPeakData)
    expect_equal(tmp@sampleData, load_xmse@sampleData)
    expect_equal(length(tmp@processHistory), length(load_xmse@processHistory))
    expect_equal(tmp@processHistory[[1L]], load_xmse@processHistory[[1L]])
    expect_equal(tmp@processHistory[[2L]], load_xmse@processHistory[[2L]])
    expect_equal(tmp@processHistory[[3L]], load_xmse@processHistory[[3L]])
    expect_equal(tmp@processHistory[[4L]], load_xmse@processHistory[[4L]])
    expect_equal(tmp@processHistory[[5L]], load_xmse@processHistory[[5L]])
    ## The 6th param object contains functions for which the comparison fails
    ## because of the name/namespace mentioned. See e.g.
    ## tmp@processHistory[[6]]@param load_xmse@processHistory[[6]]@param

    ## Check the spectraPath parameter.
    bp <- dataStorageBasePath(tmp@spectra)
    ## manually change dataStorage path of backend
    sd <- read.table(file.path(param@path, "backend_data.txt"), header = TRUE)
    sd$dataStorage <- sub("faahKO", "other", sd$dataStorage)
    write.table(sd, file = file.path(param@path, "backend_data.txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    expect_error(loadResults(XcmsExperiment(), param), "invalid class")
    A <- loadResults(XcmsExperiment(), param, spectraPath = bp)
    expect_true(validObject(spectra(A)@backend))

    param <- PlainTextParam(tempdir())
    expect_error(loadResults(MsExperiment(), param), "No 'sample_data")
})
