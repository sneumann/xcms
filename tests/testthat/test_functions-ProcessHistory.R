test_that("ProcessHistory constructor and class works", {
    skip_on_os(os = "windows", arch = "i386")

    ph <- ProcessHistory()
    expect_true(inherits(ph, "ProcessHistory"))

    ph@type <- "FAILURE"
    expect_error(validObject(ph))
    expect_error(ProcessHistory(type = "any"))

    ## Accessor methods.
    expect_error(processType(ph) <- "OOOO")
    processType(ph) <- .PROCSTEP.UNKNOWN
    expect_equal(processType(ph), .PROCSTEP.UNKNOWN)
    expect_error(processDate(ph) <- c("a", "b"))
    processDate(ph) <- "A"
    expect_equal(processDate(ph), "A")
    expect_error(processInfo(ph) <- c("a", "b"))
    processInfo(ph) <- "B"
    expect_equal(processInfo(ph), "B")
    fileIndex(ph) <- 1:3
    expect_equal(fileIndex(ph), 1:3)
})

test_that("XProcessHistory works", {
    skip_on_os(os = "windows", arch = "i386")

    ph <- XProcessHistory()
    expect_true(is(ph, "XProcessHistory"))
    expect_true(inherits(ph, "ProcessHistory"))

    ph <- XProcessHistory(info = "some info",
                          type = .PROCSTEP.PEAK.DETECTION)
    expect_equal(ph@info, "some info")
    expect_equal(ph@type, .PROCSTEP.PEAK.DETECTION)

    ph@type <- "other"
    expect_error(validObject(ph))

    ph <- XProcessHistory(info = "some info",
                          type = .PROCSTEP.PEAK.DETECTION,
                          param = CentWaveParam())

    expect_true(is(ph@param, "CentWaveParam"))
    expect_true(is(processParam(ph), "CentWaveParam"))
})

test_that("GenericProcessHistory works", {
    skip_on_os(os = "windows", arch = "i386")

    xs <- list()
    xs <- c(xs, GenericProcessHistory(fun = "mean"))
    xs <- c(xs, GenericProcessHistory(fun = "median"))
    xs <- c(xs, GenericProcessHistory(fun = "mean"))
    expect_true(length(xs) == 3)
    expect_equal(unlist(lapply(xs, function(z) processParam(z)@fun)),
                 c("mean", "median", "mean"))
    xs <- dropGenericProcessHistoryList(xs, fun = "mean")
    expect_true(length(xs) == 1)
    expect_equal(processParam(xs[[1]])@fun, "median")
})

test_that(".process_history_subset_samples works", {
    skip_on_os(os = "windows", arch = "i386")

    ph <- list(XProcessHistory(CentWaveParam()))
    expect_equal(.process_history_subset_samples(ph), ph)
    ph <- list(XProcessHistory(CentWaveParam(), fileIndex = 1:10))
    res <- .process_history_subset_samples(ph, 1:4)
    expect_equal(res[[1]]@fileIndex, 1:4)
    expect_true(is(res[[1]]@param, "CentWaveParam"))
    ph <- list(XProcessHistory(CentWaveParam(), fileIndex = 1:7),
               XProcessHistory(PeakDensityParam(
                          sampleGroups = c("a", "b", "a", "c", "c", "a", "b"))))
    res <- .process_history_subset_samples(ph, c(2, 4, 1, 2))
    expect_equal(res[[1]]@fileIndex, 1:4)
    expect_equal(res[[2]]@fileIndex, 1:4)
    expect_equal(res[[2]]@param@sampleGroups, c("b", "c", "a", "b"))
    res <- .process_history_subset_samples(ph, 5)
    expect_equal(res[[1]]@fileIndex, 1)
    expect_equal(res[[2]]@fileIndex, 1)
    expect_equal(res[[2]]@param@sampleGroups, c("c"))
})
