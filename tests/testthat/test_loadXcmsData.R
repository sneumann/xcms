test_that("loadXcmsData works", {
    expect_error(loadXcmsData("some"), "one of")

    a <- loadXcmsData("xdata")
    expect_s4_class(a, "XCMSnExp")
    expect_true(validObject(a))

    b <- loadXcmsData("xmse")
    expect_s4_class(b, "XcmsExperiment")
    expect_true(validObject(b))

    res <- loadXcmsData("faahko_sub")
    expect_s4_class(res, "XCMSnExp")
    expect_true(validObject(res))
    expect_false(hasAdjustedRtime(res))
    expect_false(hasFeatures(res))

    res <- loadXcmsData("faahko_sub2")
    expect_s4_class(res, "XcmsExperiment")
    expect_true(validObject(res))
    expect_true(hasChromPeaks(res))
    expect_false(hasAdjustedRtime(res))
    expect_false(hasFeatures(res))
})
