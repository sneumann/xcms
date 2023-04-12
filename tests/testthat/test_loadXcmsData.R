test_that("loadXcmsData works", {
    expect_error(loadXcmsData("some"), "one of")

    a <- loadXcmsData("xdata")
    expect_s4_class(a, "XCMSnExp")
    expect_true(validObject(a))

    b <- loadXcmsData("xmse")
    expect_s4_class(b, "XcmsExperiment")
    expect_true(validObject(b))
})
