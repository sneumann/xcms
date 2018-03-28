test_that("fitModel works", {
    vals <- featureValues(xod_xgrg)
    dat <- data.frame(injection_idx = 1:length(fileNames(xod_xgrg)))
    fits <- xcms:::fitModel(formula = y ~ injection_idx, y = vals, minVals = 3,
                            data = dat)
    ## Check that we've got NULL for features with less thanb 3 values.
    nulls <- apply(vals, MARGIN = 1, function(z) any(is.na(z)))
    expect_equal(nulls, lengths(fits) == 0)

    ## Check that robustbase would work
    if (requireNamespace("robustbase", quietly = TRUE))
        fits <- xcms:::fitModel(formula = y ~ injection_idx, y = vals,
                                minVals = 3, data = dat, method = "lmrob")
})
