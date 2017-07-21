## Unit tests for functions in functions-normalization.R

test.fitModel <- function() {
    vals <- featureValues(xod_xgrg)
    dat <- data.frame(injection_idx = 1:length(fileNames(xod_xgrg)))
    fits <- xcms:::fitModel(formula = y ~ injection_idx, y = vals, minVals = 3,
                            data = dat)
    ## Check that we've got NULL for features with less thanb 3 values.
    nulls <- apply(vals, MARGIN = 1, function(z) any(is.na(z)))
    checkEquals(nulls, lengths(fits) == 0)
}
