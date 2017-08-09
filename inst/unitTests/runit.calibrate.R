test_calibrate_XCMSnExp <- function() {
    do_plot <- FALSE
    
    tmp <- filterFile(faahko_xod, file = 1)
    
    ## Check shift calibration.
    mzs <- chromPeaks(tmp)[c(3, 6, 7, 13, 17, 32, 45)]
    mzs_shift <- mzs + 0.0001
    prm <- CalibrantMassParam(mz = mzs_shift, method = "shift")
    res <- calibrate(tmp, prm)
    checkTrue(isCalibrated(res))
    checkEquals(chromPeaks(tmp)[, -1], chromPeaks(res)[, -1])
    checkEquals(chromPeaks(tmp)[, 1] + 0.0001, chromPeaks(res)[, 1])
    diffs <- chromPeaks(res)[, "mz"] - chromPeaks(tmp)[, "mz"]
    X <- chromPeaks(res)[, "mz"]
    if (do_plot)
        plot(X, diffs)

    ## Check linear.
    mzs_lin <- mzs + 0.00005 + mzs * 0.000002
    max_dif <- max(mzs_lin - mzs)
    prm <- CalibrantMassParam(mz = mzs_lin, method = "linear", mzabs = max_dif)
    res <- calibrate(tmp, prm)
    checkTrue(isCalibrated(res))
    diffs <- chromPeaks(res)[, "mz"] - chromPeaks(tmp)[, "mz"]
    X <- chromPeaks(res)[, "mz"]
    if (do_plot)
        plot(X, diffs)
    res_lm <- lm(diffs ~ X)
    checkEquals(unname(coefficients(res_lm)[1]), 0.00005, tolerance = 1e-5)
    checkEquals(unname(coefficients(res_lm)[2]), 0.000002, tolerance = 1e-5)

    ## edgeshift
    prm <- CalibrantMassParam(mz = mzs_lin, method = "edgeshift",
                              mzabs = max_dif)    
    res <- calibrate(tmp, prm)
    checkTrue(isCalibrated(res))
    diffs <- chromPeaks(res)[, "mz"] - chromPeaks(tmp)[, "mz"]
    X <- chromPeaks(res)[, "mz"]
    if (do_plot)
        plot(X, diffs)
    mz_sorted <- chromPeaks(tmp)[, "mz"]
    ## Diff has to be constant before and after the linear range.
    lower_idx <- which(chromPeaks(tmp)[, "mz"] < min(mzs))
    checkTrue(all(diffs[lower_idx] == diffs[lower_idx][1]))
    upper_idx <- which(chromPeaks(tmp)[, "mz"] > max(mzs))
    checkTrue(all(diffs[upper_idx] == diffs[upper_idx][1]))
    lin_idx <- 1:length(diffs)
    lin_idx <- lin_idx[!(lin_idx %in% lower_idx)]
    lin_idx <- lin_idx[!(lin_idx %in% upper_idx)]
    lin_mod <- lm(diffs[lin_idx] ~ X[lin_idx])
    checkEquals(unname(coefficients(lin_mod)[1]), 0.00005, tolerance = 1e-5)
    checkEquals(unname(coefficients(lin_mod)[2]), 0.000002, tolerance = 1e-5)
    
    ## Test with a single mass, fall back to shift.
    prm <- CalibrantMassParam(mz = mzs_lin[1], method = "edgeshift",
                              mzabs = max_dif)    
    res <- calibrate(tmp, prm)
    diffs <- chromPeaks(res)[, "mz"] - chromPeaks(tmp)[, "mz"]
    min_diff <- min(abs(chromPeaks(tmp)[, "mz"] - mzs_lin[1]))
    checkEquals(diffs, rep(min_diff, length(diffs)))

    ## Check errors.
    checkException(calibrate(tmp, 4))
    checkException(calibrate(tmp, CalibrantMassParam(mz = list(mzs, mzs))))    
}

