test_that("fitModel works", {
    vals <- featureValues(xod_xgrg)
    dat <- data.frame(injection_idx = 1:length(fileNames(xod_xgrg)))
    fits <- xcms:::rowFitModel(formula = y ~ injection_idx, y = vals,
                               minVals = 3, data = dat)
    ## Check that we've got NA for features with less than 3 values.
    nas <- apply(vals, MARGIN = 1, function(z) any(is.na(z)))
    expect_equal(nas, vapply(fits, function(z) any(is.na(z)), logical(1)))

    ## Check that robustbase would work
    if (requireNamespace("robustbase", quietly = TRUE))
        fits <- xcms:::rowFitModel(formula = y ~ injection_idx, y = vals,
                                   minVals = 3, data = dat, method = "lmrob")
})


test_that("model adjustment with batch works", {
    ## Here we test that linear model based adjustment with a batch is
    ## working.
    y <- c(2, 3, 2.7, 3.5, 3.8, 4.6, 5.9, 8, 4, 5.1, 5.6, 6.8, 7.1, 8.1, 8.9,
           9.3)
    inj_idx <- 1:length(y)
    btch <- c(rep("a", 8), rep("b", 8))
    dta <- data.frame(inj_idx = inj_idx, btch = btch)
    plot(inj_idx, y)
    lmod <- lm(y ~ inj_idx + btch)
    prd <- predict(lmod, newdata = data.frame(y = y, inj_idx = inj_idx,
                                              btch = btch))
    ## y_new <- y - prd + mean(y)
    y_new <- y - prd + mean(lmod$fitted.values + lmod$residuals)
    
    points(inj_idx, y_new, pch = 16, col = "grey")

    expect_equal(mean(y), mean(y_new))

    mdl <- fitModel(y ~ inj_idx + btch, y = y,
                    data = data.frame(inj_idx = inj_idx, btch = btch))
    expect_equal(mdl$coefficients, lmod$coefficients)
    res <- applyModelAdjustment(y, data = data.frame(inj_idx = inj_idx,
                                                     btch = btch),
                                lmod = mdl)
    expect_equal(res, unname(y_new))

    ## Check if we have only NAs in one batch
    y[btch == "a"] <- NA
    mdl <- xcms:::fitModel(y ~ inj_idx + btch, y = y, data = dta)
    expect_true(is.na(mdl))
    ## plot(x = dta$inj_idx, y_2)
    ## abline(sum(mdl$coefficients[c(1, 3)]), mdl$coefficients[2])
    ## ## Adjust.
    ## res <- xcms:::applyModelAdjustment(y = y_2, data = dta, lmod = mdl)
    ## points(dta$inj_idx, res, pch = 16)
    ## mean(res, na.rm = TRUE)
    ## mean(y_2, na.rm = TRUE)
})

test_that("fitModel works on matrix and vector", {
    y <- c(2, 3, 2.7, 3.5, 3.8, 4.6, 5.9, 8, 4, 5.1, 5.6, 6.8, 7.1)
    inj_idx <- 1:length(y)
    dta <- data.frame(inj_idx = inj_idx)

    expect_error(xcms:::fitModel())
    expect_error(xcms:::fitModel(formula = ~ y))
    expect_error(xcms:::fitModel(formula = ~ y, data = dta))
    expect_error(xcms:::fitModel(formula = ~ inj_idx, data = dta, y = y))
    expect_error(xcms:::fitModel(y ~ inj_idx, data = dta, y = y, method = "adfd"))
    expect_error(xcms:::fitModel(y ~ inj_idx, data = dta, y = y, weights = 3))
    
    res <- xcms:::fitModel(y ~ inj_idx, data = dta, y = y)
    expect_equal(res$coefficients, lm(y ~ inj_idx)$coefficients)
    rres <- xcms:::fitModel(y ~ inj_idx, data = dta, y = y, method = "lmrob")
    expect_true(all(res$coefficients != rres$coefficients))

    ## Test with weights
    res2 <- xcms:::fitModel(y ~ inj_idx, data = dta, y = y,
                            weights = abs(rnorm(length(y))))
    expect_true(all(res$coefficients != res2$coefficients))
    rres2 <- xcms:::fitModel(y ~ inj_idx, data = dta, y = y, method = "lmrob",
                             weights = abs(rnorm(length(y))))
    expect_true(all(rres$coefficients != rres2$coefficients))
    
    ymat <- matrix(rep(y, 5), nrow = 5, byrow = TRUE)
    res_3 <- fitModel(y ~ inj_idx, data = dta, y = ymat)
    expect_equal(res_3$coefficients, res$coefficients)

    ymat[2, ] <- ymat[2, ] + 3
    res <- rowFitModel(y ~ inj_idx, data = dta, y = ymat)
    expect_true(length(res) == nrow(ymat))
    expect_equal(res[[1]]$coefficients, res_3$coefficients)
    slps <- vapply(res, function(z) z$coefficients[2], numeric(1))
    ints <- vapply(res, function(z) z$coefficients[1], numeric(1))
    expect_equal(unname(slps[1]), unname(slps[2]))
    expect_equal(unname(slps[1]), unname(slps[3]))
    expect_equal(unname(slps[1]), unname(slps[4]))
    expect_equal(unname(slps[1]), unname(slps[5]))
    expect_equal(unname(ints[1]) + 3, unname(ints[2]))

    ## y being a matrix, weights a vector.LLLLLLL

    ## y being a matrix, weights a matrix.
})

test_that("applyModelAdjustment works", {
    y <- c(2, 3, 2.7, 3.5, 3.8, 4.6, 5.9, 8, 4, 5.1, 5.6, 6.8, 7.1)
    inj_idx <- 1:length(y)
    btch <- c(rep("a", 8), rep("b", 5))
    dta <- data.frame(inj_idx = inj_idx, batch = btch)

    mdl <- xcms:::fitModel(y ~ ii, data = data.frame(ii = inj_idx), y = y)

    ## A single model, single row.
    res <- xcms:::applyModelAdjustment(y, data.frame(ii = inj_idx), lmod = mdl)
    plot(inj_idx, y)
    abline(mdl)
    points(inj_idx, res, col = "grey", pch = 16)
    expect_equal(mean(res), mean(y))
    expect_true(lm(res ~ inj_idx)$coefficients[2] < 1e-7)

    ## Model with only batch
    mdl <- xcms:::fitModel(y ~ batch, data = dta, y = y)
    res <- xcms:::applyModelAdjustment(y, dta, lmod = mdl)
    plot(inj_idx, y)
    abline(mdl)
    points(inj_idx, res, col = "grey", pch = 16)
    expect_equal(mean(res), mean(y))
    expect_equal(mean(res[dta$batch == "a"]), mean(res[dta$batch == "b"]))
    
    ## Model with batch-specific slope
    mdl <- xcms:::fitModel(y ~ inj_idx * batch, data = dta, y = y)
    res <- xcms:::applyModelAdjustment(y, dta, lmod = mdl)
    plot(inj_idx, y)
    points(inj_idx, res, col = "grey", pch = 16)
    expect_equal(mean(res), mean(y))
    expect_equal(mean(res[dta$batch == "a"]), mean(res[dta$batch == "b"]))

    ## A single model on a matrix.
    ymat <- matrix(rep(y, 5), nrow = 5, byrow = TRUE)
    ymat[2, ] <- ymat[2, ] + 3
    mdl <- xcms:::fitModel(y ~ ii, data = data.frame(ii = inj_idx), y = y)
    res <- xcms:::applyModelAdjustment(y, data.frame(ii = inj_idx), lmod = mdl)
    resm <- xcms:::applyModelAdjustment(ymat, data.frame(ii = inj_idx),
                                        lmod = mdl)
    expect_equal(resm[1, ], res)
    expect_true(lm(resm[2, ] ~ inj_idx)$coefficients[2] < 1e-7)

    ## multiple models with a matrix.
    mdls <- xcms:::rowFitModel(y ~ ii, data = data.frame(ii = inj_idx),
                               y = ymat)
    resm <- xcms:::applyModelAdjustment(ymat, data.frame(ii = inj_idx),
                                        lmod = mdls)
    expect_equal(resm[1, ], res)
    expect_equal(resm[2, ], res + 3)

    ## Use only some for estimation
    idx <- 1:8

    mdl <- fitModel(y ~ ii, data = data.frame(ii = inj_idx[idx]), y = y[idx])
    res <- applyModelAdjustment(y, data.frame(ii = inj_idx), lmod = mdl)
    plot(inj_idx, y, ylim = c(0, 8))
    abline(mdl)
    points(inj_idx, res, col = "grey", pch = 16)
})

test_that("replaceNaOnEnds works", {
    x <- c(NA, 3, 4, 6, 4, 2, NA, 3, NA, 4, 5, 6, NA)
    expect_equal(replaceNaOnEnds(x), c(3, 3, 4, 6, 4, 2, NA, 3, NA, 4, 5, 6, 6))

    expect_error(replaceNaOnEnds(x, batch = 1:4))

    batch <- c(rep("b", 7), rep("a", 6))
    inji <- c(1:7, 1:6)
    res <- replaceNaOnEnds(x, batch = batch, injIndex = inji)
    expect_equal(res, c(3, 3, 4, 6, 4, 2, 2, 3, NA, 4, 5, 6, 6))

    ## Should also work if randomized.
    idx <- sample(1:length(x))
    x[idx]
    res <- replaceNaOnEnds(x[idx], batch = batch[idx], injIndex = inji[idx])
    expect_equal(res[order(idx)], c(3, 3, 4, 6, 4, 2, 2, 3, NA, 4, 5, 6, 6))

    ## With a matrix...
    xmat <- matrix(rep(x, each = 100), nrow = 100)
    res <- replaceNaOnEnds(xmat, batch = batch, injIndex = inji)
    expect_equal(res[1, ], c(3, 3, 4, 6, 4, 2, 2, 3, NA, 4, 5, 6, 6))
    expect_true(nrow(unique(res)) == 1)
    expect_equal(unique(res)[1, ], c(3, 3, 4, 6, 4, 2, 2, 3, NA, 4, 5, 6, 6))
})
