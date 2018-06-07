#' @include DataClasses.R functions-utils.R

#' @title Fit linear model(s) to data
#'
#' @description
#'
#' `fitModel` fits a single linear model to the data which can be passed as a
#' `numeric` vector or a `matrix`.
#'
#' `rowFitModel` fits row-wise linear models to the `matrix` submitted with
#' argument `y`.
#' 
#' @details
#'
#' For `method = "lmrob"`, `fitModel` and `rowFitModel` perform robust
#' regression using the [lmrob()] function with `settings = "KS2014"` and
#' `method = "SMDB"`.
#'
#' `rowFitModel` performs by default parallel fitting of the models
#' based on the global parallel processing settings.
#' 
#' @note
#'
#' Be aware that fits can be unstable if they are based on only few
#' measurements. While `minVals` allows to set a lower threshold for the
#' number of non-NA values in `y` there is no check whether e.g. enough
#' values are available per batch for a model in the form `y ~ idx * batch`.
#' The user is advised to flag or exclude such cases before or after fitting.
#'
#' @param formula `formula` representing the model that should be fitted.
#'
#' @param data `data.frame` containing the explanatory variables in `formula`.
#'
#' @param y `numeric` or `matrix` representing the response variable y in
#'     `formula`. This can be e.g the matrix of feature abundances returned by
#'     [featureValues()].
#'
#' @param minVals `integer(1)` defining the minimum number of values to be
#'     used for the fitting. Model fitting is skipped if less than `minVals`
#'     non missing values are available in `y`, in which case `NA` is returned.
#'
#' @param method `character` defining the method/function to be used for
#'     model fitting. Allowed values are `"lm"` for least squares
#'     regression and `"lmrob"` for robust regression using the
#'     [lmrob()] function.
#'
#' @param BPPARAM optional parameter specifying parallel processing settings.
#'
#' @return `fitModel` returns the fitted linear model. `rowFitModel` a `list`
#'     of linear models, one for each row in `y`, or `NA` for rows with too
#'     few data points or for which the fitting failed.
#'
#' @author Johannes Rainer
#'
#' @noRd
#'
#' @md
#'
#' @seealso [applyModelAdjustment()] for a function to perform linear model
#'     based abundance adjustment.
#' 
#' @rdname model-fitting
#'
#' @examples
#'
#' ## Fitting a simple model to a data vector.
#' y <- c(2, 3, 2.7, 3.5, 3.8, 4.6, 5.9, 8, 4, 5.1, 5.6, 6.8, 7.1)
#' inj_idx <- 1:length(y)
#' dta <- data.frame(inj_idx = inj_idx)
#'
#' plot(inj_idx, y)
#' res <- fitModel(y ~ inj_idx, data = dta, y = y)
#' abline(res)
fitModel <- function(formula, data, y, method = c("lm", "lmrob"), control,
                     minVals = 4, weights = rep(1, length(y))) {
    method <- match.arg(method, c("lm", "lmrob"))
    if (method == "lmrob" & missing(control)) {
        ## Force use of the KS2014 settings in lmrob and increase the
        ## scale-finding iterations to avoid some of the warnings.
        control <- robustbase::lmrob.control("KS2014")
        control$maxit.scale <- 10000
        control$k.max <- 10000
        control$refine.tol <- 1e-7
    }
    if (missing(formula) || !is(formula, "formula"))
        stop("'formula' has to be submitted and should be a formula!")
    if (missing(data) || !is(data, "data.frame"))
        stop("'data' has to be submitted and should be a 'data.frame'!")
    if (missing(y))
        stop("'y' is missing with no default.")
    if (is.matrix(y)) {
        nc <- nrow(y)
        y <- as.numeric(y)
        data <- data[rep(1:nrow(data), each = nc), , drop = FALSE]
    }
    if (length(y) != nrow(data))
        stop("length of 'y' has to match the number of rows of 'data'")
    ## Check that 'data' contains the variables we're looking for.
    vars <- all.vars(formula)
    if (vars[1] != "y")
        stop("'formula' should start with 'y ~'")
    if (!all(vars[-1] %in% colnames(data)))
        stop("All of the variables from 'formula' have to be present in 'data'")
    ## data shouldn't contain a column y.
    if (any(colnames(data) == "y"))
        stop("'data' should not contain a column named 'y'")
    if (length(weights) != length(y))
        stop("'weights' has to have the same length than 'y'")
    force(weights)
    data$y <- y
    ## Check valid measurements:
    nona <- !is.na(data$y)
    res <- NA
    if (sum(nona) >= minVals) {
        if (method == "lm")
            res <- tryCatch(
                ## Note: have to use do.call, otherwise the weights parameter
                ## is not found.
                do.call("lm", args = list(formula = formula, data = data,
                        weights = weights, model = FALSE)),
                ## lm(formula = formula, data = data, weights = weights.,
                ##    model = FALSE),
                error = function(e) {
                    paste0("Failed to fit model: ", e)
                })
        if (method == "lmrob") {
            set.seed(123)
            res <- tryCatch(
                do.call("lmrob", args = list(formula = formula, data = data,
                                             model = FALSE, control = control,
                                             weights = weights)),
                ## robustbase::lmrob(formula = formula, data = data, model = FALSE,
                ##                   control = control, w = weights),
                error = function(e) {
                    paste0("Failed to fit model: ", e)
                })
            ## if (lmeth == "rlm")
            ##     stop("Not yet implemented")
            ## ##     return(MASS::rlm(formula., data = data.))
        }
        if (is.character(res)) {
            warning(res, call. = FALSE)
            res <- NA
        }
    }
    res
}

#' @noRd
#'
#' @md
#' 
#' @rdname model-fitting
rowFitModel <- function(formula, data, y, minVals = 4,
                        method = c("lm", "lmrob"), BPPARAM = bpparam(),
                        weights) {
    method <- match.arg(method, c("lm", "lmrob"))
    if (missing(formula) || !is(formula, "formula"))
        stop("'formula' has to be submitted and should be a formula!")
    if (missing(data) || !is(data, "data.frame"))
        stop("'data' has to be submitted and should be a 'data.frame'!")
    if (missing(y))
        stop("'y' is missing with no default.")
    if (ncol(y) != nrow(data))
        stop("ncol(y) has to match nrow(data)!")
    ## Check that 'data' contains the variables we're looking for.
    vars <- all.vars(formula)
    if (vars[1] != "y")
        stop("'formula' should start with 'y ~'")
    if (!all(vars[-1] %in% colnames(data)))
        stop("All of the variables from 'formula' have to be present in 'data'")
    ## data shouldn't contain a column y.
    if (any(colnames(data) == "y"))
        stop("'data' should not contain a column named 'y'")
    if (!missing(weights) && is.matrix(weights)) {
        if (nrow(weights) != nrow(y))
            stop("If 'weights' is a 'matrix' its number of rows have to match",
                 " the number of rows of 'y'")
    } else weights <- rep(1, ncol(y))
    ## Done with checking.
    ## Subset data to contain only explanatory variables
    data <- data[, vars[-1], drop = FALSE]
    if (is.null(rownames(y)))
        rownames(y) <- 1:nrow(y)
    y <- split.data.frame(y, f = factor(rownames(y), levels = rownames(y)))
    sttngs <- list()
    if (method == "lmrob") {
        ## Force use of the KS2014 settings in lmrob and increase the
        ## scale-finding iterations to avoid some of the warnings.
        sttngs <- robustbase::lmrob.control("KS2014")
        sttngs$maxit.scale <- 10000
        sttngs$k.max <- 10000
        sttngs$refine.tol <- 1e-7
    }
    if (is.matrix(weights))
        bpmapply(FUN = fitModel, y,
                 split.data.frame(weights, f = 1:nrow(weights)),
                 MoreArgs = list(formula = formula, minVals = minVals,
                                 data = data, method = method,
                                 control = sttngs, weights = weights),
                 BPPARAM = BPPARAM)
    else 
        bplapply(y, fitModel, formula = formula, minVals = minVals,
                 data = data, method = method, control = sttngs,
                 weights = weights, BPPARAM = BPPARAM)
}


#' Adjust a single data vector based on the provided model.
#'
#' @noRd
.adjust_with_linear_model <- function(y, data, model) {
    y_new <- y
    if (length(y) != nrow(data))
        stop("length of 'y' has to match the number of rows of 'data'")
    ## Catch problems predicting the value, if e.g. explanatory variables
    ## have additional factor levels in newdata
    preds <- tryCatch(predict(model, newdata = data.frame(y = y, data)),
                      error = function(e) {
                          warning("Failed to adjust value: ", e, call. = FALSE)
                      })
    ## Adjust the drift and add add the mean of the values on which the fit
    ## was performed.
    ## This is in accordance with the code from Wehrens et al 2016:
    ## https://github.com/rwehrens/BatchCorrMetabolomics
    ## I personally would add mean(y) instead to ensure that the mean before
    ## and after normalization is the same.
    if (is.numeric(preds))
        y_new <- y - preds + mean(model$fitted.values + model$residuals)
        ## y_new <- y - preds + mean(y)
    y_new
}


#' @title Model-based abundance adjustment
#'
#' @description
#' 
#' The functions listed here allow to perform linear model-based adjustment of
#' feature abundances similar to [Wehrens 2018]. This comprises simple
#' between-batch differences, injection index dependent signal drift
#' adjustment and combinations thereof.
#' Fitting of linear models is performed with the [fitModel()] or
#' [rowFitModel()] functions.
#'
#' `applyModelAdjustment` adjusts the data in `y` based on the linear
#' model(s) provided with `lmod` (estimated e.g. using [rowFitModel()]). The
#' function performs row-wise adjustment, if `y` is a `matrix`. If `lmod` is
#' a single linear model, each row in `y` is then adjusted with the same model.
#' To adjust each row with its own model, pass a `list` of linear models to
#' `lmod` with length equal to the number of rows in `y`. The list can have
#' also `NULL` elements for rows that should/could not be adjusted.
#'
#' `adjustWithModel` adjusts the provided values based on the linear model
#' specified with `formula`. This function is a convenience function that
#' calls first [rowFitModel()] to estimate the effects to adjust (e.g. on a
#' subset of samples/columns in `y`) and subsequently adjusts the data
#' using the [applyModelAdjustment()].
#'
#' @details
#'
#' For some rows/features values can become negative after adjustment.
#' To avoid this, a constant can be added to the adjusted intensities of
#' such features. Parameter `shiftNegative` allows to specify how this
#' constant is to be determined. For `shiftNegative = "min"`, if one
#' of the adjusted values of a row is `< 0`, the minimum intensity
#' is added to each intensity. This shifting is done on a per-feature basis.
#' Alternatively, the `globalMin` *globally* shifts the complete
#' matrix by the minimum value (if it is negative).
#' 
#' @note
#'
#' Rows in `y` with less than `minValues` non-NA values are returned unadjusted
#' by the `adjustWithModel` function. The function evaluates however only if
#' there are **in total** at least `minValues` available, but not e.g. within
#' a batch for models of the form `y ~ idx * batch`. The user is advised to
#' check and identify such cases.
#'
#' @param y For `applyModelAdjustment`: `numeric` or `matrix` with values that
#'     should be adjusted. For `adjustWithModel`: a `matrix` with values to be
#'     adjusted.
#'
#' @param data `data.frame` with the same explanatory variables as used by the
#'     model fitting.
#' 
#' @param lmod `list` of linear model fits such as returned by
#'     [rowFitModel()]. Can also be a single linear model, in which case
#'     the data is adjusted with the same global model. See details for more
#'     information.
#'
#' @param shiftNegative `character` specifying the method to be used to
#'     avoid adjusted values to become negative. Allowed values are
#'     `"none"` (no shift, default), `"min"` (shift intensities of rows with
#'     at least one negative value by adding this value to all intensities for
#'     that row) and `"globalMin"` (shifts the complete
#'     matrix by them smallest negative intensity). See details for more
#'     information.
#'
#' @param fitOnSubset For `adjustWithModel`: `numeric` or `logical` optionally
#'     specifying a subset of columns in `y` that should be used for the model
#'     fitting. Can be e.g. the index of quality control sample columns in `y`
#'     if the model should be fit exclusively on those while adjusting all
#'     values in `y`.
#'
#' @param minValues `numeric(1)` defining the minimum number of data points
#'     required to perform the model fitting.
#'
#' @inheritParams model-fitting
#'
#' @return
#'
#' `applyModelAdjustment` returns, depending on the input `y`, a `numeric` or a
#' `matrix` with the adjusted values.
#' `adjustWithModel` returns a `matrix`, same dimension than `y` with the
#' adjusted values.
#' 
#' @noRd
#'
#' @author Johannes Rainer
#'
#' @md
#' 
#' @rdname model-based-adjustment
#'
#' @references
#' 
#' Wehrens R, Hageman JA, van Eeuwijk F, Kooke R, Flood PJ, Wijnker E,
#' Keurentjes JJ, Lommen A, van Eekelen HD, Hall RD Mumm R and de Vos RC.
#' Improved batch correction in untargeted MS-based metabolomics.
#' \emph{Metabolomics} 2016; 12:88.
#'
#' @examples
#'
#' ## Adjusting values using a model that includes injection index and
#' ## batch:
#' y <- c(2, 3, 2.7, 3.5, 3.8, 4.6, 5.9, 8, 4, 5.1, 5.6, 6.8, 7.1, 8.1, 8.9, 9.3)
#' dta <- data.frame(inj_idx = 1:length(y),
#'     batch = c(rep("a", 8), rep("b", 8)))
#' plot(dta$inj_idx, y, col = ifelse(dta$batch == "a", "red", "blue"))
#'
#' ## Adjusting the data with a model that assumes similar injection index
#' ## dependency, but different absolute abundances between the batches.
#' mdl <- fitModel(y ~ inj_idx + batch, y = y, data = dta)
#'
#' ## Adjusting the data.
#' res <- applyModelAdjustment(y, data = dta, lmod = mdl)
#' points(dta$inj_idx, res, col = ifelse(dta$batch == "a", "red", "blue"),
#'     pch = 16)
#'
#' ## Adjusting the data with a batch-only model.
#' ## We fit a model that adjusts only differences between the batches but
#' ## not any injection index dependent effects.
#' plot(dta$inj_idx, y, col = ifelse(dta$batch == "a", "red", "blue"))
#'
#' mdl <- xcms:::fitModel(y ~ batch, y = y, data = dta)
#' res <- xcms:::applyModelAdjustment(y, data = dta, lmod = mdl)
#' points(dta$inj_idx, res, col = ifelse(dta$batch == "a", "red", "blue"),
#'     pch = 16)
#'
#' ## This did not adjust the signal drift in the data, but removed differences
#' ## between the two batches:
#' ## before normalization
#' mean(y[dta$batch == "a"])
#' mean(y[dta$batch == "b"])
#'
#' ## after normalization
#' mean(res[dta$batch == "a"])
#' mean(res[dta$batch == "b"])
#' 
#' 
#' ## Adjusting a model that assumes different slopes in each batch.
#' y <- c(2, 3, 2.7, 3.5, 3.8, 4.6, 5.9, 8, 4, 4.1, 4.3, 5.1, 5.2, 5.8, 6.3, 6.1)
#' dta <- data.frame(inj_idx = 1:length(y),
#'     batch = c(rep("a", 8), rep("b", 8)))
#' plot(dta$inj_idx, y, col = ifelse(dta$batch == "a", "red", "blue"))
#'
#' ## Fit a model assuming a batch-independent drift as before
#' mdl_1 <- xcms:::fitModel(y ~ inj_idx + batch, y = y, data = dta)
#' res_1 <- xcms:::applyModelAdjustment(y, data = dta, lmod = mdl_1)
#' points(dta$inj_idx, res_1, col = ifelse(dta$batch == "a", "red", "blue"),
#'     pch = 0)
#'
#' ## And now adjust with a batch-dependent signal drift
#' mdl_2 <- xcms:::fitModel(y ~ inj_idx * batch, y = y, data = dta)
#'
#' ## The fitted models for each batch are:
#' abline(mdl_2$coefficients[1], mdl_2$coefficients[2], col = "red")
#' abline(sum(mdl_2$coefficients[c(1, 3)]), sum(mdl_2$coefficients[c(2, 4)]),
#'     col = "blue")
#'
#' ## Adjusting the values with the model
#' res_2 <- xcms:::applyModelAdjustment(y, data = dta, lmod = mdl_2)
#' points(dta$inj_idx, res_2, col = ifelse(dta$batch == "a", "red", "blue"),
#'     pch = 16)
#'
applyModelAdjustment <- function(y, data, lmod,
                            shiftNegative = c("none", "min","globalMin")) {
    shiftNegative <- match.arg(shiftNegative)
    if (is(lmod, "lm") | is(lmod, "lmrob"))
        lmod <- list(lmod)
    if (!is.matrix(y))
        y <- matrix(y, nrow = 1)
    if (length(lmod) != nrow(y))
        lmod <- rep(lmod[1], nrow(y))
    y_new <- y
    for (i in which(lengths(lmod) > 0)) {
        y_new[i, ] <- .adjust_with_linear_model(y = y[i, ], data,
                                                model = lmod[[i]])
    }
    ## Check if we have to shift values...
    if (any(y_new < 0, na.rm = TRUE)) {
	if (shiftNegative == "none") {
	    message("Note: some adjusted values are < 0.")
	}
	if (shiftNegative == "min") {
	    ## Shift selected rows by their row min
	    mins <- apply(y_new, MARGIN = 1, min, na.rm = TRUE)
	    idx <- which(mins < 1)
	    y_new[idx, ] <- y_new[idx, ] + abs(mins[idx]) + 1e-6
	    message("Shifting ", length(idx), " of the ", nrow(y_new), " rows ",
		    "to avoid negative values.")
	}
	if (shiftNegative == "globalMin") {
	    ## Shifting ALL rows by the smallest value in the matrix.
	    shiftVal <- abs(min(y_new, na.rm = TRUE))
	    message("Shifting all values by ", format(shiftVal, digits = 3))
	    y_new <- y_new + shiftVal + 1e-6
	}
    }
    if (nrow(y_new) == 1)
        y_new[1, ]
    else y_new
}

#' @noRd
#'
#' @param imputeEnds `logical(1)` whether the first and last measurement used
#'     for for the model fitting in each row should be replaced by the closest
#'     non-missing value if it is `NA`. See details for more information.
#'
#' The `imputeEnds` parameter in `adjustWithModel` allows to perform a more
#' conservative model fit by replacing `NA`s at the first and last value in
#' a row that is used for fitting (i.e. `y[x, fitOnSubset][1]` and
#' `y[x, fitOnSubset][ncol(y)]` for any row `x`) with the closest non-missing
#' value in that row. Be aware that this assumes the values to be ordered by
#' injection index and makes only sense if the model specified in `formula`
#' is defined to adjust an injection order-dependent signal drift.
#' `imputeEnds = TRUE` avoids thus exxagerated adjustments at the beginning and
#' end of a injection series, if there are no valid measurements available
#' at the ends for model fitting.
#' 
#' @md
#' 
#' @rdname model-based-adjustment
adjustWithModel <- function(y, data = NULL, model = y ~ injection_idx,
                            fitOnSubset = 1:ncol(y), minVals = 4,
                            method = "lm",
                            shiftNegative = c("none", "min","globalMin")) {
    ## Input argument checking...
    if (!is.matrix(y))
        stop("'y' is supposed to be a matrix")
    shiftNegative <- match.arg(shiftNegative)
    if (is.logical(fitOnSubset))
	fitOnSubset <- which(fitOnSubset)
    if (!all(fitOnSubset %in% 1:ncol(y)))
	stop("'fitOnSubset' should contain indices between 1 and 'ncol(y)'")
    data_fit <- data
    if (!is.null(data_fit))
	data_fit <- data_fit[fitOnSubset, , drop = FALSE]
    ## First fitting the model.
    message("Fitting the model to the features ... ", appendLF = FALSE)
    y_fit <- y[, fitOnSubset, drop = FALSE]
    lms <- rowFitModel(formula = model, data = data_fit,
                       y = y[, fitOnSubset, drop = FALSE],
                       minVals = minVals, method = method)
    message("OK")
    message("Applying models to adjust values ... ", appendLF = FALSE)
    y_new <- applyModelAdjustment(y = y, data = data, lmod = lms,
                                  shiftNegative = shiftNegative)
    message("OK")
    if (sum(lengths(lms) == 0))
	message("Did not correct ", sum(lengths(lms) == 0), " of the ",
		length(y), " rows because of too few data points to fit the ",
		"model.")
    y_new
}

#' Replace `NA`s at the beginning and end of an injection series in each batch
#' with the closest non-NA value. This aims to avoid exagerated fits estimating
#' the signal dependent drift.
#'
#' This function could be useful if some of the features adjusted with
#' [adjustWithModel()] show exagerated adjustment slopes.
#' 
#' @md
#'
#' @noRd
#' 
#' @examples
#' x <- c(NA, 3, 4, 6, 4, 2, NA, 3, NA, 4, 5, 6, NA)
#' injIndex <- c(1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6)
#' batch <- c(rep("b", 7), rep("a", 6))
#'
#' replaceNaOnEnds(x, injIndex, batch = batch)
replaceNaOnEnds <- function(x, injIndex, batch, minValues = 4) {
    if (is.null(dim(x))) {
        if (sum(!is.na(x)) < 4)
            return(x)
        if (missing(injIndex))
            injIndex <- seq_along(x)
        if (missing(batch))
            batch <- rep("a", length(x))
        lx <- length(x)
        if (length(injIndex) != lx | length(batch) != lx)
            stop("length of 'injIndex' and 'batch' should match length of 'x'")
        for (btch in unique(batch)) {
            cur_btch <- batch == btch
            idx <- order(injIndex[cur_btch])
            nona <- !is.na(x[cur_btch][idx])
            if (!nona[1]) {
                ## Replace first if NA
                x[cur_btch][idx][1] <- x[cur_btch][idx][min(which(nona))]
            }
            if (!nona[length(nona)]) {
                ## Replace last
                x[cur_btch][idx][length(idx)] <- x[cur_btch][idx][max(which(nona))]
            }
        }
        x
    } else {
        ## Process matrix.
        if (missing(injIndex))
            injIndex <- seq_len(ncol(x))
        if (missing(batch))
            batch <- rep("a", ncol(x))
        t(apply(x, MARGIN = 1, FUN = replaceNaOnEnds, 
                injIndex = injIndex, batch = batch, minValues = minValues))
    }
}
