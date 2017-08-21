#' @include DataClasses.R functions-utils.R


#' @title Fit linear model row-wise to a matrix or data.frame
#' 
#' @description Simple function to fit linear models row-wise to the provided
#'     data.
#'
#' @details For \code{method = "lmrob"} robust regression is performed using
#'     the \code{\link[robustbase]{lmrob}} function with settings
#'     \code{settings = "KS2014"} and \code{method = "SMDB"}.
#'     The function will perform by default parallel fitting of the models
#'     based on the global parallel processing settings.
#' 
#' @note Between batch correction in the form of \code{y ~ idx * batch} is
#'     currently problematic, because we don't yet check if there are too few
#'     values within each batch.
#'
#' @param formula \code{formula} representing the model.
#'
#' @param data \code{data.frame} containing the data to be fitted (e.g. the
#'     \code{pData} of an \code{\link{XCMSnExp}} object.
#'
#' @param y \code{matrix} or \code{data.frame} with the response variable. The
#'     model is fit to each row of this matrix (which can be e.g. the
#'     \code{\link{featureValues}} matrix).
#'
#' @param minVals \code{integer(1)} defining the minimum number of values to be
#'     used for the fitting. Model fitting is skipped for rows in \code{y} with
#'     less than \code{minVals} non-NA values.
#'
#' @param method \code{character} defining the method/function to be used for
#'     model fitting. Allowed values are \code{"lm"} for least squares
#'     regression and \code{"lmrob"} for robust regression using the
#'     \code{\link[robustbase]{lmrob}} function.
#'
#' @param BPPARAM optional parameter specifying parallel processing settings.
#'
#' @return A \code{list} with the fitted linear models or \code{NULL} for rows
#'     with too few data points.
#'
#' @noRd
#' 
#' @author Johannes Rainer
fitModel <- function(formula, data, y, minVals = 4,
                     method = c("lm", "lmrob"), BPPARAM = bpparam()) {
    method <- match.arg(method, c("lm", "lmrob"))
    if (missing(formula) || !is(formula, "formula"))
        stop("'formula' has to be submitted and should be a formula!")
    if (missing(data) || !is(data, "data.frame"))
        stop("'data' has to be a 'data.frame'!")
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
    ## Done with checking.
    force(y)
    force(data)
    force(formula)
    force(minVals)
    force(method)
    force(BPPARAM)
    ## Subset data to contain only explanatory variables
    data <- data[, vars[-1], drop = FALSE]
    if (is.null(rownames(y)))
        rownames(y) <- 1:nrow(y)
    y <- split.data.frame(y, f = rownames(y))
    ## Determine fetures we skip because of too few data points.
    do_em <- which(unlist(lapply(y, function(z) sum(!is.na(z)) >= minVals)))
    res <- vector("list", length(y))
    names(res) <- names(y)
    sttngs <- list()
    if (method == "lmrob") {
        ## Force use of the KS2014 settings in lmrob and increase the
        ## scale-finding iterations to avoid some of the warnings.
        sttngs <- lmrob.control("KS2014")
        sttngs$maxit.scale <- 10000
        sttngs$k.max <- 10000
        sttngs$refine.tol <- 1e-7
    }
    if (length(do_em)) {
        ## fit the model
        res[do_em] <- bplapply(y[do_em], FUN = function(z, formula., data.,
                                                        minVals., lmeth,
                                                        sttngs) {
            ## TODO: need to check what happens if we're also performing between
            ## batch correction and we have too few samples per batch! 
            ## ## Removing all missing values - could eventually skip that.
            ## z <- as.numeric(z)
            ## not_na <- !is.na(z)
            ## data. <- droplevels(data.frame(y = z[not_na],
            ##                                data.[not_na, , drop = FALSE]))
            data. <- data.frame(y = as.numeric(z), data.)
            if (lmeth == "lm")
                return(lm(formula., data = data., model = FALSE))
            if (lmeth == "lmrob") {
                set.seed(123)
                return(lmrob(formula., data = data., model = FALSE,
                                         setting = sttngs))
            }
            if (lmeth == "rlm")
                stop("Not yet implemented")
            ##     return(MASS::rlm(formula., data = data.))
        }, formula. = formula, data. = data, minVals. = minVals, lmeth = method,
        sttngs = sttngs, BPPARAM = BPPARAM)
    }
    res
}

## Define a simple function that does the adjustment for us.
#' @title Adjust the injection order-dependent signal drift using linear models
#'
#' @description \code{adjustDriftWithModel} first fits the specified model to
#'     each individual row in \code{y} and subsequently adjusts \code{y} based
#'     on these fitted models. This enables a signal drift and batch correction
#'     as described in [Wehrens 2016].
#'
#' @details For some rows/features values can become negative after adjustment.
#'     To avoid this, a constant can be added to the adjusted intensities of
#'     such features. Parameter \code{shiftNegative} allows to specify how this
#'     constant is to be determined. For \code{shiftNegative = "min"}, if one
#'     of the adjusted values of a row is \code{< 1}, the minimum intensity (+1)
#'     is added to each intensity. Shifting values for rows that do not only
#'     have negative values, but values \code{< 1}, ensures that adjusted values
#'     are larger 1 (which might be important if \code{y} is in log2 scale. 
#'     This shifting is done on a per-feature basis. Alternatively, the
#'     \code{globalMin} \emph{globally} shifts the complete matrix by the
#'     minimum value (if it is negative).
#'
#' @note Rows with fewer than \code{minValues} data points that can be used
#'     for the model fit are returned un-adjusted.
#'
#' @param x \code{numeric} \code{matrix} or \code{data.frame} with the values
#'     that should be corrected.
#'
#' @param data \code{data.frame} with additional variables to the model.
#'
#' @param fitOnSubset \code{numeric} or \code{logical} optionally specifying a
#'     subset of columns in \code{y} that should be used for the model fitting.
#'     Can be e.g. the index of quality control sample columns in \code{y} if
#'     the model should be fit exclusively on those while adjusting all columns
#'     of \code{y}.
#'
#' @param minValues \code{numeric(1)} defining the minimum number of data points
#'     required to perform the model fitting.
#'
#' @param method \code{character} specifying the model fitting function that
#'     should be used (\code{"lm", \code{"rlm"} or \code{"lmrob"}}.
#'
#' @param shiftNegative \code{character} specifying the method to be used to
#'     avoid adjusted values to become negative. Allowed values are
#'     \code{"none"} (no shift), \code{"min"} (shift intensities of rows with
#'     at least one negative value by adding this value +1 to all intensities)
#'     and \code{"globalMin"} (shifts the complete
#'     matrix by them smallest negative intensity). See details for more
#'     information.
#' 
#' @return A \code{list} with two elements: \code{"x"} with the adjusted input
#'     matrix \code{y} and \code{"fit"} with the fitted models. The latter can
#'     be used for quality control purposes or to e.g. identify the most
#'     adjusted rows.
#'
#' @author Johannes Rainer
#'
#' @references
#' Wehrens R, Hageman JA, van Eeuwijk F, Kooke R, Flood PJ, Wijnker E,
#' Keurentjes JJ, Lommen A, van Eekelen HD, Hall RD Mumm R and de Vos RC.
#' Improved batch correction in untargeted MS-based metabolomics.
#' \emph{Metabolomics} 2016; 12:88.
#' @noRd
adjustDriftWithModel <- function(y, data = NULL, model = y ~ injection_idx,
                                 fitOnSubset = 1:ncol(y), minVals = 4,
                                 method = "lm",
                                 shiftNegative = c("none", "min","globalMin")) {
    shiftNegative <- match.arg(shiftNegative)
    ## Input argument checking...
    if (is.logical(fitOnSubset))
        fitOnSubset <- which(fitOnSubset)
    if (!all(fitOnSubset %in% 1:ncol(y)))
        stop("'fitOnSubset' should contain indices between 1 and 'ncol(y)'")
    data_fit <- data
    if (!is.null(data_fit))
        data_fit <- data_fit[fitOnSubset, , drop = FALSE]
    ## First fitting the model.
    message("Fitting the model to the features ... ", appendLF = FALSE)
    lms <- fitModel(formula = model, data = data_fit,
                    y = y[, fitOnSubset, drop = FALSE],
                    minVals = minVals, method = method)
    message("OK")
    message("Applying models to adjust values ... ", appendLF = FALSE)
    if (is.null(rownames(y)))
        rownames(y) <- seq_len(nrow(y))
    y_cn <- colnames(y)
    y <- split.data.frame(y, f = rownames(y))
    res <- mapply(y, lms, FUN = function(z, lmod, data.) {
        z <- as.numeric(z)
        if (length(lmod) == 0)
            return(z)
        rownames(data.) <- NULL
        preds <- predict(lmod, newdata = data.frame(y = z, data.))
        z + mean(z, na.rm = TRUE) - preds
    }, MoreArgs = list(data. = data), SIMPLIFY = FALSE)
    res <- do.call(rbind, res)
    message("OK")
    if (sum(lengths(lms) == 0))
        message("Did not correct ", sum(lengths(lms) == 0), " of the ",
                length(y), " rows because of too few data points to fit the ",
                "model.")
    rm(y)
    ## Check if we have to shift values...
    if (any(res < 1, na.rm = TRUE)) {
        if (shiftNegative == "none") {
            message("Note: some adjusted values are < 1.")
        }
        if (shiftNegative == "min") {
            ## Shift selected rows by their row min + 1
            ## Include here also < 1 so that values potentially in log scale
            ## between 0 and 1 are adjusted as well.
            mins <- apply(res, MARGIN = 1, function(z) min(z, na.rm = TRUE))
            idx <- which(mins < 1)
            res[idx, ] <- res[idx, ] + abs(mins[idx]) + 1
            message("Shifting ", length(idx), " of the ", nrow(res), " rows ",
                    "to avoid negative values.")
        }
        ## if (shiftNegative == "log") {
        ##     ## Shift selected rows by the difference.
        ##     mins <- apply(res, MARGIN = 1, function(z) min(z, na.rm = TRUE))
        ##     idx <- which(mins < 1)
        ##     ## res[idx, ] <- res[idx, ] + 1
        ##     res[idx, ] <- res[idx, ] + (1 - mins)
        ##     message("Shifting ", length(idx), " of the ", nrow(res),
        ##             " rows to avoid values between 0 and 1.")
        ## }
        if (shiftNegative == "globalMin") {
            ## Shifting ALL rows by the smallest value in the matrix.
            shiftVal <- abs(min(res, na.rm = TRUE)) + 1
            message("Shifting all values by ", format(shiftVal, digits = 3))
            res <- res + shiftVal
        }
    }
    colnames(res) <- y_cn
    return(list(y = res, fit = lms))
}
