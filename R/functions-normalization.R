#' @include DataClasses.R functions-utils.R


#' @title Fit linear model row-wise to a matrix or data.frame
#' 
#' @description Simple function to fit linear models row-wise to the provided
#'     data.
#'
#' @details For \code{method = "lmrob"} robust regression is performed using
#'     the \code{\link[robustbase]{lmrob}} function with settings
#'     \code{settings = "KS2014"} and \code{method = "SMDB"}. 
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
#' @return A \code{list} with the fitted linear models or \code{NULL} for rows
#'     with too few data points.
#'
#' @noRd
#'
#' @author Johannes Rainer
fitModel <- function(formula, data, y, minVals = 4,
                     method = c("lm", "lmrob", "rlm"), BPPARAM = bpparam()) {
    method <- match.arg(method, c("lm", "lmrob", "rlm"))
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
        ## sttngs <- robustbase::lmrob.control("KS2014")
        ## sttngs$maxit.scale <- 10000
        ## sttngs$k.max <- 10000
        ## sttngs$refine.tol <- 1e-7
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
                stop("Not yet implemented")
                ## set.seed(123)
                ## return(robustbase::lmrob(formula., data = data., model = FALSE,
                ##                          setting = sttngs))
            }
            if (lmeth == "rlm")
                stop("Not yet implemented")
                ## return(MASS::rlm(formula., data = data.))
        }, formula. = formula, data. = data, minVals. = minVals, lmeth = method,
        sttngs = sttngs, BPPARAM = BPPARAM)
    }
    res
}

