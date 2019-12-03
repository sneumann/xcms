#' @title Replace missing values with a proportion of the row minimum
#'
#' @description
#'
#' `imputeRowMin` imputes missing values in `x` by replacing `NA`s in each row
#' with a proportion of the minimal value for that row (i.e.
#' `min_fraction * min(x[i, ])`).
#'
#' @param x `matrix` with abundances, rows being features/metabolites and
#'     columns samples.
#'
#' @param min_fraction `numeric(1)` with the fraction of the row minimum that
#'     should be used to replace `NA` values in that row.
#'
#' @seealso `imputeLCMD` package for more left censored imputation functions.
#'
#' @family imputation functions
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @examples
#'
#' library(faahKO)
#' data("faahko")
#'
#' xset <- group(faahko)
#' mat <- groupval(xset, value = "into")
#'
#' mat_imp <- imputeRowMin(mat)
#'
#' head(mat)
#' head(mat_imp)
#'
#' ## Replace with 1/8 of the row mimimum
#' head(imputeRowMin(mat, min_fraction = 1/8))
imputeRowMin <- function(x, min_fraction = 1/2) {
    for (i in 1:nrow(x)) {
        nas <- is.na(x[i, ])
        if (all(nas))
            next
        if (any(nas)) {
            minx <- min(x[i, !nas])
            x[i, nas] <- minx * min_fraction
        }
    }
    x
}

#' @title Impute missing values with random numbers based on the row minimum
#'
#' @description
#'
#' Replace missing values with random numbers.
#' When using the `method = "mean_sd"`, random numbers will be generated
#' from a normal distribution based
#' on (a fraction of) the row min and a standard deviation estimated from the
#' linear relationship between row standard deviation and mean of the full data
#' set. Parameter `sd_fraction` allows to further reduce the estimated
#' standard deviation.
#' When using the method `method = "from_to"`, random numbers between 2 specific values
#' will be generated.
#'
#' @details
#'
#' For method **mean_sd**, imputed
#' values are taken from a normal distribution with mean being a
#' user defined fraction of the row minimum and the standard deviation
#' estimated for that mean based on the linear relationship between row
#' standard deviations and row means in the full matrix `x`.
#'
#' To largely avoid imputed values being negative or larger than the *real*
#' values, the standard deviation for the random number generation is estimated
#' ignoring the intercept of the linear model estimating the relationship
#' between standard deviation and mean. If `abs = TRUE` `NA` values are
#' replaced with the absolute value of the random values.
#'
#' For method **from_to**, imputed values are taken between 2 user defined
#' fractions of the row minimum.
#'
#' @inheritParams imputeRowMin
#'
#' @param method method `character(1)` defining the imputation method.
#' See description for details. Defaults to `method = "mean_sd"`.
#'
#' @param min_fraction `numeric(1)` with the fraction of the row minimum that
#' should be used to replace `NA` values in that row in case that `mean_sd`
#' method is specified. When  using `from_to` method, this value will be the
#' one used to calculate the maximum value for replace `NA` values in that row.
#'
#' @param min_fraction_from `numeric(1)` with the fraction of the row minimum
#' that should be used to calculate the minimum value for replace `NA` values
#' in that row. This parameter is used only in case that `from_to` method is
#' specified.
#'
#' @param sd_fraction `numeric(1)` factor to reduce the estimated standard
#'     deviation. This parameter is used only in case that `mean_sd` method is
#'     specified.
#'
#' @param abs `logical(1)` to force imputed values to be strictly positive.
#'
#' @family imputation functions
#'
#' @seealso `imputeLCMD` package for more left censored imputation functions.
#'
#' @md
#'
#' @author Johannes Rainer, Mar Garcia-Aloy
#'
#' @examples
#'
#' library(faahKO)
#' data("faahko")
#'
#' xset <- group(faahko)
#' mat <- groupval(xset, value = "into")
#'
#' ## Estimate the relationship between row sd and mean. The standard deviation
#' ## of the random distribution is estimated on this relationship.
#' mns <- rowMeans(mat, na.rm = TRUE)
#' sds <- apply(mat, MARGIN = 1, sd, na.rm = TRUE)
#' plot(mns, sds)
#' abline(lm(sds ~ mns))
#'
#' mat_imp_meansd <- imputeRowMinRand(mat, method = "mean_sd")
#' mat_imp_fromto <- imputeRowMinRand(mat, method = "from_to")
#'
#' head(mat)
#' head(mat_imp_meansd)
#' head(mat_imp_fromto)
imputeRowMinRand <- function(x, method = c("mean_sd", "from_to"),
                             min_fraction = 1/2, min_fraction_from = 1/1000,
                             sd_fraction = 1, abs = TRUE) {
    method <- match.arg(method)
    if (method == "mean_sd") {
        row_means <- rowMeans(x, na.rm = TRUE)
        row_sds <- apply(x, MARGIN = 1, sd, na.rm = TRUE)
        sd_mean_lm <- lm(row_sds ~ row_means)
        for (i in 1:nrow(x)) {
            nas <- is.na(x[i, ])
            if (all(nas))
                next
            if (any(nas)) {
                minx <- min(x[i, !nas]) * min_fraction
                rndm <- rnorm(n = sum(nas), mean = minx,
                              sd = (minx * abs(sd_mean_lm$coefficients[2])) *
                                  sd_fraction)
                if (abs)
                    x[i, nas] <- abs(rndm)
                else x[i, nas] <- rndm
            }
        }
    }
    if (method == "from_to") {
        for (i in 1:nrow(x)) {
            x[i, is.na(x[i, ])] <-
                runif(min = min(x[i, ], na.rm = TRUE) * min_fraction_from,
                      max = min(x[i, ], na.rm = TRUE) * min_fraction,
                      n = sum(is.na(x[i, ])))
        }
    }
    x
}
