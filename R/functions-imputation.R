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
#' Replace missing values with random numbers from a normal distribution based
#' on (a fraction of) the row min and a standard deviation estimated from the
#' linear relationship between row standard deviation and mean of the full data
#' set. Parameter `sd_fraction` allows to further reduce the estimated
#' standard deviation.
#'
#' @details
#'
#' Imputed values are taken from a normal distribution with mean being a
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
#' @inheritParams imputeRowMin
#'
#' @param sd_fraction `numeric(1)` factor to reduce the estimated standard
#'     deviation.
#'
#' @param abs `logical(1)` to force imputed values to be strictly positive.
#'
#' @family imputation functions
#'
#' @seealso `imputeLCMD` package for more left censored imputation functions.
#'
#' @md
#'
#' @author Johannes Rainer
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
#' mat_imp <- imputeRowMinRand(mat)
#'
#' head(mat)
#' head(mat_imp)
imputeRowMinRand <- function(x, min_fraction = 1/2,
                             sd_fraction = 1, abs = TRUE) {
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
    x
}
