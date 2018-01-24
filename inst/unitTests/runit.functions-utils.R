test_plotMsData <- function() {
    msd <- extractMsData(faahko_od, mz = c(334.9, 335.1), rt = c(2700, 2900))
    plotMsData(msd[[1]])
}

## library(xcms)
## library(RUnit)
## Test the .grow_trues
## test_grow_trues <- function() {
##     ## Compare performance with MSnbase:::utils.clean
##     Test <- c(1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0,
##               1, 0)
##     Expect <- c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE,
##                 FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
##                 TRUE, TRUE, TRUE, TRUE)
##     res_2 <- xcms:::.grow_trues(Test > 0)
##     checkEquals(res_2, Expect)
    
##     Test <- c(0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0)
##     Expect <- c(FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE,
##                 TRUE, FALSE)
##     res_2 <- xcms:::.grow_trues(Test > 0)
##     checkEquals(res_2, Expect)

##     Test <- c(0, 1, NA, 0, 0, 1)
##     Expect <- c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE)
##     res_2 <- xcms:::.grow_trues(Test > 0)
##     checkEquals(res_2, Expect)

##     Test <- c(0, NA, 1, 0, 0, 1, 0, 0)
##     Expect <- c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE)
##     res_2 <- xcms:::.grow_trues(Test > 0)
##     checkEquals(res_2, Expect)

##     Test <- c(0, 1, 0, 0, NA, 0, 1)
##     Expect <- c(TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE)
##     res_2 <- xcms:::.grow_trues(Test > 0)
##     checkEquals(res_2, Expect)

##     Test <- c(NA, 1, NA, NA, NA, NA, 1)
##     Expect <- c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE)
##     res_2 <- xcms:::.grow_trues(Test > 0)
##     checkEquals(res_2, Expect)    
## }

## benchmark_grow_trues <- function() {
##     set.seed(123)
##     Test <- rnorm(n = 30000)
##     Test[Test < 0] <- 0
##     Test2 <- Test > 0
##     res_1 <- MSnbase:::utils.clean(Test)
##     res_2 <- .clean(Test2)

##     Test <- c(0, 0, 1, 1, 0, 0, 0, 1, 0, 0)
##     Expect <- c(FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE)
##     res_1 <- MSnbase:::utils.clean(Test)
    
## }
