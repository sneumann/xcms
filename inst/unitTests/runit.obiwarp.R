test.obiwarp.default <- function() {
    faahko_sub <- faahko[,1:2]
    xr <- retcor(faahko_sub, method="obiwarp", profStep = 10)
}

test.obiwarp.local <- function() {
    faahko_sub <- faahko[,1:2]
    xr <- retcor(faahko_sub, method="obiwarp", localAlignment=1, profStep = 10)
}

test.obiwarp.cor <- function() {
    faahko_sub <- faahko[,1:2]
    xr <- retcor(faahko_sub, method="obiwarp", distFunc="cor", profStep = 10)
}

test.obiwarp.cor_opt <- function() {
    faahko_sub <- faahko[,1:2]
    xr <- retcor(faahko_sub, method="obiwarp", distFunc="cor_opt", profStep = 10)
}

test.obiwarp.cov <- function() {
    faahko_sub <- faahko[,1:2]
    xr <- retcor(faahko_sub, method="obiwarp", distFunc="cov", profStep = 10)
}

test.obiwarp.euc <- function() {
    faahko_sub <- faahko[,1:2]
    xr <- retcor(faahko_sub, method="obiwarp", distFunc="euc", profStep = 10)
}

test.obiwarp.prd <- function() {
    faahko_sub <- faahko[,1:2]
    xr <- retcor(faahko_sub, method="obiwarp", distFunc="prd", profStep = 10)
}
