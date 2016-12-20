test.obiwarp.default <- function() {
    xr <- retcor(faahko, method="obiwarp", profStep = 10)
}

test.obiwarp.local <- function() {
    xr <- retcor(faahko, method="obiwarp", localAlignment=1, profStep = 10)
}

test.obiwarp.cor <- function() {
    xr <- retcor(faahko, method="obiwarp", distFunc="cor", profStep = 10)
}

test.obiwarp.cor_opt <- function() {
    xr <- retcor(faahko, method="obiwarp", distFunc="cor_opt", profStep = 10)
}

test.obiwarp.cov <- function() {
    xr <- retcor(faahko, method="obiwarp", distFunc="cov", profStep = 10)
}

test.obiwarp.euc <- function() {
    xr <- retcor(faahko, method="obiwarp", distFunc="euc", profStep = 10)
}

test.obiwarp.prd <- function() {
    xr <- retcor(faahko, method="obiwarp", distFunc="prd", profStep = 10)
}
