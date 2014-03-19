test.obiwarp.default <- function() {
    xr <- retcor(faahko, method="obiwarp")
}

test.obiwarp.local <- function() {
    xr <- retcor(faahko, method="obiwarp", localAlignment=1)
}

test.obiwarp.cor <- function() {
    xr <- retcor(faahko, method="obiwarp", distFunc="cor")
}

test.obiwarp.cor_opt <- function() {
    xr <- retcor(faahko, method="obiwarp", distFunc="cor_opt")
}

test.obiwarp.cov <- function() {
    xr <- retcor(faahko, method="obiwarp", distFunc="cov")
}

test.obiwarp.euc <- function() {
    xr <- retcor(faahko, method="obiwarp", distFunc="euc")
}

test.obiwarp.prd <- function() {
    xr <- retcor(faahko, method="obiwarp", distFunc="prd")
}
