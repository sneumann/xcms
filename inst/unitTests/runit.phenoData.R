testMultiFactor <- function() {
    library(faahKO)

    files <- c("./E260/const/KO/ko15.CDF", "./E260/const/KO/ko16.CDF",
               "./E260/const/KO/ko18.CDF", "./E260/const/KO/ko19.CDF",
               "./E260/const/KO/ko21.CDF", "./E260/const/KO/ko22.CDF",
               "./E261/const/WT/wt15.CDF", "./E261/const/WT/wt16.CDF",
               "./E261/const/WT/wt18.CDF", "./E261/const/WT/wt19.CDF",
               "./E261/const/WT/wt21.CDF", "./E261/const/WT/wt22.CDF")
    pd <- xcms:::phenoDataFromPaths(files)
    xs <- faahko

    ##xcms::phenoData(xs) <- pd
    ## https://stat.ethz.ch/pipermail/r-devel/2008-April/049184.html
    xs <- xcms::`phenoData<-`(xs, pd)

    xsg <- group(xs)
}

testMultiFactorDiffreport <- function() {
    library(faahKO)

    files <- c("./E260/const/KO/ko15.CDF", "./E260/const/KO/ko16.CDF",
               "./E260/const/KO/ko18.CDF", "./E260/const/KO/ko19.CDF",
               "./E260/const/KO/ko21.CDF", "./E260/const/KO/ko22.CDF",
               "./E261/const/WT/wt15.CDF", "./E261/const/WT/wt16.CDF",
               "./E261/const/WT/wt18.CDF", "./E261/const/WT/wt19.CDF",
               "./E261/const/WT/wt21.CDF", "./E261/const/WT/wt22.CDF")
    pd <- xcms:::phenoDataFromPaths(files)
    xs <- faahko

    ##xcms::phenoData(xs) <- pd
    ## https://stat.ethz.ch/pipermail/r-devel/2008-April/049184.html
    xs <- xcms::`phenoData<-`(xs, pd)

    xs <- fillPeaks(group(xs))
    dr <- diffreport(xs, class1="E260.KO", class2="E261.WT")
}
