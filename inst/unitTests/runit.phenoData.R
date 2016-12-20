library(faahKO)
data(faahko)
files <- system.file(c("cdf/KO/ko15.CDF", "cdf/KO/ko16.CDF",
                       "cdf/KO/ko18.CDF", "cdf/KO/ko19.CDF",
                       "cdf/KO/ko21.CDF", "cdf/KO/ko22.CDF",
                       "cdf/WT/wt15.CDF", "cdf/WT/wt16.CDF",
                       "cdf/WT/wt18.CDF", "cdf/WT/wt19.CDF",
                       "cdf/WT/wt21.CDF", "cdf/WT/wt22.CDF"),
                     package = "faahKO")

testMultiFactor <- function() {
    pd <- xcms:::phenoDataFromPaths(files)
    xs <- faahko

    ##xcms::phenoData(xs) <- pd
    ## https://stat.ethz.ch/pipermail/r-devel/2008-April/049184.html
    xs <- xcms::`phenoData<-`(xs, pd)

    xsg <- group(xs)
}

testMultiFactorDiffreport <- function() {
    pd <- xcms:::phenoDataFromPaths(files)
    xs <- faahko

    ##xcms::phenoData(xs) <- pd
    ## https://stat.ethz.ch/pipermail/r-devel/2008-April/049184.html
    xs <- xcms::`phenoData<-`(xs, pd)
    xs <- group(xs)
    ## Setting the filepaths again; otherwise we will have problem finding these
    ## files ... obviously.
    filepaths(xs) <- files
    xs <- fillPeaks(xs)
    dr <- diffreport(xs, class1="KO", class2="WT")
    checkTrue(nrow(dr) > 0)
}
