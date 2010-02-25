testMultiFactor <- function() {
  library(faahKO)

  files <- c("./E260/const/KO/ko15.CDF", "./E260/const/KO/ko16.CDF",
             "./E260/const/KO/ko18.CDF", "./E260/const/KO/ko19.CDF",
             "./E260/const/KO/ko21.CDF", "./E260/const/KO/ko22.CDF",
             "./E261/const/WT/wt15.CDF", "./E261/const/WT/wt16.CDF",
             "./E261/const/WT/wt18.CDF", "./E261/const/WT/wt19.CDF",
             "./E261/const/WT/wt21.CDF", "./E261/const/WT/wt22.CDF")
  pd <- xcms:::phenoDataFromPaths(files)
  
  phenoData(faahko) <- pd
  xsg <- group(faahko)
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
  
  phenoData(faahko) <- pd
  xs <- fillPeaks(group(faahko))
  dr <- diffreport(xsgf, class1="E260.KO", class2="E261.WT")
}
