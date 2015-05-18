## testing whether sampclass is working as we expect. internally, xcms does not use the
## alphanumeric ordering of the factor levels, so, we have to make sure that we're
## returning the sample classes as expected, even if they might be extracted from the
## phenodata data.frame that has alphanumeric ordered factors!
test.sampclass <- function(){
    library(faahKO)
    xset <- faahko
    ## grouping the peaks
    xset <- group(xset, method="density")
    ## reversing the order of the classes.
    xset.revorder <- xset
    sampclass(xset.revorder) <- c(rep("WT", 6), rep("KO", 6))
    xset.revorder <- group(xset.revorder, method="density")
    ## check if we get what we want:
    checkEquals(groups(xset)[, "KO"], groups(xset.revorder)[, "WT"])

    ## repeat that but submitting already a factor
    xset.revorder.f <- xset
    sampclass(xset.revorder.f) <- factor(c(rep("WT", 6), rep("KO", 6)))
    xset.revorder.f <- group(xset.revorder.f, method="density")
    checkEquals(groups(xset)[, "KO"], groups(xset.revorder.f)[, "WT"])

    ## next: pheno data contains a column class with a factor.
    pd <- data.frame(class=factor(c(rep("WT", 6), rep("KO", 6))))
    xset.pheno <- xset
    phenoData(xset.pheno) <- pd
    xset.pheno <- group(xset.pheno, method="density")
    checkEquals(groups(xset)[, "KO"], groups(xset.pheno)[, "WT"])

    ## next checking what happens if we submit a multi-column data.frame
    ## to sampclass<-
    pd <- data.frame(dummy=rep(c("a", "b"), 6), group=c(rep("KO", 6), rep("WT", 6)),
                     pair=c(1:6, 1:6))
    phenoData(xset) <- pd
    ## No class column, so we're returning the interaction.
    checkEquals(sampclass(xset), interaction(pd, drop=TRUE))
    ## now we're going to submit 2 columns of the pd
    sampclass(xset) <- pd[, c("group", "pair")]
    checkEquals(as.character(sampclass(xset)),
                as.character(interaction(pd[, c("group", "pair"), drop=TRUE])))
}
