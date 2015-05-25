testSplit <- function() {
    xsl <- split(faahko,sampclass(faahko))
    checkEqualsNumeric(length(xsl), 2)
}

testSplitAll <- function() {
    xsl <- split(faahko,sampnames(faahko))
    checkEqualsNumeric(length(xsl), length(sampnames(faahko)))
}

testSplitFirst <- function() {
    xsl <- split(faahko,c(1,2,2,2,2,2,2,2,2,2,2,2))
    checkEqualsNumeric(length(xsl), 2)
}
testSplitLast <- function() {
    xsl <- split(faahko,c(2,2,2,2,2,2,2,2,2,2,2,1))
    checkEqualsNumeric(length(xsl), 2)
}
testSplitMiddle <- function() {
    xsl <- split(faahko,c(2,2,2,2,2,1,2,2,2,2,2,2))
    checkEqualsNumeric(length(xsl), 2)
}
testSplitNone <- function() {
    xsl <- split(faahko,c(2,2,2,2,2,2,2,2,2,2,2,2))
    checkEqualsNumeric(length(xsl), 1)
}


testCombine <- function() {
    xsl <- split(faahko,sampclass(faahko))
    checkEqualsNumeric(length(sampnames(c(xsl[[1]], xsl[[2]]))), 12)
}


testSplitFactorShort = function() {
  f = c(1,2,2)
  xsl <- split(faahko, f)
  for (x in unique(f)) {
    num.samps = sum(x==f)
    checkEqualsNumeric(length(xsl[[x]]@filepaths), num.samps)
    checkEqualsNumeric(nrow(xsl[[x]]@phenoData), num.samps)
    checkEqualsNumeric(length(xsl[[x]]@rt$raw), num.samps)
  }
}

testSplitFactorLongDrop = function() {
  f = c(1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,3)
  actual.f = f[1:nrow(faahko@phenoData)]
  xsl <- split(faahko, f)
  for (x in unique(actual.f)) {
    num.samps = sum(x==actual.f)
    checkEqualsNumeric(length(xsl[[x]]@filepaths), num.samps)
    checkEqualsNumeric(nrow(xsl[[x]]@phenoData), num.samps)
    checkEqualsNumeric(length(xsl[[x]]@rt$raw), num.samps)
  }
}

testSplitPhenoData <- function(){
    xset <- faahko
    ## make a dummy phenoDate data.frame
    pd <- data.frame(class=xset$class, exp=rep(1:6, 2))
    rownames(pd) <- sampnames(xset)
    phenoData(xset) <- pd
    ## split by sampleclass
    xsetList <- split(xset, sampclass(xset))
    checkEquals(xsetList[[1]]$class, droplevels(pd[pd$class=="KO", "class"]))
    checkEquals(xsetList[[2]]$exp, pd[pd$class=="WT", "exp"])
    ## now split by the experiment
    xsetList <- split(xset, xset$exp)
    ## just check some stuff...
    checkEquals(length(xsetList[[1]]$class), 2)
    checkEquals(nrow(phenoData(xsetList[[4]])), 2)
    ## have to make sure we also do have the mslevel, profinfo,
    ## scanrange, dataCorrection and polarity set.
    checkEquals(mslevel(xset), mslevel(xsetList[[1]]))
    checkEquals(scanrange(xset), scanrange(xsetList[[1]]))
    checkEquals(xset@polarity, xsetList[[1]]@polarity)
}


testSubset <- function(){
    ## first testing just the plain xset
    ##data(faahko, package="faahKO")
    xset <- faahko
    idx <- 8
    xsub <- xset[, idx]
    .compare(xset, xsub, idx)
    ## repeating the same using different ordering.
    idx <- c(8, 1, 5)
    xsub <- xset[, idx]
    .compare(xset, xsub, idx)

    ## performing a grouping, retention time correction and second grouping
    ##data(faahko, package="faahKO")
    xset <- group(faahko, method="density")
    xset <- retcor(xset, method="loess", family="symmetric")
    xset <- group(xset)
    ##xset <- fillPeaks(xset)
    xsub <- xset[, idx]
    .compare(xset, xsub, idx)

    ## subset by names
    idx <- c("ko15", "wt19", "ko21", "wt22")
    xsub <- xset[, idx]
    intidx <- match(idx, sampnames(xset))
    .compare(xset, xsub, intidx)

    ## and at last using logical.
    idx <- rep(FALSE, 12)
    idx[c(2, 8)] <- TRUE
    xsub <- xset[, idx]
    intidx <- which(idx)
    .compare(xset, xsub, intidx)
    ## testing some errors...
    checkException(xset[1, ])
    checkException(xset[, 20])
    checkException(xset[, "not there"])
}


    ## internal testing functions
.compare <- function(xs, xsub, idx){
    ## sampclass?
    checkEquals(as.character(sampclass(xsub)),
                as.character(sampclass(xs)[idx]))
    ## sampnames?
    checkEquals(sampnames(xsub),
                sampnames(xs)[idx])
    ## filepaths?
    checkEquals(filepaths(xsub), filepaths(xs)[idx])
    ## peaks
    peaks.all <- peaks(xs)
    peaks.sub <- peaks(xsub)
    cat("Comparing peaks...")
    for(i in 1:length(idx)){
        ## peaks
        checkEquals(peaks.sub[peaks.sub[, "sample"]==i, "into"],
                    peaks.all[peaks.all[, "sample"]==idx[i], "into"])
        ## rt, raw
        checkEquals(xs@rt$raw[[idx[i]]], xsub@rt$raw[[i]])
        ## rt, corrected
        checkEquals(xs@rt$corrected[[idx[i]]], xsub@rt$corrected[[i]])
    }
    cat("OK\n")
    if(length(xs@groupidx) > 0){
        cat("Comparing peak groups...")
        ## check the groups and groupidx using the groupval method.
        into.all <- groupval(xs, value="into", method="maxint")
        into.sub <- groupval(xsub, value="into", method="maxint")
        for(i in 1:length(idx)){
            checkEquals(into.all[, idx[i]], into.sub[, i])
        }
        cat("OK\n")
    }
}
