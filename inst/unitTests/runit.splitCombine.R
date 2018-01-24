library(faahKO)
fs <- c(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
        system.file('cdf/KO/ko16.CDF', package = "faahKO"),
        system.file('cdf/KO/ko18.CDF', package = "faahKO"),
        system.file('cdf/KO/ko19.CDF', package = "faahKO"))

xs_1 <- xcmsSet(fs, profparam = list(step = 0), method = "centWave",
                noise = 10000, snthresh = 50)

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

############################################################
## Test splitting and .@processHistory slot.
test_split_ProcessHistory <- function() {
    spl <- split(faahko, sampclass(faahko))
    checkTrue(length(spl[[1]]@.processHistory) == 0)

    spl <- split(xs_1, f = c(1, 2, 1, 2))
    ph <- xs_1@.processHistory[c(1, 3)]
    ph[[1]]@fileIndex <- 1L
    ph[[2]]@fileIndex <- 2L
    checkEquals(spl[[1]]@.processHistory, ph)

    ph <- xs_1@.processHistory[c(2, 4)]
    ph[[1]]@fileIndex <- 1L
    ph[[2]]@fileIndex <- 2L
    checkEquals(spl[[2]]@.processHistory, ph)

    ## Add fake ProcessHistory steps.
    ph <- xs_1@.processHistory
    ph <- c(list(xcms:::ProcessHistory(fileIndex = 1:4)), ph,
            list(xcms:::ProcessHistory()))
    xs_2 <- xs_1
    xs_2@.processHistory <- ph
    ##
    spl <- split(xs_2, f = c(2, 1, 1, 1))
    ph <- xs_2@.processHistory[c(1, 3, 4, 5)]
    ph[[1]]@fileIndex <- 1L:3L
    ph[[2]]@fileIndex <- 1L
    ph[[3]]@fileIndex <- 2L
    ph[[4]]@fileIndex <- 3L
    checkEquals(spl[[1]]@.processHistory, ph)

    ph <- xs_2@.processHistory[c(1, 2)]
    ph[[1]]@fileIndex <- 1L
    ph[[2]]@fileIndex <- 1L
    checkEquals(spl[[2]]@.processHistory, ph)
}

test_c_ProcessHistory <- function() {
    spl <- split(faahko, sampclass(faahko))
    conc <- do.call("c", spl)
    checkEquals(length(conc@.processHistory), 0)

    spl <- split(xs_1, c(1, 1, 2, 2))
    conc <- c(spl[[1]], spl[[2]])
    checkEquals(xs_1@.processHistory, conc@.processHistory)

    ## Different ordering
    xs_1.1 <- xs_1[, c(1, 3)]
    xs_1.2 <- xs_1[, 2]
    xs_1.3 <- xs_1[, 4]
    ## Add a fake processing for the second one.
    ph <- xs_1.2@.processHistory
    ph <- c(list(xcms:::ProcessHistory(fileIndex = 1)), ph)
    xs_1.2@.processHistory <- ph
    checkTrue(xcms:::.validProcessHistory(xs_1.2))
    ## Combine them.
    conc <- c(xs_1.1, xs_1.2, xs_1.3)
    ## 1st
    checkEquals(conc@.processHistory[[1]], xs_1@.processHistory[[1]])
    ## 2nd
    ph <- xs_1@.processHistory[[3]]
    ph@fileIndex <- 2L
    checkEquals(conc@.processHistory[[2]], ph)
    ## 3rd
    ph <- xs_1.2@.processHistory[[1]]
    ph@fileIndex <- 3L
    checkEquals(conc@.processHistory[[3]], ph)
    ## 4th
    ph <- xs_1@.processHistory[[2]]
    ph@fileIndex <- 3L
    checkEquals(conc@.processHistory[[4]], ph)
    ## 5th
    ph <- xs_1@.processHistory[[4]]
    ph@fileIndex <- 4L
    checkEquals(conc@.processHistory[[5]], ph)
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

## Issue #133
test_c_empty <- function() {
    library(msdata)
    suppressWarnings(
        xs <- xcmsSet(system.file("microtofq/MM8.mzML", package="msdata"),
                      method="centWave", ppm=25, peakwidth=c(20, 50))
    )
    xs2 <- xcmsSet(system.file("microtofq/MM14.mzML", package="msdata"),
                   method="centWave", ppm=25, peakwidth=c(20, 50))
    comb <- c(xs, xs2)
    checkTrue(nrow(peaks(comb)) == 0)
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

    ## Testing subsetting with .processHistory:
    xsub <- xs_1[, c(2, 3)]
    ph <- xs_1@.processHistory[c(2, 3)]
    ph[[1]]@fileIndex <- 1L
    ph[[2]]@fileIndex <- 2L
    checkEquals(xsub@.processHistory, ph)
    ## Reverse ordering:
    xsub <- xs_1[, c(3, 2)]
    ph <- xs_1@.processHistory[[3]]
    ph@fileIndex <- 1L
    checkEquals(xcms:::.getProcessHistory(xsub, fileIndex = 1), list(ph))
    ph <- xs_1@.processHistory[[2]]
    ph@fileIndex <- 2L
    checkEquals(xcms:::.getProcessHistory(xsub, fileIndex = 2), list(ph))

    ## Add fake ProcessHistory before and after the real ones.
    ph <- xs_1@.processHistory
    ph <- c(list(xcms:::ProcessHistory(fileIndex = 1:4)), ph,
            list(xcms:::ProcessHistory(fileIndex = 1:4)))
    xs_1@.processHistory <- ph
    xsub <- xs_1[, c(2, 3)]
    ph <- xs_1@.processHistory[c(1, 3, 4, 6)]
    ph[[1]]@fileIndex <- 1L:2L
    ph[[2]]@fileIndex <- 1L
    ph[[3]]@fileIndex <- 2L
    ph[[4]]@fileIndex <- 1L:2L
    checkEquals(xsub@.processHistory, ph)
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

############################################################
## Test the internal .getProcessHistory function to retrieve
## specific processing history steps.
test_getProcessHistory <- function() {

    checkEquals(xcms:::.getProcessHistory(xs_1, fileIndex = 2:3),
                xs_1@.processHistory[2:3])
}
