## Unit tests related to the Chromatogram class.
library(xcms)
library(RUnit)

test_Chromatogram_class <- function() {
    ch <- new("Chromatogram")
    ch@mz <- 3
    checkException(validObject(ch))
    ch@mz <- c(1, 3)
    ch@precursorMz <- 4
    checkException(validObject(ch))
    ch@precursorMz <- c(4, 4)
    ch@productMz <- 5
    checkException(validObject(ch))
    ##
    int <- rnorm(100, mean = 200, sd = 2)
    rt <- rnorm(100, mean = 300, sd = 3)
    ## check exceptions:
    checkException(xcms:::Chromatogram(intensity = int))
    chr <- Chromatogram()
    chr@rtime <- rt
    chr@intensity <- int
    checkException(validObject(chr))
    ## issue #145: values are ordered based on rtime
    chr <- Chromatogram(intensity = int, rtime = rt)
    checkEquals(rtime(chr), sort(rt))
    checkEquals(intensity(chr), int[order(rt)])
    rt <- sort(rt)
    ch <- xcms:::Chromatogram(intensity = int, rtime = rt)
    checkEquals(rtime(ch), rt)
    checkEquals(intensity(ch), int)
    checkException(xcms:::Chromatogram(aggregationFun = "other"))
    ch@aggregationFun <- "max"
    checkTrue(validObject(ch))
    checkEquals(aggregationFun(ch), "max")
    ch@aggregationFun <- "sum"
    checkTrue(validObject(ch))
    checkEquals(aggregationFun(ch), "sum")
    ch@aggregationFun <- "mean"
    checkTrue(validObject(ch))
    checkEquals(aggregationFun(ch), "mean")
    ch@aggregationFun <- "min"
    checkTrue(validObject(ch))
    checkEquals(aggregationFun(ch), "min")
    ch@fromFile <- 3L
    checkTrue(validObject(ch))
    checkEquals(fromFile(ch), 3L)
    checkEquals(length(ch), length(rt))
    ## as.data.frame
    df <- as.data.frame(ch)
    checkEquals(df, data.frame(rtime = rt, intensity = int))
    ch <- xcms:::Chromatogram(mz = c(1, 3))
    checkEquals(ch@mz, c(1, 3))
    checkEquals(mz(ch), c(1, 3))
    checkEquals(mz(ch, filter = TRUE), c(0, 0))
    ch <- xcms:::Chromatogram(filterMz = c(1, 3))
    checkEquals(ch@filterMz, c(1, 3))
    checkEquals(mz(ch, filter = TRUE), c(1, 3))
    checkEquals(mz(ch, filter = FALSE), c(0, 0))
    ch <- xcms:::Chromatogram(precursorMz = 123)
    checkEquals(ch@precursorMz, c(123, 123))
    checkEquals(precursorMz(ch), c(123, 123))
    ch <- xcms:::Chromatogram(productMz = 123)
    checkEquals(ch@productMz, c(123, 123))
    checkEquals(productMz(ch), c(123, 123))
}

test_extractChromatograms <- function() {
    ## OnDiskMSnExp
    ## TIC
    chrs <- extractChromatograms(filterFile(od_x, file = 2))
    spctr <- spectra(filterFile(od_x, file = 2))
    ints <- unlist(lapply(spctr, function(z)
        return(sum(intensity(z)))))
    checkEquals(intensity(chrs[[1]]), ints)
    checkEquals(rtime(chrs[[1]]), unlist(lapply(spctr, rtime)))
    ## BPC
    chrs <- extractChromatograms(filterFile(od_x, file = 2),
                                 aggregationFun = "max")
    ints <- unlist(lapply(spctr, function(z)
        return(max(intensity(z)))))
    checkEquals(intensity(chrs[[1]]), ints)
    checkEquals(rtime(chrs[[1]]), unlist(lapply(spctr, rtime)))
    ## XCMSnExp
    xod_x <- faahko_xod
    chrs <- extractChromatograms(filterFile(xod_x, file = 2))
    ints <- unlist(lapply(spctr, function(z)
        return(sum(intensity(z)))))
    checkEquals(intensity(chrs[[1]]), ints)
    checkEquals(rtime(chrs[[1]]), unlist(lapply(spctr, rtime)))
    ## BPC
    chrs <- extractChromatograms(filterFile(xod_x, file = 2),
                                 aggregationFun = "max")
    ints <- unlist(lapply(spctr, function(z)
        return(max(intensity(z)))))
    checkEquals(intensity(chrs[[1]]), ints)
    checkEquals(rtime(chrs[[1]]), unlist(lapply(spctr, rtime)))
    ## with adjusted retention times.
    chrs <- extractChromatograms(filterFile(xod_xgr, file = 2),
                                 adjustedRtime = FALSE, aggregationFun = "max")
    checkEquals(intensity(chrs[[1]]), ints)
    checkEquals(rtime(chrs[[1]]), unlist(lapply(spctr, rtime)))
    chrs <- extractChromatograms(filterFile(xod_xgr, file = 2,
                                            keepAdjustedRtime = TRUE),
                                 aggregationFun = "max")
    checkEquals(intensity(chrs[[1]]), ints)
    checkEquals(rtime(chrs[[1]]), rtime(xod_xgr, bySample = TRUE)[[2]])

    ## Now subsetting for mz:
    tmp <- filterFile(od_x, file = 2)
    chrs <- extractChromatograms(tmp, mz = c(300, 400))
    checkEquals(mz(chrs[[1]], filter = TRUE), c(300, 400))
    suppressWarnings(spctr <- spectra(filterMz(tmp, mz = c(300, 400))))
    ints <- unlist(lapply(spctr, function(z)
        return(sum(intensity(z)))))
    ints2 <- intensity(chrs[[1]])
    ints2[is.na(ints2)] <- 0
    checkEquals(ints2, ints)
    checkEquals(rtime(chrs[[1]]), unlist(lapply(spctr, rtime)))
    ## with adjusted retention times
    chrs <- extractChromatograms(filterFile(xod_xgr, file = 2,
                                            keepAdjustedRtime = TRUE),
                                 mz = c(300, 400))
    ints <- unlist(lapply(spctr, function(z)
        return(sum(intensity(z)))))
    ints2 <- intensity(chrs[[1]])
    ints2[is.na(ints2)] <- 0
    checkEquals(ints2, ints)
    checkEquals(rtime(chrs[[1]]), rtime(xod_xgr, bySample = TRUE)[[2]])
    
    ## Now subsetting for rt:
    chrs <- extractChromatograms(od_x, rt = c(2700, 2900))
    checkTrue(all(rtime(chrs[[1]]) >= 2700 & rtime(chrs[[1]]) <= 2900))
    checkTrue(all(rtime(chrs[[2]]) >= 2700 & rtime(chrs[[2]]) <= 2900))
    checkTrue(all(rtime(chrs[[3]]) >= 2700 & rtime(chrs[[3]]) <= 2900))
    spctr <- spectra(filterRt(od_x, rt = c(2700, 2900)))
    ints <- split(unlist(lapply(spctr, function(z) sum(intensity(z)))),
                  f = unlist(lapply(spctr, fromFile)))
    checkEquals(ints[[1]], intensity(chrs[[1]]))
    checkEquals(ints[[2]], intensity(chrs[[2]]))
    checkEquals(ints[[3]], intensity(chrs[[3]]))
    ## Using adjusted rt:
    chrs2 <- extractChromatograms(xod_xgr, rt = c(2700, 2900))
    checkTrue(all(rtime(chrs2[[1]]) >= 2700 & rtime(chrs2[[1]]) <= 2900))
    checkTrue(all(rtime(chrs2[[2]]) >= 2700 & rtime(chrs2[[2]]) <= 2900))
    checkTrue(all(rtime(chrs2[[3]]) >= 2700 & rtime(chrs2[[3]]) <= 2900))
    checkTrue(length(chrs[[1]]) != length(chrs2[[1]]))
    checkTrue(length(chrs[[2]]) == length(chrs2[[2]]))
    checkTrue(length(chrs[[3]]) != length(chrs2[[3]]))
    tmp <- filterRt(xod_xgr, rt = c(2700, 2900))
    checkEquals(rtime(chrs2[[1]]), rtime(tmp, bySample = TRUE)[[1]])
    checkEquals(rtime(chrs2[[2]]), rtime(tmp, bySample = TRUE)[[2]])
    checkEquals(rtime(chrs2[[3]]), rtime(tmp, bySample = TRUE)[[3]])
    ## Check the values...
    keepSp <- which(adjustedRtime(xod_xgr) >= 2700 &
                    adjustedRtime(xod_xgr) <= 2900)
    tmp <- xod_xgr[keepSp]
    ints <- unlist(lapply(spectra(tmp), function(z) sum(intensity(z))))
    intsL <- split(ints, fromFile(tmp))
    checkEquals(intensity(chrs2[[1]]), intsL[[1]])
    checkEquals(intensity(chrs2[[2]]), intsL[[2]])
    checkEquals(intensity(chrs2[[3]]), intsL[[3]])
    
    ## Now subsetting for rt and mz:
    chrs <- extractChromatograms(od_x, rt = c(2700, 2900), mz = 335)
    checkTrue(all(rtime(chrs[[1]]) >= 2700 & rtime(chrs[[1]]) <= 2900))
    checkTrue(all(rtime(chrs[[2]]) >= 2700 & rtime(chrs[[2]]) <= 2900))
    checkTrue(all(rtime(chrs[[3]]) >= 2700 & rtime(chrs[[3]]) <= 2900))
    spctr <- spectra(filterMz(filterRt(od_x, rt = c(2700, 2900)), mz = 335))
    ints <- split(unlist(lapply(spctr, function(z) {
        if (z@peaksCount)
            return(sum(intensity(z)))
        else return(NA)
    })), f = unlist(lapply(spctr, fromFile)))
    checkEquals(ints[[1]], intensity(chrs[[1]]))
    checkEquals(ints[[2]], intensity(chrs[[2]]))
    checkEquals(ints[[3]], intensity(chrs[[3]]))
    ## Using adjusted rt: LLLL
    chrs <- extractChromatograms(xod_xgr, rt = c(2700, 2900), mz = 335)
    checkTrue(all(rtime(chrs[[1]]) >= 2700 & rtime(chrs[[1]]) <= 2900))
    checkTrue(all(rtime(chrs[[2]]) >= 2700 & rtime(chrs[[2]]) <= 2900))
    checkTrue(all(rtime(chrs[[3]]) >= 2700 & rtime(chrs[[3]]) <= 2900))
    spctr <- spectra(filterMz(filterRt(xod_xgr, rt = c(2700, 2900)), mz = 335))
    ints <- split(unlist(lapply(spctr, function(z) {
        if (z@peaksCount)
            return(sum(intensity(z)))
        else return(NA)
    })), f = unlist(lapply(spctr, fromFile)))
    checkEquals(ints[[1]], intensity(chrs[[1]]))
    checkEquals(ints[[2]], intensity(chrs[[2]]))
    checkEquals(ints[[3]], intensity(chrs[[3]]))
    ## Check the rtime.
    tmp <- filterRt(xod_xgr, rt = c(2700, 2900))
    checkEquals(rtime(chrs[[1]]), rtime(tmp, bySample = TRUE)[[1]])
    checkEquals(rtime(chrs[[2]]), rtime(tmp, bySample = TRUE)[[2]])
    checkEquals(rtime(chrs[[3]]), rtime(tmp, bySample = TRUE)[[3]])
    
    ## What if we're completely off?
    chrs <- extractChromatograms(od_x, rt = c(5000, 5500))
    checkTrue(length(chrs) == 0)
    chrs <- extractChromatograms(od_x, rt = c(2600, 2700), mz = 12000)
    checkTrue(length(chrs) == 0)
}

dontrun_test_with_MRM <- function() {
    ## Test how we could read the data.
    ## chromatogramsInfo
    library(msdata)
    fls <- proteomics(full.names = TRUE)

    library(mzR)
    msf <- mzR::openMSfile(fls[2], "pwiz")
    chrs <- chromatograms(msf)
    chrsI <- chromatogram(msf)
    ## The same essentially.
    nChrom(msf)
    length(chrs)
    nrow(chrs[[1]])
    mzR::close(msf)
    ## 
    msf <- mzR::openMSfile(fls[1], "pwiz")
    chrs <- chromatograms(msf)
    chrs <- chromatograms(msf)
    nChrom(msf)
    length(chrs)
    nrow(chrs[[1]])

    ## Now, we've got the following info: cvParam
    ## accession="MS:1000235" name="total ion current chromatogram" value=""
    ## Check http://proteowizard.sourceforge.net/dox/namespacepwiz_1_1msdata.html
    ## Potentially interesting:
    ## o ChromatogramIdentity nope, no header info.
    ## OK, have to look for chromatogram with index="1", then within <precursor>
    ## for cvParam accession="MS:1000827" and its value -> Q1 or precursorMz
    ## then within <product> for cvParam accession="MS:1000827" and its value
    ## -> Q3.

    ## https://sourceforge.net/p/proteowizard/mailman/message/27571266/
}
