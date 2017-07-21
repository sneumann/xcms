## Unit tests related to the Chromatogram class.
test_deprecated_Chromatogram <- function() {
    chrs <- extractChromatograms(filterFile(od_x, file = 2))
    plotChromatogram(chrs)
}


test_extractChromatograms <- function() {
    ## OnDiskMSnExp
    ## TIC
    chrs <- chromatogram(filterFile(od_x, file = 2))
    plot(chrs)
    spctr <- spectra(filterFile(od_x, file = 2))
    ints <- unlist(lapply(spctr, function(z)
        return(sum(intensity(z)))))
    checkEquals(intensity(chrs[1, 1]), ints)
    checkEquals(rtime(chrs[1, 1]), unlist(lapply(spctr, rtime)))
    ## BPC
    chrs <- chromatogram(filterFile(od_x, file = 2),
                         aggregationFun = "max")
    ints <- unlist(lapply(spctr, function(z)
        return(max(intensity(z)))))
    checkEquals(intensity(chrs[1, 1]), ints)
    checkEquals(rtime(chrs[1, 1]), unlist(lapply(spctr, rtime)))
    ## XCMSnExp
    xod_x <- faahko_xod
    chrs <- chromatogram(filterFile(xod_x, file = 2))
    ints <- unlist(lapply(spctr, function(z)
        return(sum(intensity(z)))))
    checkEquals(intensity(chrs[1, 1]), ints)
    checkEquals(rtime(chrs[1, 1]), unlist(lapply(spctr, rtime)))
    ## BPC
    chrs <- chromatogram(filterFile(xod_x, file = 2),
                         aggregationFun = "max")
    ints <- unlist(lapply(spctr, function(z)
        return(max(intensity(z)))))    
    checkEquals(intensity(chrs[1, 1]), ints)
    checkEquals(rtime(chrs[1, 1]), unlist(lapply(spctr, rtime)))
    ## with adjusted retention times.
    chrs <- chromatogram(filterFile(xod_xgr, file = 2),
                         adjustedRtime = FALSE, aggregationFun = "max")
    checkEquals(intensity(chrs[1, 1]), ints)
    checkEquals(rtime(chrs[1, 1]), unlist(lapply(spctr, rtime)))
    chrs <- chromatogram(filterFile(xod_xgr, file = 2,
                                    keepAdjustedRtime = TRUE),
                         aggregationFun = "max")
    checkEquals(intensity(chrs[[1]]), ints)
    checkEquals(rtime(chrs[1, 1]), rtime(xod_xgr, bySample = TRUE,
                                         adjusted = TRUE)[[2]])
    ## Subset to certain mz range in all files.
    chrs_adj <- chromatogram(xod_xgr, mz = c(300, 330))
    chrs_raw <- chromatogram(xod_x, mz = c(300, 330))
    checkTrue(sum(rtime(chrs_adj[1, 1]) != rtime(chrs_raw[1, 1])) >
              length(chrs_raw[1, 1]) / 2)
    checkEquals(rtime(chrs_adj[1, 1]), rtime(xod_xgr, bySample = TRUE)[[1]])
    checkEquals(rtime(chrs_adj[1, 2]), rtime(xod_xgr, bySample = TRUE)[[2]])
    checkEquals(rtime(chrs_adj[1, 3]), rtime(xod_xgr, bySample = TRUE)[[3]])
    
    ## Now subsetting for mz:
    tmp <- filterFile(od_x, file = 2)
    chrs <- chromatogram(tmp, mz = c(300, 400))
    checkEquals(mz(chrs[1, 1], filter = TRUE), c(300, 400))
    suppressWarnings(spctr <- spectra(filterMz(tmp, mz = c(300, 400))))
    ints <- unlist(lapply(spctr, function(z)
        return(sum(intensity(z)))))
    ints2 <- intensity(chrs[1, 1])
    ints2[is.na(ints2)] <- 0
    checkEquals(ints2, ints)
    checkEquals(rtime(chrs[1, 1]), unlist(lapply(spctr, rtime)))
    ## with adjusted retention times
    chrs <- chromatogram(filterFile(xod_xgr, file = 2,
                                    keepAdjustedRtime = TRUE),
                         mz = c(300, 400))
    ints <- unlist(lapply(spctr, function(z)
        return(sum(intensity(z)))))
    ints2 <- intensity(chrs[1, 1])
    ints2[is.na(ints2)] <- 0
    checkEquals(ints2, ints)
    checkEquals(rtime(chrs[1, 1]), rtime(xod_xgr, bySample = TRUE)[[2]])
    
    ## Now subsetting for rt:
    chrs <- chromatogram(od_x, rt = c(2700, 2900))
    checkTrue(all(rtime(chrs[1, 1]) >= 2700 & rtime(chrs[1, 1]) <= 2900))
    checkTrue(all(rtime(chrs[1, 2]) >= 2700 & rtime(chrs[1, 2]) <= 2900))
    checkTrue(all(rtime(chrs[1, 3]) >= 2700 & rtime(chrs[1, 3]) <= 2900))
    spctr <- spectra(filterRt(od_x, rt = c(2700, 2900)))
    ints <- split(unlist(lapply(spctr, function(z) sum(intensity(z)))),
                  f = unlist(lapply(spctr, fromFile)))
    checkEquals(ints[[1]], intensity(chrs[1, 1]))
    checkEquals(ints[[2]], intensity(chrs[1, 2]))
    checkEquals(ints[[3]], intensity(chrs[1, 3]))
    ## Using adjusted rt:
    chrs2 <- chromatogram(xod_xgr, rt = c(2700, 2900))
    checkTrue(all(rtime(chrs2[1, 1]) >= 2700 & rtime(chrs2[1, 1]) <= 2900))
    checkTrue(all(rtime(chrs2[1, 2]) >= 2700 & rtime(chrs2[1, 2]) <= 2900))
    checkTrue(all(rtime(chrs2[1, 3]) >= 2700 & rtime(chrs2[1, 3]) <= 2900))
    checkTrue(length(chrs[1, 1]) != length(chrs2[1, 1]))
    checkTrue(length(chrs[1, 2]) == length(chrs2[1, 2]))
    checkTrue(length(chrs[1, 3]) != length(chrs2[1, 3]))
    tmp <- filterRt(xod_xgr, rt = c(2700, 2900))
    checkEquals(rtime(chrs2[1, 1]), rtime(tmp, bySample = TRUE)[[1]])
    checkEquals(rtime(chrs2[1, 2]), rtime(tmp, bySample = TRUE)[[2]])
    checkEquals(rtime(chrs2[1, 3]), rtime(tmp, bySample = TRUE)[[3]])
    ## Check the values...
    keepSp <- which(adjustedRtime(xod_xgr) >= 2700 &
                    adjustedRtime(xod_xgr) <= 2900)
    tmp <- xod_xgr[keepSp]
    ints <- unlist(lapply(spectra(tmp), function(z) sum(intensity(z))))
    intsL <- split(ints, fromFile(tmp))
    checkEquals(intensity(chrs2[1, 1]), intsL[[1]])
    checkEquals(intensity(chrs2[1, 2]), intsL[[2]])
    checkEquals(intensity(chrs2[1, 3]), intsL[[3]])
    
    ## Now subsetting for rt and mz:
    chrs <- chromatogram(od_x, rt = c(2700, 2900), mz = 335)
    checkTrue(all(rtime(chrs[1, 1]) >= 2700 & rtime(chrs[1, 1]) <= 2900))
    checkTrue(all(rtime(chrs[1, 2]) >= 2700 & rtime(chrs[1, 2]) <= 2900))
    checkTrue(all(rtime(chrs[1, 3]) >= 2700 & rtime(chrs[1, 3]) <= 2900))
    spctr <- spectra(filterMz(filterRt(od_x, rt = c(2700, 2900)), mz = 335))
    ints <- split(unlist(lapply(spctr, function(z) {
        if (z@peaksCount)
            return(sum(intensity(z)))
        else return(NA)
    })), f = unlist(lapply(spctr, fromFile)))
    checkEquals(ints[[1]], intensity(chrs[1, 1]))
    checkEquals(ints[[2]], intensity(chrs[1, 2]))
    checkEquals(ints[[3]], intensity(chrs[1, 3]))
    ## Using adjusted rt: LLLL
    chrs <- chromatogram(xod_xgr, rt = c(2700, 2900), mz = 335)
    checkTrue(all(rtime(chrs[1, 1]) >= 2700 & rtime(chrs[1, 1]) <= 2900))
    checkTrue(all(rtime(chrs[1, 2]) >= 2700 & rtime(chrs[1, 2]) <= 2900))
    checkTrue(all(rtime(chrs[1, 3]) >= 2700 & rtime(chrs[1, 3]) <= 2900))
    spctr <- spectra(filterMz(filterRt(xod_xgr, rt = c(2700, 2900)), mz = 335))
    ints <- split(unlist(lapply(spctr, function(z) {
        if (z@peaksCount)
            return(sum(intensity(z)))
        else return(NA)
    })), f = unlist(lapply(spctr, fromFile)))
    checkEquals(ints[[1]], intensity(chrs[1, 1]))
    checkEquals(ints[[2]], intensity(chrs[1, 2]))
    checkEquals(ints[[3]], intensity(chrs[1, 3]))
    ## Check the rtime.
    tmp <- filterRt(xod_xgr, rt = c(2700, 2900))
    checkEquals(rtime(chrs[1, 1]), rtime(tmp, bySample = TRUE)[[1]])
    checkEquals(rtime(chrs[1, 2]), rtime(tmp, bySample = TRUE)[[2]])
    checkEquals(rtime(chrs[1, 3]), rtime(tmp, bySample = TRUE)[[3]])
    
    ## What if we're completely off?
    chrs <- chromatogram(od_x, rt = c(5000, 5500))
    checkTrue(nrow(chrs) == 0)
    ## Now rt is within range, but mz is completely off. We expect Chromatograms
    ## with same length than there are spectra in the rt range, but all NA
    ## values.
    chrs <- chromatogram(od_x, rt = c(2600, 2700), mz = 12000)
    rts <- split(rtime(od_x), f = fromFile(od_x))
    rts <- lapply(rts, function(z) z[z >= 2600 & z <= 2700])
    checkEquals(unname(lengths(chrs[1, , drop = TRUE])), unname(lengths(rts)))
    ## All have to be NA.
    checkTrue(all(unlist(lapply(chrs[1, ], function(z) is.na(intensity(z))))))

    ## Multiple ranges.
    rtr <- matrix(c(2700, 2900, 2600, 2800), ncol = 2, byrow = TRUE)
    mzr <- matrix(c(355, 355, 344, 344), ncol = 2, byrow = TRUE)
    chrs <- chromatogram(od_x, rt = rtr, mz = mzr)
    
    checkTrue(all(rtime(chrs[1, 1]) >= 2700 & rtime(chrs[1, 1]) <= 2900))
    checkTrue(all(rtime(chrs[1, 2]) >= 2700 & rtime(chrs[1, 2]) <= 2900))
    checkTrue(all(rtime(chrs[1, 3]) >= 2700 & rtime(chrs[1, 3]) <= 2900))
    checkTrue(all(rtime(chrs[2, 1]) >= 2600 & rtime(chrs[2, 1]) <= 2800))
    checkTrue(all(rtime(chrs[2, 2]) >= 2600 & rtime(chrs[2, 2]) <= 2800))
    checkTrue(all(rtime(chrs[2, 3]) >= 2600 & rtime(chrs[2, 3]) <= 2800))
    spctr <- spectra(filterMz(filterRt(od_x, rt = rtr[1, ]),
                              mz = mzr[1, ]))
    ints <- split(unlist(lapply(spctr, function(z) {
        if (z@peaksCount)
            return(sum(intensity(z)))
        else return(NA)
    })), f = unlist(lapply(spctr, fromFile)))
    checkEquals(ints[[1]], intensity(chrs[1, 1]))
    checkEquals(ints[[2]], intensity(chrs[1, 2]))
    checkEquals(ints[[3]], intensity(chrs[1, 3]))
    spctr <- spectra(filterMz(filterRt(od_x, rt = rtr[2, ]),
                              mz = mzr[2, ]))
    ints <- split(unlist(lapply(spctr, function(z) {
        if (z@peaksCount)
            return(sum(intensity(z)))
        else return(NA)
    })), f = unlist(lapply(spctr, fromFile)))
    checkEquals(ints[[1]], intensity(chrs[2, 1]))
    checkEquals(ints[[2]], intensity(chrs[2, 2]))
    checkEquals(ints[[3]], intensity(chrs[2, 3]))

    ## Multiple ranges with complete off ranges.
    rtr <- matrix(c(2700, 2900, 5000, 5500, 2600, 2800), ncol = 2, byrow = TRUE)
    mzr <- matrix(c(355, 355, 500, 500, 344, 344), ncol = 2, byrow = TRUE)
    chrs <- chromatogram(od_x, rt = rtr, mz = mzr)
    checkTrue(nrow(chrs) == 3)
    checkTrue(all(lengths(chrs[2, ]) == 0))
    
    rtr <- matrix(c(2700, 2900, 2700, 2900, 2600, 2800), ncol = 2, byrow = TRUE)
    mzr <- matrix(c(355, 355, 100000, 100000, 344, 344), ncol = 2, byrow = TRUE)
    chrs <- chromatogram(od_x, rt = rtr, mz = mzr)
    checkTrue(nrow(chrs) == 3)
    ## All values in the 2nd Chromosome object have to be NA.
    checkTrue(all(unlist(lapply(chrs[2, ], function(z) is.na(intensity(z))))))
}

## dontrun_test_with_MRM <- function() {
##     ## Test how we could read the data.
##     ## chromatogramsInfo
##     library(msdata)
##     fls <- proteomics(full.names = TRUE)

##     library(mzR)
##     msf <- mzR::openMSfile(fls[2], "pwiz")
##     chrs <- chromatograms(msf)
##     chrsI <- chromatogram(msf)
##     ## The same essentially.
##     nChrom(msf)
##     length(chrs)
##     nrow(chrs[[1]])
##     mzR::close(msf)
##     ## 
##     msf <- mzR::openMSfile(fls[1], "pwiz")
##     chrs <- chromatograms(msf)
##     chrs <- chromatograms(msf)
##     nChrom(msf)
##     length(chrs)
##     nrow(chrs[[1]])

##     ## Now, we've got the following info: cvParam
##     ## accession="MS:1000235" name="total ion current chromatogram" value=""
##     ## Check http://proteowizard.sourceforge.net/dox/namespacepwiz_1_1msdata.html
##     ## Potentially interesting:
##     ## o ChromatogramIdentity nope, no header info.
##     ## OK, have to look for chromatogram with index="1", then within <precursor>
##     ## for cvParam accession="MS:1000827" and its value -> Q1 or precursorMz
##     ## then within <product> for cvParam accession="MS:1000827" and its value
##     ## -> Q3.

##     ## https://sourceforge.net/p/proteowizard/mailman/message/27571266/
## }
