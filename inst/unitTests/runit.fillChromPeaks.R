test_fillChromPeaks <- function() {
    ## No adjusted retention times
    checkTrue(!xcms:::.hasFilledPeaks(xod_xg))
    res <- fillChromPeaks(xod_xg)
    checkTrue(xcms:::.hasFilledPeaks(res))
    ph <- processHistory(res, type = xcms:::.PROCSTEP.PEAK.FILLING)
    checkTrue(length(ph) == 1)
    checkEquals(ph[[1]]@param, FillChromPeaksParam())
    ## Check if the signal corresponds to what we expect for some peaks.
    fp <- chromPeaks(res)
    fp <- fp[fp[, "is_filled"] == 1, ]
    idxs <- sample(1:nrow(fp), 5)
    for (i in idxs) {
        cfp <- fp[i, , drop = FALSE]
        tmp <- filterFile(xod_xg, file = cfp[1, "sample"])
        chr <- extractChromatograms(tmp, rt = cfp[1, c("rtmin", "rtmax")],
                                    mz = cfp[1, c("mzmin", "mzmax")])[[1]]
        into <- sum(intensity(chr), na.rm = TRUE) *
            (cfp[1, "rtmax"] - cfp[1, "rtmin"]) / (length(chr) - 1)
        checkEquals(unname(into), unname(cfp[1, "into"]))
    }
    ## Plot the data for some...
    if (FALSE) {
        pk_idx <- groupval(res)[1, ]
        pks <- chromPeaks(res)[pk_idx, ]
        rtr <- c(min(pks[, "rtmin"]), max(pks[, "rtmax"]))
        rtr[1] <- rtr[1] - 10
        rtr[2] <- rtr[2] + 10
        chrs <- extractChromatograms(res, rt = rtr, mz = c(min(pks[, "mzmin"]),
                                                           max(pks[, "mzmax"])))
        plot(3, 3, pch = NA, xlim = range(lapply(chrs, rtime), na.rm = TRUE),
             ylim = range(lapply(chrs, intensity), na.rm = TRUE), xlab = "rt",
             ylab = "int")
        for (i in 1:length(chrs)) {
            points(rtime(chrs[[i]]), intensity(chrs[[i]]), type = "l",
                   col = ifelse(pks[i, "is_filled"], yes = "red", no = "black"))
            abline(v = pks[i, c("rtmin", "rtmax")],
                   col = ifelse(pks[i, "is_filled"], yes = "red", no = "black"))
        }
    }

    ## Check for the NAs if there is really no signal
    gv <- groupval(res)
    feat_i <- which(is.na(gv[, 1]))
    tmp <- chromPeaks(res)[featureDefinitions(res)$peakidx[[feat_i]],
                           c("rtmin", "rtmax", "mzmin", "mzmax")]
    ## Get the intensities for the first one.
    pkArea <- apply(tmp, median, MARGIN = 2)
    chr <- extractChromatograms(res, rt = pkArea[1:2], mz = pkArea[3:4])
    checkTrue(length(chr) == 0)
    ## Get also the spectra:
    spctr <- spectra(filterRt(filterFile(xod_xg, file = 1), rt = pkArea[1:2]))
    mzs <- unlist(lapply(spctr, mz))
    ## No spectra for the fiven mz:
    checkEquals(sum(mzs >= pkArea[3] & mzs <= pkArea[4]), 0)
    
    ## Check increasing the expandRt and expandMz to see whether we get rid of
    ## the NA.
    res_2 <- fillChromPeaks(xod_xg, param = FillChromPeaksParam(expandMz = 1))
    ## Check if the mzrange is now indeed broader for the integrated ones.
    fp <- chromPeaks(res)
    fp <- fp[fp[, "is_filled"] == 1, ]
    fp2 <- chromPeaks(res_2)
    fp2 <- fp2[fp2[, "is_filled"] == 1, ]
    checkEquals(fp2[, "mzmax"] - fp2[, "mzmin"],
                2 * (fp[, "mzmax"] - fp[, "mzmin"]))
    
    res_2 <- fillChromPeaks(xod_xg, param = FillChromPeaksParam(expandRt = 1))
    ## Check if the mzrange is now indeed broader for the integrated ones.
    fp <- chromPeaks(res)
    fp <- fp[fp[, "is_filled"] == 1, ]
    fp2 <- chromPeaks(res_2)
    fp2 <- fp2[fp2[, "is_filled"] == 1, ]
    checkEquals(fp2[, "rtmax"] - fp2[, "rtmin"],
                2 * (fp[, "rtmax"] - fp[, "rtmin"]))
    ## Check using ppm
    res_2 <- fillChromPeaks(xod_xg, param = FillChromPeaksParam(ppm = 40,
                                                                expandMz = 5,
                                                                expandRt = 2))
    checkTrue(all(!is.na(rowSums(groupval(res_2)))))
    ## Drop them.
    res_rem <- dropFilledChromPeaks(res)
    checkTrue(!xcms:::.hasFilledPeaks(res_rem))
    checkEquals(res_rem, xod_xg)
    ## Drop feature definitions from res -> also filled peaks should be dropped.
    res_rem <- dropFeatureDefinitions(res)
    checkTrue(!xcms:::.hasFilledPeaks(res_rem))
    checkTrue(!any(chromPeaks(res_rem)[, "is_filled"] == 1))
    checkEquals(res_rem, xod_x)
    
    ## With adjusted rtime.
    res_2 <- fillChromPeaks(xod_xgrg)
    ## Check if the signal corresponds to what we expect for some peaks.
    fp <- chromPeaks(res_2)
    fp <- fp[fp[, "is_filled"] == 1, ]
    ## These have to be different from before!
    fp_raw <- chromPeaks(res)
    fp_raw <- fp_raw[fp_raw[, "is_filled"] == 1, ]
    checkTrue(all(fp_raw[, "rt"] != fp[, "rt"]))
    checkTrue(all(fp_raw[, "rtmin"] != fp[, "rtmin"]))
    checkTrue(all(fp_raw[, "rtmax"] != fp[, "rtmax"]))
    checkEquals(fp_raw[, "mz"], fp[, "mz"])
    checkEquals(fp_raw[, "mzmin"], fp[, "mzmin"])
    checkEquals(fp_raw[, "mzmax"], fp[, "mzmax"])
    ## Values are expected to be different, but still correlated!
    checkTrue(all(fp_raw[, "into"] != fp[, "into"]))
    checkTrue(cor(fp_raw[, "into"], fp[, "into"]) > 0.99)
    ## Check if we can get the same data using the provided range.
    ## Use the .rawMat function
    first <- filterFile(xod_xgrg, file = 1, keepAdjustedRtime = TRUE)
    spctr <- spectra(first)
    mzs <- lapply(spctr, mz)
    vps <- lengths(mzs)
    ints <- unlist(lapply(spctr, intensity), use.names = FALSE)
    mzs <- unlist(mzs, use.names = FALSE)
    rtim <- rtime(first)
    idx <- which(fp[, "sample"] == 1)
    for (i in idx) {
        mtx <- xcms:::.rawMat(mz = mzs, int = ints, scantime = rtim,
                              valsPerSpect = vps,
                              rtrange = fp[i, c("rtmin", "rtmax")],
                              mzrange = fp[i, c("mzmin", "mzmax")])
        into <- sum(mtx[, 3], na.rm = TRUE) *
            ((fp[i, "rtmax"] - fp[i, "rtmin"]) /
             (sum(rtim >= fp[i, "rtmin"] & rtim <= fp[i, "rtmax"]) - 1))
        checkEquals(unname(into), unname(fp[i, "into"]))
    }

    ## Drop them.
    res_rem <- dropFilledChromPeaks(res_2)
    checkTrue(!xcms:::.hasFilledPeaks(res_rem))
    checkEquals(res_rem, xod_xgrg)

    ## ## Alternative without rawMat:
    ## for (i in idx) {
    ##     chrs <- extractChromatograms(first, rt = fp[i, c("rtmin", "rtmax")],
    ##                                  mz = fp[i, c("mzmin", "mzmax")])[[1]]
    ##     into <- sum(intensity(chrs), na.rm = TRUE) / (length(chrs) - 1) *
    ##         (fp[i, "rtmax"] - fp[i, "rtmin"])
    ##     checkEquals(unname(into), unname(fp[i, "into"]))
    ## }
}

dontrun_exhaustive_fillChromPeaks_test <- function() {
    fls <- c("/Users/jo/data/2016/2016-11/NoSN/190516_POOL_N_POS_15.mzML",
             "/Users/jo/data/2016/2016-11/NoSN/190516_POOL_N_POS_19.mzML",
             "/Users/jo/data/2016/2016-11/NoSN/190516_POOL_N_POS_11.mzML")
    raw <- readMSData2(fls)
    pks <- findChromPeaks(raw, param = CentWaveParam(peakwidth = c(0.8, 20),
                                                     ppm = 40))
    pks_noRt <- groupChromPeaks(pks, param = PeakDensityParam(minFraction = 0.6))
    filled <- fillChromPeaks(pks_noRt)
    fp <- chromPeaks(filled)
    fp <- fp[fp[, "is_filled"] == 1, ]
    idxs <- sample(1:nrow(fp), 20)
    for (i in idxs) {
        cfp <- fp[i, , drop = FALSE]
        tmp <- filterFile(pks_noRt, file = cfp[1, "sample"],
                          keepAdjustedRtime = TRUE)
        chr <- extractChromatograms(tmp, rt = cfp[1, c("rtmin", "rtmax")],
                                    mz = cfp[1, c("mzmin", "mzmax")])[[1]]
        ## into <- sum(intensity(chr), na.rm = TRUE) *
        ##     (cfp[1, "rtmax"] - cfp[1, "rtmin"]) /
        ##     (sum(rtime(tmp) >= cfp[1, "rtmin"] & rtime(tmp) <= cfp[1, "rtmax"]) - 1)
        into <- sum(intensity(chr), na.rm = TRUE) *
            (cfp[1, "rtmax"] - cfp[1, "rtmin"]) / (length(chr) - 1)
        checkEquals(unname(into), unname(cfp[1, "into"]))
    }
    ## Check those with an NA.
    gv <- groupval(filled)
    with_na <- is.na(rowSums(gv))
    idxs <- sample(which(with_na), 20)
    for (i in idxs) {
        tmp <- chromPeaks(pks_noRt)[featureDefinitions(pks_noRt)$peakidx[[i]],
                                    c("rtmin", "rtmax", "mzmin", "mzmax")]
        ## Get the intensities for the first one.
        pkArea <- apply(tmp, median, MARGIN = 2)
        smpl <- which(is.na(gv[i, ]))
        tmp <- filterFile(pks_noRt, file = smpl, keepAdjustedRtime = TRUE)
        chr <- extractChromatograms(tmp, rt = pkArea[1:2], mz = pkArea[3:4])
        checkTrue(length(chr) == 0)
    }
    still_missing <- is.na(rowSums(groupval(filled)))
    ## Try using ppm:
    filled <- fillChromPeaks(pks_noRt, param = FillChromPeaksParam(ppm = 40))
    checkTrue(sum(still_missing) > sum(is.na(rowSums(groupval(filled)))))
    filled <- fillChromPeaks(pks_noRt, param = FillChromPeaksParam(ppm = 40,
                                                                   expandMz = 2))
    checkTrue(sum(still_missing) > sum(is.na(rowSums(groupval(filled)))))
    
    ## With adjusted retention times.
    pks <- adjustRtime(pks, param = ObiwarpParam())
    pks <- groupChromPeaks(pks, param = PeakDensityParam(minFraction = 0.6))
    filled <- fillChromPeaks(pks)
    ##
    fp <- chromPeaks(filled)
    fp <- fp[fp[, "is_filled"] == 1, ]
    idxs <- sample(1:nrow(fp), 20)
    for (i in idxs) {
        cfp <- fp[i, , drop = FALSE]
        tmp <- filterFile(pks, file = cfp[1, "sample"], keepAdjustedRtime = TRUE)
        chr <- extractChromatograms(tmp, rt = cfp[1, c("rtmin", "rtmax")],
                                    mz = cfp[1, c("mzmin", "mzmax")])[[1]]
        ## into <- sum(intensity(chr), na.rm = TRUE) *
        ##     (cfp[1, "rtmax"] - cfp[1, "rtmin"]) /
        ##     (sum(rtime(tmp) >= cfp[1, "rtmin"] & rtime(tmp) <= cfp[1, "rtmax"]) - 1)
        into <- sum(intensity(chr), na.rm = TRUE) *
            (cfp[1, "rtmax"] - cfp[1, "rtmin"]) / (length(chr) - 1)
        checkEquals(unname(into), unname(cfp[1, "into"]))
    }
    ## Check those with an NA.
    gv <- groupval(filled)
    with_na <- is.na(rowSums(gv))
    idxs <- sample(which(with_na), 20)
    for (i in idxs) {
        tmp <- chromPeaks(pks)[featureDefinitions(pks)$peakidx[[i]],
                               c("rtmin", "rtmax", "mzmin", "mzmax")]
        ## Get the intensities for the first one.
        pkArea <- apply(tmp, median, MARGIN = 2)
        smpl <- which(is.na(gv[i, ]))
        tmp <- filterFile(pks, file = smpl, keepAdjustedRtime = TRUE)
        chr <- extractChromatograms(tmp, rt = pkArea[1:2], mz = pkArea[3:4])
        checkTrue(length(chr) == 0)
    }
    
    ## what is the mz for those with NA against those without?
    boxplot(list(filled = featureDefinitions(filled)[!with_na, "mzmed"],
                 failed = featureDefinitions(filled)[with_na, "mzmed"]),
            varwidth = TRUE)
    ## They are about the same, no difference.
    
    ## what is the intensity for those with NA against those without?
    maxo <- rowMeans(groupval(filled, value = "maxo"), na.rm = TRUE)
    boxplot(list(filled = log2(maxo[!with_na]), failed = log2(maxo[with_na])),
            varwidth = TRUE)
    ## No difference.
}

dontrun_test_getPeakInt_functions <- function() {
    ## Testing whether the .getPeakInt functions are correct and which one is
    ## more performant.
    tmp <- filterFile(xod_xgrg, file = 3)
    pkInt <- xcms:::.getPeakInt(tmp, chromPeaks(tmp))
    pkInt2 <- xcms:::.getPeakInt2(tmp, chromPeaks(tmp))
    pkInt3 <- xcms:::.getPeakInt3(tmp, chromPeaks(tmp))
    checkEquals(pkInt, pkInt2)
    checkEquals(pkInt, pkInt3)
    checkEquals(pkInt, chromPeaks(tmp)[, "into"])
    checkEquals(pkInt2, chromPeaks(tmp)[, "into"])
    checkEquals(pkInt3, chromPeaks(tmp)[, "into"])

    library(microbenchmark)
    microbenchmark(xcms:::.getPeakInt(tmp, chromPeaks(tmp)[1, , drop = FALSE]),
                   xcms:::.getPeakInt2(tmp, chromPeaks(tmp)[1, , drop = FALSE]),
                   xcms:::.getPeakInt3(tmp, chromPeaks(tmp)[1, , drop = FALSE]),
                   times = 10)
    ## 258 ms vs 503 ms vs 499 ms 
    microbenchmark(xcms:::.getPeakInt(tmp, chromPeaks(tmp)[1:5, ]),
                   xcms:::.getPeakInt2(tmp, chromPeaks(tmp)[1:5, ]),
                   xcms:::.getPeakInt3(tmp, chromPeaks(tmp)[1:5, ]),
                   times = 10)
    ## 1269 ms vs 586 ms vs 587 ms
    microbenchmark(xcms:::.getPeakInt(tmp, chromPeaks(tmp)[1:10, ]),
                   xcms:::.getPeakInt2(tmp, chromPeaks(tmp)[1:10, ]),
                   xcms:::.getPeakInt3(tmp, chromPeaks(tmp)[1:10, ]),
                   times = 10)
    ## 2577 ms vs 594 ms vs 562 ms
    microbenchmark(xcms:::.getPeakInt(tmp, chromPeaks(tmp)),
                   xcms:::.getPeakInt2(tmp, chromPeaks(tmp)),
                   xcms:::.getPeakInt3(tmp, chromPeaks(tmp)),
                   times = 10)
    ## 16447 ms vs 676 ms vs 556 ms
    ## Well. getPeakInt2 and getPeakInt3 are considerably faster!
}

## This provides some extensive tests.
dontrun_getPeakInt_validity <- function() {
    ## Do extensive tests on the .getPeakInt3 function to ensure it is really
    ## returning what it should.
    ## faahKO centWave peaks
    tmp <- filterFile(xod_xgrg, file = 2)
    pkInt2 <- xcms:::.getPeakInt2(tmp, chromPeaks(tmp))
    checkEquals(pkInt2, chromPeaks(tmp)[, "into"])

    ## faahKO matchedFilter peaks
    tmp <- findChromPeaks(od_x, param = MatchedFilterParam())
    tmp <- filterFile(tmp, file = 1)
    pkInt2 <- xcms:::.getPeakInt2(tmp, chromPeaks(tmp))
    checkEquals(pkInt2, chromPeaks(tmp)[, "into"])

    ## own file centWave peaks
    fl <- "/Users/jo/data/2016/2016-11/NoSN/190516_POOL_N_POS_19.mzML"
    raw_x <- readMSData2(fl)
    ## Default centWave - completely off.
    tmp <- findChromPeaks(raw_x, param = CentWaveParam(verboseColumns = TRUE))
    pkInt3 <- xcms:::.getPeakInt3(tmp, chromPeaks(tmp))
    checkEquals(pkInt3, chromPeaks(tmp)[, "into"])

    ## With OK settings
    options(originalCentWave = FALSE)
    tmp <- findChromPeaks(raw_x, param = CentWaveParam(peakwidth = c(0.8, 40)))
    pkInt3 <- xcms:::.getPeakInt3(tmp, chromPeaks(tmp))
    checkEquals(pkInt3, chromPeaks(tmp)[, "into"])

    library(microbenchmark)
    microbenchmark(xcms:::.getPeakInt3(tmp, chromPeaks(tmp)),
                   xcms:::.getPeakInt2(tmp, chromPeaks(tmp)), times = 10)
    ## 
    
    ## Reproduce with msdata files:
    fl <- system.file("microtofq/MM14.mzML", package = "msdata")
    raw <- readMSData2(fl)
    tmp <- findChromPeaks(raw, param = CentWaveParam(peakwidth = c(1, 20)))
    pkInt3 <- xcms:::.getPeakInt3(tmp, chromPeaks(tmp))
    checkEquals(pkInt3, chromPeaks(tmp)[, "into"])
    
    ## own file matchedFilter peaks
    ## TODO: continue here - but be aware, for MatchedFilterParam we might have
    ## to do the stuff on the profile matrix!
    tmp <- findChromPeaks(raw_x, param = MatchedFilterParam())
}
