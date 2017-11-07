test_fillChromPeaks <- function() {
    ## No adjusted retention times
    checkTrue(!xcms:::.hasFilledPeaks(xod_xg))
    res <- fillChromPeaks(xod_xg)
    checkTrue(xcms:::.hasFilledPeaks(res))
    ph <- processHistory(res, type = xcms:::.PROCSTEP.PEAK.FILLING)
    checkTrue(length(ph) == 1)
    checkEquals(ph[[1]]@param, FillChromPeaksParam())
    ## Check parameter filled in featureValues (issue #157)
    checkEquals(featureValues(res, filled = FALSE), featureValues(xod_xg))
    
    ## Check if the signal corresponds to what we expect for some peaks.
    fp <- chromPeaks(res)
    fp <- fp[fp[, "is_filled"] == 1, ]
    idxs <- sample(1:nrow(fp), 5)
    for (i in idxs) {
        cfp <- fp[i, , drop = FALSE]
        tmp <- filterFile(xod_xg, file = cfp[1, "sample"])
        chr <- chromatogram(tmp, rt = cfp[1, c("rtmin", "rtmax")],
                            mz = cfp[1, c("mzmin", "mzmax")])[1, 1]
        into <- sum(intensity(chr), na.rm = TRUE) *
            (cfp[1, "rtmax"] - cfp[1, "rtmin"]) / (length(chr) - 1)
        checkEquals(unname(into), unname(cfp[1, "into"]))
    }
    ## Plot the data for some...
    if (FALSE) {
        pk_idx <- featureValues(res)[1, ]
        pks <- chromPeaks(res)[pk_idx, ]
        rtr <- c(min(pks[, "rtmin"]), max(pks[, "rtmax"]))
        rtr[1] <- rtr[1] - 10
        rtr[2] <- rtr[2] + 10
        chrs <- chromatogram(res, rt = rtr, mz = c(min(pks[, "mzmin"]),
                                                   max(pks[, "mzmax"])))[1, ]
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

    ## Check if the results are similar that we get with findChromPeaks
    for (i in 1:length(fileNames(xod_xg))) {
        fnd_pks <- chromPeaks(xod_xg)[chromPeaks(xod_xg)[, "sample"] == i, ]
        prm <- processHistory(tmp, type ="Peak detection")[[1]]@param
        ## Extract the data for these using the internal function.
        fld_pks <- xcms:::.getChromPeakData(filterFile(xod_xg, i),
                                            peakArea = fnd_pks,
                                            sample_idx = i,
                                            cn = colnames(fnd_pks))
        ## rt
        checkTrue(cor(fnd_pks[, "rt"], fld_pks[, "rt"]) > 0.99)
        ## mz
        checkTrue(cor(fnd_pks[, "mz"], fld_pks[, "mz"]) > 0.99)
        checkEquals(fnd_pks[, "mz"], fld_pks[, "mz"])
        ## into
        checkTrue(cor(fnd_pks[, "into"], fld_pks[, "into"]) > 0.99)
        checkEquals(fnd_pks[, "into"], fld_pks[, "into"])
        ## checkEquals(fnd_pks[, "into"], fld_pks[, "into"])    
        ## maxo
        checkEquals(fnd_pks[, "maxo"], fld_pks[, "maxo"])    
        checkEquals(fnd_pks[, "maxo"], fld_pks[, "maxo"])
    }
    
    ## Check for the NAs if there is really no signal
    gv <- featureValues(res)
    feat_i <- which(is.na(gv[, 1]))
    tmp <- chromPeaks(res)[featureDefinitions(res)$peakidx[[feat_i]],
                           c("rtmin", "rtmax", "mzmin", "mzmax")]
    ## Get the intensities for the first one.
    pkArea <- apply(tmp, median, MARGIN = 2)
    chr <- chromatogram(res, rt = pkArea[1:2], mz = pkArea[3:4])[1, ]
    checkTrue(all(unlist(lapply(chr, function(z) is.na(intensity(z))))))
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
    checkTrue(all(!is.na(rowSums(featureValues(res_2)))))
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

## fillChromPeaks for MSW peak detection.
test_fillChromPeaks_MSW <- function() {
    p <- MzClustParam()
    fticr_xodg <- groupChromPeaks(fticr_xod, param = p)
    checkException(res <- fillChromPeaks(fticr_xod))
    res <- fillChromPeaks(fticr_xodg)
    
    ## Got a signal for all of em.
    checkTrue(!any(is.na(featureValues(res))))
    ## 1) Compare with what I get for xcmsSet.
    tmp_x <- fticr_xs
    tmp_x <- group(tmp_x, method = "mzClust")
    tmp_x <- fillPeaks(tmp_x, method = "MSW")
    ## Compare
    checkEquals(unname(groupval(tmp_x)), unname(featureValues(res)))
    checkEquals(unname(groupval(tmp_x, value = "maxo")),
                unname(featureValues(res, value = "maxo")))
    checkEquals(unname(groupval(tmp_x, value = "into")),
                unname(featureValues(res, value = "into")))
    checkEquals(unname(groupval(tmp_x, value = "mz")),
                unname(featureValues(res, value = "mz")))
    checkEquals(unname(groupval(tmp_x, value = "mzmin")),
                unname(featureValues(res, value = "mzmin")))
    checkEquals(unname(groupval(tmp_x, value = "mzmax")),
                unname(featureValues(res, value = "mzmax")))
    ## OK
    ## 2) Check if the fillChromPeaks returns same/similar data than the
    ##    findChromPeaks does:
    fdef <- featureDefinitions(fticr_xodg)
    pkArea <- do.call(
        rbind,
        lapply(
            fdef$peakidx, function(z) {
                tmp <- chromPeaks(fticr_xodg)[z, c("rtmin", "rtmax",
                                                   "mzmin", "mzmax"),
                                              drop = FALSE]
                pa <- c(median(tmp[, 1]), median(tmp[, 2]),
                        median(tmp[, 3]), median(tmp[, 4]))
                return(pa)
            }
        ))
    colnames(pkArea) <- c("rtmin", "rtmax", "mzmin", "mzmax")
    pkArea <- cbind(group_idx = 1:nrow(pkArea), pkArea,
                    mzmed = fdef$mzmed)
    ## Get peak data for all peaks in the first file
    allPks <- xcms:::.getMSWPeakData(filterFile(fticr_xodg, file = 1),
                                     peakArea = pkArea,
                                     sample_idx = 1,
                                     cn = colnames(chromPeaks(fticr_xodg)))
    curP <- chromPeaks(res)[chromPeaks(res)[, "sample"] == 1, ]
    curP <- curP[order(curP[, "mz"]), ]
    checkEquals(allPks[, "mz"], curP[, "mz"])
    checkEquals(allPks[, "maxo"], curP[, "maxo"])
    checkTrue(cor(allPks[, "into"], curP[, "into"]) > 0.99) ## Not exactly the
    ## same but highly similar.
}

test_fillChromPeaks_matchedFilter <- function() {
    tmp <- findChromPeaks(faahko_od, param = MatchedFilterParam())
    sg <- rep(1, length(fileNames(tmp)))
    tmp <- groupChromPeaks(tmp, param = PeakDensityParam(sampleGroups = sg))

    tmp_filled <- fillChromPeaks(tmp)
    checkTrue(sum(is.na(featureValues(tmp_filled))) <
              sum(is.na(featureValues(tmp))))
    nas <- is.na(featureValues(tmp)[, 1]) | is.na(featureValues(tmp)[, 2])
    checkTrue(cor(featureValues(tmp, value = "into")[!nas, 1],
                  featureValues(tmp, value = "into")[!nas, 2]) > 0.97)
    ## plot(featureValues(tmp, value = "into")[!nas, 1],
    ##      featureValues(tmp, value = "into")[!nas, 2])
    checkTrue(cor(featureValues(tmp_filled, value = "into")[, 1],
                  featureValues(tmp_filled, value = "into")[, 2],
                  use = "complete.obs") > 0.97)
    ## plot(featureValues(tmp_filled, value = "into")[, 1],
    ##      featureValues(tmp_filled, value = "into")[, 2])

    ## Check signal generation for already found peaks.
    for (i in 1:length(fileNames(tmp))) {
        fnd_pks <- chromPeaks(tmp)[chromPeaks(tmp)[, "sample"] == i, ]
        prm <- processHistory(tmp, type ="Peak detection")[[1]]@param
        ## Extract the data for these using the internal function.
        fld_pks <- xcms:::.getChromPeakData_matchedFilter(filterFile(tmp, i),
                                                          peakArea = fnd_pks,
                                                          sample_idx = i,
                                                          param = prm,
                                                          cn = colnames(fnd_pks))
        ## rt can not be the same, since for fillChromPeaks it is the rt of the
        ## maximum signal and for findChromPeaks it is the rt of the apex of the
        ## filtered/fitted peak.
        checkTrue(cor(fnd_pks[, "rt"], fld_pks[, "rt"]) > 0.99)
        ## mz: also not the same; most likely due to slightly different binning.
        diffs <- fnd_pks[, "mz"] - fld_pks[, "mz"]
        checkTrue(max(diffs) < 1e-4)
        ## into
        checkEquals(fnd_pks[, "into"], fld_pks[, "into"])    
        ## maxo
        checkEquals(fnd_pks[, "maxo"], fld_pks[, "maxo"])    
    }

    ## modify fillChromPeaks settings.
    tmp_fld_2 <- fillChromPeaks(
        tmp, param = FillChromPeaksParam(ppm = 40, expandRt = 1))
    checkTrue(sum(is.na(featureValues(tmp_filled))) <
              sum(is.na(featureValues(tmp))))
    checkTrue(sum(is.na(featureValues(tmp_fld_2))) <
              sum(is.na(featureValues(tmp_filled))))
    nas <- is.na(featureValues(tmp)[, 1]) | is.na(featureValues(tmp)[, 2])
    checkTrue(cor(featureValues(tmp_fld_2, value = "into")[, 1],
                  featureValues(tmp_fld_2, value = "into")[, 2],
                  use = "complete.obs") > 0.97)
}

dontrun_exhaustive_fillChromPeaks_matchedFilter <- function() {
    ## Different step sizes.
    prm <- MatchedFilterParam(binSize = 0.6)
    tmp <- findChromPeaks(faahko_od, param = prm)
    sg <- rep(1, length(fileNames(tmp)))
    tmp <- groupChromPeaks(tmp, param = PeakDensityParam(sampleGroups = sg))
    tmp_fld <- fillChromPeaks(tmp)
    checkTrue(sum(is.na(featureValues(tmp_fld))) <
              sum(is.na(featureValues(tmp))))
    nas <- is.na(featureValues(tmp)[, 1]) | is.na(featureValues(tmp)[, 2])
    checkTrue(cor(featureValues(tmp, value = "into")[!nas, 1],
                  featureValues(tmp, value = "into")[!nas, 2]) > 0.97)
    checkTrue(cor(featureValues(tmp_fld, value = "into")[, 1],
                  featureValues(tmp_fld, value = "into")[, 2],
                  use = "complete.obs") > 0.97)
    ## Check signal generation for already found peaks.
    for (i in 1:length(fileNames(tmp))) {
        fnd_pks <- chromPeaks(tmp)[chromPeaks(tmp)[, "sample"] == i, ]
        prm <- processHistory(tmp, type ="Peak detection")[[1]]@param
        ## Extract the data for these using the internal function.
        fld_pks <- xcms:::.getChromPeakData_matchedFilter(filterFile(tmp, i),
                                                          peakArea = fnd_pks,
                                                          sample_idx = i,
                                                          param = prm,
                                                          cn = colnames(fnd_pks))
        ## rt can not be the same, since for fillChromPeaks it is the rt of the
        ## maximum signal and for findChromPeaks it is the rt of the apex of the
        ## filtered/fitted peak.
        checkTrue(cor(fnd_pks[, "rt"], fld_pks[, "rt"]) > 0.99)
        ## mz: also not the same; in here we're weighting by sum of intensities
        ## per mz and in findChromPeaks by individual intensities.
        ## diffs <- fnd_pks[, "mz"] - fld_pks[, "mz"]
        ## checkTrue(max(diffs) < 1e-4)
        checkTrue(cor(fnd_pks[, "mz"], fld_pks[, "mz"]) > 0.99)
        ## into
        checkTrue(cor(fnd_pks[, "into"], fld_pks[, "into"]) > 0.99)
        ## checkEquals(fnd_pks[, "into"], fld_pks[, "into"])    
        ## maxo
        checkEquals(fnd_pks[, "maxo"], fld_pks[, "maxo"])    
    }

    ## Imputation.
    prm <- MatchedFilterParam(binSize = 0.2, impute = "lin")
    tmp <- findChromPeaks(faahko_od, param = prm)
    tmp <- groupChromPeaks(tmp, param = PeakDensityParam(sampleGroups = sg))
    tmp_fld <- fillChromPeaks(tmp)
    checkTrue(sum(is.na(featureValues(tmp_fld))) <
              sum(is.na(featureValues(tmp))))
    nas <- is.na(featureValues(tmp)[, 1]) | is.na(featureValues(tmp)[, 2])
    checkTrue(cor(featureValues(tmp, value = "into")[!nas, 1],
                  featureValues(tmp, value = "into")[!nas, 2]) > 0.9)
    ## Check signal generation for already found peaks.
    for (i in 1:length(fileNames(tmp))) {
        fnd_pks <- chromPeaks(tmp)[chromPeaks(tmp)[, "sample"] == i, ]
        prm <- processHistory(tmp, type ="Peak detection")[[1]]@param
        ## Extract the data for these using the internal function.
        fld_pks <- xcms:::.getChromPeakData_matchedFilter(filterFile(tmp, i),
                                                          peakArea = fnd_pks,
                                                          sample_idx = i,
                                                          param = prm,
                                                          cn = colnames(fnd_pks))
        ## rt can not be the same, since for fillChromPeaks it is the rt of the
        ## maximum signal and for findChromPeaks it is the rt of the apex of the
        ## filtered/fitted peak.
        checkTrue(cor(fnd_pks[, "rt"], fld_pks[, "rt"]) > 0.99)
        ## mz: also not the same; in here we're weighting by sum of intensities
        ## per mz and in findChromPeaks by individual intensities.
        ## diffs <- fnd_pks[, "mz"] - fld_pks[, "mz"]
        ## checkTrue(max(diffs) < 1e-4)
        checkTrue(cor(fnd_pks[, "mz"], fld_pks[, "mz"]) > 0.99)
        ## into
        checkTrue(cor(fnd_pks[, "into"], fld_pks[, "into"]) > 0.99)
        ## checkEquals(fnd_pks[, "into"], fld_pks[, "into"])    
        ## maxo
        checkTrue(cor(fnd_pks[, "maxo"], fld_pks[, "maxo"]) > 0.99)    
    }
    
    ## Own files.
    fls <- c("/Users/jo/data/2016/2016-11/NoSN/190516_POOL_N_POS_15.mzML",
             "/Users/jo/data/2016/2016-11/NoSN/190516_POOL_N_POS_19.mzML",
             "/Users/jo/data/2016/2016-11/NoSN/190516_POOL_N_POS_11.mzML")
    raw <- readMSData(fls, mode = "onDisk")
    pks <- findChromPeaks(raw, param = MatchedFilterParam(binSize = 0.05))
    sg <- rep(1, length(fileNames(raw)))
    pks <- groupChromPeaks(pks, param = PeakDensityParam(sampleGroups = sg))
    tmp_fld <- fillChromPeaks(pks)
    nas <- is.na(featureValues(pks)[, 1]) | is.na(featureValues(pks)[, 2])
    checkTrue(cor(featureValues(pks, value = "into")[!nas, 1],
                  featureValues(pks, value = "into")[!nas, 2]) > 0.97)
    checkTrue(cor(featureValues(tmp_fld, value = "into")[, 1],
                  featureValues(tmp_fld, value = "into")[, 2],
                  use = "complete.obs") > 0.97)
    ## Check signal generation for already found peaks.
    for (i in 1:length(fileNames(pks))) {
        fnd_pks <- chromPeaks(pks)[chromPeaks(pks)[, "sample"] == i, ]
        prm <- processHistory(pks, type ="Peak detection")[[1]]@param
        ## Extract the data for these using the internal function.
        fld_pks <- xcms:::.getChromPeakData_matchedFilter(filterFile(pks, i),
                                                          peakArea = fnd_pks,
                                                          sample_idx = i,
                                                          param = prm,
                                                          cn = colnames(fnd_pks))
        ## rt can not be the same, since for fillChromPeaks it is the rt of the
        ## maximum signal and for findChromPeaks it is the rt of the apex of the
        ## filtered/fitted peak.
        checkTrue(cor(fnd_pks[, "rt"], fld_pks[, "rt"]) > 0.99)
        ## mz: also not the same; in here we're weighting by sum of intensities
        ## per mz and in findChromPeaks by individual intensities.
        ## diffs <- fnd_pks[, "mz"] - fld_pks[, "mz"]
        ## checkTrue(max(diffs) < 1e-4)
        checkTrue(cor(fnd_pks[, "mz"], fld_pks[, "mz"]) > 0.99)
        ## into
        checkTrue(cor(fnd_pks[, "into"], fld_pks[, "into"]) > 0.99)
        checkEquals(fnd_pks[, "into"], fld_pks[, "into"])    
        ## maxo
        checkTrue(cor(fnd_pks[, "maxo"], fld_pks[, "maxo"]) > 0.99)    
        checkEquals(fnd_pks[, "maxo"], fld_pks[, "maxo"])    
    }
}

dontrun_exhaustive_fillChromPeaks_MSW <- function() {
    library(xcms)
    library(RUnit)
    library(msdata)
    fticrf <- list.files(system.file("fticr", package = "msdata"),
                         recursive = TRUE, full.names = TRUE)
    fticr <- readMSData(fticrf, msLevel. = 1, mode = "onDisk")
    p <- MSWParam(scales = c(1, 7), ampTh = 0.005,
                  SNR.method = "data.mean", winSize.noise = 500)
    fticr <- findChromPeaks(fticr, param = p)
    ## Now create the MzClustParam parameter object: we're assuming here that
    ## both samples are from the same sample group.
    p <- MzClustParam()
    fticr <- groupChromPeaks(fticr, param = p)
    res <- fillChromPeaks(fticr)
    ## Got a signal for all of em.
    checkTrue(!any(is.na(featureValues(res))))
    
    ## 1) Compare with what I get for xcmsSet.
    tmp_x <- xcmsSet(fticrf, method = "MSW", SNR.method = "data.mean",
                     winSize.noise = 500, scales = c(1, 7),
                     amp.Th = 0.005)
    sampclass(tmp_x) <- rep(1, length(sampnames(tmp_x)))
    tmp_x <- group(tmp_x, method = "mzClust")
    checkEquals(unname(groupval(tmp_x)), unname(featureValues(fticr)))
    checkEquals(unname(groupval(tmp_x, value = "maxo")),
                unname(featureValues(fticr, value = "maxo")))
    checkEquals(unname(groupval(tmp_x, value = "into")),
                unname(featureValues(fticr, value = "into")))
    checkEquals(unname(groupval(tmp_x, value = "mz")),
                unname(featureValues(fticr, value = "mz")))
    checkEquals(unname(groupval(tmp_x, value = "mzmin")),
                unname(featureValues(fticr, value = "mzmin")))
    checkEquals(unname(groupval(tmp_x, value = "mzmax")),
                unname(featureValues(fticr, value = "mzmax")))
    ## Fill peaks
    tmp_x <- fillPeaks(tmp_x, method = "MSW")
    checkTrue(!any(is.na(groupval(tmp_x))))
    checkEquals(unname(groupval(tmp_x)), unname(featureValues(res)))
    checkEquals(unname(groupval(tmp_x, value = "maxo")),
                unname(featureValues(res, value = "maxo")))
    ## This below could be made equal if we used the same approach to define
    ## which values to use. See comments on issue #130.
    ## checkEquals(unname(groupval(tmp_x, value = "into")),
    ##             unname(featureValues(res, value = "into")))
    ## plot(groupval(tmp_x, value = "into"), featureValues(res, value = "into"))
    checkTrue(cor(as.numeric(groupval(tmp_x, value = "into")),
                  as.numeric(featureValues(res, value = "into"))) > 0.999)
    ## plot(groupval(tmp_x, value = "mz"), featureValues(res, value = "mz"))
    checkTrue(cor(as.numeric(groupval(tmp_x, value = "mz")),
                  as.numeric(featureValues(res, value = "mz"))) > 0.999)
    checkEquals(unname(groupval(tmp_x, value = "mzmin")),
                unname(featureValues(res, value = "mzmin")))
    checkEquals(unname(groupval(tmp_x, value = "mzmax")),
                unname(featureValues(res, value = "mzmax")))

    ## Check if I could get what MSW gets.
    fdef <- featureDefinitions(fticr)
    pkArea <- do.call(
        rbind,
        lapply(
            fdef$peakidx, function(z) {
                tmp <- chromPeaks(fticr)[z, c("rtmin", "rtmax",
                                               "mzmin", "mzmax"),
                                          drop = FALSE]
                pa <- c(median(tmp[, 1]), median(tmp[, 2]),
                        median(tmp[, 3]), median(tmp[, 4]))
                return(pa)
            }
        ))
    colnames(pkArea) <- c("rtmin", "rtmax", "mzmin", "mzmax")
    pkArea <- cbind(group_idx = 1:nrow(pkArea), pkArea,
                    mzmed = fdef$mzmed)
    allPks <- xcms:::.getMSWPeakData(filterFile(fticr, file = 1),
                                     peakArea = pkArea,
                                     sample_idx = 1,
                                     cn = colnames(chromPeaks(fticr)))
    ## Get all chrom peaks part of a feature
    curP <- chromPeaks(res)[unique(unlist(featureDefinitions(res)$peakidx)), ]
    curP <- curP[curP[, "sample"] == 1, ]
    curP <- curP[order(curP[, "mz"]), ]
    ## checkEquals(allPks[, "mz"], curP[, "mz"])
    ## checkEquals(allPks[, "maxo"], curP[, "maxo"])

    ## ## INVESTIGATE FURTHER!
    fld <- curP[, "is_filled"] == 1
    plot(allPks[fld, "mz"], curP[fld, "mz"], main = "filled")
    abline(0, 1, col = "grey") ## same
    checkEquals(allPks[fld, "mz"], curP[fld, "mz"])
    plot(allPks[!fld, "mz"], curP[!fld, "mz"], main = "not filled")
    abline(0, 1, col = "grey") ## highly similar!!!    
    ## checkEquals(allPks[!fld, "mz"], curP[!fld, "mz"]) ## HIGHLY similar though
    ## into:
    plot(allPks[fld, "into"], curP[fld, "into"], main = "filled")
    checkEquals(allPks[fld, "into"], curP[fld, "into"])
    abline(0, 1, col = "grey") ## same
    plot(allPks[!fld, "into"], curP[!fld, "into"], main = "not filled")
    abline(0, 1, col = "grey") ## somewhat different!!!

    plot(allPks[fld, "maxo"], curP[fld, "maxo"], main = "filled")
    checkEquals(allPks[fld, "maxo"], curP[fld, "maxo"])
    abline(0, 1, col = "grey") ## same
    plot(allPks[!fld, "maxo"], curP[!fld, "maxo"], main = "not filled")
    abline(0, 1, col = "grey") ## somewhat different!!!
}


dontrun_exhaustive_fillChromPeaks_test <- function() {
    fls <- c("/Users/jo/data/2016/2016-11/NoSN/190516_POOL_N_POS_15.mzML",
             "/Users/jo/data/2016/2016-11/NoSN/190516_POOL_N_POS_19.mzML",
             "/Users/jo/data/2016/2016-11/NoSN/190516_POOL_N_POS_11.mzML")
    raw <- readMSData(fls, mode = "onDisk")
    pks <- findChromPeaks(raw, param = CentWaveParam(peakwidth = c(0.8, 20),
                                                     ppm = 40))
    sg <- rep(1, length(fileNames(raw)))
    pks_noRt <- groupChromPeaks(pks, param = PeakDensityParam(minFraction = 0.6,
                                                              sampleGroups = sg))
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
    gv <- featureValues(filled)
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
    still_missing <- is.na(rowSums(featureValues(filled)))
    ## Try using ppm:
    filled <- fillChromPeaks(pks_noRt, param = FillChromPeaksParam(ppm = 40))
    checkTrue(sum(still_missing) > sum(is.na(rowSums(featureValues(filled)))))
    filled <- fillChromPeaks(pks_noRt, param = FillChromPeaksParam(ppm = 40,
                                                                   expandMz = 2))
    checkTrue(sum(still_missing) > sum(is.na(rowSums(featureValues(filled)))))
    ## Check that the mz and rt are all within the mzmin-mzmax and rtmin-rtmax
    ## of the features.
    fts <- featureDefinitions(filled)
    for (i in 1:nrow(fts)) {
        pks <- chromPeaks(filled)[fts[i, "peakidx"][[1]], ]
        checkTrue(all(pks[, "mz"] >= fts[i, "mzmin"] &
                      pks[, "mz"] <= fts[i, "mzmax"]))
        checkTrue(all(pks[, "rt"] >= fts[i, "rtmin"] &
                      pks[, "rt"] <= fts[i, "rtmax"]))
    }
    checkTrue(all(chromPeaks(filled)[, "rt"] >= chromPeaks(filled)[, "rtmin"] &
                  chromPeaks(filled)[, "rt"] <= chromPeaks(filled)[, "rtmax"]))
    to_test <- chromPeaks(filled)[, "mz"] >= chromPeaks(filled)[, "mzmin"] &
        chromPeaks(filled)[, "mz"] <= chromPeaks(filled)[, "mzmax"]
    chromPeaks(filled)[!to_test, ]
    ## checkTrue(all(to_test))
    
    ## With adjusted retention times.
    pks <- adjustRtime(pks, param = ObiwarpParam())
    sg <- rep(1, length(fileNames(pks)))
    pks <- groupChromPeaks(pks, param = PeakDensityParam(minFraction = 0.6,
                                                         sampleGroups = sg))
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
    gv <- featureValues(filled)
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
    maxo <- rowMeans(featureValues(filled, value = "maxo"), na.rm = TRUE)
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
    raw_x <- readMSData(fl, mode = "onDisk")
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
    raw <- readMSData(fl, mode = "onDisk")
    tmp <- findChromPeaks(raw, param = CentWaveParam(peakwidth = c(1, 20)))
    pkInt3 <- xcms:::.getPeakInt3(tmp, chromPeaks(tmp))
    checkEquals(pkInt3, chromPeaks(tmp)[, "into"])
    
    ## own file matchedFilter peaks
    ## TODO: continue here - but be aware, for MatchedFilterParam we might have
    ## to do the stuff on the profile matrix!
    tmp <- findChromPeaks(raw_x, param = MatchedFilterParam())
}
