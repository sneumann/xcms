## Test findChromPeaks centWave

## library(faahKO)
fs <- c(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
        system.file('cdf/KO/ko16.CDF', package = "faahKO"),
        system.file('cdf/KO/ko18.CDF', package = "faahKO"),
        system.file('cdf/KO/ko19.CDF', package = "faahKO"))

## library(msdata)
## mzf <- c(system.file("microtofq/MM14.mzML", package = "msdata"),
##          system.file("microtofq/MM8.mzML", package = "msdata"))
## f <- msdata::proteomics(full.names = TRUE, pattern = "TMT_Erwinia")
xr <- deepCopy(faahko_xr_1)
onDisk <- filterFile(faahko_od, file = 1)

## Run peak detection in MS1 and ensure that the resulting object contains
## still all MS levels.
test_findChromPeaks_multiple_MS <- function() {
    msn_file <- system.file(package = "msdata", 
                            "proteomics/MS3TMT10_01022016_32917-33481.mzML.gz")
    msn_file <- system.file(
        package = "msdata",
        "proteomics/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzML.gz")
    msn_data <- readMSData(msn_file, mode = "onDisk")
    msn_xdata <- findChromPeaks(msn_data, param = CentWaveParam())
    checkEquals(msLevel(msn_data), msLevel(msn_xdata))
}

## Ensure that the reported peak integrated signal corresponds to the correct
## data. This is to ensure that issue #135 was fixed correctly.
test_findChromPeaks_centWave_peakIntensity <- function() {
    ## Reproduce with msdata files:
    fl <- system.file("microtofq/MM14.mzML", package = "msdata")
    raw <- readMSData(fl, mode = "onDisk")
    options(originalCentWave = TRUE)
    tmp <- findChromPeaks(raw, param = CentWaveParam(peakwidth = c(2, 10)))
    ## ## Use the getPeakInt2 which uses the rawMat function.
    ## pkI2 <- xcms:::.getPeakInt2(tmp, chromPeaks(tmp))
    ## ## Use the getPeakInt3 which uses the getEIC C function.
    ## pkI3 <- xcms:::.getPeakInt3(tmp, chromPeaks(tmp))
    ## ## These fail for the original centWave code.
    ## checkTrue(sum(pkI2 != chromPeaks(tmp)[, "into"]) > length(pkI2) / 2)
    ## ## checkEquals(unname(pkI2), unname(chromPeaks(tmp)[, "into"]))
    ## ## checkEquals(unname(pkI3), unname(chromPeaks(tmp)[, "into"]))
    ## checkEquals(pkI2, pkI3)
    ## Try with new implementation.
    options(originalCentWave = FALSE)
    tmp2 <- findChromPeaks(raw, param = CentWaveParam(peakwidth = c(2, 10)))
    ## Find different number of peaks:
    checkTrue(nrow(chromPeaks(tmp2)) != nrow(chromPeaks(tmp)))
    ## Are the peaks similar?
    id_1 <- paste(chromPeaks(tmp)[, "mz"], chromPeaks(tmp)[, "rt"])
    id_2 <- paste(chromPeaks(tmp2)[, "mz"], chromPeaks(tmp2)[, "rt"])
    ## But all of the ones from the old are ALSO in the new one.
    checkTrue(all(id_1 %in% id_2))
    ## Are the peaks the same?
    cp2 <- chromPeaks(tmp2)[id_2 %in% id_1, ]
    cn <- colnames(cp2)
    cn <- cn[!(cn %in% c("intb", "into", "rtmin", "rtmax"))]
    checkEquals(cp2[, cn], chromPeaks(tmp)[, cn])
    ## Are the values related?
    plot(cp2[, "into"], chromPeaks(tmp)[, "into"])   ## Very similar
    plot(cp2[, "intb"], chromPeaks(tmp)[, "intb"])   ## Very similar
    plot(cp2[, "rtmin"], chromPeaks(tmp)[, "rtmin"])   ## Very similar
    plot(cp2[, "rtmax"], chromPeaks(tmp)[, "rtmax"])   ## Very similar
    ## Use the getPeakInt3 which uses the getEIC C function.
    ## pkI2_2 <- xcms:::.getPeakInt2(tmp2, chromPeaks(tmp2))
    ## pkI3_2 <- xcms:::.getPeakInt3(tmp2, chromPeaks(tmp2))
    ## ## These fail for the original centWave code.
    ## checkEquals(unname(pkI2_2), unname(chromPeaks(tmp2)[, "into"]))
    ## checkEquals(unname(pkI3_2), unname(chromPeaks(tmp2)[, "into"]))
    ## checkEquals(pkI2_2, pkI3_2)

    
    ## The same for one of the test files; this works even with the original
    ## centWave code
    options(originalCentWave = TRUE)
    tmp <- filterFile(xod_xgrg, file = 3)
    ## ## Use the getPeakInt2 which uses the rawMat function.
    ## pkI2 <- xcms:::.getPeakInt2(tmp, chromPeaks(tmp))
    ## ## Use the getPeakInt3 which uses the getEIC C function.
    ## pkI3 <- xcms:::.getPeakInt3(tmp, chromPeaks(tmp))
    ## checkEquals(pkI2, pkI3)
    ## checkEquals(unname(pkI2), unname(chromPeaks(tmp)[, "into"]))
    ## checkEquals(unname(pkI3), unname(chromPeaks(tmp)[, "into"]))
    ## New modified centWave.
    options(originalCentWave = FALSE)
    tmp2 <- findChromPeaks(filterFile(faahko_od, file = 3),
                           CentWaveParam(noise = 10000, snthresh = 40))
    ## Even the identified peaks are identical!
    checkEquals(chromPeaks(tmp), chromPeaks(tmp2))
    ## Use the getPeakInt2 which uses the rawMat function.
    ## pkI2 <- xcms:::.getPeakInt2(tmp2, chromPeaks(tmp2))
    ## ## Use the getPeakInt3 which uses the getEIC C function.
    ## pkI3 <- xcms:::.getPeakInt3(tmp2, chromPeaks(tmp2))
    ## checkEquals(pkI2, pkI3)
    ## checkEquals(unname(pkI2), unname(chromPeaks(tmp2)[, "into"]))
    ## checkEquals(unname(pkI3), unname(chromPeaks(tmp2)[, "into"]))
    options(originalCentWave = TRUE)
}

## Exhaustive comparison of the original and the modified centWave. We don't
## expect everything to be the same, but the results should be ideally comparable
dontrun_exhaustive_original_new_centWave_comparison <- function() {
    ## faahKO test files.
    fl <- system.file("cdf/ko15.CDF", package = "msdata")
    raw <- readMSData(fl, mode = "onDisk")
    ## Default settings
    cwp <- CentWaveParam()
    options(originalCentWave = TRUE)
    orig <- findChromPeaks(raw, param = cwp)
    options(originalCentWave = FALSE)
    modi <- findChromPeaks(raw, param = cwp)
    ## Same number of peaks?
    checkTrue(nrow(chromPeaks(orig)) == nrow(chromPeaks(modi)))  ## YES
    ## Peaks the same?
    checkEquals(chromPeaks(orig), chromPeaks(modi))  ## YES
    ## modified settings:
    cwp <- CentWaveParam(ppm = 40, peakwidth = c(4, 40))
    options(originalCentWave = TRUE)
    orig <- findChromPeaks(raw, param = cwp)
    options(originalCentWave = FALSE)
    modi <- findChromPeaks(raw, param = cwp)
    ## Same number of peaks?
    checkTrue(nrow(chromPeaks(orig)) == nrow(chromPeaks(modi)))  ## YES
    ## Peaks the same?
    checkEquals(chromPeaks(orig), chromPeaks(modi))  ## YES
    ## Another round
    cwp <- CentWaveParam(ppm = 10, peakwidth = c(1, 10))
    options(originalCentWave = TRUE)
    orig <- findChromPeaks(raw, param = cwp)
    options(originalCentWave = FALSE)
    modi <- findChromPeaks(raw, param = cwp)
    ## Same number of peaks?
    checkTrue(nrow(chromPeaks(orig)) == nrow(chromPeaks(modi)))  ## YES
    ## Peaks the same?
    checkEquals(chromPeaks(orig), chromPeaks(modi))  ## YES
    
    ## msdata test files.
    fl <- system.file("microtofq/MM14.mzML", package = "msdata")
    raw <- readMSData(fl, mode = "onDisk")
    cwp <- CentWaveParam(ppm = 10, peakwidth = c(1, 10))
    options(originalCentWave = TRUE)
    orig <- findChromPeaks(raw, param = cwp)
    options(originalCentWave = FALSE)
    modi <- findChromPeaks(raw, param = cwp)
    ## Same number of peaks?
    checkTrue(nrow(chromPeaks(orig)) == nrow(chromPeaks(modi)))  ## YES
    ## Peaks the same? NO
    orig_id <- paste(chromPeaks(orig)[, "rt"], chromPeaks(orig)[, "mz"])
    modi_id <- paste(chromPeaks(modi)[, "rt"], chromPeaks(modi)[, "mz"])
    checkTrue(all(orig_id %in% modi_id))
    ## Check peaks and colnames.
    cn <- colnames(chromPeaks(orig))
    cn_diff <- c("into", "intb", "rtmin", "rtmax")
    cn <- cn[!(cn %in% cn_diff)]
    checkEquals(chromPeaks(orig)[, cn], chromPeaks(modi)[, cn])  ## YES
    ## Check those that are different.
    for (i in cn_diff) {
        plot(chromPeaks(orig)[, i], chromPeaks(modi)[, i])
        checkTrue(cor(chromPeaks(orig)[, i], chromPeaks(modi)[, i]) > 0.99)
    }
    ## Different settings.
    cwp <- CentWaveParam(ppm = 40, peakwidth = c(1, 20))
    options(originalCentWave = TRUE)
    orig <- findChromPeaks(raw, param = cwp)
    options(originalCentWave = FALSE)
    modi <- findChromPeaks(raw, param = cwp)
    ## Same number of peaks? NO
    checkTrue(nrow(chromPeaks(orig)) < nrow(chromPeaks(modi)))  ## YES
    ## Peaks the same? NO
    orig_id <- paste(chromPeaks(orig)[, "rt"], chromPeaks(orig)[, "mz"])
    modi_id <- paste(chromPeaks(modi)[, "rt"], chromPeaks(modi)[, "mz"])
    ## Are all from orig also in modi?
    checkTrue(all(orig_id %in% modi_id))  ## YES
    modi_pks <- chromPeaks(modi)[match(orig_id, modi_id), ]
    ## Compare the peaks.
    cn <- colnames(chromPeaks(orig))
    cn_diff <- c("into", "intb", "rtmin", "rtmax")
    cn <- cn[!(cn %in% cn_diff)]
    checkEquals(chromPeaks(orig)[, cn], modi_pks[, cn])  ## YES
    for (i in cn_diff) {
        plot(chromPeaks(orig)[, i], modi_pks[, i])
        checkTrue(cor(chromPeaks(orig)[, i], modi_pks[, i]) > 0.99)
    }
    ## Different settings.
    cwp <- CentWaveParam(ppm = 40, peakwidth = c(1, 10))
    options(originalCentWave = TRUE)
    orig <- findChromPeaks(raw, param = cwp)
    options(originalCentWave = FALSE)
    modi <- findChromPeaks(raw, param = cwp)
    ## Same number of peaks?
    checkTrue(nrow(chromPeaks(orig)) == nrow(chromPeaks(modi)))  ## YES
    ## Peaks the same?
    orig_id <- paste(chromPeaks(orig)[, "rt"], chromPeaks(orig)[, "mz"])
    modi_id <- paste(chromPeaks(modi)[, "rt"], chromPeaks(modi)[, "mz"])
    ## Are all from orig also in modi?
    checkTrue(all(orig_id %in% modi_id))  ## YES
    ## Compare the peaks.
    cn <- colnames(chromPeaks(orig))
    cn_diff <- c("into", "intb", "rtmin", "rtmax")
    cn <- cn[!(cn %in% cn_diff)]
    checkEquals(chromPeaks(orig)[, cn], chromPeaks(modi)[, cn])  ## YES
    for (i in cn_diff) {
        plot(chromPeaks(orig)[, i], chromPeaks(modi)[, i])
        checkTrue(cor(chromPeaks(orig)[, i], chromPeaks(modi)[, i]) > 0.99)
    }
    
    ## Other file
    fl <- system.file("microtofq/MM8.mzML", package = "msdata")
    raw <- readMSData(fl, mode = "onDisk")
    cwp <- CentWaveParam(ppm = 30, peakwidth = c(1, 10))
    options(originalCentWave = TRUE)
    orig <- findChromPeaks(raw, param = cwp)
    options(originalCentWave = FALSE)
    modi <- findChromPeaks(raw, param = cwp)
    ## Same number of peaks?
    checkTrue(nrow(chromPeaks(orig)) == nrow(chromPeaks(modi)))  ## YES
    ## Peaks the same? NO, but all in orig are in modi
    orig_id <- paste(chromPeaks(orig)[, "rt"], chromPeaks(orig)[, "mz"])
    modi_id <- paste(chromPeaks(modi)[, "rt"], chromPeaks(modi)[, "mz"])
    checkTrue(all(orig_id %in% modi_id))
    ## Check peaks and colnames.
    cn <- colnames(chromPeaks(orig))
    cn_diff <- c("into", "intb", "rtmin", "rtmax")
    cn <- cn[!(cn %in% cn_diff)]
    checkEquals(chromPeaks(orig)[, cn], chromPeaks(modi)[, cn])  ## YES
    ## Check those that are different.
    for (i in cn_diff) {
        plot(chromPeaks(orig)[, i], chromPeaks(modi)[, i])
        checkTrue(cor(chromPeaks(orig)[, i], chromPeaks(modi)[, i]) > 0.99)
    }
    ## Different settings.
    cwp <- CentWaveParam(ppm = 40, peakwidth = c(3, 30))
    options(originalCentWave = TRUE)
    orig <- findChromPeaks(raw, param = cwp)
    options(originalCentWave = FALSE)
    modi <- findChromPeaks(raw, param = cwp)
    ## Same number of peaks? NO
    checkTrue(nrow(chromPeaks(orig)) < nrow(chromPeaks(modi)))
    ## Peaks the same? NO, but all in orig are in modi
    orig_id <- paste(chromPeaks(orig)[, "rt"], chromPeaks(orig)[, "mz"])
    modi_id <- paste(chromPeaks(modi)[, "rt"], chromPeaks(modi)[, "mz"])
    checkTrue(all(orig_id %in% modi_id))
    ## Check peaks and colnames.
    cmn <- match(orig_id, modi_id)
    cn <- colnames(chromPeaks(orig))
    cn_diff <- c("into", "intb", "rtmin", "rtmax")
    cn <- cn[!(cn %in% cn_diff)]
    checkEquals(chromPeaks(orig)[, cn], chromPeaks(modi)[cmn, cn])  ## YES
    ## Check those that are different.
    for (i in cn_diff) {
        plot(chromPeaks(orig)[, i], chromPeaks(modi)[cmn, i])
        checkTrue(cor(chromPeaks(orig)[, i], chromPeaks(modi)[cmn, i]) > 0.99)
    }
    
    ## own files.
    fl <- "/Users/jo/data/2016/2016-11/NoSN/190516_POOL_N_POS_12.mzML"
    raw <- readMSData(fl, mode = "onDisk")
    ## Default settings
    cwp <- CentWaveParam()
    options(originalCentWave = TRUE)
    orig <- findChromPeaks(raw, param = cwp)
    options(originalCentWave = FALSE)
    modi <- findChromPeaks(raw, param = cwp)
    ## Same number of peaks? NO
    checkTrue(nrow(chromPeaks(orig)) < nrow(chromPeaks(modi)))
    ## Peaks the same? NO, but all in orig are in modi
    orig_id <- paste(chromPeaks(orig)[, "rt"], chromPeaks(orig)[, "mz"])
    modi_id <- paste(chromPeaks(modi)[, "rt"], chromPeaks(modi)[, "mz"])
    checkTrue(all(orig_id %in% modi_id))
    ## Check peaks and colnames.
    cmn <- match(orig_id, modi_id)
    cn <- colnames(chromPeaks(orig))
    cn_diff <- c("into", "intb", "rtmin", "rtmax", "maxo")
    cn <- cn[!(cn %in% cn_diff)]
    checkEquals(chromPeaks(orig)[, cn], chromPeaks(modi)[cmn, cn])  ## YES
    ## Check those that are different.
    for (i in cn_diff) {
        plot(chromPeaks(orig)[, i], chromPeaks(modi)[cmn, i])
        checkTrue(cor(chromPeaks(orig)[, i], chromPeaks(modi)[cmn, i]) > 0.99)
    }
    ## Different settings
    cwp <- CentWaveParam(ppm = 30, peakwidth = c(1, 10))
    options(originalCentWave = TRUE)
    orig <- findChromPeaks(raw, param = cwp)
    options(originalCentWave = FALSE)
    modi <- findChromPeaks(raw, param = cwp)
    ## Same number of peaks? NO, more with modified versions
    checkTrue(nrow(chromPeaks(orig)) < nrow(chromPeaks(modi)))
    ## Peaks the same? NO, but all in orig are in modi
    orig_id <- paste(chromPeaks(orig)[, "rt"], chromPeaks(orig)[, "mz"])
    modi_id <- paste(chromPeaks(modi)[, "rt"], chromPeaks(modi)[, "mz"])
    checkTrue(all(orig_id %in% modi_id))
    ## Check peaks and colnames.
    cmn <- match(orig_id, modi_id)
    cn <- colnames(chromPeaks(orig))
    cn_diff <- c("into", "intb", "rtmin", "rtmax", "maxo")
    cn <- cn[!(cn %in% cn_diff)]
    checkEquals(chromPeaks(orig)[, cn], chromPeaks(modi)[cmn, cn])  ## YES
    ## Check those that are different.
    for (i in cn_diff) {
        plot(chromPeaks(orig)[, i], chromPeaks(modi)[cmn, i])
        checkTrue(cor(chromPeaks(orig)[, i], chromPeaks(modi)[cmn, i]) > 0.98)
    }
    ## Different settings
    cwp <- CentWaveParam(ppm = 40, peakwidth = c(1, 60))
    options(originalCentWave = TRUE)
    orig <- findChromPeaks(raw, param = cwp)
    options(originalCentWave = FALSE)
    modi <- findChromPeaks(raw, param = cwp)
    ## Same number of peaks? NO, more with modified versions
    checkTrue(nrow(chromPeaks(orig)) < nrow(chromPeaks(modi)))
    ## Peaks the same? NO, but all in orig are in modi
    orig_id <- paste(chromPeaks(orig)[, "rt"], chromPeaks(orig)[, "mz"])
    modi_id <- paste(chromPeaks(modi)[, "rt"], chromPeaks(modi)[, "mz"])
    checkTrue(all(orig_id %in% modi_id))
    ## Check peaks and colnames.
    cmn <- match(orig_id, modi_id)
    cn <- colnames(chromPeaks(orig))
    cn_diff <- c("into", "intb", "rtmin", "rtmax", "maxo")
    cn <- cn[!(cn %in% cn_diff)]
    checkEquals(chromPeaks(orig)[, cn], chromPeaks(modi)[cmn, cn])  ## YES
    ## Check those that are different.
    for (i in cn_diff) {
        plot(chromPeaks(orig)[, i], chromPeaks(modi)[cmn, i])
        checkTrue(cor(chromPeaks(orig)[, i], chromPeaks(modi)[cmn, i]) > 0.98)
    }

    ## Is there something common to the peaks found only by the modified version?
    common_pks <- chromPeaks(modi)[cmn, ]
    unique_pks <- chromPeaks(modi)[-cmn, ]
    cn <- colnames(common_pks)
    ## mz
    boxplot(list(common = common_pks[, "mz"], unique = unique_pks[, "mz"]),
            varwidth = TRUE, main = "mz")
    ## OK, average mz of unique peaks is smaller.

    ## mzrange
    boxplot(list(common = common_pks[, "mzmax"] - common_pks[, "mzmin"],
                 unique = unique_pks[, "mzmax"] - unique_pks[, "mzmin"]),
            varwidth = TRUE, main = "mz range")
    ## Seems to be smaller too.

    ## rt
    boxplot(list(common = common_pks[, "rt"], unique = unique_pks[, "rt"]),
            varwidth = TRUE, main = "rt")
    ## hm, rt larger in unique

    ## rtrange
    boxplot(list(common = log2(common_pks[, "rtmax"] - common_pks[, "rtmin"]),
                 unique = log2(unique_pks[, "rtmax"] - unique_pks[, "rtmin"])),
            varwidth = TRUE, main = "rt range")
    ## rtrange is same.

    ## into
    boxplot(list(common = log2(common_pks[, "into"]),
                 unique = log2(unique_pks[, "into"])),
            varwidth = TRUE, main = "into")
    ## same.

    ## sn
    boxplot(list(common = log2(common_pks[, "sn"]),
                 unique = log2(unique_pks[, "sn"])),
            varwidth = TRUE, main = "sn")
    ## sn is higher in unique

    ## Check chromatogram for some:
    chrPlot <- function(i, raw) {
        rtr <- common_pks[i, c("rtmin", "rtmax")]
        rtr[1] <- rtr[1] - 2
        rtr[2] <- rtr[2] + 2
        chr_cmn <- chromatogram(raw, rt = rtr,
                                mz = common_pks[i, c("mzmin", "mzmax")])
        rtr <- unique_pks[i, c("rtmin", "rtmax")]
        rtr[1] <- rtr[1] - 2
        rtr[2] <- rtr[2] + 2
        chr_unq <- chromatogram(raw, rt = rtr,
                                mz = unique_pks[i, c("mzmin", "mzmax")])
        par(mfrow = c(1, 2))
        plot(rtime(chr_cmn[[1]]), intensity(chr_cmn[1, 1]), main = "common peak",
             type = "l")
        abline(v = common_pks[i, c("rtmin", "rtmax")], col = "grey")
        plot(rtime(chr_unq[[1]]), intensity(chr_unq[1, 1]), main = "unique peak",
             type = "l")
        abline(v = unique_pks[i, c("rtmin", "rtmax")], col = "grey")
    }
    chrPlot(1, raw = raw)

    i <- sample(1:nrow(unique_pks), 10)
    chrPlot(i[1], raw = raw)
    chrPlot(i[2], raw = raw)
    chrPlot(i[3], raw = raw)
    chrPlot(i[4], raw = raw)
    chrPlot(i[5], raw = raw)
    chrPlot(i[6], raw = raw)
    chrPlot(i[7], raw = raw)
    chrPlot(i[8], raw = raw)
    chrPlot(i[9], raw = raw)
    chrPlot(i[10], raw = raw)

    ## Summary:
    ## Same peaks are identified by both methods, but the modified centWave finds
    ## eventually more peaks.
    ## The into, intb, rtmin and rtmax differ for some peaks found by both
    ## methods but values are highly correlated (R > 0.99).
}


## Compare what's the difference between the original centWave and the new one
## fixing the peak integration problem (issue #135)
dontrun_compare_orig_new_centWave <- function() {
    mzVals <- xr@env$mz
    intVals <- xr@env$intensity
    ## Define the values per spectrum:
    valsPerSpect <- diff(c(xr@scanindex, length(mzVals)))
    res_orig <- xcms:::.centWave_orig(mz = mzVals,
                                      int = intVals,
                                      scantime = xr@scantime,
                                      valsPerSpect,
                                      noise = 4000, verboseColumns = TRUE)
    res_new <- xcms:::.centWave_new(mz = mzVals,
                                    int = intVals,
                                    scantime = xr@scantime,
                                    valsPerSpect,
                                    noise = 4000, verboseColumns = TRUE)
    checkEquals(res_orig, res_new) ## That should always work.
    
    ## Use the previously defined ones as ROIs - this should simulate the
    ## addIsoROIs
    ## Slightly increase the mz like in addIsoROIs.
    rL <- as.data.frame(res_orig)
    ## Extend the mzmin and mzmax if needed.
    tittle <- rL[, "mz"] * (25 / 2) / 1E6
    expand_mz <- (rL[, "mzmax"] - rL[, "mzmin"]) < (tittle * 2)
    if (any(expand_mz)) {
        rL[expand_mz, "mzmin"] <- rL[expand_mz, "mz"] -
            tittle[expand_mz]
        rL[expand_mz, "mzmax"] <- rL[expand_mz, "mz"] + tittle[expand_mz]
    }
    rL <- split(rL, f = 1:nrow(rL))
    res_orig2 <- xcms:::.centWave_orig(mz = mzVals,
                                       int = intVals,
                                       scantime = xr@scantime,
                                       valsPerSpect,
                                       noise = 4000, verboseColumns = TRUE,
                                       roiList = rL,
                                       firstBaselineCheck = FALSE)
    res_new2 <- xcms:::.centWave_new(mz = mzVals,
                                     int = intVals,
                                     scantime = xr@scantime,
                                     valsPerSpect,
                                     noise = 4000, verboseColumns = TRUE,
                                     roiList = rL,
                                     firstBaselineCheck = FALSE)
    ## I get more peaks with the modified version:
    checkTrue(nrow(res_orig2) < nrow(res_new2))
    ## Sort them by rt and mz
    idx <- order(res_orig2[, "rt"], res_orig2[, "mz"])
    res_orig2 <- res_orig2[idx, ]
    idx <- order(res_new2[, "rt"], res_new2[, "mz"])
    res_new2 <- res_new2[idx, ]
    ## Are identified peaks similar?
    id_1 <- paste(res_orig2[, "rt"], res_orig2[, "mz"])
    id_2 <- paste(res_new2[, "rt"], res_new2[, "mz"])
    ## Are all from res_orig2 in res_new2?
    checkTrue(length(id_1) == sum(id_1 %in% id_2)) ## YES
    ## Are the values for these the same?
    same_new2 <- res_new2[match(id_1, id_2), ]
    ## check columns:
    cn <- colnames(same_new2)
    cn <- cn[!(cn %in% c("mzmin"))]
    for (i in cn) {
        checkEquals(res_orig2[, i], same_new2[, i])
    }
    ## So, everything is identical, EXCEPT mzmin:
    idx <- which(same_new2[, "mzmin"] != res_orig2[, "mzmin"])
    res_orig2[idx, "mzmin"]
    same_new2[idx, "mzmin"]
    ## Are the different ones 0?
    checkTrue(all(res_orig2[idx, "mzmin"] == 0))
}

test_do_findChromPeaks_centWave <- function() {
    ## xr <- xcmsRaw(fs[1], profstep = 0)
    ## We expect that changing a parameter has an influence on the result.
    mzVals <- xr@env$mz
    intVals <- xr@env$intensity
    ## Define the values per spectrum:
    valsPerSpect <- diff(c(xr@scanindex, length(mzVals)))
    res1 <- do_findChromPeaks_centWave(mz = mzVals,
                                       int = intVals,
                                       scantime = xr@scantime,
                                       valsPerSpect,
                                       snthresh = 200,
                                       noise = 4000)
    ## Eventually disable the sleep option to improve speed!
    res2 <- do_findChromPeaks_centWave(mz = mzVals,
                                       int = intVals,
                                       scantime = xr@scantime,
                                       valsPerSpect,
                                       snthresh = 500,
                                       noise = 4000, sleep = 0.01)
    checkTrue(nrow(res1) > nrow(res2))

    ## Check scanrange on findPeaks.centWave.
    res_1 <- findPeaks.centWave(xr, scanrange = c(90, 345), noise = 2000)
    xr <- xr[90:345]
    mzVals <- xr@env$mz
    intVals <- xr@env$intensity
    ## Define the values per spectrum:
    valsPerSpect <- diff(c(xr@scanindex, length(mzVals)))
    res_2 <- do_findChromPeaks_centWave(mz = mzVals, int = intVals,
                                        scantime = xr@scantime, valsPerSpect,
                                        noise = 2000)
    checkEquals(res_1@.Data, res_2)
}

## Evaluate the peak detection method using the centWave method on
## OnDiskMSnExp and on MSnExp objects.
test_findChromPeaks_centWave <- function() {
    ## Control
    library(MSnbase)
    ## xr <- xcmsRaw(fs[1], profstep = 0)
    ppm <- 40
    snthresh <- 40
    res_x <- findPeaks.centWave(xr, ppm = ppm, snthresh = snthresh,
                                noise = 100000)@.Data
    ## Bypass xcmsRaw
    xs <- xcmsSet(fs[1], profparam = list(profstep = 0), ppm = ppm,
                  snthresh = snthresh, method = "centWave",
                  noise = 100000)
    checkEquals(xs@peaks[, colnames(res_x)], res_x)
    ## OnDiskMSnExp
    ## onDisk <- readMSData(fs[1], msLevel. = 1, mode = "onDisk")
    cwp <- CentWaveParam(ppm = ppm, snthresh = snthresh, noise = 100000)
    res <- findChromPeaks(onDisk, param = cwp, return.type = "list")
    checkEquals(res[[1]], peaks(xs)@.Data)

    checkException(findChromPeaks(onDisk, param = cwp, msLevel = 2))
    ## ## MSnExp
    ## inMem <- readMSData(f[1], msLevel. = 1)
    ## suppressWarnings(
    ##     res_2 <- findChromPeaks(inMem, param = cwp, return.type = "list")
    ## )
    ## checkEquals(res_2[[1]], peaks(xs)@.Data)

    ## returning an xcmsSet
    res <- findChromPeaks(onDisk, param = cwp, return.type = "xcmsSet")
    checkEquals(peaks(res), peaks(xs))
    ## suppressWarnings(
    ##     res <- findChromPeaks(inMem, param = cwp, return.type = "xcmsSet")
    ## )
    ## checkEquals(peaks(res), peaks(xs))

    ## Return type XCMSnExp
    res <- findChromPeaks(onDisk, param = cwp)
    checkTrue(hasChromPeaks(res))
    checkTrue(!hasAdjustedRtime(res))
    checkTrue(!hasFeatures(res))
    checkEquals(peaks(xs)@.Data, chromPeaks(res)[, -ncol(chromPeaks(res))])
}

dontrun_test_benchmark_centWaves <- function() {
    library(msdata)
    f <- msdata::proteomics(full.names = TRUE, pattern = "TMT_Erwinia")
    library(microbenchmark)
    library(MSnbase)
    library(xcms)
    ##
    ## xr <- xcmsRaw(f[1], profstep = 0)
    ppm <- 40
    snthresh <- 40

    cwp <- CentWaveParam(ppm = ppm, snthresh = snthresh)
    ## onDisk <- readMSData(f[1], msLevel. = 1, mode = "onDisk")
    register(SerialParam())
    system.time(
        tmp <- findChromPeaks(onDisk, param = cwp)
    ) ## 9.7sec
    system.time(
        tmp <- findChromPeaks(onDisk, param = cwp, return.type = "xcmsSet")
    ) ## 12sec
    system.time(
        tmp <- xcmsSet(f[1], profparam = list(profstep = 0), ppm = ppm,
                       snthresh = snthresh, method = "centWave")
    ) ## 11.99sec

    inMem <- readMSData(f[1], msLevel. = 1)
    register(SerialParam())

    ## findChromPeaks,MSnExp and findPeaks.centWave should be about similar.
    microbenchmark(findPeaks.centWave(xr, ppm = ppm, snthresh = snthresh),
                   findChromPeaks(inMem, param = cwp), times = 3)
    ## findPeaks.centWave is about 1 second faster.

    ## findChromPeaks,OnDiskMSnExp and xcmsSet should be about similar.
    microbenchmark(xcmsSet(f[1], profparam = list(profstep = 0), ppm = ppm,
                           snthresh = snthresh, method = "centWave"),
                   findChromPeaks(onDisk, param = cwp),
                   findChromPeaks(inMem, param = cwp),
                   times = 3)
}

dontrun_test_benchmark_centWaves <- function() {
    library(msdata)
    f <- msdata::proteomics(full.names = TRUE, pattern = "TMT_Erwinia")
    library(microbenchmark)
    library(MSnbase)
    library(xcms)
    ##
    ## xr <- xcmsRaw(f[1], profstep = 0)
    ppm <- 40
    snthresh <- 40

    cwp <- CentWaveParam(ppm = ppm, snthresh = snthresh)
    ## onDisk <- readMSData(f[1], msLevel. = 1, mode = "onDisk")
    register(SerialParam())
    system.time(
        tmp <- findChromPeaks(onDisk, param = cwp)
    ) ## 9.7sec
    system.time(
        tmp <- findChromPeaks(onDisk, param = cwp, return.type = "xcmsSet")
    ) ## 12sec
    system.time(
        tmp <- xcmsSet(f[1], profparam = list(profstep = 0), ppm = ppm,
                       snthresh = snthresh, method = "centWave")
    ) ## 11.99sec

    inMem <- readMSData(f[1], msLevel. = 1)
    register(SerialParam())

    ## findChromPeaks,MSnExp and findPeaks.centWave should be about similar.
    microbenchmark(findPeaks.centWave(xr, ppm = ppm, snthresh = snthresh),
                   findChromPeaks(inMem, param = cwp), times = 3)
    ## findPeaks.centWave is about 1 second faster.

    ## findChromPeaks,OnDiskMSnExp and xcmsSet should be about similar.
    microbenchmark(xcmsSet(f[1], profparam = list(profstep = 0), ppm = ppm,
                           snthresh = snthresh, method = "centWave"),
                   findChromPeaks(onDisk, param = cwp),
                   findChromPeaks(inMem, param = cwp),
                   times = 3)
}




############################################################
## This is only relevant during development of the do_ function
## to evaluate that results are identical.
dontrun_test_do_findChromPeaks_centWave_impl <- function() {

    for (i in 1:length(fs)) {
        ppm = 25
        peakwidth = c(20, 50)
        snthresh = 10
        prefilter = c(3, 100)
        mzCenterFun = "wMean"
        integrate = 1
        mzdiff = -0.001
        fitgauss = FALSE
        noise = 0
        verboseColumns = FALSE

        xr <- xcmsRaw(fs[i])

        ## Default settings
        .runAndCompare(xr, ppm, peakwidth, snthresh, prefilter, mzCenterFun,
                       integrate, mzdiff, fitgauss, noise, verboseColumns)
        ## xcms: 14.6 sec
        ## do_ : 13 sec

        ppm <- 10
        .runAndCompare(xr, ppm, peakwidth, snthresh, prefilter, mzCenterFun,
                   integrate, mzdiff, fitgauss, noise, verboseColumns)
        ## xcms: 15 sec
        ## do_ : 13.3 sec

        peakwidth <- c(3, 30)
        .runAndCompare(xr, ppm, peakwidth, snthresh, prefilter, mzCenterFun,
                       integrate, mzdiff, fitgauss, noise, verboseColumns)
        ## xcms: 11.4 sec
        ## do_ :  9.5 sec

        snthresh <- 15
        .runAndCompare(xr, ppm, peakwidth, snthresh, prefilter, mzCenterFun,
                       integrate, mzdiff, fitgauss, noise, verboseColumns)
        ## xcms: 10.6 sec
        ## do_ :  8.8 sec

        fitgauss <- TRUE
        .runAndCompare(xr, ppm, peakwidth, snthresh, prefilter, mzCenterFun,
                       integrate, mzdiff, fitgauss, noise, verboseColumns)
        ## xcms: 12.5 sec
        ## do_ : 10.7 sec

        verboseColumns <- TRUE
        .runAndCompare(xr, ppm, peakwidth, snthresh, prefilter, mzCenterFun,
                       integrate, mzdiff, fitgauss, noise, verboseColumns)
        ## xcms: 12.2 sec
        ## do_ : 10.6 sec
    }
}

## That's to compare the functions in version 1.49.7.
.runAndCompare <- function(xr, ppm, peakwidth, snthresh, prefilter, mzCenterFun,
                           integrate, mzdiff, fitgauss, noise, verboseColumns) {
    require(RUnit)
    mz <- xr@env$mz
    int <- xr@env$intensity
    scantime <- xr@scantime
    scanindex <- xr@scanindex
    a <- system.time(
        ## That's the method called inside do_...
        xrDo <- xcms:::.centWave_orig(mz = mz, int = int, scantime = scantime,
                                      valsPerSpect = diff(c(scanindex, length(mz))),
                                      ppm = ppm, peakwidth = peakwidth,
                                      snthresh = snthresh,
                                      prefilter = prefilter,
                                      mzCenterFun = mzCenterFun,
                                      integrate = integrate,
                                      mzdiff = mzdiff,
                                      fitgauss = fitgauss,
                                      noise = noise,
                                      verboseColumns = verboseColumns)
    ) ## 12.7
    ## Run the original centWave code on xcmsRaw:
    b <- system.time(
        xrPeaks <- xcms:::.findPeaks.centWave_orig(xr,
                                                   ppm = ppm,
                                                   peakwidth = peakwidth,
                                                   snthresh = snthresh,
                                                   prefilter = prefilter,
                                                   mzCenterFun = mzCenterFun,
                                                   integrate = integrate,
                                                   mzdiff = mzdiff,
                                                   fitgauss = fitgauss,
                                                   noise = noise,
                                                   verbose.columns = verboseColumns)
    )  ## 15.4
    ## Compare.
    cat("DO: ", a, "\n")
    cat("XCMS: ", b, "\n")
    if (!checkEquals(new("xcmsPeaks", xrDo), xrPeaks))
        stop("do_ and xcms yield different results!")
}


## Some speed tests.
.otherTest <- function() {
    Testv <- c(2, 4.2, 34.1, 34.5, 6.4, 6.3, 1.2)
    RforM <- matrix(nrow = 0, ncol = length(Testv))
    system.time(
        for(i in 1:5000){
            RforM <- rbind(RforM, Testv)
        }
    ) ## 1.27
    ## with append to list.
    RforL <- vector("list", 0)
    system.time(
        for(i in 1:5000){
            RforL <- c(RforL, Testv)
        }
    ) ## 1.12
    system.time(
        RapplyL <- lapply(1:5000, function(z) {return(Testv)})
    ) ## 0.003
    RM <- matrix(nrow=5000, ncol = length(Testv))
    system.time(
        for (i in 1:5000) {
            RM[i, ] <- Testv
        }
    ) ## 0.006

    ## Compare adding to list instead of adding to existing. [[]]
    RexL <- vector("list", 5000)
    system.time(
        for (i in 1:5000){
            RexL[[i]] <- Testv
        }
    ) ## 0.005
    ## Dynamically...
    RexL <- list()
    system.time(
        for (i in 1:5000){
            RexL[[i]] <- Testv
        }
    ) ## 0.005

}
