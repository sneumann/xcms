############################################################
## IO functions.
#' @include functions-utils.R

############################################################
## readRawData
##
##' Function to read the raw data from an cdf, mzml file. Might eventually
##' replace the loadRaw methods. The function returns list with rt, tic,
##' scanindex, mz and intensity. For more details on the motivation to implement
##' this function see issue #65 on github.
##' @param x The file name.
##' @param includeMSn logical(1) indicating whether MS level > 1 should be loaded
##' too. Only supported for mzML files.
##' @param dropEmptyScans Scans/spectra without peaks are not returned if
##' \code{dropEmptyScans = TRUE}. If \code{FALSE} all spectra from the input
##' file are returned. This is to be consistent with the code before
##' xcms version 1.51.1 (see issue #67
##' https://github.com/sneumann/xcms/issues/67).
##'
##' @param backend \code{character} allowing to manually specify the mzR
##'     backend. If \code{NULL}, it uses the automatic backend determination
##'     from mzR.
##' 
##' @return A \code{list} with rt, tic, scanindex, mz and intensity.
##' 
##' @noRd
readRawData <- function(x, includeMSn = FALSE, dropEmptyScans = TRUE,
                        backend = NULL) {
    ## def_backend <- "Ramp"  ## Eventually use pwiz...
    header_cols <- c("retentionTime", "acquisitionNum", "totIonCurrent")
    msd <- mzR::openMSfile(x, backend = backend)
    on.exit(mzR::close(msd))
    ## That's due to issue https://github.com/lgatto/MSnbase/issues/151
    on.exit(rm(msd), add = TRUE)
    on.exit(gc(), add = TRUE)
    hdr <- mzR::header(msd)
    if (any(colnames(hdr) == "polarity"))
        header_cols <- c(header_cols, "polarity")
    idx_ms1 <- which(hdr$msLevel == 1)
    ## Drop empty spectra; see https://github.com/sneumann/xcms/issues/67
    if (dropEmptyScans & length(idx_ms1) > 0) {
        idx_ms1 <- idx_ms1[hdr[idx_ms1, "peaksCount"] > 0]
    }
    ## Fix issue #174 in RMassBank.
    if (length(idx_ms1) == 0 & !includeMSn)
        stop("No MS1 data found in file ", x, "!")
    if (length(idx_ms1)) {
        pks <- mzR::peaks(msd, idx_ms1)
        ## Fix problem with single spectrum files (issue #66)
        if (is(pks, "matrix"))
            pks <- list(pks)
        valsPerSpect <- lengths(pks) / 2
        pks <- do.call(rbind, pks)
        hdr_ms1 <- hdr[idx_ms1, header_cols,
                       drop = FALSE]
        resList <- list(rt = hdr_ms1$retentionTime,
                        acquisitionNum = hdr_ms1$acquisitionNum,
                        tic = hdr_ms1$totIonCurrent,
                        scanindex = valueCount2ScanIndex(valsPerSpect),
                        mz = pks[, 1],
                        intensity = pks[, 2],
                        polarity = hdr_ms1$polarity)
    } else {
        warning("No MS1 spectra available in file ", basename(x))
        resList <- list(rt = numeric(),
                        acquisitionNum = integer(),
                        tic = numeric(),
                        scanindex = integer(),
                        mz = numeric(),
                        intensity = numeric(),
                        polarity = numeric())
    }
    if (includeMSn) {
        ## Read MSn data
        idx_ms2 <- which(hdr$msLevel > 1)
        ## Drop empty spectra; see https://github.com/sneumann/xcms/issues/67
        if (dropEmptyScans & length(idx_ms2) > 0) {
            idx_ms2 <- idx_ms2[hdr[idx_ms2, "peaksCount"] > 0]
        }
        if (length(idx_ms2) > 0) {
            hdr_ms2 <- hdr[idx_ms2, , drop = FALSE]
            pks <- mzR::peaks(msd, idx_ms2)
            if (is(pks, "matrix"))
                pks <- list(pks)
            valsPerSpect <- lengths(pks) / 2
            pks <- do.call(rbind, pks)
            resList$MSn <- list(rt = hdr_ms2$retentionTime,
                                acquisitionNum = hdr_ms2$acquisitionNum,
                                precursorNum = hdr_ms2$precursorScanNum,
                                precursorMZ = hdr_ms2$precursorMZ,
                                precursorIntensity = hdr_ms2$precursorIntensity,
                                peaksCount = hdr_ms2$peaksCount,
                                msLevel = hdr_ms2$msLevel,
                                precursorCharge = hdr_ms2$precursorCharge,
                                scanindex = valueCount2ScanIndex(valsPerSpect),
                                collisionEnergy = hdr_ms2$collisionEnergy,
                                mz = pks[, 1],
                                intensity = pks[, 2])
        } else {
            warning("MSn spectra requested but none present in the file.")
        }
    }
    resList
}
