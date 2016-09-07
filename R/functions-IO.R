############################################################
## IO functions.
#' @include functions-utils.R

############################################################
## isCdfFile
##
## Just guessing whether the file is a CDF file based on its ending.
isCdfFile <- function(x) {
    fileEnds <- c("cdf", "nc")
    ## check for endings and and ending followed by a . (e.g. cdf.gz)
    patts <- paste0("\\.", fileEnds, "($|\\.)")
    res <- sapply(patts, function(z) {
        grep(z, x, ignore.case = TRUE)
    })
    return(any(unlist(res)))
}

############################################################
## isMzMLFile
##
## Just guessing whether the file is a mzML file based on its ending.
isMzMLFile <- function(x) {
    fileEnds <- c("mzxml", "mzml", "mzdata")
    ## check for endings and and ending followed by a . (e.g. mzML.gz)
    patts <- paste0("\\.", fileEnds, "($|\\.)")
    res <- sapply(patts, function(z) {
        grep(z, x, ignore.case = TRUE)
    })
    return(any(unlist(res)))
}

############################################################
## readRawData
##
## Function to read the raw data from an cdf, mzml file. Might eventually
## replace the loadRaw methods.
## returns list with rt, tic, scanindex, mz and intensity
readRawData <- function(x, includeMSn = FALSE) {
    def_backend <- "Ramp"  ## Eventually use pwiz...
    header_cols <- c("retentionTime", "acquisitionNum", "totIonCurrent")
    if (isCdfFile(x)) {
        backend <- "netCDF"
    } else {
        if (isMzMLFile(x)) {
            backend <- def_backend
            header_cols <- c(header_cols, "polarity")
        } else {
            stop("Unknown file type.")
        }
    }
    msd <- mzR::openMSfile(x, backend = backend)
    on.exit(mzR::close(msd))
    ## That's due to issue https://github.com/lgatto/MSnbase/issues/151
    on.exit(gc(), add = TRUE)
    hdr <- mzR::header(msd)
    idx_ms1 <- which(hdr$msLevel == 1)
    if (length(idx_ms1) == 0)
        stop("No MS1 data found in file ", x, "!")
    pks <- mzR::peaks(msd, idx_ms1)
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

    if (includeMSn) {
        if (backend == "netCDF") {
            warning("Reading of MSn spectra for NetCDF not supported.")
        } else {
            ## Read MSn data
            idx_ms2 <- which(hdr$msLevel > 1)
            if (length(idx_ms2) > 0) {
                hdr_ms2 <- hdr[idx_ms2, , drop = FALSE]
                pks <- mzR::peaks(msd, idx_ms2)
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
                warning("MSn spectra requested but not found.")
            }
        }
    }
    resList
}
