## Unsorted utility functions.
#' @include DataClasses.R

############################################################
## valueCount2ScanIndex
##
## @description Simple helper function that converts the number of values
## per scan/spectrum to an integer vector that can be passed to the base
## xcms functions/downstream C functions.
##
## @title Create index vector for internal C calls
## @param valCount Numeric vector representing the number of values per
## spectrum.
## @return An integer vector with the index (0-based) in the mz or intensity
## vectors indicating the start of a spectrum.
## @author Johannes Rainer
valueCount2ScanIndex <- function(valCount){
    ## Convert into 0 based.
    valCount <- cumsum(valCount)
    return(as.integer(c(0, valCount[-length(valCount)])))
}

############################################################
## useOriginalCode
##
## Simple function allowing the user to enable using the orignal
## code instead of the new implementations.
## This sets options.
##
##' @title Enable usage of old xcms code
##'
##' @description This function allows to enable the usage of old, partially
##' deprecated code from xcms by setting a corresponding global option. See
##' details for functions affected.
##'
##' @note Usage of old code is strongly dicouraged. This function is thought
##' to be used mainly in the transition phase from xcms to xcms version 3.
##' @details The functions/methods that will be affected by this are:
##' \itemize{
##' \item \code{\link{do_detectFeatures_matchedFilter}}
##' }
##' @param x Logical of length one to specify whether or not original
##' old code should be used in corresponding functions. If not provided the
##' function simply returns the value of the global option.
##' @return A logical of length one indicating whether old code is being
##' used.
##' @author Johannes Rainer
useOriginalCode <- function(x) {
    if (missing(x)) {
        res <- options()$BioC$xcms$useOriginalCode
        if (is.null(res))
            return(FALSE)
        return(res)
    }
    if (!is.logical(x))
        stop("'x' has to be logical.")
    b_opts <- getOption("BioC")
    x_opts <- b_opts$xcms
    x_opts$useOriginalCode <- x[1]
    b_opts$xcms <-  x_opts
    options(BioC = b_opts)
    return(options()$BioC$xcms$useOriginalCode)
}

## .getOriginalFunction <- function(x) {
##     if (!any(names(.ORIGINAL_FUNCTIONS)))
## }
## .ORIGINAL_FUNCTIONS <- c(
##     matchedFilter = ".matchedFilter_orig"
## )


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
