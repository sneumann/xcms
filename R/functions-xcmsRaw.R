## All functions for xcmsRaw methods should be put here.
#' @include DataClasses.R c.R

xcmsRaw <- function(filename, profstep = 1, profmethod = "bin",
                    profparam = list(),
                    includeMSn = FALSE, mslevel = NULL,
                    scanrange = NULL) {

    object <- new("xcmsRaw")
    object@env <- new.env(parent=.GlobalEnv)

    ## Change between old and new code; new code uses the "newer" mzR approach
    ## to read the file.
    ## if (useOriginalCode()) {
    object@filepath <- xcmsSource(filename)
    rawdata <- loadRaw(object@filepath, includeMSn = includeMSn)
    ## } else {
    ##     object@filepath <- new("xcmsFileSource", filename)
    ##     rawdata <- readRawData(filename, includeMSn = includeMSn)
    ## }

    rtdiff <- diff(rawdata$rt)
    if (any(rtdiff == 0))
        warning("There are identical scantimes.")

    if (any(rtdiff < 0)) {
        badtimes <- which(rtdiff < 0)
        stop(paste("Time for scan ", badtimes[1], " (",
                   rawdata$rt[[badtimes[1]]], ") greater than scan ",
                   badtimes[1]+1, " (", rawdata$rt[[badtimes[1]+1]], ")",
                   sep = ""))
    }

    object@scantime <- rawdata$rt
    object@tic <- rawdata$tic
    object@scanindex <- rawdata$scanindex
    object@env$mz <- rawdata$mz
    object@env$intensity <- rawdata$intensity

    ## Doing first an eventual scanrange subsetting so that we don't have to
    ## re-calculate the profile matrix later.
    if (length(scanrange) < 2) {
        scanrange <- c(1, length(object@scantime))
    } else {
        scanrange <- range(scanrange)
    }
    if (min(scanrange) < 1 | max(scanrange) > length(object@scantime)) {
        scanrange[1] <- max(1, scanrange[1])
        scanrange[2] <- min(length(object@scantime), scanrange[2])
        message("Provided scanrange was adjusted to ", scanrange[1]," - ", scanrange[2])
    }
    if (!is.null(rawdata$acquisitionNum)) {
        ## defined only for mzData and mzXML
        object@acquisitionNum <- rawdata$acquisitionNum
    }
    if (!is.null(rawdata$polarity)) {
        object@polarity <- factor(rawdata$polarity,
                                  levels = c(0, 1, -1),
                                  labels = c("negative", "positive", "unknown"))
    }
    ##
    ## After the MS1 data, take care of MSn
    ##
    if(!is.null(rawdata$MSn) ) {
        object@env$msnMz <- rawdata$MSn$mz
        object@env$msnIntensity <- rawdata$MSn$intensity
        object@msnScanindex <- rawdata$MSn$scanindex
        object@msnAcquisitionNum <- rawdata$MSn$acquisitionNum
        object@msnLevel <- rawdata$MSn$msLevel
        object@msnRt <- rawdata$MSn$rt
        object@msnPrecursorScan <- match(rawdata$MSn$precursorNum,
                                         object@acquisitionNum)
        object@msnPrecursorMz <- rawdata$MSn$precursorMZ
        object@msnPrecursorIntensity <- rawdata$MSn$precursorIntensity
        object@msnPrecursorCharge <- rawdata$MSn$precursorCharge
        object@msnCollisionEnergy <- rawdata$MSn$collisionEnergy
    }
    ## setting the scanrange.
    scanrange(object) <- as.numeric(scanrange)

    ## Subset by scanrange
    object <- object[scanrange[1]:scanrange[2]]

    mslevel(object) <- as.numeric(mslevel)
    object@mzrange <- range(object@env$mz, na.rm = TRUE)

    object@profmethod <- profmethod
    object@profparam <- profparam
    ## Creating profile matrix if profstep > 0
    if (profstep)
        profStep(object) <- profstep

    ## if (!is.null(scanrange)) {
    ##     ## Scanrange filtering
    ##     keepidx <- seq.int(1, length(object@scantime)) %in%
    ##         seq.int(scanrange[1], scanrange[2])
    ##     object <- split(object, f=keepidx)[["TRUE"]]
    ## }

    if (!missing(mslevel) & !is.null(mslevel)) {
        ## Issue #101
        ## Copy the ms2 level data only if mslevel > 1
        if (max(mslevel) > 1) {
            object <- msn2ms(object)
            object <- split(object, f=object@msnLevel==mslevel)$"TRUE"
        }
        ## fix xcmsRaw metadata, or always calculate later than here ?
    }
    return(object)
}

############################################################
## specNoise
specNoise <- function(spec, gap = quantile(diff(spec[,"mz"]), .9)) {

    ## In a spectrum with just one raw peak we can't calculate noise
    if (nrow(spec) < 2) {
        return(0)
    }

    intmean <- mean(spec[,"intensity"])

    mzlen <- diff(range(spec[,"mz"]))
    mzdiff <- diff(spec[,"mz"])
    gaplen <- sum(mzdiff[mzdiff > gap])

    weighted.mean(c(intmean, min(spec[,"intensity"])/2), c(1 - gaplen/mzlen,
                                                           gaplen/mzlen))
}

############################################################
## specPeaks
specPeaks <- function(spec, sn = 20, mzgap = .2) {

    noise <- specNoise(spec)

    spectab <- matrix(nrow = 0, ncol = 3)
    colnames(spectab) <- c("mz", "intensity", "fwhm")

    while (spec[i <- which.max(spec[,"intensity"]), "intensity"] > noise*sn) {

        mz <- spec[i,"mz"]
        intensity <- spec[i,"intensity"]
        fwhmrange <- descendValue(spec[,"intensity"], spec[i,"intensity"]/2, i)

        if (fwhmrange[1] > 1 && fwhmrange[2] < nrow(spec)) {
            fwhm1 <- spec[fwhmrange[1],"mz"] - (spec[fwhmrange[1],"intensity"]-intensity/2)*diff(spec[fwhmrange[1]-1:0,"mz"])/diff(spec[fwhmrange[1]-1:0,"intensity"])
            fwhm2 <- spec[fwhmrange[2],"mz"] - (spec[fwhmrange[2],"intensity"]-intensity/2)*diff(spec[fwhmrange[2]+1:0,"mz"])/diff(spec[fwhmrange[2]+1:0,"intensity"])

            fwhm <- fwhm2-fwhm1

            if (!any(abs(spectab[,"mz"] - mz) <= mzgap))
                spectab <- rbind(spectab, c(mz, intensity, fwhm))
        }

        peakrange <- descendValue(spec[,"intensity"], min(noise*sn, spec[i,"intensity"]/4), i)
        spec[seq(peakrange[1], peakrange[2]),"intensity"] <- 0
    }

    spectab
}


## issue #74
##' @title Extract an EIC from the profile matrix
##'
##' @description The \code{profEIC} does extract the EIC not from the raw data,
##' but from the profile matrix. To get the EIC from the raw data use the
##' \code{\link{rawEIC}} method. The \code{profEIC} is a replacement of the
##' old \code{getEIC} method implementation (both of functions \code{getEICold}
##' and \code{getEICnew}) supporting the same input arguments and returning the
##' same result object, but with more sanity checks and using the newer binning
##' and interpolation functions.
##'
##' @note This method uses the new binning and linear interpolation functionality
##' i.e. the \code{\link{binYonX}} and \code{\link{imputeLinInterpol}}. In
##' contrast to the old \code{getEIC} implementation (pre xcms 1.51.1), this
##' method performs also considerably more input parameter validations.
##' @noRd
profEIC <- function(object, mzrange, rtrange = NULL, step = 0.1) {
    ## Input argument checking:
    if (missing(mzrange)) {
        mzrange <- matrix(object@mzrange, nrow = 1)
    } else {
        if (length(mzrange) == 2)
            mzrange <- matrix(as.numeric(mzrange), nrow = 1)
        if (!is.matrix(mzrange))
            stop("'mzrange' is supposed to be a two-column matrix with each row",
                 " containing the min and max value specifying the mz range.")
    }
    if (all(c("mzmin", "mzmax") %in% colnames(mzrange)))
        mzrange <- mzrange[, c("mzmin", "mzmax"), drop = FALSE]
    ## rtrange
    if (is.null(rtrange)) {
        rtrange <- matrix(range(object@scantime), nrow = 1)
    } else {
        if (length(rtrange) == 2)
            rtrange <- matrix(as.numeric(rtrange), nrow = 1)
        if (!is.matrix(rtrange))
            stop("'rtrange' is supposed to be a two-column matrix with each row",
                 " containing the min and max value specifying the rt range.")
    }
    ## Check sizes...
    if (nrow(rtrange) != nrow(mzrange)) {
        ## Try to fix:
        if (nrow(rtrange) == 1)
            rtrange <- matrix(rep(rtrange[1, ], each = nrow(mzrange)),
                              nrow = nrow(mzrange))
        if (nrow(mzrange) == 1)
            mzrange <- matrix(rep(mzrange[1, ], each = nrow(rtrange)),
                              nrow = nrow(rtrange))
        if (nrow(rtrange) != nrow(mzrange))
            stop("'rtrange' and 'mzrange' have to have the same number of rows!")
    }
    if (ncol(rtrange) != 2 | ncol(mzrange) != 2)
        stop("Number of columns of 'rtrange' and 'mzrange' have to be 2!")
    ## Profile generation settings:
    profmat <- NULL
    pi <- profinfo(object)
    method <- pi$method
    baselevel <- pi$baselevel
    if (is.null(baselevel)) {
        baseValue <- min(object@env$intensity, na.rm = TRUE) / 2
    } else {
        baseValue <- baselevel
    }
    basespace <- pi$basespace
    impMeths <- c("none", "lin", "linbase", "intlin")
    names(impMeths) <- c("bin", "binlin", "binlinbase", "intlin")
    impute <- impMeths[method]
    if (impute == "intlin") {
        ## Have to calculate the full profile matrix...
        suppressMessages(
            profmat <- profMat(object, method = "intlin", step = step)
        )
    }
    ## Re-use the existing profile matrix?
    if (length(object@env$profile) > 0) {
        ## If the settings are the same:
        if (profStep(object) == step)
            profmat <- object@env$profile
    }
    if (is.null(profmat)) {
        valsPerSpect <- diff(c(object@scanindex, length(object@env$mz)))
        toIdx <- cumsum(valsPerSpect)
        fromIdx <- c(1L, toIdx[-length(toIdx)] + 1L)
    } else {
        mass <- seq(floor(min(object@env$mz) / step) * step,
                    ceiling(max(object@env$mz) / step) * step,
                    by = step)
    }
    ## Get the global ranges; individual ranges have to be within these!
    object_rtrange <- range(object@scantime)
    object_mzrange <- object@mzrange
    ## Initialize object
    eic <- vector("list", nrow(rtrange))
    for (i in 1:nrow(rtrange)) {
        rtr <- range(rtrange[i, ])
        mzr <- range(mzrange[i, ])
        ## Check the parameters!
        if (!(rtr[2] <= object_rtrange[2] & rtr[1] >= object_rtrange[1]))
            stop("'rtrange' number ", i, " (", paste(rtr, collapse = ", "), ") ",
                 "is outside of the retention time range of 'object'")
        if (!(mzr[2] <= object_mzrange[2] & mzr[1] >= object_mzrange[1]))
            stop("'mzrange' number ", i, " (", paste(mzr, collapse = ", "), ") ",
                 "is outside of the mz value range of 'object'")
        if (!is.null(profmat)) {
            ## Re-use the existing profile matrix to calculate.
            imz <- findRange(mass, c(mzr[1]-.5*step, mzr[2]+0.5*step), TRUE)
            irt <- which(object@scantime >= rtr[1]
                         & object@scantime <= rtr[2])
            if (length(imz) == 0)
                stop("Specified mz range ", mzr, " outside of the measured",
                     " mz range!")
            if (length(irt) == 0)
                stop("Specified retention time range ", rtr, " outside of ",
                     "the measured retention time range!")
            eic[[i]] <- cbind(rt = object@scantime[irt],
                              intensity = colMax(profmat[imz[1]:imz[2], irt,
                                                         drop=FALSE]))
        } else {
            ## 1) Determine which spectra to consider
            inSpectra <- which(object@scantime >= rtr[1]
                               & object@scantime <= rtr[2])
            if (length(inSpectra) == 0)
                stop("Specified retention time range ", rtr, " outside of ",
                     "the measured retention time range!")
            mass <- seq(floor(mzr[1]/step)*step,
                        ceiling(mzr[2]/step)*step, by = step)
            nBins <- length(mass)
            binRes <- binYonX(x = object@env$mz,
                              y = object@env$intensity,
                              nBins = nBins,
                              binFromX = mass[1],
                              binToX = mass[nBins],
                              fromIdx = fromIdx[inSpectra],
                              toIdx = toIdx[inSpectra],
                              baseValue = ifelse(impute == "none",
                                                 yes = 0, no = NA),
                              shiftByHalfBinSize = TRUE,
                              sortedX = TRUE)
            if (length(binRes) == 1)
                binRes <- list(binRes)
            bin_size <- binRes[[1]]$x[2] - binRes[[1]]$x[1]
            if (is.null(basespace)) {
                distance <- floor(0.075 / bin_size)
            } else {
                distance <- floor(basespace[1] / bin_size)
            }
            binVals <- lapply(binRes, function(z) {
                return(max(imputeLinInterpol(z$y, method = impute,
                                             distance = distance,
                                             noInterpolAtEnds = TRUE,
                                             baseValue = baseValue),
                           na.rm = TRUE))
            })
            eic[[i]] <- cbind(rt = object@scantime[inSpectra],
                              intensity = unlist(binVals, use.names = FALSE))
        }
    }
    new("xcmsEIC", eic = list(xcmsRaw=eic), mzrange = mzrange, rtrange = rtrange,
        rt = "raw", groupnames = character(0))
}

profEIC2 <- function(object, mzrange, rtrange = NULL, step = 0.1) {
    ## This version does guess the indices to be passed directly to binYonX.
    ## binFromX binToX represents the mz range.
    ## fromIdx toIdx can be used to specify the

}


############################################################
## getEICOld
## that's the original getEIC version.
## We've got a problem if step = 0! (relates to issue #39)
getEICOld <- function(object, mzrange, rtrange = NULL, step = 0.1) {
    ## if mzrange and rtrange is not provided use the full range.
    if(missing(mzrange)){
        mzrange <- matrix(object@mzrange, nrow=1)
        colnames(mzrange) <- c("mzmin", "mzmax")
    }
    if(is.null(rtrange)){
        rtrange <- matrix(range(object@scantime), nrow=1)
        colnames(rtrange) <- c("rtmin", "rtmax")
    }
    profFun <- match.profFun(object)
    if (all(c("mzmin","mzmax") %in% colnames(mzrange)))
        mzrange <- mzrange[,c("mzmin", "mzmax"),drop=FALSE]

### Create EIC buffer
    mrange <- range(mzrange)
    mass <- seq(floor(mrange[1]/step)*step, ceiling(mrange[2]/step)*step, by = step)
    bufsize <- min(100, length(mass))
    buf <- profFun(object@env$mz, object@env$intensity, object@scanindex,
                   bufsize, mass[1], mass[bufsize], TRUE, object@profparam)
    bufidx <- integer(length(mass))
    idxrange <- c(1, bufsize)
    bufidx[idxrange[1]:idxrange[2]] <- 1:bufsize

    if (missing(rtrange))
        eic <- matrix(nrow = nrow(mzrange), ncol = ncol(buf))
    else
        eic <- vector("list", nrow(rtrange))

    for (i in order(mzrange[,1])) {
        imz <- findRange(mass, c(mzrange[i,1]-.5*step, mzrange[i,2]+.5*step), TRUE)
### Update EIC buffer if necessary
        if (bufidx[imz[2]] == 0) {
            bufidx[idxrange[1]:idxrange[2]] <- 0
            idxrange <- c(max(1, min(imz[1], length(mass)-bufsize+1)), min(bufsize+imz[1]-1, length(mass)))
            bufidx[idxrange[1]:idxrange[2]] <- 1:(diff(idxrange)+1)
            buf <- profFun(object@env$mz, object@env$intensity, object@scanindex,
                           diff(idxrange)+1, mass[idxrange[1]], mass[idxrange[2]],
                           TRUE, object@profparam)
        }
        if (missing(rtrange)) {
            eic[i,] <- colMax(buf[bufidx[imz[1]:imz[2]],,drop=FALSE])
        } else {
            eic[[i]] <- matrix(c(object@scantime, colMax(buf[bufidx[imz[1]:imz[2]],,drop=FALSE])),
                               ncol = 2)[object@scantime >= rtrange[i,1] & object@scantime <= rtrange[i,2],,drop=FALSE]
            colnames(eic[[i]]) <- c("rt", "intensity")
        }
    }

    invisible(new("xcmsEIC", eic = list(xcmsRaw=eic), mzrange = mzrange, rtrange = rtrange,
                  rt = "raw", groupnames = character(0)))

}

############################################################
## getEICNew
## what's different in this method?
## 1) we're not (re-) calculating the profile matrix if it already exists and if the step argument
##    is the same.
## 2) by not using the buffer with the fixed (max) size of 100 we're no longer limited to small m/z
##    ranges, thus we can use the method to extract the EIC for the full m/z range (i.e. the base
##    peak chromatogram BPC).
## 3) the method might be slower.
## We've got a problem if step = 0! (relates to issue #39)
getEICNew <- function(object, mzrange, rtrange = NULL,
                      step = 0.1, BPPARAM = bpparam()) {
    ## if mzrange and rtrange is not provided use the full range.
    if(missing(mzrange)){
        if(length(object@mzrange) == 2){
            mzrange <- matrix(object@mzrange, nrow=1)
        }else{
            mzrange <- matrix(c(min(object@env$mz), max(object@env$mz)), nrow=1)
        }
        colnames(mzrange) <- c("mzmin", "mzmax")
    }
    if(is.null(rtrange)){
        rtrange <- matrix(range(object@scantime), nrow=1)
        colnames(rtrange) <- c("rtmin", "rtmax")
    }
    ## rtrange and mzrange have to have the same number of rows!
    if(nrow(rtrange)!=nrow(mzrange)){
        stop("rtrange and mzrange have to have the same number of rows!")
    }
    profFun <- match.profFun(object)
    if (all(c("mzmin","mzmax") %in% colnames(mzrange)))
        mzrange <- mzrange[,c("mzmin", "mzmax"),drop=FALSE]

    ## check if we have the profile and if, if the profile step fits the step...
    if(any(names(object@env) == "profile" )){
        pStep <- profStep(object)
        if (length(pStep) == 0)
            pStep <- step
        if(pStep != step){
            ## delete that profile matrix since the step differs.
            rm(list="profile", envir=object@env)
        }
    }

    mass <- seq(floor(min(object@env$mz)/step)*step,
                ceiling(max(object@env$mz)/step)*step, by = step)
    ## check if we've got already the profile matrix available, if yes, we don't have to
    ## re-calculate anything.
    if(!any(names(object@env) == "profile")){
        ## calculate the profile matrix.
        object@env$profile <- profFun(object@env$mz, object@env$intensity,
                                      object@scanindex, length(mass), mass[1],
                                      mass[length(mass)], TRUE, object@profparam)
    }

    ## once we've got the full profile matrix we go on and extract the EICs.
    parms <- vector("list", length=nrow(rtrange))
    for(i in 1:length(parms)){
        parms[[i]] <- list( mzrange=mzrange[i, ], rtrange=rtrange[i, ] )
    }
    ## check if we could run the code on multiple cpus...
    eic <- bplapply(parms, FUN=function(z){
                      imz <- findRange(mass, c(z$mzrange[1]-.5*step, z$mzrange[2]+0.5*step), TRUE)
                      irt <- which(object@scantime >= z$rtrange[1] & object@scantime <= z$rtrange[2])
                      e <- matrix(c(object@scantime[irt],
                                    colMax(object@env$profile[imz[1]:imz[2], irt, drop=FALSE])), ncol=2)
                      colnames(e) <- c("rt", "intensity")
                      return(e)
                  }, BPPARAM=BPPARAM)

    invisible(new("xcmsEIC", eic = list(xcmsRaw=eic), mzrange = mzrange, rtrange = rtrange,
                  rt = "raw", groupnames = character(0)))
}

############################################################
## split.xcmsRaw
split.xcmsRaw <- function(x, f, drop = TRUE, ...)
{
    if (length(x@msnLevel)>0)
        warning ("MSn information will be dropped")

    if (!is.factor(f))
        f <- factor(f)

    scanidx <- unclass(f)

    lcsets <- vector("list", length(levels(f)))
    names(lcsets) <- levels(f)

    for (i in unique(scanidx)) {
        lcsets[[i]] <- x

        lcsets[[i]]@env <- new.env(parent=.GlobalEnv)

        lcsets[[i]]@tic = x@tic[scanidx == i]
        lcsets[[i]]@scantime = x@scantime[scanidx == i]
        lcsets[[i]]@polarity = x@polarity[scanidx == i]
        lcsets[[i]]@acquisitionNum = x@acquisitionNum[scanidx == i]
        lcsets[[i]]@mzrange = x@mzrange[scanidx == i]

        startindex = x@scanindex[which(scanidx == i)]+1

        endindex = x@scanindex[which(scanidx == i) +1]
        endindex[which(is.na(endindex))] <- length(x@env$mz)

        if (length(endindex) > 1) {

            scanlength <- endindex-startindex+1

            lcsets[[i]]@scanindex <- as.integer(c(0, cumsum(scanlength[1:length(scanlength)-1])))
            ptidx <- unlist(sequences(cbind(startindex, endindex)))
        } else {
            ## Single Scan
            ptidx <- 0:endindex
            lcsets[[i]]@scanindex <- as.integer(0)
        }

        lcsets[[i]]@env$mz <- x@env$mz[ptidx]
        lcsets[[i]]@env$intensity <- x@env$intensity[ptidx]

        profStep(lcsets[[i]]) <- profStep(x)
    }

    if (drop)
        lcsets <- lcsets[seq(along = lcsets) %in% scanidx]

    lcsets
}

############################################################
## sequences
sequences <- function(seqs) {
    apply(seqs, 1, FUN=function(x) {x[1]:x[2]})
}

############################################################
## match.profFun
match.profFun <- function(object) {
    match.fun(.profFunctions[[profMethod(object)]])
}

############################################################
## remakeTIC
remakeTIC<-function(object){
    for(i in 1:length(object@scanindex)){
        object@tic[i]<-sum(getScan(object, i)[,"intensity"])
    }
    return(object)
}

############################################################
## getPeaks
#' @description Replacement function for the original getPeaks method/function
#'     that does no longer use the deprecated \code{profFun} functions. This
#'     function uses the \code{binYonX} and \code{imputeLinInterpol} to perform
#'     the binning (and missing value imputation).
#'
#' @param object An \code{xcmsRaw} object.
#'
#' @param peakrange \code{matrix} with 4 required columns \code{"mzmin"},
#'     \code{"mzmax"}, \code{"rtmin"} and \code{"rtmax"}.
#'
#' @param step \code{numeric(1)} defining the bin size for the profile matrix
#'     generation.
#'
#' @author Johannes Rainer
#'
#' @noRd
.getPeaks_new <- function(object, peakrange, step = 0.1) {
    ## Here we're avoiding the profFun call.
    if (all(c("mzmin","mzmax","rtmin","rtmax") %in% colnames(peakrange)))
        peakrange <- peakrange[,c("mzmin","mzmax","rtmin","rtmax"),drop=FALSE]
    stime <- object@scantime

    pi <- profinfo(object)
    method <- pi$method
    if (missing(step))
        step <- pi$step
    if (step == 0)
        step <- 0.1
    baselevel <- pi$baselevel
    basespace <- pi$basespace
    vps <- diff(c(object@scanindex, length(object@env$mz)))
    
    cat("method: ", method, "\n")
    cat("step: ", step, "\n")
    ## Create the profile matrix:
    pMat <- .createProfileMatrix(mz = object@env$mz, int = object@env$intensity,
                                 valsPerSpect = vps,
                                 method = method,
                                 step = step,
                                 baselevel = baselevel,
                                 basespace = basespace,
                                 returnBreaks = TRUE,
                                 baseValue = 0,
                                 mzrange. = NULL)
    brks <- pMat$breaks
    pMat <- pMat$profMat  ## rows are masses, cols are retention times/scans.
    bin_size <- diff(brks[1:2])
    bin_half <- bin_size / 2
    ## Calculate the mean mass per bin using the breaks used for the binning.
    ## Note: these define the real mass breaks as they have been used for the
    ## binning. Simply using seq(floor...) as in the original code is wrong
    ## because the mass bins are calculated wrongly. The bin size is != step,
    ## bin size is marginally smaller and, for larger mz the correct mass
    ## bin will be wrongly identified.
    mass <- brks[-length(brks)] + bin_half ## midpoint for the breaks
    mass_range <- range(mass)

    ## Prepare the result matrix.
    cnames <- c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "into", "maxo")
    rmat <- matrix(nrow = nrow(peakrange), ncol = length(cnames))
    colnames(rmat) <- cnames

    for (i in order(peakrange[, 1])) {
        imz <- findRange(mass, c(peakrange[i, 1] - bin_half,
                                 peakrange[i, 2] + bin_half), TRUE)
        iret <- findRange(stime, peakrange[i, 3:4], TRUE)
        idx_imz <- imz[1]:imz[2]
        idx_iret <- iret[1]:iret[2]
        ## Extract the intensity matrix for the mz-rt range: rows are mz, cols
        ## rt values.
        ymat <- pMat[idx_imz, idx_iret, drop = FALSE]
        ## Define the maximum intensity, is one value per mz.
        ymax <- colMax(ymat)
        iymax <- which.max(ymax)

        ## The width in rt.
        pwid <- diff(stime[iret])/diff(iret)

        ## Calculate sum across rt. For each mz we get one value.
        rosm <- rowSums(ymat)
        limz <- length(idx_imz)
        if (length(rosm) != limz) { ## that happens for some reason
            warning("weighted.mean  : x and w must have the same length \n")
            rosm <- rep(1, limz)  ## fallback to mean
        }
        ## mean mz:
        rmat[i, 1] <- weighted.mean(mass[idx_imz], rosm) ## mz; its not the
        ## position of the largest intensity!
        if (is.nan(rmat[i,1]) || is.na(rmat[i,1])) ##  R2.11 :  weighted.mean()
            ## results in NA (not NaN) for zero weights
            rmat[i, 1] <- mean(peakrange[i, 1:2])

        rmat[i, 2:3] <- peakrange[i, 1:2]            ## mzmin, mzmax
        rmat[i, 4] <- stime[idx_iret][iymax] ## rt
        rmat[i, 5:6] <- peakrange[i, 3:4]            ## rtmin, rtmax

        if (peakrange[i, 3] <  stime[1] ||
            peakrange[i, 4] > stime[length(stime)] ||
            is.nan(pwid)) {
            warning("getPeaks: Peak  m/z:", peakrange[i, 1], "-",
                    peakrange[i, 2], ",  RT:", peakrange[i, 3], "-",
                    peakrange[i, 4], "is out of retention time range for ",
                    "this sample (", object@filepath,
                    "), using zero intensity value.\n")
            rmat[i, 7:8] <- 0
        } else {
            rmat[i, 7] <- pwid * sum(ymax)  ## into
            rmat[i, 8] <- ymax[iymax]       ## maxo
        }
    }
    invisible(rmat)
}

#' @description Original getPeaks function. This should be removed at some point
#'     as it uses deprecated API.
#' @noRd
.getPeaks_orig <- function(object, peakrange, step = 0.1) {
    profFun <- match.profFun(object)
    if (all(c("mzmin","mzmax","rtmin","rtmax") %in% colnames(peakrange)))
        peakrange <- peakrange[,c("mzmin","mzmax","rtmin","rtmax"),drop=FALSE]
    stime <- object@scantime

### Create EIC buffer
    ## This is NOT calculated for the full file.
    mrange <- range(peakrange[,1:2])
    ## These mass bins are slightly different from the ones that are used
    ## by the binning function, since within the binning function the step/bin
    ## size is recalculated!
    mass <- seq(floor(mrange[1]/step)*step, ceiling(mrange[2]/step)*step, by = step)
    bufsize <- min(100, length(mass))
    buf <- profFun(object@env$mz, object@env$intensity, object@scanindex,
                   bufsize, mass[1], mass[bufsize], TRUE, object@profparam)
    bufidx <- integer(length(mass))
    idxrange <- c(1, bufsize)
    bufidx[idxrange[1]:idxrange[2]] <- 1:bufsize

    cnames <- c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "into", "maxo")
    rmat <- matrix(nrow = nrow(peakrange), ncol = length(cnames))
    colnames(rmat) <- cnames

    for (i in order(peakrange[,1])) {
        imz <- findRange(mass, c(peakrange[i,1]-.5*step, peakrange[i,2]+.5*step), TRUE)
        iret <- findRange(stime, peakrange[i,3:4], TRUE)

### Update EIC buffer if necessary
        if (bufidx[imz[2]] == 0) {
            bufidx[idxrange[1]:idxrange[2]] <- 0
            idxrange <- c(max(1, imz[1]), min(bufsize+imz[1]-1, length(mass)))
            bufidx[idxrange[1]:idxrange[2]] <- 1:(diff(idxrange)+1)
            buf <- profFun(object@env$mz, object@env$intensity, object@scanindex,
                           diff(idxrange)+1, mass[idxrange[1]], mass[idxrange[2]],
                           TRUE, object@profparam)
        }
        ## Extract the intensity matrix for the mz-rt range: rows are mz, cols
        ## rt values.
        ymat <- buf[bufidx[imz[1]:imz[2]],iret[1]:iret[2],drop=FALSE]
        ## Define the maximum intensity, is one value per mz.
        ymax <- colMax(ymat)
        iymax <- which.max(ymax)

        ## The width in rt.
        pwid <- diff(stime[iret])/diff(iret)

        ## Calculate sum across rt. For each mz we get one value.
        rosm <- rowSums(ymat)
        limz <- length(imz[1]:imz[2])
        if (length(rosm) != limz) { ## that happens for some reason
            warning("weighted.mean  : x and w must have the same length \n")
            rosm <- rep(1, limz)  ## fallback to mean
        }
        ## mean mz:
        rmat[i,1] <- weighted.mean(mass[imz[1]:imz[2]], rosm) ## mz; its not the
        ## position of the largest intensity!
        if (is.nan(rmat[i,1]) || is.na(rmat[i,1])) ##  R2.11 :  weighted.mean()  results in NA (not NaN) for zero weights
            rmat[i,1] <- mean(peakrange[i,1:2])

        rmat[i,2:3] <- peakrange[i,1:2]            ## mzmin, mzmax
        rmat[i,4] <- stime[iret[1]:iret[2]][iymax] ## rt
        rmat[i,5:6] <- peakrange[i,3:4]            ## rtmin, rtmax

        if (peakrange[i,3] <  stime[1] || peakrange[i,4] > stime[length(stime)] || is.nan(pwid)) {
            warning("getPeaks: Peak  m/z:",peakrange[i,1],"-",peakrange[i,2], ",  RT:",peakrange[i,3],"-",peakrange[i,4],
                    "is out of retention time range for this sample (",object@filepath,"), using zero intensity value.\n")
            rmat[i,7:8] <- 0
        } else {
            rmat[i,7] <- pwid*sum(ymax)  ## into
            rmat[i,8] <- ymax[iymax]     ## maxo
        }
    }
    invisible(rmat)

}

msn2xcmsRaw <- function(xmsn) {
    x <- deepCopy(xmsn)

    x@tic <- x@msnAcquisitionNum
                                        
    x@scantime <- x@msnRt          # Fake time in secs
    x@acquisitionNum <- x@msnAcquisitionNum
    x@scanindex <- x@msnScanindex

    x@env$mz <- x@env$msnMz
    x@env$intensity <- x@env$msnIntensity
    invisible(x)
}

