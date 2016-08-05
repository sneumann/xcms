## All functions for xcmsRaw methods should be put here.
#' @include DataClasses.R c.R

xcmsRaw <- function(filename, profstep = 1, profmethod = "bin",
                    profparam = list(),
                    includeMSn = FALSE, mslevel=NULL,
                    scanrange=NULL) {

    object <- new("xcmsRaw")
    object@env <- new.env(parent=.GlobalEnv)
    object@filepath <- xcmsSource(filename)
    rawdata <- loadRaw(object@filepath, includeMSn)

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
    ## setting the scanrange.
    scanrange(object) <- as.numeric(scanrange)
    mslevel(object) <- as.numeric(mslevel)

    object@profmethod <- profmethod
    object@profparam <- profparam
    if (profstep)
        profStep(object) <- profstep

    if (!is.null(rawdata$acquisitionNum)) {
        ## defined only for mzData and mzXML
        object@acquisitionNum <- rawdata$acquisitionNum
    }

    if (!is.null(rawdata$polarity)) {
        object@polarity <- factor(rawdata$polarity,
                                  levels=c(0,1,-1),
                                  labels=c("negative", "positive", "unknown"));
    }

    if (!is.null(scanrange)) {
        ## Scanrange filtering
        keepidx <- seq.int(1, length(object@scantime)) %in% seq.int(scanrange[1], scanrange[2])
        object <- split(object, f=keepidx)[["TRUE"]]
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
        object@msnPrecursorScan <- match(rawdata$MSn$precursorNum, object@acquisitionNum)
        object@msnPrecursorMz <- rawdata$MSn$precursorMZ
        object@msnPrecursorIntensity <- rawdata$MSn$precursorIntensity
        object@msnPrecursorCharge <- rawdata$MSn$precursorCharge
        object@msnCollisionEnergy <- rawdata$MSn$collisionEnergy
    }

    if (!missing(mslevel) & !is.null(mslevel)) {
        object <- msn2ms(object)
        object <- split(object, f=object@msnLevel==mslevel)$"TRUE"
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
