## All Methods for xcmsRaw should be here.
#' @include functions-xcmsRaw.R

############################################################
## show
setMethod("show", "xcmsRaw", function(object) {

    cat("An \"xcmsRaw\" object with", length(object@scantime), "mass spectra\n\n")

    if (length(object@scantime)>0) {
        cat("Time range: ", paste(round(range(object@scantime), 1), collapse = "-"),
            " seconds (", paste(round(range(object@scantime)/60, 1), collapse = "-"),
            " minutes)\n", sep = "")
        cat("Mass range:", paste(round(range(object@env$mz), 4), collapse = "-"),
            "m/z\n")
        cat("Intensity range:", paste(signif(range(object@env$intensity), 6), collapse = "-"),
            "\n\n")
    }

    ## summary MSn data
    if (!is.null(object@msnLevel)) {
        cat("MSn data on ", length(unique(object@msnPrecursorMz)), " mass(es)\n")
        cat("\twith ", length(object@msnPrecursorMz)," MSn spectra\n")
    }

    cat("Profile method:", object@profmethod, "\n")
    cat("Profile step: ")

    if (is.null(object@env$profile))
        cat("no profile data\n")
    else {
        profmz <- profMz(object)
        cat(profStep(object), " m/z (", length(profmz), " grid points from ",
            paste(object@mzrange, collapse = " to "), " m/z)\n", sep = "")
    }
    if (length(object@profparam)) {
        cat("Profile parameters: ")
        for (i in seq(along = object@profparam)) {
            if (i != 1) cat("                    ")
            cat(names(object@profparam)[i], " = ", object@profparam[[i]], "\n", sep = "")
        }
    }

    memsize <- object.size(object)
    for (key in ls(object@env))
        memsize <- memsize + object.size(object@env[[key]])
    cat("\nMemory usage:", signif(memsize/2^20, 3), "MB\n")
})

############################################################
## sortMz
setMethod("revMz", "xcmsRaw", function(object) {

    for (i in 1:length(object@scanindex)) {
        idx <- (object@scanindex[i]+1):min(object@scanindex[i+1],
                                           length(object@env$mz), na.rm=TRUE)
        object@env$mz[idx] <- rev(object@env$mz[idx])
        object@env$intensity[idx] <- rev(object@env$intensity[idx])
    }
})

############################################################
## sortMz
setMethod("sortMz", "xcmsRaw", function(object) {

    for (i in 1:length(object@scanindex)) {
        idx <- (object@scanindex[i]+1):min(object@scanindex[i+1],
                                           length(object@env$mz), na.rm=TRUE)
        ord <- order(object@env$mz[idx])
        object@env$mz[idx] <- object@env$mz[idx[ord]]
        object@env$intensity[idx] <- object@env$intensity[idx[ord]]
    }
})

############################################################
## plotTIC
setMethod("plotTIC", "xcmsRaw", function(object, ident = FALSE, msident = FALSE) {

    if (all(object@tic == 0))
        points <- cbind(object@scantime, rawEIC(object,mzrange=range(object@env$mz))$intensity)  else
    points <- cbind(object@scantime, object@tic)

    plot(points, type="l", main="TIC Chromatogram", xlab="Seconds",
         ylab="Intensity")

    if (ident) {
        idx <- integer(0)
        ticdev <- dev.cur()
        if ((dev.cur()+1) %in% dev.list())
            msdev <- dev.cur()+1
        else
            msdev <- integer(0)
        while(length(id <- identify(points, labels = round(points[,1], 1), n = 1))) {
            idx <- c(idx, id)
            if (!length(msdev)) {
                options("device")$device()
                msdev <- dev.cur()
            }
            dev.set(msdev)
            plotScan(object, id, ident = msident)
            dev.set(ticdev)
        }
        return(idx)
    }

    invisible(points)
})

############################################################
## getScan
setMethod("getScan", "xcmsRaw", function(object, scan, mzrange = numeric()) {

    if (scan < 0)
        scan <- length(object@scantime) + 1 + scan

    idx <- seq(object@scanindex[scan]+1, min(object@scanindex[scan+1],
                                             length(object@env$mz), na.rm=TRUE))

    if (length(mzrange) >= 2) {
        mzrange <- range(mzrange)
        idx <- idx[object@env$mz[idx] >= mzrange[1] & object@env$mz[idx] <= mzrange[2]]
    }

    points <- cbind(mz = object@env$mz[idx], intensity = object@env$intensity[idx])

    invisible(points)
})

############################################################
## getSpec
setMethod("getSpec", "xcmsRaw", function(object, ...) {

    ## FIXME: unnecessary dependency on profile matrix?
    sel <- profRange(object, ...)

    scans <- list(length(sel$scanidx))
    uniquemz <- numeric()
    for (i in seq(along = sel$scanidx)) {
        scans[[i]] <- getScan(object, sel$scanidx[i], sel$mzrange)
        uniquemz <- unique(c(uniquemz, scans[[i]][,"mz"]))
    }
    uniquemz <- sort(uniquemz)

    intmat <- matrix(nrow = length(uniquemz), ncol = length(sel$scanidx))
    for (i in seq(along = sel$scanidx)) {
        scan <- getScan(object, sel$scanidx[i], sel$mzrange)
        intmat[,i] <- approx(scan, xout = uniquemz)$y
    }

    points <- cbind(mz = uniquemz, intensity = rowMeans(intmat))

    invisible(points)
})

############################################################
## findPeaks.matchedFilter
setMethod("findPeaks.matchedFilter", "xcmsRaw",
          function(object, fwhm = 30, sigma = fwhm/2.3548, max = 5,
                   snthresh = 10, step = 0.1, steps = 2,
                   mzdiff = 0.8 - step*steps, index = FALSE, sleep = 0,
                   verbose.columns = FALSE, scanrange= numeric()) {

    profFun <- match.profFun(object)

    scanrange.old <- scanrange
    ## sanitize if too few or too many scanrange is given
    if (length(scanrange) < 2)
        scanrange <- c(1, length(object@scantime))
    else
        scanrange <- range(scanrange)

    ## restrict and sanitize scanrange to maximally cover all scans
    scanrange[1] <- max(1,scanrange[1])
    scanrange[2] <- min(length(object@scantime),scanrange[2])

    ## Mild warning if the actual scanrange doesn't match the scanrange argument
    if (!(identical(scanrange.old,scanrange)) && (length(scanrange.old) >0)) {
        cat("Warning: scanrange was adjusted to ",scanrange,"\n")

        ## Scanrange filtering
        keepidx <- seq.int(1, length(object@scantime)) %in% seq.int(scanrange[1], scanrange[2])
        object <- split(object, f=keepidx)[["TRUE"]]
    }


### Create EIC buffer
    mrange <- range(object@env$mz)
    mass <- seq(floor(mrange[1]/step)*step, ceiling(mrange[2]/step)*step, by = step)
    bufsize <- min(100, length(mass))
    buf <- profFun(object@env$mz, object@env$intensity, object@scanindex,
                   bufsize, mass[1], mass[bufsize], TRUE, object@profparam)
    bufMax <- profMaxIdxM(object@env$mz, object@env$intensity, object@scanindex,
                          bufsize, mass[1], mass[bufsize], TRUE,
                          object@profparam)
    bufidx <- integer(length(mass))
    idxrange <- c(1, bufsize)
    bufidx[idxrange[1]:idxrange[2]] <- 1:bufsize
    lookahead <- steps-1
    lookbehind <- 1

    scantime <- object@scantime
    N <- nextn(length(scantime))
    xrange <- range(scantime)
    x <- c(0:(N/2), -(ceiling(N/2-1)):-1)*(xrange[2]-xrange[1])/(length(scantime)-1)

    filt <- -attr(eval(deriv3(~ 1/(sigma*sqrt(2*pi))*exp(-x^2/(2*sigma^2)), "x")), "hessian")
    filt <- filt/sqrt(sum(filt^2))
    filt <- fft(filt, inverse = TRUE)/length(filt)

    cnames <- c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "into", "intf", "maxo", "maxf", "i", "sn")
    rmat <- matrix(nrow = 2048, ncol = length(cnames))
    num <- 0

    for (i in seq(length = length(mass)-steps+1)) {
        if (i %% 500 == 0) {
            cat(round(mass[i]), ":", num, " ", sep = "")
            flush.console()
        }
### Update EIC buffer if necessary
        if (bufidx[i+lookahead] == 0) {
            bufidx[idxrange[1]:idxrange[2]] <- 0
            idxrange <- c(max(1, i - lookbehind), min(bufsize+i-1-lookbehind, length(mass)))
            bufidx[idxrange[1]:idxrange[2]] <- 1:(diff(idxrange)+1)
            buf <- profFun(object@env$mz, object@env$intensity, object@scanindex,
                           diff(idxrange)+1, mass[idxrange[1]], mass[idxrange[2]],
                           TRUE, object@profparam)
            bufMax <- profMaxIdxM(object@env$mz, object@env$intensity, object@scanindex,
                                  diff(idxrange)+1, mass[idxrange[1]], mass[idxrange[2]],
                                  TRUE, object@profparam)
        }
        ymat <- buf[bufidx[i:(i+steps-1)],,drop=FALSE]
        ysums <- colMax(ymat)
        yfilt <- filtfft(ysums, filt)
        gmax <- max(yfilt)
        for (j in seq(length = max)) {
            maxy <- which.max(yfilt)
            noise <- mean(ysums[ysums > 0])
            ##noise <- mean(yfilt[yfilt >= 0])
            sn <- yfilt[maxy]/noise
            if (yfilt[maxy] > 0 && yfilt[maxy] > snthresh*noise && ysums[maxy] > 0) {
                peakrange <- descendZero(yfilt, maxy)
                intmat <- ymat[,peakrange[1]:peakrange[2],drop=FALSE]
                mzmat <- matrix(object@env$mz[bufMax[bufidx[i:(i+steps-1)],
                                                     peakrange[1]:peakrange[2]]],
                                nrow = steps)
                which.intMax <- which.colMax(intmat)
                mzmat <- mzmat[which.intMax]
                if (all(is.na(mzmat))) {
                    yfilt[peakrange[1]:peakrange[2]] <- 0
                    next
                }
                mzrange <- range(mzmat, na.rm = TRUE)
                massmean <- weighted.mean(mzmat, intmat[which.intMax], na.rm = TRUE)
                ## This case (the only non-na m/z had intensity 0) was reported
                ## by Gregory Alan Barding "binlin processing"
                if(any(is.na(massmean))) {
                    massmean <- mean(mzmat, na.rm = TRUE)
                }

                pwid <- (scantime[peakrange[2]] - scantime[peakrange[1]])/(peakrange[2] - peakrange[1])
                into <- pwid*sum(ysums[peakrange[1]:peakrange[2]])
                intf <- pwid*sum(yfilt[peakrange[1]:peakrange[2]])
                maxo <- max(ysums[peakrange[1]:peakrange[2]])
                maxf <- yfilt[maxy]
                if (sleep > 0) {
                    plot(scantime, yfilt, type = "l", main = paste(mass[i], "-", mass[i+1]), ylim=c(-gmax/3, gmax))
                    points(cbind(scantime, yfilt)[peakrange[1]:peakrange[2],], type = "l", col = "red")
                    points(scantime, colSums(ymat), type = "l", col = "blue", lty = "dashed")
                    abline(h = snthresh*noise, col = "red")
                    Sys.sleep(sleep)
                }
                yfilt[peakrange[1]:peakrange[2]] <- 0
                num <- num + 1
### Double the size of the output matrix if it's full
                if (num > nrow(rmat)) {
                    nrmat <- matrix(nrow = 2*nrow(rmat), ncol = ncol(rmat))
                    nrmat[seq(length = nrow(rmat)),] = rmat
                    rmat <- nrmat
                }
                rmat[num,] <- c(massmean, mzrange[1], mzrange[2], maxy, peakrange, into, intf, maxo, maxf, j, sn)
            } else
                break
        }
    }
    cat("\n")
    colnames(rmat) <- cnames
    rmat <- rmat[seq(length = num),]
    max <- max-1 + max*(steps-1) + max*ceiling(mzdiff/step)
    if (index)
        mzdiff <- mzdiff/step
    else {
        rmat[,"rt"] <- scantime[rmat[,"rt"]]
        rmat[,"rtmin"] <- scantime[rmat[,"rtmin"]]
        rmat[,"rtmax"] <- scantime[rmat[,"rtmax"]]
    }
    uorder <- order(rmat[,"into"], decreasing=TRUE)
    uindex <- rectUnique(rmat[,c("mzmin","mzmax","rtmin","rtmax"),drop=FALSE],
                         uorder, mzdiff)
    rmat <- rmat[uindex,,drop=FALSE]
    invisible(new("xcmsPeaks", rmat))
})

############################################################
## findPeaks.centWave
setMethod("findPeaks.centWave", "xcmsRaw", function(object, ppm=25, peakwidth=c(20,50), snthresh=10,
                                                    prefilter=c(3,100), mzCenterFun="wMean", integrate=1, mzdiff=-0.001,
                                                    fitgauss=FALSE, scanrange= numeric(), noise=0, ## noise.local=TRUE,
                                                    sleep=0, verbose.columns=FALSE, ROI.list=list(), 
                                                    firstBaselineCheck=TRUE, roiScales=NULL) {
    if (!isCentroided(object))
        warning("It looks like this file is in profile mode. centWave can process only centroid mode data !\n")

    mzCenterFun <- paste("mzCenter", mzCenterFun, sep=".")
    if (!exists(mzCenterFun, mode="function"))
        stop("Error: >",mzCenterFun,"< not defined ! \n")

    if (!is.logical(firstBaselineCheck))
      stop("Error: parameter >firstBaselineCheck< is not a vector ! \n")
    if (length(firstBaselineCheck) != 1)
      stop("Error: parameter >firstBaselineCheck< is not a single logical ! \n")
    if (!is.null(roiScales)){
      if (!is.vector(roiScales))
        stop("Error: parameter >roiScales< is not a vector ! \n")
      if(!is.numeric(roiScales))
        stop("Error: parameter >roiScales< is not a vector of type numeric ! \n")
      if(!length(roiScales) == length(ROI.list))
        stop("Error: length of parameter >roiScales< is not equal to the length of parameter >ROI.list< ! \n")
    }
    
    scanrange.old <- scanrange
    if (length(scanrange) < 2)
        scanrange <- c(1, length(object@scantime))
    else
        scanrange <- range(scanrange)

    scanrange[1] <- max(1,scanrange[1])
    scanrange[2] <- min(length(object@scantime),scanrange[2])

    if (!(identical(scanrange.old,scanrange)) && (length(scanrange.old) >0))
        cat("Warning: scanrange was adjusted to ",scanrange,"\n")

    basenames <- c("mz","mzmin","mzmax","rt","rtmin","rtmax","into","intb","maxo","sn")
    verbosenames <- c("egauss","mu","sigma","h","f", "dppm", "scale","scpos","scmin","scmax","lmin","lmax")

    ## Peak width: seconds to scales
    scalerange <- round((peakwidth / mean(diff(object@scantime))) / 2)

    if (length(z <- which(scalerange==0)))
        scalerange <- scalerange[-z]

    if (length(scalerange) < 1)
        stop("No scales ? Please check peak width!\n")

    if (length(scalerange) > 1)
        scales <- seq(from=scalerange[1], to=scalerange[2], by=2)  else
    scales <- scalerange;

    minPeakWidth <-  scales[1];
    noiserange <- c(minPeakWidth*3, max(scales)*3);
    maxGaussOverlap <- 0.5;
    minPtsAboveBaseLine <- max(4,minPeakWidth-2);
    minCentroids <- minPtsAboveBaseLine ;
    scRangeTol <-  maxDescOutlier <- floor(minPeakWidth/2);

    ## If no ROIs are supplied then search for them.
    if (length(ROI.list) == 0) {
        cat("\n Detecting mass traces at",ppm,"ppm ... \n"); flush.console();
        ROI.list <- findmzROI(object,scanrange=scanrange,dev=ppm * 1e-6,minCentroids=minCentroids, prefilter=prefilter, noise=noise)
        if (length(ROI.list) == 0) {
            cat("No ROIs found ! \n")

            if (verbose.columns) {
                nopeaks <- new("xcmsPeaks", matrix(nrow=0, ncol=length(basenames)+length(verbosenames)))
                colnames(nopeaks) <- c(basenames, verbosenames)
            } else {
                nopeaks <- new("xcmsPeaks", matrix(nrow=0, ncol=length(basenames)))
                colnames(nopeaks) <- c(basenames)
            }

            return(invisible(nopeaks))
        }
    }

    peaklist <- list()
    scantime <- object@scantime
    Nscantime <- length(scantime)
    lf <- length(ROI.list)

    cat('\n Detecting chromatographic peaks ... \n % finished: '); lp <- -1;

    for (f in  1:lf) {

        ## Show progress
        perc <- round((f/lf) * 100)
        if ((perc %% 10 == 0) && (perc != lp))
        {
            cat(perc," ",sep="");
            lp <- perc;
        }
        flush.console()

        feat <- ROI.list[[f]]
        N <- feat$scmax - feat$scmin + 1

        peaks <- peakinfo <- NULL

        mzrange <- c(feat$mzmin,feat$mzmax)
        sccenter <- feat$scmin[1] + floor(N/2) - 1
        scrange <- c(feat$scmin,feat$scmax)
        ## scrange + noiserange, used for baseline detection and wavelet analysis
        sr <- c(max(scanrange[1],scrange[1] - max(noiserange)),min(scanrange[2],scrange[2] + max(noiserange)))
        eic <- rawEIC(object,mzrange=mzrange,scanrange=sr)
        d <- eic$intensity
        td <- sr[1]:sr[2]
        scan.range <- c(sr[1],sr[2])
        ## original mzROI range
        mzROI.EIC <- rawEIC(object,mzrange=mzrange,scanrange=scrange)
        omz <- rawMZ(object,mzrange=mzrange,scanrange=scrange)
        if (all(omz == 0))
            stop("centWave: debug me: (omz == 0)?\n")
        od  <- mzROI.EIC$intensity
        otd <- mzROI.EIC$scan
        if (all(od == 0))
            stop("centWave: debug me: (all(od == 0))?\n")

        ##  scrange + scRangeTol, used for gauss fitting and continuous data above 1st baseline detection
        ftd <- max(td[1], scrange[1] - scRangeTol) : min(td[length(td)], scrange[2] + scRangeTol)
        fd <- d[match(ftd,td)]

        ## 1st type of baseline: statistic approach
        if (N >= 10*minPeakWidth)  ## in case of very long mass trace use full scan range for baseline detection
            noised <- rawEIC(object,mzrange=mzrange,scanrange=scanrange)$intensity else
        noised <- d;
        ## 90% trimmed mean as first baseline guess
        noise <- estimateChromNoise(noised, trim=0.05, minPts=3*minPeakWidth)

        ## any continuous data above 1st baseline ?
        if (firstBaselineCheck & !continuousPtsAboveThreshold(fd,threshold=noise,num=minPtsAboveBaseLine))
            next;

        ## 2nd baseline estimate using not-peak-range
        lnoise <- getLocalNoiseEstimate(d,td,ftd,noiserange,Nscantime, threshold=noise,num=minPtsAboveBaseLine)

        ## Final baseline & Noise estimate
        baseline <- max(1,min(lnoise[1],noise))
        sdnoise <- max(1,lnoise[2])
        sdthr <-  sdnoise * snthresh

        ## is there any data above S/N * threshold ?
        if (!(any(fd - baseline >= sdthr)))
            next;

        wCoefs <- MSW.cwt(d, scales=scales, wavelet='mexh')
        if  (!(!is.null(dim(wCoefs)) && any(wCoefs- baseline >= sdthr)))
            next;

        if (td[length(td)] == Nscantime) ## workaround, localMax fails otherwise
            wCoefs[nrow(wCoefs),] <- wCoefs[nrow(wCoefs)-1,] * 0.99
        localMax <- MSW.getLocalMaximumCWT(wCoefs)
        rL <- MSW.getRidge(localMax)
        wpeaks <- sapply(rL,
                         function(x) {
                             w <- min(1:length(x),ncol(wCoefs))
                             any(wCoefs[x,w]- baseline >= sdthr)
                         })
        if (any(wpeaks)) {
            wpeaksidx <- which(wpeaks)
            ## check each peak in ridgeList
            for (p in 1:length(wpeaksidx)) {
                opp <- rL[[wpeaksidx[p]]]
                pp <- unique(opp)
                if (length(pp) >= 1) {
                    dv <- td[pp] %in% ftd
                    if (any(dv)) { ## peaks in orig. data range
                        ## Final S/N check
                        if (any(d[pp[dv]]- baseline >= sdthr)) {
                            if(!is.null(roiScales)){
                                ## use given scale
                                best.scale.nr <- which(scales == roiScales[[f]])
                                if(best.scale.nr > length(opp))
                                    best.scale.nr <- length(opp)
                            } else {
                                ## try to decide which scale describes the peak best
                                inti <- numeric(length(opp))
                                irange = rep(ceiling(scales[1]/2),length(opp))
                                for (k in 1:length(opp)) {
                                    kpos <- opp[k]
                                    r1 <- ifelse(kpos-irange[k] > 1,kpos-irange[k],1)
                                    r2 <- ifelse(kpos+irange[k] < length(d),kpos+irange[k],length(d))
                                    inti[k] <- sum(d[r1:r2])
                                }
                                maxpi <- which.max(inti)
                                if (length(maxpi) > 1) {
                                    m <- wCoefs[opp[maxpi],maxpi]
                                    bestcol <- which(m == max(m),arr.ind=T)[2]
                                    best.scale.nr <- maxpi[bestcol]
                                } else  best.scale.nr <- maxpi
                            }

                            best.scale <-  scales[best.scale.nr]
                            best.scale.pos <- opp[best.scale.nr]

                            pprange <- min(pp):max(pp)
                            ## maxint <- max(d[pprange])
                            lwpos <- max(1,best.scale.pos - best.scale)
                            rwpos <- min(best.scale.pos + best.scale,length(td))
                            p1 <- match(td[lwpos],otd)[1]
                            p2 <- match(td[rwpos],otd); p2 <- p2[length(p2)]
                            if (is.na(p1)) p1<-1
                            if (is.na(p2)) p2<-N
                            mz.value <- omz[p1:p2]
                            mz.int <- od[p1:p2]
                            maxint <- max(mz.int)

                            ## re-calculate m/z value for peak range
                            mzrange <- range(mz.value)
                            mzmean <- do.call(mzCenterFun,list(mz=mz.value,intensity=mz.int))

                            ## Compute dppm only if needed
                            dppm <- NA
                            if (verbose.columns)
                                if (length(mz.value) >= (minCentroids+1))
                                    dppm <- round(min(running(abs(diff(mz.value)) /(mzrange[2] *  1e-6),fun=max,width=minCentroids))) else
                            dppm <- round((mzrange[2]-mzrange[1]) /  (mzrange[2] *  1e-6))

                            peaks <- rbind(peaks,
                                           c(mzmean,mzrange,           ## mz
                                             NA,NA,NA,                   ## rt, rtmin, rtmax,
                                             NA,                         ## intensity (sum)
                                             NA,                         ## intensity (-bl)
                                             maxint,                     ## max intensity
                                             round((maxint - baseline) / sdnoise),  ##  S/N Ratio
                                             NA,                         ## Gaussian RMSE
                                             NA,NA,NA,                   ## Gaussian Parameters
                                             f,                          ## ROI Position
                                             dppm,                       ## max. difference between the [minCentroids] peaks in ppm
                                             best.scale,                 ## Scale
                                             td[best.scale.pos], td[lwpos], td[rwpos],  ## Peak positions guessed from the wavelet's (scan nr)
                                             NA,NA ))                    ## Peak limits (scan nr)

                            peakinfo <- rbind(peakinfo,c(best.scale, best.scale.nr, best.scale.pos, lwpos, rwpos))  ## Peak positions guessed from the wavelet's
                        }
                    }
                }
            }  ##for
        } ## if


        ##  postprocessing
        if (!is.null(peaks)) {
            colnames(peaks) <- c(basenames, verbosenames)

            colnames(peakinfo) <- c("scale","scaleNr","scpos","scmin","scmax")

            for (p in 1:dim(peaks)[1]) {
                ## find minima, assign rt and intensity values
                if (integrate == 1) {
                    lm <- descendMin(wCoefs[,peakinfo[p,"scaleNr"]], istart= peakinfo[p,"scpos"])
                    gap <- all(d[lm[1]:lm[2]] == 0) ## looks like we got stuck in a gap right in the middle of the peak
                    if ((lm[1]==lm[2]) || gap )## fall-back
                        lm <- descendMinTol(d, startpos=c(peakinfo[p,"scmin"], peakinfo[p,"scmax"]), maxDescOutlier)
                } else
                    lm <- descendMinTol(d,startpos=c(peakinfo[p,"scmin"],peakinfo[p,"scmax"]),maxDescOutlier)

                ## narrow down peak rt boundaries by skipping zeros
                pd <- d[lm[1]:lm[2]]; np <- length(pd)
                lm.l <-  xcms:::findEqualGreaterUnsorted(pd,1)
                lm.l <- max(1, lm.l - 1)
                lm.r <- xcms:::findEqualGreaterUnsorted(rev(pd),1)
                lm.r <- max(1, lm.r - 1)
                lm <- lm + c(lm.l - 1, -(lm.r - 1) )

                peakrange <- td[lm]
                peaks[p,"rtmin"] <- scantime[peakrange[1]]
                peaks[p,"rtmax"] <- scantime[peakrange[2]]

                peaks[p,"maxo"] <- max(d[lm[1]:lm[2]])

                pwid <- (scantime[peakrange[2]] - scantime[peakrange[1]])/(peakrange[2] - peakrange[1])
                if (is.na(pwid))
                    pwid <- 1

                peaks[p,"into"] <- pwid*sum(d[lm[1]:lm[2]])

                db <-  d[lm[1]:lm[2]] - baseline
                peaks[p,"intb"] <- pwid*sum(db[db>0])

                peaks[p,"lmin"] <- lm[1];
                peaks[p,"lmax"] <- lm[2];

                if (fitgauss) {
                    ## perform gaussian fits, use wavelets for inital parameters
                    md <- max(d[lm[1]:lm[2]]);d1 <- d[lm[1]:lm[2]]/md; ## normalize data for gaussian error calc.
                    pgauss <- fitGauss(td[lm[1]:lm[2]],d[lm[1]:lm[2]],pgauss =
                                       list(mu=peaks[p,"scpos"],sigma=peaks[p,"scmax"]-peaks[p,"scmin"],h=peaks[p,"maxo"]))
                    rtime <- peaks[p,"scpos"]
                    if (!any(is.na(pgauss)) && all(pgauss > 0)) {
                        gtime <- td[match(round(pgauss$mu),td)]
                        if (!is.na(gtime)) {
                            rtime <- gtime
                            peaks[p,"mu"] <- pgauss$mu; peaks[p,"sigma"] <- pgauss$sigma; peaks[p,"h"] <- pgauss$h;
                            peaks[p,"egauss"] <- sqrt((1/length(td[lm[1]:lm[2]])) * sum(((d1-gauss(td[lm[1]:lm[2]],pgauss$h/md,pgauss$mu,pgauss$sigma))^2)))
                        }
                    }
                    peaks[p,"rt"] <- scantime[rtime]
                    ## avoid fitting side effects
                    if (peaks[p,"rt"] < peaks[p,"rtmin"])
                        peaks[p,"rt"] <- scantime[peaks[p,"scpos"]]
                } else
                    peaks[p,"rt"] <- scantime[peaks[p,"scpos"]]
            }
            peaks <- joinOverlappingPeaks(td,d,otd,omz,od,scantime,scan.range,peaks,maxGaussOverlap,mzCenterFun=mzCenterFun)
        }



        if ((sleep >0) && (!is.null(peaks))) {
            tdp <- scantime[td]; trange <- range(tdp)
            egauss <- paste(round(peaks[,"egauss"],3),collapse=", ")
            cdppm <- paste(peaks[,"dppm"],collapse=", ")
            csn <- paste(peaks[,"sn"],collapse=", ")
            par(bg = "white")
            l <- layout(matrix(c(1,2,3),nrow=3,ncol=1,byrow=T),heights=c(.5,.75,2));
            par(mar= c(2, 4, 4, 2) + 0.1)
            plotRaw(object,mzrange=mzrange,rtrange=trange,log=TRUE,title='')
            title(main=paste(f,': ', round(mzrange[1],4),' - ',round(mzrange[2],4),' m/z , dppm=',cdppm,', EGauss=',egauss ,',  S/N =',csn,sep=''))
            par(mar= c(1, 4, 1, 2) + 0.1)
            image(y=scales[1:(dim(wCoefs)[2])],z=wCoefs,col=terrain.colors(256),xaxt='n',ylab='CWT coeff.')
            par(mar= c(4, 4, 1, 2) + 0.1)
            plot(tdp,d,ylab='Intensity',xlab='Scan Time');lines(tdp,d,lty=2)
            lines(scantime[otd],od,lty=2,col='blue') ## original mzbox range
            abline(h=baseline,col='green')
            bwh <- length(sr[1]:sr[2]) - length(baseline)
            if (odd(bwh)) {bwh1 <-  floor(bwh/2); bwh2 <- bwh1+1} else {bwh1<-bwh2<-bwh/2}
            if  (any(!is.na(peaks[,"scpos"])))
            {   ## plot centers and width found through wavelet analysis
                abline(v=scantime[na.omit(peaks[(peaks[,"scpos"] >0),"scpos"])],col='red')
            }
            abline(v=na.omit(c(peaks[,"rtmin"],peaks[,"rtmax"])),col='green',lwd=1)
            if (fitgauss) {
                tdx <- seq(min(td),max(td),length.out=200)
                tdxp <- seq(trange[1],trange[2],length.out=200)
                fitted.peaks <- which(!is.na(peaks[,"mu"]))
                for (p in fitted.peaks)
                {   ## plot gaussian fits
                    yg<-gauss(tdx,peaks[p,"h"],peaks[p,"mu"],peaks[p,"sigma"])
                    lines(tdxp,yg,col='blue')
                }
            }
            Sys.sleep(sleep)
        }

        if (!is.null(peaks)) {
            peaklist[[length(peaklist)+1]] <- peaks
        }

    } ## f

    if (length(peaklist) == 0) {
        cat("\nNo peaks found !\n")

        if (verbose.columns) {
            nopeaks <- new("xcmsPeaks", matrix(nrow=0, ncol=length(basenames)+length(verbosenames)))
            colnames(nopeaks) <- c(basenames, verbosenames)
        } else {
            nopeaks <- new("xcmsPeaks", matrix(nrow=0, ncol=length(basenames)))
            colnames(nopeaks) <- c(basenames)
        }

        return(invisible(nopeaks))
    }

    p <- do.call(rbind,peaklist)

    if (!verbose.columns)
        p <- p[,basenames,drop=FALSE]

    uorder <- order(p[,"into"], decreasing=TRUE)
    pm <- as.matrix(p[,c("mzmin","mzmax","rtmin","rtmax"),drop=FALSE])
    uindex <- rectUnique(pm,uorder,mzdiff,ydiff = -0.00001) ## allow adjacent peaks
    pr <- p[uindex,,drop=FALSE]
    cat("\n",dim(pr)[1]," Peaks.\n")

    invisible(new("xcmsPeaks", pr))
})

############################################################
## findPeaks.centWaveWithPredictedIsotopeROIs
setMethod("findPeaks.centWaveWithPredictedIsotopeROIs", "xcmsRaw", function(object, ppm=25, peakwidth=c(20,50), snthresh=10, 
                                                                            prefilter=c(3,100), mzCenterFun="wMean", integrate=1, mzdiff=-0.001,
                                                                            fitgauss=FALSE, scanrange= numeric(), noise=0, ## noise.local=TRUE,
                                                                            sleep=0, verbose.columns=FALSE, ROI.list=list(), 
                                                                            firstBaselineCheck=TRUE, roiScales=NULL, snthreshIsoROIs=6.25, maxcharge=3, maxiso=5, mzIntervalExtension=TRUE) {
  ## perform tradictional peak picking
  xcmsPeaks <- findPeaks.centWave(
    object = object, ppm=ppm, peakwidth=peakwidth, snthresh=snthresh,
    prefilter=prefilter, mzCenterFun=mzCenterFun, integrate=integrate, mzdiff=mzdiff,
    fitgauss=fitgauss, scanrange=scanrange, noise=noise, ## noise.local=noise.local,
    sleep=sleep, verbose.columns=TRUE, ROI.list=ROI.list,
    firstBaselineCheck=firstBaselineCheck, roiScales=roiScales
  )
  
  xcmsPeaksWithAdditionalIsotopeFeatures <- findPeaks.addPredictedIsotopeFeatures(
    object = object, ppm=ppm, peakwidth=peakwidth, 
    prefilter=prefilter, mzCenterFun=mzCenterFun, integrate=integrate, mzdiff=mzdiff,
    fitgauss=fitgauss, scanrange=scanrange, noise=noise, ## noise.local=noise.local,
    sleep=sleep, verbose.columns=verbose.columns, 
    xcmsPeaks=xcmsPeaks, snthresh=snthreshIsoROIs, maxcharge=maxcharge, maxiso=maxiso, mzIntervalExtension=mzIntervalExtension
  )
  
  return(xcmsPeaksWithAdditionalIsotopeFeatures)
})
setMethod("findPeaks.addPredictedIsotopeFeatures", "xcmsRaw", function(object, ppm=25, peakwidth=c(20,50), 
                                                                            prefilter=c(3,100), mzCenterFun="wMean", integrate=1, mzdiff=-0.001,
                                                                            fitgauss=FALSE, scanrange= numeric(), noise=0, ## noise.local=TRUE,
                                                                            sleep=0, verbose.columns=FALSE, 
                                                                            xcmsPeaks, snthresh=6.25, maxcharge=3, maxiso=5, mzIntervalExtension=TRUE) {
  if(nrow(xcmsPeaks) == 0){
    warning("Warning: There are no features (parameter >xcmsPeaks<) for the prediction of isotope ROIs !\n")
    return(xcmsPeaks)
  }
  if(class(xcmsPeaks) != "xcmsPeaks")
    stop("Error: parameter >xcmsPeaks< is not of class 'xcmsPeaks' ! \n")
  if(any(is.na(match(x = c("scmin", "scmax"), table = colnames(xcmsPeaks)))))
    stop("Error: peak list >xcmsPeaks< is missing the columns 'scmin' and 'scmax' ! Please set parameter >verbose.columns< to TRUE for peak picking with 'centWave' and try again ! \n")
  
  addNewIsotopeROIs <- TRUE
  addNewAdductROIs  <- FALSE
  polarity <- object@polarity
  
  ####################################################################################
  ## predict new ROIs
  
  ## convert present peaks to list of lists
  presentROIs.list <- list()
  for(peakIdx in 1:nrow(xcmsPeaks)){
    presentROIs.list[[peakIdx]] <- list(
      mz        = xcmsPeaks[[peakIdx, "mz"]],## XXX not used!
      mzmin     = xcmsPeaks[[peakIdx, "mzmin"]],
      mzmax     = xcmsPeaks[[peakIdx, "mzmax"]],
      scmin     = xcmsPeaks[[peakIdx, "scmin"]],
      scmax     = xcmsPeaks[[peakIdx, "scmax"]],
      length    = -1,## XXX not used!
      intensity = xcmsPeaks[[peakIdx, "intb"]],## XXX not used!
      scale     = xcmsPeaks[[peakIdx, "scale"]]## XXX not used!
    )
    
    if(abs(xcmsPeaks[[peakIdx, "mzmax"]] - xcmsPeaks[[peakIdx, "mzmin"]]) < xcmsPeaks[[peakIdx, "mz"]] * ppm / 1E6){
      presentROIs.list[[peakIdx]]$mzmin <- xcmsPeaks[[peakIdx, "mz"]] - xcmsPeaks[[peakIdx, "mz"]] * (ppm/2) / 1E6
      presentROIs.list[[peakIdx]]$mzmax <- xcmsPeaks[[peakIdx, "mz"]] + xcmsPeaks[[peakIdx, "mz"]] * (ppm/2) / 1E6
    }
  }
  
  ## fetch predicted ROIs
  resultObj <- createAdditionalROIs(object, presentROIs.list, ppm, addNewIsotopeROIs, maxcharge, maxiso, mzIntervalExtension, addNewAdductROIs, polarity)
  newRoiCounter <- resultObj$newRoiCounter
  numberOfAdditionalIsotopeROIs <- resultObj$numberOfAdditionalIsotopeROIs
  numberOfAdditionalAdductROIs <- resultObj$numberOfAdditionalAdductROIs
  newROI.matrix <- resultObj$newROI.matrix
  
  if(nrow(newROI.matrix) == 0)
    return(xcmsPeaks)
  
  ## remove ROIs with weak signal content
  intensityThreshold <- 10
  newROI.matrix <- removeROIsWithoutSignal(object, newROI.matrix, intensityThreshold)
  
  ## convert to list of lists
  newROI.list <- list()
  for(idx in 1:nrow(newROI.matrix))
    ## c("mz", "mzmin", "mzmax", "scmin", "scmax", "length", "intensity")
    newROI.list[[length(newROI.list) + 1]] <- as.list(newROI.matrix[idx, ])
  
  cat("Predicted ROIs: ", length(newROI.list), " new ROIs (", numberOfAdditionalIsotopeROIs, " isotope ROIs, ", numberOfAdditionalAdductROIs, " adduct ROIs) for ", length(presentROIs.list)," present ROIs.", "\n")
  
  ####################################################################################
  ## perform peak picking for predicted ROIs
  roiScales <- unlist(lapply(X = newROI.list, FUN = function(x){x$scale}))
  xcmsPeaks2 <- findPeaks.centWave(
    object = object, ppm=ppm, peakwidth=peakwidth, snthresh=snthresh,
    prefilter=prefilter, mzCenterFun=mzCenterFun, integrate=integrate, mzdiff=mzdiff,
    fitgauss=fitgauss, scanrange=scanrange, noise=noise, ## noise.local=noise.local,
    sleep=sleep, verbose.columns=verbose.columns, ROI.list=newROI.list, firstBaselineCheck=FALSE, roiScales=roiScales
  )
  
  if(nrow(xcmsPeaks2) > 0){
    ## remove NaN values
    rowsWithNaN <- which(apply(X = xcmsPeaks2[, c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax")], MARGIN = 1, FUN = function(x){any(is.na(x))}))
    if(length(rowsWithNaN) > 0)
      xcmsPeaks2 <- xcmsPeaks2[-rowsWithNaN, ]
    
    noArea <- which((xcmsPeaks2[, "mzmax"] - xcmsPeaks2[, "mzmin"]) == 0 || (xcmsPeaks2[, "rtmax"] - xcmsPeaks2[, "rtmin"]) == 0)
    if(length(noArea) > 0)
      xcmsPeaks2 <- xcmsPeaks2[-noArea, ]
  }
  
  ## make present peaks and new peaks distinct by removing overlapping peaks
  if(nrow(xcmsPeaks2) > 0){
    ## remove ROIs which are already there
    overlapProportionThreshold <- 0.01
    drop <- apply(X = xcmsPeaks2, MARGIN = 1, FUN = function(x){
      roiInt  <- x[["into"]]
      peakInt  <- xcmsPeaks[, "into"]
      roiMzMin  <- x[["mzmin"]]
      roiMzMax  <- x[["mzmax"]]
      peakMzMin <- xcmsPeaks[, "mzmin"]
      peakMzMax <- xcmsPeaks[, "mzmax"]
      roiMzCenter  = (roiMzMin  + roiMzMax ) / 2;
      peakMzCenter = (peakMzMin + peakMzMax) / 2;
      roiMzRadius  = (roiMzMax  - roiMzMin ) / 2;
      peakMzRadius = (peakMzMax - peakMzMin) / 2;
      overlappingmz <- abs(peakMzCenter - roiMzCenter) <= (roiMzRadius + peakMzRadius)
      
      roiRtMin  <- x[["rtmin"]]
      roiRtMax  <- x[["rtmax"]]
      peakRtMin <- xcmsPeaks[, "rtmin"]
      peakRtMax <- xcmsPeaks[, "rtmax"]
      roiRtCenter  = (roiRtMin  + roiRtMax ) / 2;
      peakRtCenter = (peakRtMin + peakRtMax) / 2;
      roiRtRadius  = (roiRtMax  - roiRtMin ) / 2;
      peakRtRadius = (peakRtMax - peakRtMin) / 2;
      overlappingrt <- abs(peakRtCenter - roiRtCenter) <= (roiRtRadius + peakRtRadius)
      
      overlapping <- overlappingmz & overlappingrt
      
      overlappingPeaks <- which(overlapping)
      overlappingPeaksInt <- peakInt[overlappingPeaks]
      
      removeROI <- FALSE
      peaksToRemove <- NULL
      if(any(overlapping)){
        if(any(overlappingPeaksInt > roiInt))
          return(TRUE)
        else
          return(overlappingPeaks)
      } else {
        ## no overlap
        return(FALSE)
      }
      
      return(isOverlap)
    })
    
    removeROI <- unlist(lapply(X = drop, FUN = function(x){
      if(is.logical(x)){
        return(x)
      } else {
        return(FALSE)
      }
    }))
    removePeaks <- unique(unlist(lapply(X = drop, FUN = function(x){
      if(is.logical(x)){
        return(NULL)
      } else {
        return(x)
      }
    })))
    
    if(length(removePeaks) > 0)
      xcmsPeaks <- xcmsPeaks[-removePeaks, ]
    xcmsPeaks2 <- xcmsPeaks2[!removeROI, ]
  }
  
  ## merge result with present results
  if(!verbose.columns)
    xcmsPeaks <- xcmsPeaks[, c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "into", "intb", "maxo", "sn")]
  
  xcmsPeaks <- rbind(xcmsPeaks, xcmsPeaks2)
  
  return(xcmsPeaks)
})

removeROIsOutOfRange <- function(object, roi.matrix){
  ## c("mz", "mzmin", "mzmax", "scmin", "scmax", "length", "intensity")
  numberOfROIs <- nrow(roi.matrix)
  
  minMz <- min(object@env$mz)
  maxMz <- max(object@env$mz)
  minScanRange <- 1
  maxScanRange <- length(object@scantime)
  #minScanRange <- min(object@scantime)
  #maxScanRange <- max(object@scantime)
  
  roiWithinRange <- rep(x = TRUE, times = numberOfROIs)
  roiWithinRange <- roiWithinRange & (roi.matrix[, "mzmin"] >= minMz)
  roiWithinRange <- roiWithinRange & (roi.matrix[, "mzmax"] <= maxMz)
  roiWithinRange <- roiWithinRange & (roi.matrix[, "scmin"] >= minScanRange)
  roiWithinRange <- roiWithinRange & (roi.matrix[, "scmax"] <= maxScanRange)
  
  roi.matrix <- roi.matrix[roiWithinRange, ]
  
  return(roi.matrix)
}
removeROIsWithoutSignal <- function(object, roi.matrix, intensityThreshold){
  ## c("mz", "mzmin", "mzmax", "scmin", "scmax", "length", "intensity")
  numberOfROIs <- nrow(roi.matrix)
  sufficientSignalThere <- rep(x = TRUE, times = numberOfROIs)
  for(roiIdx in seq_len(numberOfROIs)){
    mzrange   <- c(roi.matrix[[roiIdx, "mzmin"]], roi.matrix[[roiIdx, "mzmax"]])
    scanrange <- c(roi.matrix[[roiIdx, "scmin"]], roi.matrix[[roiIdx, "scmax"]])
    mzROI.EIC <- rawEIC(object, mzrange=mzrange, scanrange=scanrange)
    sumOfIntensities <- sum(mzROI.EIC$intensity)
    
    if(sumOfIntensities < intensityThreshold)
      sufficientSignalThere[[roiIdx]] <- FALSE
  }
  roi.matrix <- roi.matrix[sufficientSignalThere, ]
  
  return(roi.matrix)
}

createAdditionalROIs <- function(object, ROI.list, ppm, addNewIsotopeROIs, maxcharge, maxiso, mzIntervalExtension, addNewAdductROIs, polarity){
  ###############################################################################################
  ## isotope ROIs
  if(addNewIsotopeROIs){
    ## init
    isotopeDistance <- 1.0033548378
    charges <- 1:maxcharge
    isos <- 1:maxiso
    
    isotopeStepSizesForCharge <- list()
    for(charge in charges)
      isotopeStepSizesForCharge[[charge]] <- isotopeDistance / charge
    
    isotopeStepSizes <- list()
    for(charge in charges)
      isotopeStepSizes[[charge]] <- list()
    
    for(charge in charges)
      for(iso in isos)
        isotopeStepSizes[[charge]][[iso]] <- isotopeStepSizesForCharge[[charge]] * iso
    
    isotopePopulationMz <- list()
    for(charge in charges)
      for(iso in isos)
        isotopePopulationMz[[length(isotopePopulationMz) + 1]] <- isotopeStepSizes[[charge]][[iso]]
    isotopePopulationMz <- unlist(unique(isotopePopulationMz))
    
    numberOfIsotopeROIs <- length(ROI.list) * length(isotopePopulationMz)
    isotopeROIs.matrix <- matrix(nrow = numberOfIsotopeROIs, ncol = 8)
    colnames(isotopeROIs.matrix) <- c("mz", "mzmin", "mzmax", "scmin", "scmax", "length", "intensity", "scale")
    
    ## complement found ROIs
    for(roiIdx in 1:(length(ROI.list))){
      for(mzIdx in 1:length(isotopePopulationMz)){
        ## create new ROI!
        mzDifference <- isotopePopulationMz[[mzIdx]]
        if(mzIntervalExtension)
          ## extend m/z interval for weak peaks
          #mzIntervalExtension <- ROI.list[[roiIdx]]$mz * ppm / 1E6
          mzIntervalExtension <- (ROI.list[[roiIdx]]$mzmax - ROI.list[[roiIdx]]$mzmin) * 2
        else
          mzIntervalExtension <- 0
        
        idx <- (roiIdx - 1) * length(isotopePopulationMz) + mzIdx
        isotopeROIs.matrix[idx, ] <- c(
          ROI.list[[roiIdx]]$mz + mzDifference,## XXX not used!
          ROI.list[[roiIdx]]$mzmin + mzDifference - mzIntervalExtension,
          ROI.list[[roiIdx]]$mzmax + mzDifference + mzIntervalExtension,
          ROI.list[[roiIdx]]$scmin,
          ROI.list[[roiIdx]]$scmax,
          ROI.list[[roiIdx]]$length,## XXX not used!
          -1,  #ROI.list[[roiIdx]]$intensity ## XXX not used!
          ROI.list[[roiIdx]]$scale
        )
      }
    }
  } else {
    ## no isotope ROIs
    isotopeROIs.matrix <- matrix(nrow = 0, ncol = 8)
    colnames(isotopeROIs.matrix) <- c("mz", "mzmin", "mzmax", "scmin", "scmax", "length", "intensity", "scale")
  }
  ###############################################################################################
  ## adduct ROIs
  if(addNewAdductROIs){
    ## considered adduct distances
    ## reference: Huang N.; Siegel M.M.1; Kruppa G.H.; Laukien F.H.; J Am Soc Mass Spectrom 1999, 10, 1166–1173; Automation of a Fourier transform ion cyclotron resonance mass spectrometer for acquisition, analysis, and e-mailing of high-resolution exact-mass electrospray ionization mass spectral data
    ## see also for contaminants: Interferences and contaminants encountered in modern mass spectrometry (Bernd O. Keller, Jie Sui, Alex B. Young and Randy M. Whittal, ANALYTICA CHIMICA ACTA, 627 (1): 71-81)
    
    mH  <-  1.0078250322
    mNa <- 22.98976928
    mK  <- 38.96370649
    mC  <- 12
    mN  <- 14.003074004
    mO  <- 15.994914620
    mS  <- 31.972071174
    mCl <- 34.9688527
    mBr <- 78.918338
    mF  <- 18.998403163
    mDMSO    <- mC*2+mH*6+mS+mO     # dimethylsulfoxid
    mACN     <- mC*2+mH*3+mN        # acetonitril
    mIsoProp <- mC*3+mH*8+mO        # isopropanol
    mNH4     <- mN+mH*4             # ammonium
    mCH3OH   <- mC+mH*3+mO+mH       # methanol
    mH2O     <- mH*2+mO             # water
    mFA      <- mC+mH*2+mO*2        # formic acid
    mHAc     <- mC+mH*3+mC+mO+mO+mH # acetic acid
    mTFA     <- mC+mF*3+mC+mO+mO+mH # trifluoroacetic acid
    
    switch(polarity,
           "positive"={
             adductPopulationMz <- unlist(c(
               ## [M+H]+ to [M+H]+  (Reference)
               function(mass){ mass-mH+mNH4 },               ## [M+H]+ to [M+NH4]+
               function(mass){ mass-mH+mNa },                ## [M+H]+ to [M+Na]+
               function(mass){ mass+mCH3OH },                ## [M+H]+ to [M+CH3OH+H]+
               function(mass){ mass-mH+mK },                 ## [M+H]+ to [M+K]+
               function(mass){ mass+mACN },                  ## [M+H]+ to [M+ACN+H]+
               function(mass){ mass-2*mH+2*mNa },            ## [M+H]+ to [M+2Na-H]+
               function(mass){ mass+mIsoProp },              ## [M+H]+ to [M+IsoProp+H]+
               function(mass){ mass-mH+mACN+mNa },           ## [M+H]+ to [M+ACN+Na]+
               function(mass){ mass-2*mH+2*mK },             ## [M+H]+ to [M+2K-H]+
               function(mass){ mass+mDMSO },                 ## [M+H]+ to [M+DMSO+H]+
               function(mass){ mass+2*mACN },                ## [M+H]+ to [M+2*ACN+H]+
               function(mass){ mass+mIsoProp+mNa },          ## [M+H]+ to [M+IsoProp+Na+H]+ TODO double-charged?
               function(mass){ (mass-mH)*2+mH },             ## [M+H]+ to [2M+H]+
               function(mass){ (mass-mH)*2+mNH4 },           ## [M+H]+ to [2M+NH4]+
               function(mass){ (mass-mH)*2+mNa },            ## [M+H]+ to [2M+Na]+
               function(mass){ (mass-mH)*2+mK },             ## [M+H]+ to [2M+K]+
               function(mass){ (mass-mH)*2+mACN+mH },        ## [M+H]+ to [2M+ACN+H]+
               function(mass){ (mass-mH)*2+mACN+mNa },       ## [M+H]+ to [2M+ACN+Na]+
               function(mass){((mass-mH)*2+3*mH2O+2*mH)/2 }, ## [M+H]+ to [2M+3*H2O+2*H]2+
               function(mass){ (mass+mH)/2 },                ## [M+H]+ to [M+2*H]2+
               function(mass){ (mass+mNH4)/2 },              ## [M+H]+ to [M+H+NH4]2+
               function(mass){ (mass+mNa)/2 },               ## [M+H]+ to [M+H+Na]2+
               function(mass){ (mass+mK)/2 },                ## [M+H]+ to [M+H+K]2+
               function(mass){ (mass+mACN+mH)/2 },           ## [M+H]+ to [M+ACN+2*H]2+
               function(mass){ (mass-mH+2*mNa)/2 },          ## [M+H]+ to [M+2*Na]2+
               function(mass){ (mass+2*mACN+mH)/2 },         ## [M+H]+ to [M+2*ACN+2*H]2+
               function(mass){ (mass+3*mACN+mH)/2 },         ## [M+H]+ to [M+3*ACN+2*H]2+
               function(mass){ (mass+2*mH)/3 },              ## [M+H]+ to [M+3*H]3+
               function(mass){ (mass+mH+mNa)/3 },            ## [M+H]+ to [M+2*H+Na]3+
               function(mass){ (mass+2*mNa)/3 },             ## [M+H]+ to [M+H+2*Na]3+
               function(mass){ (mass-mH+3*mNa)/3 }           ## [M+H]+ to [M+3*Na]3+
             ))
           },
           "negative"={
             adductPopulationMz <- unlist(c(
               ## [M-H]+ to [M-H]+  (Reference)
               function(mass){ mass-mH2O },             ## [M-H]+ to [M-H2O-H]+
               function(mass){ mass-mH+mNa },           ## [M-H]+ to [M+Na-2*H]+
               function(mass){ mass+mH+mCl },           ## [M-H]+ to [M+Cl]+
               function(mass){ mass-mH+mK },            ## [M-H]+ to [M+K-2*H]+
               function(mass){ mass+mFA },              ## [M-H]+ to [M+FA-H]+
               function(mass){ mass+mHAc },             ## [M-H]+ to [M+HAc-H]+
               function(mass){ mass+mH+mBr },           ## [M-H]+ to [M+Br]+
               function(mass){ mass+mTFA },             ## [M-H]+ to [M+TFA-H]+
               function(mass){ (mass+mH)*2-mH },        ## [M-H]+ to [2M-H]+
               function(mass){ (mass+mH)*2+mFA-mH },    ## [M-H]+ to [2M+FA-H]+
               function(mass){ (mass+mH)*2+mHAc-mH },   ## [M-H]+ to [2M+HAc-H]+
               function(mass){ (mass+mH)*3-mH },        ## [M-H]+ to [3M-H]+
               function(mass){ (mass-mH)/2 },           ## [M-H]+ to [M-2*H]2+
               function(mass){ (mass-2*mH)/3 }          ## [M-H]+ to [M-3*H]3+
             ))
           },
           "unknown"={
             warning(paste("Unknown polarity! No adduct ROIs have been added.", sep = ""))
           },
           stop(paste("Unknown polarity (", polarity, ")!", sep = ""))
    )
    
    numberOfAdductROIs <- length(ROI.list) * length(adductPopulationMz)
    adductROIs.matrix <- matrix(nrow = numberOfAdductROIs, ncol = 8)
    colnames(adductROIs.matrix) <- c("mz", "mzmin", "mzmax", "scmin", "scmax", "length", "intensity", "scale")
    
    for(roiIdx in 1:(length(ROI.list))){
      for(mzIdx in 1:length(adductPopulationMz)){
        ## create new ROI!
        mzDifference <- adductPopulationMz[[mzIdx]](ROI.list[[roiIdx]]$mz)
        idx <- (roiIdx - 1) * length(adductPopulationMz) + mzIdx
        if(ROI.list[[roiIdx]]$mzmin + mzDifference > 0){
          adductROIs.matrix[idx, ] <- c(
            ROI.list[[roiIdx]]$mz + mzDifference,## XXX not used!
            ROI.list[[roiIdx]]$mzmin + mzDifference,
            ROI.list[[roiIdx]]$mzmax + mzDifference,
            ROI.list[[roiIdx]]$scmin,
            ROI.list[[roiIdx]]$scmax,
            ROI.list[[roiIdx]]$length,## XXX not used!
            -1,  #ROI.list[[roiIdx]]$intensity ## XXX not used!
            ROI.list[[roiIdx]]$scale
          )
        }
      }
    }
  } else {
    ## no adduct ROIs
    adductROIs.matrix <- matrix(nrow = 0, ncol = 8)
    colnames(adductROIs.matrix) <- c("mz", "mzmin", "mzmax", "scmin", "scmax", "length", "intensity", "scale")
  }
  
  numberOfAdditionalIsotopeROIsUnfiltered <- nrow(isotopeROIs.matrix)
  numberOfAdditionalAdductROIsUnfiltered  <- nrow(adductROIs.matrix )
  numberOfAdditionalROIsUnfiltered        <- numberOfAdditionalIsotopeROIsUnfiltered + numberOfAdditionalAdductROIsUnfiltered
  newROI.matrixUnfiltered <- rbind(isotopeROIs.matrix, adductROIs.matrix)
  
  ###############################################################################################
  ## filter out m/z's out of range and without sufficient intensity
  intensityThreshold <- 10
  
  if(addNewIsotopeROIs) isotopeROIs.matrix <- removeROIsOutOfRange(object, isotopeROIs.matrix)
  if(addNewAdductROIs)  adductROIs.matrix  <- removeROIsOutOfRange(object, adductROIs.matrix)
  if(addNewIsotopeROIs) isotopeROIs.matrix <- removeROIsWithoutSignal(object, isotopeROIs.matrix, intensityThreshold)
  if(addNewAdductROIs)  adductROIs.matrix  <- removeROIsWithoutSignal(object, adductROIs.matrix, intensityThreshold)
  
  numberOfAdditionalIsotopeROIs <- nrow(isotopeROIs.matrix)
  numberOfAdditionalAdductROIs  <- nrow(adductROIs.matrix )
  numberOfAdditionalROIs        <- numberOfAdditionalIsotopeROIs + numberOfAdditionalAdductROIs
  
  ###############################################################################################
  ## box
  newROI.matrix <- rbind(isotopeROIs.matrix, adductROIs.matrix)
  
  resultObj <- list()
  ## unfiltered
  resultObj$newROI.matrixUnfiltered <- newROI.matrixUnfiltered
  resultObj$numberOfAdditionalROIsUnfiltered        <- numberOfAdditionalROIsUnfiltered
  resultObj$numberOfAdditionalIsotopeROIsUnfiltered <- numberOfAdditionalIsotopeROIsUnfiltered
  resultObj$numberOfAdditionalAdductROIsUnfiltered  <- numberOfAdditionalAdductROIsUnfiltered
  ## filtered
  resultObj$newROI.matrix <- newROI.matrix
  resultObj$numberOfAdditionalROIs        <- numberOfAdditionalROIs
  resultObj$numberOfAdditionalIsotopeROIs <- numberOfAdditionalIsotopeROIs
  resultObj$numberOfAdditionalAdductROIs  <- numberOfAdditionalAdductROIs
  
  return(resultObj)
}

############################################################
## findPeaks.MSW
setMethod("findPeaks.MSW", "xcmsRaw", function(object, snthresh=3, verbose.columns = FALSE, ...)
{
    ## Should consider to put MassSpecWavelet into Imports instead of Suggests.
          ## require(MassSpecWavelet) || stop("Couldn't load MassSpecWavelet")

          ## MassSpecWavelet Calls
          peakInfo <- peakDetectionCWT(object@env$intensity, SNR.Th=snthresh, ...)
          majorPeakInfo <- peakInfo$majorPeakInfo

          sumIntos <- function(into, inpos, scale){
              scale=floor(scale)
              sum(into[(inpos-scale):(inpos+scale)])
          }

          maxIntos <- function(into, inpos, scale){
              scale=floor(scale)
              max(into[(inpos-scale):(inpos+scale)])
          }


          betterPeakInfo <- tuneInPeakInfo(object@env$intensity,
                                           majorPeakInfo)

          peakIndex <- betterPeakInfo$peakIndex

          ## sum and max of raw values, sum and max of filter-response
          rints<-NA;fints<-NA
          maxRints<-NA;maxFints<-NA

          for (a in 1:length(peakIndex)){
              rints[a]<-    sumIntos(object@env$intensity,peakIndex[a],
                                     betterPeakInfo$peakScale[a])
              maxRints[a]<- maxIntos(object@env$intensity,peakIndex[a],
                                     betterPeakInfo$peakScale[a])
          }
          ## filter-response is not summed here, the maxF-value is the one
          ## which was "xcmsRaw$into" earlier


          ## Assemble result

          basenames <- c("mz","mzmin","mzmax","rt","rtmin","rtmax",
                         "into","maxo","sn","intf","maxf")

          peaklist <- matrix(-1, nrow = length(peakIndex), ncol = length(basenames))

          colnames(peaklist) <- c(basenames)

          peaklist[,"mz"] <- object@env$mz[peakIndex]
          peaklist[,"mzmin"] <- object@env$mz[(peakIndex-betterPeakInfo$peakScale)]
          peaklist[,"mzmax"] <- object@env$mz[(peakIndex+betterPeakInfo$peakScale)]

          peaklist[,"rt"]    <- rep(-1, length(peakIndex))
          peaklist[,"rtmin"] <- rep(-1, length(peakIndex))
          peaklist[,"rtmax"] <- rep(-1, length(peakIndex))

          peaklist[,"into"] <- rints ## sum of raw-intensities
          peaklist[,"maxo"] <- maxRints
          peaklist[,"intf"] <- rep(NA,length(peakIndex))
          peaklist[,"maxf"] <- betterPeakInfo$peakValue

          peaklist[,"sn"]   <- betterPeakInfo$peakSNR

          cat('\n')

          ## Filter additional (verbose) columns
          if (!verbose.columns)
              peaklist <- peaklist[,basenames,drop=FALSE]

          invisible(new("xcmsPeaks", peaklist))
      }
      )

############################################################
## findPeaks.MS1
setMethod("findPeaks.MS1", "xcmsRaw", function(object)
      {
          if (is.null(object@msnLevel)) {
              stop("xcmsRaw contains no MS2 spectra\n")
          }

          ## Select all MS2 scans, they have an MS1 parent defined
          peakIndex <- object@msnLevel == 2

          ## (empty) return object
          basenames <- c("mz","mzmin","mzmax",
                         "rt","rtmin","rtmax",
                         "into","maxo","sn")
          peaklist <- matrix(-1, nrow = length(which(peakIndex)),
                             ncol = length(basenames))
          colnames(peaklist) <- c(basenames)

          ## Assemble result

          peaklist[,"mz"] <- object@msnPrecursorMz[peakIndex]
          peaklist[,"mzmin"] <- object@msnPrecursorMz[peakIndex]
          peaklist[,"mzmax"] <- object@msnPrecursorMz[peakIndex]


          if (any(!is.na(object@msnPrecursorScan))&&any(object@msnPrecursorScan!=0)) {
              peaklist[,"rt"] <- peaklist[,"rtmin"] <- peaklist[,"rtmax"] <- object@scantime[object@msnPrecursorScan[peakIndex]]
          } else {
              ## This happened with ReAdW mzxml
              cat("MS2 spectra without precursorScan references, using estimation")
              ## which object@Scantime are the biggest wich are smaller than the current object@msnRt[peaklist]?
              ms1Rts<-rep(0,length(which(peakIndex)))
              i<-1
              for (a in which(peakIndex)){
                  ms1Rts[i] <- object@scantime[max(which(object@scantime<object@msnRt[a]))]
                  i<-i+1
              }
              peaklist[,"rt"] <-  ms1Rts
              peaklist[,"rtmin"] <-  ms1Rts
              peaklist[,"rtmax"] <- ms1Rts
          }

          if (any(object@msnPrecursorIntensity!=0)) {
              peaklist[,"into"] <- peaklist[,"maxo"] <- peaklist[,"sn"] <- object@msnPrecursorIntensity[peakIndex]
          } else {
              ## This happened with Agilent MzDataExport 1.0.98.2
              warning("MS2 spectra without precursorIntensity, setting to zero")
              peaklist[,"into"] <- peaklist[,"maxo"] <- peaklist[,"sn"] <- 0
          }

          cat('\n')

          invisible(new("xcmsPeaks", peaklist))
      })

############################################################
## findPeaks
setMethod("findPeaks", "xcmsRaw", function(object, method=getOption("BioC")$xcms$findPeaks.method,
                                           ...) {

    method <- match.arg(method, getOption("BioC")$xcms$findPeaks.methods)
    if (is.na(method))
        stop("unknown method : ", method)
    method <- paste("findPeaks", method, sep=".")
    invisible(do.call(method, list(object, ...)))
})

############################################################
## getPeaks
setMethod("getPeaks", "xcmsRaw", function(object, peakrange, step = 0.1) {

    profFun <- match.profFun(object)
    if (all(c("mzmin","mzmax","rtmin","rtmax") %in% colnames(peakrange)))
        peakrange <- peakrange[,c("mzmin","mzmax","rtmin","rtmax"),drop=FALSE]
    stime <- object@scantime

### Create EIC buffer
    mrange <- range(peakrange[,1:2])
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
        ymat <- buf[bufidx[imz[1]:imz[2]],iret[1]:iret[2],drop=FALSE]
        ymax <- colMax(ymat)
        iymax <- which.max(ymax)

        pwid <- diff(stime[iret])/diff(iret)

        rosm <- rowSums(ymat)
        limz <- length(imz[1]:imz[2])
        if (length(rosm) != limz) { ## that happens for some reason
            warning("weighted.mean  : x and w must have the same length \n")
            rosm <- rep(1, limz)  ## fallback to mean
        }
        rmat[i,1] <- weighted.mean(mass[imz[1]:imz[2]], rosm)
        if (is.nan(rmat[i,1]) || is.na(rmat[i,1])) ##  R2.11 :  weighted.mean()  results in NA (not NaN) for zero weights
            rmat[i,1] <- mean(peakrange[i,1:2])

        rmat[i,2:3] <- peakrange[i,1:2]
        rmat[i,4] <- stime[iret[1]:iret[2]][iymax]
        rmat[i,5:6] <- peakrange[i,3:4]

        if (peakrange[i,3] <  stime[1] || peakrange[i,4] > stime[length(stime)] || is.nan(pwid)) {
            warning("getPeaks: Peak  m/z:",peakrange[i,1],"-",peakrange[i,2], ",  RT:",peakrange[i,3],"-",peakrange[i,4],
                    "is out of retention time range for this sample (",object@filepath,"), using zero intensity value.\n")
            rmat[i,7:8] <- 0
        } else {
            rmat[i,7] <- pwid*sum(ymax)
            rmat[i,8] <- ymax[iymax]
        }
    }
    invisible(rmat)
})

############################################################
## plotPeaks
setMethod("plotPeaks", "xcmsRaw", function(object, peaks, figs, width = 200) {

    if (missing(figs)) {
        figs <- c(floor(sqrt(nrow(peaks))), ceiling(sqrt(nrow(peaks))))
        if (prod(figs) < nrow(peaks))
            figs <- rep(ceiling(sqrt(nrow(peaks))), 2)
    }

    mzi <- round((peaks[,c("mzmin","mzmax")]-object@mzrange[1])/profStep(object) + 1)

    screens <- split.screen(figs)
    on.exit(close.screen(all.screens = TRUE))

    for (i in seq(length = min(nrow(peaks), prod(figs)))) {
        screen(screens[i])
        par(cex.main = 1, font.main = 1, mar = c(0, 0, 1, 0) + 0.1)
        xlim <- c(-width/2, width/2) + peaks[i,"rt"]
        ##main <- paste(peaks[i,"i"], " ", round(peaks[i,"mz"]),
        main <- paste(round(peaks[i,"mz"]),
                      " ", round(peaks[i,"rt"]), sep = "")
        plot(object@scantime, colMax(object@env$profile[mzi[i,],,drop=FALSE]),
             type = "l", xlim = xlim, ylim = c(0, peaks[i,"maxo"]), main = main,
             xlab = "", ylab = "", xaxt = "n", yaxt = "n")
        abline(v = peaks[i,c("rtmin","rtmax")], col = "grey")
    }
})

############################################################
## getEIC
setMethod("getEIC", "xcmsRaw", function(object, mzrange, rtrange = NULL, step = 0.1) {
              FUN <- getOption("BioC")$xcms$getEIC.method
              if(FUN == "getEICOld"){
                  return(getEICOld(object=object, mzrange=mzrange, rtrange=rtrange, step=step))
              }else if(FUN == "getEICNew"){
                  return(getEICNew(object=object, mzrange=mzrange, rtrange=rtrange, step=step))
              }else{
                  stop("Method ", FUN, " not known! getEIC.method should be either getEICOld or getEICnew!")
              }
          })

############################################################
## rawMat
setMethod("rawMat", "xcmsRaw", function(object,
                                        mzrange = numeric(),
                                        rtrange = numeric(),
                                        scanrange = numeric(),
                                        log=FALSE) {

    if (length(rtrange) >= 2) {
        rtrange <- range(rtrange)
        scanidx <- (object@scantime >= rtrange[1]) & (object@scantime <= rtrange[2])
        scanrange <- c(match(TRUE, scanidx),
                       length(scanidx) - match(TRUE, rev(scanidx)))
    }
    else if (length(scanrange) < 2)
        scanrange <- c(1, length(object@scantime))
    else scanrange <- range(scanrange)
    startidx <- object@scanindex[scanrange[1]] + 1
    endidx <- length(object@env$mz)
    if (scanrange[2] < length(object@scanindex))
        endidx <- object@scanindex[scanrange[2] + 1]
    ##scans <- integer(endidx - startidx + 1)
    scans <- rep(scanrange[1]:scanrange[2],
                 diff(c(object@scanindex[scanrange[1]:scanrange[2]], endidx)))
    ##for (i in scanrange[1]:scanrange[2]) {
    ##    idx <- (object@scanindex[i] + 1):min(object@scanindex[i +
    ##        1], length(object@env$mz), na.rm = TRUE)
    ##    scans[idx - startidx + 1] <- i
    ##}
    rtrange <- c(object@scantime[scanrange[1]], object@scantime[scanrange[2]])
    masses <- object@env$mz[startidx:endidx]
    int <- object@env$intensity[startidx:endidx]
    massidx <- 1:length(masses)
    if (length(mzrange) >= 2) {
        mzrange <- range(mzrange)
        massidx <- massidx[(masses >= mzrange[1]) & (masses <= mzrange[2])]
    }
    else mzrange <- range(masses)

    y <- int[massidx]
    if (log && (length(y)>0))
        y <- log(y + max(1 - min(y), 0))

    cbind(time = object@scantime[scans[massidx]], mz = masses[massidx],
          intensity = y)
})

############################################################
## plotRaw
setMethod("plotRaw", "xcmsRaw", function(object,
                                         mzrange = numeric(),
                                         rtrange = numeric(),
                                         scanrange = numeric(),
                                         log=FALSE,title='Raw Data' ) {

    raw <- rawMat(object, mzrange, rtrange, scanrange, log)

    if (nrow(raw) > 0) {
        y <- raw[,"intensity"]
        ylim <- range(y)
        y <- y/ylim[2]
        colorlut <- terrain.colors(16)
        col <- colorlut[y*15+1]
        plot(cbind(raw[,"time"], raw[,"mz"]), pch=20, cex=.5,
             main = title, xlab="Seconds", ylab="m/z", col=col,
             xlim=range(raw[,"time"]), ylim=range(raw[,"mz"]))
    } else {
        if (length(rtrange) >= 2) {
            rtrange <- range(rtrange)
            scanidx <- (object@scantime >= rtrange[1]) & (object@scantime <= rtrange[2])
            scanrange <- c(match(TRUE, scanidx),
                           length(scanidx) - match(TRUE, rev(scanidx)))
        } else if (length(scanrange) < 2)
            scanrange <- c(1, length(object@scantime)) else
        scanrange <- range(scanrange)
        plot(c(NA,NA), main = title, xlab="Seconds", ylab="m/z",
             xlim=c(object@scantime[scanrange[1]],object@scantime[scanrange[2]]), ylim=mzrange)
    }

    invisible(raw)
})

############################################################
## profMz
setMethod("profMz", "xcmsRaw", function(object) {
    object@mzrange[1]+profStep(object)*(0:(dim(object@env$profile)[1]-1))
})

############################################################
## profMethods
setMethod("profMethod", "xcmsRaw", function(object) {
    object@profmethod
})
setReplaceMethod("profMethod", "xcmsRaw", function(object, value) {

    if (! (value %in% names(.profFunctions)))
        stop("Invalid profile method")
    object@profmethod <- value
    profStep(object) <- profStep(object)
    object
})

############################################################
## profStep
setMethod("profStep", "xcmsRaw", function(object) {
    if (is.null(object@env$profile))
        0
    else
        diff(object@mzrange)/(nrow(object@env$profile)-1)
})
setReplaceMethod("profStep", "xcmsRaw", function(object, value) {

    if ("profile" %in% ls(object@env))
        rm("profile", envir = object@env)
    if (!value)
        return(object)
    if (length(object@env$mz)==0) {
        warning("MS1 scans empty. Skipping profile matrix calculation.")
        return(object)
    }
    minmass <- round(min(object@env$mz)/value)*value
    maxmass <- round(max(object@env$mz)/value)*value

    num <- (maxmass - minmass)/value + 1
    profFun <- match.profFun(object)
    object@env$profile <- profFun(object@env$mz, object@env$intensity,
                                  object@scanindex, num, minmass, maxmass,
                                  FALSE, object@profparam)
    object@mzrange <- c(minmass, maxmass)
    return(object)
})

############################################################
## profStepPad
setReplaceMethod("profStepPad", "xcmsRaw", function(object, value) {

    if ("profile" %in% ls(object@env))
        rm("profile", envir = object@env)
    if (!value)
        return(object)

    if (length(object@env$mz)==0) {
        warning("MS1 scans empty. Skipping profile matrix calculation.")
        return(object)
    }
    mzrange <- range(object@env$mz)
    ## calculate the profile matrix with whole-number limits
    minmass <- floor(mzrange[1])
    maxmass <- ceiling(mzrange[2])

    num <- (maxmass - minmass)/value + 1
    profFun <- match.profFun(object)
    object@env$profile <- profFun(object@env$mz, object@env$intensity,
                                  object@scanindex, num, minmass, maxmass,
                                  FALSE, object@profparam)
    object@mzrange <- c(minmass, maxmass)
    return(object)
})

############################################################
## profMedFilt
setMethod("profMedFilt", "xcmsRaw", function(object, massrad = 0, scanrad = 0) {

    contdim <- dim(object@env$profile)
    object@env$profile <- medianFilter(object@env$profile, massrad, scanrad)
})

############################################################
## profRange
setMethod("profRange", "xcmsRaw", function(object,
                                           mzrange = numeric(),
                                           rtrange = numeric(),
                                           scanrange = numeric(), ...) {

    if (length(object@env$profile)) {
        contmass <- profMz(object)
        if (length(mzrange) == 0) {
            mzrange <- c(min(contmass), max(contmass))
        } else if (length(mzrange) == 1) {
            closemass <- contmass[which.min(abs(contmass-mzrange))]
            mzrange <- c(closemass, closemass)
        } else if (length(mzrange) > 2) {
            mzrange <- c(min(mzrange), max(mzrange))
        }
        massidx <- which((contmass >= mzrange[1]) & (contmass <= mzrange[2]))
    } else {
        if (length(mzrange) == 0) {
            mzrange <- range(object@env$mz)
        } else {
            mzrange <- c(min(mzrange), max(mzrange))
        }
        massidx <- integer()
    }
    if (mzrange[1] == mzrange[2])
        masslab <- paste(mzrange[1], "m/z")
    else
        masslab <- paste(mzrange[1], "-", mzrange[2], " m/z", sep="")


    if (length(rtrange) == 0) {
        if (length(scanrange) == 0)
            scanrange <- c(1, length(object@scanindex))
        else if (length(scanrange) == 1)
            scanrange <- c(scanrange, scanrange)
        else if (length(scanrange) > 2)
            scanrange <- c(max(1, min(scanrange)), min(max(scanrange), length(object@scantime)))
        rtrange <- c(object@scantime[scanrange[1]], object@scantime[scanrange[2]])
    } else if (length(rtrange) == 1) {
        closetime <- object@scantime[which.min(abs(object@scantime-rtrange))]
        rtrange <- c(closetime, closetime)
    } else if (length(rtrange) > 2) {
        rtrange <- c(min(rtrange), max(rtrange))
    }

    if (rtrange[1] == rtrange[2])
        timelab <- paste(round(rtrange[1],1), "seconds")
    else
        timelab <- paste(round(rtrange[1],1), "-", round(rtrange[2],1), " seconds", sep="")


    if (length(scanrange) == 0) {
        scanidx <- which((object@scantime >= rtrange[1]) & (object@scantime <= rtrange[2]))
        scanrange <- c(min(scanidx), max(scanidx))
    } else {
        scanidx <- scanrange[1]:scanrange[2]
    }

    if (scanrange[1] == scanrange[2])
        scanlab <- paste("scan", scanrange[1])
    else
        scanlab <- paste("scans ", scanrange[1], "-", scanrange[2], sep="")

    list(mzrange = mzrange, masslab = masslab, massidx = massidx,
         scanrange = scanrange, scanlab = scanlab, scanidx = scanidx,
         rtrange = rtrange, timelab = timelab)
})

############################################################
## rawEIC
setMethod("rawEIC", "xcmsRaw", function(object,
                                        mzrange = numeric(),
                                        rtrange = numeric(),
                                        scanrange = numeric())  {

    if (length(rtrange) >= 2) {
        rtrange <- range(rtrange)
        scanidx <- (object@scantime >= rtrange[1]) & (object@scantime <= rtrange[2])
        scanrange <- c(match(TRUE, scanidx), length(scanidx) - match(TRUE, rev(scanidx)))
    }  else if (length(scanrange) < 2)
        scanrange <- c(1, length(object@scantime))
    else
        scanrange <- range(scanrange)

    scanrange[1] <- max(1,scanrange[1])
    scanrange[2] <- min(length(object@scantime),scanrange[2])

    if (!is.double(object@env$mz))  object@env$mz <- as.double(object@env$mz)
    if (!is.double(object@env$intensity)) object@env$intensity <- as.double(object@env$intensity)
    if (!is.integer(object@scanindex)) object@scanindex <- as.integer(object@scanindex)

    .Call("getEIC",object@env$mz,object@env$intensity,object@scanindex,as.double(mzrange),as.integer(scanrange),as.integer(length(object@scantime)), PACKAGE ='xcms' )
})

############################################################
## plotEIC
setMethod("plotEIC", "xcmsRaw", function(object,
                                         mzrange = numeric(),
                                         rtrange = numeric(),
                                         scanrange = numeric(),
                                         type="l", add=FALSE, ...)  {
              if(length(mzrange)==0)
                  mzrange <- range(object@env$mz)
              if(length(rtrange)==0)
                  rtrange <- range(object@scantime)
    EIC <-  rawEIC(object,mzrange=mzrange, rtrange=rtrange, scanrange=scanrange)
    points <- cbind(object@scantime[EIC$scan], EIC$intensity)
    if(add){
        points(points, type=type, ...)
    }else{
        plot(points, type=type, main=paste("Extracted Ion Chromatogram  m/z  ",mzrange[1]," - ",mzrange[2],sep=""), xlab="Seconds",
             ylab="Intensity", xlim=rtrange, ...)
    }
    invisible(points)
})

############################################################
## rawMZ
setMethod("rawMZ", "xcmsRaw", function(object,
                                       mzrange = numeric(),
                                       rtrange = numeric(),
                                       scanrange = numeric())  {

    if (length(rtrange) >= 2) {
        rtrange <- range(rtrange)
        scanidx <- (object@scantime >= rtrange[1]) & (object@scantime <= rtrange[2])
        scanrange <- c(match(TRUE, scanidx), length(scanidx) - match(TRUE, rev(scanidx)))
    }  else if (length(scanrange) < 2)
        scanrange <- c(1, length(object@scantime))
    else
        scanrange <- range(scanrange)

    scanrange[1] <- max(1,scanrange[1])
    scanrange[2] <- min(length(object@scantime),scanrange[2])

    if (!is.double(object@env$mz))  object@env$mz <- as.double(object@env$mz)
    if (!is.double(object@env$intensity)) object@env$intensity <- as.double(object@env$intensity)
    if (!is.integer(object@scanindex)) object@scanindex <- as.integer(object@scanindex)

    .Call("getMZ",object@env$mz,object@env$intensity,object@scanindex,as.double(mzrange),as.integer(scanrange),as.integer(length(object@scantime)), PACKAGE ='xcms' )
})

############################################################
## findmzROI
setMethod("findmzROI", "xcmsRaw", function(object, mzrange=c(0.0,0.0), scanrange=c(1,length(object@scantime)),dev, minCentroids, prefilter=c(0,0), noise=0){

    scanrange[1] <- max(1,scanrange[1])
    scanrange[2] <- min(length(object@scantime),scanrange[2])

    ## mzrange not implemented yet
    if (!is.double(object@env$mz))  object@env$mz <- as.double(object@env$mz)
    if (!is.double(object@env$intensity)) object@env$intensity <- as.double(object@env$intensity)
    if (!is.integer(object@scanindex)) object@scanindex <- as.integer(object@scanindex)

    ## .Call("findmzROI", object@env$mz,object@env$intensity,object@scanindex, as.double(mzrange),
    ##       as.integer(scanrange), as.integer(length(object@scantime)),
    ##       as.double(dev), as.integer(minCentroids), as.integer(prefilter), as.integer(noise), PACKAGE ='xcms' )

    ROIs <- NULL
    withRestarts(
        tryCatch({
            ROIs <- .Call("findmzROI",
                          object@env$mz,
                          object@env$intensity,
                          object@scanindex,
                          as.double(mzrange),
                          as.integer(scanrange),
                          as.integer(length(object@scantime)),
                          as.double(dev),
                          as.integer(minCentroids),
                          as.integer(prefilter),
                          as.integer(noise),
                          PACKAGE ='xcms' )
        },
        error=function(e) {if (grepl("m/z sort assumption violated !", e$message))
                           {invokeRestart("fixSort")} else {simpleError(e)}}),
        fixSort = function() {
            ## Check and fix "m/z sort assumption violated !"
            for(i in 1:length(object@scanindex)){
                scan <- getScan(object, scan=i)
                if(is.unsorted(scan[,"mz"])){
                    message("Scan ", i, " is unsorted. Fixing.")
                    o <- order(scan[,"mz"])
                    start <- object@scanindex[i] + 1
                    end <- start+nrow(scan) - 1
                    object@env$mz[start:end] <- scan[o, "mz"]
                    object@env$intensity[start:end] <- scan[o, "intensity"]
                }
            }

            ## Re-run now with fixed m/z order
            ROIs <<- .Call("findmzROI",
                           object@env$mz,
                           object@env$intensity,
                           object@scanindex,
                           as.double(mzrange),
                           as.integer(scanrange),
                           as.integer(length(object@scantime)),
                           as.double(dev),
                           as.integer(minCentroids),
                           as.integer(prefilter),
                           as.integer(noise),
                           PACKAGE ='xcms' )
        }
    )
    return(ROIs)
})

############################################################
## findKalmanROI
setMethod("findKalmanROI", "xcmsRaw", function(object, mzrange=c(0.0,0.0),
                                               scanrange=c(1,length(object@scantime)), minIntensity,
                                               minCentroids, consecMissedLim, criticalVal, ppm,  segs, scanBack){

    scanrange[1] <- max(1,scanrange[1])
    scanrange[2] <- min(length(object@scantime),scanrange[2])

    ## var type checking
    if (!is.double(object@env$mz))  object@env$mz <- as.double(object@env$mz)
    if (!is.double(object@env$intensity)) object@env$intensity <- as.double(object@env$intensity)
    if (!is.integer(object@scanindex)) object@scanindex <- as.integer(object@scanindex)
    if (!is.double(object@scantime)) object@scantime <- as.double(object@scantime)


    .Call("massifquant", object@env$mz,object@env$intensity,object@scanindex, object@scantime,
          as.double(mzrange), as.integer(scanrange), as.integer(length(object@scantime)),
          as.double(minIntensity),as.integer(minCentroids),as.double(consecMissedLim),
          as.double(ppm), as.double(criticalVal), as.integer(segs), as.integer(scanBack), PACKAGE ='xcms' )
})

############################################################
## findPeaks.massifquant
setMethod("findPeaks.massifquant", "xcmsRaw", function(object, ppm=10, peakwidth=c(20,50), snthresh=10,
                                                       prefilter=c(3,100), mzCenterFun="wMean", integrate=1, mzdiff=-0.001,
                                                       fitgauss=FALSE, scanrange= numeric(), noise=0, ## noise.local=TRUE,
                                                       sleep=0, verbose.columns=FALSE, criticalValue = 1.125, consecMissedLimit = 2,
                                                       unions = 1, checkBack = 0, withWave = 0) {

    cat("\n Massifquant, Copyright (C) 2013 Brigham Young University.");
    cat("\n Massifquant comes with ABSOLUTELY NO WARRANTY. See LICENSE for details.\n");
    flush.console();

    ##keeep this check since massifquant doesn't check internally
    if (!isCentroided(object))
        warning("It looks like this file is in profile mode. massifquant can process only centroid mode data !\n")

    cat("\n Detecting  mass traces at",ppm,"ppm ... \n"); flush.console();
    massifquantROIs = findKalmanROI(object, minIntensity = prefilter[2], minCentroids = peakwidth[1], criticalVal = criticalValue,
    consecMissedLim = consecMissedLimit, segs = unions, scanBack = checkBack,ppm=ppm);

    if (withWave == 1) {
        featlist = findPeaks.centWave(object, ppm, peakwidth, snthresh,
        prefilter, mzCenterFun, integrate, mzdiff, fitgauss,
        scanrange, noise, sleep, verbose.columns, ROI.list= massifquantROIs);
    }
    else {
        basenames <- c("mz","mzmin","mzmax","rtmin","rtmax","rt", "into")
        if (length(massifquantROIs) == 0) {
            cat("\nNo peaks found !\n");
            nopeaks <- new("xcmsPeaks", matrix(nrow=0, ncol=length(basenames)));
            colnames(nopeaks) <- basenames;
            return(invisible(nopeaks));
        }

        p <- t(sapply(massifquantROIs, unlist));
        colnames(p) <- basenames;

        #get the max intensity for each feature
        maxo <- sapply(seq_len(nrow(p)), function(i) {
            raw <- rawMat(object, mzrange = p[i,c("mzmin", "mzmax")],
                          scanrange = p[i,c("rtmin", "rtmax")])
            max(raw[,3])
        })
        p <- cbind(p, maxo)

        #calculate median index
        p[,"rt"] = as.integer(p[,"rtmin"] + ( (p[,"rt"] + 1) / 2 ) - 1);
        #convert from index into actual time
        p[,"rtmin"] = object@scantime[p[,"rtmin"]];
        p[,"rtmax"] = object@scantime[p[,"rtmax"]];
        p[,"rt"] = object@scantime[p[,"rt"]];

        uorder <- order(p[,"into"], decreasing=TRUE);
        pm <- as.matrix(p[,c("mzmin","mzmax","rtmin","rtmax"),drop=FALSE]);

        uindex <- rectUnique(pm,uorder,mzdiff,ydiff = -0.00001) ## allow adjacent peaks;
        featlist <- p[uindex,,drop=FALSE];
        cat("\n",dim(featlist)[1]," Peaks.\n");
        invisible(new("xcmsPeaks", featlist));
    }
    return(invisible(featlist));
})

############################################################
## isCentroided
setMethod("isCentroided", "xcmsRaw", function(object){
    if (length(getScan(object,length(object@scantime) / 2)) >2 ) {
        quantile(diff(getScan(object,length(object@scantime) / 2)[,"mz"]),.25)  > 0.025
    } else {
        TRUE
    }
})

############################################################
## msnparent2ms
setMethod("msnparent2ms", "xcmsRaw", function(object) {
    xr <- new("xcmsRaw")

    xr@env$mz=object@msnPrecursorMz
    xr@env$intensity=object@msnPrecursorIntensity
    xr@scantime = object@msnRt
    xr@scanindex = seq(1,length(object@msnRt))
    xr@acquisitionNum = seq(1,length(object@msnRt))
    xr@mzrange = range(object@msnPrecursorMz)

    xr
})

############################################################
## msn2ms
setMethod("msn2ms", "xcmsRaw", function(object) {

    object@tic <- rep(0, length(object@msnAcquisitionNum)) ##

    object@scantime <- object@msnRt
    object@acquisitionNum <- object@msnAcquisitionNum
    object@scanindex <- object@msnScanindex

    object@env$mz <- object@env$msnMz
    object@env$intensity <- object@env$msnIntensity
    invisible(object)

})

############################################################
## deepCopy
setMethod("deepCopy", "xcmsRaw", function(object) {

    x <- object
    x@env <- new.env(parent=.GlobalEnv)

    for (variable in ls(object@env)) {
        eval(parse(text=paste("x@env$",variable," <- object@env$",variable,sep="")))
    }

    invisible(x)
})

############################################################
## levelplot
## levelplot for xcmsRaw objects; contains code from the image method, but uses the levelplot
## from the lattice package.
setMethod("levelplot", "xcmsRaw", function(x, log=TRUE,
                                           col.regions=colorRampPalette(brewer.pal(9, "YlOrRd"))(256), ...){
    ## some code taken from plotSurf...
    sel <- profRange(x, ...)
    zvals <- x@env$profile[sel$massidx, sel$scanidx]
    if(log){
        zvals <- log(zvals+max(c(-min(zvals), 1)))
    }
    ## y axis is time
    yvals <- x@scantime[sel$scanidx]
    yrange <- range(yvals)
    ## y has to be sequentially increasing!
    if(length(unique(yvals))!=length(yvals))
        yvals <- 1:length(yvals)
    ## x is m/z
    xvals <- profMz(x)[sel$massidx]
    ## that's much slower...
    ## grid <- expand.grid(x=xvals, y=yvals)
    ## grid <- cbind(grid, z=zvals)
    ## grid$z <- zvals
    ######
    ## now i have to match the x and y to the z.
    ## as.numeric of z returns values by column(!), i.e. the first nrow(z) correspond to
    ## the x of 1.
    xvals <- rep(xvals, ncol(zvals))
    yvals <- rep(yvals, each=nrow(zvals))
    zvals <- as.numeric(zvals)
    ## get the file name
    fileNpath <- x@filepath[1]
    fileName <- unlist(strsplit(fileNpath, split=.Platform$file.sep))
    fileName <- fileName[length(fileName)]
    plt <- levelplot(zvals~xvals*yvals,
                     xlab="m/z", ylab="Time",
                     colorkey=list(height=1, width=0.7),
              main=list(fileName, side=1, line=0.5), col.regions=col.regions)
    ## trellis.focus("legend", side="right", clipp.off=TRUE, highlight=FALSE)
    ## grid.text("Int.", 0.2, 0, hjust=0.5, vjust=1)
    ## trellis.unfocus()
    plt
})

############################################################
## profinfo
setMethod("profinfo", "xcmsRaw", function(object) {
              pinfo <- object@profparam
              ## fill with additional method and step.
              pinfo$method <- profMethod(object)
              pinfo$step <- profStep(object)
              return(pinfo)
          })

############################################################
## scanrange
setMethod("scanrange", "xcmsRaw", function(object) {
              if(.hasSlot(object, "scanrange")){
                  srange <- object@scanrange
                  if(length(srange) == 0){
                      ## for compatibility with the xcmsSet and xcmsRaw functions,
                      ## which default the scanrange argument to NULL.
                      return(NULL)
                  }
                  return(srange)
              }else{
                  warning("No slot scanrange available, consider updating the",
                          " object using the 'updateObject' method.")
                  return(NULL)
              }
          })
setReplaceMethod("scanrange", "xcmsRaw", function(object, value) {
                     if(.hasSlot(object, "scanrange")){
                         object@scanrange <- value
                     }else{
                         warning("Object has no slot scanrange, condider updating",
                                 " the object using the 'updateObject' method.")
                     }
                     object
                 })

############################################################
## mslevel
setMethod("mslevel", "xcmsRaw", function(object){
              if(.hasSlot(object, "mslevel")){
                  mlevel <- object@mslevel
                  if(length(mlevel) == 0){
                      ## for compatibility with the xcmsSet and xcmsRaw functions,
                      ## which default the mslevel argument to NULL.
                      return(NULL)
                  }
                  return(object@mslevel)
              }else{
                  warning("No slot mslevel available, consider updating the",
                          " object using the 'updateObject' method.")
                  return(NULL)
              }
          })
setReplaceMethod("mslevel", "xcmsRaw", function(object, value){
                     if(.hasSlot(object, "mslevel")){
                         object@mslevel <- value
                     }else{
                         warning("Object has no slot mslevel, consider updating",
                                 " the object using the 'updateObject' method.")
                     }
                     object
                 })

############################################################
## plotScan
setMethod("plotScan", "xcmsRaw", function(object, scan, mzrange = numeric(),
                                          ident = FALSE)
      {
          if (scan<1 || scan>length(object@scanindex) ) {
              warning("scan out of range")
              return()
          }

          ## handle last spectrum
          if (scan == length(object@scanindex)) {
              followingScanIndex <- length(object@env$mz)
          } else {
              followingScanIndex <- object@scanindex[scan+1]
          }

          ## hendle empty spectra
          if (object@scanindex[scan] == length(object@env$mz) ||
              object@scanindex[scan] == followingScanIndex) {
              warning("empty scan")
              return()
          }

          idx <- (object@scanindex[scan]+1):min(followingScanIndex,
                                                length(object@env$mz), na.rm=TRUE)
          if (length(mzrange) >= 2) {
              mzrange <- range(mzrange)
              idx <- idx[object@env$mz[idx] >= mzrange[1] & object@env$mz[idx] <= mzrange[2]]
          }
          points <- cbind(object@env$mz[idx], object@env$intensity[idx])
          title = paste("Mass Spectrum: ", round(object@scantime[scan], 1),
          " seconds (scan ", scan, ")", sep = "")
          plot(points, type="h", main = title, xlab="m/z", ylab="Intensity")

          if (ident)
              return(identify(points, labels = round(points[,1], 1)))

          invisible(points)
      })

############################################################
## plotSpec
setMethod("plotSpec", "xcmsRaw", function(object, ident = FALSE,
                                          vline = numeric(0), ...) {

    sel <- profRange(object, ...)

    title = paste("Averaged Mass Spectrum: ", sel$timelab, " (",
    sel$scanlab, ")",  sep = "")
    points <- cbind(profMz(object)[sel$massidx],
                    rowMeans(object@env$profile[sel$massidx,sel$scanidx,drop=FALSE]))
    plot(points, type="l", main = title, xlab="m/z", ylab="Intensity")
    if (length(vline))
        abline(v = vline, col = "red")

    if (ident)
        return(identify(points, labels = round(points[,1], 1)))

    invisible(points)
})

############################################################
## plotChrom
setMethod("plotChrom", "xcmsRaw", function(object, base = FALSE, ident = FALSE,
                                           fitgauss = FALSE, vline = numeric(0), ...) {

    sel <- profRange(object, ...)


    if (base) {
        title = paste("Base Peak Chromatogram: ", sel$masslab, sep = "")
        pts <- cbind(object@scantime[sel$scanidx],
                     colMax(object@env$profile[sel$massidx,sel$scanidx,drop=FALSE]))
    }
    else {
        title = paste("Averaged Ion Chromatogram: ", sel$masslab, sep = "")
        pts <- cbind(object@scantime[sel$scanidx],
                     colMeans(object@env$profile[sel$massidx,sel$scanidx,drop=FALSE]))
    }
    plot(pts, type="l", main = title, xlab="Seconds", ylab="Intensity")
    if (length(vline))
        abline(v = vline, col = "red")

    if (fitgauss) {
        fit <- nls(y ~ SSgauss(x, mu, sigma, h), data.frame(x = pts[,1], y = pts[,2]))
        points(pts[,1], fitted(fit), type = "l", col = "red", lwd = 2)
        return(fit)
    }

    if (ident)
        return(identify(pts, labels = round(pts[,1], 1)))

    invisible(pts)
})

############################################################
## image
setMethod("image", "xcmsRaw", function(x, col = rainbow(256), ...) {
    sel <- profRange(x, ...)

    zlim <- log(range(x@env$intensity))

    title <- paste("XC/MS Log Intensity Image (Profile Method: ",
                   x@profmethod, ")", sep = "")
    if (zlim[1] < 0) {
        zlim <- log(exp(zlim)+1)
        image(profMz(x)[sel$massidx], x@scantime[sel$scanidx],
              log(x@env$profile[sel$massidx, sel$scanidx]+1),
              col = col, zlim = zlim, main = title, xlab="m/z", ylab="Seconds")
    } else
        image(profMz(x)[sel$massidx], x@scantime[sel$scanidx],
              log(x@env$profile[sel$massidx, sel$scanidx]),
              col = col, zlim = zlim, main = title, xlab="m/z", ylab="Seconds")
})

############################################################
## plotSurf
setMethod("plotSurf", "xcmsRaw", function(object, log = FALSE,
                                          aspect = c(1, 1, .5), ...) {

    require(rgl) || stop("Couldn't load package rgl")

    sel <- profRange(object, ...)

    y <- object@env$profile[sel$massidx, sel$scanidx]
    if (log)
        y <- log(y+max(1-min(y), 0))
    ylim <- range(y)

    x <- seq(0, aspect[1], length=length(sel$massidx))
    z <- seq(0, aspect[2], length=length(sel$scanidx))
    y <- y/ylim[2]*aspect[3]

    colorlut <- terrain.colors(256)
    col <- colorlut[y/aspect[3]*255+1]

    rgl.clear("shapes")
    rgl.clear("bbox")
    rgl.surface(x, z, y, color = col, shininess = 128)
    rgl.points(0, 0, 0, alpha = 0)

    mztics <- pretty(sel$mzrange, n = 5*aspect[1])
    rttics <- pretty(sel$rtrange, n = 5*aspect[2])
    inttics <- pretty(c(0,ylim), n = 10*aspect[3])
    inttics <- inttics[inttics > 0]

    rgl.bbox(xat = (mztics - sel$mzrange[1])/diff(sel$mzrange)*aspect[1],
             xlab = as.character(mztics),
             yat = inttics/ylim[2]*aspect[3],
             ylab = as.character(inttics),
             zat = (rttics - sel$rtrange[1])/diff(sel$rtrange)*aspect[2],
             zlab = as.character(rttics),
             ylen = 0, alpha=0.5)
})

############################################################
## getMsnScan
setMethod("getMsnScan", "xcmsRaw", function(object, scan, mzrange = numeric()) {

    if (scan < 0)
        scan <- length(object@msnRt) + 1 + scan

    idx <- seq(object@msnScanindex[scan]+1, min(object@msnScanindex[scan+1],
                                                length(object@env$msnMz), na.rm=TRUE))

    if (length(mzrange) >= 2) {
        mzrange <- range(mzrange)
        idx <- idx[object@env$msnMz[idx] >= mzrange[1] & object@env$msnMz[idx] <= mzrange[2]]
    }

    points <- cbind(mz = object@env$msnMz[idx], intensity = object@env$msnIntensity[idx])

    invisible(points)
})

############################################################
## AutoLockMass
setMethod("AutoLockMass", "xcmsRaw", function(object) {
    if(length(grep("xml|mzData|mzXML|mzML", object@filepath, ignore.case=TRUE)) >= 1){
        tempFreq<-diff(which(diff(object@scantime) == 0))-1
        idx <- which(tempFreq != floor(mean(tempFreq))) ## only needed for newer lockmass signal
        if(is.nan(mean(tempFreq)) ){
            dn<-density(diff(object@scantime))
            lockMassScans <- quantile(dn$x, .75) ## hopefully always correct (?)
            inx<-which(diff(object@scantime) >= lockMassScans) ## these seems to be some of the new files
            return(inx)
        }else if(all(tempFreq == mean(tempFreq)) ){
            freqLock<-mean(tempFreq)
        } else if(all(idx == which(tempFreq != floor(mean(tempFreq) )) )){
            ## for the newer mzML and mzXML not sure why the change?
            ## This means that there is only one gap :( ??
            stop("This file is different from the normally seen files and requires special programming\n
                        This functionality has not been implemented yet\n ")
            ## these files seem to come either from newer MS units or/and msconvert ....
        } else {
            freqLock<-mean(tempFreq)
            warning("\nLock mass frequency wasn't detected correctly", immediate.=TRUE)
        }

        if(diff(object@scantime[1:5])[1] == 0 ){
            start<-1
        } else{
            start<-freqLock
        }
        return(makeacqNum(object, freqLock, start))

    } else if(length(grep("cdf", object@filepath, ignore.case=TRUE)) >= 1){
        ## check to see if we have the X02.CDF files around
        ## These files should be the lock mass channel
        file02<-list.files(gsub("01.CDF", "02.CDF", object@filepath), recursive=T)
        if(length(file02)> 0){
            xr<-xcmsRaw(file02)
            lockMass<-sapply(xr@scantime, function(x, object){
                which.min(abs(object@scantime - x))
            }, object)
            return(lockMass)
        } else {
            ## we couldn't find the files so lets try to find them automatically
            hr <- hist(diff(object@scantime), breaks=4, plot=FALSE)
            if(length(hr$counts) > 2){
                idx<-which(hr$counts == 0)
                ## could have something here about which way the plot is ie cor R is - or +
                                        # if(cor(hr$mids, hr$counts) < 0){
                inx<-which(diff(object@scantime) >= hr$mids[(max(idx))])
                                        # } else {
                                        #       inx<-which()
                                        # }
            }else if(length(hr$counts) == 2){
                inx<-which(diff(object@scantime) >= hr$mids[2])
            } else {
                stop("File appears to have been run without lock mass\n ")
            }
            if(length(inx) <= 1){
                warning("\nLock mass frequency wasn't detected", immediate.=TRUE)
                return(0)
            }
            ## above we're looking for scantimes that are much longer than the normal scan times
            tempFreq<-diff(inx)-1
            if(all(tempFreq == median(tempFreq)) ){
                freqLock<-median(tempFreq)
            }else{
                freqLock<-median(tempFreq)
                warning("Lock mass frequency wasn't detected correctly\n", immediate.=TRUE)
            }

            if(inx[1] == 0 || inx[1] == 1){
                start<-1
            }else{
                start<-freqLock
            }
                                        #return(inx)
            return(makeacqNum(object, freqLock, start))
        }
    } else{
        stop("Couldn't detect file type\n")
    }
})

############################################################
## makeacqNum
setMethod("makeacqNum", "xcmsRaw", function(object, freq, start=1) {

    freq<-freq+1 ##nessary for the start at +1 and others since 1st scan is +1

    acqNum<-numeric()
    fo<-seq(from=start, to=length(object@scanindex), by=freq)
    for(i in fo){
        acqNum<-c(acqNum, i,i+1)
    }
    return(acqNum)
})

############################################################
## stitch
setMethod("stitch", "xcmsRaw", function(object, lockMass) {
    if(length(grep("xml|mzData", object@filepath, ignore.case=TRUE)) >= 1){
        type<-stitch.xml
    } else if(length(grep("cdf", object@filepath, ignore.case=TRUE)) >= 1){
        ## lets check to see if lockMass is one scan or two
        if(any(diff(lockMass) == 0)){
            type<-stitch.netCDF.new
        }else {
            type<-stitch.netCDF
        }
    } else{
        stop("Unknown stitch method \n")
    }

    invisible(do.call(type, list(object, lockMass)))
})

############################################################
## stitch.xml
setMethod("stitch.xml", "xcmsRaw", function(object, lockMass) {

    ob<-new("xcmsRaw")
    ob@env$mz<-object@env$mz
    ob@env$intensity<-object@env$intensity
    ob@scanindex<-object@scanindex
    ob@scantime<-object@scantime

    ob@acquisitionNum<-1:length(ob@scanindex)
    ob@filepath<-object@filepath
    ob@mzrange<-range(ob@env$mz)
    ob@profmethod<-object@profmethod
    ob@tic<-object@tic
    ob@profparam<-list()

    arr<-array(dim=c(2,max(diff(ob@scanindex)), length(ob@scanindex)))
    if(lockMass[1] == 1){
        lockMass<-lockMass[3:length(lockMass)]
    }
    lockMass<-matrix(lockMass, ncol=2, byrow=TRUE)
    if((lockMass[nrow(lockMass),2]+2) > length(ob@scanindex)){
        lockMass<-lockMass[1:(nrow(lockMass)-1),]
    }

    for(i in 1:(length(ob@scanindex)-1)){
        if(any(i == lockMass[,1])){
            arr[1,,i] <-c(object@env$mz[(object@scanindex[(i-1)]+1):object@scanindex[i]],
                          rep(NA, (max(diff(object@scanindex))-
                                   length((object@scanindex[(i-1)]+1):object@scanindex[i])) ))

            arr[2,,i] <-c(object@env$intensity[(object@scanindex[(i-1)]+1):object@scanindex[i]],
                          rep(NA, (max(diff(object@scanindex)) -
                                   length((object@scanindex[(i-1)]+1):object@scanindex[i])) ))

        } else if(any(i == lockMass[,2])){
            arr[1,,i] <-c(object@env$mz[(object@scanindex[i+1]+1):object@scanindex[(i+2)]],
                          rep(NA, (max(diff(object@scanindex)) -
                                   length((object@scanindex[i+1]+1):object@scanindex[(i+2)])) ))

            arr[2,,i] <-c(object@env$intensity[(object@scanindex[i+1]+1):object@scanindex[(i+2)]],
                          rep(NA, (max(diff(object@scanindex)) -
                                   length((object@scanindex[i+1]+1):object@scanindex[(i+2)])) ))

        } else{
            arr[1,,i] <-c(object@env$mz[(object@scanindex[i]+1):object@scanindex[i+1]],
                          rep(NA, (max(diff(object@scanindex))-
                                   length((object@scanindex[i]+1):object@scanindex[i+1])) ))

            arr[2,,i] <-c(object@env$intensity[(object@scanindex[i]+1):object@scanindex[i+1]],
                          rep(NA, (max(diff(object@scanindex)) -
                                   length((object@scanindex[i]+1):object@scanindex[i+1])) ))
        }
        ## mz is in 1; Intensity is in 2
        ##remake scanindex
        if(i == 1){
            ob@scanindex[i]<-as.integer(0)

        }else if(i == length(ob@scanindex)-1){
            ob@scanindex[i]<-as.integer(length(na.omit(arr[1,,(i-1)]))+ob@scanindex[(i-1)])
            ob@scanindex[i+1]<-as.integer(length(na.omit(arr[1,,i]))+ob@scanindex[i])
                                        #			ob@scanindex[i+1]<-as.integer(length(ob@env$mz))
        }else{
            ob@scanindex[i]<-as.integer(length(na.omit(arr[1,,(i-1)]))+ob@scanindex[(i-1)])
        }
    }

    NAidx<-is.na(arr[1,,])
    ob@env$mz<-as.numeric(arr[1,,][!NAidx])
    ob@env$intensity<-as.numeric(arr[2,,][!NAidx])
    ob<-remakeTIC(ob) ## remake TIC

    return(ob)
})

############################################################
## stitch.netCDF
setMethod("stitch.netCDF", "xcmsRaw", function(object, lockMass) {
    if(length(lockMass) == 0 | all(lockMass == 0)){
        return(object)
    }

    ob<-new("xcmsRaw")

    ob@filepath<-object@filepath
    ob@mzrange<-range(ob@env$mz)
    ob@profmethod<-object@profmethod
    ob@profparam<-list()

    arr<-array(dim=c(2,max(diff(object@scanindex)), (length(object@scanindex)+length(lockMass)) ))
    ob@scanindex <- integer(length=length(arr[1,1,]))
    ob@acquisitionNum<-1:length(ob@scanindex)

    if(lockMass[1] == 1){
        lockMass<-lockMass[3:length(lockMass)]
    }
    lockMass<-matrix(lockMass, ncol=2, byrow=TRUE)
    if((lockMass[nrow(lockMass),2]+2) > length(ob@scanindex)){
        lockMass<-lockMass[1:(nrow(lockMass)-1),]
    } ## remove the last lock mass scan if it's at the end of the run

    add<-0
    arrMax<-length(arr[1,,1])
    scanIx<-integer(length(arr[1,1,]))
    for(i in 1:length(object@scanindex)){
                                        #		if((i+add) > length(object@scanindex)){
                                        #			break
                                        #		}
        scan<-getScan(object, i)
        arr[1,,i+add]<- c(scan[,"mz"], rep(NA, (arrMax-length(scan[,"mz"])) ))
        arr[2,,i+add]<- c(scan[,"intensity"], rep(NA, (arrMax-length(scan[,"intensity"])) ))
        scanIx[(i+add)+1]<- (scanIx[(i+add)])+nrow(scan)
        ##proably going to need a cut at the end of scanIx +1 problem

        if(any(i == lockMass[,1])){
            arr[1,,i+1+add]<- c(scan[,"mz"], rep(NA, (arrMax-length(scan[,"mz"])) ))
            arr[2,,i+1+add]<- c(scan[,"intensity"], rep(NA, (arrMax-length(scan[,"intensity"])) ))
            scanIx[(i+add)+2]<- (scanIx[1+i+add])+nrow(scan)

            scan<-getScan(object, i+1)
            arr[1,,i+2+add]<- c(scan[,"mz"], rep(NA, (arrMax-length(scan[,"mz"])) ))
            arr[2,,i+2+add]<- c(scan[,"intensity"], rep(NA, (arrMax-length(scan[,"intensity"])) ))
            scanIx[(i+add)+3]<- (scanIx[(i+2)+add])+nrow(scan)

            add<-add+2
        }
    }

    NAidx<-is.na(arr[1,,])
    ob@env$mz<-as.numeric(arr[1,,][!NAidx])
    ob@env$intensity<-as.numeric(arr[2,,][!NAidx])
    ##remake scanindex
                                        #	scanInx<- as.integer(apply(arr[1,,], 2, function(x){
                                        #		inx<-is.na(x)
                                        #		length(x[!inx]) ## need to add these length together
                                        #	}))
    ob@scanindex<-as.integer(scanIx)
    ob@scantime <- sapply(1:length(ob@scanindex), function(x, time){
        time*x
    }, mean(diff(object@scantime)))
    ob<-remakeTIC(ob) ## remake TIC

    return(ob)
})

############################################################
## stitch.netCDF.new
setMethod("stitch.netCDF.new", "xcmsRaw", function(object, lockMass) {
    if(length(lockMass) == 0 | all(lockMass == 0)){
        return(object)
    }

    ob<-new("xcmsRaw")

    ob@filepath<-object@filepath
    ob@mzrange<-range(ob@env$mz)
    ob@profmethod<-object@profmethod
    ob@profparam<-list()

    arr<-array(dim=c(2,max(diff(object@scanindex)), (length(object@scanindex)+length(lockMass)) ))
    ob@scanindex <- integer(length=length(arr[1,1,]))
    ob@acquisitionNum<-1:length(ob@scanindex)

    if(lockMass[1] == 1){
        lockMass<-lockMass[2:length(lockMass)]
    }
                                        # lockMass<-matrix(lockMass, ncol=2, byrow=TRUE)

    add<-0
    arrMax<-length(arr[1,,1])
    scanIx<-integer(length(arr[1,1,]))
    for(i in 1:length(object@scanindex)){
        scan<-getScan(object, i)
        arr[1,,i+add]<- c(scan[,"mz"], rep(NA, (arrMax-length(scan[,"mz"])) ))
        arr[2,,i+add]<- c(scan[,"intensity"], rep(NA, (arrMax-length(scan[,"intensity"])) ))
        scanIx[(i+add)+1]<- (scanIx[(i+add)])+nrow(scan)

        if(any(i == lockMass)){
            arr[1,,i+1+add]  <- c(scan[,"mz"], rep(NA, (arrMax-length(scan[,"mz"])) ))
            arr[2,,i+1+add]  <- c(scan[,"intensity"], rep(NA, (arrMax-length(scan[,"intensity"])) ))
            scanIx[(i+add)+2]<- (scanIx[1+i+add])+nrow(scan)

            add<-add+1
            ## for the moment lets be dirty and add the scan before
            ## upgrade later to 1/2 and 1/2 from each scan
        }
    }

    NAidx<-is.na(arr[1,,])
    ob@env$mz<-as.numeric(arr[1,,][!NAidx])
    ob@env$intensity<-as.numeric(arr[2,,][!NAidx])
    ## above is to remove any NA buffers from the array

    ob@scanindex<-as.integer(scanIx)
    ob@scantime <- sapply(1:length(ob@scanindex), function(x, time){
        time*x
    }, mean(diff(object@scantime))) ## remake the scantime vector
    ob<-remakeTIC(ob) ## remake TIC

    return(ob)
})

