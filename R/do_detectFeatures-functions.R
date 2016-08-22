## All low level (API) analysis functions for feature detection should go in here.
#' @include c.R functions-binning.R cwTools.R

############################################################
## centWave
##
## Some notes on a potential speed up:
## Tried:
## o initialize peaks matrix in the inner loop instead of rbind: slower.
## o pre-initialize the peaks list; slower.
## o keep the peaks matrix small, add additional columns only if fitgauss or
##   verboseColumns is TRUE.
## Conclusion:
## o speed improvement can only come from internal methods called withihn.
##
##
do_detectFeatures_centWave <- function(mz, int, scantime, valsPerSpect,
                                       ppm = 25,
                                       peakwidth = c(20, 50),
                                       snthresh = 10,
                                       prefilter = c(3, 100),
                                       mzCenterFun = "wMean",
                                       integrate = 1,
                                       mzdiff = -0.001,
                                       fitgauss = FALSE,
                                       noise = 0,
                                       verboseColumns = FALSE,
                                       ROIs = list()) {
    ## TODO @jo Ensure in upstream method that data is in centroided mode!
    ## Input argument checking.
    if (missing(mz) | missing(int) | missing(scantime) | missing(valsPerSpect))
        stop("Arguments 'mz', 'int', 'scantime' and 'valsPerSpect'",
             " are required!")
    if (length(mz) != length(int) | length(valsPerSpect) != length(scantime)
        | length(mz) != sum(valsPerSpect))
        stop("Lengths of 'mz', 'int' and of 'scantime','valsPerSpect'",
             " have to much. Also, 'length(mz)' should be equal to",
             " 'sum(valsPerSpect)'.")
    scanindex <- valueCount2ScanIndex(valsPerSpect) ## Get index vector for C calls
    mz <- as.double(mz)
    int <- as.double(int)
    ## Fix the mzCenterFun
    mzCenterFun <- paste("mzCenter",
                         gsub(mzCenterFun, pattern = "mzCenter.",
                              replacement = "", fixed = TRUE), sep=".")
    if (!exists(mzCenterFun, mode="function"))
        stop("Error: >", mzCenterFun, "< not defined !")

    ## Define the result column names.
    basenames <- c("mz","mzmin","mzmax","rt","rtmin","rtmax",
                   "into","intb","maxo","sn")
    verbosenames <- c("egauss","mu","sigma","h","f", "dppm", "scale","scpos",
                      "scmin","scmax","lmin","lmax")

    ## Peak width: seconds to scales
    scalerange <- round((peakwidth / mean(diff(scantime))) / 2)
    if(length(z <- which(scalerange==0)))
        scalerange <- scalerange[-z]
    if(length(scalerange) < 1)
        stop("No scales ? Please check peak width!\n")
    if(length(scalerange) > 1){
        scales <- seq(from=scalerange[1], to=scalerange[2], by=2)
    }else{
        scales <- scalerange
    }

    ## Define more variables.
    minPeakWidth <-  scales[1]
    noiserange <- c(minPeakWidth*3, max(scales)*3)
    maxGaussOverlap <- 0.5
    minPtsAboveBaseLine <- max(4, minPeakWidth - 2)
    minCentroids <- minPtsAboveBaseLine
    scRangeTol <-  maxDescOutlier <- floor(minPeakWidth/2)
    scanrange <- c(1, length(scantime))

    ## Search for potential peaks in M/Z direction.
    ## If no ROIs are supplied then search for them.
    if (length(ROIs) == 0) {
        message("Detecting mass traces at ", ppm, "ppm")
        ## We're including the findmzROI code in this function to reduce the need to copy
        ## objects etc.
        ## We could also sort the data by M/Z anyway; wouldn't need that much time.
        withRestarts(
            tryCatch({
                ROIs <- .Call("findmzROI",
                              mz, int, scanindex,
                              as.double(c(0.0, 0.0)),
                              as.integer(scanrange),
                              as.integer(length(scantime)),
                              as.double(ppm * 1e-6),
                              as.integer(minCentroids),
                              as.integer(prefilter),
                              as.integer(noise),
                              PACKAGE ='xcms' )
            },
            error=function(e){if (grepl("m/z sort assumption violated !", e$message))
                              {invokeRestart("fixSort")} else {simpleError(e)}}),
            fixSort = function() {
                ## Force ordering of values within spectrum by mz:
                ##  o split values into a list -> mz per spectrum, intensity per spectrum.
                ##  o define the ordering.
                ##  o re-order the mz and intensity and unlist again.
                ## Note: the Rle split is faster than the "conventional" factor split.
                splitF <- Rle(1:length(valsPerSpect), valsPerSpect)
                mzl <- as.list(S4Vectors::split(mz, f = splitF))
                oidx <- lapply(mzl, order)
                mz <<- unlist(mapply(mzl, oidx, FUN = function(y, z) {
                    return(y[z])
                }, SIMPLIFY = FALSE, USE.NAMES = FALSE), use.names = FALSE)
                int <<- unlist(mapply(as.list(split(int, f = splitF)), oidx,
                                      FUN=function(y, z) {
                                          return(y[z])
                                      }, SIMPLIFY = FALSE, USE.NAMES = FALSE),
                               use.names = FALSE)
                rm(mzl)
                rm(splitF)
                ROIs <<- .Call("findmzROI",
                                   mz, int, scanindex,
                                   as.double(c(0.0, 0.0)),
                                   as.integer(scanrange),
                                   as.integer(length(scantime)),
                                   as.double(ppm * 1e-6),
                                   as.integer(minCentroids),
                                   as.integer(prefilter),
                                   as.integer(noise),
                                   PACKAGE ='xcms' )
            }
        )

        ## ROI.list <- findmzROI(object, scanrange=scanrange, dev=ppm * 1e-6,
        ##                       minCentroids=minCentroids, prefilter=prefilter,
        ##                       noise=noise)
        if (length(ROIs) == 0) {
            warning("No ROIs found!")
            if (verboseColumns) {
                nopeaks <- matrix(nrow=0, ncol=length(basenames)+length(verbosenames))
                colnames(nopeaks) <- c(basenames, verbosenames)
            } else {
                nopeaks <- matrix(nrow=0, ncol=length(basenames))
                colnames(nopeaks) <- c(basenames)
            }
            return(nopeaks)
        }
    }

    ## Process the ROIs
    peaklist <- list()
    Nscantime <- length(scantime)
    lf <- length(ROIs)

    ## We might want to replace that with a lapply!
    for (f in  1:lf) {
        feat <- ROIs[[f]]
        N <- feat$scmax - feat$scmin + 1

        peaks <- peakinfo <- NULL
        mzrange <- c(feat$mzmin, feat$mzmax)
        sccenter <- feat$scmin[1] + floor(N/2) - 1
        scrange <- c(feat$scmin, feat$scmax)
        ## scrange + noiserange, used for baseline detection and wavelet analysis
        sr <- c(max(scanrange[1], scrange[1] - max(noiserange)),
                min(scanrange[2], scrange[2] + max(noiserange)))

        ## Directly call the C:
        eic <- .Call("getEIC", mz, int, scanindex, as.double(mzrange), as.integer(sr),
                     as.integer(length(scanindex)), PACKAGE = "xcms")

        d <- eic$intensity
        td <- sr[1]:sr[2]
        scan.range <- c(sr[1],sr[2])
        ## original mzROI range; can't we extract that directly from the eic???
        ## mzROI.EIC <- .Call("getEIC", mz, int, scanindex, as.double(mzrange), as.integer(scrange),
        ##                    as.integer(length(scanindex)), PACKAGE="xcms")
        idxs <- which(eic$scan %in% seq(scrange[1], scrange[2]))
        mzROI.EIC <- list(scan=eic$scan[idxs], intensity=eic$intensity[idxs])
        ## Get the actual M/Z matching these values.
        omz <- .Call("getMZ",mz, int, scanindex, as.double(mzrange), as.integer(scrange),
                     as.integer(length(scantime)), PACKAGE = 'xcms')

        if (all(omz == 0))
            stop("centWave: debug me: (omz == 0)?\n")
        od  <- mzROI.EIC$intensity
        otd <- mzROI.EIC$scan
        if (all(od == 0))
            stop("centWave: debug me: (all(od == 0))?\n")

        ## scrange + scRangeTol, used for gauss fitting and continuous
        ## data above 1st baseline detection
        ftd <- max(td[1], scrange[1] - scRangeTol) : min(td[length(td)],
                                                         scrange[2] + scRangeTol)
        fd <- d[match(ftd,td)]

        ## 1st type of baseline: statistic approach
        if (N >= 10 * minPeakWidth) {
            ## in case of very long mass trace use full scan range for
            ## baseline detection
            noised <- .Call("getEIC", mz, int, scanindex, as.double(mzrange),
                            as.integer(scanrange), as.integer(length(scanindex)),
                            PACKAGE="xcms")$intensity
        }else{
            noised <- d
        }
        ## 90% trimmed mean as first baseline guess
        noise <- estimateChromNoise(noised, trim = 0.05, minPts = 3 * minPeakWidth)

        ## any continuous data above 1st baseline ?
        if(!continuousPtsAboveThreshold(fd, threshold = noise,
                                        num = minPtsAboveBaseLine))
            next

        ## 2nd baseline estimate using not-peak-range
        lnoise <- getLocalNoiseEstimate(d, td, ftd, noiserange, Nscantime,
                                        threshold = noise,
                                        num = minPtsAboveBaseLine)

        ## Final baseline & Noise estimate
        baseline <- max(1, min(lnoise[1], noise))
        sdnoise <- max(1, lnoise[2])
        sdthr <-  sdnoise * snthresh

        ## is there any data above S/N * threshold ?
        if (!(any(fd - baseline >= sdthr)))
            next

        wCoefs <- MSW.cwt(d, scales = scales, wavelet = 'mexh')
        if (!(!is.null(dim(wCoefs)) && any(wCoefs- baseline >= sdthr)))
            next

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
                            ## try to decide which scale describes the peak best
                            inti <- numeric(length(opp))
                            irange = rep(ceiling(scales[1]/2), length(opp))
                            for (k in 1:length(opp)) {
                                kpos <- opp[k]
                                r1 <- ifelse(kpos-irange[k] > 1, kpos-irange[k], 1)
                                r2 <- ifelse(kpos+irange[k] < length(d),
                                             kpos+irange[k], length(d))
                                inti[k] <- sum(d[r1:r2])
                            }
                            maxpi <- which.max(inti)
                            if (length(maxpi) > 1) {
                                m <- wCoefs[opp[maxpi], maxpi]
                                bestcol <- which(m == max(m), arr.ind = TRUE)[2]
                                best.scale.nr <- maxpi[bestcol]
                            } else best.scale.nr <- maxpi

                            best.scale <-  scales[best.scale.nr]
                            best.scale.pos <- opp[best.scale.nr]

                            pprange <- min(pp):max(pp)
                            ## maxint <- max(d[pprange])
                            lwpos <- max(1, best.scale.pos - best.scale)
                            rwpos <- min(best.scale.pos + best.scale, length(td))
                            p1 <- match(td[lwpos], otd)[1]
                            p2 <- match(td[rwpos], otd); p2 <- p2[length(p2)]
                            if (is.na(p1)) p1 <- 1
                            if (is.na(p2)) p2 <- N
                            mz.value <- omz[p1:p2]
                            mz.int <- od[p1:p2]
                            maxint <- max(mz.int)

                            ## re-calculate m/z value for peak range
                            mzrange <- range(mz.value)
                            mzmean <- do.call(mzCenterFun,
                                              list(mz = mz.value,
                                                   intensity = mz.int))

                            ## Compute dppm only if needed
                            dppm <- NA
                            if (verboseColumns) {
                                if (length(mz.value) >= (minCentroids+1)) {
                                    dppm <- round(min(running(abs(diff(mz.value)) /(mzrange[2] * 1e-6),
                                                              fun=max,width=minCentroids)))
                                } else {
                                    dppm <- round((mzrange[2]-mzrange[1]) / (mzrange[2] * 1e-6))
                                }
                            }
                            peaks <- rbind(peaks,
                                           c(mzmean, mzrange,           ## mz
                                             NA, NA, NA,                ## rt, rtmin, rtmax,
                                             NA,                        ## intensity (sum)
                                             NA,                        ## intensity (-bl)
                                             maxint,                    ## max intensity
                                             round((maxint - baseline) / sdnoise),  ##  S/N Ratio
                                             NA,                        ## Gaussian RMSE
                                             NA,NA,NA,                  ## Gaussian Parameters
                                             f,                         ## ROI Position
                                             dppm,                      ## max. difference between the [minCentroids] peaks in ppm
                                             best.scale,                ## Scale
                                             td[best.scale.pos], td[lwpos], td[rwpos],  ## Peak positions guessed from the wavelet's (scan nr)
                                             NA,NA ))                   ## Peak limits (scan nr)

                            peakinfo <- rbind(peakinfo,
                                              c(best.scale, best.scale.nr,
                                                best.scale.pos, lwpos, rwpos))
                            ## Peak positions guessed from the wavelet's
                        }
                    }
                }
            }  ##for
        } ## if

        ##  postprocessing
        if (!is.null(peaks)) {
            colnames(peaks) <- c(basenames, verbosenames)

            colnames(peakinfo) <- c("scale","scaleNr","scpos","scmin","scmax")

            ## Exchange with lapply?
            for (p in 1:nrow(peaks)) {
                ## find minima, assign rt and intensity values
                if (integrate == 1) {
                    lm <- descendMin(wCoefs[, peakinfo[p, "scaleNr"]],
                                     istart = peakinfo[p, "scpos"])
                    gap <- all(d[lm[1]:lm[2]] == 0) ## looks like we got stuck in a gap right in the middle of the peak
                    if ((lm[1] == lm[2]) || gap )     ## fall-back
                        lm <- descendMinTol(d,
                                            startpos = c(peakinfo[p, "scmin"],
                                                         peakinfo[p, "scmax"]),
                                            maxDescOutlier)
                } else
                    lm <- descendMinTol(d, startpos = c(peakinfo[p,"scmin"],
                                                        peakinfo[p,"scmax"]),
                                        maxDescOutlier)

                ## narrow down peak rt boundaries by skipping zeros
                pd <- d[lm[1]:lm[2]]
                np <- length(pd)
                lm.l <-  xcms:::findEqualGreaterUnsorted(pd, 1)
                lm.l <- max(1, lm.l - 1)
                lm.r <- xcms:::findEqualGreaterUnsorted(rev(pd), 1)
                lm.r <- max(1, lm.r - 1)
                lm <- lm + c(lm.l - 1, -(lm.r - 1))

                peakrange <- td[lm]
                peaks[p, "rtmin"] <- scantime[peakrange[1]]
                peaks[p, "rtmax"] <- scantime[peakrange[2]]

                peaks[p, "maxo"] <- max(d[lm[1]:lm[2]])

                pwid <- (scantime[peakrange[2]] - scantime[peakrange[1]])/(peakrange[2] - peakrange[1])
                if (is.na(pwid))
                    pwid <- 1

                peaks[p, "into"] <- pwid*sum(d[lm[1]:lm[2]])

                db <-  d[lm[1]:lm[2]] - baseline
                peaks[p, "intb"] <- pwid * sum(db[db>0])

                peaks[p, "lmin"] <- lm[1]
                peaks[p, "lmax"] <- lm[2]

                if (fitgauss) {
                    ## perform gaussian fits, use wavelets for inital parameters
                    md <- max(d[lm[1]:lm[2]])
                    d1 <- d[lm[1]:lm[2]]/md ## normalize data for gaussian error calc.
                    pgauss <- fitGauss(td[lm[1]:lm[2]], d[lm[1]:lm[2]],
                                       pgauss = list(mu = peaks[p,"scpos"],
                                                     sigma=peaks[p,"scmax"] -
                                                         peaks[p,"scmin"],
                                                     h = peaks[p,"maxo"]))
                    rtime <- peaks[p, "scpos"]
                    if (!any(is.na(pgauss)) && all(pgauss > 0)) {
                        gtime <- td[match(round(pgauss$mu), td)]
                        if (!is.na(gtime)) {
                            rtime <- gtime
                            peaks[p, "mu"] <- pgauss$mu
                            peaks[p, "sigma"] <- pgauss$sigma
                            peaks[p, "h"] <- pgauss$h;
                            peaks[p, "egauss"] <- sqrt((1/length(td[lm[1]:lm[2]])) *
                                                       sum(((d1 - gauss(td[lm[1]:lm[2]],
                                                                        pgauss$h/md,
                                                                        pgauss$mu,
                                                                        pgauss$sigma))^2)))
                        }
                    }
                    peaks[p, "rt"] <- scantime[rtime]
                    ## avoid fitting side effects
                    if (peaks[p, "rt"] < peaks[p, "rtmin"])
                        peaks[p, "rt"] <- scantime[peaks[p, "scpos"]]
                } else
                    peaks[p, "rt"] <- scantime[peaks[p, "scpos"]]
            }
            peaks <- joinOverlappingPeaks(td, d, otd, omz, od, scantime,
                                          scan.range, peaks, maxGaussOverlap,
                                          mzCenterFun = mzCenterFun)
        }

        if (!is.null(peaks)) {
            peaklist[[length(peaklist)+1]] <- peaks
        }

    } ## f

    if (length(peaklist) == 0) {
        warning("No peaks found!")

        if (verboseColumns) {
            nopeaks <- matrix(nrow=0, ncol=length(basenames)+length(verbosenames))
            colnames(nopeaks) <- c(basenames, verbosenames)
        } else {
            nopeaks <- matrix(nrow=0, ncol=length(basenames))
            colnames(nopeaks) <- c(basenames)
        }

        return(invisible(nopeaks))
    }

    p <- do.call(rbind, peaklist)

    if (!verboseColumns)
        p <- p[, basenames, drop=FALSE]

    uorder <- order(p[,"into"], decreasing=TRUE)
    pm <- as.matrix(p[,c("mzmin","mzmax","rtmin","rtmax"),drop=FALSE])
    uindex <- rectUnique(pm,uorder,mzdiff,ydiff = -0.00001) ## allow adjacent peaks
    pr <- p[uindex,,drop=FALSE]

    return(pr)

}




############################################################
## massifquant
##
do_detectFeatures_massifquant <- function() {
}

############################################################
## matchedFilter
##
##  That's the function that matches the code from the
##  findPeaks.matchedFilter method from the xcms package.
##  This function takes basic R-objects and might thus be used as the base analysis
##  method for a future xcms API.
##  mz is a numeric vector with all M/Z values.
##  int is a numeric vector with the intensities.
##  valsPerSpect: is an integer vector with the number of values per spectrum. This will
##     be converted to what xcms calls the scanindex.
##  TODO: in the long run it would be better to avoid buffer creation, extending, filling
##     and all this stuff being done in a for loop.
## profFun: bin, binlin, binlinbase, intlin
do_detectFeatures_matchedFilter <- function(mz,
                                            int,
                                            scantime,
                                            valsPerSpect,
                                            profFun = "bin",
                                            profparam=list(),
                                            fwhm = 30,
                                            sigma = fwhm/2.3548,
                                            max = 5,
                                            snthresh = 10,
                                            step = 0.1,
                                            steps = 2,
                                            mzdiff = 0.8 - step * steps,
                                            index = FALSE,
                                            verboseColumns = FALSE){

    ## Input argument checking.
    if (missing(mz) | missing(int) | missing(scantime) | missing(valsPerSpect))
        stop("Arguments 'mz', 'int', 'scantime' and 'valsPerSpect'",
             " are required!")
    if (length(mz) != length(int) | length(valsPerSpect) != length(scantime)
        | length(mz) != sum(valsPerSpect))
        stop("Lengths of 'mz', 'int' and of 'scantime','valsPerSpect'",
             " have to much. Also, 'length(mz)' should be equal to",
             " 'sum(valsPerSpect)'.")
    ## get the profile/binning function:
    profFun <- match.arg(profFun, names(.profFunctions))
    profFun <- .profFunctions[[profFun]]
    ## Calculate a the "scanindex" from the number of values per spectrum:
    scanindex <- valueCount2ScanIndex(valsPerSpect)

    ## Create EIC buffer
    mrange <- range(mz)
    mass <- seq(floor(mrange[1]/step)*step, ceiling(mrange[2]/step)*step, by = step)
    bufsize <- min(100, length(mass))
    ## This returns a matrix, ncol equals the number of spectra, nrow the bufsize.
    buf <- do.call(profFun, args = list(mz, int, scanindex, bufsize, mass[1],
                                        mass[bufsize], TRUE, profparam))
    ## buf <- profFun(mz, int, scanindex, bufsize, mass[1], mass[bufsize],
    ##                TRUE, profparam)
    bufMax <- profMaxIdxM(mz, int, scanindex, bufsize, mass[1], mass[bufsize],
                          TRUE, profparam)
    bufidx <- integer(length(mass))
    idxrange <- c(1, bufsize)
    bufidx[idxrange[1]:idxrange[2]] <- 1:bufsize
    lookahead <- steps-1
    lookbehind <- 1

    N <- nextn(length(scantime))
    xrange <- range(scantime)
    x <- c(0:(N/2), -(ceiling(N/2-1)):-1)*(xrange[2]-xrange[1])/(length(scantime)-1)

    filt <- -attr(eval(deriv3(~ 1/(sigma*sqrt(2*pi))*exp(-x^2/(2*sigma^2)), "x")), "hessian")
    filt <- filt/sqrt(sum(filt^2))
    filt <- fft(filt, inverse = TRUE)/length(filt)

    cnames <- c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "into", "intf",
                "maxo", "maxf", "i", "sn")
    rmat <- matrix(nrow = 2048, ncol = length(cnames))
    num <- 0

    for (i in seq(length = length(mass)-steps+1)) {
        ## Update EIC buffer if necessary
        if (bufidx[i+lookahead] == 0) {
            bufidx[idxrange[1]:idxrange[2]] <- 0
            idxrange <- c(max(1, i - lookbehind), min(bufsize+i-1-lookbehind, length(mass)))
            bufidx[idxrange[1]:idxrange[2]] <- 1:(diff(idxrange)+1)
            buf <- do.call(profFun, args = list(mz, int, scanindex, diff(idxrange)+1,
                                                mass[idxrange[1]], mass[idxrange[2]],
                                                TRUE, profparam))
            ## buf <- profFun(mz, int, scanindex, diff(idxrange)+1, mass[idxrange[1]],
            ##                mass[idxrange[2]], TRUE, profparam)
            bufMax <- profMaxIdxM(mz, int, scanindex, diff(idxrange)+1, mass[idxrange[1]],
                                  mass[idxrange[2]], TRUE, profparam)
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
                intmat <- ymat[, peakrange[1]:peakrange[2], drop = FALSE]
                mzmat <- matrix(mz[bufMax[bufidx[i:(i+steps-1)],
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

                pwid <- (scantime[peakrange[2]] - scantime[peakrange[1]]) /
                    (peakrange[2] - peakrange[1])
                into <- pwid*sum(ysums[peakrange[1]:peakrange[2]])
                intf <- pwid*sum(yfilt[peakrange[1]:peakrange[2]])
                maxo <- max(ysums[peakrange[1]:peakrange[2]])
                maxf <- yfilt[maxy]
                yfilt[peakrange[1]:peakrange[2]] <- 0
                num <- num + 1
                ## Double the size of the output matrix if it's full
                if (num > nrow(rmat)) {
                    nrmat <- matrix(nrow = 2*nrow(rmat), ncol = ncol(rmat))
                    nrmat[seq(length = nrow(rmat)),] = rmat
                    rmat <- nrmat
                }
                rmat[num,] <- c(massmean, mzrange[1], mzrange[2], maxy, peakrange,
                                into, intf, maxo, maxf, j, sn)
            } else
                break
        }
    }
    ## cat("\n")
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
    ## Select for each unique mzmin, mzmax, rtmin, rtmax the largest peak and report that.
    uorder <- order(rmat[,"into"], decreasing=TRUE)
    uindex <- rectUnique(rmat[,c("mzmin","mzmax","rtmin","rtmax"),drop=FALSE],
                                uorder, mzdiff)
    rmat <- rmat[uindex,,drop=FALSE]
    return(rmat)
}

############################################################
## The code of this function is basically the same than of the original
## findPeaks.matchedFilter method in xcms with the following differences:
##  o Create the full 'profile matrix' (i.e. the M/Z binned matrix) once instead of
##    repeatedly creating a "buffer" of 100 M/Z values.
##  o Append the identified peaks to a list instead of generating a matrix with a fixed
##    set of rows which is doubled in its size each time more peaks are identified than
##    there are rows in the matrix.
do_detectFeatures_matchedFilter_new <- function(mz,
                                                int,
                                                scantime,
                                                valsPerSpect,
                                                profFun = "bin",
                                                profparam = list(),
                                                fwhm = 30,
                                                sigma = fwhm/2.3548,
                                                max = 5,
                                                snthresh = 10,
                                                step = 0.1,
                                                steps = 2,
                                                mzdiff = 0.8 - step * steps,
                                                index = FALSE,
                                                verboseColumns = FALSE){

    ## Input argument checking.
    if (missing(mz) | missing(int) | missing(scantime) | missing(valsPerSpect))
        stop("Arguments 'mz', 'int', 'scantime' and 'valsPerSpect'",
             " are required!")
    if (length(mz) != length(int) | length(valsPerSpect) != length(scantime)
        | length(mz) != sum(valsPerSpect))
        stop("Lengths of 'mz', 'int' and of 'scantime','valsPerSpect'",
             " have to much. Also, 'length(mz)' should be equal to",
             " 'sum(valsPerSpect)'.")
    ## get the profile/binning function:
    profFun <- match.arg(profFun, names(.profFunctions))
    profFun <- .profFunctions[[profFun]]
    ## Calculate a the "scanindex" from the number of values per spectrum:
    scanindex <- valueCount2ScanIndex(valsPerSpect)

    ## Create the full profile matrix.
    mrange <- range(mz)
    mass <- seq(floor(mrange[1]/step)*step, ceiling(mrange[2]/step)*step, by = step)
    ## bufsize <- min(100, length(mass))
    bufsize <- length(mass)
    ## This returns a matrix, ncol equals the number of spectra, nrow the bufsize.
    ## buf <- profFun(mz, int, scanindex, bufsize, mass[1], mass[bufsize],
    ##                TRUE, profparam)
    buf <- do.call(profFun, args = list(mz, int, scanindex, bufsize, mass[1], mass[bufsize],
                   TRUE, profparam))

    ## The full matrix, nrow is the total number of (binned) M/Z values.
    bufMax <- profMaxIdxM(mz, int, scanindex, bufsize, mass[1], mass[bufsize],
                          TRUE, profparam)
    ## bufidx <- integer(length(mass))
    ## idxrange <- c(1, bufsize)
    ## bufidx[idxrange[1]:idxrange[2]] <- 1:bufsize
    bufidx <- 1L:length(mass)
    lookahead <- steps-1
    lookbehind <- 1

    N <- nextn(length(scantime))
    xrange <- range(scantime)
    x <- c(0:(N/2), -(ceiling(N/2-1)):-1)*(xrange[2]-xrange[1])/(length(scantime)-1)

    filt <- -attr(eval(deriv3(~ 1/(sigma*sqrt(2*pi))*exp(-x^2/(2*sigma^2)), "x")), "hessian")
    filt <- filt/sqrt(sum(filt^2))
    filt <- fft(filt, inverse = TRUE)/length(filt)

    cnames <- c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "into", "intf",
                "maxo", "maxf", "i", "sn")
    num <- 0

    ResList <- list()

    ## Can not do much here, lapply/apply won't work because of the 'steps' parameter.
    ## That's looping through the masses, i.e. rows of the profile matrix.
    for (i in seq(length = length(mass)-steps+1)) {

        ymat <- buf[bufidx[i:(i+steps-1)], , drop = FALSE]
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
                intmat <- ymat[, peakrange[1]:peakrange[2], drop = FALSE]
                mzmat <- matrix(mz[bufMax[bufidx[i:(i+steps-1)],
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

                pwid <- (scantime[peakrange[2]] - scantime[peakrange[1]]) /
                    (peakrange[2] - peakrange[1])
                into <- pwid*sum(ysums[peakrange[1]:peakrange[2]])
                intf <- pwid*sum(yfilt[peakrange[1]:peakrange[2]])
                maxo <- max(ysums[peakrange[1]:peakrange[2]])
                maxf <- yfilt[maxy]
                yfilt[peakrange[1]:peakrange[2]] <- 0
                num <- num + 1
                ResList[[num]] <- c(massmean, mzrange[1], mzrange[2], maxy, peakrange,
                                    into, intf, maxo, maxf, j, sn)
            } else
                break
        }
    }
    rmat <- do.call(rbind, ResList)
    colnames(rmat) <- cnames
    max <- max-1 + max*(steps-1) + max*ceiling(mzdiff/step)
    if (index)
        mzdiff <- mzdiff/step
    else {
        rmat[, "rt"] <- scantime[rmat[, "rt"]]
        rmat[, "rtmin"] <- scantime[rmat[, "rtmin"]]
        rmat[, "rtmax"] <- scantime[rmat[, "rtmax"]]
    }
    ## Select for each unique mzmin, mzmax, rtmin, rtmax the largest peak and report that.
    uorder <- order(rmat[, "into"], decreasing = TRUE)
    uindex <- rectUnique(rmat[, c("mzmin", "mzmax", "rtmin", "rtmax"),
                              drop = FALSE],
                         uorder, mzdiff)
    rmat <- rmat[uindex,,drop = FALSE]
    return(rmat)
}

############################################################
## The code of this function is basically the same than of the original
## findPeaks.matchedFilter method in xcms with the following differences:
##  o Create the full 'profile matrix' (i.e. the M/Z binned matrix) once instead of
##    repeatedly creating a "buffer" of 100 M/Z values.
##  o Append the identified peaks to a list instead of generating a matrix with a fixed
##    set of rows which is doubled in its size each time more peaks are identified than
##    there are rows in the matrix.
##  o Use binYonX and imputeLinInterpol instead of the profBin... methods.
do_detectFeatures_matchedFilter_newer <- function(mz,
                                                  int,
                                                  scantime,
                                                  valsPerSpect,
                                                  profFun = "bin",
                                                  fwhm = 30,
                                                  sigma = fwhm/2.3548,
                                                  max = 5,
                                                  snthresh = 10,
                                                  step = 0.1,
                                                  steps = 2,
                                                  mzdiff = 0.8 - step * steps,
                                                  index = FALSE,
                                                  verboseColumns = FALSE,
                                                  ...){
    ## ... arguments are passed down to the binning function.

    ## Input argument checking.
    if (missing(mz) | missing(int) | missing(scantime) | missing(valsPerSpect))
        stop("Arguments 'mz', 'int', 'scantime' and 'valsPerSpect'",
             " are required!")
    if (length(mz) != length(int) | length(valsPerSpect) != length(scantime)
        | length(mz) != sum(valsPerSpect))
        stop("Lengths of 'mz', 'int' and of 'scantime','valsPerSpect'",
             " have to much. Also, 'length(mz)' should be equal to",
             " 'sum(valsPerSpect)'.")

    ## Generate the 'profile' matrix, i.e. perform the binning:
    mrange <- range(mz)
    mass <- seq(floor(mrange[1]/step)*step, ceiling(mrange[2]/step)*step, by = step)
    ## Get the profile/binning function: allowed: bin, binlin, binlinbase and intlin
    profFun <- match.arg(profFun, names(.profFunctions))
    ## Select the profFun and the settings for it...
    if (profFun == "intlin") {
        profFun = "profIntLinM"
        ## Calculate a the "scanindex" from the number of values per spectrum:
        scanindex <- valueCount2ScanIndex(valsPerSpect)
        bufsize <- length(mass)
        ## This returns a matrix, ncol equals the number of spectra, nrow the bufsize.
        ## buf <- profFun(mz, int, scanindex, bufsize, mass[1], mass[bufsize],
        ##                TRUE, profparam)
        buf <- do.call(profFun, args = list(mz, int, scanindex, bufsize, mass[1], mass[bufsize],
                                            TRUE, profparam))

        ## The full matrix, nrow is the total number of (binned) M/Z values.
        bufMax <- profMaxIdxM(mz, int, scanindex, bufsize, mass[1], mass[bufsize],
                              TRUE, profparam)
    } else {
        cat("Using binYonX\n")
        ## Define argument imputeMethod
        if (profFun == "bin") {
            imputeMethod = "no"
        } else if (profFun == "binlin") {
            imputeMethod = "lin"
        } else if (profFun == "binlinbase") {
            imputeMethod = "linbase"
        } else {
            stop("Don't know profile function '", profFun, "'!")
        }
        profFun = "binYonX"
        ## Create and translate settings for binYonX
        toX <- cumsum(valsPerSpect)
        fromX <- c(1L, toX[-length(toX)] + 1L)
        shiftBy = TRUE
        binRes <- binYonX(mz, int, fromIdx = fromX, toIdx = toX, binFromX = mass[1],
                          binToX = mass[length(mass)], shiftByHalfBinSize = shiftBy,
                          impute = imputeMethod, sortedY = TRUE, binSize = step,
                          returnIndex = TRUE)
        bufMax <- unlist(lapply(binRes, function(z) return(z$index)))
        buf <- do.call(cbind, lapply(binRes, function(z) return(z$y)))
    }

    bufidx <- 1L:length(mass)
    lookahead <- steps-1
    lookbehind <- 1

    N <- nextn(length(scantime))
    xrange <- range(scantime)
    x <- c(0:(N/2), -(ceiling(N/2-1)):-1)*(xrange[2]-xrange[1])/(length(scantime)-1)

    filt <- -attr(eval(deriv3(~ 1/(sigma*sqrt(2*pi))*exp(-x^2/(2*sigma^2)), "x")), "hessian")
    filt <- filt/sqrt(sum(filt^2))
    filt <- fft(filt, inverse = TRUE)/length(filt)

    cnames <- c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "into", "intf",
                "maxo", "maxf", "i", "sn")
    num <- 0

    ResList <- list()

    ## Can not do much here, lapply/apply won't work because of the 'steps' parameter.
    ## That's looping through the masses, i.e. rows of the profile matrix.
    for (i in seq(length = length(mass)-steps+1)) {

        ymat <- buf[bufidx[i:(i+steps-1)], , drop = FALSE]
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
                intmat <- ymat[, peakrange[1]:peakrange[2], drop = FALSE]
                mzmat <- matrix(mz[bufMax[bufidx[i:(i+steps-1)],
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

                pwid <- (scantime[peakrange[2]] - scantime[peakrange[1]]) /
                    (peakrange[2] - peakrange[1])
                into <- pwid*sum(ysums[peakrange[1]:peakrange[2]])
                intf <- pwid*sum(yfilt[peakrange[1]:peakrange[2]])
                maxo <- max(ysums[peakrange[1]:peakrange[2]])
                maxf <- yfilt[maxy]
                yfilt[peakrange[1]:peakrange[2]] <- 0
                num <- num + 1
                ResList[[num]] <- c(massmean, mzrange[1], mzrange[2], maxy, peakrange,
                                    into, intf, maxo, maxf, j, sn)
            } else
                break
        }
    }
    rmat <- do.call(rbind, ResList)
    colnames(rmat) <- cnames
    max <- max-1 + max*(steps-1) + max*ceiling(mzdiff/step)
    if (index)
        mzdiff <- mzdiff/step
    else {
        rmat[, "rt"] <- scantime[rmat[, "rt"]]
        rmat[, "rtmin"] <- scantime[rmat[, "rtmin"]]
        rmat[, "rtmax"] <- scantime[rmat[, "rtmax"]]
    }
    ## Select for each unique mzmin, mzmax, rtmin, rtmax the largest peak and report that.
    uorder <- order(rmat[, "into"], decreasing = TRUE)
    uindex <- rectUnique(rmat[, c("mzmin", "mzmax", "rtmin", "rtmax"),
                              drop = FALSE],
                         uorder, mzdiff)
    rmat <- rmat[uindex,,drop = FALSE]
    return(rmat)
}


############################################################
## MSW
##
do_detectFeatures_MSW <- function() {
}

############################################################
## MS1
##
do_detectFeatures_MS1 <- function() {
}


