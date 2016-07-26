## All low level (API) analysis functions for feature detection should go in here.
#' @include c.R cwTools.R

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
do_detectFeatures_matchedFilter <- function() {
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


