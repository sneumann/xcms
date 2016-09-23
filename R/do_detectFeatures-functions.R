## All low level (API) analysis functions for feature detection should go in here.
#' @include c.R functions-binning.R cwTools.R

############################################################
## centWave
##
## TODO LLLL @jo update this method with the new one in methods-xcmsRaw!
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
##' @title Feature detection using the centWave method
##'
##' @description This function performs peak density and wavelet based feature
##' detection for high resulution LC/MS data in centroid mode [Tautenhahn 2008].
##'
##' @details This algorithm is most suitable for high resolution
##' LC/\{TOF,OrbiTrap,FTICR\}-MS data in centroid mode. In the first phase the
##' method identifies \emph{regions of interest} (ROIs) representing mass traces
##' that are characterized as regions with less than \code{ppm} m/z deviation in
##' consecutive scans in the LC/MS map. These ROIs are then subsequently
##' analyzed using continuous wavelet transform (CWT) to locate chromatographic
##' peaks on different scales. The first analysis step is skipped, if regions
##' of interest are passed with the \code{ROIs} parameter.
##'
##' @note The \emph{centWave} was designed to work on centroided mode, thus it
##' is expected that such data is presented to the function.
##'
##' This function exposes core feature detection functionality of
##' the \emph{centWave} method. While this function can be called directly,
##' users will generally call the corresponding method for the data object
##' instead.
##'
##' @param mz Numeric vector with the individual m/z values from all scans/
##' spectra of one file/sample.
##' @param int Numeric vector with the individual intensity values from all
##' scans/spectra of one file/sample.
##' @param scantime Numeric vector of length equal to the number of
##' spectra/scans of the data representing the retention time of each scan.
##' @param valsPerSpect Numeric vector with the number of values for each
##' spectrum.
##' @param ppm Maximal tolerated m/z deviation in consecutive scans in parts
##' per million (ppm).
##' @param peakwidth Numeric of length 2 with the expected approximate
##' feature/peak width in chromatographic space. Given as a range (min, max)
##' in seconds.
##' @param snthresh Signal to noise ratio cutoff.
##' @param prefilter Numeric of length 2: \code{c(k, I)} specifying the prefilter
##' step for the first analysis step (ROI detection). Mass traces are only
##' retained if they contain at least \code{k} peaks with intensity >= \code{I}.
##' @param mzCenterFun Name of the function to calculate the m/z center of the
##' feature. Allowed are: \code{"wMean"}: intensity weighted mean of the feature's
##' m/z values, \code{"mean"}: mean of the feature's m/z values, \code{"apex"}:
##' use the m/z value at the peak apex, \code{"wMeanApex3"}: intensity weighted
##' mean of the m/z value at the peak apex and the m/z values left and right of
##' it and \code{"meanApex3"}: mean of the m/z value of the peak apex and the
##' m/z values left and right of it.
##' @param integrate Integration method. For \code{integrate = 1} peak limits
##' are found through descent on the mexican hat filtered data, for
##' \code{integrate = 2} the descent is done on the real data. The latter method
##' is more accurate but prone to noise, while the former is more robust, but
##' less exact.
##' @param mzdiff Numeric representing the minimum difference in m/z dimension
##' for peaks with overlapping retention times; can be negatove to allow overlap.
##' @param fitgauss Logical whether or not a Gaussian should be fitted to each
##' peak.
##' @param noise Numeric of length 1 allowing to set a minimum intensity required
##' for centroids to be considered in the first analysis step (centroids with
##' intensity \code{< noise} are omitted from ROI detection).
##' @param verboseColumns Logical whether additional feature meta data columns
##' should be returned.
##' @param ROIs An optional list of regions-of-interest (ROI) representing
##' detected mass traces. If ROIs are submitted the first analysis step is
##' omitted and feature detection is performed on the submitted ROIs. Each
##' ROI object in the list is expected to have the following slots specified:
##' \code{scmin} (start scan index), \code{scmax} (end scan index),
##' \code{mzmin} (minimum m/z), \code{mzmax} (maximum m/z), \code{length}
##' (number of scans), \code{intensity} (summed intensity).
##'
##' @family core feature detection functions
##' @references
##' Ralf Tautenhahn, Christoph B\"{o}ttcher, and Steffen Neumann "Highly
##' sensitive feature detection for high resolution LC/MS" \emph{BMC Bioinformatics}
##' 2008, 9:504
##' @return
##' A matrix, each row representing an intentified feature, with columns:
##' \describe{
##' \item{mz}{Intensity weighted mean of m/z values of the feature across scans.}
##' \item{mzmin}{Minimum m/z of the feature.}
##' \item{mzmax}{Maximum m/z of the feature.}
##' \item{rt}{Retention time of the feature's midpoint.}
##' \item{rtmin}{Minimum retention time of the feature.}
##' \item{rtmax}{Maximum retention time of the feature.}
##' \item{into}{Integrated (original) intensity of the feature.}
##' \item{intb}{Per-feature baseline corrected integrated feature intensity.}
##' \item{maxo}{Maximum intensity of the feature.}
##' \item{sn}{Signal to noise ratio, defined as \code{(maxo - baseline)/sd},
##' \code{sd} being the standard deviation of local chromatographic noise.}
##' \item{egauss}{RMSE of Gaussian fit.}
##' }
##' Additional columns for \code{verboseColumns = TRUE}:
##' \describe{
##' \item{mu}{Gaussian parameter mu.}
##' \item{sigma}{Gaussian parameter sigma.}
##' \item{h}{Gaussian parameter h.}
##' \item{f}{Region number of the m/z ROI where the peak was localized.}
##' \item{dppm}{m/z deviation of mass trace across scanns in ppk.}
##' \item{scale}{Scale on which the feature was localized.}
##' \item{scpos}{Peak position found by wavelet analysis (scan number).}
##' \item{scmin}{Left peak limit found by wavelet analysis (scan number).}
##' \item{scmax}{Right peak limit found by wavelet analysis (scan numer).}
##' }
##' @author Ralf Tautenhahn, Johannes Rainer
##' @examples
##' ## Load the test file
##' library(faahKO)
##' fs <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
##' xr <- xcmsRaw(fs)
##'
##' ## Extracting the data from the xcmsRaw for do_detectFeatures_centWave
##' mzVals <- xr@env$mz
##' intVals <- xr@env$intensity
##' ## Define the values per spectrum:
##' valsPerSpect <- diff(c(xr@scanindex, length(mzVals)))
##'
##' res <- do_detectFeatures_centWave(mz = mzVals, int = intVals,
##' scantime = xr@scantime, valsPerSpect = valsPerSpect)
##' head(res)
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
    .centWave_old(mz = mz, int = int, scantime = scantime,
                  valsPerSpect = valsPerSpect, ppm = ppm, peakwidth = peakwidth,
                  snthresh = snthresh, prefilter = prefilter,
                  mzCenterFun = mzCenterFun, integrate = integrate,
                  mzdiff = mzdiff, fitgauss = fitgauss, noise = noise,
                  verboseColumns = verboseColumns, ROIs = ROIs)
}
.centWave_old <- function(mz, int, scantime, valsPerSpect,
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
    ## TODO @jo Ensure in upstream method that data is in centroided mode!
    ## TODO @jo Ensure the upstream method did eventual sub-setting on scanrange
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
    if (!is.double(mz))
        mz <- as.double(mz)
    if (!is.double(int))
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

    ## Search for potential peaks in m/z direction.
    ## If no ROIs are supplied then search for them.
    if (length(ROIs) == 0) {
        message("Detecting mass traces at ", ppm, "ppm")
        ## We're including the findmzROI code in this function to reduce
        ## the need to copy objects etc.
        ## We could also sort the data by m/z anyway; wouldn't need that
        ## much time. Once we're using classes from MSnbase we can be
        ## sure that values are correctly sorted.
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
                ##  o split values into a list -> mz per spectrum, intensity per
                ##    spectrum.
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
        eic <- .Call("getEIC", mz, int, scanindex, as.double(mzrange),
                     as.integer(sr), as.integer(length(scanindex)),
                     PACKAGE = "xcms")

        d <- eic$intensity
        td <- sr[1]:sr[2]
        scan.range <- c(sr[1],sr[2])
        ## original mzROI range; can't we extract that directly from the eic???
        ## mzROI.EIC <- .Call("getEIC", mz, int, scanindex, as.double(mzrange), as.integer(scrange),
        ##                    as.integer(length(scanindex)), PACKAGE="xcms")
        idxs <- which(eic$scan %in% seq(scrange[1], scrange[2]))
        mzROI.EIC <- list(scan=eic$scan[idxs], intensity=eic$intensity[idxs])
        ## Get the actual m/z matching these values.
        omz <- .Call("getMZ",mz, int, scanindex, as.double(mzrange),
                     as.integer(scrange),
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
                                    dppm <- round(min(running(abs(diff(mz.value)) /
                                                              (mzrange[2] * 1e-6),
                                                              fun=max,width=minCentroids)))
                                } else {
                                    dppm <- round((mzrange[2]-mzrange[1]) /
                                                  (mzrange[2] * 1e-6))
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

                pwid <- (scantime[peakrange[2]] - scantime[peakrange[1]]) /
                    (peakrange[2] - peakrange[1])
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
do_detectFeatures_massifquant <- function(mz,
                                          int,
                                          scantime,
                                          valsPerSpect,
                                          ppm = 10,
                                          peakwidth = c(20, 50),
                                          snthresh = 10,
                                          prefilter = c(3, 100),
                                          mzCenterFun = "wMean",
                                          integrate = 1,
                                          mzdiff = -0.001,
                                          fitgauss = FALSE,
                                          noise = 0,
                                          verboseColumns = FALSE,
                                          criticalValue = 1.125,
                                          consecMissedLimit = 2,
                                          unions = 1,
                                          checkBack = 0,
                                          withWave = 0) {
}
## The original code.
## Not much to speed up here; it's more code tidying.
## LLLL Test this function
.massifquant <- function(mz,
                         int,
                         scantime,
                         valsPerSpect,
                         ppm = 10,
                         peakwidth = c(20,50),
                         snthresh = 10,
                         prefilter = c(3,100),
                         mzCenterFun = "wMean",
                         integrate = 1,
                         mzdiff = -0.001,
                         fitgauss = FALSE,
                         noise = 0, ## noise.local=TRUE,
                         verboseColumns = FALSE,
                         criticalValue = 1.125,
                         consecMissedLimit = 2,
                         unions = 1,
                         checkBack = 0,
                         withWave = 0) {
    cat("\n Massifquant, Copyright (C) 2013 Brigham Young University.")
    cat("\n Massifquant comes with ABSOLUTELY NO WARRANTY.",
        " See LICENSE for details.\n", sep ="")
    flush.console()

    ## TODO @jo Ensure in upstream method that data is in centroided mode!
    ## TODO @jo Ensure the upstream method did eventual sub-setting on scanrange
    ## Input argument checking.
    if (missing(mz) | missing(int) | missing(scantime) | missing(valsPerSpect))
        stop("Arguments 'mz', 'int', 'scantime' and 'valsPerSpect'",
             " are required!")
    if (length(mz) != length(int) | length(valsPerSpect) != length(scantime)
        | length(mz) != sum(valsPerSpect))
        stop("Lengths of 'mz', 'int' and of 'scantime','valsPerSpect'",
             " have to much. Also, 'length(mz)' should be equal to",
             " 'sum(valsPerSpect)'.")
    if (!is.double(mz))
        mz <- as.double(mz)
    if (!is.double(int))
        int <- as.double(int)
    ## Fix the mzCenterFun
    mzCenterFun <- paste("mzCenter",
                         gsub(mzCenterFun, pattern = "mzCenter.",
                              replacement = "", fixed = TRUE), sep=".")
    if (!exists(mzCenterFun, mode="function"))
        stop("Error: >", mzCenterFun, "< not defined !")

    cat("\n Detecting  mass traces at ",ppm,"ppm ... \n")
    flush.console()
    massifquantROIs <- do_findKalmanROI(mz = mz, int = int, scantime = scantime,
                                        valsPerSpect = valsPerSpect,
                                        minIntensity = prefilter[2],
                                        minCentroids = peakwidth[1],
                                        criticalVal = criticalValue,
                                        consecMissedLim = consecMissedLimit,
                                        segs = unions, scanBack = checkBack,
                                        ppm = ppm)
    if (withWave == 1) {
        featlist <- do_detectFeatures_centWave(mz = mz, int = int,
                                               scantime = scantime,
                                               valsPerSpect = valsPerSpect,
                                               ppm = ppm, peakwidth = peakwidth,
                                               snthresh = snthresh,
                                               prefilter = prefilter,
                                               mzCenterFun = mzCenterFun,
                                               integrate = integrate,
                                               mzdiff = mzdiff,
                                               fitgauss = fitgauss,
                                               noise = noise,
                                               verboseColumns = verboseColumns,
                                               ROIs = massifquantROIs)
    }
    else {
        ## Get index vector for C calls
        scanindex <- valueCount2ScanIndex(valsPerSpect)
        basenames <- c("mz","mzmin","mzmax","rtmin","rtmax","rt", "into")
        if (length(massifquantROIs) == 0) {
            cat("\nNo peaks found !\n")
            nopeaks <- matrix(nrow=0, ncol=length(basenames))
            colnames(nopeaks) <- basenames
            return(nopeaks)
        }

        ## Get the max intensity for each feature.
        maxo <- lapply(massifquantROIs, function(z) {
            raw <- xcms:::.rawMat(mz = mz, int = int, scantime = scantime,
                           valsPerSpect = valsPerSpect,
                           mzrange = c(z$mzmin, z$mzmax),
                           scanrange = c(z$scmin, z$scmax))
            max(raw[, 3])
        })

        ## p <- t(sapply(massifquantROIs, unlist))
        p <- do.call(rbind, lapply(massifquantROIs, unlist, use.names = FALSE))
        colnames(p) <- basenames
        p <- cbind(p, maxo = unlist(maxo))

        #calculate median index
        p[, "rt"] <- as.integer(p[, "rtmin"] + ( (p[, "rt"] + 1) / 2 ) - 1)
        #convert from index into actual time
        p[, "rtmin"] <- scantime[p[, "rtmin"]]
        p[, "rtmax"] <- scantime[p[, "rtmax"]]
        p[, "rt"] <- scantime[p[, "rt"]]

        uorder <- order(p[, "into"], decreasing = TRUE)
        pm <- as.matrix(p[, c("mzmin", "mzmax", "rtmin", "rtmax"),
                          drop = FALSE])
        uindex <- rectUnique(pm, uorder, mzdiff, ydiff = -0.00001) ## allow adjacent peaks;
        featlist <- p[uindex, , drop = FALSE]
        cat("\n", dim(featlist)[1]," Peaks.\n");
        return(featlist)
    }
    return(featlist)
}

## The version of matchedFilter:
## .matchedFilter_orig: original code, iterative buffer creation.
## .matchedFilter_binYonX_iter: iterative buffer creation but using our binning function.
## .matchedFilter_no_iter: original binning, but a single binning call.
## .matchedFilter_binYonX_no_iter: single binning call using our binning function.

############################################################
## matchedFilter
##
##  That's the function that matches the code from the
##  findPeaks.matchedFilter method from the xcms package.
##  The peak detection is performed on the binned data. Depending on the variable
##  `bufsize` the function iteratively bins the intensity values on m/z dimension into
##  bins of size `step` always binning into `bufsize` bins. While ensuring low memory
##  usage, this iterative buffering is actually quite time consuming.
##  The loop runs over the variable `mass` which corresponds to the midpoints of the
##  bins. Peak detection is performed for bin `i` considering also `steps` neighboring
##  bins. As detailed above, if the index i is outside of the buffer size, the binnin
##  is performed for the next chunk of `bufsize` bins.
##
##  This function takes basic R-objects and might thus be used as the base analysis
##  method for a future xcms API.
##  mz is a numeric vector with all m/z values.
##  int is a numeric vector with the intensities.
##  valsPerSpect: is an integer vector with the number of values per spectrum.
##  This will be converted to what xcms calls the scanindex.
##  TODO: in the long run it would be better to avoid buffer creation, extending,
##  filling and all this stuff being done in a for loop.
## impute: none (=bin), binlin, binlinbase, intlin
## baseValue default: min(int)/2 (smallest value in the whole data set).
##
##' @title Feature detection in the chromatographic time domain
##'
##' @description This function identifies features in the chromatographic
##' time domain as described in [Smith 2006]. The intensity values are
##' binned by cutting The LC/MS data into slices (bins) of a mass unit
##' (\code{binSize} m/z) wide. Within each bin the maximal intensity is
##' selected. The feature detection is then performed in each bin by
##' extending it based on the \code{steps} parameter to generate slices
##' comprising bins \code{current_bin - steps +1} to \code{current_bin + steps - 1}.
##' Each of these slices is then filtered with matched filtration using
##' a second-derative Gaussian as the model feature/peak shape. After filtration
##' features are detected using a signal-to-ration cut-off. For more details
##' and illustrations see [Smith 2006].
##'
##' @details The intensities are binned by the provided m/z values within each
##' spectrum (scan). Binning is performed such that the bins are centered around
##' the m/z values (i.e. the first bin includes all m/z values between
##' \code{min(mz) - bin_size/2} and \code{min(mz) + bin_size/2}).
##'
##' For more details on binning and missing value imputation see
##' \code{\link{binYonX}} and \code{\link{imputeLinInterpol}} methods.
##'
##' @note
##' This function exposes core feature detection functionality of
##' the \emph{matchedFilter} method. While this function can be called directly,
##' users will generally call the corresponding method for the data object
##' instead (e.g. the \code{link{findPeaks.matchedFilter}} method).
##'
##' @inheritParams do_detectFeatures_centWave
##' @inheritParams imputeLinInterpol
##' @param binSize Numeric of length one specifying the width of the
##' bins/slices in m/z dimension.
##' @param impute Character string specifying the method to be used for missing
##' value imputation. Allowed values are \code{"none"} (no linear interpolation),
##' \code{"lin"} (linear interpolation), \code{"linbase"} (linear interpolation
##' within a certain bin-neighborhood) and \code{"intlin"}. See
##' \code{\link{imputeLinInterpol}} for more details.
##' @param fwhm Numeric of length one specifying the full width at half maximum
##' of matched filtration gaussian model peak. Only used to calculate the actual
##' sigma, see below.
##' @param sigma Numeric of length one specifying the standard deviation (width)
##' of the matched filtration model peak.
##' @param max Numeric of length one representing the maximum number of peaks
##' that are expected/will be identified per slice.
##' @param snthresh Numeric of length one defining the signal to noise cutoff
##' to be used in the feature detection step.
##' @param steps Numeric of length one defining the number of bins to be
##' merged before filtration (i.e. the number of neighboring bins that will be
##' joined to the slice in which filtration and peak detection will be
##' performed).
##' @param mzdiff Numeric of length one defining the minimum difference
##' in m/z for peaks with overlapping retention times
##' @param index Logical specifying whether indicies should be returned instead
##' of values for m/z and retention times.
##' @return A matrix, each row representing an intentified feature, with columns:
##' \describe{
##' \item{mz}{Intensity weighted mean of m/z values of the feature across scans.}
##' \item{mzmin}{Minimum m/z of the feature.}
##' \item{mzmax}{Maximum m/z of the feature.}
##' \item{rt}{Retention time of the feature's midpoint.}
##' \item{rtmin}{Minimum retention time of the feature.}
##' \item{rtmax}{Maximum retention time of the feature.}
##' \item{into}{Integrated (original) intensity of the feature.}
##' \item{intf}{Integrated intensity of the filtered peak.}
##' \item{maxo}{Maximum intensity of the feature.}
##' \item{maxf}{Maximum intensity of the filtered peak.}
##' \item{i}{Rank of feature in merged EIC (\code{<= max}).}
##' \item{sn}{Signal to noise ratio of the feature}
##' }
##' @references
##' Colin A. Smith, Elizabeth J. Want, Grace O'Maille, Ruben Abagyan and
##' Gary Siuzdak. "XCMS: Processing Mass Spectrometry Data for Metabolite
##' Profiling Using Nonlinear Peak Alignment, Matching, and Identification"
##' \emph{Anal. Chem.} 2006, 78:779-787.
##' @author Colin A Smith, Johannes Rainer
##' @family core feature detection functions
##' @seealso \code{\link{binYonX}} for a binning function,
##' \code{\link{imputeLinInterpol}} for the interpolation of missing values.
##' @examples
##' ## Load the test file
##' library(faahKO)
##' fs <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
##' xr <- xcmsRaw(fs)
##'
##' ## Extracting the data from the xcmsRaw for do_detectFeatures_centWave
##' mzVals <- xr@env$mz
##' intVals <- xr@env$intensity
##' ## Define the values per spectrum:
##' valsPerSpect <- diff(c(xr@scanindex, length(mzVals)))
##'
##' res <- do_detectFeatures_matchedFilter(mz = mzVals, int = intVals,
##' scantime = xr@scantime, valsPerSpect = valsPerSpect)
##' head(res)
do_detectFeatures_matchedFilter <- function(mz,
                                            int,
                                            scantime,
                                            valsPerSpect,
                                            binSize = 0.1,
                                            impute = "none",
                                            baseValue,
                                            distance,
                                            fwhm = 30,
                                            sigma = fwhm/2.3548,
                                            max = 5,
                                            snthresh = 10,
                                            steps = 2,
                                            mzdiff = 0.8 - binSize * steps,
                                            index = FALSE
                                            ){
    ## Use original code
    if (useOriginalCode()) {
        ## warning("Old xcms code was used; be aware that this code",
        ##         " may contain bugs.")
        return(.matchedFilter_orig(mz, int, scantime, valsPerSpect,
                                   binSize, impute, baseValue, distance,
                                   fwhm, sigma, max, snthresh,
                                   steps, mzdiff, index))
    } else {
        return(.matchedFilter_binYonX_no_iter(mz, int, scantime, valsPerSpect,
                                              binSize, impute, baseValue,
                                              distance, fwhm, sigma, max,
                                              snthresh, steps, mzdiff, index
                                              ))
    }
}
.matchedFilter_orig <- function(mz,
                                int,
                                scantime,
                                valsPerSpect,
                                binSize = 0.1,
                                impute = "none",
                                baseValue,
                                distance,
                                fwhm = 30,
                                sigma = fwhm/2.3548,
                                max = 5,
                                snthresh = 10,
                                steps = 2,
                                mzdiff = 0.8 - binSize * steps,
                                index = FALSE
                                ){
    ## Map arguments to findPeaks.matchedFilter arguments.
    step <- binSize
    profMeths <- c("profBinM", "profBinLinM", "profBinLinBaseM", "profIntLinM")
    names(profMeths) <- c("none", "lin", "linbase", "intlin")
    impute <- match.arg(impute, names(profMeths))
    profFun <- profMeths[impute]

    ## Input argument checking.
    if (missing(mz) | missing(int) | missing(scantime) | missing(valsPerSpect))
        stop("Arguments 'mz', 'int', 'scantime' and 'valsPerSpect'",
             " are required!")
    if (length(mz) != length(int) | length(valsPerSpect) != length(scantime)
        | length(mz) != sum(valsPerSpect))
        stop("Lengths of 'mz', 'int' and of 'scantime','valsPerSpect'",
             " have to much. Also, 'length(mz)' should be equal to",
             " 'sum(valsPerSpect)'.")
    ## Calculate a the "scanindex" from the number of values per spectrum:
    scanindex <- valueCount2ScanIndex(valsPerSpect)

    ## Create EIC buffer
    mrange <- range(mz)
    ## Create a numeric vector of masses; these will be the mid-points of the bins.
    mass <- seq(floor(mrange[1]/step)*step, ceiling(mrange[2]/step)*step, by = step)
    ## Calculate the /real/ bin size (as in xcms.c code).
    bin_size <- (max(mass) - min(mass)) / (length(mass) - 1)
    bufsize <- min(100, length(mass))
    ## Define profparam:
    profp <- list()
    if (missing(baseValue))
        baseValue <- min(int, na.rm = TRUE) / 2
    profp$baselevel <- baseValue
    if (!missing(distance)) {
        profp$basespace <- distance * bin_size
    } else {
        profp$basespace <- 0.075
        distance <- floor(0.075 / bin_size)
    }
    ## This returns a matrix, ncol equals the number of spectra, nrow the bufsize.
    buf <- do.call(profFun, args = list(mz, int, scanindex, bufsize, mass[1],
                                        mass[bufsize], TRUE, profp))
    bufMax <- profMaxIdxM(mz, int, scanindex, bufsize, mass[1], mass[bufsize],
                          TRUE, profp)
    bufidx <- integer(length(mass))
    idxrange <- c(1, bufsize)
    bufidx[idxrange[1]:idxrange[2]] <- 1:bufsize
    lookahead <- steps-1
    lookbehind <- 1  ## always 1?

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

    for (i in seq(length = (length(mass)-steps+1))) {
        ## Update EIC buffer if necessary
        if (bufidx[i+lookahead] == 0) {
            bufidx[idxrange[1]:idxrange[2]] <- 0
            idxrange <- c(max(1, i - lookbehind), min(bufsize+i-1-lookbehind,
                                                      length(mass)))
            bufidx[idxrange[1]:idxrange[2]] <- 1:(diff(idxrange)+1)
            buf <- do.call(profFun, args = list(mz, int, scanindex, diff(idxrange)+1,
                                                mass[idxrange[1]], mass[idxrange[2]],
                                                TRUE, profp))
            bufMax <- profMaxIdxM(mz, int, scanindex, diff(idxrange)+1,
                                  mass[idxrange[1]],
                                  mass[idxrange[2]], TRUE, profp)
        }
        ymat <- buf[bufidx[i:(i+steps-1)],,drop=FALSE]
        ysums <- colMax(ymat)
        yfilt <- filtfft(ysums, filt)
        gmax <- max(yfilt)
        ## Just look for 'max' number of peaks within the bin/slice.
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
    ## Select for each unique mzmin, mzmax, rtmin, rtmax the largest peak
    ## and report that.
    uorder <- order(rmat[,"into"], decreasing=TRUE)
    uindex <- rectUnique(rmat[,c("mzmin","mzmax","rtmin","rtmax"),drop=FALSE],
                                uorder, mzdiff)
    rmat <- rmat[uindex,,drop=FALSE]
    return(rmat)
}

############################################################
## Same as do_detectFeatures_matchedFilter except:
##
## o Using the binYtoX and imputeLinInterpol instead of the
##   profBin* methods.
## THIS IS MATTER OF REMOVAL
.matchedFilter_binYonX_iter <- function(mz,
                                        int,
                                        scantime,
                                        valsPerSpect,
                                        binSize = 0.1,
                                        impute = "none",
                                        baseValue,
                                        distance,
                                        fwhm = 30,
                                        sigma = fwhm/2.3548,
                                        max = 5,
                                        snthresh = 10,
                                        steps = 2,
                                        mzdiff = 0.8 - binSize * steps,
                                        index = FALSE
                                        ){

    ## Input argument checking.
    if (missing(mz) | missing(int) | missing(scantime) | missing(valsPerSpect))
        stop("Arguments 'mz', 'int', 'scantime' and 'valsPerSpect'",
             " are required!")
    if (length(mz) != length(int) | length(valsPerSpect) != length(scantime)
        | length(mz) != sum(valsPerSpect))
        stop("Lengths of 'mz', 'int' and of 'scantime','valsPerSpect'",
             " have to much. Also, 'length(mz)' should be equal to",
             " 'sum(valsPerSpect)'.")
    ## Get the profile/binning function: allowed: bin, lin, linbase and intlin
    impute <- match.arg(impute, c("none", "lin", "linbase", "intlin"))
    if (impute == "intlin")
        stop("Not yet implemented!")
    toIdx <- cumsum(valsPerSpect)
    fromIdx <- c(1L, toIdx[-length(toIdx)] + 1L)

    ## Create EIC buffer
    mrange <- range(mz)
    mass <- seq(floor(mrange[1]/binSize)*binSize,
                ceiling(mrange[2]/binSize)*binSize, by = binSize)
    bufsize <- min(100, length(mass))

    ## Calculate the breaks; we will re-use these in all calls.
    ## Calculate breaks and "correct" binSize; using seq ensures we're closer
    ## to the xcms profBin* results.
    binFromX <- min(mass)
    binToX <- max(mass)
    bin_size <- (binToX - binFromX) / (length(mass) - 1)
    brks <- seq(binFromX - bin_size/2, binToX + bin_size/2, by = bin_size)

    ## Problem with sequential binning is that we don't want to have the condition
    ## <= last_break in each iteration since this would cause some values being
    ## considered for multiple bins. Thus we add an additional last bin, for which
    ## we however want to get rid of later.
    binRes <- binYonX(mz, int,
                      breaks = brks[1:(bufsize+2)],
                      fromIdx = fromIdx,
                      toIdx = toIdx,
                      baseValue = ifelse(impute == "none", yes = 0, no = NA),
                      sortedX = TRUE,
                      returnIndex = TRUE
                      )
    if (length(toIdx) == 1)
        binRes <- list(binRes)
    ## Remove the last bin value unless bufsize + 2 is equal to the length of brks
    if (length(brks) > (bufsize + 2)) {
        binRes <- lapply(binRes, function(z) {
            len <- length(z$x)
            return(list(x = z$x[-len], y = z$y[-len], index = z$index[-len]))
        })
    }
    bufMax <- do.call(cbind, lapply(binRes, function(z) return(z$index)))
    bin_size <- binRes[[1]]$x[2] - binRes[[1]]$x[1]
    if (missing(baseValue))
        baseValue <- min(int, na.rm = TRUE) / 2
    if (missing(distance))
        distance <- floor(0.075 / bin_size)
    binVals <- lapply(binRes, function(z) {
        return(imputeLinInterpol(z$y, method = impute,
                                 noInterpolAtEnds = TRUE,
                                 distance = distance,
                                 baseValue = baseValue))
    })
    buf <- do.call(cbind, binVals)

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
            idxrange <- c(max(1, i - lookbehind), min(bufsize+i-1-lookbehind,
                                                      length(mass)))
            bufidx[idxrange[1]:idxrange[2]] <- 1:(diff(idxrange)+1)
            ## Avoid the problem reported above for the sequential buffering:
            ## add an additional bin for which we remove the value afterwards.
            additionalBin <- 0
            if ((idxrange[2] + 1) < length(brks)) {
                additionalBin <- 1
            }
            ## Re-fill buffer.
            binRes <- binYonX(mz, int,
                              breaks = brks[idxrange[1]:(idxrange[2] + 1 +
                                                         additionalBin)],
                              fromIdx = fromIdx,
                              toIdx = toIdx,
                              baseValue = ifelse(impute == "none", yes = 0,
                                                 no = NA),
                              sortedX = TRUE,
                              returnIndex = TRUE
                              )
            if (length(toIdx) == 1)
                binRes <- list(binRes)
            if (additionalBin == 1) {
                binRes <- lapply(binRes, function(z) {
                    len <- length(z$x)
                    return(list(x = z$x[-len], y = z$y[-len],
                                index = z$index[-len]))
                })
            }
            bufMax <- do.call(cbind, lapply(binRes, function(z) return(z$index)))
            binVals <- lapply(binRes, function(z) {
                return(imputeLinInterpol(z$y, method = impute,
                                         noInterpolAtEnds = TRUE,
                                         distance = distance,
                                         baseValue = baseValue))
            })
            buf <- do.call(cbind, binVals)
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
    colnames(rmat) <- cnames
    rmat <- rmat[seq(length = num),]
    max <- max-1 + max*(steps-1) + max*ceiling(mzdiff/binSize)
    if (index)
        mzdiff <- mzdiff/binSize
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
##  o Create the full 'profile matrix' (i.e. the m/z binned matrix) once
##    instead of repeatedly creating a "buffer" of 100 m/z values.
##  o Append the identified peaks to a list instead of generating a matrix
##    with a fixed set of rows which is doubled in its size each time more
##    peaks are identified than there are rows in the matrix.
.matchedFilter_no_iter <- function(mz,
                                   int,
                                   scantime,
                                   valsPerSpect,
                                   binSize = 0.1,
                                   impute = "none",
                                   baseValue,
                                   distance,
                                   fwhm = 30,
                                   sigma = fwhm/2.3548,
                                   max = 5,
                                   snthresh = 10,
                                   steps = 2,
                                   mzdiff = 0.8 - binSize * steps,
                                   index = FALSE
                                   ){
    ## Map arguments to findPeaks.matchedFilter arguments.
    step <- binSize
    profMeths <- c("profBinM", "profBinLinM", "profBinLinBaseM", "profIntLinM")
    names(profMeths) <- c("none", "lin", "linbase", "intlin")
    impute <- match.arg(impute, names(profMeths))
    profFun <- profMeths[impute]

    ## Input argument checking.
    if (missing(mz) | missing(int) | missing(scantime) | missing(valsPerSpect))
        stop("Arguments 'mz', 'int', 'scantime' and 'valsPerSpect'",
             " are required!")
    if (length(mz) != length(int) | length(valsPerSpect) != length(scantime)
        | length(mz) != sum(valsPerSpect))
        stop("Lengths of 'mz', 'int' and of 'scantime','valsPerSpect'",
             " have to much. Also, 'length(mz)' should be equal to",
             " 'sum(valsPerSpect)'.")
    ## Calculate a the "scanindex" from the number of values per spectrum:
    scanindex <- valueCount2ScanIndex(valsPerSpect)

    ## Create the full profile matrix.
    mrange <- range(mz)
    mass <- seq(floor(mrange[1]/step)*step, ceiling(mrange[2]/step)*step, by = step)
    ## Calculate the /real/ bin size (as in xcms.c code).
    bin_size <- (max(mass) - min(mass)) / (length(mass) - 1)
    ## bufsize <- min(100, length(mass))
    bufsize <- length(mass)
    ## Define profparam:
    profp <- list()
    if (!missing(baseValue))
        profp$baselevel <- baseValue
    if (!missing(distance))
        profp$basespace <- distance * bin_size
    ## This returns a matrix, ncol equals the number of spectra, nrow the bufsize.
    buf <- do.call(profFun, args = list(mz, int, scanindex, bufsize, mass[1],
                                        mass[bufsize], TRUE, profp))

    ## The full matrix, nrow is the total number of (binned) m/z values.
    bufMax <- profMaxIdxM(mz, int, scanindex, bufsize, mass[1], mass[bufsize],
                          TRUE, profp)
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
    for (i in seq(length = (length(mass)-steps+1))) {

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
    if (length(ResList) == 0) {
        rmat <- matrix(nrow = 0, ncol = length(cnames))
        colnames(rmat) <- cnames
        return(rmat)
    }
    rmat <- do.call(rbind, ResList)
    if (is.null(dim(rmat))) {
        rmat <- matrix(rmat, nrow = 1)
    }
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
##  o Create the full 'profile matrix' (i.e. the m/z binned matrix) once
##    instead of repeatedly creating a "buffer" of 100 m/z values.
##  o Append the identified peaks to a list instead of generating a matrix
##    with a fixed set of rows which is doubled in its size each time more
##    peaks are identified than there are rows in the matrix.
##  o Use binYonX and imputeLinInterpol instead of the profBin... methods.
.matchedFilter_binYonX_no_iter <- function(mz,
                                           int,
                                           scantime,
                                           valsPerSpect,
                                           binSize = 0.1,
                                           impute = "none",
                                           baseValue,
                                           distance,
                                           fwhm = 30,
                                           sigma = fwhm/2.3548,
                                           max = 5,
                                           snthresh = 10,
                                           steps = 2,
                                           mzdiff = 0.8 - binSize * steps,
                                           index = FALSE
                                           ){
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
    mass <- seq(floor(mrange[1] / binSize) * binSize,
                ceiling(mrange[2] / binSize) * binSize,
                by = binSize)
    ## Get the imputation method: allowed: none, lin, linbase and intlin
    impute <- match.arg(impute, c("none", "lin", "linbase", "intlin"))
    ## Select the profFun and the settings for it...
    if (impute == "intlin") {
        ## intlin not yet implemented...
        profFun = "profIntLinM"
        profp <- list()
        ## Calculate a the "scanindex" from the number of values per spectrum:
        scanindex <- valueCount2ScanIndex(valsPerSpect)
        bufsize <- length(mass)
        buf <- do.call(profFun, args = list(mz, int, scanindex, bufsize,
                                            mass[1], mass[bufsize],
                                            TRUE))
        ## The full matrix, nrow is the total number of (binned) m/z values.
        bufMax <- profMaxIdxM(mz, int, scanindex, bufsize, mass[1],
                              mass[bufsize], TRUE, profp)
    } else {
        ## Binning the data.
        ## Create and translate settings for binYonX
        toIdx <- cumsum(valsPerSpect)
        fromIdx <- c(1L, toIdx[-length(toIdx)] + 1L)
        shiftBy = TRUE
        binFromX <- min(mass)
        binToX <- max(mass)
        ## binSize <- (binToX - binFromX) / (length(mass) - 1)
        ## brks <- seq(binFromX - binSize/2, binToX + binSize/2, by = binSize)
        brks <- breaks_on_nBins(fromX = binFromX, toX = binToX,
                                nBins = length(mass), shiftByHalfBinSize = TRUE)
        binRes <- binYonX(mz, int,
                          breaks = brks,
                          fromIdx = fromIdx,
                          toIdx = toIdx,
                          baseValue = ifelse(impute == "none", yes = 0, no = NA),
                          sortedX = TRUE,
                          returnIndex = TRUE
                          )
        if (length(toIdx) == 1)
            binRes <- list(binRes)
        bufMax <- do.call(cbind, lapply(binRes, function(z) return(z$index)))
        bin_size <- binRes[[1]]$x[2] - binRes[[1]]$x[1]
        ## Missing value imputation
        if (missing(baseValue))
            baseValue <- min(int, na.rm = TRUE) / 2
        if (missing(distance))
            distance <- floor(0.075 / bin_size)
        binVals <- lapply(binRes, function(z) {
            return(imputeLinInterpol(z$y, method = impute, distance = distance,
                                     noInterpolAtEnds = TRUE,
                                     baseValue = baseValue))
        })
        buf <- do.call(cbind, binVals)
    }

    bufidx <- 1L:length(mass)
    lookahead <- steps-1
    lookbehind <- 1

    N <- nextn(length(scantime))
    xrange <- range(scantime)
    x <- c(0:(N/2), -(ceiling(N/2-1)):-1) *
        (xrange[2]-xrange[1])/(length(scantime)-1)

    filt <- -attr(eval(deriv3(~ 1/(sigma*sqrt(2*pi)) *
                                  exp(-x^2/(2*sigma^2)), "x")), "hessian")
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
                massmean <- weighted.mean(mzmat, intmat[which.intMax],
                                          na.rm = TRUE)
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
                ResList[[num]] <- c(massmean, mzrange[1], mzrange[2], maxy,
                                    peakrange, into, intf, maxo, maxf, j, sn)
            } else
                break
        }
    }
    if (length(ResList) == 0) {
        rmat <- matrix(nrow = 0, ncol = length(cnames))
        colnames(rmat) <- cnames
        return(rmat)
    }
    rmat <- do.call(rbind, ResList)
    if (is.null(dim(rmat))) {
        rmat <- matrix(rmat, nrow = 1)
    }
    colnames(rmat) <- cnames
    max <- max-1 + max*(steps-1) + max*ceiling(mzdiff/binSize)
    if (index)
        mzdiff <- mzdiff/binSize
    else {
        rmat[, "rt"] <- scantime[rmat[, "rt"]]
        rmat[, "rtmin"] <- scantime[rmat[, "rtmin"]]
        rmat[, "rtmax"] <- scantime[rmat[, "rtmax"]]
    }
    ## Select for each unique mzmin, mzmax, rtmin, rtmax the largest peak
    ## and report that.
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
##' @title Feature detection for single-spectrum non-chromatography MS data
##'
##' @description This function performs feature detection in mass spectrometry
##' direct injection spectrum using a wavelet based algorithm.
##'
##' @details This is a wrapper around the peak picker in Bioconductor's
##' \code{MassSpecWavelet} package calling
##' \code{\link[MassSpecWavelet]{peakDetectionCWT}} and
##' \code{\link[MassSpecWavelet]{tuneInPeakInfo}} functions.
##'
##' @inheritParams do_detectFeatures_centWave
##' @param ... Additional parameters to be passed to the
##' \code{\link[MassSpecWavelet]{peakDetectionCWT}} function.
##'
##' @return
##' A matrix, each row representing an intentified feature, with columns:
##' \describe{
##' \item{mz}{m/z value of the feature at the centroid position.}
##' \item{mzmin}{Minimum m/z of the feature.}
##' \item{mzmax}{Maximum m/z of the feature.}
##' \item{rt}{Always \code{-1}.}
##' \item{rtmin}{Always \code{-1}.}
##' \item{rtmax}{Always \code{-1}.}
##' \item{into}{Integrated (original) intensity of the feature.}
##' \item{maxo}{Maximum intensity of the feature.}
##' \item{intf}{Always \code{NA}.}
##' \item{maxf}{Maximum MSW-filter response of the feature.}
##' \item{sn}{Signal to noise ratio.}
##' }
##'
##' @family core feature detection functions
##' @seealso \code{\link[MassSpecWavelet]{peakDetectionCWT}} from the \code{MassSpecWavelet}.
##' @author Joachim Kutzera, Steffen Neumann, Johannes Rainer
do_detectFeatures_MSW <- function(mz, int, snthresh = 3,
                                  verboseColumns = FALSE, ...) {
    ## Input argument checking.
    if (missing(int))
        stop("Argument 'int' is missing!")
    if (missing(mz))
        stop("Argument 'mz' is missing!")
    if (length(int) != length(mz))
        stop("Length of 'int' and 'mz' do not match!")
    if (!is.numeric(int) | !is.numeric(mz))
        stop("'int' and 'mz' are supposed to be numeric vectors!")

    .MSW(int = int, mz = mz, snthresh = snthresh,
         verboseColumns = verboseColumns, ...)
}
############################################################
## The original code
## This should be removed at some point.
.MSW_orig <- function(mz, int, snthresh = 3, verboseColumns = FALSE, ...) {

    ## MassSpecWavelet Calls
    peakInfo <- peakDetectionCWT(int, SNR.Th=snthresh, ...)
    majorPeakInfo <- peakInfo$majorPeakInfo

    sumIntos <- function(into, inpos, scale){
        scale=floor(scale)
        sum(into[(inpos-scale):(inpos+scale)])
    }

    maxIntos <- function(into, inpos, scale){
        scale=floor(scale)
              max(into[(inpos-scale):(inpos+scale)])
    }

    betterPeakInfo <- tuneInPeakInfo(int,
                                     majorPeakInfo)

    peakIndex <- betterPeakInfo$peakIndex

    ## sum and max of raw values, sum and max of filter-response
    rints<-NA;fints<-NA
    maxRints<-NA;maxFints<-NA

    for (a in 1:length(peakIndex)) {
        rints[a] <- sumIntos(int,peakIndex[a],
                             betterPeakInfo$peakScale[a])
        maxRints[a] <- maxIntos(int,peakIndex[a],
                                betterPeakInfo$peakScale[a])
    }
    ## filter-response is not summed here, the maxF-value is the one
    ## which was "xcmsRaw$into" earlier

    ## Assemble result

    basenames <- c("mz","mzmin","mzmax","rt","rtmin","rtmax",
                   "into","maxo","sn","intf","maxf")

    peaklist <- matrix(-1, nrow = length(peakIndex), ncol = length(basenames))

    colnames(peaklist) <- c(basenames)

    peaklist[,"mz"] <- mz[peakIndex]
    peaklist[,"mzmin"] <- mz[(peakIndex-betterPeakInfo$peakScale)]
    peaklist[,"mzmax"] <- mz[(peakIndex+betterPeakInfo$peakScale)]

    peaklist[,"rt"]    <- rep(-1, length(peakIndex))
    peaklist[,"rtmin"] <- rep(-1, length(peakIndex))
    peaklist[,"rtmax"] <- rep(-1, length(peakIndex))

    peaklist[,"into"] <- rints ## sum of raw-intensities
    peaklist[,"maxo"] <- maxRints
    peaklist[,"intf"] <- rep(NA, length(peakIndex))
    peaklist[,"maxf"] <- betterPeakInfo$peakValue

    peaklist[,"sn"]   <- betterPeakInfo$peakSNR

    cat('\n')

    ## Filter additional (verbose) columns
    if (!verboseColumns)
        peaklist <- peaklist[,basenames,drop=FALSE]

    peaklist
}
############################################################
## Slightly modified and tuned original code
.MSW <- function(mz, int, snthresh = 3, verboseColumns = FALSE, ...) {

    ## MassSpecWavelet Calls
    peakInfo <- peakDetectionCWT(int, SNR.Th = snthresh, ...)
    majorPeakInfo <- peakInfo$majorPeakInfo

    sumIntos <- function(into, inpos, scale){
        scale = floor(scale)
        sum(into[(inpos-scale):(inpos+scale)])
    }
    maxIntos <- function(into, inpos, scale){
        scale = floor(scale)
        max(into[(inpos-scale):(inpos+scale)])
    }
    betterPeakInfo <- tuneInPeakInfo(int,
                                     majorPeakInfo)
    peakIndex <- betterPeakInfo$peakIndex
    nPeaks <- length(peakIndex)

    ## sum and max of raw values, sum and max of filter-response
    rints <- numeric(nPeaks)
    fints <- NA
    maxRints <- numeric(nPeaks)
    maxFints <- NA

    for (a in 1:nPeaks) {
        rints[a] <- sumIntos(int, peakIndex[a],
                             betterPeakInfo$peakScale[a])
        maxRints[a] <- maxIntos(int, peakIndex[a],
                                betterPeakInfo$peakScale[a])
    }
    ## filter-response is not summed here, the maxF-value is the one
    ## which was "xcmsRaw$into" earlier

    ## Assemble result
    basenames <- c("mz","mzmin","mzmax","rt","rtmin","rtmax",
                   "into","maxo","sn","intf","maxf")

    peaklist <- matrix(-1, nrow = nPeaks, ncol = length(basenames))
    colnames(peaklist) <- c(basenames)

    peaklist[,"mz"] <- mz[peakIndex]
    peaklist[,"mzmin"] <- mz[(peakIndex - betterPeakInfo$peakScale)]
    peaklist[,"mzmax"] <- mz[(peakIndex + betterPeakInfo$peakScale)]

    ## peaklist[,"rt"]    <- rep(-1, length(peakIndex))
    ## peaklist[,"rtmin"] <- rep(-1, length(peakIndex))
    ## peaklist[,"rtmax"] <- rep(-1, length(peakIndex))

    peaklist[,"into"] <- rints ## sum of raw-intensities
    peaklist[,"maxo"] <- maxRints
    peaklist[,"intf"] <- rep(NA, nPeaks)
    peaklist[,"maxf"] <- betterPeakInfo$peakValue

    peaklist[,"sn"]   <- betterPeakInfo$peakSNR

    cat('\n')

    ## Filter additional (verbose) columns
    if (!verboseColumns)
        peaklist <- peaklist[,basenames,drop=FALSE]

    peaklist
}



############################################################
## MS1
##
do_detectFeatures_MS1 <- function() {
}


############################################################
## do_findKalmanROI
do_findKalmanROI <- function(mz, int, scantime, valsPerSpect,
                             mzrange = c(0.0, 0.0),
                             scanrange = c(1, length(scantime)),
                             minIntensity, minCentroids, consecMissedLim,
                             criticalVal, ppm, segs, scanBack) {
    if (missing(mz) | missing(int) | missing(scantime) | missing(valsPerSpect))
        stop("Arguments 'mz', 'int', 'scantime' and 'valsPerSpect'",
             " are required!")
    if (length(mz) != length(int) | length(valsPerSpect) != length(scantime)
        | length(mz) != sum(valsPerSpect))
        stop("Lengths of 'mz', 'int' and of 'scantime','valsPerSpect'",
             " have to much. Also, 'length(mz)' should be equal to",
             " 'sum(valsPerSpect)'.")
    scanindex <- valueCount2ScanIndex(valsPerSpect) ## Get index vector for C calls
    ## Call the C function.
    if (!is.double(mz))
        mz <- as.double(mz)
    if (!is.double(int))
        int <- as.double(int)
    if (!is.integer(scanindex))
        scanindex <- as.integer(scanindex)
    if (!is.double(scantime))
        scantime <- as.double(scantime)
    .Call("massifquant", mz, int, scanindex, scantime, as.double(mzrange),
          as.integer(scanrange), as.integer(length(scantime)),
          as.double(minIntensity),as.integer(minCentroids),as.double(consecMissedLim),
          as.double(ppm), as.double(criticalVal), as.integer(segs), as.integer(scanBack), PACKAGE ='xcms' )
}
