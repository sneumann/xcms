## All low level (API) analysis functions for chromatographic peak detection
## should go in here.
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
#' @title Core API function for centWave peak detection
#'
#' @description This function performs peak density and wavelet based
#'     chromatographic peak detection for high resolution LC/MS data in centroid
#'     mode [Tautenhahn 2008].
#'
#' @details
#'
#' This algorithm is most suitable for high resolution
#' LC/\{TOF,OrbiTrap,FTICR\}-MS data in centroid mode. In the first phase
#' the method identifies \emph{regions of interest} (ROIs) representing
#' mass traces that are characterized as regions with less than \code{ppm}
#' m/z deviation in consecutive scans in the LC/MS map. In detail, starting
#' with a single m/z, a ROI is extended if a m/z can be found in the next scan
#' (spectrum) for which the difference to the mean m/z of the ROI is smaller
#' than the user defined \code{ppm} of the m/z. The mean m/z of the ROI is then
#' updated considering also the newly included m/z value.
#'
#' These ROIs are then, after some cleanup, analyzed using continuous wavelet
#' transform (CWT) to locate chromatographic peaks on different scales. The
#' first analysis step is skipped, if regions of interest are passed with
#' the \code{roiList} parameter.
#'
#' @note The \emph{centWave} was designed to work on centroided mode, thus it
#'     is expected that such data is presented to the function.
#'
#'     This function exposes core chromatographic peak detection functionality
#'     of the \emph{centWave} method. While this function can be called
#'     directly, users will generally call the corresponding method for the
#'     data object instead.
#'
#' @param mz Numeric vector with the individual m/z values from all scans/
#'     spectra of one file/sample.
#' 
#' @param int Numeric vector with the individual intensity values from all
#'     scans/spectra of one file/sample.
#' 
#' @param scantime Numeric vector of length equal to the number of
#'     spectra/scans of the data representing the retention time of each scan.
#' 
#' @param valsPerSpect Numeric vector with the number of values for each
#'     spectrum.
#'
#' @param sleep \code{numeric(1)} defining the number of seconds to wait between
#'     iterations. Defaults to \code{sleep = 0}. If \code{> 0} a plot is
#'     generated visualizing the identified chromatographic peak. Note: this
#'     argument is for backward compatibility only and will be removed in
#'     future.
#' 
#' @inheritParams findChromPeaks-centWave
#'
#' @family core peak detection functions
#'
#' @references
#' Ralf Tautenhahn, Christoph B\"{o}ttcher, and Steffen Neumann "Highly
#'     sensitive feature detection for high resolution LC/MS"
#'     \emph{BMC Bioinformatics} 2008, 9:504
#'
#' @return
#'     A matrix, each row representing an identified chromatographic peak,
#'     with columns:
#'     \describe{
#' 
#'     \item{mz}{Intensity weighted mean of m/z values of the peak across
#'     scans.}
#'     \item{mzmin}{Minimum m/z of the peak.}
#'     \item{mzmax}{Maximum m/z of the peak.}
#'     \item{rt}{Retention time of the peak's midpoint.}
#'     \item{rtmin}{Minimum retention time of the peak.}
#'     \item{rtmax}{Maximum retention time of the peak.}
#'     \item{into}{Integrated (original) intensity of the peak.}
#'     \item{intb}{Per-peak baseline corrected integrated peak intensity.}
#'     \item{maxo}{Maximum intensity of the peak.}
#'     \item{sn}{Signal to noise ratio, defined as \code{(maxo - baseline)/sd},
#'     \code{sd} being the standard deviation of local chromatographic noise.}
#'     \item{egauss}{RMSE of Gaussian fit.}
#'     }
#'     Additional columns for \code{verboseColumns = TRUE}:
#'     \describe{
#' 
#'     \item{mu}{Gaussian parameter mu.}
#'     \item{sigma}{Gaussian parameter sigma.}
#'     \item{h}{Gaussian parameter h.}
#'     \item{f}{Region number of the m/z ROI where the peak was localized.}
#'     \item{dppm}{m/z deviation of mass trace across scans in ppm.}
#'     \item{scale}{Scale on which the peak was localized.}
#'     \item{scpos}{Peak position found by wavelet analysis (scan number).}
#'     \item{scmin}{Left peak limit found by wavelet analysis (scan number).}
#'     \item{scmax}{Right peak limit found by wavelet analysis (scan numer).}
#'     }
#' 
#' @author Ralf Tautenhahn, Johannes Rainer
#'
#' @seealso \code{\link{centWave}} for the standard user interface method.
#'
#' @examples
#' ## Load the test file
#' library(faahKO)
#' fs <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
#' xr <- xcmsRaw(fs, profstep = 0)
#'
#' ## Extracting the data from the xcmsRaw for do_findChromPeaks_centWave
#' mzVals <- xr@env$mz
#' intVals <- xr@env$intensity
#' ## Define the values per spectrum:
#' valsPerSpect <- diff(c(xr@scanindex, length(mzVals)))
#'
#' ## Calling the function. We're using a large value for noise to speed up
#' ## the call in the example performance - in a real use case we would either
#' ## set the value to a reasonable value or use the default value.
#' res <- do_findChromPeaks_centWave(mz = mzVals, int = intVals,
#' scantime = xr@scantime, valsPerSpect = valsPerSpect, noise = 10000)
#' head(res)
do_findChromPeaks_centWave <- function(mz, int, scantime, valsPerSpect,
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
                                       roiList = list(),
                                       firstBaselineCheck = TRUE,
                                       roiScales = NULL,
                                       sleep = 0) {
    if (getOption("originalCentWave", default = TRUE)) {
        ## message("DEBUG: using original centWave.")
        .centWave_orig(mz = mz, int = int, scantime = scantime,
                       valsPerSpect = valsPerSpect, ppm = ppm, peakwidth = peakwidth,
                       snthresh = snthresh, prefilter = prefilter,
                       mzCenterFun = mzCenterFun, integrate = integrate,
                       mzdiff = mzdiff, fitgauss = fitgauss, noise = noise,
                       verboseColumns = verboseColumns, roiList = roiList,
                       firstBaselineCheck = firstBaselineCheck,
                       roiScales = roiScales, sleep = sleep)
    } else {
        ## message("DEBUG: using modified centWave.")
        .centWave_new(mz = mz, int = int, scantime = scantime,
                      valsPerSpect = valsPerSpect, ppm = ppm, peakwidth = peakwidth,
                      snthresh = snthresh, prefilter = prefilter,
                      mzCenterFun = mzCenterFun, integrate = integrate,
                      mzdiff = mzdiff, fitgauss = fitgauss, noise = noise,
                      verboseColumns = verboseColumns, roiList = roiList,
                      firstBaselineCheck = firstBaselineCheck,
                      roiScales = roiScales, sleep = sleep)
    }
}
############################################################
## ORIGINAL code from xcms_1.49.7
.centWave_orig <- function(mz, int, scantime, valsPerSpect,
                           ppm = 25, peakwidth = c(20,50), snthresh = 10,
                           prefilter = c(3,100), mzCenterFun = "wMean",
                           integrate = 1, mzdiff = -0.001, fitgauss = FALSE,
                           noise = 0, ## noise.local=TRUE,
                           sleep = 0, verboseColumns = FALSE, roiList = list(),
                           firstBaselineCheck = TRUE, roiScales = NULL) {
    ## TODO @jo Ensure in upstream method that data is in centroided mode!
    ## TODO @jo Ensure the upstream method did eventual sub-setting on scanrange
    ## Input argument checking.
    if (missing(mz) | missing(int) | missing(scantime) | missing(valsPerSpect))
        stop("Arguments 'mz', 'int', 'scantime' and 'valsPerSpect'",
             " are required!")
    if (length(mz) != length(int) | length(valsPerSpect) != length(scantime)
        | length(mz) != sum(valsPerSpect))
        stop("Lengths of 'mz', 'int' and of 'scantime','valsPerSpect'",
             " have to match. Also, 'length(mz)' should be equal to",
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
        stop("Function '", mzCenterFun, "' not defined !")

    if (!is.logical(firstBaselineCheck))
        stop("Parameter 'firstBaselineCheck' should be logical!")
    if (length(firstBaselineCheck) != 1)
        stop("Parameter 'firstBaselineCheck' should be a single logical !")
    if (length(roiScales) > 0)
        if (length(roiScales) != length(roiList) | !is.numeric(roiScales))
            stop("If provided, parameter 'roiScales' has to be a numeric with",
                 " length equal to the length of 'roiList'!")
    ## if (!is.null(roiScales)) {
    ##     if (!is.numeric(roiScales) | length(roiScales) != length(roiList))
    ##         stop("Parameter 'roiScales' has to be a numeric of length equal to",
    ##              " parameter 'roiList'!")
    ##}

    basenames <- c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax",
                   "into", "intb", "maxo", "sn")
    verbosenames <- c("egauss", "mu", "sigma", "h", "f", "dppm", "scale",
                      "scpos", "scmin", "scmax", "lmin", "lmax")

    ## Peak width: seconds to scales
    scalerange <- round((peakwidth / mean(diff(scantime))) / 2)

    if (length(z <- which(scalerange == 0)))
        scalerange <- scalerange[-z]
    if (length(scalerange) < 1) {
        warning("No scales? Please check peak width!")
        if (verboseColumns) {
            nopeaks <- matrix(nrow = 0, ncol = length(basenames) +
                                            length(verbosenames))
            colnames(nopeaks) <- c(basenames, verbosenames)
        } else {
            nopeaks <- matrix(nrow = 0, ncol = length(basenames))
            colnames(nopeaks) <- c(basenames)
        }
        return(invisible(nopeaks))
    }

    if (length(scalerange) > 1)
        scales <- seq(from = scalerange[1], to = scalerange[2], by = 2)
    else
        scales <- scalerange

    minPeakWidth <-  scales[1]
    noiserange <- c(minPeakWidth * 3, max(scales) * 3)
    maxGaussOverlap <- 0.5
    minPtsAboveBaseLine <- max(4, minPeakWidth - 2)
    minCentroids <- minPtsAboveBaseLine
    scRangeTol <-  maxDescOutlier <- floor(minPeakWidth / 2)
    scanrange <- c(1, length(scantime))

    ## If no ROIs are supplied then search for them.
    if (length(roiList) == 0) {
        message("Detecting mass traces at ", ppm, " ppm ... ", appendLF = FALSE)
        ## flush.console();
        ## We're including the findmzROI code in this function to reduce
        ## the need to copy objects etc.
        ## We could also sort the data by m/z anyway; wouldn't need that
        ## much time. Once we're using classes from MSnbase we can be
        ## sure that values are correctly sorted.
        withRestarts(
            tryCatch({
                tmp <- capture.output(
                    roiList <- .Call("findmzROI",
                                     mz, int, scanindex,
                                     as.double(c(0.0, 0.0)),
                                     as.integer(scanrange),
                                     as.integer(length(scantime)),
                                     as.double(ppm * 1e-6),
                                     as.integer(minCentroids),
                                     as.integer(prefilter),
                                     as.integer(noise),
                                     PACKAGE ='xcms' )
                )
            },
            error = function(e){
                if (grepl("m/z sort assumption violated !", e$message)) {
                    invokeRestart("fixSort")
                } else {
                    simpleError(e)
                }
            }),
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
                tmp <- capture.output(
                    roiList <<- .Call("findmzROI",
                                      mz, int, scanindex,
                                      as.double(c(0.0, 0.0)),
                                      as.integer(scanrange),
                                      as.integer(length(scantime)),
                                      as.double(ppm * 1e-6),
                                      as.integer(minCentroids),
                                      as.integer(prefilter),
                                      as.integer(noise),
                                      PACKAGE ='xcms' )
                )
            }
        )
        message("OK")
        ## ROI.list <- findmzROI(object,scanrange=scanrange,dev=ppm * 1e-6,minCentroids=minCentroids, prefilter=prefilter, noise=noise)
        if (length(roiList) == 0) {
            warning("No ROIs found! \n")
            if (verboseColumns) {
                nopeaks <- matrix(nrow = 0, ncol = length(basenames) +
                                                length(verbosenames))
                colnames(nopeaks) <- c(basenames, verbosenames)
            } else {
                nopeaks <- matrix(nrow = 0, ncol = length(basenames))
                colnames(nopeaks) <- c(basenames)
            }
            return(invisible(nopeaks))
        }
    }

    ## Second stage: process the ROIs
    peaklist <- list()
    Nscantime <- length(scantime)
    lf <- length(roiList)

    ## cat('\n Detecting chromatographic peaks ... \n % finished: ')
    ## lp <- -1
    message("Detecting chromatographic peaks in ", length(roiList),
            " regions of interest ...", appendLF = FALSE)

    for (f in  1:lf) {

        ## ## Show progress
        ## perc <- round((f/lf) * 100)
        ## if ((perc %% 10 == 0) && (perc != lp))
        ## {
        ##     cat(perc," ",sep="");
        ##     lp <- perc;
        ## }
        ## flush.console()

        feat <- roiList[[f]]
        N <- feat$scmax - feat$scmin + 1
        peaks <- peakinfo <- NULL
        mzrange <- c(feat$mzmin, feat$mzmax)
        sccenter <- feat$scmin[1] + floor(N/2) - 1
        scrange <- c(feat$scmin, feat$scmax)
        ## scrange + noiserange, used for baseline detection and wavelet analysis
        sr <- c(max(scanrange[1], scrange[1] - max(noiserange)),
                min(scanrange[2], scrange[2] + max(noiserange)))
        eic <- .Call("getEIC", mz, int, scanindex, as.double(mzrange),
                     as.integer(sr), as.integer(length(scanindex)),
                     PACKAGE = "xcms")
        ## eic <- rawEIC(object,mzrange=mzrange,scanrange=sr)
        d <- eic$intensity
        td <- sr[1]:sr[2]
        scan.range <- c(sr[1], sr[2])
        ## original mzROI range
        idxs <- which(eic$scan %in% seq(scrange[1], scrange[2]))
        mzROI.EIC <- list(scan=eic$scan[idxs], intensity=eic$intensity[idxs])
        ## mzROI.EIC <- rawEIC(object,mzrange=mzrange,scanrange=scrange)
        omz <- .Call("getMZ", mz, int, scanindex, as.double(mzrange),
                     as.integer(scrange), as.integer(length(scantime)),
                     PACKAGE = 'xcms')
        ## omz <- rawMZ(object,mzrange=mzrange,scanrange=scrange)
        if (all(omz == 0)) {
            warning("centWave: no peaks found in ROI.")
            next
        }
        od  <- mzROI.EIC$intensity
        otd <- mzROI.EIC$scan
        if (all(od == 0)) {
            warning("centWave: no peaks found in ROI.")
            next
        }

        ## scrange + scRangeTol, used for gauss fitting and continuous
        ## data above 1st baseline detection
        ftd <- max(td[1], scrange[1] - scRangeTol) : min(td[length(td)],
                                                         scrange[2] + scRangeTol)
        fd <- d[match(ftd, td)]

        ## 1st type of baseline: statistic approach
        if (N >= 10*minPeakWidth) {
            ## in case of very long mass trace use full scan range
            ## for baseline detection
            noised <- .Call("getEIC", mz, int, scanindex, as.double(mzrange),
                            as.integer(scanrange), as.integer(length(scanindex)),
                            PACKAGE="xcms")$intensity
            ## noised <- rawEIC(object,mzrange=mzrange,scanrange=scanrange)$intensity
        } else {
            noised <- d
        }
        ## 90% trimmed mean as first baseline guess
        noise <- estimateChromNoise(noised, trim = 0.05,
                                    minPts = 3 * minPeakWidth)
        ## any continuous data above 1st baseline ?
        if (firstBaselineCheck &
            !continuousPtsAboveThreshold(fd, threshold = noise,
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
            wCoefs[nrow(wCoefs),] <- wCoefs[nrow(wCoefs) - 1, ] * 0.99
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
                            ## if(!is.null(roiScales)) {
                            ## allow roiScales to be a numeric of length 0
                            if(length(roiScales) > 0) {
                                ## use given scale
                                best.scale.nr <- which(scales == roiScales[[f]])
                                if(best.scale.nr > length(opp))
                                    best.scale.nr <- length(opp)
                            } else {
                                ## try to decide which scale describes the peak best
                                inti <- numeric(length(opp))
                                irange <- rep(ceiling(scales[1]/2), length(opp))
                                for (k in 1:length(opp)) {
                                    kpos <- opp[k]
                                    r1 <- ifelse(kpos - irange[k] > 1,
                                                 kpos-irange[k], 1)
                                    r2 <- ifelse(kpos + irange[k] < length(d),
                                                 kpos + irange[k], length(d))
                                    inti[k] <- sum(d[r1:r2])
                                }
                                maxpi <- which.max(inti)
                                if (length(maxpi) > 1) {
                                    m <- wCoefs[opp[maxpi], maxpi]
                                    bestcol <- which(m == max(m),
                                                     arr.ind = TRUE)[2]
                                    best.scale.nr <- maxpi[bestcol]
                                } else  best.scale.nr <- maxpi
                            }

                            best.scale <-  scales[best.scale.nr]
                            best.scale.pos <- opp[best.scale.nr]

                            pprange <- min(pp):max(pp)
                            ## maxint <- max(d[pprange])
                            lwpos <- max(1,best.scale.pos - best.scale)
                            rwpos <- min(best.scale.pos + best.scale, length(td))
                            p1 <- match(td[lwpos], otd)[1]
                            p2 <- match(td[rwpos], otd)
                            p2 <- p2[length(p2)]
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
                                if (length(mz.value) >= (minCentroids + 1)) {
                                    dppm <- round(min(running(abs(diff(mz.value)) /
                                                              (mzrange[2] *  1e-6),
                                                              fun = max,
                                                              width = minCentroids)))
                                } else {
                                    dppm <- round((mzrange[2] - mzrange[1]) /
                                                  (mzrange[2] * 1e-6))
                                }
                            }
                            peaks <- rbind(peaks,
                                           c(mzmean,mzrange, ## mz
                                             NA, NA, NA,     ## rt, rtmin, rtmax,
                                             NA,             ## intensity (sum)
                                             NA,             ## intensity (-bl)
                                             maxint,         ## max intensity
                                             round((maxint - baseline) / sdnoise),  ##  S/N Ratio
                                             NA,             ## Gaussian RMSE
                                             NA,NA,NA,       ## Gaussian Parameters
                                             f,              ## ROI Position
                                             dppm,           ## max. difference between the [minCentroids] peaks in ppm
                                             best.scale,     ## Scale
                                             td[best.scale.pos],
                                             td[lwpos],
                                             td[rwpos],  ## Peak positions guessed from the wavelet's (scan nr)
                                             NA, NA))                    ## Peak limits (scan nr)
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
            colnames(peakinfo) <- c("scale", "scaleNr", "scpos",
                                    "scmin", "scmax")
            for (p in 1:dim(peaks)[1]) {
                ## find minima, assign rt and intensity values
                if (integrate == 1) {
                    lm <- descendMin(wCoefs[, peakinfo[p,"scaleNr"]],
                                     istart = peakinfo[p,"scpos"])
                    gap <- all(d[lm[1]:lm[2]] == 0) ## looks like we got stuck in a gap right in the middle of the peak
                    if ((lm[1] == lm[2]) || gap )## fall-back
                        lm <- descendMinTol(d,
                                            startpos = c(peakinfo[p, "scmin"],
                                                         peakinfo[p, "scmax"]),
                                            maxDescOutlier)
                } else {
                    lm <- descendMinTol(d, startpos = c(peakinfo[p, "scmin"],
                                                        peakinfo[p, "scmax"]),
                                        maxDescOutlier)
                }
                ## narrow down peak rt boundaries by skipping zeros
                pd <- d[lm[1]:lm[2]]
                np <- length(pd)
                lm.l <-  findEqualGreaterUnsorted(pd, 1)
                lm.l <- max(1, lm.l - 1)
                lm.r <- findEqualGreaterUnsorted(rev(pd), 1)
                lm.r <- max(1, lm.r - 1)
                lm <- lm + c(lm.l - 1, -(lm.r - 1) )

                peakrange <- td[lm]
                peaks[p, "rtmin"] <- scantime[peakrange[1]]
                peaks[p, "rtmax"] <- scantime[peakrange[2]]
                peaks[p, "maxo"] <- max(d[lm[1]:lm[2]])
                pwid <- (scantime[peakrange[2]] - scantime[peakrange[1]]) /
                    (peakrange[2] - peakrange[1])
                if (is.na(pwid))
                    pwid <- 1
                peaks[p, "into"] <- pwid * sum(d[lm[1]:lm[2]])
                db <-  d[lm[1]:lm[2]] - baseline
                peaks[p, "intb"] <- pwid * sum(db[db>0])
                peaks[p, "lmin"] <- lm[1]
                peaks[p, "lmax"] <- lm[2]

                if (fitgauss) {
                    ## perform gaussian fits, use wavelets for inital parameters
                    md <- max(d[lm[1]:lm[2]])
                    d1 <- d[lm[1]:lm[2]] / md ## normalize data for gaussian error calc.
                    pgauss <- fitGauss(td[lm[1]:lm[2]], d[lm[1]:lm[2]],
                                       pgauss = list(mu = peaks[p, "scpos"],
                                                     sigma = peaks[p, "scmax"] -
                                                         peaks[p, "scmin"],
                                                     h = peaks[p, "maxo"]))
                    rtime <- peaks[p, "scpos"]
                    if (!any(is.na(pgauss)) && all(pgauss > 0)) {
                        gtime <- td[match(round(pgauss$mu), td)]
                        if (!is.na(gtime)) {
                            rtime <- gtime
                            peaks[p, "mu"] <- pgauss$mu
                            peaks[p, "sigma"] <- pgauss$sigma
                            peaks[p, "h"] <- pgauss$h
                            peaks[p,"egauss"] <- sqrt((1 / length(td[lm[1]:lm[2]])) *
                                                      sum(((d1-gauss(td[lm[1]:lm[2]],
                                                                     pgauss$h / md,
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

        ## BEGIN - plotting/sleep
        if ((sleep >0) && (!is.null(peaks))) {
            tdp <- scantime[td]; trange <- range(tdp)
            egauss <- paste(round(peaks[,"egauss"],3),collapse=", ")
            cdppm <- paste(peaks[,"dppm"],collapse=", ")
            csn <- paste(peaks[,"sn"],collapse=", ")
            par(bg = "white")
            l <- layout(matrix(c(1,2,3),nrow=3,ncol=1,byrow=T),heights=c(.5,.75,2));
            par(mar= c(2, 4, 4, 2) + 0.1)
            ## plotRaw(object,mzrange=mzrange,rtrange=trange,log=TRUE,title='')
            ## Do plotRaw manually.
            raw_mat <- .rawMat(mz = mz, int = int, scantime = scantime,
                               valsPerSpect = valsPerSpect, mzrange = mzrange,
                               rtrange = rtrange, log = TRUE)
            if (nrow(raw_mat) > 0) {
                y <- raw_mat[, "intensity"]
                ylim <- range(y)
                y <- y / ylim[2]
                colorlut <- terrain.colors(16)
                col <- colorlut[y * 15 + 1]
                plot(raw_mat[, "time"], raw_mat[, "mz"], pch = 20, cex = .5,
                     main = "", xlab = "Seconds", ylab = "m/z", col = col,
                     xlim = trange)
            } else {
                plot(c(NA, NA), main = "", xlab = "Seconds", ylab = "m/z",
                     xlim = trange, ylim = mzrange)
            }
            ## done
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
        ## -- END plotting/sleep

        if (!is.null(peaks)) {
            peaklist[[length(peaklist) + 1]] <- peaks
        }
    } ## f

    if (length(peaklist) == 0) {
        warning("No peaks found!")

        if (verboseColumns) {
            nopeaks <- matrix(nrow = 0, ncol = length(basenames) +
                                            length(verbosenames))
            colnames(nopeaks) <- c(basenames, verbosenames)
        } else {
            nopeaks <- matrix(nrow = 0, ncol = length(basenames))
            colnames(nopeaks) <- c(basenames)
        }
        message(" FAIL: none found!")
        return(nopeaks)
    }
    p <- do.call(rbind, peaklist)
    if (!verboseColumns)
        p <- p[, basenames, drop = FALSE]

    uorder <- order(p[, "into"], decreasing = TRUE)
    pm <- as.matrix(p[,c("mzmin", "mzmax", "rtmin", "rtmax"), drop = FALSE])
    uindex <- rectUnique(pm, uorder, mzdiff, ydiff = -0.00001) ## allow adjacent peaks
    pr <- p[uindex, , drop = FALSE]
    message(" OK: ", nrow(pr), " found.")

    return(pr)
}
## This version fixes issue #135, i.e. that the peak signal is integrated based
## on the mzrange of the ROI and not of the actually reported peak.
## Issue #136.
##
## What's different to the original version?
##
## 1) The mz range of the peaks is calculated only using mz values with a
##    measured intensity. This avoids mz ranges from 0 to max mz of the peak,
##    with the mz=0 corresponding actually to scans in which no intensity was
##    measured. Search for "@MOD1" to jump to the respective code.
## 
## 2) The intensities for the peak are reloaded with the refined mz range during
##    the postprocessing. Search for "@MOD2" to jump to the respective code.
##
## What I don't like:
## o Might be better if the getEIC and getMZ C functions returned NA instead of 0
##   if nothing was measured.
## o The joinOverlappingPeaks is still calculated using the variable "d" which
##   contains all intensities from the ROI - might actually not be too bad
##   though.
.centWave_new <- function(mz, int, scantime, valsPerSpect,
                          ppm = 25, peakwidth = c(20,50), snthresh = 10,
                          prefilter = c(3,100), mzCenterFun = "wMean",
                          integrate = 1, mzdiff = -0.001, fitgauss = FALSE,
                          noise = 0, ## noise.local=TRUE,
                          sleep = 0, verboseColumns = FALSE, roiList = list(),
                          firstBaselineCheck = TRUE, roiScales = NULL) {
    if (sleep)
        warning("Parameter 'sleep' is defunct")
    ## TODO @jo Ensure in upstream method that data is in centroided mode!
    ## TODO @jo Ensure the upstream method did eventual sub-setting on scanrange
    ## Input argument checking.
    if (missing(mz) | missing(int) | missing(scantime) | missing(valsPerSpect))
        stop("Arguments 'mz', 'int', 'scantime' and 'valsPerSpect'",
             " are required!")
    if (length(mz) != length(int) | length(valsPerSpect) != length(scantime)
        | length(mz) != sum(valsPerSpect))
        stop("Lengths of 'mz', 'int' and of 'scantime','valsPerSpect'",
             " have to match. Also, 'length(mz)' should be equal to",
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
        stop("Function '", mzCenterFun, "' not defined !")

    if (!is.logical(firstBaselineCheck))
        stop("Parameter 'firstBaselineCheck' should be logical!")
    if (length(firstBaselineCheck) != 1)
        stop("Parameter 'firstBaselineCheck' should be a single logical !")
    if (length(roiScales) > 0)
        if (length(roiScales) != length(roiList) | !is.numeric(roiScales))
            stop("If provided, parameter 'roiScales' has to be a numeric with",
                 " length equal to the length of 'roiList'!")

    basenames <- c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax",
                   "into", "intb", "maxo", "sn")
    verbosenames <- c("egauss", "mu", "sigma", "h", "f", "dppm", "scale",
                      "scpos", "scmin", "scmax", "lmin", "lmax")

    ## Peak width: seconds to scales
    scalerange <- round((peakwidth / mean(diff(scantime))) / 2)

    if (length(z <- which(scalerange == 0)))
        scalerange <- scalerange[-z]
    if (length(scalerange) < 1) {
        warning("No scales? Please check peak width!")
        if (verboseColumns) {
            nopeaks <- matrix(nrow = 0, ncol = length(basenames) +
                                            length(verbosenames))
            colnames(nopeaks) <- c(basenames, verbosenames)
        } else {
            nopeaks <- matrix(nrow = 0, ncol = length(basenames))
            colnames(nopeaks) <- c(basenames)
        }
        return(invisible(nopeaks))
    }

    if (length(scalerange) > 1)
        scales <- seq(from = scalerange[1], to = scalerange[2], by = 2)
    else
        scales <- scalerange

    minPeakWidth <-  scales[1]
    noiserange <- c(minPeakWidth * 3, max(scales) * 3)
    maxGaussOverlap <- 0.5
    minPtsAboveBaseLine <- max(4, minPeakWidth - 2)
    minCentroids <- minPtsAboveBaseLine
    scRangeTol <-  maxDescOutlier <- floor(minPeakWidth / 2)
    scanrange <- c(1, length(scantime))

    ## If no ROIs are supplied then search for them.
    if (length(roiList) == 0) {
        message("Detecting mass traces at ", ppm, " ppm ... ", appendLF = FALSE)
        ## flush.console();
        ## We're including the findmzROI code in this function to reduce
        ## the need to copy objects etc.
        ## We could also sort the data by m/z anyway; wouldn't need that
        ## much time. Once we're using classes from MSnbase we can be
        ## sure that values are correctly sorted.
        withRestarts(
            tryCatch({
                tmp <- capture.output(
                    roiList <- .Call("findmzROI",
                                     mz, int, scanindex,
                                     as.double(c(0.0, 0.0)),
                                     as.integer(scanrange),
                                     as.integer(length(scantime)),
                                     as.double(ppm * 1e-6),
                                     as.integer(minCentroids),
                                     as.integer(prefilter),
                                     as.integer(noise),
                                     PACKAGE ='xcms' )
                )
            },
            error = function(e){
                if (grepl("m/z sort assumption violated !", e$message)) {
                    invokeRestart("fixSort")
                } else {
                    simpleError(e)
                }
            }),
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
                tmp <- capture.output(
                    roiList <<- .Call("findmzROI",
                                      mz, int, scanindex,
                                      as.double(c(0.0, 0.0)),
                                      as.integer(scanrange),
                                      as.integer(length(scantime)),
                                      as.double(ppm * 1e-6),
                                      as.integer(minCentroids),
                                      as.integer(prefilter),
                                      as.integer(noise),
                                      PACKAGE ='xcms' )
                )
            }
        )
        message("OK")
        if (length(roiList) == 0) {
            warning("No ROIs found! \n")
            if (verboseColumns) {
                nopeaks <- matrix(nrow = 0, ncol = length(basenames) +
                                                length(verbosenames))
                colnames(nopeaks) <- c(basenames, verbosenames)
            } else {
                nopeaks <- matrix(nrow = 0, ncol = length(basenames))
                colnames(nopeaks) <- c(basenames)
            }
            return(invisible(nopeaks))
        }
    }

    ## Second stage: process the ROIs
    peaklist <- list()
    Nscantime <- length(scantime)
    lf <- length(roiList)

    ## cat('\n Detecting chromatographic peaks ... \n % finished: ')
    ## lp <- -1
    message("Detecting chromatographic peaks in ", length(roiList),
            " regions of interest ...", appendLF = FALSE)

    for (f in  1:lf) {

        ## cat("\nProcess roi ", f, "\n")
        feat <- roiList[[f]]
        N <- feat$scmax - feat$scmin + 1
        peaks <- peakinfo <- NULL
        mzrange <- c(feat$mzmin, feat$mzmax)
        mzrange_ROI <- mzrange
        sccenter <- feat$scmin[1] + floor(N/2) - 1
        scrange <- c(feat$scmin, feat$scmax)
        ## scrange + noiserange, used for baseline detection and wavelet analysis
        sr <- c(max(scanrange[1], scrange[1] - max(noiserange)),
                min(scanrange[2], scrange[2] + max(noiserange)))
        eic <- .Call("getEIC", mz, int, scanindex, as.double(mzrange),
                     as.integer(sr), as.integer(length(scanindex)),
                     PACKAGE = "xcms")
        d <- eic$intensity
        td <- sr[1]:sr[2]
        scan.range <- c(sr[1], sr[2])
        ## original mzROI range
        idxs <- which(eic$scan %in% seq(scrange[1], scrange[2]))
        mzROI.EIC <- list(scan=eic$scan[idxs], intensity=eic$intensity[idxs])
        omz <- .Call("getMZ", mz, int, scanindex, as.double(mzrange),
                     as.integer(scrange), as.integer(length(scantime)),
                     PACKAGE = 'xcms')
        if (all(omz == 0)) {
            warning("centWave: no peaks found in ROI.")
            next
        }
        od  <- mzROI.EIC$intensity
        otd <- mzROI.EIC$scan
        if (all(od == 0)) {
            warning("centWave: no peaks found in ROI.")
            next
        }
        ## scrange + scRangeTol, used for gauss fitting and continuous
        ## data above 1st baseline detection
        ftd <- max(td[1], scrange[1] - scRangeTol) : min(td[length(td)],
                                                         scrange[2] + scRangeTol)
        fd <- d[match(ftd, td)]

        ## 1st type of baseline: statistic approach
        if (N >= 10*minPeakWidth) {
            ## in case of very long mass trace use full scan range
            ## for baseline detection
            noised <- .Call("getEIC", mz, int, scanindex, as.double(mzrange),
                            as.integer(scanrange), as.integer(length(scanindex)),
                            PACKAGE="xcms")$intensity
        } else {
            noised <- d
        }
        ## 90% trimmed mean as first baseline guess
        noise <- estimateChromNoise(noised, trim = 0.05,
                                    minPts = 3 * minPeakWidth)
        ## any continuous data above 1st baseline ?
        if (firstBaselineCheck &
            !continuousPtsAboveThreshold(fd, threshold = noise,
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
            wCoefs[nrow(wCoefs),] <- wCoefs[nrow(wCoefs) - 1, ] * 0.99
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
                            ## if(!is.null(roiScales)) {
                            ## allow roiScales to be a numeric of length 0
                            if(length(roiScales) > 0) {
                                ## use given scale
                                best.scale.nr <- which(scales == roiScales[[f]])
                                if(best.scale.nr > length(opp))
                                    best.scale.nr <- length(opp)
                            } else {
                                ## try to decide which scale describes the peak best
                                inti <- numeric(length(opp))
                                irange <- rep(ceiling(scales[1]/2), length(opp))
                                for (k in 1:length(opp)) {
                                    kpos <- opp[k]
                                    r1 <- ifelse(kpos - irange[k] > 1,
                                                 kpos-irange[k], 1)
                                    r2 <- ifelse(kpos + irange[k] < length(d),
                                                 kpos + irange[k], length(d))
                                    inti[k] <- sum(d[r1:r2])
                                }
                                maxpi <- which.max(inti)
                                if (length(maxpi) > 1) {
                                    m <- wCoefs[opp[maxpi], maxpi]
                                    bestcol <- which(m == max(m),
                                                     arr.ind = TRUE)[2]
                                    best.scale.nr <- maxpi[bestcol]
                                } else  best.scale.nr <- maxpi
                            }

                            best.scale <-  scales[best.scale.nr]
                            best.scale.pos <- opp[best.scale.nr]

                            pprange <- min(pp):max(pp)
                            ## maxint <- max(d[pprange])
                            lwpos <- max(1,best.scale.pos - best.scale)
                            rwpos <- min(best.scale.pos + best.scale, length(td))
                            p1 <- match(td[lwpos], otd)[1]
                            p2 <- match(td[rwpos], otd)
                            p2 <- p2[length(p2)]
                            ## cat("p1: ", p1, " p2: ", p2, "\n")
                            if (is.na(p1)) p1 <- 1
                            if (is.na(p2)) p2 <- N
                            mz.value <- omz[p1:p2]
                            ## cat("mz.value: ", paste0(mz.value, collapse = ", "),
                            ##     "\n")
                            mz.int <- od[p1:p2]
                            maxint <- max(mz.int)
                            ## @MOD1: Remove mz values for which no intensity was
                            ## measured. Would be better if getEIC returned NA
                            ## if nothing was measured.
                            mzorig <- mz.value
                            mz.value <- mz.value[mz.int > 0]
                            mz.int <- mz.int[mz.int > 0]
                            ## Call next to avoid reporting peaks without mz
                            ## values (issue #165).
                            if (length(mz.value) == 0)
                                next
                            ## cat("mz.value: ", paste0(mz.value, collapse = ", "),
                            ##     "\n")

                            ## re-calculate m/z value for peak range
                            ## cat("mzrange refined: [",
                            ##     paste0(mzrange, collapse = ", "), "]")
                            ## hm, shouldn't we get rid of the mz = 0 here?
                            mzrange <- range(mz.value)
                            ## cat(" -> [",
                            ##     paste0(mzrange, collapse = ", "), "]\n")
                            mzmean <- do.call(mzCenterFun,
                                              list(mz = mz.value,
                                                   intensity = mz.int))

                            ## Compute dppm only if needed
                            dppm <- NA
                            if (verboseColumns) {
                                if (length(mz.value) >= (minCentroids + 1)) {
                                    dppm <- round(min(running(abs(diff(mz.value)) /
                                                              (mzrange[2] *  1e-6),
                                                              fun = max,
                                                              width = minCentroids)))
                                } else {
                                    dppm <- round((mzrange[2] - mzrange[1]) /
                                                  (mzrange[2] * 1e-6))
                                }
                            }
                            peaks <- rbind(peaks,
                                           c(mzmean,mzrange, ## mz
                                             NA, NA, NA,     ## rt, rtmin, rtmax,
                                             NA,             ## intensity (sum)
                                             NA,             ## intensity (-bl)
                                             maxint,         ## max intensity
                                             round((maxint - baseline) / sdnoise),  ##  S/N Ratio
                                             NA,             ## Gaussian RMSE
                                             NA,NA,NA,       ## Gaussian Parameters
                                             f,              ## ROI Position
                                             dppm,           ## max. difference between the [minCentroids] peaks in ppm
                                             best.scale,     ## Scale
                                             td[best.scale.pos],
                                             td[lwpos],
                                             td[rwpos],  ## Peak positions guessed from the wavelet's (scan nr)
                                             NA, NA))                    ## Peak limits (scan nr)
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
            colnames(peakinfo) <- c("scale", "scaleNr", "scpos",
                                    "scmin", "scmax")
            for (p in 1:dim(peaks)[1]) {
                ## @MOD2
                ## Fix for issue #135: reload the EIC data if the
                ## mzrange differs from that of the ROI, but only if the mz
                ## range of the peak is different from the one of the ROI.
                mzr <- peaks[p, c("mzmin", "mzmax")]
                if (any(mzr != mzrange_ROI)) {
                    eic <- .Call("getEIC", mz, int, scanindex,
                                 as.double(mzr), as.integer(sr),
                                 as.integer(length(scanindex)),
                                 PACKAGE = "xcms")
                    current_ints <- eic$intensity
                    ## Force re-loading also of a potential additional peak in
                    ## the same ROI.
                    mzrange_ROI <- c(0, 0)
                } else {
                    current_ints <- d
                }
                ## find minima, assign rt and intensity values
                if (integrate == 1) {
                    lm <- descendMin(wCoefs[, peakinfo[p,"scaleNr"]],
                                     istart = peakinfo[p,"scpos"])
                    gap <- all(current_ints[lm[1]:lm[2]] == 0) ## looks like we got stuck in a gap right in the middle of the peak
                    if ((lm[1] == lm[2]) || gap )## fall-back
                        lm <- descendMinTol(current_ints,
                                            startpos = c(peakinfo[p, "scmin"],
                                                         peakinfo[p, "scmax"]),
                                            maxDescOutlier)
                } else {
                    lm <- descendMinTol(current_ints,
                                        startpos = c(peakinfo[p, "scmin"],
                                                     peakinfo[p, "scmax"]),
                                        maxDescOutlier)
                }
                ## narrow down peak rt boundaries by skipping zeros
                pd <- current_ints[lm[1]:lm[2]]
                np <- length(pd)
                lm.l <-  findEqualGreaterUnsorted(pd, 1)
                lm.l <- max(1, lm.l - 1)
                lm.r <- findEqualGreaterUnsorted(rev(pd), 1)
                lm.r <- max(1, lm.r - 1)
                lm <- lm + c(lm.l - 1, -(lm.r - 1) )

                peakrange <- td[lm]
                peaks[p, "rtmin"] <- scantime[peakrange[1]]
                peaks[p, "rtmax"] <- scantime[peakrange[2]]
                peaks[p, "maxo"] <- max(current_ints[lm[1]:lm[2]])
                pwid <- (scantime[peakrange[2]] - scantime[peakrange[1]]) /
                    (peakrange[2] - peakrange[1])
                if (is.na(pwid))
                    pwid <- 1
                peaks[p, "into"] <- pwid * sum(current_ints[lm[1]:lm[2]])
                db <-  current_ints[lm[1]:lm[2]] - baseline
                peaks[p, "intb"] <- pwid * sum(db[db>0])
                peaks[p, "lmin"] <- lm[1]
                peaks[p, "lmax"] <- lm[2]
                ## cat("[", paste0(peaks[p, c("rtmin", "rtmax")], collapse = ", "),
                ##     "] into ", peaks[p, "into"], "\n")

                if (fitgauss) {
                    ## perform gaussian fits, use wavelets for inital parameters
                    md <- max(current_ints[lm[1]:lm[2]])
                    d1 <- current_ints[lm[1]:lm[2]] / md ## normalize data for gaussian error calc.
                    pgauss <- fitGauss(td[lm[1]:lm[2]], current_ints[lm[1]:lm[2]],
                                       pgauss = list(mu = peaks[p, "scpos"],
                                                     sigma = peaks[p, "scmax"] -
                                                         peaks[p, "scmin"],
                                                     h = peaks[p, "maxo"]))
                    rtime <- peaks[p, "scpos"]
                    if (!any(is.na(pgauss)) && all(pgauss > 0)) {
                        gtime <- td[match(round(pgauss$mu), td)]
                        if (!is.na(gtime)) {
                            rtime <- gtime
                            peaks[p, "mu"] <- pgauss$mu
                            peaks[p, "sigma"] <- pgauss$sigma
                            peaks[p, "h"] <- pgauss$h
                            peaks[p,"egauss"] <- sqrt((1 / length(td[lm[1]:lm[2]])) *
                                                      sum(((d1-gauss(td[lm[1]:lm[2]],
                                                                     pgauss$h / md,
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
            ## Use d here instead of current_ints
            peaks <- joinOverlappingPeaks(td, d, otd, omz, od,
                                          scantime, scan.range, peaks,
                                          maxGaussOverlap,
                                          mzCenterFun = mzCenterFun)
        }
        if (!is.null(peaks)) {
            peaklist[[length(peaklist) + 1]] <- peaks
        }
    } ## f

    if (length(peaklist) == 0) {
        warning("No peaks found!")

        if (verboseColumns) {
            nopeaks <- matrix(nrow = 0, ncol = length(basenames) +
                                            length(verbosenames))
            colnames(nopeaks) <- c(basenames, verbosenames)
        } else {
            nopeaks <- matrix(nrow = 0, ncol = length(basenames))
            colnames(nopeaks) <- c(basenames)
        }
        message(" FAIL: none found!")
        return(nopeaks)
    }

    ## cat("length peaklist: ", length(peaklist), "\n")
    p <- do.call(rbind, peaklist)
    if (!verboseColumns)
        p <- p[, basenames, drop = FALSE]
    return(p)
    
    uorder <- order(p[, "into"], decreasing = TRUE)
    pm <- as.matrix(p[,c("mzmin", "mzmax", "rtmin", "rtmax"), drop = FALSE])
    uindex <- rectUnique(pm, uorder, mzdiff, ydiff = -0.00001) ## allow adjacent peaks
    pr <- p[uindex, , drop = FALSE]
    message(" OK: ", nrow(pr), " found.")

    return(pr)
}




############################################################
## massifquant
##
#' @title Core API function for massifquant peak detection
#'
#' @description Massifquant is a Kalman filter (KF)-based chromatographic peak
#'     detection for XC-MS data in centroid mode. The identified peaks
#'     can be further refined with the \emph{centWave} method (see
#'     \code{\link{do_findChromPeaks_centWave}} for details on centWave)
#'     by specifying \code{withWave = TRUE}.
#'
#' @details This algorithm's performance has been tested rigorously
#'     on high resolution LC/{OrbiTrap, TOF}-MS data in centroid mode.
#'     Simultaneous kalman filters identify peaks and calculate their
#'     area under the curve. The default parameters are set to operate on
#'     a complex LC-MS Orbitrap sample. Users will find it useful to do some
#'     simple exploratory data analysis to find out where to set a minimum
#'     intensity, and identify how many scans an average peak spans. The
#'     \code{consecMissedLimit} parameter has yielded good performance on
#'     Orbitrap data when set to (\code{2}) and on TOF data it was found best
#'     to be at (\code{1}). This may change as the algorithm has yet to be
#'     tested on many samples. The \code{criticalValue} parameter is perhaps
#'     most dificult to dial in appropriately and visual inspection of peak
#'     identification is the best suggested tool for quick optimization.
#'     The \code{ppm} and \code{checkBack} parameters have shown less influence
#'     than the other parameters and exist to give users flexibility and
#'     better accuracy.
#' 
#' @inheritParams do_findChromPeaks_centWave
#'
#' @inheritParams findChromPeaks-centWave
#'
#' @inheritParams findChromPeaks-massifquant
#'
#' @return
#' A matrix, each row representing an identified chromatographic peak,
#'     with columns:
#'     \describe{
#'     \item{mz}{Intensity weighted mean of m/z values of the peaks across
#'     scans.}
#'     \item{mzmin}{Minumum m/z of the peak.}
#'     \item{mzmax}{Maximum m/z of the peak.}
#'     \item{rtmin}{Minimum retention time of the peak.}
#'     \item{rtmax}{Maximum retention time of the peak.}
#'     \item{rt}{Retention time of the peak's midpoint.}
#'     \item{into}{Integrated (original) intensity of the peak.}
#'     \item{maxo}{Maximum intensity of the peak.}
#'     }
#' 
#'     If \code{withWave} is set to \code{TRUE}, the result is the same as
#'     returned by the \code{\link{do_findChromPeaks_centWave}} method.
#' 
#' @family core peak detection functions
#' 
#' @seealso \code{\link{massifquant}} for the standard user interface method.
#'
#' @references
#' Conley CJ, Smith R, Torgrip RJ, Taylor RM, Tautenhahn R and Prince JT
#' "Massifquant: open-source Kalman filter-based XC-MS isotope trace feature
#' detection" \emph{Bioinformatics} 2014, 30(18):2636-43.
#'
#' @author Christopher Conley
#'
#' @examples
#' library(faahKO)
#' library(xcms)
#' cdfpath <- system.file("cdf", package = "faahKO")
#' cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)
#'
#' ## Read the first file
#' xraw <- xcmsRaw(cdffiles[1])
#' ## Extract the required data
#' mzVals <- xraw@env$mz
#' intVals <- xraw@env$intensity
#' ## Define the values per spectrum:
#' valsPerSpect <- diff(c(xraw@scanindex, length(mzVals)))
#'
#' ## Perform the peak detection using massifquant
#' res <- do_findChromPeaks_massifquant(mz = mzVals, int = intVals,
#' scantime = xraw@scantime, valsPerSpect = valsPerSpect)
#' head(res)
do_findChromPeaks_massifquant <- function(mz,
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
                                          withWave = FALSE) {
    message("\n Massifquant, Copyright (C) 2013 Brigham Young University.")
    message(" Massifquant comes with ABSOLUTELY NO WARRANTY.",
            " See LICENSE for details.", sep ="")
    ## flush.console()

    ## TODO @jo Ensure in upstream method that data is in centroided mode!
    ## TODO @jo Ensure the upstream method did eventual sub-setting on scanrange
    ## Input argument checking.
    if (missing(mz) | missing(int) | missing(scantime) | missing(valsPerSpect))
        stop("Arguments 'mz', 'int', 'scantime' and 'valsPerSpect'",
             " are required!")
    if ((length(mz) != length(int)) | (length(valsPerSpect) != length(scantime))
        | (length(mz) != sum(valsPerSpect)))
        stop("Lengths of 'mz', 'int' and of 'scantime','valsPerSpect'",
             " have to match. Also, 'length(mz)' should be equal to",
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

    message("\n Detecting  mass traces at ",ppm,"ppm ... ", appendLF = FALSE)
    flush.console()
    massifquantROIs <- do_findKalmanROI(mz = mz, int = int, scantime = scantime,
                                        valsPerSpect = valsPerSpect,
                                        minIntensity = prefilter[2],
                                        minCentroids = peakwidth[1],
                                        criticalVal = criticalValue,
                                        consecMissedLim = consecMissedLimit,
                                        segs = unions, scanBack = checkBack,
                                        ppm = ppm)
    message("OK")
    if (withWave) {
        featlist <- do_findChromPeaks_centWave(mz = mz, int = int,
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
                                               roiList = massifquantROIs)
    }
    else {
        ## Get index vector for C calls
        scanindex <- valueCount2ScanIndex(valsPerSpect)
        basenames <- c("mz","mzmin","mzmax","rtmin","rtmax","rt", "into")
        if (length(massifquantROIs) == 0) {
            warning("\nNo peaks found!")
            nopeaks <- matrix(nrow=0, ncol=length(basenames))
            colnames(nopeaks) <- basenames
            return(nopeaks)
        }

        ## Get the max intensity for each peak.
        maxo <- lapply(massifquantROIs, function(z) {
            raw <- .rawMat(mz = mz, int = int, scantime = scantime,
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
        message(" ", dim(featlist)[1]," Peaks.");
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
#' @title Core API function for matchedFilter peak detection
#'
#' @description This function identifies peaks in the chromatographic
#'     time domain as described in [Smith 2006]. The intensity values are
#'     binned by cutting The LC/MS data into slices (bins) of a mass unit
#'     (\code{binSize} m/z) wide. Within each bin the maximal intensity is
#'     selected. The peak detection is then performed in each bin by
#'     extending it based on the \code{steps} parameter to generate slices
#'     comprising bins \code{current_bin - steps +1} to
#'     \code{current_bin + steps - 1}.
#'     Each of these slices is then filtered with matched filtration using
#'     a second-derative Gaussian as the model peak shape. After filtration
#'     peaks are detected using a signal-to-ration cut-off. For more details
#'     and illustrations see [Smith 2006].
#'
#' @details The intensities are binned by the provided m/z values within each
#'     spectrum (scan). Binning is performed such that the bins are centered
#'     around the m/z values (i.e. the first bin includes all m/z values between
#'     \code{min(mz) - bin_size/2} and \code{min(mz) + bin_size/2}).
#'
#'     For more details on binning and missing value imputation see
#'     \code{\link{binYonX}} and \code{\link{imputeLinInterpol}} methods.
#'
#' @note This function exposes core peak detection functionality of
#'     the \emph{matchedFilter} method. While this function can be called
#'     directly, users will generally call the corresponding method for the
#'     data object instead (e.g. the \code{link{findPeaks.matchedFilter}}
#'     method).
#'
#' @inheritParams do_findChromPeaks_centWave
#' 
#' @inheritParams findChromPeaks-centWave
#' 
#' @inheritParams imputeLinInterpol
#' 
#' @inheritParams findChromPeaks-matchedFilter
#'
#' @return A matrix, each row representing an identified chromatographic peak,
#'     with columns:
#'     \describe{
#'     \item{mz}{Intensity weighted mean of m/z values of the peak across scans.}
#'     \item{mzmin}{Minimum m/z of the peak.}
#'     \item{mzmax}{Maximum m/z of the peak.}
#'     \item{rt}{Retention time of the peak's midpoint.}
#'     \item{rtmin}{Minimum retention time of the peak.}
#'     \item{rtmax}{Maximum retention time of the peak.}
#'     \item{into}{Integrated (original) intensity of the peak.}
#'     \item{intf}{Integrated intensity of the filtered peak.}
#'     \item{maxo}{Maximum intensity of the peak.}
#'     \item{maxf}{Maximum intensity of the filtered peak.}
#'     \item{i}{Rank of peak in merged EIC (\code{<= max}).}
#'     \item{sn}{Signal to noise ratio of the peak}
#'     }
#' 
#' @references
#' Colin A. Smith, Elizabeth J. Want, Grace O'Maille, Ruben Abagyan and
#' Gary Siuzdak. "XCMS: Processing Mass Spectrometry Data for Metabolite
#' Profiling Using Nonlinear Peak Alignment, Matching, and Identification"
#' \emph{Anal. Chem.} 2006, 78:779-787.
#'
#' @author Colin A Smith, Johannes Rainer
#'
#' @family core peak detection functions
#'
#' @seealso \code{\link{binYonX}} for a binning function,
#'     \code{\link{imputeLinInterpol}} for the interpolation of missing values.
#'     \code{\link{matchedFilter}} for the standard user interface method.
#' 
#' @examples
#' ## Load the test file
#' library(faahKO)
#' fs <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
#' xr <- xcmsRaw(fs)
#'
#' ## Extracting the data from the xcmsRaw for do_findChromPeaks_centWave
#' mzVals <- xr@env$mz
#' intVals <- xr@env$intensity
#' ## Define the values per spectrum:
#' valsPerSpect <- diff(c(xr@scanindex, length(mzVals)))
#'
#' res <- do_findChromPeaks_matchedFilter(mz = mzVals, int = intVals,
#' scantime = xr@scantime, valsPerSpect = valsPerSpect)
#' head(res)
do_findChromPeaks_matchedFilter <- function(mz,
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
                                            index = FALSE,
                                            sleep = 0
                                            ){
    ## Use original code
    if (useOriginalCode()) {
        ## warning("Old xcms code was used; be aware that this code",
        ##         " may contain bugs.")
        return(.matchedFilter_orig(mz, int, scantime, valsPerSpect,
                                   binSize, impute, baseValue, distance,
                                   fwhm, sigma, max, snthresh,
                                   steps, mzdiff, index, sleep = sleep))
    } else {
        return(.matchedFilter_binYonX_no_iter(mz, int, scantime, valsPerSpect,
                                              binSize, impute, baseValue,
                                              distance, fwhm, sigma, max,
                                              snthresh, steps, mzdiff, index,
                                              sleep = sleep
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
                                index = FALSE,
                                sleep = 0
                                ){
    .Deprecated(msg = paste0("Use of the original code with iterative binning",
                             " is discouraged!"))
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
             " have to match. Also, 'length(mz)' should be equal to",
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
        baseValue <- numeric()
    if (length(baseValue) == 0)
        baseValue <- min(int, na.rm = TRUE) / 2
    profp$baselevel <- baseValue
    if (missing(distance))
        distance <- numeric()
    if (length(distance) != 0) {
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

                ## -- begin sleep/plot
                if (sleep > 0) {
                    plot(scantime, yfilt, type = "l",
                         main = paste(mass[i], "-", mass[i+1]),
                         ylim = c(-gmax/3, gmax))
                    points(cbind(scantime, yfilt)[peakrange[1]:peakrange[2],],
                           type = "l", col = "red")
                    points(scantime, colSums(ymat), type = "l", col = "blue",
                           lty = "dashed")
                    abline(h = snthresh*noise, col = "red")
                    Sys.sleep(sleep)
                }
                ## -- end sleep plot
                
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
                                           index = FALSE,
                                           sleep = 0
                                           ){
    ## Input argument checking.
    if (missing(mz) | missing(int) | missing(scantime) | missing(valsPerSpect))
        stop("Arguments 'mz', 'int', 'scantime' and 'valsPerSpect'",
             " are required!")
    if (length(mz) != length(int) | length(valsPerSpect) != length(scantime)
        | length(mz) != sum(valsPerSpect))
        stop("Lengths of 'mz', 'int' and of 'scantime','valsPerSpect'",
             " have to match. Also, 'length(mz)' should be equal to",
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
        ## Calculate the "scanindex" from the number of values per spectrum:
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
            baseValue <- numeric()
        if (length(baseValue) == 0)
            baseValue <- min(int, na.rm = TRUE) / 2
        if (missing(distance))
            distance <- numeric()
        if (length(distance) == 0)
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

                ## begin sleep/plot
                if (sleep > 0) {
                    plot(scantime, yfilt, type = "l",
                         main = paste(mass[i], "-", mass[i+1]),
                         ylim=c(-gmax/3, gmax))
                    points(cbind(scantime, yfilt)[peakrange[1]:peakrange[2],],
                           type = "l", col = "red")
                    points(scantime, colSums(ymat), type = "l", col = "blue",
                           lty = "dashed")
                    abline(h = snthresh*noise, col = "red")
                    Sys.sleep(sleep)
                }
                ## end sleep/plot
                
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
#' @title Core API function for single-spectrum non-chromatography MS data
#'     peak detection
#'
#' @description This function performs peak detection in mass spectrometry
#'     direct injection spectrum using a wavelet based algorithm.
#'
#' @details This is a wrapper around the peak picker in Bioconductor's
#'     \code{MassSpecWavelet} package calling
#'     \code{\link{peakDetectionCWT}} and
#'     \code{\link{tuneInPeakInfo}} functions. See the
#'     \emph{xcmsDirect} vignette for more information.
#'
#' @inheritParams do_findChromPeaks_centWave
#' 
#' @inheritParams findChromPeaks-centWave
#'
#' @param ... Additional parameters to be passed to the
#'     \code{\link{peakDetectionCWT}} function.
#'
#' @return A matrix, each row representing an identified peak, with columns:
#'     \describe{
#'     \item{mz}{m/z value of the peak at the centroid position.}
#'     \item{mzmin}{Minimum m/z of the peak.}
#'     \item{mzmax}{Maximum m/z of the peak.}
#'     \item{rt}{Always \code{-1}.}
#'     \item{rtmin}{Always \code{-1}.}
#'     \item{rtmax}{Always \code{-1}.}
#'     \item{into}{Integrated (original) intensity of the peak.}
#'     \item{maxo}{Maximum intensity of the peak.}
#'     \item{intf}{Always \code{NA}.}
#'     \item{maxf}{Maximum MSW-filter response of the peak.}
#'     \item{sn}{Signal to noise ratio.}
#'     }
#'
#' @family core peak detection functions
#'
#' @seealso \code{\link{MSW}} for the standard user interface
#'     method. \code{\link{peakDetectionCWT}} from the
#'     \code{MassSpecWavelet} package.
#' 
#' @author Joachim Kutzera, Steffen Neumann, Johannes Rainer
do_findPeaks_MSW <- function(mz, int, snthresh = 3,
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

    ## cat('\n')

    ## Filter additional (verbose) columns
    if (!verboseColumns)
        peaklist <- peaklist[,basenames,drop=FALSE]

    peaklist
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
## MS1
## This one might be too cumbersome to do it for plain vectors. It would be ideal
## for MSnExp objects though.
##
## do_findChromPeaks_MS1 <- function(mz, int, scantime, valsPerSpect) {
##     ## Checks: do I have
## }


## ## Original code: TODO REMOVE ME once method is validated.
## do_predictIsotopeROIs <- function(object,
##                                   xcmsPeaks, ppm=25,
##                                   maxcharge=3, maxiso=5, mzIntervalExtension=TRUE) {
##   if(nrow(xcmsPeaks) == 0){
##     warning("Warning: There are no peaks (parameter >xcmsPeaks<) for the prediction of isotope ROIs !\n")
##     return(list())
##   }
##   if(class(xcmsPeaks) != "xcmsPeaks")
##     stop("Error: parameter >xcmsPeaks< is not of class 'xcmsPeaks' ! \n")
##   if(any(is.na(match(x = c("scmin", "scmax"), table = colnames(xcmsPeaks)))))
##     stop("Error: peak list >xcmsPeaks< is missing the columns 'scmin' and 'scmax' ! Please set parameter >verbose.columns< to TRUE for peak picking with 'centWave' and try again ! \n")

##   addNewIsotopeROIs <- TRUE
##   addNewAdductROIs  <- FALSE
##   polarity <- NA

##   ## convert present peaks to list of lists
##   presentROIs.list <- list()
##   for(peakIdx in 1:nrow(xcmsPeaks)){
##     presentROIs.list[[peakIdx]] <- list(
##       mz        = xcmsPeaks[[peakIdx, "mz"]],## XXX not used!
##       mzmin     = xcmsPeaks[[peakIdx, "mzmin"]],
##       mzmax     = xcmsPeaks[[peakIdx, "mzmax"]],
##       scmin     = xcmsPeaks[[peakIdx, "scmin"]],
##       scmax     = xcmsPeaks[[peakIdx, "scmax"]],
##       length    = -1,## XXX not used!
##       intensity = xcmsPeaks[[peakIdx, "intb"]],## XXX not used!
##       scale     = xcmsPeaks[[peakIdx, "scale"]]## XXX not used!
##     )

##     if(abs(xcmsPeaks[[peakIdx, "mzmax"]] - xcmsPeaks[[peakIdx, "mzmin"]]) < xcmsPeaks[[peakIdx, "mz"]] * ppm / 1E6){
##       presentROIs.list[[peakIdx]]$mzmin <- xcmsPeaks[[peakIdx, "mz"]] - xcmsPeaks[[peakIdx, "mz"]] * (ppm/2) / 1E6
##       presentROIs.list[[peakIdx]]$mzmax <- xcmsPeaks[[peakIdx, "mz"]] + xcmsPeaks[[peakIdx, "mz"]] * (ppm/2) / 1E6
##     }
##   }

##   ## fetch predicted ROIs
##   resultObj <- createAdditionalROIs(object, presentROIs.list, ppm, addNewIsotopeROIs, maxcharge, maxiso, mzIntervalExtension, addNewAdductROIs, polarity)
##   newRoiCounter <- resultObj$newRoiCounter
##   numberOfAdditionalIsotopeROIs <- resultObj$numberOfAdditionalIsotopeROIs
##   numberOfAdditionalAdductROIs <- resultObj$numberOfAdditionalAdductROIs
##   newROI.matrix <- resultObj$newROI.matrix

##   if(nrow(newROI.matrix) == 0)
##     return(list())

##   ## This should not be needed, as it has already been performed above.
##   ## remove ROIs with weak signal content
##   intensityThreshold <- 10
##   newROI.matrix <- removeROIsWithoutSignal(object, newROI.matrix, intensityThreshold)

##   ## convert to list of lists
##   newROI.list <- list()
##   for(idx in 1:nrow(newROI.matrix))
##     ## c("mz", "mzmin", "mzmax", "scmin", "scmax", "length", "intensity")
##     newROI.list[[length(newROI.list) + 1]] <- as.list(newROI.matrix[idx, ])

##   cat("Predicted ROIs: ", length(newROI.list), " new ROIs (", numberOfAdditionalIsotopeROIs, " isotope ROIs, ", numberOfAdditionalAdductROIs, " adduct ROIs) for ", length(presentROIs.list)," present ROIs.", "\n")

##   return(newROI.list)
## }

## Tuned from the original code.
#' @param peaks. \code{matrix} or \code{data.frame} with peaks for which
#' isotopes should be predicted. Required columns are \code{"mz"},
#' \code{"mzmin"}, \code{"mzmax"}, \code{"scmin"}, \code{"scmax"},
#' \code{"intb"} and \code{"scale"}.
#'
#' @return a \code{matrix} with columns \code{"mz"}, \code{"mzmin"},
#' \code{"mzmax"}, \code{"scmin"}, \code{"scmax"}, \code{"length"} (always -1),
#' \code{"intensity"} (always -1) and \code{"scale"}.
#' @noRd
do_define_isotopes <- function(peaks., maxCharge = 3, maxIso = 5,
                               mzIntervalExtension = TRUE) {
    req_cols <- c("mz", "mzmin", "mzmax", "scmin", "scmax", "scale")
    if (is.null(dim(peaks.)))
        stop("'peaks.' has to be a matrix or data.frame!")
    if (!all(req_cols %in% colnames(peaks.))) {
        not_there <- req_cols[!(req_cols %in% colnames(peaks.))]
        stop("'peaks.' lacks required columns ",
             paste0("'", not_there, "'", collapse = ","), "!")
    }
    if (is.data.frame(peaks.))
        peaks. <- as.matrix(peaks.)

    isotopeDistance <- 1.0033548378
    charges <- 1:maxCharge
    isos <- 1:maxIso

    isotopePopulationMz <- unique(as.numeric(matrix(isos, ncol = 1) %*%
                                             (isotopeDistance / charges)))
    ## split the peaks into a list.
    roiL <- split(peaks.[, req_cols, drop = FALSE], f = 1:nrow(peaks.))

    newRois <- lapply(roiL, function(z) {
        if (mzIntervalExtension)
            mz_ext <- (z[3] - z[2]) * 2
        else
            mz_ext <- 0
        return(cbind(mz = z[1] + isotopePopulationMz,
                     mzmin = z[2] + isotopePopulationMz - mz_ext,
                     mzmax = z[3] + isotopePopulationMz + mz_ext,
                     scmin = z[4],
                     scmax = z[5],
                     length = -1,
                     intensity = -1,
                     scale = z[6])
               )
    })
    return(do.call(rbind, newRois))
}

#' @param peaks. see do_define_isotopes
#' 
#' @param polarity character(1) defining the polarity, either \code{"positive"}
#'     or \code{"negative"}.
#' 
#' @return see do_define_isotopes.
#'
#' @noRd
do_define_adducts <- function(peaks., polarity = "positive") {
    req_cols <- c("mz", "mzmin", "mzmax", "scmin", "scmax", "scale")
    if (is.null(dim(peaks.)))
        stop("'peaks.' has to be a matrix or data.frame!")
    if (!all(req_cols %in% colnames(peaks.))) {
        not_there <- req_cols[!(req_cols %in% colnames(peaks.))]
        stop("'peaks.' lacks required columns ",
             paste0("'", not_there, "'", collapse = ","), "!")
    }
    if (is.data.frame(peaks.))
        peaks. <- as.matrix(peaks.)
    ## considered adduct distances
    ## reference: Huang N.; Siegel M.M.1; Kruppa G.H.; Laukien F.H.; J Am Soc
    ## Mass Spectrom 1999, 10, 1166–1173; Automation of a Fourier transform ion
    ## cyclotron resonance mass spectrometer for acquisition, analysis, and
    ## e-mailing of high-resolution exact-mass electrospray ionization mass
    ## spectral data
    ## see also for contaminants: Interferences and contaminants encountered
    ## in modern mass spectrometry (Bernd O. Keller, Jie Sui, Alex B. Young
    ## and Randy M. Whittal, ANALYTICA CHIMICA ACTA, 627 (1): 71-81)

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
    mDMSO    <- mC * 2 + mH * 6 + mS + mO        # dimethylsulfoxid
    mACN     <- mC * 2 + mH * 3 + mN             # acetonitril
    mIsoProp <- mC * 3 + mH * 8 + mO             # isopropanol
    mNH4     <- mN + mH * 4                      # ammonium
    mCH3OH   <- mC + mH * 3 + mO + mH            # methanol
    mH2O     <- mH * 2 + mO                      # water
    mFA      <- mC + mH * 2 + mO * 2             # formic acid
    mHAc     <- mC + mH * 3 + mC + mO + mO + mH  # acetic acid
    mTFA     <- mC + mF * 3 + mC + mO + mO + mH  # trifluoroacetic acid

    switch(polarity,
           "positive"={
               adductPopulationMz <- unlist(c(
                   ## [M+H]+ to [M+H]+  (Reference)
                   ## [M+H]+ to [M+NH4]+
                   function(mass){ mass - mH + mNH4 },
                   ## [M+H]+ to [M+Na]+
                   function(mass){ mass - mH + mNa },
                   ## [M+H]+ to [M+CH3OH+H]+
                   function(mass){ mass + mCH3OH },
                   ## [M+H]+ to [M+K]+
                   function(mass){ mass - mH + mK },
                   ## [M+H]+ to [M+ACN+H]+
                   function(mass){ mass + mACN },
                   ## [M+H]+ to [M+2Na-H]+
                   function(mass){ mass - 2 * mH + 2 * mNa },
                   ## [M+H]+ to [M+IsoProp+H]+
                   function(mass){ mass + mIsoProp },
                   ## [M+H]+ to [M+ACN+Na]+
                   function(mass){ mass - mH + mACN + mNa },
                   ## [M+H]+ to [M+2K-H]+
                   function(mass){ mass - 2 * mH + 2 * mK },
                   ## [M+H]+ to [M+DMSO+H]+
                   function(mass){ mass + mDMSO },
                   ## [M+H]+ to [M+2*ACN+H]+
                   function(mass){ mass + 2 * mACN },
                   ## [M+H]+ to [M+IsoProp+Na+H]+ TODO double-charged?
                   function(mass){ mass + mIsoProp + mNa },
                   ## [M+H]+ to [2M+H]+
                   function(mass){ (mass - mH) * 2 + mH },
                   ## [M+H]+ to [2M+NH4]+
                   function(mass){ (mass - mH) * 2 + mNH4 },
                   ## [M+H]+ to [2M+Na]+
                   function(mass){ (mass - mH) * 2 + mNa },
                   ## [M+H]+ to [2M+K]+
                   function(mass){ (mass - mH) * 2 + mK },
                   ## [M+H]+ to [2M+ACN+H]+
                   function(mass){ (mass - mH) * 2 + mACN + mH },
                   ## [M+H]+ to [2M+ACN+Na]+
                   function(mass){ (mass - mH) * 2 + mACN + mNa },
                   ## [M+H]+ to [2M+3*H2O+2*H]2+
                   function(mass){((mass - mH) * 2 + 3 * mH2O + 2 * mH) / 2 },
                   ## [M+H]+ to [M+2*H]2+
                   function(mass){ (mass + mH) / 2 },
                   ## [M+H]+ to [M+H+NH4]2+
                   function(mass){ (mass + mNH4) / 2 },
                   ## [M+H]+ to [M+H+Na]2+
                   function(mass){ (mass + mNa) / 2 },
                   ## [M+H]+ to [M+H+K]2+
                   function(mass){ (mass + mK) / 2 },
                   ## [M+H]+ to [M+ACN+2*H]2+
                   function(mass){ (mass + mACN + mH) / 2 },
                   ## [M+H]+ to [M+2*Na]2+
                   function(mass){ (mass - mH + 2 * mNa) / 2 },
                   ## [M+H]+ to [M+2*ACN+2*H]2+
                   function(mass){ (mass + 2 * mACN + mH) / 2 },
                   ## [M+H]+ to [M+3*ACN+2*H]2+
                   function(mass){ (mass + 3 * mACN + mH) / 2 },
                   ## [M+H]+ to [M+3*H]3+
                   function(mass){ (mass + 2 * mH) / 3 },
                   ## [M+H]+ to [M+2*H+Na]3+
                   function(mass){ (mass + mH + mNa) / 3 },
                   ## [M+H]+ to [M+H+2*Na]3+
                   function(mass){ (mass + 2 * mNa) / 3 },
                   ## [M+H]+ to [M+3*Na]3+
                   function(mass){ (mass - mH + 3 * mNa) / 3 }
               ))
           },
           "negative" = {
               adductPopulationMz <- unlist(c(
                   ## [M-H]+ to [M-H]+  (Reference)
                   ## [M-H]+ to [M-H2O-H]+
                   function(mass){ mass - mH2O },
                   ## [M-H]+ to [M+Na-2*H]+
                   function(mass){ mass - mH + mNa },
                   ## [M-H]+ to [M+Cl]+
                   function(mass){ mass + mH + mCl },
                   ## [M-H]+ to [M+K-2*H]+
                   function(mass){ mass - mH + mK },
                   ## [M-H]+ to [M+FA-H]+
                   function(mass){ mass + mFA },
                   ## [M-H]+ to [M+HAc-H]+
                   function(mass){ mass + mHAc },
                   ## [M-H]+ to [M+Br]+
                   function(mass){ mass + mH + mBr },
                   ## [M-H]+ to [M+TFA-H]+
                   function(mass){ mass + mTFA },
                   ## [M-H]+ to [2M-H]+
                   function(mass){ (mass + mH) * 2 - mH },
                   ## [M-H]+ to [2M+FA-H]+
                   function(mass){ (mass + mH) * 2 + mFA - mH },
                   ## [M-H]+ to [2M+HAc-H]+
                   function(mass){ (mass + mH) * 2 + mHAc - mH },
                   ## [M-H]+ to [3M-H]+
                   function(mass){ (mass + mH) * 3 - mH },
                   ## [M-H]+ to [M-2*H]2+
                   function(mass){ (mass - mH) / 2 },
                   ## [M-H]+ to [M-3*H]3+
                   function(mass){ (mass - 2 * mH) / 3 }
               ))
           },
           "unknown"={
               warning("Unknown polarity! No adduct ROIs have been added.")
           },
           stop("Unknown polarity (", polarity, ")!")
           )

    req_cols <- c("mz", "mzmin", "mzmax", "scmin", "scmax", "scale")
    roiL <- split(peaks.[, req_cols, drop = FALSE], f = 1:nrow(peaks.))

    newRois <- lapply(roiL, function(z) {
        mzDiff <- unlist(lapply(adductPopulationMz, function(x) {
            do.call(x, list(mass = z[1]))
        }))
        return(cbind(mz = z[1] + mzDiff,
                     mzmin = z[2] + mzDiff,
                     mzmax = z[3] + mzDiff,
                     scmin = z[4],
                     scmax = z[5],
                     length = -1,
                     intensity = -1,
                     scale = z[6]))
    })

    newRois <- do.call(rbind, newRois)
    ## Remove ROIs with negative or zero mzmin.
    return(newRois[newRois[, "mzmin"] > 0, , drop = FALSE])
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
             " have to match. Also, 'length(mz)' should be equal to",
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
    tmp <- capture.output(
        res <- .Call("massifquant", mz, int, scanindex, scantime,
                     as.double(mzrange), as.integer(scanrange),
                     as.integer(length(scantime)), as.double(minIntensity),
                     as.integer(minCentroids), as.double(consecMissedLim),
                     as.double(ppm), as.double(criticalVal), as.integer(segs),
                     as.integer(scanBack), PACKAGE ='xcms' )
    )
    res
}

############################################################
## do_findChromPeaks_centWaveWithPredIsoROIs
## 1) Run a centWave.
## 2) Predict isotope ROIs for the identified peaks.
## 3) centWave on the predicted isotope ROIs.
## 4) combine both lists of identified peaks removing overlapping ones by
##    keeping the peak with the largest signal intensity.
#' @title Core API function for two-step centWave peak detection with isotopes
#'
#' @description The \code{do_findChromPeaks_centWaveWithPredIsoROIs} performs a
#'     two-step centWave based peak detection: chromatographic peaks are
#'     identified using centWave followed by a prediction of the location of
#'     the identified peaks' isotopes in the mz-retention time space. These
#'     locations are fed as \emph{regions of interest} (ROIs) to a subsequent
#'     centWave run. All non overlapping peaks from these two peak detection
#'     runs are reported as the final list of identified peaks.
#'
#' @details For more details on the centWave algorithm see
#'     \code{\link{centWave}}.
#'
#' @inheritParams findChromPeaks-centWave
#' 
#' @inheritParams findChromPeaks-centWaveWithPredIsoROIs
#' 
#' @inheritParams do_findChromPeaks_centWave
#'
#' @family core peak detection functions
#' 
#' @return A matrix, each row representing an identified chromatographic peak.
#'     All non-overlapping peaks identified in both centWave runs are reported.
#'     The matrix columns are:
#'     \describe{
#'     \item{mz}{Intensity weighted mean of m/z values of the peaks across scans.}
#'     \item{mzmin}{Minimum m/z of the peaks.}
#'     \item{mzmax}{Maximum m/z of the peaks.}
#'     \item{rt}{Retention time of the peak's midpoint.}
#'     \item{rtmin}{Minimum retention time of the peak.}
#'     \item{rtmax}{Maximum retention time of the peak.}
#'     \item{into}{Integrated (original) intensity of the peak.}
#'     \item{intb}{Per-peak baseline corrected integrated peak intensity.}
#'     \item{maxo}{Maximum intensity of the peak.}
#'     \item{sn}{Signal to noise ratio, defined as \code{(maxo - baseline)/sd},
#'     \code{sd} being the standard deviation of local chromatographic noise.}
#'     \item{egauss}{RMSE of Gaussian fit.}
#'     }
#'     Additional columns for \code{verboseColumns = TRUE}:
#'     \describe{
#'     \item{mu}{Gaussian parameter mu.}
#'     \item{sigma}{Gaussian parameter sigma.}
#'     \item{h}{Gaussian parameter h.}
#'     \item{f}{Region number of the m/z ROI where the peak was localized.}
#'     \item{dppm}{m/z deviation of mass trace across scans in ppm.}
#'     \item{scale}{Scale on which the peak was localized.}
#'     \item{scpos}{Peak position found by wavelet analysis (scan number).}
#'     \item{scmin}{Left peak limit found by wavelet analysis (scan number).}
#'     \item{scmax}{Right peak limit found by wavelet analysis (scan numer).}
#'     }
#' 
#' @rdname do_findChromPeaks_centWaveWithPredIsoROIs
#'
#' @author Hendrik Treutler, Johannes Rainer
do_findChromPeaks_centWaveWithPredIsoROIs <-
    function(mz, int, scantime, valsPerSpect, ppm = 25, peakwidth = c(20, 50),
             snthresh = 10, prefilter = c(3, 100), mzCenterFun = "wMean",
             integrate = 1, mzdiff = -0.001, fitgauss = FALSE, noise = 0,
             verboseColumns = FALSE, roiList = list(),
             firstBaselineCheck = TRUE, roiScales = NULL, snthreshIsoROIs = 6.25,
             maxCharge = 3, maxIso = 5, mzIntervalExtension = TRUE,
             polarity = "unknown") {
        ## Input argument checking: most of it will be done in
        ## do_findChromPeaks_centWave
        polarity <- match.arg(polarity, c("positive", "negative", "unknown"))

        ## 1) First centWave
        feats_1 <- do_findChromPeaks_centWave(mz = mz, int = int,
                                              scantime = scantime,
                                              valsPerSpect = valsPerSpect,
                                              ppm = ppm,
                                              peakwidth = peakwidth,
                                              snthresh = snthresh,
                                              prefilter = prefilter,
                                              mzCenterFun = mzCenterFun,
                                              integrate = integrate,
                                              mzdiff = mzdiff, fitgauss = fitgauss,
                                              noise = noise,
                                              verboseColumns = TRUE,
                                              roiList = roiList,
                                              firstBaselineCheck = firstBaselineCheck,
                                              roiScales = roiScales)
        return(do_findChromPeaks_addPredIsoROIs(mz = mz, int = int,
                                                scantime = scantime,
                                                valsPerSpect = valsPerSpect,
                                                ppm = ppm,
                                                peakwidth = peakwidth,
                                                snthresh = snthreshIsoROIs,
                                                prefilter = prefilter,
                                                mzCenterFun = mzCenterFun,
                                                integrate = integrate,
                                                mzdiff = mzdiff,
                                                fitgauss = fitgauss,
                                                noise = noise,
                                                verboseColumns = verboseColumns,
                                                peaks. = feats_1,
                                                maxCharge = maxCharge,
                                                maxIso = maxIso,
                                                mzIntervalExtension = mzIntervalExtension,
                                                polarity = polarity))
    }
#' @description The \code{do_findChromPeaks_centWaveAddPredIsoROIs} performs
#'     centWave based peak detection based in regions of interest (ROIs)
#'     representing predicted isotopes for the peaks submitted with argument
#'     \code{peaks.}. The function returns a matrix with the identified peaks
#'     consisting of all input peaks and peaks representing predicted isotopes
#'     of these (if found by the centWave algorithm).
#'
#' @param peaks. A matrix or \code{xcmsPeaks} object such as one returned by
#'     a call to \code{link{do_findChromPeaks_centWave}} or
#'     \code{link{findPeaks.centWave}} (both with \code{verboseColumns = TRUE})
#'     with the peaks for which isotopes should be predicted and used for an
#'     additional peak detectoin using the centWave method. Required columns
#'     are: \code{"mz"}, \code{"mzmin"}, \code{"mzmax"}, \code{"scmin"},
#'     \code{"scmax"}, \code{"scale"} and \code{"into"}.
#'
#' @param snthresh For \code{do_findChromPeaks_addPredIsoROIs}:
#'     numeric(1) defining the signal to noise threshold for the centWave
#'     algorithm. For \code{do_findChromPeaks_centWaveWithPredIsoROIs}:
#'     numeric(1) defining the signal to noise threshold for the initial
#'     (first) centWave run.
#'
#' @inheritParams findChromPeaks-centWave
#' 
#' @inheritParams do_findChromPeaks_centWave
#'
#' @rdname do_findChromPeaks_centWaveWithPredIsoROIs
do_findChromPeaks_addPredIsoROIs <-
    function(mz, int, scantime, valsPerSpect, ppm = 25, peakwidth = c(20, 50),
             snthresh = 6.25, prefilter = c(3, 100), mzCenterFun = "wMean",
             integrate = 1, mzdiff = -0.001, fitgauss = FALSE, noise = 0,
             verboseColumns = FALSE, peaks. = NULL,
             maxCharge = 3, maxIso = 5, mzIntervalExtension = TRUE,
             polarity = "unknown") {
        ## Input argument checking: most of it will be done in
        ## do_findChromPeaks_centWave
        polarity <- match.arg(polarity, c("positive", "negative", "unknown"))

        ## These variables might at some point be added as function args.
        addNewIsotopeROIs <- TRUE
        addNewAdductROIs <- FALSE
        ## 2) predict isotope and/or adduct ROIs
        f_mod <- peaks.
        ## Extend the mzmin and mzmax if needed.
        tittle <- peaks.[, "mz"] * (ppm / 2) / 1E6
        expand_mz <- (peaks.[, "mzmax"] - peaks.[, "mzmin"]) < (tittle * 2)
        if (any(expand_mz)) {
            f_mod[expand_mz, "mzmin"] <- peaks.[expand_mz, "mz"] -
                tittle[expand_mz]
            f_mod[expand_mz, "mzmax"] <- peaks.[expand_mz, "mz"] + tittle[expand_mz]
        }
        ## Add predicted ROIs
        if (addNewIsotopeROIs) {
            iso_ROIs <- do_define_isotopes(peaks. = f_mod,
                                           maxCharge = maxCharge,
                                           maxIso = maxIso,
                                           mzIntervalExtension = mzIntervalExtension)
        } else {
            iso_ROIs <- matrix(nrow = 0, ncol = 8)
            colnames(iso_ROIs) <- c("mz", "mzmin", "mzmax", "scmin", "scmax",
                                    "length", "intensity", "scale")
        }
        if (addNewAdductROIs) {
            add_ROIs <- do_define_adducts(peaks. = f_mod, polarity = polarity)
        } else {
            add_ROIs <- matrix(nrow = 0, ncol = 8)
            colnames(iso_ROIs) <- c("mz", "mzmin", "mzmax", "scmin", "scmax",
                                    "length", "intensity", "scale")
        }
        newROIs <- rbind(iso_ROIs, add_ROIs)
        rm(f_mod)
        if (nrow(newROIs) == 0)
            return(peaks.)
        ## Remove ROIs that are out of mz range:
        mz_range <- range(mz)
        newROIs <- newROIs[newROIs[, "mzmin"] >= mz_range[1] &
                           newROIs[, "mzmax"] <= mz_range[2], , drop = FALSE]
        ## Remove ROIs with too low signal:
        keep_me <- logical(nrow(newROIs))
        scanindex <- as.integer(valueCount2ScanIndex(valsPerSpect))
        for (i in 1:nrow(newROIs)) {
            vals <- .Call("getEIC", mz, int, scanindex,
                          as.double(newROIs[i, c("mzmin", "mzmax")]),
                          as.integer(newROIs[i, c("scmin", "scmax")]),
                          as.integer(length(scantime)), PACKAGE ='xcms' )
            keep_me[i] <- sum(vals$intensity, na.rm = TRUE) >= 10
        }
        newROIs <- newROIs[keep_me, , drop = FALSE]

        if (nrow(newROIs) == 0) {
            warning("No isotope or adduct ROIs for the identified peaks with a ",
                    "valid signal found!")
            return(peaks.)
        }
        
        ## 3) centWave using the identified ROIs.
        roiL <- split(as.data.frame(newROIs), f = 1:nrow(newROIs))
        feats_2 <- do_findChromPeaks_centWave(mz = mz, int = int,
                                              scantime = scantime,
                                              valsPerSpect = valsPerSpect,
                                              ppm = ppm, peakwidth = peakwidth,
                                              snthresh = snthresh,
                                              prefilter = prefilter,
                                              mzCenterFun = mzCenterFun,
                                              integrate = integrate,
                                              mzdiff = mzdiff, fitgauss = fitgauss,
                                              noise = noise,
                                              verboseColumns = verboseColumns,
                                              roiList = roiL,
                                              firstBaselineCheck = FALSE,
                                              roiScales = newROIs[, "scale"])
        ## Clean up of the results:
        if (nrow(feats_2) > 0) {
            ## remove NaNs
            any_na <- is.na(rowSums(feats_2[, c("mz", "mzmin", "mzmax", "rt",
                                                "rtmin", "rtmax")]))
            if (any(any_na))
                feats_2 <- feats_2[!any_na, , drop = FALSE]
            no_mz_width <- (feats_2[, "mzmax"] - feats_2[, "mzmin"]) == 0
            no_rt_width <- (feats_2[, "rtmax"] - feats_2[, "rtmin"]) == 0
            ## remove empty area
            ## no_area <- (feats_2[, "mzmax"] - feats_2[, "mzmin"]) == 0 ||
            ##     (feats_2[, "rtmax"] - feats_2[, "rtmin"]) == 0
            no_area <- no_mz_width || no_rt_width
            if (any(no_area))
                feats_2 <- feats_2[!no_area, , drop = FALSE]
        }
        ## 4) Check and remove ROIs overlapping with peaks.
        if (nrow(feats_2) > 0) {
            ## Comparing each ROI with each peak; slightly modified from the original
            ## code in which we prevent calling apply followed by two lapply.
            removeROIs <- rep(FALSE, nrow(feats_2))
            removeFeats <- rep(FALSE, nrow(peaks.))
            overlapProportionThreshold <- 0.01
            for (i in 1:nrow(feats_2)) {
                ## Compare ROI i with all peaks (peaks) and check if its
                ## overlapping
                ## mz
                roiMzCenter <- (feats_2[i, "mzmin"] + feats_2[i, "mzmax"]) / 2
                peakMzCenter <- (peaks.[, "mzmin"] + peaks.[, "mzmax"]) / 2
                roiMzRadius <- (feats_2[i, "mzmax"] - feats_2[i, "mzmin"]) / 2
                peakMzRadius <- (peaks.[, "mzmax"] - peaks.[, "mzmin"]) / 2
                overlappingMz <- abs(peakMzCenter - roiMzCenter) <=
                    (roiMzRadius + peakMzRadius)
                ## rt
                roiRtCenter <- (feats_2[i, "rtmin"] + feats_2[i, "rtmax"]) / 2
                peakRtCenter <- (peaks.[, "rtmin"] + peaks.[, "rtmax"]) / 2
                roiRtRadius <- (feats_2[i, "rtmax"] - feats_2[i, "rtmin"]) / 2
                peakRtRadius <- (peaks.[, "rtmax"] - peaks.[, "rtmin"]) / 2
                overlappingRt <- abs(peakRtCenter - roiRtCenter) <=
                    (roiRtRadius + peakRtRadius)
                is_overlapping <- overlappingMz & overlappingRt
                ## Now determine whether we remove the ROI or the peak, depending
                ## on the raw signal intensity.
                if (any(is_overlapping)) {
                    if (any(peaks.[is_overlapping, "into"] > feats_2[i, "into"])) {
                        removeROIs[i] <- TRUE
                    } else {
                        removeFeats[is_overlapping] <- TRUE
                    }
                }
            }
            feats_2 <- feats_2[!removeROIs, , drop = FALSE]
            peaks. <- peaks.[!removeFeats, , drop = FALSE]
        }
        if (!verboseColumns)
            peaks. <- peaks.[ , c("mz", "mzmin", "mzmax", "rt", "rtmin",
                                  "rtmax", "into", "intb", "maxo", "sn")]
        if (nrow(feats_2) == 0)
            return(peaks.)
        else
            return(rbind(peaks., feats_2))
    }

do_findChromPeaks_addPredIsoROIs_mod <-
    function(mz, int, scantime, valsPerSpect, ppm = 25, peakwidth = c(20, 50),
             snthresh = 6.25, prefilter = c(3, 100), mzCenterFun = "wMean",
             integrate = 1, mzdiff = -0.001, fitgauss = FALSE, noise = 0,
             verboseColumns = FALSE, peaks. = NULL,
             maxCharge = 3, maxIso = 5, mzIntervalExtension = TRUE,
             polarity = "unknown") {
        ## Input argument checking: most of it will be done in
        ## do_findChromPeaks_centWave
        polarity <- match.arg(polarity, c("positive", "negative", "unknown"))

        ## These variables might at some point be added as function args.
        addNewIsotopeROIs <- TRUE
        addNewAdductROIs <- FALSE
        ## 2) predict isotope and/or adduct ROIs
        f_mod <- peaks.
        ## Extend the mzmin and mzmax if needed.
        tittle <- peaks.[, "mz"] * (ppm / 2) / 1E6
        expand_mz <- (peaks.[, "mzmax"] - peaks.[, "mzmin"]) < (tittle * 2)
        if (any(expand_mz)) {
            f_mod[expand_mz, "mzmin"] <- peaks.[expand_mz, "mz"] -
                tittle[expand_mz]
            f_mod[expand_mz, "mzmax"] <- peaks.[expand_mz, "mz"] + tittle[expand_mz]
        }
        ## Add predicted ROIs
        if (addNewIsotopeROIs) {
            iso_ROIs <- do_define_isotopes(peaks. = f_mod,
                                           maxCharge = maxCharge,
                                           maxIso = maxIso,
                                           mzIntervalExtension = mzIntervalExtension)
        } else {
            iso_ROIs <- matrix(nrow = 0, ncol = 8)
            colnames(iso_ROIs) <- c("mz", "mzmin", "mzmax", "scmin", "scmax",
                                    "length", "intensity", "scale")
        }
        if (addNewAdductROIs) {
            add_ROIs <- do_define_adducts(peaks. = f_mod, polarity = polarity)
        } else {
            add_ROIs <- matrix(nrow = 0, ncol = 8)
            colnames(iso_ROIs) <- c("mz", "mzmin", "mzmax", "scmin", "scmax",
                                    "length", "intensity", "scale")
        }
        newROIs <- rbind(iso_ROIs, add_ROIs)
        rm(f_mod)
        if (nrow(newROIs) == 0)
            return(peaks.)
        ## Remove ROIs that are out of mz range:
        mz_range <- range(mz)
        newROIs <- newROIs[newROIs[, "mzmin"] >= mz_range[1] &
                           newROIs[, "mzmax"] <= mz_range[2], , drop = FALSE]
        ## Remove ROIs with too low signal:
        keep_me <- logical(nrow(newROIs))
        scanindex <- as.integer(valueCount2ScanIndex(valsPerSpect))
        for (i in 1:nrow(newROIs)) {
            vals <- .Call("getEIC", mz, int, scanindex,
                          as.double(newROIs[i, c("mzmin", "mzmax")]),
                          as.integer(newROIs[i, c("scmin", "scmax")]),
                          as.integer(length(scantime)), PACKAGE ='xcms' )
            keep_me[i] <- sum(vals$intensity, na.rm = TRUE) >= 10
        }
        newROIs <- newROIs[keep_me, , drop = FALSE]

        if (nrow(newROIs) == 0) {
            warning("No isotope or adduct ROIs for the identified peaks with a ",
                    "valid signal found!")
            return(peaks.)
        }
        cat("No. of input peaks: ", nrow(peaks.), "\n")
        
        ## 3) centWave using the identified ROIs.
        roiL <- split(as.data.frame(newROIs), f = 1:nrow(newROIs))
        cat("Identified iso ROIs: ", length(roiL), "\n")
        feats_2 <- do_findChromPeaks_centWave(mz = mz, int = int,
                                              scantime = scantime,
                                              valsPerSpect = valsPerSpect,
                                              ppm = ppm, peakwidth = peakwidth,
                                              snthresh = snthresh,
                                              prefilter = prefilter,
                                              mzCenterFun = mzCenterFun,
                                              integrate = integrate,
                                              mzdiff = mzdiff, fitgauss = fitgauss,
                                              noise = noise,
                                              verboseColumns = verboseColumns,
                                              roiList = roiL,
                                              firstBaselineCheck = FALSE,
                                              roiScales = newROIs[, "scale"])
        cat("No. of chrom. peaks found in ROIs: ", nrow(feats_2), "\n")
        ## Clean up of the results:
        if (nrow(feats_2) > 0) {
            ## remove NaNs
            any_na <- is.na(rowSums(feats_2[, c("mz", "mzmin", "mzmax", "rt",
                                                "rtmin", "rtmax")]))
            if (any(any_na))
                feats_2 <- feats_2[!any_na, , drop = FALSE]
            no_mz_width <- (feats_2[, "mzmax"] - feats_2[, "mzmin"]) == 0
            no_rt_width <- (feats_2[, "rtmax"] - feats_2[, "rtmin"]) == 0

            cat("No. of peaks with NA values: ", sum(any_na), "\n")
            cat("No. of peaks without mz width: ", sum(no_mz_width), "\n")
            cat("No. of peaks without rt width: ", sum(no_rt_width), "\n")
            ## remove empty area
            ## no_area <- (feats_2[, "mzmax"] - feats_2[, "mzmin"]) == 0 ||
            ##     (feats_2[, "rtmax"] - feats_2[, "rtmin"]) == 0
            ## no_area <- no_mz_width || no_rt_width
            no_area <- no_mz_width
            if (any(no_area))
                feats_2 <- feats_2[!no_area, , drop = FALSE]
        }
        cat("After removing NAs or empty are peaks: ", nrow(feats_2), "\n")
        ## 4) Check and remove ROIs overlapping with peaks.
        if (nrow(feats_2) > 0) {
            ## Comparing each ROI with each peak; slightly modified from the original
            ## code in which we prevent calling apply followed by two lapply.
            removeROIs <- rep(FALSE, nrow(feats_2))
            removeFeats <- rep(FALSE, nrow(peaks.))
            overlapProportionThreshold <- 0.01
            for (i in 1:nrow(feats_2)) {
                ## Compare ROI i with all peaks (peaks) and check if its
                ## overlapping
                ## mz
                roiMzCenter <- (feats_2[i, "mzmin"] + feats_2[i, "mzmax"]) / 2
                peakMzCenter <- (peaks.[, "mzmin"] + peaks.[, "mzmax"]) / 2
                roiMzRadius <- (feats_2[i, "mzmax"] - feats_2[i, "mzmin"]) / 2
                peakMzRadius <- (peaks.[, "mzmax"] - peaks.[, "mzmin"]) / 2
                overlappingMz <- abs(peakMzCenter - roiMzCenter) <=
                    (roiMzRadius + peakMzRadius)
                ## rt
                roiRtCenter <- (feats_2[i, "rtmin"] + feats_2[i, "rtmax"]) / 2
                peakRtCenter <- (peaks.[, "rtmin"] + peaks.[, "rtmax"]) / 2
                roiRtRadius <- (feats_2[i, "rtmax"] - feats_2[i, "rtmin"]) / 2
                peakRtRadius <- (peaks.[, "rtmax"] - peaks.[, "rtmin"]) / 2
                overlappingRt <- abs(peakRtCenter - roiRtCenter) <=
                    (roiRtRadius + peakRtRadius)
                is_overlapping <- overlappingMz & overlappingRt
                ## Now determine whether we remove the ROI or the peak, depending
                ## on the raw signal intensity.
                if (any(is_overlapping)) {
                    if (any(peaks.[is_overlapping, "into"] > feats_2[i, "into"])) {
                        removeROIs[i] <- TRUE
                    } else {
                        removeFeats[is_overlapping] <- TRUE
                    }
                }
            }
            feats_2 <- feats_2[!removeROIs, , drop = FALSE]
            peaks. <- peaks.[!removeFeats, , drop = FALSE]
        }
        cat("After removing overlapping peaks: ", nrow(feats_2), "\n")
        cat("After removing overlapping peaks (peaks.): ", nrow(peaks.), "\n")
        if (!verboseColumns)
            peaks. <- peaks.[ , c("mz", "mzmin", "mzmax", "rt", "rtmin",
                                  "rtmax", "into", "intb", "maxo", "sn")]
        ## For now just return the new ones.
        return(feats_2)
        if (nrow(feats_2) == 0)
            return(peaks.)
        else
            return(rbind(peaks., feats_2))
    }

#' @title Identify peaks in chromatographic data using matchedFilter
#'
#' @description
#'
#' The function performs peak detection using the [matchedFilter] algorithm
#' on chromatographic data (i.e. with only intensities and retention time).
#'
#' @param int `numeric` with intensity values.
#'
#' @param rt `numeric` with the retention time for the intensities. Length has
#'     to be equal to `length(int)`.
#'
#' @param fwhm `numeric(1)` specifying the full width at half maximum
#'     of matched filtration gaussian model peak. Only used to calculate the
#'     actual sigma, see below.
#'
#' @param sigma `numeric(1)` specifying the standard deviation (width)
#'     of the matched filtration model peak.
#'
#' @param max `numeric(1)` with the maximal number of peaks that are expected/
#'     will bbe detected in the data
#'
#' @param snthresh `numeric(1)` defining the signal to noise cut-off to be used
#'     in the peak detection step.
#'
#' @param ... currently ignored.
#' 
#' @family peak detection functions for chromatographic data
#'
#' @seealso [matchedFilter] for a detailed description of the peak detection
#'     method.
#'
#' @author Johannes Rainer
#'
#' @return
#' 
#' A matrix, each row representing an identified chromatographic peak, with
#' columns:
#' 
#' - `"rt"`: retention time of the peak's midpoint (time of the maximum signal).
#' - `"rtmin"`: minimum retention time of the peak.
#' - `"rtmax"`: maximum retention time of the peak.
#' - `"into"`: integrated (original) intensity of the peak.
#' - `"intf"`: integrated intensity of the filtered peak.
#' - `"maxo"`: maximum (original) intensity of the peak.
#' - `"maxf"`" maximum intensity of the filtered peak.
#' - `"sn"`: signal to noise ratio of the peak.
#'
#' @md
#'
#' @examples
#'
#' ## Read one file from the faahKO package
#' od <- readMSData(system.file("cdf/KO/ko15.CDF", package = "faahKO"),
#'     mode = "onDisk")
#'
#' ## Extract chromatographic data for a small m/z range
#' chr <- chromatogram(od, mz = c(272.1, 272.3))[1, 1]
#'
#' pks <- peaksWithMatchedFilter(intensity(chr), rtime(chr))
#' pks
#'
#' ## Plotting the data
#' plot(rtime(chr), intensity(chr), type = "h")
#' rect(xleft = pks[, "rtmin"], xright = pks[, "rtmax"], ybottom = c(0, 0),
#'     ytop = pks[, "maxo"], border = "red")
peaksWithMatchedFilter <- function(int, rt, fwhm = 30, sigma = fwhm / 2.3548,
                                   max = 20, snthresh = 10, ...) {
    if (missing(int) | missing(rt))
        stop("Arguments 'int' and 'rt' are required")
    if (length(int) != length(rt))
        stop("'int' and 'rt' must have the same length")

    ## Replace NAs with 0 - that's how the original code handled it.
    nas <- is.na(int)
    int[nas] <- 0
    
    n_vals <- length(int)
    N <- nextn(n_vals)
    rtrange <- range(rt)
    x <- c(0:(N / 2), -(ceiling(N / 2 - 1)):-1) *
        (rtrange[2] - rtrange[1]) / (n_vals - 1)
    filt <- -attr(eval(deriv3(~ 1 / (sigma * sqrt(2 * pi)) *
                                  exp(-x^2 / (2 * sigma^2)), "x")),
                  "hessian")
    filt <- filt / sqrt(sum(filt^2))
    filt <- fft(filt, inverse = TRUE) / length(filt)

    cnames <- c("rt", "rtmin", "rtmax", "into", "intf", "maxo", "maxf", "sn")
    num <- 0
    ResList <- list()
    int_filt <- filtfft(int, filt)
    glob_max <- max(int_filt)
    for (j in seq(length = max)) {
        max_idx <- which.max(int_filt)
        noise <- mean(int[int > 0])
        ##noise <- mean(yfilt[yfilt >= 0])
        sn <- int_filt[max_idx] / noise
        if (int_filt[max_idx] > 0 && int_filt[max_idx] > snthresh * noise &&
            int[max_idx] > 0) {
            peak_range <- descendZero(int_filt, max_idx)
            peak_range_idx <- peak_range[1]:peak_range[2]
            ## Define into and intf, i.e. the integrated signal.
            peak_width <- diff(rt[peak_range]) / diff(peak_range)
            into <- peak_width * sum(int[peak_range_idx])
            intf <- peak_width * sum(int_filt[peak_range_idx])
            num <- num + 1
            ResList[[num]] <- c(rt[max_idx], rt[peak_range], into, intf,
                                max(int[peak_range_idx]), int_filt[max_idx], sn)
            ## "remove" current peak from the data
            int_filt[peak_range_idx] <- 0
        } else
            break
    }
    if (length(ResList)) {
        res <- do.call(rbind, ResList)
        colnames(res) <- cnames
    } else {
        warning("No peaks found with current settings")
        res <- matrix(nrow = 0, ncol = length(cnames),
                      dimnames = list(character(), cnames))
    }
    res
}

#' @title Identify peaks in chromatographic data using centWave
#'
#' @description
#' 
#' `peaksWithCentWave` identifies (chromatographic) peaks in purely
#' chromatographic data, i.e. based on intensity and retention time values
#' without m/z values.
#'
#' @details
#'
#' The method uses the same algorithm for the peak detection than [centWave],
#' employs however a different approach to identify the initial regions in
#' which the peak detection is performed (i.e. the *regions of interest* ROI).
#' The method first identifies all local maxima in the chromatographic data and
#' defines the corresponding positions +/- `peakwidth[2]` as the ROIs. Noise
#' estimation bases also on these ROIs and can thus be different from [centWave]
#' resulting in different signal to noise ratios.
#' 
#' @param int `numeric` with intensity values.
#'
#' @param rt `numeric` with the retention time for the intensities. Length has
#'     to be equal to `length(int)`.
#'
#' @param peakwidth `numeric(2)` with the lower and upper bound of the
#'     expected peak width.
#'
#' @param snthresh `numeric(1)` defining the signal to noise ratio cutoff.
#'     Peaks with a signal to noise ratio < `snthresh` are omitted.
#' 
#' @param prefilter `numeric(2)` (`c(k, I)`): only regions of interest with at
#'     least `k` centroids with signal `>= I` are returned in the first
#'     step.
#'
#' @param integrate `numeric(1)`, integration method. For `integrate = 1` peak
#'     limits are found through descending on the mexican hat filtered data,
#'     for `integrate = 2` the descend is done on the real data. The latter
#'     method is more accurate but prone to noise, while the former is more
#'     robust, but less exact.
#'
#' @param fitgauss `logical(1)` whether or not a Gaussian should be fitted
#'     to each peak.
#' 
#' @param noise `numeric(1)` defining the minimum required intensity for
#'     centroids to be considered in the first analysis step (definition of
#'     the *regions of interest*).
#'
#' @param verboseColumns `logical(1)`: whether additional peak meta data
#'     columns should be returned.
#'
#' @param firstBaselineCheck `logical(1)`. If `TRUE` continuous data within
#'     regions of interest is checked to be above the first baseline.
#'
#' @param ... currently ignored.
#' 
#' @family peak detection functions for chromatographic data
#'
#' @seealso [centWave] for a detailed description of the peak detection
#'     method.
#'
#' @author Johannes Rainer
#'
#' @return
#' 
#' A matrix, each row representing an identified chromatographic peak, with
#' columns:
#' 
#' - `"rt"`: retention time of the peak's midpoint (time of the maximum signal).
#' - `"rtmin"`: minimum retention time of the peak.
#' - `"rtmax"`: maximum retention time of the peak.
#' - `"into"`: integrated (original) intensity of the peak.
#' - `"intb"`: per-peak baseline corrected integrated peak intensity.
#' - `"maxo"`: maximum (original) intensity of the peak.
#' - `"sn"`: signal to noise ratio of the peak defined as
#'   `(maxo - baseline)/sd` with `sd` being the standard defiatio of the local
#'   chromatographic noise.
#'
#' Additional columns for `verboseColumns = TRUE`:
#' 
#' - `"mu"`: gaussian parameter mu.
#' - `"sigma"`: gaussian parameter sigma.
#' - `"h"`: gaussian parameter h.
#' - `"f"`: region number of the m/z ROI where the peak was localized.
#' - `"dppm"`: m/z deviation of mass trace across scans in ppm (always `NA`).
#' - `"scale"`: scale on which the peak was localized.
#' - `"scpos"`: peak position found by wavelet analysis (index in `int`).
#' - `"scmin"`: left peak limit found by wavelet analysis (index in `int`).
#' - `"scmax"`: right peak limit found by wavelet analysis (index in `int`).
#'
#' @md
#'
#' @examples
#'
#' ## Reading a file
#' library(xcms)
#' od <- readMSData(system.file("cdf/KO/ko15.CDF", package = "faahKO"),
#'     mode = "onDisk")
#'
#' ## Extract chromatographic data for a small m/z range
#' mzr <- c(272.1, 272.2)
#' chr <- chromatogram(od, mz = mzr)[1, 1]
#'
#' int <- intensity(chr)
#' rt <- rtime(chr)
#'
#' ## Plot the region
#' plot(chr, type = "h")
#' 
#' ## Identify peaks in the chromatographic data
#' pks <- peaksWithCentWave(intensity(chr), rtime(chr))
#' pks
#'
#' ## Highlight the peaks
#' rect(xleft = pks[, "rtmin"], xright = pks[, "rtmax"],
#'     ybottom = rep(0, nrow(pks)), ytop = pks[, "maxo"], col = "#ff000040",
#'     border = "#00000040")
peaksWithCentWave <- function(int, rt,
                              peakwidth = c(20, 50),
                              snthresh = 10,
                              prefilter = c(3, 100),
                              integrate = 1,
                              fitgauss = FALSE,
                              noise = 0, ## noise.local=TRUE,
                              verboseColumns = FALSE,
                              firstBaselineCheck = TRUE,
                              ...
                              ) {
    if (length(peakwidth) != 2)
        stop("'peakwidth' has to be a numeric of length 2")
    ## Avoid NAs
    int[is.na(int)] <- 0
    rois <- .getRtROI(int, rt, peakwidth = peakwidth, noise = noise,
                      prefilter = prefilter)
    
    basenames <- c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "into",
                   "intb", "maxo", "sn")
    verbosenames <- c("egauss", "mu", "sigma", "h", "f", "dppm", "scale",
                      "scpos", "scmin", "scmax", "lmin", "lmax")
    peaks_names <- c(basenames, verbosenames)
    peaks_ncols <- length(peaks_names)
    peakinfo_names <- c("scale", "scaleNr", "scpos", "scmin", "scmax")

    ## Peak width: seconds to scales
    scalerange <- round((peakwidth / mean(diff(rt))) / 2)

    if (length(z <- which(scalerange == 0)))
        scalerange <- scalerange[-z]
    if (length(scalerange) < 1) {
        warning("No scales? Please check peak width!")
        if (verboseColumns)
            basenames <- c(basenames, verbosenames)
        return(invisible(matrix(nrow = 0, ncol = length(basenames),
                                dimnames = list(character(), basenames))))
    }
    if (length(scalerange) > 1)
        scales <- seq(from = scalerange[1], to = scalerange[2], by = 2)
    else
        scales <- scalerange

    minPeakWidth <-  scales[1]
    noiserange <- c(minPeakWidth * 3, max(scales) * 3)
    maxGaussOverlap <- 0.5
    minPtsAboveBaseLine <- max(4, minPeakWidth - 2)
    minCentroids <- minPtsAboveBaseLine
    scRangeTol <-  maxDescOutlier <- floor(minPeakWidth / 2)
    scanrange <- c(1, length(rt))
    Nscantime <- length(int)
    ## mzdiff <- -0.001
    mzdiff <- 0

    peaklist <- NULL

    for (i in seq_len(nrow(rois))) {
        scmin <- rois[i, "scmin"]
        scmax <- rois[i, "scmax"]
        
        N <- scmax - scmin + 1
        peaks <- matrix(ncol = peaks_ncols, nrow = 0, dimnames = list(character(), peaks_names))
        peakinfo <- matrix(ncol = 5, nrow = 0, dimnames = list(character(), peakinfo_names))
        ## Could also return the "correct one..."
        sccenter <- scmin + floor(N/2) - 1
        ## sccenter <- rois[i, "sccent"]
        scrange <- c(scmin, scmax)

        ## scrange + noiserange, used for baseline detection and wavelet analysis
        sr <- c(max(scanrange[1], scrange[1] - max(noiserange)),
                min(scanrange[2], scrange[2] + max(noiserange)))
        ## Intensities for baseline detection and wavelet analysis
        td <- sr[1]:sr[2]
        d <- int[td]
        scan.range <- c(sr[1], sr[2])
        ## the ROI
        otd <- scmin:scmax
        od <- int[otd]
        if (all(od == 0)) {
            warning("centWave: no peaks found in ROI.")
            next
        }
        ## scrange + scRangeTol, used for gauss fitting and continuous
        ## data above 1st baseline detection
        ftd <- max(td[1], scrange[1] - scRangeTol):min(td[length(td)],
                                                       scrange[2] + scRangeTol)
        fd <- int[ftd]

        ## 1st type of baseline: statistic approach
        ## in case of very long mass trace use full scan range
        ## for baseline detection
        if (N >= 10 * minPeakWidth)
            noised <- int
        else
            noised <- d
        ## 90% trimmed mean as first baseline guess
        noise <- estimateChromNoise(noised, trim = 0.05,
                                    minPts = 3 * minPeakWidth)
        ## any continuous data above 1st baseline ?
        if (firstBaselineCheck &
            !continuousPtsAboveThreshold(fd, threshold = noise,
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
        if (!(!is.null(dim(wCoefs)) && any((wCoefs - baseline) >= sdthr)))
            next
        if (td[length(td)] == Nscantime) ## workaround, localMax fails otherwise
            wCoefs[nrow(wCoefs),] <- wCoefs[nrow(wCoefs) - 1, ] * 0.99
        localMax <- MSW.getLocalMaximumCWT(wCoefs)
        rL <- MSW.getRidge(localMax)
        wpeaks <- sapply(rL,
                         function(x) {
                             w <- min(1:length(x), ncol(wCoefs))
                             any((wCoefs[x,w] - baseline) >= sdthr)
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
                        if (any(d[pp[dv]] - baseline >= sdthr)) {
                            ## ## allow roiScales to be a numeric of length 0
                            ## if(length(roiScales) > 0) {
                            ##     ## use given scale
                            ##     best.scale.nr <- which(scales == roiScales[[f]])
                            ##     if(best.scale.nr > length(opp))
                            ##         best.scale.nr <- length(opp)
                            ## } else {
                            ## try to decide which scale describes the peak best
                            inti <- numeric(length(opp))
                            irange <- rep(ceiling(scales[1]/2), length(opp))
                            for (k in 1:length(opp)) {
                                kpos <- opp[k]
                                r1 <- ifelse(kpos - irange[k] > 1,
                                             kpos - irange[k], 1)
                                r2 <- ifelse(kpos + irange[k] < length(d),
                                             kpos + irange[k], length(d))
                                inti[k] <- sum(d[r1:r2])
                            }
                            maxpi <- which.max(inti)
                            if (length(maxpi) > 1) {
                                m <- wCoefs[opp[maxpi], maxpi]
                                bestcol <- which(m == max(m),
                                                 arr.ind = TRUE)[2]
                                best.scale.nr <- maxpi[bestcol]
                            } else best.scale.nr <- maxpi

                            best.scale <-  scales[best.scale.nr]
                            best.scale.pos <- opp[best.scale.nr]

                            pprange <- min(pp):max(pp)
                            lwpos <- max(1, best.scale.pos - best.scale)
                            rwpos <- min(best.scale.pos + best.scale, length(td))
                            p1 <- match(td[lwpos], otd)[1]
                            p2 <- match(td[rwpos], otd)
                            p2 <- p2[length(p2)]
                            if (is.na(p1)) p1 <- 1
                            if (is.na(p2)) p2 <- N
                            maxint <- max(od[p1:p2])
                            
                            peaks <- rbind(
                                peaks,
                                c(1, 1, 1,    # mz, mzmin, mzmax,
                                  NA, NA, NA, # rt, rtmin, rtmax,
                                  NA,         # intensity (sum)
                                  NA,         # intensity (-bl)
                                  maxint,     # max intensity
                                  round((maxint - baseline) / sdnoise), # S/N Ratio
                                  NA,       # Gaussian RMSE
                                  NA,NA,NA, # Gaussian Parameters
                                  i,        # ROI Position
                                  NA, # max. difference between the [minCentroids] peaks in ppm
                                  best.scale, # Scale
                                  td[best.scale.pos],
                                  td[lwpos],
                                  td[rwpos], # Peak positions guessed from the wavelet's (scan nr)
                                  NA, NA))   # Peak limits (scan nr)
                            peakinfo <- rbind(
                                peakinfo,
                                c(best.scale, best.scale.nr,
                                  best.scale.pos, lwpos, rwpos))
                            ## Peak positions guessed from the wavelet's
                        }
                    }
                }
            }                           # for (p in 1:length(wpeaksidx))
        }                               # if (any(wpeaks))

        ##  postprocessing
        for (p in seq_len(nrow(peaks))) {
            ## find minima, assign rt and intensity values
            if (integrate == 1) {
                lm <- descendMin(wCoefs[, peakinfo[p, "scaleNr"]],
                                 istart = peakinfo[p, "scpos"])
                gap <- all(d[lm[1]:lm[2]] == 0) # looks like we got stuck in a gap right in the middle of the peak
                if ((lm[1] == lm[2]) || gap )   # fall-back
                    lm <- descendMinTol(d, startpos = c(peakinfo[p, "scmin"],
                                                        peakinfo[p, "scmax"]),
                                        maxDescOutlier)
            } else {
                lm <- descendMinTol(d, startpos = c(peakinfo[p, "scmin"],
                                                    peakinfo[p, "scmax"]),
                                    maxDescOutlier)
            }
            ## narrow down peak rt boundaries by skipping zeros
            lm_range <- lm[1]:lm[2]
            pd <- d[lm_range]
            np <- length(pd)
            lm.l <- findEqualGreaterUnsorted(pd, 1)
            lm.l <- max(1, lm.l - 1)
            lm.r <- findEqualGreaterUnsorted(rev(pd), 1)
            lm.r <- max(1, lm.r - 1)
            lm <- lm + c(lm.l - 1, -(lm.r - 1) )
            
            peakrange <- td[lm]
            peaks[p, "rtmin"] <- rt[peakrange[1]]
            peaks[p, "rtmax"] <- rt[peakrange[2]]
            peaks[p, "maxo"] <- max(pd)
            pwid <- (rt[peakrange[2]] - rt[peakrange[1]]) /
                (peakrange[2] - peakrange[1])
            if (is.na(pwid))
                pwid <- 1
            peaks[p, "into"] <- pwid * sum(pd)
            db <- pd - baseline
            peaks[p, "intb"] <- pwid * sum(db[db>0])
            peaks[p, "lmin"] <- lm[1]
            peaks[p, "lmax"] <- lm[2]
            
            if (fitgauss) {
                ## perform gaussian fits, use wavelets for inital parameters
                td_lm <- td[lm_range]
                md <- max(pd)
                d1 <- pd / md ## normalize data for gaussian error calc.
                pgauss <- fitGauss(td_lm, pd,
                                   pgauss = list(mu = peaks[p, "scpos"],
                                                 sigma = peaks[p, "scmax"] -
                                                     peaks[p, "scmin"],
                                                 h = peaks[p, "maxo"]))
                rtime <- peaks[p, "scpos"]
                if (!any(is.na(pgauss)) && all(pgauss > 0)) {
                    gtime <- td[match(round(pgauss$mu), td)]
                    if (!is.na(gtime)) {
                        rtime <- gtime
                        peaks[p, "mu"] <- pgauss$mu
                        peaks[p, "sigma"] <- pgauss$sigma
                        peaks[p, "h"] <- pgauss$h
                        peaks[p,"egauss"] <- sqrt(
                        (1 / length(td_lm)) *
                        sum(((d1 - gauss(td_lm, pgauss$h / md,
                                         pgauss$mu, pgauss$sigma))^2)))
                    }
                }
                peaks[p, "rt"] <- rt[rtime]
                ## avoid fitting side effects
                if (peaks[p, "rt"] < peaks[p, "rtmin"])
                    peaks[p, "rt"] <- rt[peaks[p, "scpos"]]
            } else
                peaks[p, "rt"] <- rt[peaks[p, "scpos"]]
        }   # end for (p in seq_len(nrow(peaks)))
        peaks <- joinOverlappingPeaks(td, d, otd, rep(1:length(otd)), od, rt,
                                      scan.range, peaks, maxGaussOverlap,
                                      mzCenterFun = mzCenter.wMean)
        
        if (!is.null(peaks))
            peaklist[[length(peaklist) + 1]] <- peaks
    }                                   # end of for (i in seq_len(nrow(rois)))
    
    if (length(peaklist) == 0) {
        warning("No peaks found!")
        if (verboseColumns)
            nopeaks <- matrix(nrow = 0, ncol = peans_ncols,
                              dimnames = list(character(), peaks_names))
        else
            nopeaks <- matrix(nrow = 0, ncol = length(basenames),
                              dimnames = list(character(), basenames))
        return(nopeaks)
    }
    p <- do.call(rbind, peaklist)
    if (!verboseColumns)
        p <- p[, basenames, drop = FALSE]

    uorder <- order(p[, "into"], decreasing = TRUE)
    pm <- p[, c("mzmin", "mzmax", "rtmin", "rtmax"), drop = FALSE]
    uindex <- rectUnique(pm, uorder, mzdiff, ydiff = -0.00001) # allow adjacent peaks
    p[uindex, -(1:3), drop = FALSE]
}


#' @description
#'
#' Identify regions in the intensity vector in which we might expect a peak.
#' In contrast to the original `.Call("findmzROI")` we base the search
#' exclusively on the intensity values, first identifying local maxima in a
#' moving window approach based on the lower expected peak width
#' (`peakwidth[1]`) and subsequently defining regions with width equal to
#' `peakwidth[2]` centered around the local maxima.
#'
#' @param int `numeric` with the intensities. Should **not** contain `NA`s!
#'
#' @param rt `numeric` with the retention times.
#'
#' @param peakwidth `numeric(2)` with the lowe and upper bound for the expected
#'     peak widths.
#' 
#' @param prefilter `numeric(2)` (`c(k, I)`): only regions of interest with at
#'     least `k` centroids with signal `>= I` are returned.
#'
#' @param noise `numeric(1)` defining the minimum required intensity for
#'     centroids to be considered in the first analysis step.
#'
#' @return `matrix` with two columns `"scmin"`, `"scmax"` and `"sccent"` with
#'     the index of (lower and upper) bound defining the region of interest and
#'     the position of the center.
#' 
#' @md
#'
#' @noRd
#'
#' @examples
#'
#' library(xcms)
#' od <- readMSData(system.file("cdf/KO/ko15.CDF", package = "faahKO"),
#'     mode = "onDisk")
#'
#' ## Extract chromatographic data for a small m/z range
#' chr <- chromatogram(od, mz = c(272.1, 272.3))[1, 1]
#'
#' int <- intensity(chr)
#' int[is.na(int)] <- 0
#' rt <- rtime(chr)
#' .getRtROI(int, rt)
.getRtROI <- function(int, rt, peakwidth = c(20, 50), noise = 0,
                      prefilter = c(3, 100)) {
    peakwidth <- range(peakwidth)
    if (length(prefilter) != 2)
        stop("'prefilter' has to be a 'numeric' of length 2")
    int_len <- length(int)
    if (int_len != length(rt))
        stop("lengths of 'int' and 'rt' have to match")
    ## rt to halfWindowSize:
    rt_step <- mean(diff(rt), na.rm = TRUE)
    up_bound <- ceiling(peakwidth[2] / rt_step)
    pk_idx <- which(MALDIquant:::.localMaxima(int, floor(peakwidth[1] /
                                                         rt_step)))
    ## First filter: int > noise
    pk_idx <- pk_idx[int[pk_idx] >= noise]
    if (!length(pk_idx))
        return(matrix(ncol = 2, nrow = 0))
    ## Define the ROIs
    scmin <- sapply(pk_idx - up_bound, max, y = 1)
    scmax <- sapply(pk_idx + up_bound, min, y = int_len)
    ## Second filter: at least k values larger I
    roi_idxs <- mapply(scmin, scmax, FUN = seq)
    ok <- vapply(roi_idxs,
                 FUN = function(x, k, I) {
                     sum(int[x] >= I) >= k
                 },
                 FUN.VALUE = logical(1),
                 k = prefilter[1], I = prefilter[2],
                 USE.NAMES = FALSE)
    if (any(ok))
        cbind(scmin = scmin[ok], scmax = scmax[ok], sccent = pk_idx[ok])
    else
        matrix(ncol = 3, nrow = 0)
}
