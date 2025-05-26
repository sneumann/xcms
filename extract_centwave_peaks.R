# This script contains functions extracted from the xcms package
# for performing centWave peak detection. It is intended for use
# in environments where the full xcms package installation is not
# feasible or desired.
#
# NOTE: This script relies on C functions that are part of the xcms
# package. An xcms installation is therefore required for these
# functions to work correctly. Some R utility functions might also
# be implicitly expected by the copied functions.
#
# Original source: xcms package (https://github.com/sneumann/xcms)
#
# Please ensure that the xcms package is installed and accessible in
# your R environment when using this script.
#
# Functions included:
# - do_findChromPeaks_centWave
# - .centWave_orig
# - .narrow_rt_boundaries
# - .get_beta_values
# - .scale_zero_one
# - MSW.cwt
# - MSW.extendNBase
# - MSW.extendLength
# - MSW.getLocalMaximumCWT
# - MSW.localMaximum
# - MSW.getRidge
# - running
# - invalid
# - odd
# - gauss
# - fitGauss
# - joinOverlappingPeaks
# - descendMinTol
# - gaussCoverage
# - mzCenter.wMean
# - mzCenter.mean
# - mzCenter.apex
# - mzCenter.wMeanApex3
# - mzCenter.meanApex3
# - trimm
# - estimateChromNoise
# - getLocalNoiseEstimate
# - continuousPtsAboveThresholdIdx
# - valueCount2ScanIndex
#
# Potential TODOs for missing functions if not found in source files:

#
# Compiling the accompanying C code:
#
# The C functions necessary for this script to run independently of a full
# xcms installation need to be compiled into a shared library.
#
# Prerequisites:
# 1. R installed.
# 2. A C compiler toolchain:
#    - Windows: Rtools.
#    - macOS: Xcode Command Line Tools.
#    - Linux: GCC and R development headers (e.g., r-base-dev).
#
# Steps:
# 1. Create a directory (e.g., xcms_c_src).
# 2. Place the following C source and header files (which should be
#    provided alongside this R script) into that directory:
#    - mzROI_xcms_extracted.c
#    - util_xcms_extracted.c
#    - init_xcms_extracted.c
#    - util_xcms_extracted.h
# 3. Open a terminal or command prompt.
# 4. Navigate into the directory created in step 1 (e.g., cd path/to/xcms_c_src).
# 5. Run the R compilation command:
#    R CMD SHLIB mzROI_xcms_extracted.c util_xcms_extracted.c init_xcms_extracted.c
#
# This will produce a shared library file (e.g., mzROI_xcms_extracted.so on
# Linux/macOS, or mzROI_xcms_extracted.dll on Windows). The exact name might
# vary slightly based on your OS and R's conventions (often named after the
# first .c file provided to R CMD SHLIB).
# This shared library needs to be loaded in the R script using dyn.load().
# See the section below on loading this library.
#

# Implementation based on common peak boundary finding logic,
# as the original xcms source for this specific helper was not located.
descendMin <- function(d, istart) {
    N <- length(d)
    
    # Ensure istart is within valid bounds
    if (istart < 1) istart <- 1
    if (istart > N) istart <- N

    l <- istart
    # Descend left: go left as long as the next point is smaller or equal
    # Or until we hit the start of the vector.
    while (l > 1 && d[l-1] <= d[l]) {
        l <- l - 1
    }

    r <- istart
    # Descend right: go right as long as the next point is smaller or equal
    # Or until we hit the end of the vector.
    while (r < N && d[r+1] <= d[r]) {
        r <- r + 1
    }
    
    return(as.integer(c(l, r)))
}

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
                                       sleep = 0,
                                       extendLengthMSW = FALSE,
                                       verboseBetaColumns = FALSE) {
    if (getOption("originalCentWave", default = TRUE)) {
        ## message("DEBUG: using original centWave.")
        .centWave_orig(mz = mz, int = int, scantime = scantime,
                       valsPerSpect = valsPerSpect, ppm = ppm, peakwidth = peakwidth,
                       snthresh = snthresh, prefilter = prefilter,
                       mzCenterFun = mzCenterFun, integrate = integrate,
                       mzdiff = mzdiff, fitgauss = fitgauss, noise = noise,
                       verboseColumns = verboseColumns, roiList = roiList,
                       firstBaselineCheck = firstBaselineCheck,
                       roiScales = roiScales, sleep = sleep,
                       extendLengthMSW = extendLengthMSW,
                       verboseBetaColumns = verboseBetaColumns)
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

.centWave_orig <- function(mz, int, scantime, valsPerSpect,
                           ppm = 25, peakwidth = c(20,50), snthresh = 10,
                           prefilter = c(3,100), mzCenterFun = "wMean",
                           integrate = 1, mzdiff = -0.001, fitgauss = FALSE,
                           noise = 0, ## noise.local=TRUE,
                           sleep = 0, verboseColumns = FALSE, roiList = list(),
                           firstBaselineCheck = TRUE, roiScales = NULL,
                           extendLengthMSW = FALSE, verboseBetaColumns = FALSE) {
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
    betanames <- c("beta_cor", "beta_snr")

    ## Peak width: seconds to scales
    scalerange <- round((peakwidth / mean(diff(scantime))) / 2)

    if (length(z <- which(scalerange == 0)))
        scalerange <- scalerange[-z]
    if (length(scalerange) < 1) {
        warning("No scales? Please check peak width!")
      matrix_length <- length(basenames)
      matrix_names <- basenames
      if (verboseColumns) {
        matrix_length <- matrix_length + length(verbosenames)
        matrix_names <- c(matrix_names, verbosenames)
      }
      if (verboseBetaColumns) {
        matrix_length <- matrix_length + length(betanames)
        matrix_names <- c(matrix_names, betanames)
      }
      nopeaks <- matrix(nrow = 0, ncol = matrix_length)
      colnames(nopeaks) <- matrix_names
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
                    roiList <- xcms:::.Call("findmzROI",
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
                    roiList <<- xcms:::.Call("findmzROI",
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
          matrix_length <- length(basenames)
          matrix_names <- basenames
          if (verboseColumns) {
            matrix_length <- matrix_length + length(verbosenames)
            matrix_names <- c(matrix_names, verbosenames)
          }
          if (verboseBetaColumns) {
            matrix_length <- matrix_length + length(betanames)
            matrix_names <- c(matrix_names, betanames)
          }
          nopeaks <- matrix(nrow = 0, ncol = matrix_length)
          colnames(nopeaks) <- matrix_names
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
        eic <- xcms:::.Call("getEIC", mz, int, scanindex, as.double(mzrange),
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
        omz <- xcms:::.Call("getWeightedMZ", mz, int, scanindex, as.double(mzrange),
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
            noised <- xcms:::.Call("getEIC", mz, int, scanindex, as.double(mzrange),
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
        if (firstBaselineCheck &&
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
        wCoefs <- MSW.cwt(d, scales = scales, wavelet = 'mexh',
                          extendLengthMSW = extendLengthMSW)
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
                                             NA, NA,     ## Peak limits (scan nr)
                                             NA, NA))    ## Beta fitting values
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
            colnames(peaks) <- c(basenames, verbosenames, betanames)
            colnames(peakinfo) <- c("scale", "scaleNr", "scpos",
                                    "scmin", "scmax")
            for (p in 1:dim(peaks)[1]) {
                ## find minima (peak boundaries), assign rt and intensity values
                if (integrate == 1) {
                    lm <- descendMin(wCoefs[, peakinfo[p, "scaleNr"]],
                                     istart = peakinfo[p, "scpos"])
                    gap <- all(d[lm[1]:lm[2]] == 0) # looks like we got stuck in a gap right in the middle of the peak
                    if ((lm[1] == lm[2]) || gap)   # fall-back
                        lm <- descendMinTol(
                            d, startpos = c(peakinfo[p, "scmin"],
                                            peakinfo[p, "scmax"]),
                            maxDescOutlier)
                } else {
                    lm <- descendMinTol(d, startpos = c(peakinfo[p, "scmin"],
                                                        peakinfo[p, "scmax"]),
                                        maxDescOutlier)
                }
                ## Narrow peak rt boundaries by removing values below threshold
                lm <- .narrow_rt_boundaries(lm, d)
                lm_seq <- lm[1]:lm[2]
                pd <- d[lm_seq]

                # Implement a fit of a skewed gaussian (beta distribution)
                # for peak shape and within-peak signal-to-noise ratio
                # See https://doi.org/10.1186/s12859-023-05533-4 and
                # https://github.com/sneumann/xcms/pull/685
                if(verboseBetaColumns){
                  peaks[p, c("beta_cor", "beta_snr")] <- .get_beta_values(pd)
                }

                peakrange <- td[lm]
                peaks[p, "rtmin"] <- scantime[peakrange[1]]
                peaks[p, "rtmax"] <- scantime[peakrange[2]]
                peaks[p, "maxo"] <- max(pd)
                pwid <- (scantime[peakrange[2]] - scantime[peakrange[1]]) /
                    (peakrange[2] - peakrange[1])
                if (is.na(pwid))
                    pwid <- 1
                peaks[p, "into"] <- pwid * sum(pd)
                db <- pd - baseline
                peaks[p, "intb"] <- pwid * sum(db[db > 0])
                peaks[p, "lmin"] <- lm[1]
                peaks[p, "lmax"] <- lm[2]

                if (fitgauss) {
                    ## perform gaussian fits, use wavelets for inital parameters
                    td_lm <- td[lm_seq]
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
                               rtrange = trange, log = TRUE)
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
      matrix_length <- length(basenames)
      matrix_names <- basenames
      if (verboseColumns) {
        matrix_length <- matrix_length + length(verbosenames)
        matrix_names <- c(matrix_names, verbosenames)
      }
      if (verboseBetaColumns) {
        matrix_length <- matrix_length + length(betanames)
        matrix_names <- c(matrix_names, betanames)
      }
      nopeaks <- matrix(nrow = 0, ncol = matrix_length)
      colnames(nopeaks) <- matrix_names
      message(" FAIL: none found!")
      return(nopeaks)
    }
    p <- do.call(rbind, peaklist)
    keepcols <- basenames
    if (verboseColumns){
      keepcols <- c(keepcols, verbosenames)
    }
    if(verboseBetaColumns){
      keepcols <- c(keepcols, betanames)
    }
    p <- p[, keepcols, drop = FALSE]
    uorder <- order(p[, "into"], decreasing = TRUE)
    pm <- as.matrix(p[,c("mzmin", "mzmax", "rtmin", "rtmax"), drop = FALSE])
    uindex <- xcms:::.Call("rectUnique", pm, uorder, mzdiff, -0.00001, PACKAGE = "xcms") ## allow adjacent peaks
    pr <- p[uindex, , drop = FALSE]
    message(" OK: ", nrow(pr), " found.")

    return(pr)
}

.narrow_rt_boundaries <- function(lm, d, thresh = 1) {
    lm_seq <- lm[1]:lm[2]
    above_thresh <- d[lm_seq] >= thresh
    if (any(above_thresh)) {
        ## Expand by one on each side to be consistent with old code.
        above_thresh <- above_thresh | c(above_thresh[-1], FALSE) |
            c(FALSE, above_thresh[-length(above_thresh)])
        lm <- range(lm_seq[above_thresh], na.rm = TRUE)
    }
    lm
}

.get_beta_values <- function(intensity, rtime = seq_along(intensity),
                             skews=c(3, 3.5, 4, 4.5, 5), zero.rm = TRUE){
  if (zero.rm) {
    ## remove 0 or NA intensities
    keep <- which(intensity > 0)
    rtime <- rtime[keep]
    intensity <- intensity[keep]
  }
  if (length(intensity) < 5) {
    best_cor <- NA
    beta_snr <- NA
  } else {
    beta_sequence <- rep(.scale_zero_one(rtime), each = length(skews))
    beta_vals <- t(matrix(dbeta(beta_sequence, shape1 = skews, shape2 = 5),
                          nrow = length(skews)))
    # matplot(beta_vals)
    beta_cors <- cor(intensity, beta_vals)
    best_cor <- max(beta_cors)
    best_curve <- beta_vals[, which.max(beta_cors)]
    noise_level <- sd(diff(.scale_zero_one(best_curve) -
                           .scale_zero_one(intensity)))
    beta_snr <- log10(max(intensity) / noise_level)
  }
  c(best_cor = best_cor, beta_snr = beta_snr)
}

.scale_zero_one <- function(num_vec){
  (num_vec-min(num_vec)) / (max(num_vec) - min(num_vec))
}

MSW.cwt <- function (ms, scales = 1, wavelet = "mexh", extendLengthMSW = FALSE)
{ ## modified from package MassSpecWavelet
    if (wavelet == "mexh") {
        psi_xval <- seq(-6, 6, length = 256)
        psi <- (2/sqrt(3) * pi^(-0.25)) * (1 - psi_xval^2) *
            exp(-psi_xval^2/2)
    }
    else if (is.matrix(wavelet)) {
        if (nrow(wavelet) == 2) {
            psi_xval <- wavelet[1, ]
            psi <- wavelet[2, ]
        }
        else if (ncol(wavelet) == 2) {
            psi_xval <- wavelet[, 1]
            psi <- wavelet[2, ]
        }
        else {
            stop("Unsupported wavelet format!")
        }
    }
    else {
        stop("Unsupported wavelet!")
    }
    oldLen <- length(ms)
    # IF extendLengthMSW is TRUE:
    # The new length is determined by the scales argument, so a larger peakwidth
    # will ensure more scales are run, but may slow it down. See 
    # https://github.com/sneumann/xcms/issues/445 for more information about
    # a change from using extendNBase to extendLength.
    if(extendLengthMSW){
        newLen <- 2^(ceiling(log2(max(scales)*12)))
        ms <- MSW.extendLength(x = ms, addLength = (newLen-length(ms)), 
                               method = "open")
    } else {
        ms <- MSW.extendNBase(ms, nLevel = NULL, base = 2)
    }
    
    
    len <- length(ms)
    nbscales <- length(scales)
    wCoefs <- NULL
    psi_xval <- psi_xval - psi_xval[1]
    dxval <- psi_xval[2]
    xmax <- psi_xval[length(psi_xval)]
    for (i in 1:length(scales)) {
        scale.i <- scales[i]
        f <- rep(0, len)
        j <- 1 + floor((0:(scale.i * xmax))/(scale.i * dxval))
        if (length(j) == 1)
            j <- c(1, 1)
        lenWave <- length(j)
        f[1:lenWave] <- rev(psi[j]) - mean(psi[j])
        if (length(f) > len)
        {i<-i-1;break;}   ##  stop(paste("scale", scale.i, "is too large!"))
        wCoefs.i <- 1/sqrt(scale.i) * convolve(ms, f)
        wCoefs.i <- c(wCoefs.i[(len - floor(lenWave/2) + 1):len],
                      wCoefs.i[1:(len - floor(lenWave/2))])
        wCoefs <- cbind(wCoefs, wCoefs.i)
    }
    if (i < 1) return(NA)
    scales <- scales[1:i]
    if (length(scales) == 1)
        wCoefs <- matrix(wCoefs, ncol = 1)
    colnames(wCoefs) <- scales
    wCoefs <- wCoefs[1:oldLen, , drop = FALSE]
    wCoefs
}

MSW.extendNBase <- function(x, nLevel=1, base=2, ...)
{ ## from package MassSpecWavelet
    if (!is.matrix(x)) x <- matrix(x, ncol=1)

    nR <- nrow(x)
    if (is.null(nLevel)) {
        nR1 <- nextn(nR, base)
    } else {
        nR1 <- ceiling(nR / base^nLevel) * base^nLevel
    }
    if (nR != nR1) {
        x <- MSW.extendLength(x, addLength=nR1-nR, ...)
    }
    x
}

MSW.extendLength <-
    function(x, addLength=NULL, method=c('reflection', 'open', 'circular'), direction=c('right', 'left', 'both'))
{       ## from package MassSpecWavelet
    if (is.null(addLength)) stop('Please provide the length to be added!')
    if (!is.matrix(x)) x <- matrix(x, ncol=1)
    method <- match.arg(method)
    direction <- match.arg(direction)

    nR <- nrow(x)
    nR1 <- nR + addLength
    if (direction == 'both') {
        left <- right <- addLength
    } else if (direction == 'right') {
        left <- 0
        right <- addLength
    } else if (direction == 'left') {
        left <- addLength
        right <- 0
    }

    if (right > 0) {
        x <- switch(method,
                    reflection =rbind(x, x[nR:(2 * nR - nR1 + 1), , drop=FALSE]),
                    open = rbind(x, matrix(rep(x[nR,], addLength), ncol=ncol(x), byrow=TRUE)),
                    circular = rbind(x, x[1:(nR1 - nR),, drop=FALSE]))
    }

    if (left > 0) {
        x <- switch(method,
                    reflection =rbind(x[addLength:1, , drop=FALSE], x),
                    open = rbind(matrix(rep(x[1,], addLength), ncol=ncol(x), byrow=TRUE), x),
                    circular = rbind(x[(2 * nR - nR1 + 1):nR,, drop=FALSE], x))
    }
    if (ncol(x) == 1)  x <- as.vector(x)

    x
}

MSW.getLocalMaximumCWT <-
    function(wCoefs, minWinSize=5, amp.Th=0)
{        ## from package MassSpecWavelet
    localMax <- NULL
    scales <- as.numeric(colnames(wCoefs))

    for (i in 1:length(scales)) {
        scale.i <- scales[i]
        winSize.i <- scale.i * 2 + 1
        if (winSize.i < minWinSize) {
            winSize.i <- minWinSize
        }
        temp <- MSW.localMaximum(wCoefs[,i], winSize.i)
        localMax <- cbind(localMax, temp)
    }
    ## Set the values less than peak threshold as 0
    localMax[wCoefs < amp.Th] <- 0
    colnames(localMax) <- colnames(wCoefs)
    rownames(localMax) <- rownames(wCoefs)
    localMax
}

MSW.localMaximum <-
    function (x, winSize = 5)
{   ## from package MassSpecWavelet
    len <- length(x)
    rNum <- ceiling(len/winSize)

    ## Transform the vector as a matrix with column length equals winSize
    ##		and find the maximum position at each row.
    y <- matrix(c(x, rep(x[len], rNum * winSize - len)), nrow=winSize)
    y.maxInd <- apply(y, 2, which.max)
    ## Only keep the maximum value larger than the boundary values
    selInd <- which(apply(y, 2, function(x) max(x) > x[1] & max(x) > x[winSize]))

    ## keep the result
    localMax <- rep(0, len)
    localMax[(selInd-1) * winSize + y.maxInd[selInd]] <- 1

    ## Shift the vector with winSize/2 and do the same operation
    shift <- floor(winSize/2)
    rNum <- ceiling((len + shift)/winSize)
    y <- matrix(c(rep(x[1], shift), x, rep(x[len], rNum * winSize - len - shift)), nrow=winSize)
    y.maxInd <- apply(y, 2, which.max)
    ## Only keep the maximum value larger than the boundary values
    selInd <- which(apply(y, 2, function(x) max(x) > x[1] & max(x) > x[winSize]))
    localMax[(selInd-1) * winSize + y.maxInd[selInd] - shift] <- 1

    ## Check whether there is some local maxima have in between distance less than winSize
    maxInd <- which(localMax > 0)
    selInd <- which(diff(maxInd) < winSize)
    if (length(selInd) > 0) {
        selMaxInd1 <- maxInd[selInd]
        selMaxInd2 <- maxInd[selInd + 1]
        temp <- x[selMaxInd1] - x[selMaxInd2]
        localMax[selMaxInd1[temp <= 0]] <- 0
        localMax[selMaxInd2[temp > 0]] <- 0
    }

    localMax
}

MSW.getRidge <-
    function(localMax, iInit=ncol(localMax), step=-1, iFinal=1, minWinSize=3, gapTh=3, skip=NULL)
{  ## modified from package MassSpecWavelet

    scales <- as.numeric(colnames(localMax))
    if (is.null(scales))  scales <- 1:ncol(localMax)

    maxInd_curr <- which(localMax[, iInit] > 0)
    nMz <- nrow(localMax)

    if (is.null(skip))	{
        skip <- iInit + 1
    }

    ## Identify all the peak pathes from the coarse level to detail levels (high column to low column)
    ## Only consider the shortest path
    if ( ncol(localMax) > 1 ) colInd <- seq(iInit+step, iFinal, step)
    else colInd <- 1
    ridgeList <- as.list(maxInd_curr)
    names(ridgeList) <- maxInd_curr
    peakStatus <- as.list(rep(0, length(maxInd_curr)))
    names(peakStatus) <- maxInd_curr

    ## orphanRidgeList keep the ridges disconnected at certain scale level
    ## Changed by Pan Du 05/11/06
    orphanRidgeList <- NULL
    orphanRidgeName <- NULL
    nLevel <- length(colInd)

    for (j in 1:nLevel) {
        col.j <- colInd[j]
        scale.j <- scales[col.j]

        if (colInd[j] == skip) {
            oldname <- names(ridgeList)
            ridgeList <- lapply(ridgeList, function(x) c(x, x[length(x)]))
            ##peakStatus <- lapply(peakStatus, function(x) c(x, x[length(x)]))
            names(ridgeList) <- oldname
            ##names(peakStatus) <- oldname
            next
        }

        if (length(maxInd_curr) == 0) {
            maxInd_curr <- which(localMax[, col.j] > 0)
            next
        }

        ## The slide window size is proportional to the CWT scale
        ## winSize.j <- scale.j / 2 + 1
        winSize.j <- floor(scale.j/2)
        if (winSize.j < minWinSize) {
            winSize.j <- minWinSize
        }

        selPeak.j <- NULL
        remove.j <- NULL
        for (k in 1:length(maxInd_curr)) {
            ind.k <- maxInd_curr[k]
            start.k <- ifelse(ind.k-winSize.j < 1, 1, ind.k-winSize.j)
            end.k <- ifelse(ind.k+winSize.j > nMz, nMz, ind.k+winSize.j)
            ind.curr <- which(localMax[start.k:end.k, col.j] > 0) + start.k - 1
            ##ind.curr <- which(localMax[, col.j] > 0)
            if (length(ind.curr) == 0) {
                status.k <- peakStatus[[as.character(ind.k)]]
                ## bug  work-around
                if (is.null(status.k)) status.k <- gapTh +1
                ##
                if (status.k > gapTh & scale.j >= 2) {
                    temp <- ridgeList[[as.character(ind.k)]]
                    orphanRidgeList <- c(orphanRidgeList, list(temp[1:(length(temp)-status.k)]))
                    orphanRidgeName <- c(orphanRidgeName, paste(col.j + status.k + 1, ind.k, sep='_'))
                    remove.j <- c(remove.j, as.character(ind.k))
                    next
                } else {
                    ind.curr <- ind.k
                    peakStatus[[as.character(ind.k)]] <- status.k + 1
                }
            } else {
                peakStatus[[as.character(ind.k)]] <- 0
                if (length(ind.curr) >= 2)  ind.curr <- ind.curr[which.min(abs(ind.curr - ind.k))]
            }
            ridgeList[[as.character(ind.k)]] <- c(ridgeList[[as.character(ind.k)]], ind.curr)
            selPeak.j <- c(selPeak.j, ind.curr)
        }
        ## Remove the disconnected lines from the currrent list
        if (length(remove.j) > 0) {
            removeInd <- which(names(ridgeList) %in% remove.j)
            ridgeList <- ridgeList[-removeInd]
            peakStatus <- peakStatus[-removeInd]
        }

        ## Check for duplicated selected peaks and only keep the one with the longest path.
        dupPeak.j <- unique(selPeak.j[duplicated(selPeak.j)])
        if (length(dupPeak.j) > 0) {
            removeInd <- NULL
            for (dupPeak.jk in dupPeak.j) {
                selInd <- which(selPeak.j == dupPeak.jk)
                selLen <- sapply(ridgeList[selInd], length)
                removeInd.jk <- which.max(selLen)
                removeInd <- c(removeInd, selInd[-removeInd.jk])
                orphanRidgeList <- c(orphanRidgeList, ridgeList[removeInd.jk])
                orphanRidgeName <- c(orphanRidgeName, paste(col.j, selPeak.j[removeInd.jk], sep='_'))
            }
            selPeak.j <- selPeak.j[-removeInd]
            ridgeList <- ridgeList[-removeInd]
            peakStatus <- peakStatus[-removeInd]
        }

        ## Update the names of the ridgeList as the new selected peaks
        ##if (scale.j >= 2) {
        if (length(ridgeList) > 0) names(ridgeList) <- selPeak.j
        if (length(peakStatus) > 0) names(peakStatus) <- selPeak.j
        ##}

        ## If the level is larger than 3, expand the peak list by including other unselected peaks at that level
        if (scale.j >= 2) {
            maxInd_next <- which(localMax[, col.j] > 0)
            unSelPeak.j <- maxInd_next[!(maxInd_next %in% selPeak.j)]
            newPeak.j <- as.list(unSelPeak.j)
            names(newPeak.j) <- unSelPeak.j
            ## Update ridgeList
            ridgeList <- c(ridgeList, newPeak.j)
            maxInd_curr <- c(selPeak.j, unSelPeak.j)
            ## Update peakStatus
            newPeakStatus <- as.list(rep(0, length(newPeak.j)))
            names(newPeakStatus) <- newPeak.j
            peakStatus <- c(peakStatus, newPeakStatus)
        } else {
            maxInd_curr <- selPeak.j
        }
    }

    ## Attach the peak level at the beginning of the ridge names
    if (length(ridgeList) > 0) names(ridgeList) <- paste(1, names(ridgeList), sep='_')
    if (length(orphanRidgeList) > 0) names(orphanRidgeList) <- orphanRidgeName
    ## Combine ridgeList and orphanRidgeList
    ridgeList <- c(ridgeList, orphanRidgeList)
    if (length(ridgeList) == 0) return(NULL)

    ## Reverse the order as from the low level to high level.
    ridgeList <- lapply(ridgeList, rev)
    ## order the ridgeList in increasing order
    ##ord <- order(selPeak.j)
    ##ridgeList <- ridgeList[ord]

    ## Remove possible duplicated ridges
    ridgeList <- ridgeList[!duplicated(names(ridgeList))]

    attr(ridgeList, 'class') <- 'ridgeList'
    attr(ridgeList, 'scales') <- scales
    return(ridgeList)
}

running <- function (X, Y = NULL, fun = mean, width = min(length(X), 20),
                     allow.fewer = FALSE, pad = FALSE, align = c("right", "center",
                                                       "left"), simplify = TRUE, by, ...)
{   ## from package gtools
    align = match.arg(align)
    n <- length(X)
    if (align == "left") {
        from <- 1:n
        to <- pmin((1:n) + width - 1, n)
    }
    else if (align == "right") {
        from <- pmax((1:n) - width + 1, 1)
        to <- 1:n
    }
    else {
        from <- pmax((2 - width):n, 1)
        to <- pmin(1:(n + width - 1), n)
        if (!odd(width))
            stop("width must be odd for center alignment")
    }
    elements <- apply(cbind(from, to), 1, function(x) seq(x[1],
                                                          x[2]))
    if (is.matrix(elements))
        elements <- as.data.frame(elements)
    names(elements) <- paste(from, to, sep = ":")
    if (!allow.fewer) {
        len <- sapply(elements, length)
        skip <- (len < width)
    }
    else {
        skip <- 0
    }
    run.elements <- elements[!skip]
    if (!invalid(by))
        run.elements <- run.elements[seq(from = 1, to = length(run.elements),
                                         by = by)]
    if (is.null(Y)) {
        funct <- function(which, what, fun, ...) fun(what[which],
                                                     ...)
        if (simplify)
            Xvar <- sapply(run.elements, funct, what = X, fun = fun,
                           ...)
        else Xvar <- lapply(run.elements, funct, what = X, fun = fun,
                            ...)
    } else {
        funct <- function(which, XX, YY, fun, ...) fun(XX[which],
                                                       YY[which], ...)
        if (simplify)
            Xvar <- sapply(run.elements, funct, XX = X, YY = Y,
                           fun = fun, ...)
        else Xvar <- lapply(run.elements, funct, XX = X, YY = Y,
                            fun = fun, ...)
    }
    if (allow.fewer || !pad)
        return(Xvar)
    if (simplify)
        if (is.matrix(Xvar)) {
            wholemat <- matrix(new(class(Xvar[1, 1]), NA), ncol = length(to),
                               nrow = nrow(Xvar))
            colnames(wholemat) <- paste(from, to, sep = ":")
            wholemat[, -skip] <- Xvar
            Xvar <- wholemat
        }
        else {
            wholelist <- rep(new(class(Xvar[1]), NA), length(from))
            names(wholelist) <- names(elements)
            wholelist[names(Xvar)] <- Xvar
            Xvar <- wholelist
        }
    return(Xvar)
}

invalid <- function (x)
{   ## from package gtools
    if (missing(x) || is.null(x) || length(x) == 0)
        return(TRUE)
    if (is.list(x))
        return(all(sapply(x, invalid)))
    else if (is.vector(x))
        return(all(is.na(x)))
    else return(FALSE)
}

odd <- function (x) x != as.integer(x/2) * 2;

gauss <- function(x, h, mu, sigma){
    h*exp(-(x-mu)^2/(2*sigma^2))
}

fitGauss <- function(td, d, pgauss = NA) {
    if (length(d) < 3) return(rep(NA,3))
    if (!any(is.na(pgauss))) { mu <- pgauss$mu; sigma <- pgauss$sigma;h <- pgauss$h }
    fit <- try(nls(d ~ SSgauss(td,mu,sigma,h)), silent = TRUE)
    if (class(fit) == "try-error")
        fit <- try(nls(d ~ SSgauss(td, mu, sigma, h), algorithm = 'port'),
                   silent = TRUE)
    if (class(fit) == "try-error")  return(rep(NA, 3))

    as.data.frame(t(fit$m$getPars()))
}

joinOverlappingPeaks <- function(td, d, otd, omz, od, scantime, scan.range,
                                 peaks, maxGaussOverlap=0.5, mzCenterFun) {
    ## Fix issue #284: avoid having identical peaks multiple times in this
    ## matrix.
    peaks <- unique(peaks)
    gausspeaksidx <- which(!is.na(peaks[,"mu"]))
    Ngp <- length(gausspeaksidx)
    if (Ngp == 0)
        return(peaks)

    newpeaks <- NULL

    gpeaks <- peaks[gausspeaksidx, , drop = FALSE]
    if (nrow(peaks) - Ngp > 0)
        notgausspeaks <- peaks[-gausspeaksidx, , drop = FALSE]

    if (Ngp > 1) {
        comb <- which(upper.tri(matrix(0, Ngp, Ngp)), arr.ind = TRUE)
        overlap <- logical(nrow(comb))
        overlap <- rep(FALSE, dim(comb)[1])
        for (k in seq_len(nrow(comb))) {
            p1 <- comb[k, 1]
            p2 <- comb[k, 2]
            overlap[k] <- gaussCoverage(xlim = scan.range,
                                        h1 = gpeaks[p1, "h"],
                                        mu1 = gpeaks[p1, "mu"],
                                        s1 = gpeaks[p1, "sigma"],
                                        h2 = gpeaks[p2, "h"],
                                        mu2 = gpeaks[p2, "mu"],
                                        s2 = gpeaks[p2, "sigma"]) >=
                maxGaussOverlap
        }
    } else overlap <- FALSE
    
    if (any(overlap) && (Ngp > 1)) {
        jlist <- list()
        if (length(which(overlap)) > 1) {
            gm <- comb[overlap, ]
            ## create list of connected components
            cc <- list()
            cc[[1]] <- gm[1,] ## copy first entry
            for (j in 2:dim(gm)[1]) { ## search for connections
                ccl <- unlist(cc)
                nl <- sapply(cc, function(x) length(x))
                ccidx <- rep(1:length(nl),nl)
                idx <- match(gm[j,],ccl)
                if (any(!is.na(idx))) { ## connection found, add to list
                    pos <- ccidx[idx[which(!is.na(idx))[1]]]
                    cc[[pos]] <- c(cc[[pos]],gm[j,])
                } else  ## create new list element
                    cc[[length(cc) + 1]] <- gm[j,]

            }
            ccn <- list()
            lcc <- length(cc)
            ins <- rep(FALSE,lcc)
            if (lcc > 1) {
                jcomb <- which(upper.tri(matrix(0,lcc,lcc)),arr.ind = TRUE)
                for (j in 1:dim(jcomb)[1]) {
                    j1 <- jcomb[j,1]; j2 <- jcomb[j,2]
                    if (any(cc[[j1]] %in% cc[[j2]]))
                        ccn[[length(ccn) +1]] <- unique(c(cc[[j1]],cc[[j2]]))
                    else {
                        if (!ins[j1]) {
                            ccn[[length(ccn) +1]] <- unique(cc[[j1]])
                            ins[j1] <- TRUE
                        }
                        if (!ins[j2]) {
                            ccn[[length(ccn) +1]] <- unique(cc[[j2]])
                            ins[j2] <- TRUE
                        }
                    }
                }
            } else ccn <- cc;

            size <- sapply(ccn, function(x) length(x))
            s2idx <- which(size >= 2)

            if (length(s2idx) > 0) {
                for (j in 1:length(s2idx)) {
                    pgroup <- unique(ccn[[ s2idx[j] ]])
                    jlist[[j]] <- pgroup
                }
            } else stop('(length(s2idx) = 0) ?!?')
        } else jlist[[1]] <- comb[overlap, ]

        ## join all peaks belonging to one cc
        for (j in seq_along(jlist)) {
            jidx <- jlist[[j]]
            newpeak <- gpeaks[jidx[1], , drop = FALSE]
            newmin <- min(gpeaks[jidx, "lmin"])
            newmax <- max(gpeaks[jidx, "lmax"])
            newpeak[1, "scpos"] <- -1 ## not defined after join
            newpeak[1, "scmin"] <- -1 ##    ..
            newpeak[1, "scmax"] <- -1 ##    ..
            newpeak[1, "scale"] <- -1 ##    ..

            newpeak[1, "maxo"] <- max(gpeaks[jidx, "maxo"])
            newpeak[1, "sn"]   <- max(gpeaks[jidx, "sn"])
            newpeak[1, "lmin"] <- newmin
            newpeak[1, "lmax"] <- newmax
            newpeak[1, "rtmin"] <- scantime[td[newmin]]
            newpeak[1, "rtmax"] <- scantime[td[newmax]]
            newpeak[1,"rt"] <- weighted.mean(gpeaks[jidx, "rt"],
                                             w = gpeaks[jidx, "maxo"])

            ## Re-assign m/z values
            p1 <- match(td[newmin], otd)[1]
            p2 <- match(td[newmax], otd)
            p2 <- p2[length(p2)]
            if (is.na(p1)) p1 <- 1
            if (is.na(p2)) p2 <- length(omz)
            mz.value <- omz[p1:p2]
            mz.int <- od[p1:p2]

            ## re-calculate m/z value for peak range
            mzmean <- do.call(mzCenterFun, list(mz = mz.value,
                                                intensity = mz.int))
            mzrange <- range(mz.value)
            newpeak[1, "mz"] <- mzmean
            newpeak[1, c("mzmin","mzmax")] <- mzrange

            ## re-fit gaussian
            md <- max(d[newmin:newmax])
            d1 <- d[newmin:newmax] / md
            pgauss <- fitGauss(td[newmin:newmax],
                               d[newmin:newmax],
                               pgauss = list(mu = td[newmin] +
                                                 (td[newmax] - td[newmin])/2,
                                             sigma = td[newmax] - td[newmin],
                                             h = max(gpeaks[jidx, "h"])))
            if (!any(is.na(pgauss)) && all(pgauss > 0)) {
                newpeak[1, "mu"]    <- pgauss$mu
                newpeak[1, "sigma"] <- pgauss$sigma
                newpeak[1, "h"]     <- pgauss$h
                newpeak[1, "egauss"]<- sqrt((1/length(td[newmin:newmax])) *
                                            sum(((d1 - gauss(td[newmin:newmax],
                                                           pgauss$h/md,
                                                           pgauss$mu,
                                                           pgauss$sigma))^2)))
            } else { ## re-fit after join failed
                newpeak[1, "mu"]       <- NA
                newpeak[1, "sigma"]    <- NA
                newpeak[1, "h"]        <- NA
                newpeak[1, "egauss"]   <- NA
            }

            newpeaks <- rbind(newpeaks, newpeak)
        }
        ## add the remaining peaks
        jp <- unique(unlist(jlist))
        if (dim(peaks)[1] - length(jp) > 0)
            newpeaks <- rbind(newpeaks, gpeaks[-jp, ])

    } else
        newpeaks <- gpeaks

    grt.min <- newpeaks[, "rtmin"]
    grt.max <- newpeaks[, "rtmax"]

    if (nrow(peaks) - Ngp > 0) { ## notgausspeaks
        for (k in 1:nrow(notgausspeaks)) {
            ## here we can only check if they are completely overlapped
            ## by other peaks
            if (!any((notgausspeaks[k, "rtmin"] >= grt.min) &
                     (notgausspeaks[k,"rtmax"] <= grt.max)))
                newpeaks <- rbind(newpeaks,notgausspeaks[k,])
        }
    }

    rownames(newpeaks) <- NULL
    newpeaks
}

descendMinTol <- function(d,startpos,maxDescOutlier) {
    l <- startpos[1]; r <- startpos[2]; outl <- 0; N <- length(d)
    ## left
    while ((l > 1) && (d[l] > 0) && outl <= maxDescOutlier) {
        if (outl > 0) vpos <- opos else vpos <- l
        if (d[l-1] > d[vpos]) outl <- outl + 1 else outl <- 0
        if (outl == 1) opos <- l
        l <- l -1
    }
    if (outl > 0) l <- l + outl
    ## right
    outl <- 0;
    while ((r < N) && (d[r] > 0) && outl <= maxDescOutlier) {
        if (outl > 0) vpos <- opos else vpos <- r
        if (d[r+1] > d[vpos]) outl <- outl + 1 else outl <- 0
        if (outl == 1) opos <- r
        r <- r + 1
    }
    if (outl > 0) r <- r - outl
    c(l,r)
}

gaussCoverage <- function(xlim,h1,mu1,s1,h2,mu2,s2) {
    overlap <- NA
    by = 0.05
    ## Calculate points of intersection
    a <- s2^2 - s1^2
    cc <- -( 2 * s1^2 * s2^2 * (log(h1) - log(h2)) + (s1^2 * mu2^2) - (s2^2 * mu1^2) )
    b <- ((2 * s1^2 *mu2) - (2 * s2^2 * mu1))
    D <- b^2 - (a*cc)
    if (a==0) {
        S1 <- -cc/b
        S2 <- NA
    } else if ((D < 0) || ((b^2 - (4*a*cc)) < 0)) {
        S1 <- S2 <- NA
    } else {
        S1 <- (-b + sqrt(b^2 - (4*a*cc))) / (2*a)
        S2 <- (-b - sqrt(b^2 - (4*a*cc))) / (2*a)
        if (S2 < S1)
        {
            tmp <- S1
            S1 <- S2
            S2 <- tmp
        }
    }
    if (!is.na(S1)) if (S1 < xlim[1] || S1 > xlim[2]) S1 <- NA
    if (!is.na(S2)) if (S2 < xlim[1] || S2 > xlim[2]) S2 <- NA

    x <- seq(xlim[1],xlim[2],by=by)
    vsmall <- min(sum(gauss(x,h1,mu1,s1)), sum(gauss(x,h2,mu2,s2)))

    if (!is.na(S1) && !is.na(S2)) {
        x0 <- seq(xlim[1],S1,by=by)
        xo <- seq(S1,S2,by=by)
        x1 <- seq(S2,xlim[2],by=by)
        if (gauss(x0[cent(x0)],h1,mu1,s1) < gauss(x0[cent(x0)],h2,mu2,s2)) {
            ov1 <- sum(gauss(x0,h1,mu1,s1))
        } else {
            ov1 <- sum(gauss(x0,h2,mu2,s2))
        }
        if (gauss(xo[cent(xo)],h1,mu1,s1) < gauss(xo[cent(xo)],h2,mu2,s2)) {
            ov <- sum(gauss(xo,h1,mu1,s1))
        } else {
            ov <- sum(gauss(xo,h2,mu2,s2))
        }
        if (gauss(x1[cent(x1)],h1,mu1,s1) < gauss(x1[cent(x1)],h2,mu2,s2)) {
            ov2 <- sum(gauss(x1,h1,mu1,s1))
        } else {
            ov2 <- sum(gauss(x1,h2,mu2,s2))
        }
        overlap <- ov1 + ov + ov2
    } else
        if (is.na(S1) && is.na(S2)) { ## no overlap -> intergrate smaller function
            if (gauss(x[cent(x)],h1,mu1,s1) < gauss(x[cent(x)],h2,mu2,s2)) {
                overlap <- sum(gauss(x,h1,mu1,s1))
            } else {
                overlap <- sum(gauss(x,h2,mu2,s2))
            }
        } else
            if (!is.na(S1) || !is.na(S2)) {
                if (is.na(S1)) S0 <- S2 else S0 <- S1
                x0 <- seq(xlim[1],S0,by=by)
                x1 <- seq(S0,xlim[2],by=by)
                g01 <- gauss(x0[cent(x0)],h1,mu1,s1)
                g02 <- gauss(x0[cent(x0)],h2,mu2,s2)
                g11 <- gauss(x1[cent(x1)],h1,mu1,s1)
                g12 <- gauss(x1[cent(x1)],h2,mu2,s2)
                if (g01 < g02) ov1 <- sum(gauss(x0,h1,mu1,s1)) else ov1 <- sum(gauss(x0,h2,mu2,s2))
                if (g11 < g12) ov2 <- sum(gauss(x1,h1,mu1,s1)) else ov2 <- sum(gauss(x1,h2,mu2,s2))
                if ((g01 == g02) && (g01==0)) ov1 <- 0
                if ((g11 == g12) && (g11==0)) ov2 <- 0
                overlap <- ov1 + ov2
            }

    overlap / vsmall
}

mzCenter.wMean <- function(mz,intensity) {
    weighted.mean(mz, intensity)
}

mzCenter.mean <- function(mz,intensity) {
    mean(mz)
}

mzCenter.apex <- function(mz,intensity) {
    mz[which.max(intensity)]
}

mzCenter.wMeanApex3 <- function(mz,intensity) {
    iap <- which.max(intensity)
    st <- max(1,iap-1)
    en <- min(iap+1,length(mz))
    weighted.mean(mz[st:en], intensity[st:en])
}

mzCenter.meanApex3 <- function(mz,intensity) {
    iap <- which.max(intensity)
    st <- max(1,iap-1)
    en <- min(iap+1,length(mz))
    mean(mz[st:en])
}

trimm <- function(x, trim=c(0.05,0.95)) {
    a <- sort(x[x>0])
    Na <- length(a)
    quant <- round((Na*trim[1])+1):round(Na*trim[2])
    a[quant]
}

estimateChromNoise <- function(x, trim=0.05, minPts=20) {
    gz <- which(x > 0)
    if (length(gz) < minPts)
        return(mean(x))

    mean(x[gz], trim=trim)
}

getLocalNoiseEstimate <- function(d, td, ftd, noiserange, Nscantime, threshold, num) {

    if (length(d) < Nscantime) {

        ## noiserange[2] is full d-range
        drange <- which(td %in% ftd)
        n1 <- d[-drange] ## region outside the detected ROI (wide)
        n1.cp <- continuousPtsAboveThresholdIdx(n1, threshold=threshold,num=num) ## continousPtsAboveThreshold (probably peak) are subtracted from data for local noise estimation
        n1 <- n1[!n1.cp]
        if (length(n1) > 1)  {
            baseline1 <- mean(n1)
            sdnoise1 <- sd(n1)
        } else
            baseline1 <- sdnoise1 <- 1

        ## noiserange[1]
        d1 <- drange[1]
        d2 <- drange[length(drange)]
        nrange2 <- c(max(1,d1 - noiserange[1]) : d1, d2 : min(length(d),d2 + noiserange[1]))
        n2 <- d[nrange2] ## region outside the detected ROI (narrow)
        n2.cp <- continuousPtsAboveThresholdIdx(n2, threshold=threshold,num=num) ## continousPtsAboveThreshold (probably peak) are subtracted from data for local noise estimation
        n2 <- n2[!n2.cp]
        if (length(n2) > 1)  {
            baseline2 <- mean(n2)
            sdnoise2 <- sd(n2)
        } else
            baseline2 <- sdnoise2 <- 1

    } else {
        trimmed <- trimm(d,c(0.05,0.95))
        baseline1 <- baseline2 <- mean(trimmed)
        sdnoise1 <- sdnoise2 <- sd(trimmed)
    }

    c(min(baseline1,baseline2),min(sdnoise1,sdnoise2))
}

continuousPtsAboveThresholdIdx <- function(d, threshold, num) {
    above <- d >= threshold
    result <- logical(length(d))
    for (i in 1:(length(d) - num + 1)) {
        if (all(above[i:(i + num - 1)])) {
            result[i:(i + num - 1)] <- TRUE
        }
    }
    result
}

continuousPtsAboveThreshold <- function(x, threshold, num) {
    # Assumes continuousPtsAboveThresholdIdx is defined and available.
    # Returns TRUE if any sequence of 'num' points in 'x' is >= 'threshold'.
    return(sum(continuousPtsAboveThresholdIdx(x, threshold, num)) > 0)
}

# Helper function cent() needed by gaussCoverage, but not in provided list
cent <- function(x) {
    N <- length(x)
    if (N == 1) return(1)
    floor(N/2)
}

valueCount2ScanIndex <- function(valCount){
    ## Convert into 0 based.
    valCount <- cumsum(valCount)
    return(as.integer(c(0, valCount[-length(valCount)])))
}

#' @title Perform Chromatographic Peak Detection using CentWave Algorithm
#'
#' @description
#' This function is the main user-facing entry point for performing centWave-based
#' peak detection on mass spectrometry data. It processes input spectra,
#' prepares them for the core centWave algorithm, and returns the detected peaks.
#'
#' @param spectra_data The input mass spectrometry data. This can be:
#'   \itemize{
#'     \item A `data.frame` (or `data.table`, `tibble`) containing columns:
#'       \itemize{
#'         \item `"mz"`: Numeric m/z values.
#'         \item `"intensity"`: Numeric intensity values.
#'         \item `"rt"`: Numeric retention time values (in seconds).
#'         \item `"spectrum_id"`: A unique identifier for each spectrum/scan.
#'           If this column is missing, the function will attempt to use `"rt"`
#'           to group data points, assuming retention times are unique per
#'           spectrum; a warning will be issued in this case.
#'       }
#'       The data.frame should be ordered by retention time, and ideally by m/z
#'       within each spectrum, for optimal performance.
#'     \item A `Spectra` object (from the `Spectra` package). Data will be
#'       extracted using accessor functions.
#'   }
#' @param centwave_params A `list` of parameters to be passed to the
#'   underlying `do_findChromPeaks_centWave` function, which implements the
#'   centWave algorithm. These parameters control aspects like PPM tolerance,
#'   peak width, S/N threshold, etc. For a detailed list of available
#'   parameters and their default values, refer to the documentation of
#'   `xcms::findChromPeaks.centWave`.
#'
#' @return A `data.frame` where each row represents a detected chromatographic
#'   peak. Columns include peak properties such as "mz", "mzmin", "mzmax",
#'   "rt", "rtmin", "rtmax", "into" (integrated intensity), "maxo" (maximum
#'   intensity), "sn" (signal-to-noise ratio), and potentially others if
#'   `verboseColumns = TRUE` is passed in `centwave_params`.
#'
#' @section Dependencies:
#' This script and its functions rely heavily on underlying C functions and R
#' utility functions that are normally part of the `xcms` package. For this
#' script to function correctly, it is essential that the `xcms` package is
#' installed in your R environment.
#'
#' @section TODO:
#' While `continuousPtsAboveThreshold` and `descendMin` are now defined in this
#' script, their implementations are based on common logic or related functions
#' (`continuousPtsAboveThresholdIdx`) and might differ from the exact original
#' internal `xcms` implementations if those were highly specialized. Testing
#' with known datasets is recommended.
#'
#' @examples
#' # 1. Example with a sample data.frame
#' # Create a sample data.frame mimicking MS data
#' set.seed(123)
#' n_spectra <- 10
#' n_peaks_per_spectrum <- 50
#' df_data <- data.frame(
#'   rt = sort(rep(seq(100, 100 + (n_spectra - 1) * 10, by = 10), n_peaks_per_spectrum)),
#'   mz = unlist(lapply(1:n_spectra, function(x) rnorm(n_peaks_per_spectrum, 300 + x*5, 0.01))),
#'   intensity = runif(n_spectra * n_peaks_per_spectrum, 1000, 100000),
#'   spectrum_id = sort(rep(1:n_spectra, n_peaks_per_spectrum))
#' )
#' # Add a more prominent peak for detection
#' peak_rt_center <- 130
#' peak_mz_center <- 315.05
#' peak_indices <- which(df_data$rt > peak_rt_center - 5 & df_data$rt < peak_rt_center + 5 &
#'                       df_data$mz > peak_mz_center - 0.05 & df_data$mz < peak_mz_center + 0.05)
#' df_data$intensity[peak_indices] <- df_data$intensity[peak_indices] *
#'                                    dnorm(df_data$rt[peak_indices], peak_rt_center, 2) *
#'                                    dnorm(df_data$mz[peak_indices], peak_mz_center, 0.005) * 5000
#' df_data$intensity <- pmax(0, df_data$intensity) # Ensure non-negative
#'
#' # Define centWave parameters
#' cwp_params <- list(
#'   ppm = 25,
#'   peakwidth = c(5, 12), # Shorter peakwidth for this example data
#'   snthresh = 5,
#'   prefilter = c(3, 1000), # Adjusted prefilter
#'   verboseColumns = TRUE
#' )
#'
#' # Process the data.frame
#' detected_peaks_df <- process_spectra_with_centwave(df_data, cwp_params)
#' print(head(detected_peaks_df))
#'
#' # 2. Conceptual example with a Spectra object (requires Spectra package)
#' # if (requireNamespace("Spectra", quietly = TRUE) && requireNamespace("MsBackendDataFrame", quietly = TRUE)) {
#' #   library(Spectra)
#' #   # Assuming df_data from above is available
#' #   # Create a Spectra object (requires MsBackendDataFrame or similar)
#' #   # spectra_object <- Spectra(df_data, backend = MsBackendDataFrame())
#' #
#' #   # Define centWave parameters (can be the same or different)
#' #   cwp_params_spectra <- list(ppm = 20, peakwidth = c(10, 60), snthresh = 10)
#' #
#' #   # Process the Spectra object
#' #   # detected_peaks_spectra <- process_spectra_with_centwave(spectra_object, cwp_params_spectra)
#' #   # print(head(detected_peaks_spectra))
#' # } else {
#' #   print("Spectra or MsBackendDataFrame package not available, skipping Spectra object example.")
#' # }
process_spectra_with_centwave <- function(spectra_data, centwave_params = list()) {

    mz_vec <- NULL
    int_vec <- NULL
    rt_vec <- NULL
    vps_vec <- NULL # valsPerSpect

    if (inherits(spectra_data, "data.frame")) {
        # Expected columns for data.frame: "mz", "intensity", "rt", "spectrum_id"
        # "spectrum_id" is used to group data points belonging to the same scan.
        # If "spectrum_id" is not present, "rt" will be used with a warning,
        # assuming retention times are unique identifiers for spectra.
        # Data should ideally be ordered by spectrum_id (or rt if used as id)
        # and then by mz within each spectrum.

        required_cols <- c("mz", "intensity", "rt")
        if (!all(required_cols %in% colnames(spectra_data))) {
            missing_cols <- required_cols[!required_cols %in% colnames(spectra_data)]
            stop("The input data.frame is missing the following required columns: ",
                 paste(missing_cols, collapse = ", "))
        }

        # Order by retention time first to ensure correct processing order
        spectra_data <- spectra_data[order(spectra_data$rt), ]

        if ("spectrum_id" %in% colnames(spectra_data)) {
            spectrum_ids <- spectra_data$spectrum_id
        } else {
            warning("Column 'spectrum_id' not found in data.frame. Using 'rt' to ",
                    "group data points per spectrum. Ensure 'rt' is unique per spectrum.")
            spectrum_ids <- spectra_data$rt
        }

        # Ensure rt_vec contains unique retention times, one per spectrum
        # The values in vps_vec are the counts of mz/intensity pairs per spectrum
        # The scantime (rt_vec here) must have the same length as vps_vec
        
        # Get unique retention times corresponding to each spectrum
        unique_rts_table <- table(spectrum_ids)
        rt_vec <- as.numeric(names(unique_rts_table)) # These are the unique spectrum identifiers (rt or spectrum_id)
        
        # If spectrum_id was not present, rt_vec is already correct (unique rts).
        # If spectrum_id was present, we need to get the actual rt values for these unique spectra.
        # This assumes that within each spectrum_id group, all rt values are the same.
        # We take the first rt for each unique spectrum_id.
        if ("spectrum_id" %in% colnames(spectra_data)) {
             rt_map <- unique(spectra_data[, c("spectrum_id", "rt")])
             rt_lookup <- rt_map$rt
             names(rt_lookup) <- rt_map$spectrum_id
             rt_vec <- rt_lookup[names(unique_rts_table)]
        }


        # Order rt_vec and then get vps_vec in that order
        rt_order <- order(rt_vec)
        rt_vec <- rt_vec[rt_order]
        
        # Calculate valsPerSpect based on the (ordered) spectrum_ids
        # The table function counts occurrences of each unique spectrum_id
        vps_table <- table(factor(spectrum_ids, levels = names(unique_rts_table)[rt_order]))
        vps_vec <- as.integer(vps_table)


        # Prepare mz_vec and int_vec, ensuring they match the order of spectra
        # defined by rt_vec and vps_vec
        # This requires splitting the data by spectrum_id, then ordering each
        # spectrum's data by mz (though .centWave_orig has a fallback),
        # and then unlisting.
        
        # Ensure spectra_data is ordered by the (now ordered) rt_vec and then mz
        # This is a bit complex if spectrum_id was used.
        # A simpler way: iterate through the ordered unique spectrum_ids (derived from rt_vec's names)
        
        mz_list <- list()
        int_list <- list()
        
        original_spectrum_ids_ordered <- names(vps_table) # These are ordered by rt_vec

        for(id in original_spectrum_ids_ordered) {
            current_spectrum_data <- spectra_data[spectrum_ids == id, ]
            # Optional: order by mz within each spectrum, though .centWave_orig handles unsorted mz.
            # current_spectrum_data <- current_spectrum_data[order(current_spectrum_data$mz), ]
            mz_list[[as.character(id)]] <- current_spectrum_data$mz
            int_list[[as.character(id)]] <- current_spectrum_data$intensity
        }
        
        mz_vec <- unlist(mz_list, use.names = FALSE)
        int_vec <- unlist(int_list, use.names = FALSE)


    } else if (inherits(spectra_data, "Spectra")) {
        if (!requireNamespace("Spectra", quietly = TRUE)) {
            stop("The 'Spectra' package is required to process Spectra objects. ",
                 "Please install it using BiocManager::install('Spectra').")
        }
        # Extract data using Spectra accessors
        # Ensure data is ordered by acquisitionNum or rtime by default by Spectra
        rt_vec <- Spectra::rtime(spectra_data)
        
        # Sort by retention time if not already sorted (Spectra objects usually are)
        spec_order <- order(rt_vec)
        spectra_data_ordered <- spectra_data[spec_order]
        
        rt_vec <- Spectra::rtime(spectra_data_ordered) # Re-extract ordered rtime

        # Get mz and intensity values as lists, then unlist
        mz_list <- Spectra::mz(spectra_data_ordered)
        int_list <- Spectra::intensity(spectra_data_ordered)

        mz_vec <- unlist(mz_list, use.names = FALSE)
        int_vec <- unlist(int_list, use.names = FALSE)
        
        # Calculate valsPerSpect (number of peaks per spectrum)
        vps_vec <- Spectra::peaksCount(spectra_data_ordered)
        # Alternative if peaksCount is not available or desired:
        # vps_vec <- as.integer(lengths(mz_list))

    } else {
        stop("Unsupported input type for 'spectra_data'. Please provide a ",
             "data.frame or a Spectra object.")
    }

    # Ensure all vectors are not NULL before calling centWave
    if (is.null(mz_vec) || is.null(int_vec) || is.null(rt_vec) || is.null(vps_vec)) {
        stop("Failed to prepare necessary data vectors from the input.")
    }
    
    if(length(rt_vec) != length(vps_vec)){
        stop("Internal error: Length of scantime vector (rt_vec) must match length of values-per-spectrum vector (vps_vec).")
    }
    if(length(mz_vec) != sum(vps_vec) || length(int_vec) != sum(vps_vec)){
        stop("Internal error: Total number of m/z or intensity values does not match sum of values-per-spectrum.")
    }


    # Call do_findChromPeaks_centWave
    # Note: Relies on xcms being installed for C functions.
    # Note: continuousPtsAboveThreshold and descendMin are TODOs if not defined.
    detected_peaks_matrix <- do.call(do_findChromPeaks_centWave, c(
        list(
            mz = mz_vec,
            int = int_vec,
            scantime = rt_vec,
            valsPerSpect = vps_vec
        ),
        centwave_params
    ))

    result_df <- as.data.frame(detected_peaks_matrix)
    # Column names are preserved by as.data.frame for matrices.

    return(result_df)
}
