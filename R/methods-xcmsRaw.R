## All Methods for xcmsRaw should be here.
#' @include functions-xcmsRaw.R functions-utils.R

############################################################
## show
setMethod("show", "xcmsRaw", function(object) {

    cat("An \"xcmsRaw\" object with", length(object@scantime),
        "mass spectra\n\n")
    if (length(object@scantime)>0) {
        cat("Time range: ", paste(round(range(object@scantime), 1),
                                  collapse = "-"),
            " seconds (", paste(round(range(object@scantime)/60, 1),
                                collapse = "-"),
            " minutes)\n", sep = "")
        cat("Mass range:", paste(round(range(object@env$mz), 4),
                                 collapse = "-"), "m/z\n")
        cat("Intensity range:", paste(signif(range(object@env$intensity), 6),
                                      collapse = "-"), "\n\n")
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
            cat(names(object@profparam)[i], " = ", object@profparam[[i]], "\n",
                sep = "")
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
## findPeaks.matchedFilter_orig
## We're keeping this one only for comparison in unit tests; this
## should be removed soon!
setGeneric("findPeaks.matchedFilter_orig", function(object, ...)
    standardGeneric("findPeaks.matchedFilter_orig"))
setMethod("findPeaks.matchedFilter_orig", "xcmsRaw",
          function(object, fwhm = 30, sigma = fwhm/2.3548, max = 5,
                   snthresh = 10, step = 0.1, steps = 2,
                   mzdiff = 0.8 - step*steps, index = FALSE, sleep = 0,
                   scanrange= numeric()) {

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
## findPeaks.matchedFilter
#' @title Peak detection in the chromatographic time domain
#'
#' @aliases findPeaks.matchedFilter
#' 
#' @description Find peaks in the chromatographic time domain of the
#'     profile matrix. For more details see
#'     \code{\link{do_findChromPeaks_matchedFilter}}.
#' 
#' @param object The \code{\linkS4class{xcmsRaw}} object on which peak detection
#'     should be performed.
#' 
#' @inheritParams findChromPeaks-matchedFilter
#' 
#' @param step numeric(1) specifying the width of the bins/slices in m/z
#'     dimension.
#' 
#' @param sleep (DEPRECATED). The use of this parameter is highly discouraged,
#'     as it could cause problems in parallel processing mode.
#' 
#' @param scanrange Numeric vector defining the range of scans to which the
#'     original \code{object} should be sub-setted before peak detection.
#' 
#' @author Colin A. Smith
#' 
#' @return A matrix, each row representing an intentified chromatographic peak,
#'     with columns:
#'     \describe{
#'     \item{mz}{Intensity weighted mean of m/z values of the peak across
#'     scans.}
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
#'     \item{sn}{Signal to noise ratio of the peak.}
#'     }
#' 
#' @references
#'     Colin A. Smith, Elizabeth J. Want, Grace O'Maille, Ruben Abagyan and
#'     Gary Siuzdak. "XCMS: Processing Mass Spectrometry Data for Metabolite
#'     Profiling Using Nonlinear Peak Alignment, Matching, and Identification"
#'     \emph{Anal. Chem.} 2006, 78:779-787.
#'     @family Old peak detection methods
#' 
#' @seealso \code{\link{matchedFilter}} for the new user interface.
#'     \code{\linkS4class{xcmsRaw}},
#'     \code{\link{do_findChromPeaks_matchedFilter}} for the core function
#'     performing the peak detection.
setMethod("findPeaks.matchedFilter", "xcmsRaw",
          function(object, fwhm = 30, sigma = fwhm/2.3548, max = 5,
                   snthresh = 10, step = 0.1, steps = 2,
                   mzdiff = 0.8 - step*steps, index = FALSE, sleep = 0,
                   scanrange = numeric()) {

              ## Fix issue #63:
              ## Sub-set the xcmsRaw based on scanrange
              if (length(scanrange) < 2) {
                  scanrange <- c(1, length(object@scantime))
              } else {
                  scanrange <- range(scanrange)
              }
              if (min(scanrange) < 1 | max(scanrange) > length(object@scantime)) {
                  scanrange[1] <- max(1, scanrange[1])
                  scanrange[2] <- min(length(object@scantime), scanrange[2])
                  message("Provided scanrange was adjusted to ", scanrange)
              }
              object <- object[scanrange[1]:scanrange[2]]
              scanrange <- c(1, length(object@scantime))
              ## ## Sub-set the xcmsRaw baesd on scanrange
              ## scanrange.old <- scanrange
              ## ## sanitize if too few or too many scanrange is given
              ## if (length(scanrange) < 2)
              ##     scanrange <- c(1, length(object@scantime))
              ## else
              ##     scanrange <- range(scanrange)
              ## ## restrict and sanitize scanrange to maximally cover all scans
              ## scanrange[1] <- max(1,scanrange[1])
              ## scanrange[2] <- min(length(object@scantime),scanrange[2])
              ## ## Mild warning if the actual scanrange doesn't match the scanrange
              ## ## argument
              ## if (!(identical(scanrange.old, scanrange)) &&
              ##     (length(scanrange.old) > 0)) {
              ##     cat("Warning: scanrange was adjusted to ",scanrange,"\n")
              ##     ## Scanrange filtering
              ##     keepidx <- seq.int(1, length(object@scantime))
              ##     %in% seq.int(scanrange[1], scanrange[2])
              ##     object <- split(object, f=keepidx)[["TRUE"]]
              ## }
              ## Determine the impute method:
              imputeMeths <- c("none", "lin", "linbase", "intlin")
              names(imputeMeths) <- c("bin", "binlin",
                                      "binlinbase", "intlin")
              profFun <- profMethod(object)
              profFun <- match.arg(profFun, names(imputeMeths))
              imputeMeth <- imputeMeths[profFun]
              if (imputeMeth == "linbase") {
                  profp <- object@profparam
                  if (length(profp) == 0)
                      profp <- list()
                  ## Determine the settings for this:
                  ## o distance
                  ##   Define the distance argument; that's tricky, as it
                  ##   requires the bin_size, not the step.
                  mrange <- range(object@env$mz)
                  mass <- seq(floor(mrange[1]/step)*step,
                              ceiling(mrange[2]/step)*step, by = step)
                  mlength <- length(mass)
                  bin_size <- (mass[mlength] - mass[1]) / (mlength - 1)
                  rm(mass)
                  if (length(profp$basespace) > 0) {
                      if (!is.numeric(profp$basespace))
                          stop("Profile parameter 'basespace' has to be numeric!")
                      distance <- floor(profp$basespace[1] / bin_size)
                  } else {
                      distance <- floor(0.075 / bin_size)
                  }
                  ## o baseValue
                  if (length(profp$baseleve) > 0) {
                      if (!is.numeric(profp$baselevel))
                          stop("Profile parameter 'baselevel' has to be numeric!")
                      baseValue <- profp$baselevel[1]
                  } else {
                      baseValue <- min(object@env$intensity) / 2
                  }
              } else {
                  ## For other methods these are not used anyway.
                  distance <- 0
                  baseValue <- 0
              }
              res <- do_findChromPeaks_matchedFilter(mz = object@env$mz,
                                                     int = object@env$intensity,
                                                     scantime = object@scantime,
                                                     valsPerSpect = diff(c(object@scanindex,
                                                                           length(object@env$mz))),
                                                     binSize = step,
                                                     impute = imputeMeth,
                                                     baseValue = baseValue,
                                                     distance = distance,
                                                     fwhm = fwhm,
                                                     sigma = sigma,
                                                     max = max,
                                                     snthresh = snthresh,
                                                     steps = steps,
                                                     mzdiff = mzdiff,
                                                     index = index,
                                                     sleep = sleep
                                                     )
              invisible(new("xcmsPeaks", res))
})


############################################################
## findPeaks.centWave
setMethod("findPeaks.centWave", "xcmsRaw", function(object, ppm=25,
                                                    peakwidth=c(20,50),
                                                    snthresh=10,
                                                    prefilter=c(3,100),
                                                    mzCenterFun="wMean",
                                                    integrate=1, mzdiff=-0.001,
                                                    fitgauss=FALSE,
                                                    scanrange = numeric(),
                                                    noise=0, ## noise.local=TRUE,
                                                    sleep=0,
                                                    verbose.columns=FALSE,
                                                    ROI.list=list(),
                                                    firstBaselineCheck=TRUE,
                                                    roiScales=NULL) {
    if (!isCentroided(object))
        warning("It looks like this file is in profile mode. centWave can",
                " process only centroid mode data !\n")

    ## Fix issue #64:
    ## Sub-set the xcmsRaw based on scanrange
    if (length(scanrange) < 2) {
        scanrange <- c(1, length(object@scantime))
    } else {
        scanrange <- range(scanrange)
    }
    if (min(scanrange) < 1 | max(scanrange) > length(object@scantime)) {
        scanrange[1] <- max(1, scanrange[1])
        scanrange[2] <- min(length(object@scantime), scanrange[2])
        message("Provided scanrange was adjusted to ", scanrange)
    }
    object <- object[scanrange[1]:scanrange[2]]

    vps <- diff(c(object@scanindex, length(object@env$mz)))
    res <- do_findChromPeaks_centWave(mz = object@env$mz,
                                      int = object@env$intensity,
                                      scantime = object@scantime,
                                      valsPerSpect = vps,
                                      ppm = ppm, peakwidth = peakwidth,
                                      snthresh = snthresh,
                                      prefilter = prefilter,
                                      mzCenterFun = mzCenterFun,
                                      integrate = integrate,
                                      mzdiff = mzdiff, fitgauss = fitgauss,
                                      noise = noise,
                                      verboseColumns = verbose.columns,
                                      roiList = ROI.list,
                                      firstBaselineCheck = firstBaselineCheck,
                                      roiScales = roiScales,
                                      sleep = sleep
                                      )
    invisible(new("xcmsPeaks", res))
})

############################################################
## findPeaks.centWaveWithPredictedIsotopeROIs
## Performs first a centWave analysis and based on the identified peaks
## defines ROIs for a second centWave run to check for presence of
## predicted isotopes for the first peaks.
setMethod("findPeaks.centWaveWithPredictedIsotopeROIs", "xcmsRaw",
          function(object, ppm = 25, peakwidth = c(20,50), snthresh = 10,
                   prefilter = c(3,100), mzCenterFun = "wMean", integrate = 1,
                   mzdiff = -0.001, fitgauss = FALSE, scanrange = numeric(),
                   noise = 0, sleep = 0, verbose.columns = FALSE,
                   ROI.list = list(), firstBaselineCheck = TRUE,
                   roiScales = NULL, snthreshIsoROIs = 6.25, maxcharge = 3,
                   maxiso = 5, mzIntervalExtension = TRUE) {
              if (!isCentroided(object))
                  warning("It looks like this file is in profile mode. centWave",
                          " can process only centroid mode data !\n")

              ## Sub-set the xcmsRaw based on scanrange
              if (length(scanrange) < 2) {
                  scanrange <- c(1, length(object@scantime))
              } else {
                  scanrange <- range(scanrange)
              }
              if (min(scanrange) < 1 | max(scanrange) > length(object@scantime)) {
                  scanrange[1] <- max(1, scanrange[1])
                  scanrange[2] <- min(length(object@scantime), scanrange[2])
                  message("Provided scanrange was adjusted to ", scanrange)
              }
              object <- object[scanrange[1]:scanrange[2]]

              vps <- diff(c(object@scanindex, length(object@env$mz)))
              res <- do_findChromPeaks_centWaveWithPredIsoROIs(mz = object@env$mz,
                                                               int = object@env$intensity,
                                                               scantime = object@scantime,
                                                               valsPerSpect = vps,
                                                               ppm = ppm,
                                                               peakwidth = peakwidth,
                                                               snthresh = snthresh,
                                                               prefilter = prefilter,
                                                               mzCenterFun = mzCenterFun,
                                                               integrate = integrate,
                                                               mzdiff = mzdiff,
                                                               fitgauss = fitgauss,
                                                               noise = noise,
                                                               verboseColumns = verbose.columns,
                                                               roiList = ROI.list,
                                                               firstBaselineCheck = firstBaselineCheck,
                                                               roiScales = roiScales,
                                                               snthreshIsoROIs = snthreshIsoROIs,
                                                      maxCharge = maxcharge,
                                                      maxIso = maxiso,
                                                      mzIntervalExtension = mzIntervalExtension
                                                      )
              invisible(new("xcmsPeaks", res))
          })

setMethod("findPeaks.addPredictedIsotopeFeatures",
          "xcmsRaw", function(object, ppm = 25, peakwidth = c(20,50),
                              prefilter = c(3,100), mzCenterFun = "wMean",
                              integrate = 1, mzdiff = -0.001, fitgauss = FALSE,
                              scanrange = numeric(), noise=0, ## noise.local=TRUE,
                              sleep = 0, verbose.columns = FALSE,
                              xcmsPeaks, snthresh = 6.25, maxcharge = 3,
                              maxiso = 5, mzIntervalExtension = TRUE) {
              if (!isCentroided(object))
                  warning("It looks like this file is in profile mode. centWave",
                          " can process only centroid mode data !\n")

              ## Sub-set the xcmsRaw based on scanrange
              if (length(scanrange) < 2) {
                  scanrange <- c(1, length(object@scantime))
              } else {
                  scanrange <- range(scanrange)
              }
              if (min(scanrange) < 1 | max(scanrange) > length(object@scantime)) {
                  scanrange[1] <- max(1, scanrange[1])
                  scanrange[2] <- min(length(object@scantime), scanrange[2])
                  message("Provided scanrange was adjusted to ", scanrange)
              }
              object <- object[scanrange[1]:scanrange[2]]
              if(class(xcmsPeaks) != "xcmsPeaks")
                  stop("Pparameter >xcmsPeaks< is not of class 'xcmsPeaks'!\n")

              vps <- diff(c(object@scanindex, length(object@env$mz)))
              res <- do_findChromPeaks_addPredIsoROIs(mz = object@env$mz,
                                                      int = object@env$intensity,
                                                      scantime = object@scantime,
                                                      valsPerSpect = vps,
                                                      ppm = ppm,
                                                      peakwidth = peakwidth,
                                                      snthresh = snthresh,
                                                      prefilter = prefilter,
                                                      mzCenterFun = mzCenterFun,
                                                      integrate = integrate,
                                                      mzdiff = mzdiff,
                                                      fitgauss = fitgauss,
                                                      noise = noise,
                                                      verboseColumns = verbose.columns,
                                                      peaks. = xcmsPeaks@.Data,
                                                      maxCharge = maxcharge,
                                                      maxIso = maxiso,
                                                      mzIntervalExtension = mzIntervalExtension
                                                      )
              invisible(new("xcmsPeaks", res))
          })


############################################################
## findPeaks.MSW
#' @title Peak detection for single-spectrum non-chromatography MS data
#' 
#' @aliases findPeaks.MSW
#'
#' @description This method performs peak detection in mass spectrometry
#'     direct injection spectrum using a wavelet based algorithm.
#'
#' @details This is a wrapper around the peak picker in Bioconductor's
#'     \code{MassSpecWavelet} package calling
#'     \code{\link{peakDetectionCWT}} and
#'     \code{\link{tuneInPeakInfo}} functions.
#'
#' @inheritParams findPeaks-MSW
#' 
#' @inheritParams findChromPeaks-centWave
#' 
#' @param object The \code{\linkS4class{xcmsRaw}} object on which peak
#'     detection should be performed.
#' 
#' @param verbose.columns Logical whether additional peak meta data columns
#'     should be returned.
#'
#' @return
#'     A matrix, each row representing an intentified peak, with columns:
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
#' @seealso \code{\link{MSW}} for the new user interface,
#'     \code{\link{do_findPeaks_MSW}} for the downstream analysis
#'     function or \code{\link{peakDetectionCWT}} from the
#'     \code{MassSpecWavelet} for details on the algorithm and additionally
#'     supported parameters.
#'
#' @author Joachim Kutzera, Steffen Neumann, Johannes Rainer
setMethod("findPeaks.MSW", "xcmsRaw",
          function(object, snthresh=3, verbose.columns = FALSE, ...) {
              if (length(object@scantime) > 1)
                  stop("MSW works only on single spectrum, direct injection",
                       " MS data, but 'object' has ", length(object@scantime),
                       " spectra")
              res <- do_findPeaks_MSW(mz = object@env$mz,
                                      int = object@env$intensity,
                                      snthresh = snthresh,
                                      verboseColumns = verbose.columns,
                                      ...)
              invisible(new("xcmsPeaks", res))
          })

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

setMethod("getPeaks", "xcmsRaw", function(object, peakrange, step = 0.1) {
    if (useOriginalCode())
        return(.getPeaks_orig(object, peakrange, step = step))
    else
        return(.getPeaks_new(object, peakrange, step = step))
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
## Issue #74: implement an alternative (improved) getEIC method.
setMethod("getEIC", "xcmsRaw", function(object, mzrange, rtrange = NULL,
                                        step = 0.1) {
    profEIC(object, mzrange = mzrange, rtrange = rtrange, step = step)
})

############################################################
## rawMat
#' @description Extracts a matrix with columns time (retention time), mz and
#' intensity from an xcmsRaw object.
#'
#' @noRd
setMethod("rawMat", "xcmsRaw", function(object,
                                        mzrange = numeric(),
                                        rtrange = numeric(),
                                        scanrange = numeric(),
                                        log=FALSE) {
    .rawMat(mz = object@env$mz, int = object@env$intensity,
            scantime = object@scantime,
            valsPerSpect = diff(c(object@scanindex, length(object@env$mz))),
            mzrange = mzrange, rtrange = rtrange, scanrange = scanrange,
            log = log)
})
## @jo TODO LLL replace that with an implementation in C.
## Note: this function silently drops retention times for which no intensity-mz
## pair was measured.
.rawMat <- function(mz, int, scantime, valsPerSpect, mzrange = numeric(),
                    rtrange = numeric(), scanrange = numeric,
                    log = FALSE) {
    if (length(rtrange) >= 2) {
        rtrange <- range(rtrange)
        scanrange <- range(which((scantime >= rtrange[1]) &
                                 (scantime <= rtrange[2])))
    }
    if (length(scanrange) < 2)
        scanrange <- c(1, length(valsPerSpect))
    else scanrange <- range(scanrange)
    if (scanrange[1] == 1)
        startidx <- 1
    else
        startidx <- sum(valsPerSpect[1:(scanrange[1]-1)]) + 1
    endidx <- sum(valsPerSpect[1:scanrange[2]])
    scans <- rep(scanrange[1]:scanrange[2],
                 valsPerSpect[scanrange[1]:scanrange[2]])
    masses <- mz[startidx:endidx]
    massidx <- 1:length(masses)
    if (length(mzrange) >= 2) {
        mzrange <- range(mzrange)
        massidx <- massidx[(masses >= mzrange[1] & (masses <= mzrange[2]))]
    }
    int <- int[startidx:endidx][massidx]
    if (log && (length(int) > 0))
        int <- log(int + max(1 - min(int), 0))
    cbind(time = scantime[scans[massidx]],
          mz = masses[massidx],
          intensity = int)
}

## .rawMat2 <- function(mz, int, scantime, valsPerSpect, mzrange = numeric(),
##                      rtrange = numeric(), scanrange = numeric,
##                      log = FALSE) {
##     if (length(rtrange) >= 2) {
##         rtrange <- range(rtrange)
##         scanrange <- range(which((scantime >= rtrange[1]) &
##                                  (scantime <= rtrange[2])))
##     }
##     if (length(scanrange) < 2)
##         scanrange <- c(1, length(valsPerSpect))
##     else scanrange <- range(scanrange)
##     if (scanrange[1] == 1)
##         startidx <- 1
##     else
##         startidx <- sum(valsPerSpect[1:(scanrange[1]-1)]) + 1
##     endidx <- sum(valsPerSpect[1:scanrange[2]])
##     scans <- rep(scanrange[1]:scanrange[2],
##                  valsPerSpect[scanrange[1]:scanrange[2]])
##     masses <- mz[startidx:endidx]
##     massidx <- 1:length(masses)
##     if (length(mzrange) >= 2) {
##         mzrange <- range(mzrange)
##         massidx <- massidx[(masses >= mzrange[1] & (masses <= mzrange[2]))]
##     }
##     int <- int[startidx:endidx][massidx]
##     if (log && (length(int) > 0))
##         int <- log(int + max(1 - min(int), 0))
##     cbind(time = scantime[scans[massidx]],
##           mz = masses[massidx],
##           intensity = int)
## }


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
    if (length(object@env$mz) == 0) {
        warning("MS1 scans empty. Skipping profile matrix creation.")
        return(object)
    }
    have_profmethod <- object@profmethod
    ## Re-calculate the profile matrix if method differs
    if (have_profmethod != value & profStep(object) > 0) {
        object@env$profile <- profMat(object, method = value)
    }
    object@profmethod <- value
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
## Update: related to issue #71
setReplaceMethod("profStep", "xcmsRaw", function(object, value) {
    if (!is.numeric(value) && value < 0)
        stop("'value' has to be a positive number!")
    if (length(object@env$mz) == 0) {
        warning("MS1 scans empty. Skipping profile matrix creation.")
        return(object)
    }
    if (value == 0) {
        if ("profile" %in% ls(object@env)) {
            rm("profile", envir = object@env)
            message("Removing profile matrix.")
        }
        return(object)
    }
    have_step <- profStep(object)
    ## Check if the value differs from step and only calculate if different.
    if (have_step != value) {
        ## OK, now we're re-calculating
        minmass <- round(min(object@env$mz) / value) * value
        maxmass <- round(max(object@env$mz) / value) * value
        object@mzrange <- c(minmass, maxmass)
        ## Fix for issue #98: to be in accordance with the "old" code we require
        ## that the number of rows of the profile matrix matches minmass to
        ## maxmass in steps of value
        ## To me that is somewhat problematic, as it means that the @mzrange does
        ## not correctly correspond to the range(object@env$mz)!
        tmp <- seq(minmass, maxmass, by = value)
        prf <- profMat(object, step = value)
        object@env$profile <- prf[1:min(c(length(tmp), nrow(prf))), ]
    }
    return(object)
})

############################################################
## profStepPad
## The difference to the profStep? seems only to be the way how the range
## is calculated:
## profStep: minmass <- round(min(object@env$mz)/value)*value
## profStepPad: floor(mzrange(range(object@env$mz))[1])
setReplaceMethod("profStepPad", "xcmsRaw", function(object, value) {
    if (!is.numeric(value) && value < 0)
        stop("'value' has to be a positive number!")
    if (length(object@env$mz) == 0) {
        warning("MS1 scans empty. Skipping profile matrix creation.")
        return(object)
    }
    if (value == 0) {
        if ("profile" %in% ls(object@env)) {
            rm("profile", envir = object@env)
            message("Removing profile matrix.")
        }
        return(object)
    }
    mzr <- range(object@env$mz)
    minmass <- floor(mzr[1])
    maxmass <- ceiling(mzr[2])
    object@mzrange <- c(minmass, maxmass)
    object@env$profile <- profMat(object, step = value,
                                  mzrange. = c(minmass, maxmass))
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
        scanrange <- c(match(TRUE, scanidx), length(scanidx) -
                                             match(TRUE, rev(scanidx)))
    }  else if (length(scanrange) < 2)
        scanrange <- c(1, length(object@scantime))
    else
        scanrange <- range(scanrange)

    scanrange[1] <- max(1,scanrange[1])
    scanrange[2] <- min(length(object@scantime), scanrange[2])

    if (!is.double(object@env$mz))
        object@env$mz <- as.double(object@env$mz)
    if (!is.double(object@env$intensity))
        object@env$intensity <- as.double(object@env$intensity)
    if (!is.integer(object@scanindex))
        object@scanindex <- as.integer(object@scanindex)

    .Call("getEIC", object@env$mz, object@env$intensity, object@scanindex,
          as.double(mzrange), as.integer(scanrange),
          as.integer(length(object@scantime)), PACKAGE ='xcms' )
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
                                               scanrange=c(1,length(object@scantime)),
                                               minIntensity, minCentroids,
                                               consecMissedLim, criticalVal,
                                               ppm,  segs, scanBack){

    scanrange[1] <- max(1, scanrange[1])
    scanrange[2] <- min(length(object@scantime), scanrange[2])

    ## Subset object by scanrange.

    valsPS <- diff(c(object@scanindex, length(object@env$mz)))
    do_findKalmanROI(mz = object@env$mz, int = object@env$intensity,
                     scantime = object@scantime,
                     valsPerSpect = valsPS, mzrange = mzrange,
                     scanrange = scanrange, minIntensity = minIntensity,
                     minCentroids = minCentroids,
                     consecMissedLim = consecMissedLim,
                     criticalVal = criticalVal, ppm = ppm, segs = segs,
                     scanBack = scanBack)

    ## .Call("massifquant", object@env$mz, object@env$intensity, object@scanindex,
    ##       object@scantime, as.double(mzrange), as.integer(scanrange),
    ##       as.integer(length(object@scantime)), as.double(minIntensity),
    ##       as.integer(minCentroids),as.double(consecMissedLim), as.double(ppm),
    ##       as.double(criticalVal), as.integer(segs), as.integer(scanBack),
    ##       PACKAGE ='xcms' )
})


############################################################
## findPeaks.massifquant: Note the original code returned, if withWave = 1,
## a Peaks object otherwise a matrix!
setMethod("findPeaks.massifquant", "xcmsRaw", function(object,
                                                       ppm=10,
                                                       peakwidth = c(20,50),
                                                       snthresh = 10,
                                                       prefilter = c(3,100),
                                                       mzCenterFun = "wMean",
                                                       integrate = 1,
                                                       mzdiff = -0.001,
                                                       fitgauss = FALSE,
                                                       scanrange = numeric(),
                                                       noise = 0,
                                                       sleep = 0,
                                                       verbose.columns = FALSE,
                                                       criticalValue = 1.125,
                                                       consecMissedLimit = 2,
                                                       unions = 1,
                                                       checkBack = 0,
                                                       withWave = 0) {

    if (sleep > 0)
        cat("'sleep' argument is defunct and will be ignored.")

    ## Fix issue #61
    ## Sub-set the xcmsRaw based on scanrange
    if (length(scanrange) < 2) {
        scanrange <- c(1, length(object@scantime))
    } else {
        scanrange <- range(scanrange)
    }
    if (min(scanrange) < 1 | max(scanrange) > length(object@scantime)) {
        scanrange[1] <- max(1, scanrange[1])
        scanrange[2] <- min(length(object@scantime), scanrange[2])
        message("Provided scanrange was adjusted to ", scanrange)
    }
    object <- object[scanrange[1]:scanrange[2]]
    scanrange <- c(1, length(object@scantime))
    ## scanrange.old <- scanrange
    ## ## sanitize if too few or too many scanrange is given
    ## if (length(scanrange) < 2)
    ##     scanrange <- c(1, length(object@scantime))
    ## else
    ##     scanrange <- range(scanrange)
    ## ## restrict and sanitize scanrange to maximally cover all scans
    ## scanrange[1] <- max(1,scanrange[1])
    ## scanrange[2] <- min(length(object@scantime),scanrange[2])
    ## ## Mild warning if the actual scanrange doesn't match the scanrange
    ## ## argument
    ## if (!(identical(scanrange.old, scanrange)) &&
    ##     (length(scanrange.old) > 0)) {
    ##     cat("Warning: scanrange was adjusted to ",scanrange,"\n")
    ##     ## Scanrange filtering
    ##     keepidx <- seq.int(1, length(object@scantime)) %in% seq.int(scanrange[1], scanrange[2])
    ##     object <- split(object, f=keepidx)[["TRUE"]]
    ## }
    if (!isCentroided(object))
        warning("It looks like this file is in profile mode.",
                " Massifquant can process only centroid mode data !\n")
    vps <- diff(c(object@scanindex, length(object@env$mz)))
    res <- do_findChromPeaks_massifquant(mz = object@env$mz,
                                         int = object@env$intensity,
                                         scantime = object@scantime,
                                         valsPerSpect = vps,
                                         ppm = ppm, peakwidth = peakwidth,
                                         snthresh = snthresh,
                                         prefilter = prefilter,
                                         mzCenterFun = mzCenterFun,
                                         integrate = integrate,
                                         mzdiff = mzdiff, fitgauss = fitgauss,
                                         noise = noise,
                                         verboseColumns = verbose.columns,
                                         criticalValue = criticalValue,
                                         consecMissedLimit = consecMissedLimit,
                                         unions = unions, checkBack = checkBack,
                                         withWave = as.logical(withWave))
    invisible(new("xcmsPeaks", res))
})

############################################################
## findPeaks.massifquant: Note the original code returned, if withWave = 1,
## a Peaks object otherwise a matrix!
setMethod("findPeaks.massifquant", "xcmsRaw", function(object,
                                                       ppm=10,
                                                       peakwidth = c(20,50),
                                                       snthresh = 10,
                                                       prefilter = c(3,100),
                                                       mzCenterFun = "wMean",
                                                       integrate = 1,
                                                       mzdiff = -0.001,
                                                       fitgauss = FALSE,
                                                       scanrange = numeric(),
                                                       noise = 0,
                                                       sleep = 0,
                                                       verbose.columns = FALSE,
                                                       criticalValue = 1.125,
                                                       consecMissedLimit = 2,
                                                       unions = 1,
                                                       checkBack = 0,
                                                       withWave = 0) {

    if (sleep > 0)
        cat("'sleep' argument is defunct and will be ignored.")

    ## Fix issue #61
    ## Sub-set the xcmsRaw based on scanrange
    if (length(scanrange) < 2) {
        scanrange <- c(1, length(object@scantime))
    } else {
        scanrange <- range(scanrange)
    }
    if (min(scanrange) < 1 | max(scanrange) > length(object@scantime)) {
        scanrange[1] <- max(1, scanrange[1])
        scanrange[2] <- min(length(object@scantime), scanrange[2])
        message("Provided scanrange was adjusted to ", scanrange)
    }
    object <- object[scanrange[1]:scanrange[2]]
    scanrange <- c(1, length(object@scantime))
    ## scanrange.old <- scanrange
    ## ## sanitize if too few or too many scanrange is given
    ## if (length(scanrange) < 2)
    ##     scanrange <- c(1, length(object@scantime))
    ## else
    ##     scanrange <- range(scanrange)
    ## ## restrict and sanitize scanrange to maximally cover all scans
    ## scanrange[1] <- max(1,scanrange[1])
    ## scanrange[2] <- min(length(object@scantime),scanrange[2])
    ## ## Mild warning if the actual scanrange doesn't match the scanrange
    ## ## argument
    ## if (!(identical(scanrange.old, scanrange)) &&
    ##     (length(scanrange.old) > 0)) {
    ##     cat("Warning: scanrange was adjusted to ",scanrange,"\n")
    ##     ## Scanrange filtering
    ##     keepidx <- seq.int(1, length(object@scantime)) %in% seq.int(scanrange[1], scanrange[2])
    ##     object <- split(object, f=keepidx)[["TRUE"]]
    ## }
    if (!isCentroided(object))
        warning("It looks like this file is in profile mode.",
                " Massifquant can process only centroid mode data !\n")
    vps <- diff(c(object@scanindex, length(object@env$mz)))
    res <- do_findChromPeaks_massifquant(mz = object@env$mz,
                                         int = object@env$intensity,
                                         scantime = object@scantime,
                                         valsPerSpect = vps,
                                         ppm = ppm, peakwidth = peakwidth,
                                         snthresh = snthresh,
                                         prefilter = prefilter,
                                         mzCenterFun = mzCenterFun,
                                         integrate = integrate,
                                         mzdiff = mzdiff, fitgauss = fitgauss,
                                         noise = noise,
                                         verboseColumns = verbose.columns,
                                         criticalValue = criticalValue,
                                         consecMissedLimit = consecMissedLimit,
                                         unions = unions, checkBack = checkBack,
                                         withWave = as.logical(withWave))
    invisible(new("xcmsPeaks", res))
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

    ## Array [x, y, z] with
    ## - x: mz and intensity
    ## - y: spectrum (1: max measurements within one of the spectra)
    ## - z: scans (1: number of spectra)
    arr <- array(dim = c(2, max(diff(ob@scanindex)), length(ob@scanindex)))
    if(lockMass[1] == 1){
        lockMass<-lockMass[3:length(lockMass)]
    }

    ## Remove the last lock mass if it is too close by the end
    if ((lockMass[length(lockMass)] + 2) > length(ob@scanindex))
        lockMass <- lockMass[1:(length(lockMass) - 1)]
    
    ## If the number of lockMass values is not even splitting them into a
    ## two-column matrix is not OK (causes also the first lockMass spectrum to
    ## be overwritten twice. That's to get rid of the warning in issue #173.
    if (length(lockMass) %% 2)
        lockMass <- c(lockMass, -99)
    lockMass<-matrix(lockMass, ncol=2, byrow=TRUE)
    ## if((lockMass[nrow(lockMass),2]+2) > length(ob@scanindex)){
    ##     lockMass<-lockMass[1:(nrow(lockMass)-1),]
    ## }

    ## We're looping from 1 to length - 1, thus we have to fill in the last
    ## scan later.
    for(i in 1:(length(ob@scanindex)-1)){
        if(any(i == lockMass[, 1])){
            ## Place mz and intensity values from the previous scan into the
            ## array and fill the rest with NA.
            arr[1,,i] <-c(object@env$mz[(object@scanindex[(i-1)]+1):object@scanindex[i]],
                          rep(NA, (max(diff(object@scanindex))-
                                   length((object@scanindex[(i-1)]+1):object@scanindex[i])) ))

            arr[2,,i] <-c(object@env$intensity[(object@scanindex[(i-1)]+1):object@scanindex[i]],
                          rep(NA, (max(diff(object@scanindex)) -
                                   length((object@scanindex[(i-1)]+1):object@scanindex[i])) ))

        } else if(any(i == lockMass[, 2])){
            ## Place mz and intensity values from the next scan into the array.
            arr[1,,i] <-c(object@env$mz[(object@scanindex[i+1]+1):object@scanindex[(i+2)]],
                          rep(NA, (max(diff(object@scanindex)) -
                                   length((object@scanindex[i+1]+1):object@scanindex[(i+2)])) ))

            arr[2,,i] <-c(object@env$intensity[(object@scanindex[i+1]+1):object@scanindex[(i+2)]],
                          rep(NA, (max(diff(object@scanindex)) -
                                   length((object@scanindex[i+1]+1):object@scanindex[(i+2)])) ))

        } else{
            ## Just fill with the actual values.
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
    ## Fix for #173: fill also values for the last scan.
    last_i <- length(ob@scanindex)
    fetch_idx <- (object@scanindex[last_i] + 1):length(object@env$mz)
    put_idx <- 1:length(fetch_idx)
    arr[1, put_idx, length(ob@scanindex)] <- object@env$mz[fetch_idx]
    arr[2, put_idx, length(ob@scanindex)] <- object@env$intensity[fetch_idx]

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

############################################################
## [
## Subset by scan.
#' @title Subset an xcmsRaw object by scans
#' 
#' @aliases subset-xcmsRaw
#'
#' @description Subset an \code{\linkS4class{xcmsRaw}} object by scans. The
#'     returned \code{\linkS4class{xcmsRaw}} object contains values for all
#'     scans specified with argument \code{i}. Note that the \code{scanrange}
#'     slot of the returned \code{xcmsRaw} will be
#'     \code{c(1, length(object@scantime))} and hence not \code{range(i)}.
#'
#' @details Only subsetting by scan index in increasing order or by a logical
#'     vector are supported. If not ordered, argument \code{i} is sorted
#'     automatically. Indices which are larger than the total number of scans
#'     are discarded.
#' 
#' @param x The \code{\linkS4class{xcmsRaw}} object that should be sub-setted.
#' 
#' @param i Integer or logical vector specifying the scans/spectra to which
#'     \code{x} should be sub-setted.
#' 
#' @param j Not supported.
#' 
#' @param drop Not supported.
#' 
#' @return The sub-setted \code{\linkS4class{xcmsRaw}} object.
#' 
#' @author Johannes Rainer
#' 
#' @seealso \code{\link{split.xcmsRaw}}
#' 
#' @examples
#' ## Load a test file
#' file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
#' xraw <- xcmsRaw(file)
#' ## The number of scans/spectra:
#' length(xraw@scantime)
#'
#' ## Subset the object to scans with a scan time from 3500 to 4000.
#' xsub <- xraw[xraw@scantime >= 3500 & xraw@scantime <= 4000]
#' range(xsub@scantime)
#' ## The number of scans:
#' length(xsub@scantime)
#' ## The number of values of the subset:
#' length(xsub@env$mz)
setMethod("[", signature(x = "xcmsRaw",
                         i = "logicalOrNumeric",
                         j = "missing",
                         drop = "missing"),
          function(x, i, j, drop) {
              nScans <- length(x@scantime)
              if (is.logical(i))
                  i <- which(i)
              if (length(i) == 0)
                  return(new("xcmsRaw"))
              ## Ensure i is sorted.
              if (is.unsorted(i)) {
                  warning("'i' has been sorted as subsetting supports",
                          " only sorted indices.")
                  i <- sort(i)
              }
              ## Ensure that we're not supposed to return duplicated scans!
              i <- unique(i)
              ## Check that i is within the number of scans
              outsideRange <- !(i %in% 1:nScans)
              if (any(outsideRange)) {
                  ## Throw an error or just a warning?
                  iOut <- i[outsideRange]
                  i <- i[!outsideRange]
                  warning("Some of the indices in 'i' are larger than the",
                          " total number of scans which is ", nScans)
              }
              ## If i is equal to the number of scans just return the object
              if (all(1:nScans %in% i))
                  return(x)
              have_profstep <- profStep(x)
              valsPerSpect <- diff(c(x@scanindex, length(x@env$mz)))
              ## Subset:
              ## 1) scantime
              x@scantime <- x@scantime[i]
              ## 2) scanindex
              x@scanindex <- valueCount2ScanIndex(valsPerSpect[i])
              ## 3) @env$mz
              newE <- new.env()
              valsInScn <- rep(1:nScans, valsPerSpect) %in% i
              newE$mz <- x@env$mz[valsInScn]
              ## 4) @env$intensity
              newE$intensity <- x@env$intensity[valsInScn]
              x@env <- newE
              ## 5) The remaining slots
              x@tic <- x@tic[i]
              if (length(x@polarity) == nScans)
                  x@polarity <- x@polarity[i]
              if (length(x@acquisitionNum) == nScans)
                  x@acquisitionNum <- x@acquisitionNum[i]
              x@mzrange <- range(x@env$mz)
              scanrange(x) <- c(1, length(x@scantime))
              ## Profile matrix
              if (have_profstep > 0) {
                  ## Create the profile matrix.
                  ## Call profStep<- that will also update the mzrange properly
                  profStep(x) <- have_profstep
              }
              ## TODO: what with the MSn data?
              return(x)
          })

#' @title The profile matrix
#'
#' @aliases profile-matrix profMat profMat,xcmsRaw-method
#'
#' @description The \emph{profile} matrix is an n x m matrix, n (rows)
#'     representing equally spaced m/z values (bins) and m (columns) the
#'     retention time of the corresponding scans. Each cell contains the maximum
#'     intensity measured for the specific scan and m/z values falling within
#'     the m/z bin.
#'
#'     The \code{profMat} method creates a new profile matrix or returns the
#'     profile matrix within the object's \code{@env} slot, if available.
#'     Settings for the profile matrix generation, such as \code{step} (the bin
#'     size), \code{method} or additional settings are extracted from the
#'     respective slots of the \code{\linkS4class{xcmsRaw}} object.
#'     Alternatively it is possible to specify all of the settings as
#'     additional parameters.
#'
#' @details Profile matrix generation methods:
#'     \describe{
#'     \item{bin}{The default profile matrix generation method that does a
#'     simple binning, i.e. aggregating of intensity values falling within an
#'     m/z bin.}
#'     \item{binlin}{Binning followed by linear interpolation to impute missing
#'     values. The value for m/z bins without a measured intensity are inferred
#'     by a linear interpolation between neighboring bins with a measured
#'     intensity.}
#'     \item{binlinbase}{Binning followed by a linear interpolation to impute
#'     values for empty elements (m/z bins) within a user-definable proximity to
#'     non-empty elements while stetting the element's value to the
#'     \code{baselevel} otherwise. See \code{impute = "linbase"} parameter of
#'     \code{\link{imputeLinInterpol}} for more details.}
#'     \item{intlin}{Set the elements' values to the integral of the linearly
#'     interpolated data from plus to minus half the step size.}
#'     }
#'
#' @note From \code{xcms} version 1.51.1 on only the \code{profMat} method
#'     should be used to extract the profile matrix instead of the previously
#'     default way to access it directly \emph{via} \code{object@env$profile}.
#'
#' @param object The \code{\linkS4class{xcmsRaw}} object.
#'
#' @param method The profile matrix generation method. Allowed are \code{"bin"},
#'     \code{"binlin"}, \code{"binlinbase"} and \code{"intlin"}. See details
#'     section for more information.
#'
#' @param step numeric(1) representing the m/z bin size.
#'
#' @param baselevel numeric(1) representing the base value to which
#'     empty elements (i.e. m/z bins without a measured intensity) should be
#'     set. Only considered if \code{method = "binlinbase"}. See
#'     \code{baseValue} parameter of \code{\link{imputeLinInterpol}} for more
#'     details.
#'
#' @param basespace numeric(1) representing the m/z length after
#'     which the signal will drop to the base level. Linear interpolation will
#'     be used between consecutive data points falling within
#'     \code{2 * basespace} to each other. Only considered if
#'     \code{method = "binlinbase"}. If not specified, it defaults to
#'     \code{0.075}. Internally this parameter is translated into the
#'     \code{distance} parameter of the \code{\link{imputeLinInterpol}}
#'     function by \code{distance = floor(basespace / step)}. See
#'     \code{distance} parameter of \code{\link{imputeLinInterpol}} for more
#'     details.
#'
#' @param mzrange. Optional numeric(2) manually specifying the mz value range to
#'     be used for binnind. If not provided, the whole mz value range is used.
#'
#' @seealso \code{\linkS4class{xcmsRaw}}, \code{\link{binYonX}} and
#'     \code{\link{imputeLinInterpol}} for the employed binning and
#'     missing value imputation methods, respectively.
#'     \code{\link{profMat,XCMSnExp-method}} for the method on
#'     \code{\link{XCMSnExp}} objects.
#'
#' @return \code{profMat} returns the profile matrix (rows representing scans,
#'     columns equally spaced m/z values).
#'
#' @author Johannes Rainer
#'
#' @examples
#' file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
#' ## Load the data without generating the profile matrix (profstep = 0)
#' xraw <- xcmsRaw(file, profstep = 0)
#' ## Extract the profile matrix
#' profmat <- profMat(xraw, step = 0.3)
#' dim(profmat)
#' ## If not otherwise specified, the settings from the xraw object are used:
#' profinfo(xraw)
#' ## To extract a profile matrix with linear interpolation use
#' profmat <- profMat(xraw, step = 0.3, method = "binlin")
#' ## Alternatively, the profMethod of the xraw objects could be changed
#' profMethod(xraw) <- "binlin"
#' profmat_2 <- profMat(xraw, step = 0.3)
#' all.equal(profmat, profmat_2)
#'
#' @rdname profMat-xcmsSet
#' 
#' @name profMat-xcmsSet
setMethod("profMat", signature(object = "xcmsRaw"), function(object, method,
                                                             step,
                                                             baselevel,
                                                             basespace,
                                                             mzrange.) {
    ## Call the .createProfileMatrix if the internal profile matrix is empty or
    ## if the settings differ from the internal settings.
    pi <- profinfo(object)
    if (missing(method))
        method <- pi$method
    if (missing(step))
        step <- pi$step
    if (missing(baselevel))
        baselevel <- pi$baselevel
    if (missing(basespace))
        basespace <- pi$basespace
    if (length(object@env$profile) > 0) {
        ## Check if the settings are the same...
        have_method <- pi$method
        have_step <- pi$step
        have_basespace <- pi$basespace
        have_baselevel <- pi$baselevel
        if (!is.character(all.equal(list(method, step, basespace, baselevel),
                                    list(have_method, have_step,
                                         have_basespace, have_baselevel)))) {
            return(object@env$profile)
        }
    }
    if (missing(mzrange.)) {
        mzrange. <- NULL
    } else {
        if (length(mzrange.) != 2 | !is.numeric(mzrange.))
            stop("If provided, 'mzrange.' has to be a numeric of length 2!")
    }
    ## Calculate the profile matrix
    if (step <= 0)
        stop("Can not calculate a profile matrix with step=",step, "! Either",
             " provide a positive number with argument 'step' or use the",
             " 'profStep' method to set the step parameter for the xcmsSet.")
    message("Create profile matrix with method '", method,
            "' and step ", step, " ... ", appendLF = FALSE)
    res <- .createProfileMatrix(mz = object@env$mz, int = object@env$intensity,
                                valsPerSpect = diff(c(object@scanindex,
                                                      length(object@env$mz))),
                                method = method, step = step,
                                baselevel = baselevel, basespace = basespace,
                                mzrange. = mzrange.)
    message("OK")
    res
})
