require(methods) || stop("Couldn't load package methods")

# FIXME: What about a 'clone' method? Obviously not a good idea for a huge
# profile matrix, but for smaller matrices it would be very convenient
setClass("xcmsRaw", representation(env = "environment", tic = "numeric",
                                   scantime = "numeric", scanindex = "integer",
                                   acquisitionNum = "integer",
                                   pipeline = "xcmsRawPipeline",
                                   mzrange = "numeric", gradient = "matrix",
                                   msnScanindex = "integer",
                                   msnAcquisitionNum = "integer",
                                   msnPrecursorScan = "integer",
                                   msnLevel = "integer",
                                   msnRt = "numeric",
                                   msnPrecursorMz = "numeric",
                                   msnPrecursorIntensity = "numeric",
                                   msnPrecursorCharge = "numeric",
                                   msnPeakCount = "integer",
                                   filepath = "character"),
         prototype(env = new.env(parent=.GlobalEnv), tic = numeric(0),
                   scantime = numeric(0), scanindex = integer(0),
                   acquisitionNum = integer(0),
                   pipeline = new("xcmsRawPipeline"),
                   mzrange = numeric(0),
                   gradient = matrix(nrow=0, ncol=0),
                   msnScanindex = NULL,
                   msnAcquisitionNum = integer(0),
                   msnLevel = NULL,
                   msnRt = NULL,
                   msnPrecursorScan = NULL,
                   msnPrecursorMz = NULL,
                   msnPrecursorIntensity = NULL,
                   msnPrecursorCharge = NULL,
                   msnPeakCount = NULL
                   ))

xcmsRaw <- function(filename, profstep = 1, profmethod = "intlin",
                    profparam = list(), genprof = TRUE,
                    includeMSn = FALSE, pipeline = NULL) {

    object <- new("xcmsRaw")
    object@env <- new.env(parent=.GlobalEnv)

    if (!file.exists(filename)) stop("File ",filename, " not exists. \n"   )
    if (netCDFIsFile(filename)) {
        cdf <- netCDFOpen(filename)
        if (!is.null(attr(cdf, "errortext")))
            stop(attr(cdf, "errortext"))
        on.exit(netCDFClose(cdf))
        rawdata <- netCDFRawData(cdf)
        if (includeMSn) {
            warning("Reading of MSn spectra for NetCDF not supported")
        }
    } else if (rampIsFile(filename)) {
        rampid <- rampOpen(filename)
        if (rampid < 0)
            stop("Couldn't open mzXML/mzData file")
        on.exit(rampClose(rampid))

        rawdata <- rampRawData(rampid)

        if ( includeMSn ) {
            rawdataMSn <- rampRawDataMSn(rampid)
        }
    } else
        stop("Couldn't determine file type")

    rtdiff <- diff(rawdata$rt)
    if (any(rtdiff < 0)) {
    	badtimes <- which(rtdiff < 0)
    	stop(paste("Time for scan ", badtimes[1], " (",
    	           rawdata$rt[[badtimes[1]]], ") greater than scan ",
    	           badtimes[1]+1, " (", rawdata$rt[[badtimes[1]+1]], ")",
    	           sep = ""))
    }

    object@filepath <- filename
    object@scantime <- rawdata$rt
    object@tic <- rawdata$tic
    object@scanindex <- rawdata$scanindex
    object@env$mz <- rawdata$mz
    object@env$intensity <- rawdata$intensity

    if (!is.null(rawdata$acquisitionNum)) {
      ## defined only for mzData and mzXML
      object@acquisitionNum <- rawdata$acquisitionNum
    }

    if(exists("rawdataMSn")) {
        object@env$msnMz <- rawdataMSn$mz
        object@env$msnIntensity <- rawdataMSn$intensity

        object@msnScanindex <- rawdataMSn$scanindex
        object@msnAcquisitionNum <- rawdataMSn$acquisitionNum
        object@msnLevel <- rawdataMSn$msLevel
        object@msnRt <- rawdataMSn$rt
        object@msnPrecursorScan <- match(rawdataMSn$precursorNum, object@acquisitionNum)
        object@msnPrecursorMz <- rawdataMSn$precursorMZ
        object@msnPrecursorIntensity <- rawdataMSn$precursorIntensity
        object@msnPrecursorCharge <- rawdataMSn$precursorCharge
    }

    if (is.null(pipeline)) { # create pipeline, if needed
        profargs <- c(list(profmethod = profmethod, profstep = profstep), profparam)
        genprofproto <- do.call("xcmsProtocol", c("GenProfile", profargs))
        pipeline <- new("xcmsRawPipeline", genprofproto = genprofproto)
    } else { # if given, clone, overriding with any specified parameters
        mc <- as.list(match.call())
        specified <- match(c("profstep", "profmethod"), names(mc), 0)
        if (sum(specified)) {
            genprofproto <- genProfProto(pipeline)
            genProfProto(pipeline) <- do.call("new", c(class(genprofproto),
                genprofproto, mc[specified], profparam))
        }
    }
    realprofstep <- profStep(genProfProto(pipeline))
    if (!genprof) # temporarily set step to zero to disable profile generation
        profStep(genProfProto(pipeline)) <- 0
    pipeline(object) <- pipeline
    profStep(genProfProto(object@pipeline)) <- realprofstep

    return(object)
}

setMethod("show", "xcmsRaw", function(object) {

    cat("An \"xcmsRaw\" object with", length(object@scantime), "mass spectra\n\n")

    cat("Time range: ", paste(round(range(object@scantime), 1), collapse = "-"),
        " seconds (", paste(round(range(object@scantime)/60, 1), collapse = "-"),
        " minutes)\n", sep = "")
    cat("Mass range:", paste(round(range(object@env$mz), 4), collapse = "-"),
        "m/z\n")
    cat("Intensity range:", paste(signif(range(object@env$intensity), 6), collapse = "-"),
        "\n\n")

    ## summary MSn data
    if (!is.null(object@msnLevel)) {
	cat("MSN data on ", length(unique(object@msnPrecursorMz)), " mass(es)\n")
	cat("\twith ", length(object@msnPrecursorMz)," MSn spectra\n")
    }

    cat("Profile matrix: ")
    if (is.null(object@env$profile))
        cat("none\n")
    else {
        profmz <- profMz(object)
        cat(profStep(object), " m/z step (", length(profmz), " grid points from ",
            paste(object@mzrange, collapse = " to "), " m/z)\n", sep = "")
    }
    cat("\n")

    show(object@pipeline)

    memsize <- object.size(object)
    for (key in ls(object@env))
        memsize <- memsize + object.size(object@env[[key]])
    cat("\nMemory usage:", signif(memsize/2^20, 3), "MB\n")
})

setGeneric("revMz", function(object, ...) standardGeneric("revMz"))

setMethod("revMz", "xcmsRaw", function(object) {

    for (i in 1:length(object@scanindex)) {
        idx <- (object@scanindex[i]+1):min(object@scanindex[i+1],
                                         length(object@env$mz), na.rm=TRUE)
        object@env$mz[idx] <- rev(object@env$mz[idx])
        object@env$intensity[idx] <- rev(object@env$intensity[idx])
    }
})

setGeneric("sortMz", function(object, ...) standardGeneric("sortMz"))

setMethod("sortMz", "xcmsRaw", function(object) {

    for (i in 1:length(object@scanindex)) {
        idx <- (object@scanindex[i]+1):min(object@scanindex[i+1],
                                         length(object@env$mz), na.rm=TRUE)
        ord <- order(object@env$mz[idx])
        object@env$mz[idx] <- object@env$mz[idx[ord]]
        object@env$intensity[idx] <- object@env$intensity[idx[ord]]
    }
})

setGeneric("plotTIC", function(object, ...) standardGeneric("plotTIC"))

setMethod("plotTIC", "xcmsRaw", function(object, ident = FALSE, msident = FALSE) {

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
                eval(call(options("device")$device))
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

setGeneric("getScan", function(object, ...) standardGeneric("getScan"))

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

setGeneric("getSpec", function(object, ...) standardGeneric("getSpec"))

setMethod("getSpec", "xcmsRaw", function(object, ...) {

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

specNoise <- function(spec, gap = quantile(diff(spec[,"mz"]), .9)) {

    intmean <- mean(spec[,"intensity"])

    mzlen <- diff(range(spec[,"mz"]))
    mzdiff <- diff(spec[,"mz"])
    gaplen <- sum(mzdiff[mzdiff > gap])

    weighted.mean(c(intmean, min(spec[,"intensity"])/2), c(1 - gaplen/mzlen,
                                                           gaplen/mzlen))
}

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

setGeneric("plotScan", function(object, ...) standardGeneric("plotScan"))

setMethod("plotScan", "xcmsRaw", function(object, scan, mzrange = numeric(),
                                          ident = FALSE) {

    if (object@scanindex[scan] == length(object@env$mz) ||
        object@scanindex[scan] == object@scanindex[scan+1])
        return()
    idx <- (object@scanindex[scan]+1):min(object@scanindex[scan+1],
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

setGeneric("plotSpec", function(object, ...) standardGeneric("plotSpec"))

setMethod("plotSpec", "xcmsRaw", function(object, ident = FALSE,
                                          vline = numeric(0), ...) {

    sel <- profRange(object, ...)

    title = paste("Averaged Mass Spectrum: ", sel$rtlabel, " (",
                  sel$scanlab, ")",  sep = "")
    points <- cbind(profMz(object)[sel$mzindex],
                    rowMeans(object@env$profile[sel$mzindex,sel$scanidx,drop=FALSE]))
    plot(points, type="l", main = title, xlab="m/z", ylab="Intensity")
    if (length(vline))
        abline(v = vline, col = "red")

    if (ident)
        return(identify(points, labels = round(points[,1], 1)))

    invisible(points)
})

setGeneric("plotChrom", function(object, ...) standardGeneric("plotChrom"))

setMethod("plotChrom", "xcmsRaw", function(object, base = FALSE, ident = FALSE,
                                           fitgauss = FALSE, vline = numeric(0), ...) {

    sel <- profRange(object, ...)


    if (base) {
        title = paste("Base Peak Chromatogram: ", sel$mzlabel, sep = "")
        pts <- cbind(object@scantime[sel$scanidx],
                     colMax(object@env$profile[sel$mzindex,sel$scanidx,drop=FALSE]))
    }
    else {
        title = paste("Averaged Ion Chromatogram: ", sel$mzlabel, sep = "")
        pts <- cbind(object@scantime[sel$scanidx],
                     colMeans(object@env$profile[sel$mzindex,sel$scanidx,drop=FALSE]))
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

image.xcmsRaw <- function(x, col = rainbow(256), ...) {

    sel <- profRange(x, ...)

    zlim <- log(range(x@env$intensity))

    title <- paste("XC/MS Log Intensity Image (Profile Method: ",
                   profMethod(x), ")", sep = "")
    if (zlim[1] < 0) {
        zlim <- log(exp(zlim)+1)
        image(profMz(x)[sel$mzindex], x@scantime[sel$scanidx],
              log(x@env$profile[sel$mzindex, sel$scanidx]+1),
              col = col, zlim = zlim, main = title, xlab="m/z", ylab="Seconds")
    } else
        image(profMz(x)[sel$mzindex], x@scantime[sel$scanidx],
              log(x@env$profile[sel$mzindex, sel$scanidx]),
              col = col, zlim = zlim, main = title, xlab="m/z", ylab="Seconds")
}

setGeneric("plotSurf", function(object, ...) standardGeneric("plotSurf"))

setMethod("plotSurf", "xcmsRaw", function(object, log = FALSE,
                                          aspect = c(1, 1, .5), ...) {

    require(rgl) || stop("Couldn't load package rgl")

    sel <- profRange(object, ...)

    y <- object@env$profile[sel$mzindex, sel$scanidx]
    if (log)
        y <- log(y+max(1-min(y), 0))
    ylim <- range(y)

    x <- seq(0, aspect[1], length=length(sel$mzindex))
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

filtfft <- function(y, filt) {

    yfilt <- numeric(length(filt))
    yfilt[1:length(y)] <- y
    yfilt <- fft(fft(yfilt, inverse = TRUE) * filt)

    Re(yfilt[1:length(y)])
}

.findPeaks.matchedFilter <- function(object, fwhm = 30, sigma = fwhm/2.3548,
                                     max = 5, snthresh = 10, step = 0.1,
                                     steps = 2, mzdiff = 0.8 - step*steps,
                                     index = FALSE, sleep = 0,
                                     verbose.columns = FALSE) {
    pipeline <- object@pipeline
    profStep(genProfProto(pipeline)) <- step
    maxidxproto <- genProfProto(pipeline)
    profMethod(maxidxproto) <- "maxidx"
    maxidx <- new("xcmsRawPipeline", genprofproto = maxidxproto)

    ### Create EIC buffer
    # This buffered reading should be factored into an iterator object
    mrange <- range(object@env$mz)
    mass <- seq(floor(mrange[1]/step)*step, ceiling(mrange[2]/step)*step, by = step)
    bufsize <- min(100, length(mass))
    buf <- profMat(object, pipeline, mzrange = c(mass[1],mass[bufsize]))
    bufMax <- profMat(object, maxidx, mzrange = c(mass[1],mass[bufsize]))
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
            buf <- profMat(object, pipeline,
                mzrange = c(mass[idxrange[1]], mass[idxrange[2]]))
            bufMax <- profMat(object, maxidx,
                mzrange = c(mass[idxrange[1]], mass[idxrange[2]]))
        }
        ymat <- buf[bufidx[i:(i+steps-1)],,drop=FALSE]
        ysums <- colMax(ymat)
        yfilt <- filtfft(ysums, filt)
        gmax <- max(yfilt)
        for (j in seq(length = max)) {
             maxy <- which.max(yfilt)
             noise <- mean(ysums[ysums > 0])
             #noise <- mean(yfilt[yfilt >= 0])
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
                 massrange <- range(mzmat, na.rm = TRUE)
                 massmean <- weighted.mean(mzmat, intmat[which.intMax], na.rm = TRUE)
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
                 rmat[num,] <- c(massmean, massrange[1], massrange[2], maxy, peakrange, into, intf, maxo, maxf, j, sn)
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
    invisible(rmat)
}

setProtocolClass("xcmsProtoFindPeaksMatchedFilter",
  representation(fwhm = "numeric", sigma = "numeric", max = "numeric",
    snthresh = "numeric", step = "numeric", steps = "numeric", mzdiff = "numeric",
    index = "logical"),
  c(formals(.findPeaks.matchedFilter), dispname = "Matched Filter"),
  "xcmsProtoFindPeaks")

setGeneric("findPeaks.matchedFilter", function(object, ...) standardGeneric("findPeaks.matchedFilter"))

setMethod("findPeaks.matchedFilter", "xcmsRaw", .findPeaks.matchedFilter)

.findPeaks.centWave <- function(object, scanrange=c(1,length(object@scantime)),
                                minEntries=4, dev=140e-6, snthresh=20, minPeakWidth=7,
                                noiserange=c(minPeakWidth*3,minPeakWidth*6),
                                scales=c(5,7,9,12,16,20), maxGaussOverlap = 0.5,
                                minPtsAboveBaseLine=4, scRangeTol=2,
                                maxDescOutlier=floor(minPeakWidth/2), mzdiff=-0.001,
                                rtdiff=-round(2/3 *minPeakWidth *mean(diff(object@scantime))),
                                integrate=1, sleep=0, fitgauss = FALSE, verbose.columns = FALSE)
{
    peaklist <- list()
    featlist <- findMZBoxes(object,scanrange=scanrange,dev=dev,minEntries=minEntries)
    scantime <- object@scantime
    Nscantime <- length(scantime)
    lf <- length(featlist)
    cat('\n Searching for peaks... \n % finished: '); lp <- -1;

    for (f in  1:lf) {
      feat <- featlist[[f]]
      perc <- round((f/lf) * 100)
      if ((perc %% 10 == 0) && (perc != lp)) { cat(perc,' '); lp <- perc }
      flush.console()
      N <- length(feat$mz)
      peaks <- peakinfo <- NULL
      mzrange <- range(feat$mz)
      sccenter <- feat$scan[1] + floor(N/2) - 1
      scrange <- range(feat$scan)
      ## scrange + noiserange, used for baseline detection and wavelet analysis
      sr <- c(max(scanrange[1],scrange[1] - max(noiserange)),min(scanrange[2],scrange[2] + max(noiserange)))
      eic <- rawEIC(object,mzrange=mzrange,scanrange=sr)
      d <- eic$intensity
      td <- sr[1]:sr[2]
      scan.range <- c(sr[1],sr[2])
      ## original data range (hd m/z boxes)
      omz <- feat$mz
      od <- feat$intensity
      otd <- feat$scan
      ##  scrange + scRangeTol, used for gauss fitting and continuous data above 1st baseline detection
      ftd <- max(td[1], scrange[1] - scRangeTol) : min(td[length(td)], scrange[2] + scRangeTol)
      fd <- d[match(ftd,td)]

      ## 1st type of baseline: statistic approach
      if (N >= 10*minPeakWidth)  ## in case of very long mass trace use full scan range for baseline detection
        noised <- rawEIC(object,mzrange=mzrange,scanrange=scanrange)$intensity else
            noised <- d;
      ## 90% trimmed mean as first baseline guess
      noise <- estimateChromNoise(noised,c(0.05,0.95),minPts=3*minPeakWidth)

      ## any continuous data above 1st baseline ?
      if (continuousPtsAboveThreshold(fd,threshold=noise,num=minPtsAboveBaseLine)) {
          ## 2nd baseline estimate using not-peak-range
          lnoise <- getLocalNoiseEstimate(d,td,ftd,noiserange,Nscantime)

          ## Final baseline & Noise estimate
          baseline <- max(1,min(lnoise[1],noise))
          sdnoise <- max(1,lnoise[2])
          sdnoise10 <-  sdnoise * 10^(snthresh/20)

        ## is there any data above S/N * threshold ?

        if (any(fd - baseline >= sdnoise10 )) {
            wCoefs <- MSW.cwt(d, scales=scales, wavelet='mexh')
            if (!is.null(dim(wCoefs)) && any(wCoefs- baseline >= sdnoise10)) {
                if (td[length(td)] == Nscantime) ## workaround, localMax fails otherwise
                    wCoefs[nrow(wCoefs),] <- wCoefs[nrow(wCoefs)-1,] * 0.99
                localMax <- MSW.getLocalMaximumCWT(wCoefs)
                rL <- MSW.getRidge(localMax)
                wpeaks <- sapply(rL,
                    function(x) {
                        w <- min(1:length(x),ncol(wCoefs))
                        any(wCoefs[x,w]- baseline >= sdnoise10)
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
                          if (any(d[pp[dv]]- baseline >= sdnoise10)) {
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
                                   bestcol <- which(m == max(m),arr=T)[2]
                                   best.scale.nr <- maxpi[bestcol]
                                } else  best.scale.nr <- maxpi

                              best.scale <-  scales[best.scale.nr]
                              best.scale.pos <- opp[best.scale.nr]

                              pprange <- min(pp):max(pp)
                              maxint <- max(d[pprange])
                              lwpos <- max(1,best.scale.pos - best.scale)
                              rwpos <- min(best.scale.pos + best.scale,length(td))
                              p1 <- match(td[lwpos],otd)[1]
                              p2 <- match(td[rwpos],otd); p2 <- p2[length(p2)]
                              if (is.na(p1)) p1<-1
                              if (is.na(p2)) p2<-N
                              mz.value <- omz[p1:p2]
                              mz.int <- od[p1:p2]
                              mzmean <- mzModel(mz.value,mz.int) ## re-calculate m/z value for peak range
                              mzrange <- range(mz.value)
                              if (length(mz.value) >= (minEntries+1))
                                dppm <- round(min(running(abs(diff(mz.value)) /(mzrange[2] *  1e-6),fun=max,width=minEntries))) else
                                    dppm <- round((mzrange[2]-mzrange[1]) /  (mzrange[2] *  1e-6))

                              peaks <- rbind(peaks,
                                  c(mzmean,mzrange,           ## mz
                                  NA,NA,NA,                   ## rt, rtmin, rtmax,
                                  NA,                         ## intensity (sum)
                                  NA,                         ## intensity (-bl)
                                  maxint,                     ## max intensity
                                  round(20 * log10( (maxint - baseline) / sdnoise)),  ##  S/N Ratio
                                  NA,                         ## Gaussian RMSE
                                  NA,NA,NA,                   ## Gaussian Parameters
                                  f,                          ## ROI Position
                                  dppm,                       ## max. difference between the [minEntries] peaks in ppm
                                  best.scale,                 ## Scale
                                  td[best.scale.pos], td[lwpos], td[rwpos],
                                                              ## Peak positions guessed from the wavelet's (scan nr)
                                  NA,NA ))                    ## Peak limits (scan nr)

                              peakinfo <- rbind(peakinfo,c(best.scale, best.scale.nr, best.scale.pos, lwpos, rwpos))  ## Peak positions guessed from the wavelet's
                          }
                        }
                      }
                    }  #for
                } # if
            }

            ##  postprocessing
            if (!is.null(peaks)) {
               # if (is.vector(peaks)) peaks <- data.frame(t(peaks))  else peaks <- data.frame(peaks)
                basenames <- c("mz","mzmin","mzmax","rt","rtmin","rtmax","into","intb","maxo","sn")
                colnames(peaks) <- c(basenames,"egauss","mu","sigma","h","f", "dppm", "scale","scpos","scmin","scmax","lmin","lmax")

               # if (is.vector(peakinfo)) peakinfo <- data.frame(t(peakinfo))  else peakinfo <- data.frame(peakinfo)
                colnames(peakinfo) <- c("scale","scaleNr","scpos","scmin","scmax")

                for (p in 1:dim(peaks)[1]) {
                  ## find minima, assign rt and intensity values
                  if (integrate == 1)
                    lm <- descendMin(wCoefs[,peakinfo[p,"scaleNr"]], istart= peakinfo[p,"scpos"]) else
                        lm <- descendMinTol(d,startpos=c(peakinfo[p,"scmin"],peakinfo[p,"scmax"]),maxDescOutlier)

                  peakrange <- td[lm]
                  peaks[p,"rtmin"] <- scantime[peakrange[1]]
                  peaks[p,"rtmax"] <- scantime[peakrange[2]]
                  pwid <- (scantime[peakrange[2]] - scantime[peakrange[1]])/(peakrange[2] - peakrange[1])
                  peaks[p,"into"] <- pwid*sum(d[lm[1]:lm[2]])

                  db <-  d[lm[1]:lm[2]] - baseline
                  peaks[p,"intb"] <- pwid*sum(db[db>0])
                  peaks[p,"lmin"] <- lm[1]; peaks[p,"lmax"] <- lm[2];

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
                peaks <- joinOverlappingPeaks(td,d,otd,omz,od,scantime,scan.range,peaks,maxGaussOverlap)
            }

        } ## s/n
      } ## minPtsAboveBaseLine valid

      if ((sleep >0) && (!is.null(peaks))) {
            tdp <- scantime[td]; trange <- range(tdp)
            egauss <- paste(round(peaks[,"egauss"],3),collapse=", ")
            cdppm <- paste(peaks[,"dppm"],collapse=", ")
            csn <- paste(peaks[,"sn"],collapse=", ")
            par(bg = "white")
            l <- layout(matrix(c(1,2,3),nr=3,nc=1,byrow=T),heights=c(.5,.75,2));
            par(mar= c(2, 4, 4, 2) + 0.1)
            plotRaw(object,mass=mzrange,time=trange,log=T,title='')
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
               #abline(v=scantime[na.omit(c(peaks[(peaks[,"scmin"] >0),"scmin"],peaks[(peaks[,"scmax"] >0),"scmax"]))],col='cyan')
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

      if (!is.null(peaks)) peaklist[[length(peaklist)+1]] <- peaks

    } # f
    cat('\n')
    p <- do.call("rbind",peaklist)

    if (!verbose.columns)
        p <- p[,basenames,drop=FALSE]

    uorder <- order(p[,"into"], decreasing=TRUE)
    pm <- as.matrix(p[,c("mzmin","mzmax","rtmin","rtmax"),drop=FALSE])
    uindex <- rectUnique(pm,uorder,mzdiff,rtdiff)
    pr <- p[uindex,,drop=FALSE]
    cat(dim(p)[1],' Peaks  -- rectUnique(',mzdiff,',',rtdiff,') -->  ', dim(pr)[1],' Peaks.\n',sep='')

    invisible(pr) #as.matrix(pr)
}

setProtocolClass("xcmsProtoFindPeaksCentWave",
  representation(scanrange="numeric", minEntries="numeric", dev="numeric",
    snthresh="numeric", noiserange="numeric", minPeakWidth="numeric",
    scales="numeric", maxGaussOverlap = "numeric", minPtsAboveBaseLine="numeric",
    scRangeTol="numeric", maxDescOutlier="numeric", mzdiff="numeric",
    rtdiff="numeric", integrate="numeric", fitgauss = "logical"),
  c(formals(.findPeaks.centWave), dispname = "Centroid Wavelet"),
  "xcmsProtoFindPeaks")

setGeneric("findPeaks.centWave", function(object, ...) standardGeneric("findPeaks.centWave"))

setMethod("findPeaks.centWave", "xcmsRaw", .findPeaks.centWave)

.findPeaks.mS1 <- function(object, sleep=0, verbose.columns = FALSE)
{
    if (is.null(object@msnLevel)) {
        warning("xcmsRaw contains no MS2 spectra\n");
        return (NULL);
    }

    ## Select all MS2 scans, they have an MS1 parent defined
    peakIndex <- object@msnLevel == 2

    ## (empty) return object
    basenames <- c("mz","mzmin","mzmax","rt","rtmin","rtmax",
                   "into","maxo","sn")
    peaklist <- matrix(-1, nrow = length(peakIndex), ncol = length(basenames))
    colnames(peaklist) <- c(basenames)

    ## Assemble result

    peaklist[,"mz"] <- object@msnPrecursorMz[peakIndex]
    peaklist[,"mzmin"] <- object@msnPrecursorMz[peakIndex]
    peaklist[,"mzmax"] <- object@msnPrecursorMz[peakIndex]

    if (any(object@msnPrecursorScan!=0)) {
        peaklist[,"rt"] <- peaklist[,"rtmin"] <- peaklist[,"rtmax"] <- object@scantime[object@msnPrecursorScan[peakIndex]]
    } else {
        ## This happened with ReAdW mzData
        warning("MS2 spectra without precursorScan references")
        peaklist[,"rt"] <- peaklist[,"rtmin"] <- peaklist[,"rtmax"] <- 0
    }

    if (any(object@msnPrecursorIntensity!=0)) {
        peaklist[,"into"] <- peaklist[,"maxo"] <- peaklist[,"sn"] <- object@msnPrecursorIntensity[peakIndex]
    } else {
        ## This happened with Agilent MzDataExport 1.0.98.2
        warning("MS2 spectra without precursorIntensity, setting to zero")
        peaklist[,"into"] <- peaklist[,"maxo"] <- peaklist[,"sn"] <- 0
    }

    cat('\n')

    invisible(peaklist)
}

setProtocolClass("xcmsProtoFindPeaksMS1",
  representation(includeMSn="logical"),
  c(formals(.findPeaks.mS1), dispname = "MS2 precursor Peaks"),
  "xcmsProtoFindPeaks")

setGeneric("findPeaks.mS1", function(object, ...) standardGeneric("findPeaks.mS1"))

setMethod("findPeaks.mS1", "xcmsRaw", .findPeaks.mS1)

.findPeaks.MSW <- function(object, snthresh=3, scales=seq(1,22,3), nearbyPeak=TRUE,
                           SNR.method='quantile', winSize.noise=500,
                           sleep=0, verbose.columns = FALSE)
{
  require(MassSpecWavelet) || stop("Couldn't load MassSpecWavelet")

  # MassSpecWavelet Calls
  peakInfo <- peakDetectionCWT(object@env$intensity,
                                scales=scales, SNR.Th = snthresh,
                                nearbyPeak = nearbyPeak, SNR.method=SNR.method, winSize.noise=winSize.noise)
  majorPeakInfo <- peakInfo$majorPeakInfo

  betterPeakInfo <- tuneInPeakInfo(object@env$intensity,
                                    majorPeakInfo)

  peakIndex <- betterPeakInfo$peakIndex

  # Assemble result

  basenames <- c("mz","mzmin","mzmax","rt","rtmin","rtmax",
                 "into","maxo","sn")

  peaklist <- matrix(-1, nrow = length(peakIndex), ncol = length(basenames))

  colnames(peaklist) <- c(basenames)

  peaklist[,"mz"] <- object@env$mz[peakIndex]
  peaklist[,"mzmin"] <- object@env$mz[peakIndex]
  peaklist[,"mzmax"] <- object@env$mz[peakIndex]

  peaklist[,"rt"]    <- rep(-1, length(peakIndex))
  peaklist[,"rtmin"] <- rep(-1, length(peakIndex))
  peaklist[,"rtmax"] <- rep(-1, length(peakIndex))

  peaklist[,"into"] <- betterPeakInfo$peakValue
  peaklist[,"maxo"] <- object@env$intensity[peakIndex]
  peaklist[,"sn"]   <- betterPeakInfo$peakSNR

  cat('\n')

  # Filter additional (verbose) columns
  if (!verbose.columns)
    peaklist <- peaklist[,basenames,drop=FALSE]

  invisible(peaklist)
}

setProtocolClass("xcmsProtoFindPeaksMSW",
  representation(snthresh="numeric", mzdiff="numeric", scales="numeric"),
  c(formals(.findPeaks.MSW), name = "MassSpecWavelet"),
  "xcmsProtoFindPeaks")

setGeneric("findPeaks.MSW", function(object, ...) standardGeneric("findPeaks.MSW"))

setMethod("findPeaks.MSW", "xcmsRaw", .findPeaks.MSW)

setGeneric("findPeaks", function(object, ...) standardGeneric("findPeaks"))

setMethod("findPeaks", "xcmsRaw", function(object, method=xcmsMethodDefault("findPeaks"),
                                           ...) {
    perform(xcmsProtocol("findPeaks", method, ...), object)
})

setGeneric("getPeaks", function(object, ...) standardGeneric("getPeaks"))

setMethod("getPeaks", "xcmsRaw", function(object, peakrange, step = 0.1) {

    pipeline <- object@pipeline
    profStep(genProfProto(pipeline)) <- step

    if (all(c("mzmin","mzmax","rtmin","rtmax") %in% colnames(peakrange)))
        peakrange <- peakrange[,c("mzmin","mzmax","rtmin","rtmax"),drop=FALSE]
    stime <- object@scantime

    ### Create EIC buffer
    mrange <- range(peakrange[,1:2])
    mass <- seq(floor(mrange[1]/step)*step, ceiling(mrange[2]/step)*step, by = step)
    bufsize <- min(100, length(mass))
    buf <- profMat(object, pipeline, mzrange = c(mass[1], mass[bufsize]))
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
            buf <- profMat(object, pipeline, mzrange = c(mass[idxrange[1]], mass[idxrange[2]]))
        }
        ymat <- buf[bufidx[imz[1]:imz[2]],iret[1]:iret[2],drop=FALSE]
        ymax <- colMax(ymat)
        iymax <- which.max(ymax)
        pwid <- diff(stime[iret])/diff(iret)
        rmat[i,1] <- weighted.mean(mass[imz[1]:imz[2]], rowSums(ymat))
        if (is.nan(rmat[i,1]))
            rmat[i,1] <- mean(peakrange[i,1:2])
        rmat[i,2:3] <- peakrange[i,1:2]
        rmat[i,4] <- stime[iret[1]:iret[2]][iymax]
        rmat[i,5:6] <- peakrange[i,3:4]
        rmat[i,7] <- pwid*sum(ymax)
        rmat[i,8] <- ymax[iymax]
    }

    invisible(rmat)
})

setGeneric("plotPeaks", function(object, ...) standardGeneric("plotPeaks"))

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
        #main <- paste(peaks[i,"i"], " ", round(peaks[i,"mz"]),
        main <- paste(round(peaks[i,"mz"]),
                      " ", round(peaks[i,"rt"]), sep = "")
        plot(object@scantime, colMax(object@env$profile[mzi[i,],,drop=FALSE]),
             type = "l", xlim = xlim, ylim = c(0, peaks[i,"maxo"]), main = main,
             xlab = "", ylab = "", xaxt = "n", yaxt = "n")
        abline(v = peaks[i,c("rtmin","rtmax")], col = "grey")
    }
})

setGeneric("getEIC", function(object, ...) standardGeneric("getEIC"))

setMethod("getEIC", "xcmsRaw", function(object, mzrange, rtrange = NULL, step = 0.1) {

    pipeline <- object@pipeline
    profStep(genProfProto(pipeline)) <- step

    if (all(c("mzmin","mzmax") %in% colnames(mzrange)))
        mzrange <- mzrange[,c("mzmin", "mzmax"),drop=FALSE]

    ### Create EIC buffer
    mrange <- range(mzrange)
    mass <- seq(floor(mrange[1]/step)*step, ceiling(mrange[2]/step)*step, by = step)
    bufsize <- min(100, length(mass))
    buf <- profMat(object, pipeline, mzrange = c(mass[1], mass[bufsize]))
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
            buf <- profMat(object, pipeline, mzrange = c(mass[idxrange[1]], mass[idxrange[2]]))
        }
        if (missing(rtrange))
            eic[i,] <- colMax(buf[bufidx[imz[1]:imz[2]],,drop=FALSE])
        else {
            eic[[i]] <- matrix(c(object@scantime, colMax(buf[bufidx[imz[1]:imz[2]],,drop=FALSE])),
                               ncol = 2)[object@scantime >= rtrange[i,1] & object@scantime <= rtrange[i,2],,drop=FALSE]
            colnames(eic[[i]]) <- c("rt", "intensity")
        }
    }

    invisible(eic)
})

setGeneric("rawMat", function(object, ...) standardGeneric("rawMat"))

setMethod("rawMat", "xcmsRaw", function(object,
                                         mzrange = numeric(),
                                         rtrange = numeric(),
                                         scanrange = numeric(),
                                         log=FALSE) {

    if (length(rtrange) >= 2) {
        rtrange <- range(rtrange)
        scanidx <- (object@scantime >= rtrange[1]) & (object@scantime <= rtrange[2])
        scanrange <- c(match(TRUE, scanidx), length(scanidx) - match(TRUE, rev(scanidx)))
    }
    else if (length(scanrange) < 2)
        scanrange <- c(1, length(object@scantime))
    else scanrange <- range(scanrange)
    startidx <- object@scanindex[scanrange[1]] + 1
    endidx <- length(object@env$mz)
    if (scanrange[2] < length(object@scanindex))
        endidx <- object@scanindex[scanrange[2] + 1]
    #scans <- integer(endidx - startidx + 1)
    scans <- rep(scanrange[1]:scanrange[2], diff(c(object@scanindex, length(object@env$mz)+1)))
    #for (i in scanrange[1]:scanrange[2]) {
    #    idx <- (object@scanindex[i] + 1):min(object@scanindex[i +
    #        1], length(object@env$mz), na.rm = TRUE)
    #    scans[idx - startidx + 1] <- i
    #}
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
    if (log)
        y <- log(y + max(1 - min(y), 0))

    cbind(time = object@scantime[scans[massidx]], mz = masses[massidx], intensity = y)
})

setGeneric("plotRaw", function(object, ...) standardGeneric("plotRaw"))

setMethod("plotRaw", "xcmsRaw", function(object,
                                         mzrange = numeric(),
                                         rtrange = numeric(),
                                         scanrange = numeric(),
                                         log=FALSE,title='Raw Data' ) {

    raw <- rawMat(object, mzrange, rtrange, scanrange, log)

    y <- raw[,"intensity"]
    ylim <- range(y)
    y <- y/ylim[2]
    colorlut <- terrain.colors(16)
    col <- colorlut[y*15+1]

    plot(cbind(raw[,"time"], raw[,"mz"]), pch=20, cex=.5,
        main = title, xlab="Seconds", ylab="m/z", col=col,
        xlim=range(raw[,"time"]), ylim=range(raw[,"mz"]))

    invisible(raw)
})

setGeneric("profMz", function(object) standardGeneric("profMz"))

setMethod("profMz", "xcmsRaw", function(object) {

    object@mzrange[1]+profStep(object)*(0:(dim(object@env$profile)[1]-1))
})


setMethod("profMethod", "xcmsRaw", function(object) {

    profMethod(genProfProto(object@pipeline))
})

setReplaceMethod("profMethod", "xcmsRaw", function(object, value) {

    profMethod(genProfProto(object)) <- value
    object
})

setMethod("profStep", "xcmsRaw", function(object) {

    profStep(genProfProto(object@pipeline))
})

setReplaceMethod("profStep", "xcmsRaw", function(object, value) {

    profStep(genProfProto(object)) <- value
    object
})

setGeneric("profMedFilt", function(object, ...) standardGeneric("profMedFilt"))

setMethod("profMedFilt", "xcmsRaw", function(object, massrad = 0, scanrad = 0) {

    contdim <- dim(object@env$profile)
    object@env$profile <- medianFilter(object@env$profile, massrad, scanrad)
})

setGeneric("profRange", function(object, ...) standardGeneric("profRange"))

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
        mzindex <- which((contmass >= mzrange[1]) & (contmass <= mzrange[2]))
    } else {
        if (length(mzrange) == 0) {
            mzrange <- range(object@env$mz)
        } else {
            mzrange <- c(min(mzrange), max(mzrange))
        }
        mzindex <- integer()
    }
    if (mzrange[1] == mzrange[2])
        mzlabel <- paste(mzrange[1], "m/z")
    else
        mzlabel <- paste(mzrange[1], "-", mzrange[2], " m/z", sep="")


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
        rtlabel <- paste(round(rtrange[1],1), "seconds")
    else
        rtlabel <- paste(round(rtrange[1],1), "-", round(rtrange[2],1), " seconds", sep="")


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

    list(mzrange = mzrange, mzlabel = mzlabel, mzindex = mzindex,
         scanrange = scanrange, scanlab = scanlab, scanidx = scanidx,
         rtrange = rtrange, rtlabel = rtlabel)
})

# generate a profile matrix using a given pipeline and subset
setGeneric("profMat", function(object, pipeline = object@pipeline, ...)
    standardGeneric("profMat"))
setMethod("profMat", "xcmsRaw", function(object, pipeline = object@pipeline, ...)
{
    filtprofprotos <- filtProfProtos(pipeline)

    # FIXME: handle margins

    #if (length(filtprofprotos))
    #  margins <- rowMax(sapply(filtprofprotos, profMargins,
    #    mzcount = num, scancount = length(object@scantime)))

    # FIXME: need some sort of caching mechanism
    prof <- perform(genProfProto(pipeline), object, ...)

    if (is.null(prof)) # 'NULL' returned if profile protocol is no-op
        return(NULL)

    for (proto in filtprofprotos)
        prof <- perform(proto, prof, ...)

    prof
})

setGeneric("rawEIC", function(object, ...) standardGeneric("rawEIC"))

setMethod("rawEIC", "xcmsRaw", function(object,mzrange,scanrange=c(1,length(object@scantime))){

  if (!is.double(object@env$mz))  object@env$mz <- as.double(object@env$mz)
  if (!is.double(object@env$intensity)) object@env$intensity <- as.double(object@env$intensity)
  if (!is.integer(object@scanindex)) object@scanindex <- as.integer(object@scanindex)

  .Call("getEIC",object@env$mz,object@env$intensity,object@scanindex,as.double(mzrange),as.integer(scanrange),as.integer(length(object@scantime)), PACKAGE ='xcms' )
})

setGeneric("findMZBoxes", function(object, ...) standardGeneric("findMZBoxes"))

setMethod("findMZBoxes", "xcmsRaw", function(object,mzrange=c(0.0,0.0),scanrange=c(1,length(object@scantime)),dev,minEntries,debug=0){

  ## mzrange not implemented yet
  if (!is.double(object@env$mz))  object@env$mz <- as.double(object@env$mz)
  if (!is.double(object@env$intensity)) object@env$intensity <- as.double(object@env$intensity)
  if (!is.integer(object@scanindex)) object@scanindex <- as.integer(object@scanindex)

  .Call("findmzboxes",object@env$mz,object@env$intensity,object@scanindex,as.double(mzrange),
  as.integer(scanrange),as.integer(length(object@scantime)),
  as.double(dev),as.integer(minEntries),as.integer(debug), PACKAGE ='xcms' )
})

# Profile filtering

# Add a profile filter protocol
setMethod("addFiltProfProtos", c("xcmsRaw", "list"), function(object, value)
{
  object@pipeline <- addFiltProfProtos(object@pipeline, value)
  prof <- object@env$profile
  for(proto in value) {
    prof <- perform(proto, prof)
  }
  object@env$profile <- prof
  object
})
setMethod("addFiltProfProtos", c("xcmsRaw", "xcmsProtoFilterProfile"), function(object, value)
{
  addFiltProfProtos(object, list(value))
})

setGeneric("filterProfile", function(object, ...) standardGeneric("filterProfile"))
setMethod("filterProfile", "xcmsRaw",
  function(object, method=xcmsMethodDefault("filterProfile"), subtract = FALSE, ...)
{
  proto <- xcmsProtocol("filterProfile", method, ...)
  if (subtract)
    proto <- xcmsProtocol("filterProfile", "subtract", filter = proto)
  addFiltProfProtos(object, proto)
})

# FIXME: Should be modified to accept mz/scan ranges
.filterProfile.median <- function(object, mzrad = 0, scanrad = 0, ...) {
  medianFilter(object, mzrad, scanrad)
}

setGeneric("filterProfile.median", function(object, ...)
  standardGeneric("filterProfile.median"))

setMethod("filterProfile.median", "matrix", .filterProfile.median)

setProtocolClass("xcmsProtoFilterProfileMedian",
  representation(mzrad = "numeric", scanrad = "numeric"),
  c(formals(.filterProfile.median), dispname = "Median"),
  "xcmsProtoFilterProfile")

setMethod("profMargins", "xcmsProtoFilterProfileMedian",
  function(object, mzcount, scancount)
{
  round(c(mzmargin = object@mzrad, scanmargin = object@scanrad))
})

# Pipeline

setGeneric("pipeline", function(object) standardGeneric("pipeline"))
setMethod("pipeline", "xcmsRaw", function(object) object@pipeline)

# set (and run) a new pipeline
setGeneric("pipeline<-", function(object, value) standardGeneric("pipeline<-"))
setReplaceMethod("pipeline", "xcmsRaw", function(object, value) {
    object@pipeline <- value
    # only profile matrix is managed by pipeline currently
    updateProfile(object)
})

setMethod("genProfProto", "xcmsRaw", function(object) genProfProto(object@pipeline))

setReplaceMethod("genProfProto", "xcmsRaw", function(object, value) {
    genProfProto(object@pipeline) <- value
    updateProfile(object)
})

# generate a (full) profile matrix and store in this object
setGeneric("genProfile", function(object, method = "", ...) standardGeneric("genProfile"))
setMethod("genProfile", "xcmsRaw", function(object, method = "", ...) {
    genProfProto(object) <- xcmsProtocol("GenProfile", method, ...)
    object
})

# FIXME: Need to follow 'type.method' delegation pattern of other protocols
setMethod("perform", c("xcmsProtoGenProfile", "xcmsRaw"), function(object, data, ...)
{
  performProfile(object, data@env$mz, data@env$intensity, data@scanindex, data@scantime, ...)
})

setMethod("explore", "xcmsRaw", function(object, ...)
{
  pipeline <- object@pipeline
  if (length(filtProfProtos(pipeline)))
    explore(tail(filtProfProtos(pipeline), 1)[[1]], object, ...)
  else if (length(genProfProto(pipeline)) && !is.null(object@env$profile))
    explore(genProfProto(pipeline), object, ...)
  else { # explore raw data - should raw data input be a protocol?

  }
})

# Private methods (do not export)

# update the profile by performing the appropriate protocols in the pipeline
setGeneric("updateProfile", function(object) standardGeneric("updateProfile"))
setMethod("updateProfile", "xcmsRaw", function(object) {
    # remove old profile matrix
    if ("profile" %in% ls(object@env))
        rm("profile", envir = object@env)
    # store new profile matrix if non-NULL
    profile <- profMat(object, object@pipeline)
    if (!is.null(profile)) {
      object@env$profile <- profile
      object@mzrange <- range(attr(profile, "profmz"))
    }
    object # not needed, but only we know that
})
