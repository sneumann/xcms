 
xcmsRaw <- function(filename, profstep = 1, profmethod = "intlin",
                    profparam = list(),
                    includeMSn = FALSE, mslevel=NULL) {

    object <- new("xcmsRaw")
    object@env <- new.env(parent=.GlobalEnv)

    if (!file.exists(filename)) stop("File ",filename, " not exists. \n"   )
    if (netCDFIsFile(filename)) {
        if (includeMSn) {
            warning("Reading of MSn spectra for NetCDF not supported")
        }

        cdf <- netCDFOpen(filename)
        if (!is.null(attr(cdf, "errortext")))
            stop(attr(cdf, "errortext"))
        on.exit(netCDFClose(cdf))
        rawdata <- netCDFRawData(cdf)
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
    if (any(rtdiff == 0))
       warning("There are identical scantimes.")

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
        object@msnCollisionEnergy <- rawdataMSn$collisionEnergy
    }

    if (!missing(mslevel) & !is.null(mslevel)) {
      object <- msn2ms(object)
      object <- split(object, f=object@msnLevel==mslevel)$"TRUE"
      ## fix xcmsRaw metadata, or always calculate later than here ?
    }
    
    return(object)
}

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

setGeneric("write.cdf", function(object, ...) standardGeneric("write.cdf"))

setMethod("write.cdf", "xcmsRaw", function(object, filename) {
    require(ncdf) || stop("Couldn't load package ncvar for NetCDF writing")

    scan_no <- length(object@scanindex)
    point_no <- length(object@env$mz)


    dim32bytes <- dim.def.ncdf("_32_byte_string", "", 1:32, create_dimvar=FALSE)
    dim64bytes <- dim.def.ncdf("_64_byte_string", "", 1:64, create_dimvar=FALSE)
    dimError   <- dim.def.ncdf("error_num",       "", 1:1, create_dimvar=FALSE)
    dimScans   <- dim.def.ncdf("scan_number",     "", 1:scan_no, create_dimvar=FALSE)
    dimPoints  <- dim.def.ncdf("point_number",    "", 1:point_no, create_dimvar=FALSE)

    ## Define netCDF vars
    scan_acquisition_time <- var.def.ncdf("scan_acquisition_time", "", dimScans, -1)
    total_intensity       <- var.def.ncdf("total_intensity", "", dimScans, -1)
    scan_index            <- var.def.ncdf("scan_index", "", dimScans, -1)
    total_intensity       <- var.def.ncdf("total_intensity", "", dimScans, -1)
    mass_values           <- var.def.ncdf("mass_values", "", dimPoints, -1)
    intensity_values      <- var.def.ncdf("intensity_values", "", dimPoints, -1)


    ## Define netCDF definitions

    ms <- create.ncdf(filename, list(scan_acquisition_time,
                                     scan_index, total_intensity,
                                     mass_values, intensity_values))

    ## Add values to netCDF vars
    put.var.ncdf(ms, "scan_acquisition_time", object@scantime)
    put.var.ncdf(ms, "total_intensity", object@tic)
    put.var.ncdf(ms, "scan_index", object@scanindex)
    put.var.ncdf(ms, "mass_values", object@env$mz)
    put.var.ncdf(ms, "intensity_values", object@env$intensity)


    close.ncdf(ms)

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

setGeneric("getMsnScan", function(object, ...) standardGeneric("getMsnScan"))

setMethod("getMsnScan", "xcmsRaw", function(object, scanLevel = 2, ms1Rt = -1, parentMzs = 0,
                       precision=1, userMsnIndex=NULL)
{
    if (scanLevel<1)
    {
        warning("Exit: Do you really want to have a ms ",scanLevel," Scan?")
        return(NULL)
    }

    if (is.null(userMsnIndex)) { ## if the User wants to address the data via xcms@msnScanindex
        nxcms <- new("xcmsRaw"); # creates a new empty xcmsraw-object

        nxcms@scantime <- ms1Rt
        nxcms@env$mz        <- object@env$msnMz[(object@msnScanindex[msn]+1):(object@msnScanindex[msn+1])]
        nxcms@env$intensity <- object@env$msnIntensity[(object@msnScanindex[msn]+1):(object@msnScanindex[msn+1])]

        return(nxcms);
    }

    if (parentMzs[1]==0)
        parentMzs <- rep(0,scanLevel-1)

    ## using a zero-vector if none is given
    wasonlyone=TRUE;
    if (ms1Rt < object@scantime[1]) {
        warning("Exit: ms1Rt is smaller than smallest ms1Rt in the object")
        return(NULL)
    }
    ms1ind <- max(which(object@scantime <= ms1Rt))
    if (scanLevel==1) { # in this case just the ms1schan of this rt will be returned
        nxcms <- new("xcmsRaw"); # creates a new empty xcmsraw-object

        nxcms@scantime <- ms1Rt
        nxcms@env$mz        <- object@env$mz[(object@scanindex[ms1ind]+1):(object@scanindex[ms1ind+1])]
        nxcms@env$intensity <- object@env$intensity[(object@scanindex[ms1ind]+1):(object@scanindex[ms1ind+1])]

        return(nxcms);
    }

    if (is.null(object@env$msnMz)) {
        warning("Exit: There are no MSnScans in this object.")
        return(NULL)
    }

    ##finding the peak in the s1 the user wants to have the msnpath from (searching in the ms2parentRtlist):
    ms2s <- which((object@msnRt >= ms1Rt)  &
                  (object@msnLevel == 2) &
                  (object@msnRt <= object@scantime[ms1ind+1]))
    if (length(ms2s) == 0)
    {
        warning("Exit: There is no ms2scan in this Rt-Range!")
        return(NULL)
    }
    ##cat("1> ",ms2s,"\n")
    if (length(ms2s) > 1)
    {
        if (parentMzs[1] == 0)  # more than one ms2scan aviable but no mzvalues given
            warning("More than one ms2scan available but no mz-parameters given! using first scan")
        wasonlyone=FALSE;
        diffe <- abs(object@msnPrecursorMz[ms2s] - parentMzs[1])
        msn <- ms2s[min(which(diffe == min(diffe)))] # The PArent-Rt of this ms2index ist closest to the wanted value
    } else {
        msn <- ms2s; # there is only one ms2scan in this ms1range so use this
    }
    if ((parentMzs[1] != 0) & (abs(object@msnPrecursorMz[msn] - parentMzs[1]) > 1)) {
        warning("No ms2scan parent m/z is close enought to the requested value! using closest:",object@msnPrecursorMz[msn])
    }
    msnRt <- object@msnRt[msn]
    ##cat("3> ",msnRt,"\n")
    if (scanLevel > 2) {
        for (a in 3:scanLevel) {
            msns <- which((object@msnRt >= msnRt) &
                          (object@msnLevel == a) &
                          (object@msnRt <= object@scantime[ms1ind+1]))
            ##cat("4> ",ms2s,"\n")
            if (length(msns)==0) {
                warning("Exit: There is no ms",a,"scan in this Rt-Range!")
                return(NULL)
            }
            if (length(msns)>1) {
                wasonlyone=FALSE;
                if ((length(parentMzs)< a-1) | (parentMzs[a-1] == 0)) { # more than one ms2scan aviable but no mzvalues given
                    warning("More than one ms",a,"scan available but no mzdata given! using first scan")
                    msn <- msns[1];
                } else {
                    diffe <- abs(object@msnPrecursorMz[msns] - parentMzs[a-1])
                    msn <- msns[min(which(diffe == min(diffe)))]
                }
            } else {
                msn <- msns; # there is only one ms[n-1]scan in this ms[n]ramge so use this
            }
            if (length(parentMzs)>=(a-1)&(parentMzs[1]!=0)) {
                if (abs(object@msnPrecursorMz[msn] - parentMzs[a-1]) > 1) {
                    warning("No ms",scanLevel,"scan parent m/z is close enought to the requested value! using closest: ", object@msnPrecursorMz[msn])
                }
            }
            msnRt <- object@msnRt[msn]
        }
    }
    if (wasonlyone) {
        message("Note: There was only one ms",scanLevel,"Scan for the given MS1rt.\n", sep="")
    }
    nxcms <- new("xcmsRaw"); # creates a new empty xcmsraw-object

    nxcms@scantime <- msnRt
    nxcms@env$mz        <- object@env$msnMz[(object@msnScanindex[msn]+1):(object@msnScanindex[msn+1])]
    nxcms@env$intensity <- object@env$msnIntensity[(object@msnScanindex[msn]+1):(object@msnScanindex[msn+1])]

    return(nxcms);
})

setGeneric("getSpec", function(object, ...) standardGeneric("getSpec"))

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

specNoise <- function(spec, gap = quantile(diff(spec[,"mz"]), .9)) {

    # In a spectrum with just one raw peak we can't calculate noise
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

setGeneric("plotSpec", function(object, ...) standardGeneric("plotSpec"))

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

setGeneric("plotChrom", function(object, ...) standardGeneric("plotChrom"))

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

setGeneric("image", function(x, ...) standardGeneric("image"))

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

setGeneric("plotSurf", function(object, ...) standardGeneric("plotSurf"))

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

filtfft <- function(y, filt) {

    yfilt <- numeric(length(filt))
    yfilt[1:length(y)] <- y
    yfilt <- fft(fft(yfilt, inverse = TRUE) * filt)

    Re(yfilt[1:length(y)])
}

setClass("xcmsPeaks", contains = "matrix")

setMethod("show", "xcmsPeaks", function(object) {
  cat("A matrix of", nrow(object), "peaks\n")
  cat("Column names:\n")
  print(colnames(object))
})

setGeneric("findPeaks.matchedFilter", function(object, ...) standardGeneric("findPeaks.matchedFilter"))

setMethod("findPeaks.matchedFilter", "xcmsRaw", function(object, fwhm = 30, sigma = fwhm/2.3548,
                                                         max = 5, snthresh = 10, step = 0.1,
                                                         steps = 2, mzdiff = 0.8 - step*steps,
                                                         index = FALSE, sleep = 0,
                                                         verbose.columns = FALSE) {

    profFun <- match.profFun(object)

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
                 mzrange <- range(mzmat, na.rm = TRUE)
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

setGeneric("findPeaks.centWave", function(object, ...) standardGeneric("findPeaks.centWave"))

setMethod("findPeaks.centWave", "xcmsRaw", function(object, ppm=25, peakwidth=c(20,50), snthresh=10,
                                                    prefilter=c(3,100), mzCenterFun="wMean", integrate=1, mzdiff=-0.001, 
                                                    fitgauss=FALSE, scanrange= numeric(), noise=0, # noise.local=TRUE,
                                                    sleep=0, verbose.columns=FALSE) {
    if (!isCentroided(object))
        warning("It looks like this file is in profile mode. centWave can process only centroid mode data !\n")

    mzCenterFun <- paste("mzCenter", mzCenterFun, sep=".")
    if (!exists(mzCenterFun, mode="function"))
        stop("Error: >",mzCenterFun,"< not defined ! \n")

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
    if (length(scalerange) > 1)
        scales <- seq(from=scalerange[1], to=scalerange[2], by=2)  else
            scales <- scalerange;

    dev <- ppm * 1e-6;
    minPeakWidth <-  scales[1];
    noiserange <- c(minPeakWidth*3, max(scales)*3);
    maxGaussOverlap <- 0.5;
    minPtsAboveBaseLine <- max(4,minPeakWidth-2);
    minCentroids <- minPtsAboveBaseLine ;
    scRangeTol <-  maxDescOutlier <- floor(minPeakWidth/2);

    peaklist <- list()
    cat("\n Detecting mass traces at",ppm,"ppm ... \n"); flush.console();
    featlist <- findmzROI(object,scanrange=scanrange,dev=dev,minCentroids=minCentroids, prefilter=prefilter, noise=noise)
    scantime <- object@scantime
    Nscantime <- length(scantime)
    lf <- length(featlist)
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

      feat <- featlist[[f]]
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
      if (any(omz == 0))
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
      noise <- estimateChromNoise(noised,c(0.05,0.95),minPts=3*minPeakWidth)

      ## any continuous data above 1st baseline ?
      if (!continuousPtsAboveThreshold(fd,threshold=noise,num=minPtsAboveBaseLine)) 
        next;
  
      ## 2nd baseline estimate using not-peak-range
#       if (noise.local)  ## experimental
#         lnoise <- getLocalNoiseEstimate(d,td,ftd,noiserange,Nscantime)  else  
#           lnoise <- c(noise, sd(noised))

      lnoise <- getLocalNoiseEstimate(d,td,ftd,noiserange,Nscantime)


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
                    # maxint <- max(d[pprange])
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
          }  #for
      } # if

      
      ##  postprocessing
      if (!is.null(peaks)) {
                colnames(peaks) <- c(basenames, verbosenames)

          colnames(peakinfo) <- c("scale","scaleNr","scpos","scmin","scmax")

          for (p in 1:dim(peaks)[1]) {
            ## find minima, assign rt and intensity values
            if (integrate == 1) {
                lm <- descendMin(wCoefs[,peakinfo[p,"scaleNr"]], istart= peakinfo[p,"scpos"])
                if (lm[1]==lm[2]) ## fall-back
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
          peaks <- joinOverlappingPeaks(td,d,otd,omz,od,scantime,scan.range,peaks,maxGaussOverlap,mzCenterFun=mzCenterFun)
      }



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

    } # f

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


setGeneric("findPeaks.MSW", function(object, ...) standardGeneric("findPeaks.MSW"))

setMethod("findPeaks.MSW", "xcmsRaw", function(object, snthresh=3, verbose.columns = FALSE, ...)
{
  require(MassSpecWavelet) || stop("Couldn't load MassSpecWavelet")

  # MassSpecWavelet Calls
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


  # Assemble result

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

  # Filter additional (verbose) columns
  if (!verbose.columns)
    peaklist <- peaklist[,basenames,drop=FALSE]

  invisible(new("xcmsPeaks", peaklist))
}
)

setGeneric("findPeaks.MS1", function(object, ...) standardGeneric("findPeaks.MS1"))

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



setGeneric("findPeaks", function(object, ...) standardGeneric("findPeaks"))

setMethod("findPeaks", "xcmsRaw", function(object, method=getOption("BioC")$xcms$findPeaks.method,
                                           ...) {

    method <- match.arg(method, getOption("BioC")$xcms$findPeaks.methods)
    if (is.na(method))
        stop("unknown method : ", method)
    method <- paste("findPeaks", method, sep=".")
    invisible(do.call(method, list(object, ...)))
})

setGeneric("getPeaks", function(object, ...) standardGeneric("getPeaks"))

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
        rmat[i,1] <- weighted.mean(mass[imz[1]:imz[2]], rowSums(ymat))
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
  if (log)
    y <- log(y + max(1 - min(y), 0))

  cbind(time = object@scantime[scans[massidx]], mz = masses[massidx],
        intensity = y)
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

setGeneric("profMethod", function(object) standardGeneric("profMethod"))

setMethod("profMethod", "xcmsRaw", function(object) {

    object@profmethod
})

.profFunctions <- list(intlin = "profIntLinM", binlin = "profBinLinM",
                       binlinbase = "profBinLinBaseM", bin = "profBinM")

setGeneric("profMethod<-", function(object, value) standardGeneric("profMethod<-"))

setReplaceMethod("profMethod", "xcmsRaw", function(object, value) {

    if (! (value %in% names(.profFunctions)))
        stop("Invalid profile method")

    object@profmethod <- value

    profStep(object) <- profStep(object)

    object
})

setGeneric("profStep", function(object) standardGeneric("profStep"))

setMethod("profStep", "xcmsRaw", function(object) {

    if (is.null(object@env$profile))
        0
    else
        diff(object@mzrange)/(nrow(object@env$profile)-1)
})

setGeneric("profStep<-", function(object, value) standardGeneric("profStep<-"))

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

setGeneric("rawEIC", function(object, ...) standardGeneric("rawEIC"))

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

setGeneric("plotEIC", function(object, ...) standardGeneric("plotEIC"))

setMethod("plotEIC", "xcmsRaw", function(object,
                                           mzrange = numeric(),
                                           rtrange = numeric(),
                                           scanrange = numeric())  {
   
   EIC <-  rawEIC(object,mzrange=mzrange, rtrange=rtrange, scanrange=scanrange)
   points <- cbind(object@scantime[EIC$scan], EIC$intensity)

   plot(points, type="l", main=paste("Extracted Ion Chromatogram  m/z  ",mzrange[1]," - ",mzrange[2],sep=""), xlab="Seconds",
        ylab="Intensity")

   invisible(points)
})


setGeneric("rawMZ", function(object, ...) standardGeneric("rawMZ"))

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


setGeneric("findmzROI", function(object, ...) standardGeneric("findmzROI"))

setMethod("findmzROI", "xcmsRaw", function(object, mzrange=c(0.0,0.0), scanrange=c(1,length(object@scantime)),dev, minCentroids, prefilter=c(0,0), noise=0){

  scanrange[1] <- max(1,scanrange[1])
  scanrange[2] <- min(length(object@scantime),scanrange[2])

  ## mzrange not implemented yet
  if (!is.double(object@env$mz))  object@env$mz <- as.double(object@env$mz)
  if (!is.double(object@env$intensity)) object@env$intensity <- as.double(object@env$intensity)
  if (!is.integer(object@scanindex)) object@scanindex <- as.integer(object@scanindex)

  .Call("findmzROI", object@env$mz,object@env$intensity,object@scanindex, as.double(mzrange),
  as.integer(scanrange), as.integer(length(object@scantime)),
  as.double(dev), as.integer(minCentroids), as.integer(prefilter), as.integer(noise), PACKAGE ='xcms' )
})


setGeneric("isCentroided", function(object, ...) standardGeneric("isCentroided"))

setMethod("isCentroided", "xcmsRaw", function(object){
    if (length(getScan(object,length(object@scantime) / 2)) >2 ) {
        quantile(diff(getScan(object,length(object@scantime) / 2)[,"mz"]),.25)  > 0.025
    } else {
        TRUE
    }
})

sequenceMz <- function(dat) {
    for (p in 1:dim(dat)[1] ){ # makes the index for the scan
        seq<-seq(from=dat[p,"from"], to=dat[p,"to"])
        if(p == 1){
            seqInd<-seq
        } else {
            seqInd<-c(seqInd, seq)
        }
    }
    return(seqInd)
}

if (!isGeneric("collect") )
    setGeneric("collect", function(object, ...) standardGeneric("collect"))

setMethod("collect", "xcmsRaw", function(object, rtU, mzU=0, sn=5, uCE=-1, check=FALSE, fragments=TRUE, ...) {
    for(k in 1:length(object@msnScanindex)){
        if (k ==1){
            from<-object@msnScanindex[k]
            to<-object@msnScanindex[k+1]-1
            mat<-c(from, to)
        }else if (k==length(object@msnScanindex)){
            from<-object@msnScanindex[k]
            to<-length(object@env$msnMz)
            mat<-rbind(mat, c(from, to))
        } else{
            from<-object@msnScanindex[k]
            to<-object@msnScanindex[k+1]-1
            mat<-rbind(mat, c(from, to))
        }
    }

    uniMZ<-unique(object@msnPrecursorMz)
    ref<-1

    run<-vector("list", 0)
    for(i in 1:length(uniMZ)){
	cat(paste(uniMZ[i], " ", sep=""))
	if(mzU==0){
	    mzU<-uniMZ[i]
	    check<-TRUE
	}

	uniCE<-unique(object@msnCollisionEnergy[object@msnPrecursorMz == mzU])
	for (o in 1:length(uniCE)){
	    if(uCE == -1){
	        uCE<-uniCE[o]
	        check<-TRUE
	    }
        tempIndex<-which(object@msnPrecursorMz == mzU & object@msnCollisionEnergy == uCE)
	    scanIX<-cbind(object@msnPrecursorMz[tempIndex], mat[tempIndex,2], mat[tempIndex,1], 
						object@msnRt[tempIndex], object@msnCollisionEnergy[tempIndex])
        colnames(scanIX)<-c("preMZ", "to", "from", "rt", "collisionEnergy")
	    if(!is.matrix(scanIX)){
	        tempNames<-names(scanIX)
	        scanIX<-matrix(scanIX, ncol=5)
	        colnames(scanIX)<-tempNames
	    }
	    rthold<-min(scanIX[,"rt"])

            if(dim(scanIX)[1] == 1){ ##do we need rt time window?
                dat<-scanIX
                seqInd<-sequenceMz(dat)
                scan<-matrix(c(object@env$msnMz[seqInd], object@env$msnIntensity[seqInd]), ncol=2)
                colnames(scan)<-c("mz", "intensity")
                scan<-scan[order(scan[,"mz"], scan[,"intensity"]),, drop=FALSE]
                scan[,"intensity"]<-scan[,"intensity"]/max(scan[,"intensity"])*100
                #spectab<-specPeaks(scan, sn=sn) ## used to remove noise
                run[[ref]]<-specPeaks(scan, sn=sn)
                #run[[ref]]<-deisotopeNclean(spectab) ##Just what it says. could be better

                if(!exists("runinfo") ){ # to get over the null object issue
                    runinfo<-c(mzU, min(dat[,"rt"]), max(dat[,"rt"]), ref, uCE)
                } else {
                    runinfo<-rbind(runinfo, c(mzU, min(dat[,"rt"]), max(dat[,"rt"]), ref, uCE))
                }
                ref<-ref+1
            }else if (dim(scanIX)[1] >1){
	        while(rthold <= max(scanIX[,"rt"])){
	            ahead<-scanIX[which.max(scanIX[,"rt"]>rthold),"rt"]+rtU
	            scanIndex<-scanIX[,"rt"] <= ahead & scanIX[, "rt"] > rthold
                    if (!is.matrix(scanIX[scanIndex,])){ # if dim==[1,8] it becomes a vector so check for it
                        dat<-matrix(scanIX[scanIndex,], ncol=5)
                        colnames(dat)<-names(scanIX[scanIndex,])
                    }else {
                        dat<-scanIX[scanIndex,]
                    }
                    if(dim(dat)[1] < 1 ){
                        rthold<-ahead+rthold
                        next
                    }
                    seqInd<-sequenceMz(dat)
                    scan<-matrix(c(object@env$msnMz[seqInd], object@env$msnIntensity[seqInd]), ncol=2)
                    colnames(scan)<-c("mz", "intensity")
                    scan<-scan[order(scan[,"mz"], scan[,"intensity"]),, drop=FALSE]
                    scan[,"intensity"]<-scan[,"intensity"]/max(scan[,"intensity"])*100
##normilise the intensity for signal to noise ration.
                    #spectab<-specPeaks(scan, sn=sn) # used to remove noise
                    run[[ref]]<-specPeaks(scan, sn=sn)
                    #run[[ref]]<-deisotopeNclean(spectab)
                    rthold<-max(dat[,"rt"])

                    if(!exists("runinfo") ){ # to get over the null object issue
                        runinfo<-c(mzU, min(dat[,"rt"]), max(dat[,"rt"]), ref, uCE)
                    } else {
                        runinfo<-rbind(runinfo, c(mzU, min(dat[,"rt"]), max(dat[,"rt"]), ref, uCE))
                    }
                    ref<-ref+1
                }
            }
        }

		if(!exists("runinfoFinal")){
	    	runinfoFinal<-runinfo
		} else {
	    	runinfoFinal<-rbind(runinfoFinal, runinfo)
		}
		if(check==TRUE) {
	    	mzU<-0
	    	uCE<- -1
		}
    }

    colnames(runinfoFinal)<-c("preMZ", "rtmin", "rtmax", "ref", "CollisionEnergy")
    rownames(runinfoFinal)<-rep("", dim(runinfoFinal)[1]) #just to get rid of the added row names which aren't needed.
    cat("\nCalculating accurate mass...")
    if(fragments == TRUE){
        frag<-new("xcmsFragments")
    }
    if(length(run)== dim(runinfoFinal)[1] ){## some odd stuff happens with the matrix duplicates are made
		frag@MS2spec<-run## So we check for it
        frag@specinfo<-runinfoFinal
    } else {
		dup<-duplicated(runinfoFinal)
		if(dim(runinfoFinal[!dup,])[1] == length(run)){
            frag@specinfo<-runinfoFinal[!dup,]
            frag@MS2spec<-run
		} else{
	    	cat(paste("Error: Method \'collect\' has failed. \n", sep="")) ##should never be seen :D just having fun
		}
    }

    cat(paste("\n", sep=""))
	if(length(object@env$mz) >1){
		accurateMZ<-getMZ(object, frag@specinfo)
	}else{
		accurateMZ<-as.numeric(frag@specinfo[,1])
	}
    
    frag@specinfo<-cbind(frag@specinfo[,1], accurateMZ, frag@specinfo[,2:5])
    colnames(frag@specinfo)<-c("preMZ", "AccMZ", "rtmin", "rtmax", "ref", "CollisionEnergy")
    if(fragments == TRUE){
        return(frag)
    } else {
        MSMS<-list()
        MSMS[[1]]<-frag@specinfo
        MSMS[[2]]<-frag@MS2spec
        return(MSMS)
    }
})

if (!isGeneric("getMZ") )
    setGeneric("getMZ", function(object, ...) standardGeneric("getMZ"))

setMethod("getMZ", "xcmsRaw", function(object, specinfo,  ...) {
    for(i in 1:dim(specinfo)[1]){
		A<-which(specinfo[i,"rtmin"] > object@scantime)
		if(specinfo[i,"rtmax"] > object@scantime[length(object@scantime)]){
			B<-length(object@scantime)
		} else {
				B<-which(specinfo[i,"rtmax"] < object@scantime)
		}

		massInd<-object@scanindex[max(A,na.rm=TRUE):min(B,na.rm=TRUE)]
		AccuMass<-object@env$mz[massInd[1]:massInd[length(massInd)]]
		AccuIntensity<-object@env$intensity[massInd[1]:massInd[length(massInd)]]


		index<-which(AccuMass > specinfo[i, "preMZ"]-0.25 & AccuMass < specinfo[i, "preMZ"]+0.25)
		#if (length(index) == 0){ # we should probably check for nothing found :wq

			#}
		if(!exists("accurate")){
			accurate<-weighted.mean(AccuMass[index], AccuIntensity[index], na.rm = TRUE)
		}else{
			accurate<-c(accurate, weighted.mean(AccuMass[index], AccuIntensity[index], na.rm=TRUE))
		}
    }
    return(accurate)
})


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

sequences <- function(seqs) {
   apply(seqs, 1, FUN=function(x) {x[1]:x[2]})
}

match.profFun <- function(object) {
  match.fun(.profFunctions[[profMethod(object)]])
}

setGeneric("msnparent2ms", function(object, ...) standardGeneric("msnparent2ms"))
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

setGeneric("msn2ms", function(object, ...) standardGeneric("msn2ms"))
setMethod("msn2ms", "xcmsRaw", function(object) {

    object@tic <- rep(0, length(object@msnAcquisitionNum)) ## 

    object@scantime <- object@msnRt
    object@acquisitionNum <- object@msnAcquisitionNum
    object@scanindex <- object@msnScanindex

    object@env$mz <- object@env$msnMz
    object@env$intensity <- object@env$msnIntensity
    invisible(object)

})

setGeneric("deepCopy", function(object) standardGeneric("deepCopy"))
setMethod("deepCopy", "xcmsRaw", function(object) {

    x <- object
    x@env <- new.env(parent=.GlobalEnv)

    for (variable in ls(object@env)) {
      eval(parse(text=paste("x@env$",variable," <- object@env$",variable,sep="")))
    }
    
    invisible(x)
})


  
