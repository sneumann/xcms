### THE PROFILE MATRIX

# Profile pipeline
setClass("xcmsPipelineProfile", contains = "xcmsPipeline")

# Profile matrix data class
setClass("xcmsProfile",
         representation(mzrange = "numeric", step = "numeric",
                        scantime = "numeric"),
         contains = c("xcmsData", "matrix"))

setGeneric("mzBreaks", function(object, ...) standardGeneric("mzBreaks"))
setMethod("mzBreaks", "xcmsProfile", function(object) {
  
  object@mzrange[1]+object@step*(0:(dim(object)[1]-1))
})

setMethod("show", "xcmsProfile", function(object) {

  breaks <- mzBreaks(object)
  cat(object@step, " m/z step (", length(breaks), " grid points from ",
      paste(object@mzrange, collapse = " to "), " m/z)\n", sep = "")
})

setMethod("pipeline", "xcmsProfile",
          function(object, ancestry = TRUE, local = TRUE)
          {
            pipeline <- callNextMethod(object, ancestry, local)
            if (!ancestry)
              pipeline <- new("xcmsPipelineProfile", pipeline)
            pipeline
          })

# analogous to profRange() on xcmsRaw

setGeneric("selectRange", function(object, ...) standardGeneric("selectRange"))

setMethod("selectRange", "xcmsProfile", function(object,
                                                 mzrange = numeric(),
                                                 rtrange = numeric(),
                                                 scanrange = numeric(), ...) {

    if (length(object)) {
        contmass <- mzBreaks(object)
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
            mzrange <- object@mzrange
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
            scanrange <- c(1, length(object@scantime))
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

image.xcmsProfile <- function(x, col = rainbow(256), ...) {

    sel <- selectRange(x, ...)

    zlim <- log(range(x))

    method <- method(profileMatrixProto(x@pipeline))
    title <- paste("XC/MS Log Intensity Image (Profile Method: ",
                   method, ")", sep = "")
    if (zlim[1] < 0) {
        zlim <- log(exp(zlim)+1)
        image(mzBreaks(x)[sel$mzindex], x@scantime[sel$scanidx],
              log(x[sel$mzindex, sel$scanidx]+1),
              col = col, zlim = zlim, main = title, xlab="m/z", ylab="Seconds")
    } else
        image(mzBreaks(x)[sel$mzindex], x@scantime[sel$scanidx],
              log(x[sel$mzindex, sel$scanidx]),
              col = col, zlim = zlim, main = title, xlab="m/z", ylab="Seconds")
}

setGeneric("plotSpec", function(object, ...) standardGeneric("plotSpec"))

setMethod("plotSpec", "xcmsProfile", function(object, ident = FALSE,
                                              vline = numeric(0), ...) {

    sel <- selectRange(object, ...)

    title = paste("Averaged Mass Spectrum: ", sel$rtlabel, " (",
                  sel$scanlab, ")",  sep = "")
    points <- cbind(mzBreaks(object)[sel$mzindex],
                    rowMeans(object[sel$mzindex, sel$scanidx, drop=FALSE]))
    plot(points, type="l", main = title, xlab="m/z", ylab="Intensity")
    if (length(vline))
        abline(v = vline, col = "red")

    if (ident)
        return(identify(points, labels = round(points[,1], 1)))

    invisible(points)
})

setGeneric("plotChrom", function(object, ...) standardGeneric("plotChrom"))

setMethod("plotChrom", "xcmsProfile", function(object, base = FALSE,
                                               ident = FALSE, fitgauss = FALSE,
                                               vline = numeric(0), ...) {

    sel <- selectRange(object, ...)

    if (base) {
        title = paste("Base Peak Chromatogram: ", sel$mzlabel, sep = "")
        pts <- cbind(object@scantime[sel$scanidx],
                     colMax(object[sel$mzindex, sel$scanidx, drop=FALSE]))
    }
    else {
        title = paste("Averaged Ion Chromatogram: ", sel$mzlabel, sep = "")
        pts <- cbind(object@scantime[sel$scanidx],
                     colMeans(object[sel$mzindex, sel$scanidx, drop=FALSE]))
    }
    plot(pts, type="l", main = title, xlab="Seconds", ylab="Intensity")
    if (length(vline))
        abline(v = vline, col = "red")

    if (fitgauss) {
        fit <- nls(y ~ SSgauss(x, mu, sigma, h),
                   data.frame(x = pts[,1], y = pts[,2]))
        points(pts[,1], fitted(fit), type = "l", col = "red", lwd = 2)
        return(fit)
    }

    if (ident)
        return(identify(pts, labels = round(pts[,1], 1)))

    invisible(pts)
})

setGeneric("plotSurf", function(object, ...) standardGeneric("plotSurf"))

setMethod("plotSurf", "xcmsProfile", function(object, log = FALSE,
                                           aspect = c(1, 1, .5), ...) {

  require(rgl) || stop("Couldn't load package rgl")

    sel <- selectRange(object, ...)

    y <- object[sel$mzindex, sel$scanidx]
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

# Profile filtering

# The filterProfile stage is split up into baseline subtraction and smoothing
# This avoids the need for a 'subtract' parameter in a base class and
# more importantly improves clarity.
# Fitting procedures can be shared behind the scenes.

setStage("removeProfileBaseline", "Remove baseline from profile matrix",
         "xcmsProfile")

setStage("smoothProfile", "Fit a smoother to the profile matrix",
         "xcmsProfile")

# Provide necessary margins (as list) in profile matrix for given ranges
# This is to avoid edge effects when processing subsets of the matrix
# FIXME: this is not yet used or exported
setGeneric("profMargins", function(object, ...) standardGeneric("profMargins"))

# Profile filtering

# FIXME: Should be modified to accept mz/scan ranges
.filterProfile.median <- function(object, mzrad = 0, scanrad = 0) {
  object@.Data <- object - medianFilter(object, mzrad, scanrad)
  object
}

setProtocol("median", "Median",
            representation(mzrad = "numeric", scanrad = "numeric"),
            .filterProfile.median, "removeProfileBaseline")

setMethod("profMargins", protocolClass("removeProfileBaseline", "median"),
  function(object, mzcount, scancount)
{
  round(c(mzmargin = object@mzrad, scanmargin = object@scanrad))
})
