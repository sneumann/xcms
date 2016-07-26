## Functions for xcmsEIC
#' @include DataClasses.R

plot.xcmsEIC <- function(x, y, groupidx = groupnames(x), sampleidx = sampnames(x),
                         rtrange = x@rtrange, col = rep(1, length(sampleidx)), legtext = NULL,
                         peakint = TRUE, sleep = 0, mzdec=2, ...) {

    object <- x

    if (is.numeric(groupidx))
        groupidx <- object@groupnames[groupidx]

    if (length(object@groupnames))
        grpidx <- match(groupidx, object@groupnames)
    else
        grpidx <- 1

    if (is.numeric(sampleidx))
        sampleidx <- names(object@eic)[sampleidx]
    sampidx <- match(sampleidx, names(object@eic))

    if (length(rtrange) == 1)
        rtrange <- retexp(object@rtrange, rtrange)
    else if (length(rtrange) < length(object@rtrange)) {
        tmprtrange <- matrix(rtrange, ncol = 2)
        rtrange <- matrix(nrow = nrow(object@rtrange), ncol = 2)
        rtrange[sampidx,] <- matrix(rep(t(tmprtrange),
                                        length = length(sampidx)*2), ncol = 2,
                                    byrow = TRUE)
    }

    if (!missing(col) && length(col) < length(object@eic)) {
        tmpcol <- col
        col <- vector(class(col), length(object@eic))
        col[sampidx] <- rep(tmpcol, length = length(sampidx))
    }

    if (!missing(y)) {
        xset <- y
        pks <- peaks(xset)
        pidx <- groupval(xset)
        xsgrpidx <- match(object@groupnames, groupnames(xset, template = groupidx))
        xssampidx <- match(names(object@eic), sampnames(xset))

        if (missing(col)) {
            col <- unclass(sampclass(xset))[xssampidx]
            if (length(palette()) < max(col))
                col <- rainbow(max(col), end = 0.85)[col]
        }
        lcol <- col
        for (i in seq(along = lcol)) {
            rgbvec <- pmin(col2rgb(lcol[i])+153,255)
            lcol[i] <- rgb(rgbvec[1], rgbvec[2], rgbvec[3],
                           maxColorValue = 255)
        }

        if (missing(legtext))
            legtext <- levels(sampclass(xset))
    } else {
        if (missing(col))
            col <- rep(1, length(object@eic))
    }

    for (i in grpidx) {
        maxint <- numeric(length(sampidx))
        for (j in seq(along = sampidx)) {
            eic.int <- object@eic[[sampidx[j]]][[i]][,"intensity"]
            if (length(eic.int) > 0)
                maxint[j] <- max(eic.int)
            else
                maxint[j] <- 0
        }
        plot(0, 0, type = "n", xlim = rtrange[i,], ylim = c(0, max(maxint)),
             xlab = "Retention Time (seconds)", ylab = "Intensity",
             main = paste("Extracted Ion Chromatogram:", round(object@mzrange[i,1], mzdec),
             "-", round(object@mzrange[i,2], mzdec), "m/z"))
        for (j in sampidx[order(maxint, decreasing = TRUE)]) {
            pts <- object@eic[[j]][[i]]
            if (missing(y) || !peakint)
                points(pts, type = "l", col = col[j])
            else {
                points(pts, type = "l", col = lcol[j])
                peakrange <- pks[pidx[xsgrpidx[i],xssampidx[j]], c("rtmin","rtmax")]
                if (object@rt == "raw") {
                    corrt <- xset@rt$corrected[[xssampidx[j]]]
                    cidx <- which(corrt >= peakrange[1] & corrt <= peakrange[2])
                    peakrange <- xset@rt$raw[[xssampidx[j]]][c(cidx[1],cidx[length(cidx)])]
                }
                ptsidx <- pts[,"rt"] >= peakrange[1] & pts[,"rt"] <= peakrange[2]
                points(pts[ptsidx,], type = "l", col = col[j])
            }
        }

        if (!is.null(legtext))
            if (is.numeric(col))
                legend(rtrange[i,2], max(maxint), legtext[unique(col[sampidx])],
                       col = unique(col[sampidx]), lty = 1, xjust = 1)
            else  ### THIS NEEDS A BIT OF WORK ###
                legend(rtrange[i,2], max(maxint), legtext,
                       col = unique(col[sampidx]), lty = 1, xjust = 1)


        if (sleep > 0)
            Sys.sleep(sleep)
    }
}
