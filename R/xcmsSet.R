setClass("xcmsSet", representation(peaks = "matrix", groups = "matrix",
                                   groupidx = "list", sampnames = "character",
                                   sampclass = "factor", rt = "list", 
                                   cdfpaths = "character", profinfo = "list"),
         prototype(peaks = matrix(nrow = 0, ncol = 0), 
                   groups = matrix(nrow = 0, ncol = 0), 
                   groupidx = list(), sampnames = character(0), 
                   sampclass = factor(integer(0)), rt = list(),
                   cdfpaths = character(0), profinfo = vector("list")))

xcmsSet <- function(files = list.files(pattern = ".[Cc][Dd][Ff]$", recursive = TRUE), 
                    snames = gsub("\.[^.]*$", "", basename(files)), 
                    sclass = gsub("^\.$", "sample", dirname(files)),
                    profmethod = "bin", profparam = list(), ...) {

    object <- new("xcmsSet")
    sampnames(object) <- snames
    sampclass(object) <- sclass
    cdfpaths(object) <- file.path(getwd(), files)
    
    rtlist <- list(raw = vector("list", length(snames)),
                   corrected = vector("list", length(snames)))
    
    if ("step" %in% names(list(...)))
        profstep <- list(...)$step
    else
        profstep <- 0.1
    
    profinfo(object) <- c(list(method = profmethod, step = profstep), profparam)
    
    cnames <- c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "into", "intf", "maxo", "maxf")
    
    peaklist <- vector("list", length(files))
    
    for (i in seq(along = peaklist)) {

        lcraw <- xcmsRaw(files[i], profmethod = profmethod, profparam = profparam, 
                         profstep = 0)
        cat(snames[i], ": ", sep = "")
        peaklist[[i]] <- findPeaks(lcraw, ...)[,cnames]
        peaklist[[i]] <- cbind(peaklist[[i]], sample = rep.int(i, nrow(peaklist[[i]])))
        rtlist$raw[[i]] <- lcraw@scantime
        rtlist$corrected[[i]] <- lcraw@scantime
        rm(lcraw)
        gc()
    }
    
    peaks(object) <- do.call("rbind", peaklist)
    object@rt <- rtlist
    
    object
}

if ( !isGeneric("show") )
    setGeneric("show", function(object) standardGeneric("show"))

setMethod("show", "xcmsSet", function(object) {

    cat("An \"xcmsSet\" object with", length(object@sampnames), "samples\n\n")
    
    cat("Time range: ", paste(round(range(object@peaks[,"rt"]), 1), collapse = "-"), 
        " seconds (", paste(round(range(object@peaks[,"rt"])/60, 1), collapse = "-"), 
        " minutes)\n", sep = "")
    cat("Mass range:", paste(round(range(object@peaks[,"mz"], na.rm = TRUE), 4), collapse = "-"), 
        "m/z\n")
    cat("Peaks:", nrow(object@peaks), "(about", 
        round(nrow(object@peaks)/length(object@sampnames)), "per sample)\n")
    cat("Peak Groups:", nrow(object@groups), "\n")
    cat("Sample classes:", paste(levels(object@sampclass), collapse = ", "), "\n\n")
    
    if (length(object@profinfo)) {
        cat("Profile settings: ")
        for (i in seq(along = object@profinfo)) {
            if (i != 1) cat("                  ")
            cat(names(object@profinfo)[i], " = ", object@profinfo[[i]], "\n", sep = "")
        }
        cat("\n")
    }
    
    memsize <- object.size(object)
    cat("Memory usage:", signif(memsize/2^20, 3), "MB\n")
})

c.xcmsSet <- function(...) {

    lcsets <- list(...)
    object <- new("xcmsSet")
    
    peaklist <- vector("list", length(lcsets))
    namelist <- vector("list", length(lcsets))
    classlist <- vector("list", length(lcsets))
    cdflist <- vector("list", length(lcsets))
    rtraw <- vector("list", 0)
    rtcor <- vector("list", 0)
    nsamp <- 0
    for (i in seq(along = lcsets)) {
        peaklist[[i]] <- peaks(lcsets[[i]])
        namelist[[i]] <- sampnames(lcsets[[i]])
        classlist[[i]] <- sampclass(lcsets[[i]])
        classlist[[i]] <- levels(classlist[[i]])[classlist[[i]]]
        cdflist[[i]] <- cdfpaths(lcsets[[i]])
        rtraw <- c(rtraw, lcsets[[i]]@rt$raw)
        rtcor <- c(rtcor, lcsets[[i]]@rt$corrected)
        
        sampidx <- seq(along = namelist[[i]]) + nsamp
        peaklist[[i]][,"sample"] <- sampidx[peaklist[[i]][,"sample"]]
        nsamp <- nsamp + length(namelist[[i]])
    }
    
    peaks(object) <- do.call("rbind", peaklist)
    sampnames(object) <- unlist(namelist)
    classlist <- unlist(classlist)
    sampclass(object) <- factor(classlist, unique(classlist))
    cdfpaths(object) <- unlist(cdflist)
    profinfo(object) <- profinfo(lcsets[[1]])
    object@rt <- list(raw = rtraw, corrected = rtcor)
    
    invisible(object)
}

split.xcmsSet <- function(x, f) {

    f <- factor(f)
    sampidx <- unclass(f)
    peakmat <- peaks(x)
    samples <- sampnames(x)
    classlabel <- sampclass(x)
    cdffiles <- cdfpaths(x)
    prof <- profinfo(x)
    rtraw <- x@rt$raw
    rtcor <- x@rt$corrected
    
    lcsets <- vector("list", length(levels(f)))
    names(lcsets) <- levels(f)
    
    for (i in seq(along = lcsets)) {
        lcsets[[i]] <- new("xcmsSet")
        
        samptrans <- numeric(length(f))
        samptrans[sampidx == i] <- rank(which(sampidx == i))
        samp <- samptrans[peakmat[,"sample"]]
        sidx <- which(samp != 0)
        cpeaks <- peakmat[sidx,]
        cpeaks[,"sample"] <- samp[sidx]
        peaks(lcsets[[i]]) <- cpeaks
        
        sampnames(lcsets[[i]]) <- samples[sampidx == i]
        sampclass(lcsets[[i]]) <- classlabel[sampidx == i, drop = TRUE]
        cdfpaths(lcsets[[i]]) <- cdffiles[sampidx == i]
        profinfo(lcsets[[i]]) <- prof
        lcsets[[i]]@rt$raw <- rtraw[sampidx == i]
        lcsets[[i]]@rt$corrected <- rtcor[sampidx == i]
    }
    
    lcsets
}

if( !isGeneric("peaks") )
    setGeneric("peaks", function(object) standardGeneric("peaks"))

setMethod("peaks", "xcmsSet", function(object) object@peaks)

if( !isGeneric("peaks<-") )
    setGeneric("peaks<-", function(object, value) standardGeneric("peaks<-"))

setReplaceMethod("peaks", "xcmsSet", function(object, value) {

    object@peaks <- value
    
    object
})

if( !isGeneric("groups") )
    setGeneric("groups", function(object) standardGeneric("groups"))

setMethod("groups", "xcmsSet", function(object) object@groups)

if( !isGeneric("groups<-") )
    setGeneric("groups<-", function(object, value) standardGeneric("groups<-"))

setReplaceMethod("groups", "xcmsSet", function(object, value) {

    object@groups <- value
    
    object
})

if( !isGeneric("groupidx") )
    setGeneric("groupidx", function(object) standardGeneric("groupidx"))

setMethod("groupidx", "xcmsSet", function(object) object@groupidx)

if( !isGeneric("groupidx<-") )
    setGeneric("groupidx<-", function(object, value) standardGeneric("groupidx<-"))

setReplaceMethod("groupidx", "xcmsSet", function(object, value) {

    object@groupidx <- value
    
    object
})

if( !isGeneric("sampnames") )
    setGeneric("sampnames", function(object) standardGeneric("sampnames"))

setMethod("sampnames", "xcmsSet", function(object) object@sampnames)

if( !isGeneric("sampnames<-") )
    setGeneric("sampnames<-", function(object, value) standardGeneric("sampnames<-"))

setReplaceMethod("sampnames", "xcmsSet", function(object, value) {

    object@sampnames <- value
    
    object
})

if( !isGeneric("sampclass") )
    setGeneric("sampclass", function(object) standardGeneric("sampclass"))

setMethod("sampclass", "xcmsSet", function(object) object@sampclass)

if( !isGeneric("sampclass<-") )
    setGeneric("sampclass<-", function(object, value) standardGeneric("sampclass<-"))

setReplaceMethod("sampclass", "xcmsSet", function(object, value) {

    if (is.factor(value))
        object@sampclass <- value
    else
        object@sampclass <- factor(value, unique(value))
    
    object
})

if( !isGeneric("cdfpaths") )
    setGeneric("cdfpaths", function(object) standardGeneric("cdfpaths"))

setMethod("cdfpaths", "xcmsSet", function(object) object@cdfpaths)

if( !isGeneric("cdfpaths<-") )
    setGeneric("cdfpaths<-", function(object, value) standardGeneric("cdfpaths<-"))

setReplaceMethod("cdfpaths", "xcmsSet", function(object, value) {

    object@cdfpaths <- value
    
    object
})

if( !isGeneric("profinfo") )
    setGeneric("profinfo", function(object) standardGeneric("profinfo"))

setMethod("profinfo", "xcmsSet", function(object) object@profinfo)

if( !isGeneric("profinfo<-") )
    setGeneric("profinfo<-", function(object, value) standardGeneric("profinfo<-"))

setReplaceMethod("profinfo", "xcmsSet", function(object, value) {

    object@profinfo <- value
    
    object
})

if( !isGeneric("group") )
    setGeneric("group", function(object, ...) standardGeneric("group"))

setMethod("group", "xcmsSet", function(object, bw = 30, minfrac = 0.5, minpeaks = 1,
                                       mzwid = 0.25, max = 5, sleep = 0) {

    samples <- sampnames(object)
    classlabel <- sampclass(object)
    classnames <- levels(classlabel)
    classlabel <- as.vector(unclass(classlabel))
    classnum <- integer(max(classlabel))
    for (i in seq(along = classnum))
        classnum[i] <- sum(classlabel == i)
    
    peakmat <- peaks(object)
    porder <- order(peakmat[,"mz"])
    peakmat <- peakmat[porder,]
    rownames(peakmat) <- NULL
    retrange <- range(peakmat[,"rt"])
    
    minpeakmat <- min(classnum)/2
    
    mass <- seq(peakmat[1,"mz"], peakmat[nrow(peakmat),"mz"] + mzwid, by = mzwid/2)
    masspos <- findEqualGreaterM(peakmat[,"mz"], mass)
    
    groupmat <- matrix(nrow = 512, ncol = 7 + length(classnum))
    groupindex <- vector("list", 512)
    
    endidx <- 0
    num <- 0
    gcount <- integer(length(classnum))
    for (i in seq(length = length(mass)-2)) {
        if (i %% 500 == 0) {
            cat(round(mass[i]), "")
            if (.Platform$OS.type == "windows") flush.console()
        }
        startidx <- masspos[i]
        endidx <- masspos[i+2]-1
        if (endidx - startidx + 1 < minpeakmat)
            next
        speakmat <- peakmat[startidx:endidx,,drop=FALSE]
        den <- density(speakmat[,"rt"], bw, from = retrange[1]-3*bw, to = retrange[2]+3*bw)
        maxden <- max(den$y)
        deny <- den$y
        gmat <- matrix(nrow = 5, ncol = 2+length(classnum))
        snum <- 0
        while (deny[maxy <- which.max(deny)] > maxden/20 && snum < max) {
            grange <- descendMin(deny, maxy)
            deny[grange[1]:grange[2]] <- 0
            gidx <- which(speakmat[,"rt"] >= den$x[grange[1]] & speakmat[,"rt"] <= den$x[grange[2]])
            gnum <- classlabel[unique(speakmat[gidx,"sample"])]
            for (j in seq(along = gcount))
                gcount[j] <- sum(gnum == j)
            if (! any(gcount >= classnum*minfrac & gcount >= minpeaks))
                next
            snum <- snum + 1
            num <- num + 1
            ### Double the size of the output containers if they're full
            if (num > nrow(groupmat)) {
                groupmat <- rbind(groupmat, matrix(nrow = nrow(groupmat), ncol = ncol(groupmat)))
                groupindex <- c(groupindex, vector("list", length(groupindex)))
            }
            groupmat[num, 1] <- median(speakmat[gidx, "mz"])
            groupmat[num, 2:3] <- range(speakmat[gidx, "mz"])
            groupmat[num, 4] <- median(speakmat[gidx, "rt"])
            groupmat[num, 5:6] <- range(speakmat[gidx, "rt"])
            groupmat[num, 7] <- length(gidx)
            groupmat[num, 7+seq(along = gcount)] <- gcount
            groupindex[[num]] <- sort(porder[(startidx:endidx)[gidx]])
        }
        if (sleep > 0) {
            plot(den, main = paste(round(min(speakmat[,"mz"]), 2), "-", round(max(speakmat[,"mz"]), 2)))
            for (i in seq(along = classnum)) {
                idx <- classlabel[speakmat[,"sample"]] == i
                points(speakmat[idx,"rt"], speakmat[idx,"into"]/max(speakmat[,"into"])*maxden, col = i, pch=20)
            }
            for (i in seq(length = snum))
                abline(v = groupmat[num-snum+i, 5:6], lty = "dashed", col = i)
            Sys.sleep(sleep)
        }
    }
    cat("\n")
    
    colnames(groupmat) <- c("mzmed", "mzmin", "mzmax", "rtmed", "rtmin", "rtmax",
                            "npeaks", classnames)
    
    groups(object) <- groupmat[seq(length = num),]
    groupidx(object) <- groupindex[seq(length = num)]
    
    invisible(object)
})

if( !isGeneric("groupval") )
    setGeneric("groupval", function(object, ...) standardGeneric("groupval"))

setMethod("groupval", "xcmsSet", function(object, method = c("medret", "maxint"), 
                                          value = "index", intensity = "into") {

    method <- match.arg(method)
    peakmat <- peaks(object)
    groupmat <- groups(object)
    groupindex <- groupidx(object)
    
    sampnum <- seq(length = length(sampnames(object)))
    retcol <- match("rt", colnames(peakmat))
    intcol <- match(intensity, colnames(peakmat))
    sampcol <- match("sample", colnames(peakmat))
    
    values <- matrix(nrow = length(groupindex), ncol = length(sampnum))
    
    if (method == "medret") {
        for (i in seq(along = groupindex)) {
            gidx <- groupindex[[i]][order(abs(peakmat[groupindex[[i]],retcol] - median(peakmat[groupindex[[i]],retcol])))]
            values[i,] <- gidx[match(sampnum, peakmat[gidx,sampcol])]
        }
    } else {
        for (i in seq(along = groupindex)) {
            gidx <- groupindex[[i]][order(peakmat[groupindex[[i]],intcol], decreasing = TRUE)]
            values[i,] <- gidx[match(sampnum, peakmat[gidx,sampcol])]
        }
    }
    
    if (value != "index") {
        values <- peakmat[values,value]
        dim(values) <- c(length(groupindex), length(sampnum))
    }
    colnames(values) <- sampnames(object)
    rownames(values) <- paste(round(groupmat[,"mzmed"],1), round(groupmat[,"rtmed"]), sep = "/")
    
    values
})

if( !isGeneric("retcor") )
    setGeneric("retcor", function(object, ...) standardGeneric("retcor"))

setMethod("retcor", "xcmsSet", function(object, missing = 1, extra = 1, span = .2, 
                                        family = c("gaussian", "symmetric"),
                                        plottype = c("none", "deviation", "mdevden"), 
                                        col = NULL, ty = NULL) {

    peakmat <- peaks(object)
    groupmat <- groups(object)
    samples <- sampnames(object)
    classlabel <- as.vector(unclass(sampclass(object)))
    n <- length(samples)
    corpeaks <- peakmat
    plottype <- match.arg(plottype)
    family <- match.arg(family)
    if (length(object@rt) == 2)
        rtcor <- object@rt$corrected
    else {
        fnames <- cdfpaths(object)
        rtcor <- vector("list", length(fnames))
        for (i in seq(along = fnames)) {
            cdf <- netCDFOpen(fnames[i])
            rtcor[[i]] <- netCDFVarDouble(cdf, "scan_acquisition_time")
            netCDFClose(cdf)
        }
        object@rt <- list(raw = rtcor, corrected = rtcor)
    }
    
    nsamp <- rowSums(groupmat[,match("npeaks", colnames(groupmat))+unique(classlabel),drop=FALSE])
    
    idx <- which(nsamp >= n-missing & groupmat[,"npeaks"] <= nsamp + extra)
    idx <- idx[clustunique(groupmat[idx,], order(-nsamp[idx], groupmat[idx,"npeaks"]), max=10)]
    idx <- idx[order(groupmat[idx,"rtmed"])]
    
    rt <- groupval(object, "maxint", "rt")[idx,]
    cat("Retention Time Correction Groups:", nrow(rt), "\n")
    rtdev <- rt - apply(rt, 1, median, na.rm = TRUE)
    
    rtdevsmo <- vector("list", n)
    
    for (i in 1:n) {
    
        pts <- na.omit(data.frame(rt = rt[,i], rtdev = rtdev[,i]))
        lo <- suppressWarnings(loess(rtdev ~ rt, pts, span = span, degree = 1, family = family))
        
        rtdevsmo[[i]] <- na.flatfill(predict(lo, data.frame(rt = rtcor[[i]])))
        ### Remove singularities from the loess function
        rtdevsmo[[i]][abs(rtdevsmo[[i]]) > quantile(abs(rtdevsmo[[i]]), 0.9)*2] <- NA
        
        if (length(naidx <- which(is.na(rtdevsmo[[i]]))))
            rtdevsmo[[i]][naidx] <- suppressWarnings(approx(na.omit(data.frame(rtcor[[i]], rtdevsmo[[i]])), 
                                                            xout = rtcor[[i]][naidx])$y)
        while (length(decidx <- which(diff(rtcor[[i]] - rtdevsmo[[i]]) < 0))) {
            d <- diff(rtcor[[i]] - rtdevsmo[[i]])[tail(decidx, 1)]
            rtdevsmo[[i]][tail(decidx, 1)] <- rtdevsmo[[i]][tail(decidx, 1)] - d
        }
    }
    
    if (plottype == "mdevden") {
        split.screen(matrix(c(0, 1, .3, 1, 0, 1, 0, .3), ncol = 4, byrow = TRUE))
        screen(1)
        par(mar = c(0, 4.1, 4.1, 2), xaxt = "n")
    }
    
    if (plottype %in% c("deviation", "mdevden")) {
    
        ### Set up the colors and line type
        if (missing(col)) {
            col <- integer(n)
            for (i in 1:max(classlabel))
                col[classlabel == i] <- 1:sum(classlabel == i)
        }
        if (missing(ty)) {
            ty <- integer(n)
            for (i in 1:max(col))
                ty[col == i] <- 1:sum(col == i)
        }
        if (length(palette()) < max(col))
            mypal <- rainbow(max(col), end = 0.85)
        else
            mypal <- palette()[1:max(col)]
    
        rtrange <- range(do.call("c", rtcor))
        devrange <- range(do.call("c", rtdevsmo))
        
        plot(0, 0, type="n", xlim = rtrange, ylim = devrange, main = "Retention Time Deviation vs. Retention Time", xlab = "Retention Time", ylab = "Retention Time Deviation")
        legend(rtrange[2], devrange[2], samples, col = mypal[col], lty = ty, pch = ceiling(1:n/length(mypal)), xjust = 1)
     
        for (i in 1:n) {
            points(data.frame(rt = rt[,i], rtdev = rtdev[,i]), col = mypal[col[i]], pch = ty[i], type="p")
            points(rtcor[[i]], rtdevsmo[[i]], type="l", col = mypal[col[i]], lty = ty[i])
        }
    }
    
    if (plottype == "mdevden") {
        
        screen(2)
        par(mar = c(5.1, 4.1, 0, 2), yaxt = "n")
        allden <- density(peakmat[,"rt"], bw = diff(rtrange)/200, from = rtrange[1], to = rtrange[2])[c("x","y")]
        corden <- density(rt, bw = diff(rtrange)/200, from = rtrange[1], to = rtrange[2], na.rm = TRUE)[c("x","y")]
        allden$y <- allden$y / sum(allden$y)
        corden$y <- corden$y / sum(corden$y)
        maxden <- max(allden$y, corden$y)
        plot(c(0,0), xlim = rtrange, ylim = c(0, maxden), type = "n", main = "", xlab = "Retention Time", ylab = "Peak Density")
        points(allden, type = "l", col = 1)
        points(corden, type = "l", col = 2)
        abline(h = 0, col = "grey")
        legend(rtrange[2], maxden, c("All", "Correction"), col = 1:2, lty = c(1,1), xjust = 1)
        close.screen(all = TRUE)
    }
    
    for (i in 1:n) {
    
        cfun <- stepfun(rtcor[[i]][-1] - diff(rtcor[[i]])/2, rtcor[[i]] - rtdevsmo[[i]])
        rtcor[[i]] <- rtcor[[i]] - rtdevsmo[[i]]
        
        sidx <- which(corpeaks[,"sample"] == i)
        corpeaks[sidx, c("rt", "rtmin", "rtmax")] <- cfun(corpeaks[sidx, c("rt", "rtmin", "rtmax")])
    }
    
    object@rt$corrected <- rtcor
    peaks(object) <- corpeaks
    groups(object) <- matrix(nrow = 0, ncol = 0)
    groupidx(object) <- list()
    invisible(object)
})

if( !isGeneric("plotrt") )
    setGeneric("plotrt", function(object, ...) standardGeneric("plotrt"))

setMethod("plotrt", "xcmsSet", function(object, col = NULL, ty = NULL, leg = TRUE, densplit = FALSE) {

    samples <- sampnames(object)
    classlabel <- as.vector(unclass(sampclass(object)))
    n <- length(samples)
    rtuncor <- object@rt$raw
    rtcor <- object@rt$corrected
    
    if (missing(col)) {
        col <- integer(n)
        for (i in 1:max(classlabel))
            col[classlabel == i] <- 1:sum(classlabel == i)
    }
    if (missing(ty)) {
        ty <- integer(n)
        for (i in 1:max(col))
            ty[col == i] <- 1:sum(col == i)
    }
    if (length(palette()) < max(col))
        mypal <- rainbow(max(col), end = 0.85)
    else
        mypal <- palette()[1:max(col)]
    
    rtdevsmo <- vector("list", n)
    
    for (i in 1:n)
        rtdevsmo[[i]] <- rtuncor[[i]] - rtcor[[i]]
    
    rtrange <- range(do.call("c", rtuncor))
    devrange <- range(do.call("c", rtdevsmo))
    
    if (densplit) {
        split.screen(matrix(c(0, 1, .3, 1, 0, 1, 0, .3), ncol = 4, byrow = TRUE))
        screen(1)
        par(mar = c(0, 4.1, 4.1, 2), xaxt = "n")
    }
    
    plot(0, 0, type="n", xlim = rtrange, ylim = devrange, main = "Retention Time Deviation vs. Retention Time", xlab = "Retention Time", ylab = "Retention Time Deviation")
    if (leg)
        legend(rtrange[2], devrange[2], samples, col = mypal[col], lty = ty, pch = ceiling(1:n/length(mypal)), xjust = 1)

    for (i in 1:n)
        points(rtuncor[[i]], rtdevsmo[[i]], type="l", col = mypal[col[i]], lty = ty[i])

    if (densplit) {
        
        screen(2)
        par(mar = c(5.1, 4.1, 0, 2), yaxt = "n")
        allden <- density(object@peaks[,"rt"], bw = diff(rtrange)/200, from = rtrange[1], to = rtrange[2])[c("x","y")]
        plot(allden, xlim = rtrange, type = "l", main = "", xlab = "Retention Time", ylab = "Peak Density")
        abline(h = 0, col = "grey")
        close.screen(all = TRUE)
    }
})

if( !isGeneric("fillPeaks") )
    setGeneric("fillPeaks", function(object, ...) standardGeneric("fillPeaks"))

setMethod("fillPeaks", "xcmsSet", function(object, mzrange) {

    peakmat <- peaks(object)
    groupmat <- groups(object)
    files <- cdfpaths(object)
    samp <- sampnames(object)
    classlabel <- as.vector(unclass(sampclass(object)))
    prof <- profinfo(object)
    rtcor <- object@rt$corrected
    nsamp <- rowSums(groupmat[,match("npeaks", colnames(groupmat))+unique(classlabel),drop=FALSE])
    
    gunique <- clustunique(groupmat, order(-nsamp, groupmat[,"npeaks"]))
    groupmat <- groupmat[gunique,]
    groupindex <- groupidx(object)[gunique]
    gvals <- groupval(object)[gunique,]
    
    peakrange <- matrix(nrow = nrow(gvals), ncol = 4)
    colnames(peakrange) <- c("mzmin","mzmax","rtmin","rtmax")
    
    mzmin <- peakmat[gvals,"mzmin"]
    dim(mzmin) <- c(nrow(gvals), ncol(gvals))
    peakrange[,"mzmin"] <- apply(mzmin, 1, min, na.rm = TRUE)
    mzmax <- peakmat[gvals,"mzmax"]
    dim(mzmax) <- c(nrow(gvals), ncol(gvals))
    peakrange[,"mzmax"] <- apply(mzmax, 1, min, na.rm = TRUE)
    retmin <- peakmat[gvals,"rtmin"]
    dim(retmin) <- c(nrow(gvals), ncol(gvals))
    peakrange[,"rtmin"] <- apply(retmin, 1, median, na.rm = TRUE)
    retmax <- peakmat[gvals,"rtmax"]
    dim(retmax) <- c(nrow(gvals), ncol(gvals))
    peakrange[,"rtmax"] <- apply(retmax, 1, median, na.rm = TRUE)
    
    lastpeak <- nrow(peakmat)
    peakmat <- rbind(peakmat, matrix(nrow = sum(is.na(gvals)), ncol = ncol(peakmat)))
    
    cnames <- c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "into", "intf", "maxo", "maxf")
    
    for (i in seq(along = files)) {
        
        cat(samp[i], "")
        naidx <- which(is.na(gvals[,i]))
        lcraw <- xcmsRaw(files[i], profmethod = prof$method, profstep = 0)
        if (length(prof) > 2)
            lcraw@profparam <- prof[seq(3, length(prof))]
        if (length(rtcor) == length(files))
            lcraw@scantime <- rtcor[[i]]
        newpeaks <- getPeaks(lcraw, peakrange[naidx,], step = prof$step)
        rm(lcraw)
        gc()
        newpeaks <- cbind(newpeaks[,cnames], sample = rep(i, length(naidx)))
        peakmat[lastpeak+seq(along = naidx),] <- newpeaks
        for (i in seq(along = naidx))
            groupindex[[naidx[i]]] <- c(groupindex[[naidx[i]]], lastpeak+i)
        lastpeak <- lastpeak + length(naidx)
    }
    cat("\n")
    
    peaks(object) <- peakmat
    groups(object) <- groupmat
    groupidx(object) <- groupindex
    
    invisible(object)
})

if( !isGeneric("getEIC") )
    setGeneric("getEIC", function(object, ...) standardGeneric("getEIC"))

setMethod("getEIC", "xcmsSet", function(object, mzrange, 
                                        sampidx = seq(along = sampnames(object))) {

    files <- cdfpaths(object)
    samp <- sampnames(object)
    prof <- profinfo(object)
    
    eics <- vector("list", length(sampidx))
    
    for (i in seq(along = sampidx)) {
        
        cat(samp[sampidx[i]], "")
        lcraw <- xcmsRaw(files[sampidx[i]], profmethod = prof$method, profstep = 0)
        if (length(prof) > 2)
            lcraw@profparam <- prof[seq(3, length(prof))]
        eics[[i]] <- getEIC(lcraw, mzrange, step = prof$step)
        rm(lcraw)
        gc()
    }
    cat("\n")
    
    invisible(eics)
})

if( !isGeneric("plotEIC") )
    setGeneric("plotEIC", function(object, ...) standardGeneric("plotEIC"))

setMethod("plotEIC", "xcmsSet", function(object, eics, peakrange, nums = seq(length = nrow(peakrange)), 
                                         classlabel = sampclass(object), peakindex = NULL, filebase = character(), 
                                         wh = c(640,480), sleep = 0) {

    classnames <- levels(classlabel)
    classlabel <- as.vector(unclass(classlabel))
    profstep <- profinfo(object)$step
    peakmat <- peaks(object)
    rtcor <- object@rt$corrected
    cols <- palette()[seq(length = max(classlabel))]
    lcols <- cols
    if (!is.null(peakindex)) {
        pcols <- cols
        for (i in seq(along = pcols)) {
            rgbvec <- pmin(col2rgb(cols[i])+153,255)
            cols[i] <- rgb(rgbvec[1], rgbvec[2], rgbvec[3], max = 255)
        }
    }
    
    for (i in seq(along = nums)) {
        pts <- vector("list", length(eics))
        maxint <- numeric(length(eics))
        for (j in seq(along = pts)) {
             stime <- rtcor[[j]]
             idx <- which(stime >= peakrange[nums[i],"rtmin"] & stime <= peakrange[nums[i],"rtmax"])
             pts[[j]] <- cbind(rtcor[[j]][idx], eics[[j]][nums[i],idx])
             maxint[j] <- max(pts[[j]][,2])
        }
        
        mzmin <- floor(peakrange[nums[i],"mzmin"]/profstep)*profstep
        mzmax <- ceiling(peakrange[nums[i],"mzmax"]/profstep)*profstep
        main <- paste("Extracted Ion Chromatogram:", mzmin, "-", mzmax, "m/z")
        if (!missing(filebase))
            png(paste(filebase, sprintf("%03i", as.integer(i)), ".png", sep = ""), w = wh[1], h = wh[2])
        plot(0, 0, type = "n", xlim = peakrange[nums[i],c("rtmin","rtmax")], ylim = c(0, max(maxint, na.rm = TRUE)), main = main,
             xlab = "Retention Time", ylab = "Intensity")
        for (j in order(-maxint)) {
            points(pts[[j]], type = "l", col = cols[classlabel[j]])
            if (!is.null(peakindex) && !is.na(peakindex[nums[i],j])) {
                retrange <- peakmat[peakindex[nums[i],j],c("rtmin","rtmax")]
                ptidx <- which(pts[[j]][,1] >= retrange[1] & pts[[j]][,1] <= retrange[2])
                points(pts[[j]][ptidx,], type = "l", col = pcols[classlabel[j]])
            }
        }
        legend(peakrange[nums[i],c("rtmax")], max(maxint), classnames[unique(classlabel)], col = lcols[unique(classlabel)], lty = 1, xjust = 1)
        if (!missing(filebase))
            dev.off()
        else if (sleep > 0)
            Sys.sleep(sleep)
    }
})

if( !isGeneric("diffreport") )
    setGeneric("diffreport", function(object, ...) standardGeneric("diffreport"))

setMethod("diffreport", "xcmsSet", function(object, class1 = levels(sampclass(object))[1], 
                                            class2 = levels(sampclass(object))[2],
                                            filebase = character(), eicmax = 0, 
                                            sortpval = TRUE, classeic = c(class1,class2)) {
    
    require(multtest) || stop("Couldn't load multtest")
    
    groupmat <- groups(object)
    samples <- sampnames(object)
    n <- length(samples)
    classlabel <- sampclass(object)
    classlabel <- levels(classlabel)[as.vector(unclass(classlabel))]
    
    values <- groupval(object, "medret", "into")
    indecies <- groupval(object, "medret")
    
    if (!all(c(class1,class2) %in% classlabel))
        stop("Incorrect Class Labels")
    
    c1 <- which(classlabel %in% class1)
    c2 <- which(classlabel %in% class2)
    ceic <- which(classlabel %in% classeic)
    if (length(intersect(c1, c2)) > 0)
        stop("Intersecting Classes")
    
    mean1 <- rowMeans(values[,c1], na.rm = TRUE)
    mean2 <- rowMeans(values[,c2], na.rm = TRUE)
    fold <- mean2 / mean1
    fold[!is.na(fold) & fold < 1] <- 1/fold[!is.na(fold) & fold < 1]
    
    testval <- values[,c(c1,c2)]
	testclab <- c(rep(0,length(c1)),rep(1,length(c2)))
	tstat <- mt.teststat(testval, testclab)
	pvalue <- pval(testval, testclab, tstat)
	stat <- data.frame(fold = fold, tstat = tstat, pvalue = pvalue)
	twosamp <- cbind(stat, groupmat, values)
	if (sortpval) {
	   tsidx <- order(twosamp[,"pvalue"])
	   tsidx <- tsidx[clustunique(twosamp[tsidx,c("mzmin", "mzmax", "rtmin", "rtmax")])]
	   twosamp <- twosamp[tsidx,]
	   rownames(twosamp) <- 1:nrow(twosamp)
	}
    
    if (length(filebase)) {
        write.table(twosamp, paste(filebase, ".tsv", sep = ""), quote = FALSE, sep = "\t", col.names = NA)
        if (eicmax > 0) {
            tsnum <- min(eicmax, length(tsidx))
            eicidx <- tsidx[seq(length=tsnum)]
            peakrange <- retexp(groupmat[eicidx, c("mzmin", "mzmax", "rtmin", "rtmax")])
            eics <- getEIC(object, peakrange, ceic)
            dir.create(paste(filebase, "_eic", sep=""))
            plotEIC(object, eics, peakrange, seq(length=tsnum), sampclass(object)[ceic],
                    indecies[eicidx,], file.path(paste(filebase, "_eic", sep=""), ""))
        }
    }
    
    invisible(twosamp)
})

retexp <- function(peakrange, width = 200) {

    retmean <- rowMeans(peakrange[,c("rtmin", "rtmax")])
    peakrange[,"rtmin"] <- retmean-width/2
    peakrange[,"rtmax"] <- retmean+width/2
    
    peakrange
}

#### Create a unique cluster list

clustunique <- function(clust, priority = seq(length = nrow(clust)), mzdiff = 0, max = 5) {

    use <- logical(nrow(clust))
    corder <- order(clust[,"mzmin"])
    clust <- clust[corder, c("mzmin","mzmax","rtmin","rtmax")]
    priority <- priority[corder]
    
    #mztest <- seq(min(clust[,"mzmin"]), max(clust[,"mzmax"]), by = 0.1)
    #show(range(mztest))
    #mintest <- findEqualGreaterM(clust[,"mzmin"], mztest)
    #maxtest <- findEqualGreaterM(clust[,"mzmax"], mztest)
    #show(max(maxtest[(mzdiff/.1):length(maxtest)]-mintest[1:(length(mintest)-mzdiff/.1+1)]))
    #show(max(clust[,2] - clust[,1]))
    
    use <- logical(nrow(clust))
    for (i in priority) {
        nearidx <- c((i-max):(i-1), (i+1):(i+max))
        nearidx <- nearidx[nearidx >= 1 & nearidx <= nrow(clust)]
        nearclust <- clust[nearidx,]
        if (!any(use[nearidx] & !(clust[i,"mzmin"] - nearclust[,"mzmax"] > mzdiff | nearclust[,"mzmin"] - clust[i,"mzmax"] > mzdiff) &
                 !(clust[i,"rtmin"] > nearclust[,"rtmax"] | clust[i,"rtmax"] < nearclust[,"rtmin"])))
            use[i] = TRUE
    }
    
    sort(corder[use])
}

#### Fill in NA values with the minimum and maximum value

na.flatfill <- function(x) {

    realloc <- which(!is.na(x))
    if (realloc[1] > 1)
        x[1:(realloc[1]-1)] <- x[realloc[1]]
    if (realloc[length(realloc)] < length(x))
        x[(realloc[length(realloc)]+1):length(x)] <- x[realloc[length(realloc)]]
    
    x
}

pval <- function(X, classlabel, teststat) {

    n1 <- rowSums(!is.na(X[,classlabel == 0]))
    n2 <- rowSums(!is.na(X[,classlabel == 1]))
    A <- sd(t(X[,classlabel == 0]), na.rm = TRUE)^2/n1
    B <- sd(t(X[,classlabel == 1]), na.rm = TRUE)^2/n2
    df <- (A+B)^2/(A^2/(n1-1)+B^2/(n2-1))
    
    pvalue <- 2 * (1 - pt(abs(teststat), df))
    invisible(pvalue)
}

panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y, use = "pairwise.complete.obs"))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex)
}
