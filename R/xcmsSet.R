# TODO:
# 1) phenoData -> AnnotatedDataFrame

setClass("xcmsSet", representation(peaks = "matrix", groups = "matrix",
                                   groupidx = "list", comps = "matrix",
                                   phenoData = "data.frame", rt = "list",
                                   filepaths = "character"
                                   ),
         prototype(peaks = matrix(nrow = 0, ncol = 0),
                   groups = matrix(nrow = 0, ncol = 0),
                   comps = matrix(nrow = 0, ncol = 0),
                   groupidx = list(),
                   phenoData = data.frame(), rt = list(),
                   filepaths = character(0)),
         "xcmsData")

xcmsSet <- function(files = NULL, snames = NULL, sclass = NULL,
                    profmethod = "bin", profparam = list(),
                    ..., rawpipeline = NULL, pipeline = NULL,
                    phenoData = NULL) {

    object <- new("xcmsSet")

    filepattern <- c("[Cc][Dd][Ff]", "[Nn][Cc]", "([Mm][Zz])?[Xx][Mm][Ll]",
                     "[Mm][Zz][Dd][Aa][Tt][Aa]")
    filepattern <- paste(paste("\\.", filepattern, "$", sep = ""), collapse = "|")
    if (is.null(files))
        files <- getwd()
    info <- file.info(files)
    listed <- list.files(files[info$isdir], pattern = filepattern,
                         recursive = TRUE, full.names = TRUE)
    files <- c(files[!info$isdir], listed)

    # try making paths absolute
    files_abs <- file.path(getwd(), files)
    exists <- file.exists(files_abs)
    files[exists] <- files_abs[exists]
    
    filepaths(object) <- files

    # determine experimental design
    fromPaths <- phenoDataFromPaths(files)
    if (is.null(snames))
      snames <- rownames(fromPaths)
    pdata <- phenoData
    if (is.null(pdata)) {
      pdata <- sclass
      if (is.null(pdata)) 
        pdata <- fromPaths
    }
    phenoData(object) <- pdata
    if (is.null(phenoData))
      rownames(phenoData(object)) <- snames

    # figure out 'step' argument
    mc <- as.list(match.call())
    #step <- mc$step
    #if (is.null(step))
    #  step <- 0.1

    # determine peak arguments
    peakargs <- list(...)
    #peakargs <- vargs[!(vargs %in% "step")]

    specified <- match(c("profmethod", "profparams"), names(mc), 0)
    params <- c(names(mc)[specified], peakargs)
    if (length(params) && (!is.null(pipeline) || !is.null(rawpipeline)))
      warning("Parameter(s) ",
              paste("'", params, "'", collapse = ", "),
              " overriden by the 'pipeline' and 'rawpipeline' parameters")
    if (!is.null(rawpipeline) && !is.null(pipeline))
      warning("Parameter 'rawpipeline' overriden by the 'pipeline' parameter")
    if (is.null(pipeline)) {
      if (is.null(rawpipeline)) {
        # store pipeline parameters in the raw pipeline (but don't generate)
        profargs <- c(list(method = profmethod, step = 0), profparam)
        profmatproto <- do.call("xcmsProtocol", c("profileMatrix", profargs))
        profpipe <- new("xcmsPipelineProfile", list(profmatproto))
        genprofproto <- xcmsProtocol("genProfile", "generic", pipeline=profpipe)
        peaksproto <- do.call("xcmsProtocol", c("findPeaks", peakargs))
        rawpipeline <- new("xcmsPipeline", list(genprofproto, peaksproto))
      }
      rawproto <- xcmsProtocol("processRaws", "generic", pipeline = rawpipeline)
      pipeline <- new("xcmsPipeline", list(rawproto))
    }

    perform(pipeline, object)
}

setMethod("show", "xcmsSet", function(object) {

    cat("An \"xcmsSet\" object with", nrow(object@phenoData), "samples\n\n")

    cat("Time range: ", paste(round(range(object@peaks[,"rt"]), 1), collapse = "-"),
        " seconds (", paste(round(range(object@peaks[,"rt"])/60, 1), collapse = "-"),
        " minutes)\n", sep = "")
    cat("Mass range:", paste(round(range(object@peaks[,"mz"], na.rm = TRUE), 4), collapse = "-"),
        "m/z\n")
    cat("Peaks:", nrow(object@peaks), "(about",
        round(nrow(object@peaks)/nrow(object@phenoData)), "per sample)\n")
    cat("Peak Groups:", nrow(object@groups), "\n")
    cat("Sample classes:", paste(levels(sampclass(object)), collapse = ", "), "\n\n")

    show(object@pipeline)
    cat("\n")

    memsize <- object.size(object)
    cat("Memory usage:", signif(memsize/2^20, 3), "MB\n")
})

# The c() and split() methods should be used only for interpreting results;
# not as part of the preprocessing.

# There is no attempt to ensure consistency with the pipeline.

c.xcmsSet <- function(...) {

    lcsets <- list(...)
    object <- new("xcmsSet")

    peaklist <- vector("list", length(lcsets))
    pdatalist <- vector("list", length(lcsets))
    cdflist <- vector("list", length(lcsets))
    rtraw <- vector("list", 0)
    rtcor <- vector("list", 0)
    nsamp <- 0
    for (i in seq(along = lcsets)) {
        peaklist[[i]] <- peaks(lcsets[[i]])
        pdatalist[[i]] <- phenoData(lcsets[[i]])
        cdflist[[i]] <- filepaths(lcsets[[i]])
        rtraw <- c(rtraw, lcsets[[i]]@rt$raw)
        rtcor <- c(rtcor, lcsets[[i]]@rt$corrected)

        mynsamp <- nrow(pdatalist[[i]])
        sampidx <- seq_len(mynsamp) + nsamp
        peaklist[[i]][,"sample"] <- sampidx[peaklist[[i]][,"sample"]]
        nsamp <- nsamp + mynsamp
    }

    peaks(object) <- do.call("rbind", peaklist)
    phenoData(object) <- do.call("rbind", pdatalist)
    filepaths(object) <- unlist(cdflist)
    # assumes all pipelines were the same (otherwise create a new xcmsSet)
    object@pipeline <- lcsets[[1]]@pipeline
    object@rt <- list(raw = rtraw, corrected = rtcor)

    invisible(object)
}

split.xcmsSet <- function(x, f, drop = TRUE, ...) {

    if (!is.factor(f))
        f <- factor(f)
    sampidx <- unclass(f)
    peakmat <- peaks(x)

    pdata <- phenoData(x)
    cdffiles <- filepaths(x)
    pipeline <- pipeline(x)
    rtraw <- x@rt$raw
    rtcor <- x@rt$corrected

    lcsets <- vector("list", length(levels(f)))
    names(lcsets) <- levels(f)

    for (i in unique(sampidx)) {
        lcsets[[i]] <- new("xcmsSet")

        samptrans <- numeric(length(f))
        samptrans[sampidx == i] <- rank(which(sampidx == i))
        samp <- samptrans[peakmat[,"sample"]]
        sidx <- which(samp != 0)
        cpeaks <- peakmat[sidx,]
        cpeaks[,"sample"] <- samp[sidx]
        peaks(lcsets[[i]]) <- cpeaks

        phenoData(lcsets[[i]]) <- pdata[,sampidx == i]
        filepaths(lcsets[[i]]) <- cdffiles[sampidx == i]
        lcsets[[i]]@rt$raw <- rtraw[sampidx == i]
        lcsets[[i]]@rt$corrected <- rtcor[sampidx == i]

        # Note that the results are likely different from
        # those obtained by processing ths subset of samples from the start.
        lcsets[[i]]@pipeline <- pipeline
    }

    if (drop)
        lcsets <- lcsets[seq(along = lcsets) %in% sampidx]

    lcsets
}

setGeneric("peaks", function(object) standardGeneric("peaks"))

setMethod("peaks", "xcmsSet", function(object) object@peaks)

setGeneric("peaks<-", function(object, value) standardGeneric("peaks<-"))

setReplaceMethod("peaks", "xcmsSet", function(object, value) {

    object@peaks <- value

    object
})

setGeneric("groups", function(object) standardGeneric("groups"))

setMethod("groups", "xcmsSet", function(object) object@groups)

setGeneric("groups<-", function(object, value) standardGeneric("groups<-"))

setReplaceMethod("groups", "xcmsSet", function(object, value) {

    object@groups <- value

    object
})

setGeneric("groupidx", function(object) standardGeneric("groupidx"))

setMethod("groupidx", "xcmsSet", function(object) object@groupidx)

setGeneric("groupidx<-", function(object, value) standardGeneric("groupidx<-"))

setReplaceMethod("groupidx", "xcmsSet", function(object, value) {

    object@groupidx <- value

    object
})

setGeneric("comps", function(object) standardGeneric("comps"))

setMethod("comps", "xcmsSet", function(object) object@comps)

setGeneric("comps<-", function(object, value) standardGeneric("comps<-"))

setReplaceMethod("comps", "xcmsSet", function(object, value) {

    object@comps <- value

    object
})


setMethod("sampnames", "xcmsSet", function(object) rownames(object@phenoData))

setGeneric("sampnames<-", function(object, value) standardGeneric("sampnames<-"))

setReplaceMethod("sampnames", "xcmsSet", function(object, value) {

    rownames(object@phenoData) <- value

    object
})

setGeneric("sampclass", function(object) standardGeneric("sampclass"))

setMethod("sampclass", "xcmsSet", function(object) {
    if (ncol(object@phenoData) >0) {
        interaction(object@phenoData)
    } else {
        factor()
    }
})

setGeneric("sampclass<-", function(object, value) standardGeneric("sampclass<-"))

setReplaceMethod("sampclass", "xcmsSet", function(object, value) {
  if (!is.factor(value))
    value <- factor(value, unique(value))

  object@phenoData <- data.frame(class = value)
  object
})

setGeneric("phenoData", function(object) standardGeneric("phenoData"))

setMethod("phenoData", "xcmsSet", function(object) object@phenoData)

setGeneric("phenoData<-", function(object, value) standardGeneric("phenoData<-"))

setReplaceMethod("phenoData", "xcmsSet", function(object, value) {
    if (is.matrix(value))
        value <- as.data.frame(value)
    #if (is.data.frame(value) && !("class" %in% colnames(value)))
    #    value[,"class"] <- interaction(value)
    #else
    if (!is.data.frame(value))
      value <- data.frame(class = value)
    object@phenoData <- value
    object
})

setGeneric("filepaths", function(object) standardGeneric("filepaths"))

setMethod("filepaths", "xcmsSet", function(object) object@filepaths)

setGeneric("filepaths<-", function(object, value) standardGeneric("filepaths<-"))

setReplaceMethod("filepaths", "xcmsSet", function(object, value) {

    object@filepaths <- value

    object
})

setMethod("groupnames", "xcmsSet", function(object, mzdec = 0, rtdec = 0,
                                            template = NULL) {

    if (!missing(template)) {
        tempsplit <- strsplit(template[1], "[T_]")
        tempsplit <- strsplit(unlist(tempsplit), "\\.")
        if (length(tempsplit[[1]]) > 1)
            mzdec <- nchar(tempsplit[[1]][2])
        else
            mzdec <- 0
        if (length(tempsplit[[2]]) > 1)
            rtdec <- nchar(tempsplit[[2]][2])
        else
            rtdec <- 0
    }

    mzfmt <- paste("%.", mzdec, "f", sep = "")
    rtfmt <- paste("%.", rtdec, "f", sep = "")

    gnames <- paste("M", sprintf(mzfmt, groups(object)[,"mzmed"]), "T",
                    sprintf(rtfmt, groups(object)[,"rtmed"]), sep = "")

    if (any(dup <- duplicated(gnames)))
        for (dupname in unique(gnames[dup])) {
            dupidx <- which(gnames == dupname)
            gnames[dupidx] <- paste(gnames[dupidx], seq(along = dupidx), sep = "_")
        }

    gnames
})

# Stage that performs a pipeline on each raw sample and aggregates the peaks

setStage("processRaws", "Load each sample and find its peaks")

.processRaws.generic <- function(data, pipeline)
{
  files <- filepaths(data)
  snames <- sampnames(data)

  peaklist <- vector("list", length(files))
  rtlist <- list(raw = vector("list", length(snames)),
                 corrected = vector("list", length(snames)))

  for (i in seq(along = peaklist)) {
    # FIXME: at some point the call to xcmsRaw could become a protocol
    # Till then deduce xcmsRaw parameter from findPeaks protocol

    proto <- protocol(pipeline, "findPeaks", "MS1")
    if (!is.null(proto))
      includeMSn <- TRUE
    else
      includeMSn <- FALSE

    lcraw <- xcmsRaw(files[i], genprof = FALSE, includeMSn = includeMSn)
    cat(snames[i], ": ", sep = "")
    peaklist[[i]] <- perform(pipeline, lcraw)
    peaklist[[i]] <- cbind(peaklist[[i]], sample = rep.int(i, nrow(peaklist[[i]])))
    rtlist$raw[[i]] <- lcraw@scantime
    rtlist$corrected[[i]] <- lcraw@scantime
    rm(lcraw)
    gc()
    if (nrow(peaklist[[i]]) == 0)
      warning(paste("No peaks found in sample", snames[i]))
    else if (nrow(peaklist[[i]]) == 1)
      warning(paste("Only 1 peak found in sample", snames[i]))
    else if (nrow(peaklist[[i]]) < 10)
      warning(paste("Only", nrow(peaklist[[i]]), "peaks found in sample",
                    snames[i]))
  }

  peaks(data) <- do.call("rbind", peaklist)
  data@rt <- rtlist

  data
}

setProtocol("generic", "Generic", representation(pipeline = "xcmsPipeline"),
            .processRaws.generic, "processRaws")

# Feature-level protocols

setStage("findComps", "Find Components")
setStage("retcor", "Correct Retention Time")
setStage("fillPeaks", "Impute Missing Peaks")
setStage("summarize", "Summarize Quantities")
setStage("normalize", "Normalize Quantities")
setStage("identifyCompounds", "Identify Compounds")

### Grouping

setStage("group", "Group Components")

.group.density <- function(object, bw = 30, minfrac = 0.5, minsamp = 1,
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
            flush.console()
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
            if (! any(gcount >= classnum*minfrac & gcount >= minsamp))
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

    groupmat <- groupmat[seq(length = num),]
    groupindex <- groupindex[seq(length = num)]

    # Remove groups that overlap with more "well-behaved" groups
    numsamp <- rowSums(groupmat[,(match("npeaks", colnames(groupmat))+1):ncol(groupmat),drop=FALSE])
    uorder <- order(-numsamp, groupmat[,"npeaks"])
    uindex <- rectUnique(groupmat[,c("mzmin","mzmax","rtmin","rtmax"),drop=FALSE],
                         uorder)

    groups(object) <- groupmat[uindex,]
    groupidx(object) <- groupindex[uindex]

    object
}

setProtocol("density", "Density Estimation",
            representation(bw = "numeric", minfrac = "numeric",
                           minsamp = "numeric", mzwid = "numeric",
                           max = "numeric"),
            .group.density, "group")

.group.mzClust <- function(object,
                           mzppm = 20,
                           mzabs = 0,
                           minsamp = 1,
                           minfrac=0.5)
{
    samples <- sampnames(object)
    classlabel <- sampclass(object)
    peaks <- peaks(object)
    groups <- xcms:::mzClustGeneric(peaks[,c(1,10)],
                                    sampclass=classlabel,
                                    mzppm=mzppm,mzabs=mzabs,
                                    minsamp=minsamp,
                                    minfrac=minfrac)

    if(is.null(nrow(groups$mat))) {
        matColNames <- names(groups$mat)
        groups$mat <- matrix(groups$mat,
                             ncol=length(groups$mat),byrow=F);
        colnames(groups$mat) <- matColNames
    }

    rt <- c(rep(-1,nrow(groups$mat)))

    groups(object) <- cbind(groups$mat[,c(1:3)],rt,rt,rt,groups$mat[,4:ncol(groups$mat)])
    colnames(groups(object)) <- c(colnames(groups$mat[,1:3]), "rtmed", "rtmin", "rtmax", colnames(groups$mat[,4:ncol(groups$mat)]))
    groupidx(object) <- groups$idx

    object
}

setProtocol("mzClust", "Spectrum Alignment",
            representation(mzppm = "numeric", mzabs = "numeric",
                           minsamp = "numeric",
                           minfrac = "numeric"),
            .group.mzClust, "group")

setGeneric("groupval", function(object, ...) standardGeneric("groupval"))

setMethod("groupval", "xcmsSet", function(object,
                                          method = c("medret", "maxint"),
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



.retcor.smooth <- function(object, missing = 1, extra = 1,
                                        model = c("loess", "linear"), span = .2,
                                        family = c("gaussian", "symmetric"),
                                        plottype = c("none", "deviation", "mdevden"),
                                        col = NULL, ty = NULL) {

    peakmat <- peaks(object)
    groupmat <- groups(object)
    if (length(groupmat) == 0)
        stop("No group information found")
    samples <- sampnames(object)
    classlabel <- as.vector(unclass(sampclass(object)))
    n <- length(samples)
    corpeaks <- peakmat
    method <- match.arg(model)
    plottype <- match.arg(plottype)
    family <- match.arg(family)
    if (length(object@rt) == 2)
        rtcor <- object@rt$corrected
    else {
        fnames <- filepaths(object)
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
    if (length(idx) == 0)
        stop("No peak groups found for retention time correction")
    idx <- idx[order(groupmat[idx,"rtmed"])]

    rt <- groupval(object, "maxint", "rt")[idx,]
    cat("Retention Time Correction Groups:", nrow(rt), "\n")
    rtdev <- rt - apply(rt, 1, median, na.rm = TRUE)

    if (method == "loess") {
        mingroups <- min(colSums(!is.na(rt)))
        if (mingroups < 4) {
            method <- "linear"
            warning("Too few peak groups, reverting to linear method")
        } else if (mingroups*span < 4) {
            span <- 4/mingroups
            warning("Span too small, resetting to ", round(span, 2))
        }
    }

    rtdevsmo <- vector("list", n)

    # Code for checking to see if retention time correction is overcorrecting
    rtdevrange <- range(rtdev, na.rm = TRUE)
    warn.overcorrect <- FALSE

    for (i in 1:n) {

        pts <- na.omit(data.frame(rt = rt[,i], rtdev = rtdev[,i]))

        if (method == "loess") {
            lo <- suppressWarnings(loess(rtdev ~ rt, pts, span = span, degree = 1, family = family))

            rtdevsmo[[i]] <- na.flatfill(predict(lo, data.frame(rt = rtcor[[i]])))
            ### Remove singularities from the loess function
            rtdevsmo[[i]][abs(rtdevsmo[[i]]) > quantile(abs(rtdevsmo[[i]]), 0.9)*2] <- NA

            if (length(naidx <- which(is.na(rtdevsmo[[i]]))))
                rtdevsmo[[i]][naidx] <- suppressWarnings(approx(na.omit(data.frame(rtcor[[i]], rtdevsmo[[i]])),
                                                                xout = rtcor[[i]][naidx], rule = 2)$y)
            while (length(decidx <- which(diff(rtcor[[i]] - rtdevsmo[[i]]) < 0))) {
                d <- diff(rtcor[[i]] - rtdevsmo[[i]])[tail(decidx, 1)]
                rtdevsmo[[i]][tail(decidx, 1)] <- rtdevsmo[[i]][tail(decidx, 1)] - d
            }

            rtdevsmorange <- range(rtdevsmo[[i]])
            if (any(rtdevsmorange/rtdevrange > 2)) warn.overcorrect <- TRUE
        } else {
            fit <- lsfit(pts$rt, pts$rtdev)
            rtdevsmo[[i]] <- rtcor[[i]] * fit$coef[2] + fit$coef[1]
            ptsrange <- range(pts$rt)
            minidx <- rtcor[[i]] < ptsrange[1]
            maxidx <- rtcor[[i]] > ptsrange[2]
            rtdevsmo[[i]][minidx] <- rtdevsmo[[i]][head(which(!minidx), n = 1)]
            rtdevsmo[[i]][maxidx] <- rtdevsmo[[i]][tail(which(!maxidx), n = 1)]
        }
    }

    if (warn.overcorrect) {
        warning(paste("Fitted retention time deviation curves exceed points by more than 2x.",
                      "This is dangerous and the algorithm is probably overcorrecting your data.",
                      "Consider increasing the span parameter or switching to the linear method.",
                      sep = "\n"))
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
}

setProtocol("smooth", "Smooth",
  representation(missing = "numeric", extra = "numeric",
    model = "character", span = "numeric", family = "character"),
  .retcor.smooth, "retcor")

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

.fillPeaks.extract <- function(object) {

    peakmat <- peaks(object)
    groupmat <- groups(object)
    if (length(groupmat) == 0)
        stop("No group information found")
    files <- filepaths(object)
    samp <- sampnames(object)
    classlabel <- as.vector(unclass(sampclass(object)))
    rawpipeline <- pipeline(pipeline(processRawsProto(object)),
                            outtype = "xcmsRaw")
    rtcor <- object@rt$corrected

    # Remove groups that overlap with more "well-behaved" groups
    numsamp <- rowSums(groupmat[,(match("npeaks", colnames(groupmat))+1):ncol(groupmat),drop=FALSE])
    uorder <- order(-numsamp, groupmat[,"npeaks"])
    uindex <- rectUnique(groupmat[,c("mzmin","mzmax","rtmin","rtmax"),drop=FALSE],
                         uorder)
    groupmat <- groupmat[uindex,]
    groupindex <- groupidx(object)[uindex]
    gvals <- groupval(object)[uindex,]

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

    cnames <- colnames(object@peaks)

    for (i in seq(along = files)) {

        cat(samp[i], "")
        flush.console()
        naidx <- which(is.na(gvals[,i]))
        if (length(naidx)) {
            lcraw <- xcmsRaw(files[i], pipeline = rawpipeline)
#            lcraw <- xcmsRaw(files[i])
            if (length(rtcor) == length(files))
                lcraw@scantime <- rtcor[[i]]
            newpeaks <- getPeaks(lcraw, peakrange[naidx,,drop=FALSE])
            rm(lcraw)
            gc()
            newpeaks <- cbind(newpeaks, sample = rep(i, length(naidx)))
            newcols <- colnames(newpeaks)[colnames(newpeaks) %in% cnames]
            peakmat[lastpeak+seq(along = naidx),newcols] <- newpeaks[,newcols]
            for (i in seq(along = naidx))
                groupindex[[naidx[i]]] <- c(groupindex[[naidx[i]]], lastpeak+i)
            lastpeak <- lastpeak + length(naidx)
        }
    }
    cat("\n")

    peaks(object) <- peakmat
    groups(object) <- groupmat
    groupidx(object) <- groupindex

    invisible(object)
}

setProtocol("extract", "Extract from raw data",
            fun = .fillPeaks.extract, parent = "fillPeaks")



.fillPeaks.mean <- function(object, value = "into") {
    peakmat <- peaks(object)
    groupmat <- groups(object)
    if (length(groupmat) == 0)
        stop("No group information found")

    classlabel <- as.vector(unclass(sampclass(object)))

    groupindex <- groupidx(object)
    gidxs <- groupval(object)
    gvals <- groupval(object, value = value)

    cnames <- colnames(object@peaks)

    classWise <- function(x, classes, fun) {
        unlist(tapply(x, classes, fun))
    }

    imputeMean <- function(x)  {
      isna <- is.na(x)
      x[which(isna)] <- mean(x, na.rm=TRUE)
      x[which(!isna)] <- NA
      x
    }

    gvals.new <- apply(gvals, 1, function(x) {classWise(x, classlabel, imputeMean)})

    samples.new  <- matrix(rep(1:nrow(gvals.new), times=ncol(gvals.new)),
                          ncol=ncol(gvals.new), nrow=nrow(gvals.new))

    newidx <- which(!is.na(gvals.new))

    newpeaks <- matrix(nrow = length(newidx), ncol = ncol(peakmat))
    newpeaks[, 1:6] <- groupmat[floor((newidx-1)/ncol(gvals))+1, 1:6]
    newpeaks[, which(cnames %in% c(value, "sample"))] <- cbind(gvals.new[newidx],
                                                               samples.new[newidx])
    colnames(peakmat) <- cnames
    peakmat.new <- rbind(peakmat, newpeaks)

    groupidx.new <- matrix(rep(1:ncol(gvals.new), each=nrow(gvals.new)),
                           ncol=ncol(gvals.new), nrow=nrow(gvals.new))
    groupindex.new <- tapply(seq(nrow(peakmat)+1, nrow(peakmat.new)),
                             groupidx.new[newidx], function(x){x})

    groupindex.both <- lapply(1:length(groupindex), function(x){
         w <- which(as.numeric(unlist(dimnames(groupindex.new))) %in% x)
         unlist(c(groupindex[[x]], na.omit(ifelse(length(w)==0, integer(0), groupindex.new[w]))))
     })

    peaks(object) <- peakmat.new
    groups(object) <- groupmat
    groupidx(object) <- groupindex.both

    invisible(object)
}

setProtocol("mean", "Insert mean value of non-NA peaks",
            fun = .fillPeaks.mean, parent = "fillPeaks")



setMethod("getEIC", "xcmsSet", function(object, mzrange, rtrange = 200,
                                        groupidx, sampleidx = sampnames(object),
                                        rt = c("corrected", "raw")) {

    files <- filepaths(object)
    grp <- groups(object)
    samp <- sampnames(object)
    rawpipeline <- pipeline(pipeline(processRawsProto(object)),
                            outtype = "xcmsRaw")

    rt <- match.arg(rt)

    if (is.numeric(sampleidx))
        sampleidx <- sampnames(object)[sampleidx]
    sampidx <- match(sampleidx, sampnames(object))

    if (!missing(groupidx)) {
        if (is.numeric(groupidx))
            groupidx <- groupnames(object)[unique(as.integer(groupidx))]
        grpidx <- match(groupidx, groupnames(object, template = groupidx))
    }

    if (missing(mzrange)) {
        if (missing(groupidx))
            stop("No m/z range or groups specified")
        # if (any(is.na(groupval(object, value = "mz")))) stop('Please use fillPeaks() to fill up NA values !')
        mzmin <- -rowMax(-groupval(object, value = "mzmin"))
        mzmax <- rowMax(groupval(object, value = "mzmax"))
        mzrange <- matrix(c(mzmin[grpidx], mzmax[grpidx]), ncol = 2)
    } else if (all(c("mzmin","mzmax") %in% colnames(mzrange)))
        mzrange <- mzrange[,c("mzmin", "mzmax"),drop=FALSE]
    else if (is.null(dim(mzrange)))
        stop("mzrange must be a matrix")
    colnames(mzrange) <- c("mzmin", "mzmax")

    if (length(rtrange) == 1) {
        if (missing(groupidx))
            rtrange <- matrix(rep(range(object@rt[[rt]][sampidx]), nrow(mzrange)),
                              ncol = 2, byrow = TRUE)
        else {
            rtrange <- retexp(grp[grpidx,c("rtmin","rtmax"),drop=FALSE], rtrange)
        }
    } else if (is.null(dim(rtrange)))
        stop("mzrange must be a matrix or single number")
    colnames(rtrange) <- c("rtmin", "rtmax")

    if (missing(groupidx))
        gnames <- character(0)
    else
        gnames <- groupidx

    eic <- vector("list", length(sampleidx))
    names(eic) <- sampleidx

    for (i in seq(along = sampidx)) {

        cat(sampleidx[i], "")
        flush.console()
        lcraw <- xcmsRaw(files[sampidx[i]], pipeline = rawpipeline)
        if (rt == "corrected")
            lcraw@scantime <- object@rt$corrected[[sampidx[i]]]
        eic[[i]] <- getEIC(lcraw, mzrange, rtrange)
        rm(lcraw)
        gc()
    }
    cat("\n")

    invisible(new("xcmsEIC", eic = eic, mzrange = mzrange, rtrange = rtrange,
                  rt = rt, groupnames = gnames))
})

setGeneric("diffreport", function(object, ...) standardGeneric("diffreport"))

setMethod("diffreport", "xcmsSet", function(object, class1 = levels(sampclass(object))[1],
                                            class2 = levels(sampclass(object))[2],
                                            filebase = character(), eicmax = 0, eicwidth = 200,
                                            sortpval = TRUE, classeic = c(class1,class2),
                                            metlin = FALSE) {

    require(multtest) || stop("Couldn't load multtest")

    groupmat <- groups(object)
    if (length(groupmat) == 0)
        stop("No group information found")
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
	if (metlin) {
	    neutralmass <- groupmat[,"mzmed"] + ifelse(metlin < 0, 1, -1)
	    metlin <- abs(metlin)
	    digits <- ceiling(-log10(metlin))+1
	    metlinurl <- paste("http://metlin.scripps.edu/metabo_list.php?mass_min=",
	                       round(neutralmass - metlin, digits), "&mass_max=",
	                       round(neutralmass + metlin, digits), sep="")
	    values <- cbind(metlin = metlinurl, values)
	}
	twosamp <- cbind(name = groupnames(object), stat, groupmat, values)
	if (sortpval) {
	   tsidx <- order(twosamp[,"pvalue"])
	   twosamp <- twosamp[tsidx,]
	   rownames(twosamp) <- 1:nrow(twosamp)
	}

    if (length(filebase))
        write.table(twosamp, paste(filebase, ".tsv", sep = ""), quote = FALSE, sep = "\t", col.names = NA)

    if (eicmax > 0) {
        eicmax <- min(eicmax, length(tsidx))
        eics <- getEIC(object, rtrange = eicwidth*1.1, sampleidx = ceic,
                       groupidx = tsidx[seq(length = eicmax)])
        if (length(filebase)) {
            eicdir <- paste(filebase, "_eic", sep="")
            dir.create(eicdir)
            if (capabilities("png"))
                png(file.path(eicdir, "%03d.png"), width = 640, height = 480)
            else
                pdf(file.path(eicdir, "%03d.pdf"), width = 640/72,
                    height = 480/72, onefile = FALSE)
        }
        plot(eics, object, rtrange = eicwidth)
        if (length(filebase))
            dev.off()
    }

    invisible(twosamp)
})

retexp <- function(peakrange, width = 200) {

    retmean <- rowMeans(peakrange[,c("rtmin", "rtmax"),drop=FALSE])
    peakrange[,"rtmin"] <- retmean-width/2
    peakrange[,"rtmax"] <- retmean+width/2

    peakrange
}

setMethod("explore", c("xcmsSet", "NULL"), function(object, protocol, ...)
{
  # something that allows the user to view the raw data in each sample
})

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

# derive experimental design from set of file paths
phenoDataFromPaths <- function(paths) {
  ## create factors from filesystem hierarchy
  sclass <- gsub("^\\.$", "sample", dirname(paths))
  lev <- strsplit(sclass, "/")
  levlen <- sapply(lev, length)
  if(length(lev) > 1 && !all(levlen[1] == levlen))
    stop("Directory tree must be level")
  pdata <- as.data.frame(matrix(unlist(lev), nrow=length(lev), byrow=TRUE))
  redundant <- apply(pdata, 2, function(col) length(unique(col)) == 1)
  if (!any(!redundant)) {
    redundant[length(redundant)] <- FALSE
  }
  pdata <- pdata[,!redundant,drop=FALSE]
  if (ncol(pdata) == 1) { ## if not multiple factors, behave as before
    ## Make the default group names less redundant
    scomp <- strsplit(substr(sclass, 1, min(nchar(sclass))), "")
    scomp <- matrix(c(scomp, recursive = TRUE), ncol = length(scomp))
    i <- 1
    while(all(scomp[i,1] == scomp[i,-1]) && i < nrow(scomp))
      i <- i + 1
    i <- min(i, tail(c(0, which(scomp[1:i,1] == .Platform$file.sep)), n = 1) + 1)
    if (i > 1 && i <= nrow(scomp))
      sclass <- substr(sclass, i, max(nchar(sclass)))
    pdata <- data.frame(class = sclass)
  }
  rownames(pdata) <- gsub("\\.[^.]*$", "", basename(paths))
  pdata
}
