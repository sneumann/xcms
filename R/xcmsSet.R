setClass("xcmsSet", representation(peaks = "matrix", groups = "matrix",
                                   groupidx = "list",
                                   filled="numeric",
                                   phenoData = "data.frame",
                                   rt = "list",
                                   filepaths = "character", profinfo = "list",
                                   polarity = "character"),
         prototype(peaks = matrix(nrow = 0, ncol = 0),
                   groups = matrix(nrow = 0, ncol = 0),
                   groupidx = list(),
                   filled = integer(0),
                   phenoData = data.frame(), rt = list(),
                   rt = list(),
                   filepaths = character(0), profinfo = vector("list"),
                   polarity = character(0)))

xcmsSet <- function(files = NULL, snames = NULL, sclass = NULL, phenoData = NULL,
                    profmethod = "bin", profparam = list(),
                    polarity = NULL, nSlaves=0, ...) {

    object <- new("xcmsSet")

    filepattern <- c("[Cc][Dd][Ff]", "[Nn][Cc]", "([Mm][Zz])?[Xx][Mm][Ll]",
                     "[Mm][Zz][Dd][Aa][Tt][Aa]", "[Mm][Zz][Mm][Ll]")
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

    if (length(files) == 0)
      stop("No NetCDF/mzXML/mzData/mzML files were found.\n")

    # determine experimental design
    if (is.null(snames)) {
      fromPaths <- phenoDataFromPaths(files)
      snames <- rownames(fromPaths)
    }
    pdata <- phenoData
    if (is.null(pdata)) {
      pdata <- sclass
      if (is.null(pdata))
        pdata <- fromPaths
    }
    phenoData(object) <- pdata
    if (is.null(phenoData))
      rownames(phenoData(object)) <- snames

    rtlist <- list(raw = vector("list", length(snames)),
                   corrected = vector("list", length(snames)))

    if ("step" %in% names(list(...)))
        profstep <- list(...)$step
    else
        profstep <- 0.1

    profinfo(object) <- c(list(method = profmethod, step = profstep), profparam)

    object@polarity <- as.character(polarity)
    includeMSn=FALSE
    xcmsSetArgs <- as.list(match.call())
    if (!is.null(xcmsSetArgs$method)) {
        if (xcmsSetArgs$method=="MS1") {
            includeMSn=TRUE
        }
    }

    runParallel <- 0

    if (nSlaves > 1) {
        ## If MPI is available ...
        rmpi = "Rmpi"
        if (require(rmpi,character.only=TRUE) && !is.null(nSlaves)) {
            if (is.loaded('mpi_initialize')) {

                mpi.spawn.Rslaves(nslaves=nSlaves, needlog=FALSE)

                ## If there are multiple slaves AND this process is the master,
                ## run in parallel.
                if ((mpi.comm.size() > 2)  && (mpi.comm.rank() == 0))
                    runParallel <- 1
            }
        }
    }

    if (runParallel==1) { ## ... we use MPI

        params <- list(...);
        params$profmethod <- profmethod;
        params$profparam <- profparam;
        params$includeMSn <- includeMSn;

        ft <- cbind(file=files,id=1:length(files))
        argList <- apply(ft,1,function(x) list(file=x["file"],id=as.numeric(x["id"]),params=params))

        res <- xcmsPapply(argList, findPeaksMPI)

        mpi.close.Rslaves()

        peaklist <- lapply(res, function(x) x$peaks)
        rtlist$raw <-  rtlist$corrected <-  lapply(res, function(x) x$scantime)

    } else {

      peaklist <- vector("list", length(files))

      for (i in seq(along = peaklist)) {

        lcraw <- xcmsRaw(files[i], profmethod = profmethod, profparam = profparam,
                          profstep = 0, includeMSn=includeMSn)
	## check existence of slot, absent in old xcmsSets
        if (exists("object@polarity") && length(object@polarity) >0) {
            ## Retain wanted polarity only
            lcraws <- split(lcraw, lcraw@polarity, DROP=TRUE)
            lcraw <- lcraws[[object@polarity]]
        }

        cat(snames[i], ": ", sep = "")
          peaklist[[i]] <- findPeaks(lcraw, ...)
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

    }

    peaks(object) <- do.call(rbind, peaklist)
    object@rt <- rtlist

    object
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
        cdflist[[i]] <- filepaths(lcsets[[i]])
        rtraw <- c(rtraw, lcsets[[i]]@rt$raw)
        rtcor <- c(rtcor, lcsets[[i]]@rt$corrected)

        sampidx <- seq(along = namelist[[i]]) + nsamp
        peaklist[[i]][,"sample"] <- sampidx[peaklist[[i]][,"sample"]]
        nsamp <- nsamp + length(namelist[[i]])
    }

    peaks(object) <- do.call(rbind, peaklist)
    sampnames(object) <- unlist(namelist)
    classlist <- unlist(classlist)
    sampclass(object) <- factor(classlist, unique(classlist))
    filepaths(object) <- unlist(cdflist)
    profinfo(object) <- profinfo(lcsets[[1]])
    object@rt <- list(raw = rtraw, corrected = rtcor)

    invisible(object)
}

split.xcmsSet <- function(x, f, drop = TRUE, ...) {

    if (!is.factor(f))
        f <- factor(f)
    sampidx <- unclass(f)
    peakmat <- peaks(x)
    samples <- sampnames(x)
    classlabel <- sampclass(x)
    cdffiles <- filepaths(x)
    prof <- profinfo(x)
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

        sampnames(lcsets[[i]]) <- samples[sampidx == i]
        sampclass(lcsets[[i]]) <- classlabel[sampidx == i, drop = TRUE]
        filepaths(lcsets[[i]]) <- cdffiles[sampidx == i]
        profinfo(lcsets[[i]]) <- prof
        lcsets[[i]]@rt$raw <- rtraw[sampidx == i]
        lcsets[[i]]@rt$corrected <- rtcor[sampidx == i]
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

setMethod("sampnames", "xcmsSet", function(object) rownames(object@phenoData))

setGeneric("sampnames<-", function(object, value) standardGeneric("sampnames<-"))

setReplaceMethod("sampnames", "xcmsSet", function(object, value) {

    if (length(object@phenoData)==0) {
        object@phenoData <- data.frame(class=rep("dummy", length(value)))
    }

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

  object@phenoData$class <- value
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

setGeneric("profinfo", function(object) standardGeneric("profinfo"))

setMethod("profinfo", "xcmsSet", function(object) object@profinfo)

setGeneric("profinfo<-", function(object, value) standardGeneric("profinfo<-"))

setReplaceMethod("profinfo", "xcmsSet", function(object, value) {

    object@profinfo <- value

    object
})

setGeneric("calibrate", function(object, ...) standardGeneric("calibrate"))
setMethod("calibrate", "xcmsSet", function(object,wishlist,method="linear",
                                           mzabs=0.0001, mzppm=5,
                                           neighbours=3, plotres=FALSE) {

    nsamp = length(unique(object@peaks[,"sample"]))
    if (!sum(method == c("shift","linear","edgeshift")))
        stop("unknown calibration method!")

    if (is.list(wishlist))
        if (length(wishlist) != nsamp)
            stop("Error: Number of masslists differs with number of samples")

    for (s in 1:nsamp)
    {
        peaklist = object@peaks[which(object@peaks[,"sample"]==s),]
        if (is.list(wishlist)) {
            masslist <- wishlist[s]
        }else{
            masslist <- wishlist
        }

        masses <- matchpeaks(peaklist,masslist,mzabs,mzppm,neighbours)
        if (length(masses)==0) stop("No masses close enough!")

        if (nrow(masses)==1 & method!="shift") {
            cat("Warning: only one peak found, fallback to shift.")
            method="shift"
        }

        params <-  estimate (masses, method)
        mzu <- peaklist[,"mz"]
        mposs <- masses[,"pos"]
        mdiffs <- masses[,"dif"]
        a <- params[1]
        b <- params[2]

        ## cat("a=",a,"b=",b)

        if (method != "edgeshift"){
            mzu <- mzu - (a * mzu + b)
        } else {
            mzu[c(1:(min(mposs)-1))] <- mzu[c(1:(min(mposs)-1))] - (a * mzu[min(mposs)] + b)
            mzu[c((min(mposs)):(max(mposs)))] <-
                mzu[c((min(mposs)):(max(mposs)))] - (a * mzu[c((min(mposs)):(max(mposs)))] + b)
            mzu[c((max(mposs)+1):length(mzu))] <-
                mzu[c((max(mposs)+1):length(mzu))] - (a * mzu[max(mposs)] + b)
        }

        peaklist[,"mz"] <- mzu
        object@peaks[which(object@peaks[,"sample"]==s),] <- peaklist
    }

    if (plotres) {
        plot(mzu[mposs],mdiffs, xlim=c(min(mzu),max(mzu)))
        if (method!="edgeshift") {abline(b,a)}else{
            lines(c(min(mzu),mzu[min(mposs)]),c(a * mzu[min(mposs)] + b,a * mzu[min(mposs)] + b))
            lines(c(mzu[min(mposs)],mzu[max(mposs)]),c(a * mzu[min(mposs)] + b,a * mzu[max(mposs)] + b))
            lines(c(mzu[max(mposs)],max(mzu)),c(a * mzu[max(mposs)] + b,a * mzu[max(mposs)] + b))
        }
    }

    invisible(object)
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

setGeneric("group.density", function(object, ...) standardGeneric("group.density"))

setMethod("group.density", "xcmsSet", function(object, bw = 30, minfrac = 0.5, minsamp = 1,
                                       mzwid = 0.25, max = 50, sleep = 0) {

    samples <- sampnames(object)
    classlabel <- sampclass(object)
    classnames <- as.character(unique(sampclass(object)))
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
})

setGeneric("group.mzClust", function(object, ...) standardGeneric("group.mzClust"))

setMethod("group.mzClust", "xcmsSet", function(object,
                           mzppm = 20,
                           mzabs = 0,
                           minsamp = 1,
                           minfrac=0.5)
{
    samples <- sampnames(object)
    classlabel <- sampclass(object)
    peaks <- peaks(object)
    groups <- xcms:::mzClustGeneric(peaks[,c("mz","sample")],
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
})

setGeneric("group.nearest", function(object, ...) standardGeneric("group.nearest"))

setMethod("group.nearest", "xcmsSet", function(object, mzVsRTbalance=10,
                                               mzCheck=0.2, rtCheck=15, kNN=10) {

    ## If ANN is available ...
    RANN = "RANN"
    if (!require(RANN)) {
        stop("RANN is not installed")
    }

    samples <- sampnames(object)
    peakmat <- peaks(object)
    plength <- list()
    parameters <- list(mzVsRTBalance=mzVsRTbalance, mzcheck=mzCheck, rtcheck=rtCheck,knn=kNN)

    for(i in 1:length(samples)){
        plength <- table(peaks(object)[,"sample"])
    }

    mplenv <- new.env(parent = .GlobalEnv)
    peakmat1 <- which(peakmat[,"sample"]==1)
    mplenv$mplist <- matrix(0,length(peakmat1),length(samples))
    mplenv$mplist[,1] <- which(peakmat[,"sample"]==1)
    mplenv$mplistmean = data.frame(peakmat[which(peakmat[,"sample"]==1),c("mz","rt")])
    mplenv$peakmat <- peakmat
    assign("peakmat",peakmat,env=mplenv)

    cat("sample:",basename(samples[1])," ")

    for(sample in 2:length(samples)) {
        for(mml in seq(mplenv$mplist[,1])){
            mplenv$mplistmean[mml,"mz"] <- mean(mplenv$peakmat[mplenv$mplist[mml,],"mz"])
            mplenv$mplistmean[mml,"rt"] <- mean(mplenv$peakmat[mplenv$mplist[mml,],"rt"])
        }

        cat("sample:",basename(samples[sample])," ")
        mplenv$peakIdxList <- data.frame(peakidx=which(mplenv$peakmat[,"sample"]==sample),
                                         isJoinedPeak=FALSE)
        if(length(mplenv$peakIdxList$peakidx)==0){
            cat("Warning: No peaks in sample",s,"\n")
        }
        scoreList <- data.frame(score=numeric(0),peak=integer(0),
                                mpListRow=integer(0),
                                isJoinedPeak=logical(0), isJoinedRow=logical(0))

        for(currPeak in mplenv$peakIdxList$peakidx){
            pvrScore <- patternVsRowScore(currPeak,parameters,mplenv)
            scoreList <- rbind(scoreList,pvrScore)
        }


        ## Browse scores in order of descending goodness-of-fit
        scoreListcurr <- scoreList[order(scoreList$score),]
        for (scoreIter in seq(scoreListcurr$score)) {

            iterPeak <-scoreListcurr$peak[scoreIter]
            iterRow <- scoreListcurr$mpListRow[scoreIter]

            ## Check if master list row is already assigned with peak
            if (scoreListcurr$isJoinedRow[scoreIter]==TRUE) {
                next
            }

            ## Check if peak is already assigned to some master list row
            if (scoreListcurr$isJoinedPeak[scoreIter]==TRUE) { next }

            ##  Check if score good enough
            ## Assign peak to master peak list row
            mplenv$mplist[iterRow,sample] <- iterPeak

            ## Mark peak as joined
            setTrue <- which(scoreListcurr$mpListRow==iterRow)
            scoreListcurr[setTrue,]$isJoinedRow <- TRUE
            setTrue <- which(scoreListcurr$peak==iterPeak)
            scoreListcurr[setTrue,]$isJoinedPeak <- TRUE
            mplenv$peakIdxList[which(mplenv$peakIdxList$peakidx==iterPeak),]$isJoinedPeak <- TRUE
        }

        notJoinedPeaks <- mplenv$peakIdxList[which(mplenv$peakIdxList$isJoinedPeak==FALSE),]$peakidx
        for(notJoinedPeak in notJoinedPeaks) {
            mplenv$mplist <- rbind(mplenv$mplist,matrix(0,1,dim(mplenv$mplist)[2]))
            mplenv$mplist[length(mplenv$mplist[,1]),sample] <- notJoinedPeak
        }

        ## Clear "Joined" information from all master peaklist rows
        rm(peakIdxList,envir=mplenv)
    }
    gc()

    groupmat <- matrix(0,nrow(mplenv$mplist), 7+length(levels(sampclass(object))))
    colnames(groupmat) <- c("mzmed", "mzmin", "mzmax", "rtmed", "rtmin", "rtmax",
                            "npeaks", levels(sampclass(object)))
    groupindex <- vector("list", nrow(mplenv$mplist))

    for (i in 1:nrow(mplenv$mplist)) {
        groupmat[i, "mzmed"] <- median(peakmat[mplenv$mplist[i,],"mz"])
        groupmat[i, c("mzmin", "mzmax")] <- range(peakmat[mplenv$mplist[i,],"mz"])
        groupmat[i, "rtmed"] <- median(peakmat[mplenv$mplist[i,],"rt"])
        groupmat[i, c("rtmin", "rtmax")] <- range(peakmat[mplenv$mplist[i,],"rt"])

        groupmat[i, "npeaks"] <- length(which(peakmat[mplenv$mplist[i,]]>0))

        groupindex[[i]] <- mplenv$mplist[i, (which(mplenv$mplist[i,]>0))]
    }

      groups(object) <- groupmat
    groupidx(object) <- groupindex

    invisible(object)
})


patternVsRowScore <- function(currPeak, parameters, mplenv)
{
    mplistmeanCurr <- mplenv$mplistmean[,c("mz","rt")]
    mplistmeanCurr[,"mz"] <- mplistmeanCurr[,"mz"] * parameters$mzVsRTBalance
    peakmatCurr <- mplenv$peakmat[currPeak,c("mz","rt"),drop=FALSE]
    peakmatCurr[,"mz"] <- peakmatCurr[,"mz"] * parameters$mzVsRTBalance

    nnDist <- nn2(mplistmeanCurr,peakmatCurr[,c("mz","rt"),drop=FALSE],
                  k=min(length(mplistmeanCurr[,1]),parameters$knn))

    scoreListcurr <- data.frame(score=numeric(0),peak=integer(0), mpListRow=integer(0),
                                isJoinedPeak=logical(0), isJoinedRow=logical(0))

    for(mplRow in 1:length(nnDist$nn.idx)){
        mplistMZ <- mplenv$mplistmean[nnDist$nn.idx[mplRow],"mz"]
        mplistRT <- mplenv$mplistmean[nnDist$nn.idx[mplRow],"rt"]

        ## Calculate differences between M/Z and RT values of current peak and median of the row

        diffMZ = abs(mplistMZ-mplenv$peakmat[[currPeak,"mz"]])
        diffRT = abs(mplistRT-mplenv$peakmat[[currPeak,"rt"]])
        ## What type of RT tolerance is used?

        rtTolerance=0
        rtTolerance = parameters$rtcheck


        ## Calculate if differences within tolerancdiffRT < rtTolerance)es
        if ( (diffMZ < parameters$mzcheck)& (diffRT < rtTolerance) ) {
            scoreListcurr <- rbind(scoreListcurr,
                                   data.frame(score=nnDist$nn.dists[mplRow],
                                              peak=currPeak, mpListRow=nnDist$nn.idx[mplRow],
                                              isJoinedPeak=FALSE, isJoinedRow=FALSE))
            ## goodEnough = true
            return(scoreListcurr)
        }
    }
}

setGeneric("group", function(object, ...) standardGeneric("group"))

setMethod("group", "xcmsSet", function(object, method=getOption("BioC")$xcms$group.method,
                                       ...) {

    method <- match.arg(method, getOption("BioC")$xcms$group.methods)
    if (is.na(method))
        stop("unknown method : ", method)
    method <- paste("group", method, sep=".")
    invisible(do.call(method, alist(object, ...)))
})

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

setGeneric("retcor", function(object, ...) standardGeneric("retcor"))

setMethod("retcor", "xcmsSet", function(object, method=getOption("BioC")$xcms$retcor.method,
                                        ...) {

    ## Backward compatibility for old "methods"
    if (method == "linear" || method == "loess") {
        return(invisible(do.call(retcor.peakgroups, alist(object, smooth=method, ...))))
    }

    method <- match.arg(method, getOption("BioC")$xcms$retcor.methods)
    if (is.na(method))
        stop("unknown method : ", method)
    method <- paste("retcor", method, sep=".")

    invisible(do.call(method, alist(object, ...)))
})

setGeneric("retcor.peakgroups", function(object, ...) standardGeneric("retcor.peakgroups"))
setMethod("retcor.peakgroups", "xcmsSet", function(object, missing = 1, extra = 1,
                                        smooth = c("loess", "linear"), span = .2,
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
    smooth <- match.arg(smooth)
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

    if (smooth == "loess") {
        mingroups <- min(colSums(!is.na(rt)))
        if (mingroups < 4) {
            smooth <- "linear"
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

        if (smooth == "loess") {
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
                      "Consider increasing the span parameter or switching to the linear smoothing method.",
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

        rtrange <- range(do.call(c, rtcor))
        devrange <- range(do.call(c, rtdevsmo))

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

setGeneric("retcor.obiwarp", function(object, ...) standardGeneric("retcor.obiwarp"))
setMethod("retcor.obiwarp", "xcmsSet", function(object, plottype = c("none", "deviation", "mdevden"),
                                                col = NULL, ty = NULL, profStep=1, r=NULL, g=NULL, cor = NULL, l=NULL, i_=0) {


    peakmat <- peaks(object)
    samples <- sampnames(object)
    classlabel <- as.vector(unclass(sampclass(object)))
    n <- length(samples)
    corpeaks <- peakmat
    plottype <- match.arg(plottype)

    if (length(object@rt) == 2) {
        rtcor <- object@rt$corrected
    } else {
        fnames <- filepaths(object)
        rtcor <- vector("list", length(fnames))
        for (i in seq(along = fnames)) {
            cdf <- netCDFOpen(fnames[i])
            rtcor[[i]] <- netCDFVarDouble(cdf, "scan_acquisition_time")
            netCDFClose(cdf)
        }
        object@rt <- list(raw = rtcor, corrected = rtcor)
    }
    rtimecor <- vector("list",n)
    rtdevsmo <- vector("list", n)

    if(!is.null(r)) {
        r <- c("-r",r)
    }
    if(!is.null(cor)) {
        cor <- c("-s",cor)
    }
    if(!is.null(l)){
        l <- "-l"
        i_ <- c("-i",i_)
    }
    if(!is.null(g)) {
        g <- c("-g",g)
    }

    plength <- list()
    for(i in 1:length(samples)){
        plength[i] <-length(which(peakmat[,"sample"]==i))
    }
    maxsample <- which.max(plength)
    idx <- which(seq(plength)!=maxsample)
    cat("center sample: ",samples[maxsample], "\n")

    for (s in idx) {
        cat("processing sample: ",samples[s], "\n")

        obj1 <- xcmsRaw(object@filepaths[maxsample], profmethod="bin", profstep=profStep)
        obj2 <- xcmsRaw(object@filepaths[s], profmethod="bin", profstep=profStep)

        mzmin <-  min(obj1@mzrange[1], obj2@mzrange[1])
        mzmax <-  max(obj1@mzrange[2], obj2@mzrange[2])

        mz <- seq(mzmin,mzmax, by=profStep)
        mz <- as.double(mz)
        mzval <- length(mz)

        scantime1 <- obj1@scantime
        scantime2 <- obj2@scantime

        mstdiff <- median(c(diff(scantime1), diff(scantime2)))

        rtup1 <- c(1:length(scantime1))
        rtup2 <- c(1:length(scantime2))

        mst1 <- which(diff(scantime1)>5*mstdiff)[1]
        if(!is.na(mst1)) {
            rtup1 <- which(rtup1<=mst1)
            cat("Found gaps: cut scantime-vector at ", scantime1[mst1],"seconds", "\n")
        }

        mst2 <- which(diff(scantime2)>5*mstdiff)[1]
        if(!is.na(mst2)) {
            rtup2 <- which(rtup2<=mst2)
            cat("Found gaps: cut scantime-vector at ", scantime2[mst2],"seconds", "\n")
        }

        scantime1 <- scantime1[rtup1]
        scantime2 <- scantime2[rtup2]

        rtmaxdiff <- abs(diff(c(scantime1[length(scantime1)],
                                scantime2[length(scantime2)])))
        if(rtmaxdiff>(5*mstdiff)){
            rtmax <- min(scantime1[length(scantime1)],
                         scantime2[length(scantime2)])
            rtup1 <- which(scantime1<=rtmax)
            rtup2 <- which(scantime2<=rtmax)
        }

        scantime1 <- scantime1[rtup1]
        scantime2 <- scantime2[rtup2]
        valscantime1 <- length(scantime1)
        valscantime2 <- length(scantime2)

        if(length(obj1@scantime)>valscantime1) {
            obj1@env$profile <- obj1@env$profile[,-c((valscantime1+1):length(obj1@scantime))]
        }
        if(length(obj2@scantime)>valscantime2) {
            obj2@env$profile <- obj2@env$profile[,-c((valscantime2+1):length(obj2@scantime))]
        }

        if(mzmin < obj1@mzrange[1]) {
            seqlen <- length(seq(mzmin, obj1@mzrange[1], profStep))-1
            x <- matrix(0, seqlen,dim(obj1@env$profile)[2])
            obj1@env$profile <- rbind(x, obj1@env$profile)
        }
        if(mzmax > obj1@mzrange[2]){
            seqlen <- length(seq(obj1@mzrange[2], mzmax, profStep))-1
            x <- matrix(0, seqlen, dim(obj1@env$profile)[2])
            obj1@env$profile <- rbind(obj1@env$profile, x)
        }
        if(mzmin < obj2@mzrange[1]){
            seqlen <- length(seq(mzmin, obj2@mzrange[1], profStep))-1
            x <- matrix(0, seqlen, dim(obj2@env$profile)[2])
            obj2@env$profile <- rbind(x, obj2@env$profile)
        }
        if(mzmax > obj2@mzrange[2]){
            seqlen <- length(seq(obj2@mzrange[2], mzmax, profStep))-1
            x <- matrix(0, seqlen, dim(obj2@env$profile)[2])
            obj2@env$profile <- rbind(obj2@env$profile, x)
        }

        intensity1 <- obj1@env$profile
        intensity2 <- obj2@env$profile

        rtimecor[[s]] <-.Call("R_set_from_xcms",
                              valscantime1,scantime1,mzval,mz,intensity1,
                              valscantime2,scantime2,mzval,mz,intensity2,
                              c(" ",r,cor,g,l,i_))

        if(length(obj2@scantime) > valscantime2) {
            object@rt$corrected[[s]] <- c(rtimecor[[s]],
                                          obj2@scantime[(max(rtup2)+1):length(obj2@scantime)])
        } else {
            object@rt$corrected[[s]] <- rtimecor[[s]]
        }

        rtdevsmo[[s]] <- round(rtcor[[s]]-object@rt$corrected[[s]],2)

        rm(obj1,obj2)
        gc()
    }

    rtdevsmo[[maxsample]] <- round(rtcor[[maxsample]] - object@rt$corrected[[maxsample]], 2)

    if (plottype == "mdevden") {
        split.screen(matrix(c(0, 1, .3, 1, 0, 1, 0, .3), ncol = 4, byrow = TRUE))
        screen(1)
        par(mar = c(0, 4.1, 4.1, 2), xaxt = "n")
    }

    if (plottype %in% c("deviation", "mdevden")) {

        ## Set up the colors and line type
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

        plot(0, 0, type="n", xlim = rtrange*1.4, ylim = devrange,
             main = "Retention Time Deviation vs. Retention Time",
             xlab = "Retention Time", ylab = "Retention Time Deviation")
        legend(rtrange[2], devrange[2],
               basename(samples), col = mypal[col],
               lty = ty, pch = ceiling(1:n/length(mypal)), xjust = -0.01)

        for (i in 1:n) {
            points(rtcor[[i]], rtdevsmo[[i]], type="l", col = mypal[col[i]], lty = ty[i])
        }
    }

    if (plottype == "mdevden") {
        screen(2)
        par(mar = c(5.1, 4.1, 0, 2), yaxt = "n")
        allden <- density(peakmat[,"rt"], bw = diff(rtrange)/200,
                          from = rtrange[1], to = rtrange[2])[c("x","y")]
        corden <- density(rt, bw = diff(rtrange)/200,
                          from = rtrange[1], to = rtrange[2], na.rm = TRUE)[c("x","y")]
        allden$y <- allden$y / sum(allden$y)
        corden$y <- corden$y / sum(corden$y)
        maxden <- max(allden$y, corden$y)
        plot(c(0,0), xlim = rtrange, ylim = c(0, maxden),
             type = "n", main = "",
             xlab = "Retention Time", ylab = "Peak Density")
        points(allden, type = "l", col = 1)
        points(corden, type = "l", col = 2)
        abline(h = 0, col = "grey")

        legend(rtrange[2], maxden, c("All", "Correction"),
               col = 1:2, lty = c(1,1), xjust = 1)
        close.screen(all = TRUE)
    }

    for (i in 1:n) {
        cfun <- stepfun(rtcor[[i]][-1] - diff(rtcor[[i]])/2, rtcor[[i]] - rtdevsmo[[i]])
        rtcor[[i]] <- rtcor[[i]] - rtdevsmo[[i]]

        sidx <- which(corpeaks[,"sample"] == i)
        corpeaks[sidx, c("rt", "rtmin", "rtmax")] <- cfun(corpeaks[sidx, c("rt", "rtmin", "rtmax")])
    }

    peaks(object) <- corpeaks
    groups(object) <- matrix(nrow = 0, ncol = 0)
    groupidx(object) <- list()
    invisible(object)

})


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

    rtrange <- range(do.call(c, rtuncor))
    devrange <- range(do.call(c, rtdevsmo))

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

setGeneric("fillPeaks.chrom", function(object, ...) standardGeneric("fillPeaks.chrom"))

setMethod("fillPeaks.chrom", "xcmsSet", function(object) {

    peakmat <- peaks(object)
    groupmat <- groups(object)
    if (length(groupmat) == 0)
        stop("No group information found")
    files <- filepaths(object)
    samp <- sampnames(object)
    classlabel <- as.vector(unclass(sampclass(object)))
    prof <- profinfo(object)
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
    peakrange[,"mzmin"] <- apply(mzmin, 1, median, na.rm = TRUE)
    mzmax <- peakmat[gvals,"mzmax"]
    dim(mzmax) <- c(nrow(gvals), ncol(gvals))
    peakrange[,"mzmax"] <- apply(mzmax, 1, median, na.rm = TRUE)
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
            lcraw <- xcmsRaw(files[i], profmethod = prof$method, profstep = 0)
	    ## check existence of slot, absent in old xcmsSets
	    if (exists("object@polarity") && length(object@polarity) >0) {
            ## Retain wanted polarity only
            lcraws <- split(lcraw, lcraw@polarity, DROP=TRUE)
            lcraw <- lcraws[[object@polarity]]
        }

            if (length(prof) > 2)
                lcraw@profparam <- prof[seq(3, length(prof))]
            if (length(rtcor) == length(files))
                lcraw@scantime <- rtcor[[i]]
            newpeaks <- getPeaks(lcraw, peakrange[naidx,,drop=FALSE], step = prof$step)
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
    object@filled <- seq((lastpeak+1),nrow(peakmat))
    groups(object) <- groupmat
    groupidx(object) <- groupindex

    invisible(object)
})

setGeneric("fillPeaks.MSW", function(object, ...) standardGeneric("fillPeaks.MSW"))
setMethod("fillPeaks.MSW", "xcmsSet", function(object, mrange=c(0,0)) {

    peakmat <- peaks(object)
    groupmat <- groups(object)
    if (length(groupmat) == 0)
        stop("No group information found")
    files <- filepaths(object)
    samp <- sampnames(object)
    classlabel <- as.vector(unclass(sampclass(object)))
    prof <- profinfo(object)
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
	newpeaks <- matrix(nrow=length(naidx), ncol=9)
	nppos<-0 ## line position in the newpeaks matrix
        if (length(naidx)) {
            lcraw <- xcmsRaw(files[i], profmethod = prof$method, profstep = 0)
	   ngs <- as.vector(naidx)
	   for (g in ngs)
		{
		nppos <- nppos+1
		mzpos <- which(abs(lcraw@env$mz - groupmat[g,"mzmed"]) == min(abs(lcraw@env$mz - groupmat[g,"mzmed"])))
		mmzpos <- mzpos[which(lcraw@env$intensity[mzpos] == max(lcraw@env$intensity[mzpos]))]
		mmzr <- seq((mmzpos-mrange[1]),(mmzpos+mrange[2]))
		maxo <- max(lcraw@env$intensity[mmzr])
		## this is the new one, summing the scale-range
		## calculating scale, adding intensitiesin this scale
		medMZmin <- median(peakmat[groupindex[[g]],"mzmin"])
		medMZmax <- median(peakmat[groupindex[[g]],"mzmax"])
		minMzpos <- min(which(abs(lcraw@env$mz - medMZmin) == min(abs(lcraw@env$mz - medMZmin))))
		maxMzpos <- max(which(abs(lcraw@env$mz - medMZmax) == min(abs(lcraw@env$mz - medMZmax))))
		into = sum(lcraw@env$intensity[minMzpos:maxMzpos])
		newpeaks[nppos,] <- c(groupmat[g,"mzmed"],medMZmin,medMZmax,-1,-1,-1,into,maxo,i)
		}
	    colnames(newpeaks) <- c("mz","mzmin","mzmax","rt","rtmin","rtmax","into","maxo","sample")
            rm(lcraw)
            gc()
            newcols <-colnames(newpeaks)[colnames(newpeaks) %in% cnames]
	    peakmat[lastpeak+seq(along = naidx),newcols] <- newpeaks[,newcols]
            for (i in seq(along = naidx))
                groupindex[[naidx[i]]] <- c(groupindex[[naidx[i]]], lastpeak+i)
            lastpeak <- lastpeak + length(naidx)
        }
    }
    cat("\n")
    peaks(object) <- peakmat
    object@filled <- seq((lastpeak+1),nrow(peakmat))
    groups(object) <- groupmat
    groupidx(object) <- groupindex
    invisible(object)
})

setGeneric("fillPeaks", function(object, ...) standardGeneric("fillPeaks"))
 setMethod("fillPeaks", "xcmsSet", function(object, method=getOption("BioC")$xcms$fillPeaks.method,...) {
 	method <- match.arg(method, getOption("BioC")$xcms$fillPeaks.methods)
     	if (is.na(method))
         	stop("unknown method : ", method)
     	method <- paste("fillPeaks", method, sep=".")
     	invisible(do.call(method, alist(object, ...)))
 	})

setMethod("getEIC", "xcmsSet", function(object, mzrange, rtrange = 200,
                                        groupidx, sampleidx = sampnames(object),
                                        rt = c("corrected", "raw")) {

    files <- filepaths(object)
    grp <- groups(object)
    samp <- sampnames(object)
    prof <- profinfo(object)

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
        if (any(is.na(groupval(object, value = "mz"))))
            stop('Please use fillPeaks() to fill up NA values !')
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
        lcraw <- xcmsRaw(files[sampidx[i]], profmethod = prof$method, profstep = 0)
        if (rt == "corrected")
            lcraw@scantime <- object@rt$corrected[[sampidx[i]]]
        if (length(prof) > 2)
            lcraw@profparam <- prof[seq(3, length(prof))]
        eic[[i]] <- getEIC(lcraw, mzrange, rtrange, step = prof$step)
        rm(lcraw)
        gc()
    }
    cat("\n")

    invisible(new("xcmsEIC", eic = eic, mzrange = mzrange, rtrange = rtrange,
                  rt = rt, groupnames = gnames))
})


getSpecWindow <- function(xs, gidxs, borderwidth=1)
{
    groupidx <- groupidx(xs)
    if (length(groupidx)==0) {
        stop("no groups found in xcmsSet object.")
    }
    minmaxs <- matrix(ncol=2, nrow=length(gidxs))
    cat("Processing data from sample: ")
    for (a in 1:length(gidxs)){
        ## 1st step: getting boundaries
        mzmin <- min(peaks(xs)[groupidx[[gidxs[a]]],"mzmin"]) ## lower bound
        mzmax <- max(peaks(xs)[groupidx[[gidxs[a]]],"mzmax"]) ## upper bound
        mzw <- mzmax-mzmin
        minmaxs[a,1] <- mzmin - borderwidth*mzw
        minmaxs[a,2] <- mzmax + borderwidth*mzw ## a peakwidth left and right
    }
    mzlistlist=list()
    for (s in 1:length(gidxs)) mzlistlist[[s]] <- list()
    mzlistlist$minmax <- minmaxs
    for (s in 1:length(sampnames(xs))){
        ## 2nd Step: getting each sample
        cat(" ",s)
        xr <- xcmsRaw(xs@filepaths[s])
        for (a in 1:length(gidxs)){
            pmin <- min(which(abs(xr@env$mz - minmaxs[a,1]) == min(abs(xr@env$mz - minmaxs[a,1]))))
            pmax <- min(which(abs(xr@env$mz - minmaxs[a,2]) == min(abs(xr@env$mz - minmaxs[a,2]))))
            mzlistlist[[a]][[s]] <- matrix(ncol=2, nrow=(pmax-pmin+1),
                                           data=c(xr@env$mz[pmin:pmax],xr@env$intensity[pmin:pmax]))
        }
    }
    cat("\n")
    invisible(mzlistlist)
}

plotSpecWindow <- function(xs, gidxs, borderwidth=1)
{
    if (length(groupidx(xs))==0) {
        stop("no groups found in xcmsSet object.")
    }

    mzll <- getSpecWindow(xs, gidxs, borderwidth)
    minmax <- mzll$minmax
    groupidx <- groupidx(xs)
    pcolors <- c("black","darkred", "green","blue")
    ecolors <- c("grey","red","lightgreen","lightblue")

    cat("\ngroup: ")
    for (a in 1:length(gidxs)){
        cat(" ",gidxs[a])
        nsa <- length(levels(sampclass(xs)))
        pcol <- pcolors[which(levels(sampclass(xs)) == sampclass(xs)[1])]
        ecol <- ecolors[which(levels(sampclass(xs)) == sampclass(xs)[1])]
        mzmin <- minmax[a,1]
        mzmax <- minmax[a,2]
        maxints <- NA

        for (s in 1:length(mzll[[a]])) {
            maxints[s] <- max(mzll[[a]][[s]][,2])
        }

        maxo <- max(maxints)
        for (n in 1:length(mzll[[a]])){
            par(lty=1)
            pcol <- pcolors[which(levels(sampclass(xs)) == sampclass(xs)[n])]
            ecol <- ecolors[which(levels(sampclass(xs)) == sampclass(xs)[n])]
            if(n==1) {
                plot(mzll[[a]][[n]][,1],
                     mzll[[a]][[n]][,2],
                     xlim=c(mzmin,mzmax), ylim=c(0,(maxo+10000)),
                     type='l', col=ecol, xlab="", ylab="") ## complete raw-range
            } else {
                lines(mzll[[a]][[n]][,1],
                      mzll[[a]][[n]][,2],
                      xlim=c(mzmin,mzmax), ylim=c(0,(maxo+10000)),
                      type='l', col=ecol, xlab="", ylab="") ## complete raw-range}
            }
            if (length(which(peaks(xs)[groupidx[[gidxs[a]]],"sample"] == n))>0){
                ## peak entry of the first sample
                apeak <- groupidx[[gidxs[a]]][min(which(peaks(xs)[groupidx[[gidxs[a]]],"sample"] == n))]
                if (apeak %in% xs@filled) {
                    par(lty=2)
                }else{
                    par(lty=1)
                }
                ppmin <- min(which(abs(mzll[[a]][[n]][,1] - peaks(xs)[apeak,"mzmin"]) ==
                                   min(      abs(mzll[[a]][[n]][,1] - peaks(xs)[apeak,"mzmin"]))))
                ppmax <- min(which(abs(mzll[[a]][[n]][,1] - peaks(xs)[apeak,"mzmax"]) ==
                                   min(      abs(mzll[[a]][[n]][,1] - peaks(xs)[apeak,"mzmax"]))))
                lines(mzll[[a]][[n]][ppmin:ppmax,1],mzll[[a]][[n]][ppmin:ppmax,2],
                      xlim=c(mzmin,mzmax), col=pcol, type='l')
            }
        }

        title(main=paste("m/z vs. intensity for group",gidxs[a]), xlab="m/z", ylab="intensity")
        legend("topright", as.vector(levels(sampclass(xs))),col=pcolors[1:nsa],lty=rep(1,nsa))
    }
    cat("\n")
}


setGeneric("diffreport", function(object, ...) standardGeneric("diffreport"))

setMethod("diffreport", "xcmsSet", function(object, class1 = levels(sampclass(object))[1],
                                            class2 = levels(sampclass(object))[2],
                                            filebase = character(), eicmax = 0, eicwidth = 200,
                                            sortpval = TRUE, classeic = c(class1,class2),
                                            value = c("into","maxo","intb"), metlin = FALSE, h=480,w=640, ...) {

    require(multtest) || stop("Couldn't load multtest")

    value <- match.arg(value)
    groupmat <- groups(object)
    if (length(groupmat) == 0)
        stop("No group information found")
    samples <- sampnames(object)
    n <- length(samples)
    classlabel <- sampclass(object)
    classlabel <- levels(classlabel)[as.vector(unclass(classlabel))]

    values <- groupval(object, "medret", value=value)
    indecies <- groupval(object, "medret", value = "index")

    if (!all(c(class1,class2) %in% classlabel))
        stop("Incorrect Class Labels")

    ## c1 and c2 are column indices of class1 and class2 resp.
    c1 <- which(classlabel %in% class1)
    c2 <- which(classlabel %in% class2)
    ceic <- which(classlabel %in% classeic)
    if (length(intersect(c1, c2)) > 0)
        stop("Intersecting Classes")

    mean1 <- rowMeans(values[,c1], na.rm = TRUE)
    mean2 <- rowMeans(values[,c2], na.rm = TRUE)

    ## Calculate fold change.
    ## For foldchange <1 set fold to 1/fold
    ## See tstat to check which was higher
    fold <- mean2 / mean1
    fold[!is.na(fold) & fold < 1] <- 1/fold[!is.na(fold) & fold < 1]

    testval <- values[,c(c1,c2)]
	testclab <- c(rep(0,length(c1)),rep(1,length(c2)))
	tstat <- mt.teststat(testval, testclab, ...)
	pvalue <- pval(testval, testclab, tstat)
	stat <- data.frame(fold = fold, tstat = tstat, pvalue = pvalue)
	if(length(levels(sampclass(object))) >2){
            pvalAnova<-c()
	    for(i in 1:nrow(values)){
		var<-as.numeric(values[i,])
		ano<-summary(aov(var ~ sampclass(object)) )
                pvalAnova<-append(pvalAnova, unlist(ano)["Pr(>F)1"])
	    }
	    stat<-cbind(stat, anova= pvalAnova)
	}
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
	   values<-values[tsidx,]
	} else
	   tsidx <- 1:nrow(values)

    if (length(filebase))
        write.table(twosamp, paste(filebase, ".tsv", sep = ""), quote = FALSE, sep = "\t", col.names = NA)

    if (eicmax > 0) {
        if (length(unique(peaks(object)[,"rt"])) > 1) {
            ## This looks like "normal" LC data

            eicmax <- min(eicmax, length(tsidx))
            eics <- getEIC(object, rtrange = eicwidth*1.1, sampleidx = ceic,
                           groupidx = tsidx[seq(length = eicmax)])

            if (length(filebase)) {
                eicdir <- paste(filebase, "_eic", sep="")
                boxdir <- paste(filebase, "_box", sep="")
                dir.create(eicdir)
                dir.create(boxdir)
                if (capabilities("png")){
                    xcmsBoxPlot(values[seq(length = eicmax),],
                                sampclass(object), dirpath=boxdir, pic="png",  width=w, height=h)
                    png(file.path(eicdir, "%003d.png"), width = w, height = h)
                } else {
                    xcmsBoxPlot(values[seq(length = eicmax),],
                                sampclass(object), dirpath=boxdir, pic="pdf", width=w, height=h)
                    pdf(file.path(eicdir, "%003d.pdf"), width = w/72,
                        height = h/72, onefile = FALSE)
                }
            }
            plot(eics, object, rtrange = eicwidth)

            if (length(filebase))
                dev.off()
        } else {
            ## This looks like a direct-infusion single spectrum
            if (length(filebase)) {
                specdir <- paste(filebase, "_spec", sep="")
                dir.create(specdir)
                if (capabilities("png")){
                    png(file.path(specdir, "%003d.png"), width = w, height = h)
                }else{
                    pdf(file.path(eicdir, "%003d.pdf"), width = w/72,
                        height = h/72, onefile = FALSE)
                }
            }

            plotSpecWindow(object, gidxs = tsidx[seq(length = eicmax)], borderwidth=1)

            if (length(filebase))
                dev.off()
        }
    }

    invisible(twosamp)
})

xcmsBoxPlot<-function(values, className, dirpath, pic, width=640, height=480)
{


    if (pic == "png"){
	png(file.path(dirpath, "%003d.png"), width, height)
    } else{
	pdf(file.path(dirpath, "%003d.pdf"), width = width/72, height = height/72, onefile = FALSE)
    }

    ind<-which(colnames(values) != "metlin")
    for (i in 1:nrow(values)){
	boxplot(as.numeric(values[i,ind]) ~ className, col="blue",
                outline=FALSE, main=paste("Feature ", row.names(values)[i] ))
    }
    if (length(values) > 0) {
        dev.off()
    }
}

retexp <- function(peakrange, width = 200) {

    retmean <- rowMeans(peakrange[,c("rtmin", "rtmax"),drop=FALSE])
    peakrange[,"rtmin"] <- retmean-width/2
    peakrange[,"rtmax"] <- retmean+width/2

    peakrange
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
