## All methods for xcmsSet should go here.
#' @include functions-xcmsSet.R

############################################################
## show
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

    if(.hasSlot(object, "mslevel")){
        MSn <- mslevel(object)
        if(is.null(MSn))
            MSn <- 1
        cat(paste0("Peak picking was performed on MS", MSn, ".\n"))
    }
    if(.hasSlot(object, "scanrange")){
        if(!is.null(scanrange(object))){
            cat("Scan range limited to ", scanrange(object)[1], "-", scanrange(object)[2], "\n")
        }
    }

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



#' @description This method updates an \emph{old} \code{xcmsSet} object to the latest
#' definition.
#' @title Update an \code{xcmsSet} object
#' @param object The \code{xcmsSet} object to update.
#' @param ... Optional additional arguments. Currently ignored.
#' @param verbose Currently ignored.
#' @return An updated \code{xcmsSet} containing all data from the input object.
#' @author Johannes Rainer
setMethod("updateObject", "xcmsSet", function(object, ..., verbose = FALSE) {
    ## Create a new empty xcmsSet and start filling it with the slot
    ## values (if present.)
    ## Define the default values.
    pks <- matrix(nrow = 0, ncol = 0)
    grps <- matrix(nrow = 0, ncol = 0)
    grpidx <- list()
    flld <- integer(0)
    phnD <- data.frame()
    theRt <- list()
    flpth <- character(0)
    prfnf <- vector("list")
    datCorr <- integer(0)
    pol <- character(0)
    prgInfo <- list()
    msl <- numeric(0)
    scnr <- numeric(0)
    prgCb <- function(progress) NULL

    ## Now replace the values with the slots... if present.
    if (.hasSlot(object, "peaks"))
        pks <- object@peaks
    if (.hasSlot(object, "groups"))
        grps <- object@groups
    if (.hasSlot(object, "groupidx"))
        grpidx <- object@groupidx
    if (.hasSlot(object, "filles"))
        flld <- object@filles
    if (.hasSlot(object, "phenoData"))
        phnD <- object@phenoData
    if (.hasSlot(object, "rt"))
        theRt <- object@rt
    if (.hasSlot(object, "filepaths"))
        flpth <- object@filepaths
    if (.hasSlot(object, "profinfo"))
        prfnf <- object@profinfo
    if (.hasSlot(object, "dataCorrection"))
        datCorr <- object@dataCorrection
    if (.hasSlot(object, "polarity"))
        pol <- object@polarity
    if (.hasSlot(object, "polarity"))
        prgInfo <- object@progressInfo
    if (.hasSlot(object, "mslevel"))
        msl <- object@mslevel
    if (.hasSlot(object, "scanrange"))
        scnr <- object@scanrange
    if (.hasSlot(object, "progressCallback"))
        prgCb <- object@progressCallback

    ## Generate the new object.
    newXs <- new("xcmsSet",
                 peaks = pks,
                 groups = grps,
                 groupidx = grpidx,
                 filled = flld,
                 phenoData = phnD,
                 rt = theRt,
                 filepaths = flpth,
                 profinfo = prfnf,
                 dataCorrection = datCorr,
                 polarity = pol,
                 progressInfo = prgInfo,
                 mslevel = msl,
                 scanrange = scnr,
                 progressCallback = prgCb
                 )
    return(newXs)
})

############################################################
## peaks
setMethod("peaks", "xcmsSet", function(object) object@peaks)
setReplaceMethod("peaks", "xcmsSet", function(object, value) {
    object@peaks <- value
    object
})

############################################################
## groups
setMethod("groups", "xcmsSet", function(object) object@groups)
setReplaceMethod("groups", "xcmsSet", function(object, value) {
    object@groups <- value
    object
})

############################################################
## groupidx
setMethod("groupidx", "xcmsSet", function(object) object@groupidx)
setReplaceMethod("groupidx", "xcmsSet", function(object, value) {
    object@groupidx <- value
    object
})

############################################################
## sampnames
setMethod("sampnames", "xcmsSet", function(object) rownames(object@phenoData))
setReplaceMethod("sampnames", "xcmsSet", function(object, value) {
    if (length(object@phenoData)==0) {
        object@phenoData <- data.frame(class=rep("dummy", length(value)))
    }
    rownames(object@phenoData) <- value
    object
})

############################################################
## sampclass
setMethod("sampclass", "xcmsSet", function(object) {
              if (ncol(object@phenoData) >0) {
                  if(any(colnames(object@phenoData)=="class")){
                      sclass <- object$class
                      ## in any rate: transform class to a character vector
                      ## and generate a new factor on that with the levels
                      ## being in the order of the first occurrence of the
                      ## elements (i.e. no alphanumeric ordering).
                      sclass <- as.character(sclass)
                      sclass <- factor(sclass, levels=unique(sclass))
                      return(sclass)
                  }
        interaction(object@phenoData, drop=TRUE)
    } else {
        factor()
    }
})
setReplaceMethod("sampclass", "xcmsSet", function(object, value) {
    ## if we're submitting a data.frame, we're using interaction on that.
    if(class(value)=="data.frame"){
        message("Setting the class labels as the interaction of the data.frame columns.")
        value <- as.character(interaction(value, drop=TRUE))
    }
    if (!is.factor(value))
        value <- factor(value, unique(value))
    object@phenoData$class <- value
    object
})

############################################################
## phenoData
setMethod("phenoData", "xcmsSet", function(object) object@phenoData)
setReplaceMethod("phenoData", "xcmsSet", function(object, value) {
    if (is.matrix(value))
        value <- as.data.frame(value)
    ## if (is.data.frame(value) && !("class" %in% colnames(value)))
    ##     value[,"class"] <- interaction(value)
    ## else
    if (!is.data.frame(value))
        value <- data.frame(class = value)
    object@phenoData <- value
    object
})

############################################################
## filepaths
setMethod("filepaths", "xcmsSet", function(object) object@filepaths)
setReplaceMethod("filepaths", "xcmsSet", function(object, value) {
    object@filepaths <- value
    object
})

############################################################
## profinfo
setMethod("profinfo", "xcmsSet", function(object) object@profinfo)
setReplaceMethod("profinfo", "xcmsSet", function(object, value) {
    object@profinfo <- value
    object
})

############################################################
## calibrate
setMethod("calibrate", "xcmsSet", function(object,calibrants,method="linear",
                                           mzabs=0.0001, mzppm=5,
                                           neighbours=3, plotres=FALSE) {

    nsamp = length(unique(object@peaks[,"sample"]))
    if (!sum(method == c("shift","linear","edgeshift")))
        stop("unknown calibration method!")

    if (is.list(calibrants))
        if (length(calibrants) != nsamp)
            stop("Error: Number of masslists differs with number of samples")

    for (s in 1:nsamp){
        peaklist = object@peaks[which(object@peaks[,"sample"]==s),]
        if (is.list(calibrants)) {
            masslist <- calibrants[s]
        }else{
            masslist <- calibrants
        }

        masses <- matchpeaks(peaklist,masslist,mzabs,mzppm,neighbours)
        if (length(masses)==0){
            warning("No masses close enough!\n")
            next
        }

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

############################################################
## groupnames
setMethod("groupnames", "xcmsSet", function(object, mzdec = 0, rtdec = 0,
                                            template = NULL) {

    if ( nrow(object@groups)<1 || length(object@groupidx) <1) {
        stop("No group information. Use group().")
    }

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

############################################################
## group.density
setMethod("group.density", "xcmsSet", function(object, bw = 30, minfrac = 0.5, minsamp = 1,
                                               mzwid = 0.25, max = 50, sleep = 0) {

    samples <- sampnames(object)
    classlabel <- sampclass(object)
    classnames <- as.character(unique(sampclass(object)))
    classlabel <- as.vector(unclass(classlabel))
    classnum <- table(classlabel)

    peakmat <- peaks(object)
    porder <- order(peakmat[,"mz"])
    peakmat <- peakmat[porder,, drop=FALSE]
    rownames(peakmat) <- NULL
    retrange <- range(peakmat[,"rt"])

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
        if (endidx - startidx < 0)
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

    groupmat <- groupmat[seq(length = num),,drop=FALSE]
    groupindex <- groupindex[seq(length = num)]

    ## Remove groups that overlap with more "well-behaved" groups
    numsamp <- rowSums(groupmat[,(match("npeaks", colnames(groupmat))+1):ncol(groupmat),drop=FALSE])
    uorder <- order(-numsamp, groupmat[,"npeaks"])

    uindex <- rectUnique(groupmat[,c("mzmin","mzmax","rtmin","rtmax"),drop=FALSE],
                         uorder)

    groups(object) <- groupmat[uindex,,drop=FALSE]
    groupidx(object) <- groupindex[uindex]

    object
})

############################################################
## group.mzClust
setMethod("group.mzClust", "xcmsSet", function(object,
                                               mzppm = 20,
                                               mzabs = 0,
                                               minsamp = 1,
                                               minfrac=0.5)
      {
          samples <- sampnames(object)
          classlabel <- sampclass(object)
          peaks <- peaks(object)
          groups <- mzClustGeneric(peaks[,c("mz","sample")],
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

          groups(object) <- cbind(groups$mat[,(1:3),drop=FALSE],rt,rt,rt,groups$mat[,(4:ncol(groups$mat)),drop=FALSE])
          colnames(groups(object)) <- c(colnames(groups$mat[,(1:3),drop=FALSE]), "rtmed", "rtmin", "rtmax", colnames(groups$mat[,(4:ncol(groups$mat)),drop=FALSE]))
          groupidx(object) <- groups$idx

          object
      })

############################################################
## group.nearest
setMethod("group.nearest", "xcmsSet", function(object, mzVsRTbalance=10,
                                               mzCheck=0.2, rtCheck=15, kNN=10) {

    ## If ANN is available ...
    ##RANN = "RANN"
    ##if (!require(RANN)) {
    ##    stop("RANN is not installed")
    ##}

    ## classlabel <- sampclass(object)
    classlabel <- as.vector(unclass(sampclass(object)))

    samples <- sampnames(object)
    gcount <- integer(length(unique(sampclass(object))))

    peakmat <- peaks(object)
    parameters <- list(mzVsRTBalance = mzVsRTbalance, mzcheck = mzCheck, rtcheck = rtCheck, knn = kNN)

    ptable <- table(peaks(object)[,"sample"])
    pord <- ptable[order(ptable, decreasing = TRUE)]
    sid <- as.numeric(names(pord))
    pn <- as.numeric(pord)

    samples <- sampnames(object)
    cat("sample:", basename(samples[sid[1]]), " ")

    mplenv <- new.env(parent = .GlobalEnv)
    mplenv$mplist <- matrix(0, pn[1], length(sid))
    mplenv$mplist[, sid[1]] <- which(peakmat[,"sample"] == sid[1])
    mplenv$mplistmean <- data.frame(peakmat[which(peakmat[,"sample"] == sid[1]),c("mz","rt")])
    mplenv$peakmat <- peakmat
    assign("peakmat", peakmat, envir = mplenv)

    sapply(sid[2:length(sid)], function(sample, mplenv, object){
        ## require(parallel)
        ## cl <- makeCluster(getOption("cl.cores", nSlaves))
        ## clusterEvalQ(cl, library(RANN))
        ## parSapply(cl, 2:length(samples), function(sample,mplenv, object){
        for(mml in seq(mplenv$mplist[,1])){
            mplenv$mplistmean[mml,"mz"] <- mean(mplenv$peakmat[mplenv$mplist[mml,],"mz"])
            mplenv$mplistmean[mml,"rt"] <- mean(mplenv$peakmat[mplenv$mplist[mml,],"rt"])
        }

        cat("sample:",basename(samples[sample])," ")
        mplenv$peakIdxList <- data.frame(peakidx=which(mplenv$peakmat[,"sample"]==sample),
                                         isJoinedPeak=FALSE)
        if(length(mplenv$peakIdxList$peakidx)==0){
            cat("Warning: No peaks in sample\n")
        }
        ## scoreList <- data.frame(score=numeric(0),peak=integer(0),mpListRow=integer(0),
        ##                                                       isJoinedPeak=logical(0), isJoinedRow=logical(0))
        ##
        ##       for(currPeak in mplenv$peakIdxList$peakidx){
        ##                       pvrScore <- patternVsRowScore(currPeak,parameters,mplenv) ## does the actual NN
        ##                       scoreList <- rbind(scoreList,pvrScore)
        ##       }
        ## this really doesn't take a long time not worth parallel version here.
        ## but make an apply loop now faster even with rearranging the data :D : PB
        scoreList <- sapply(mplenv$peakIdxList$peakidx, function(currPeak, para, mplenv){
            patternVsRowScore(currPeak,para,mplenv)
        }, parameters, mplenv)
        if(is.list(scoreList)){
            idx<-which(scoreList != "NULL")
            scoreList<-matrix(unlist(scoreList[idx]), ncol=5, nrow=length(idx), byrow=T)
            colnames(scoreList)<-c("score", "peak", "mpListRow", "isJoinedPeak", "isJoinedRow")
        } else {
            scoreList <- data.frame(score=unlist(scoreList["score",]), peak=unlist(scoreList["peak",]), mpListRow=
                                    unlist(scoreList["mpListRow",]), isJoinedPeak=unlist(scoreList["isJoinedPeak",]),
                                    isJoinedRow=unlist(scoreList["isJoinedRow",]))
        }

        ## Browse scores in order of descending goodness-of-fit
        scoreListcurr <- scoreList[order(scoreList[,"score"]),]
        if (nrow(scoreListcurr) > 0)
            for (scoreIter in 1:nrow(scoreListcurr)) {

                iterPeak <-scoreListcurr[scoreIter, "peak"]
                iterRow <- scoreListcurr[scoreIter, "mpListRow"]

                ## Check if master list row is already assigned with peak
                if (scoreListcurr[scoreIter, "isJoinedRow"]==TRUE) {
                    next
                }

                ## Check if peak is already assigned to some master list row
                if (scoreListcurr[scoreIter, "isJoinedPeak"]==TRUE) { next }

                ##  Check if score good enough
                ## Assign peak to master peak list row
                mplenv$mplist[iterRow,sample] <- iterPeak

                ## Mark peak as joined
                setTrue <- which(scoreListcurr[,"mpListRow"] == iterRow)
                scoreListcurr[setTrue,"isJoinedRow"] <- TRUE
                setTrue <- which(scoreListcurr[,"peak"] == iterPeak)
                scoreListcurr[setTrue, "isJoinedPeak"] <- TRUE
                mplenv$peakIdxList[which(mplenv$peakIdxList$peakidx==iterPeak),]$isJoinedPeak <- TRUE
            }

        notJoinedPeaks <- mplenv$peakIdxList[which(mplenv$peakIdxList$isJoinedPeak==FALSE),]$peakidx
        for(notJoinedPeak in notJoinedPeaks) {
            mplenv$mplist <- rbind(mplenv$mplist,matrix(0,1,dim(mplenv$mplist)[2]))
            mplenv$mplist[length(mplenv$mplist[,1]),sample] <- notJoinedPeak
        }

        ## Clear "Joined" information from all master peaklist rows
        rm(list = "peakIdxList", envir=mplenv)

        ## updateProgressInfo
        object@progressInfo$group.nearest <- (sample - 1) / (length(samples) - 1)
        progressInfoUpdate(object)
    }, mplenv, object)
    ## stopCluster(cl)
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

        gnum <- classlabel[unique(peakmat[mplenv$mplist[i,],"sample"])]
        for (j in seq(along = gcount))
            gcount[j] <- sum(gnum == j)
        groupmat[i, 7+seq(along = gcount)] <- gcount

        groupindex[[i]] <- mplenv$mplist[i, (which(mplenv$mplist[i,]>0))]
    }

    groups(object) <- groupmat
    groupidx(object) <- groupindex

    invisible(object)
})

############################################################
## group
setMethod("group", "xcmsSet", function(object, method=getOption("BioC")$xcms$group.method,
                                       ...) {
    method <- match.arg(method, getOption("BioC")$xcms$group.methods)
    if (is.na(method))
        stop("unknown method : ", method)
    method <- paste("group", method, sep=".")
    invisible(do.call(method, alist(object, ...)))
})

############################################################
## groupval
setMethod("groupval", "xcmsSet", function(object, method = c("medret", "maxint"),
                                          value = "index", intensity = "into") {

    if ( nrow(object@groups)<1 || length(object@groupidx) <1) {
        stop("xcmsSet is not been group()ed.")
    }

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

############################################################
## retcor
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

############################################################
## retcor.peakgroups
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
            xraw <- xcmsRaw(fnames[i])
            rtcor[[i]] <- xraw@scantime
        }
        object@rt <- list(raw = rtcor, corrected = rtcor)
    }

    nsamp <- rowSums(groupmat[,match("npeaks", colnames(groupmat))+unique(classlabel),drop=FALSE])

    idx <- which(nsamp >= n-missing & groupmat[,"npeaks"] <= nsamp + extra)
    if (length(idx) == 0)
        stop("No peak groups found for retention time correction")
    idx <- idx[order(groupmat[idx,"rtmed"])]

    rt <- groupval(object, "maxint", "rt")[idx,, drop=FALSE]
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

    ## Code for checking to see if retention time correction is overcorrecting
    rtdevrange <- range(rtdev, na.rm = TRUE)
    warn.overcorrect <- FALSE

    for (i in 1:n) {

        pts <- na.omit(data.frame(rt = rt[,i], rtdev = rtdev[,i]))

        if (smooth == "loess") {
            lo <- suppressWarnings(loess(rtdev ~ rt, pts, span = span, degree = 1, family = family))

            rtdevsmo[[i]] <- xcms:::na.flatfill(predict(lo, data.frame(rt = rtcor[[i]])))
### Remove singularities from the loess function
            rtdevsmo[[i]][abs(rtdevsmo[[i]]) > quantile(abs(rtdevsmo[[i]]), 0.9)*2] <- NA

            if (length(naidx <- which(is.na(rtdevsmo[[i]]))))
                rtdevsmo[[i]][naidx] <- suppressWarnings(approx(na.omit(data.frame(rtcor[[i]], rtdevsmo[[i]])),
                                                                xout = rtcor[[i]][naidx], rule = 2)$y)
            while (length(decidx <- which(diff(rtcor[[i]] - rtdevsmo[[i]]) < 0))) {
                d <- diff(rtcor[[i]] - rtdevsmo[[i]])[tail(decidx, 1)]
                rtdevsmo[[i]][tail(decidx, 1)] <- rtdevsmo[[i]][tail(decidx, 1)] - d
                if (abs(d) <= 1e-06)
                    break;
            }

            rtdevsmorange <- range(rtdevsmo[[i]])
            if (any(rtdevsmorange/rtdevrange > 2)) warn.overcorrect <- TRUE
        } else {
            if (nrow(pts) < 2) {
                stop("Not enough ``well behaved'' peak groups even for linear smoothing of retention times")
            }
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

        ## define the colors and line types and returns a list of
        ## mypal, col and ty. Uses the original code if no colors are
        ## submitted. Supports manually selected colors (e.g. in hex)
        vals <- defineColAndTy(col, ty, classlabel)
        col <- vals$col
        mypal <- vals$mypal
        ty <- vals$ty

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
        close.screen(all.screens = TRUE)
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

############################################################
## retcor.obiwarp
setMethod("retcor.obiwarp", "xcmsSet", function(object, plottype = c("none", "deviation"),
                                                profStep=1, center=NULL,
                                                col = NULL, ty = NULL,
                                                response=1, distFunc="cor_opt",
                                                gapInit=NULL, gapExtend=NULL,
                                                factorDiag=2, factorGap=1,
                                                localAlignment=0, initPenalty=0) {

    if (is.null(gapInit)) {
        if (distFunc=="cor") {gapInit=0.3}
        if (distFunc=="cor_opt") {gapInit=0.3}
        if (distFunc=="cov") {gapInit=0.0}
        if (distFunc=="euc") {gapInit=0.9}
        if (distFunc=="prd") {gapInit=0.0}
    }

    if (is.null(gapExtend)) {
        if (distFunc=="cor") {gapExtend=2.4}
        if (distFunc=="cor_opt") {gapExtend=2.4}
        if (distFunc=="cov") {gapExtend= 11.7}
        if (distFunc=="euc") {gapExtend= 1.8}
        if (distFunc=="prd") {gapExtend= 7.8}
    }

    peakmat <- peaks(object)
    samples <- sampnames(object)
    classlabel <- as.vector(unclass(sampclass(object)))
    N <- length(samples)
    corpeaks <- peakmat
    plottype <- match.arg(plottype)

    if (length(object@rt) == 2) {
        rtcor <- object@rt$corrected
    } else {
        fnames <- filepaths(object)
        rtcor <- vector("list", length(fnames))
        for (i in seq(along = fnames)) {
            xraw <- xcmsRaw(fnames[i])
            rtcor[[i]] <- xraw@scantime
        }
        object@rt <- list(raw = rtcor, corrected = rtcor)
    }

    rtimecor <- vector("list", N)
    rtdevsmo <- vector("list", N)
    plength <- rep(0, N)

    if (missing(center)) {
        for(i in 1:N){
            plength[i] <- length(which(peakmat[,"sample"]==i))
        }
        center <- which.max(plength)
    }

    cat("center sample: ", samples[center], "\nProcessing: ")
    idx <- which(seq(1,N) != center)
    obj1 <- xcmsRaw(object@filepaths[center], profmethod="bin", profstep=0)

	## added t automatically find the correct scan range from the xcmsSet object
	if(length(obj1@scantime) != length(object@rt$raw[[center]])){
		##figure out the scan time range
		scantime.start	<-object@rt$raw[[center]][1]
		scantime.end	<-object@rt$raw[[center]][length(object@rt$raw[[center]])]

		scanrange.start	<-which.min(abs(obj1@scantime - scantime.start))
		scanrange.end	<-which.min(abs(obj1@scantime - scantime.end))
		scanrange<-c(scanrange.start, scanrange.end)
		obj1 <- xcmsRaw(object@filepaths[center], profmethod="bin", profstep=0, scanrange=scanrange)
	} else{
		scanrange<-NULL
	}

    for (si in 1:length(idx)) {
        s <- idx[si]
        cat(samples[s], " ")

        profStepPad(obj1) <- profStep ## (re-)generate profile matrix, since it might have been modified during previous iteration
		if(is.null(scanrange)){
			obj2 <- xcmsRaw(object@filepaths[s], profmethod="bin", profstep=0)
		} else{
			obj2 <- xcmsRaw(object@filepaths[s], profmethod="bin", profstep=0, scanrange=scanrange)
		}
        profStepPad(obj2) <- profStep ## generate profile matrix

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

        if ((mzval * valscantime1 != length(intensity1)) ||  (mzval * valscantime2 != length(intensity2)))
            stop("Dimensions of profile matrices do not match !\n")

        rtimecor[[s]] <-.Call("R_set_from_xcms",
                              valscantime1,scantime1,mzval,mz,intensity1,
                              valscantime2,scantime2,mzval,mz,intensity2,
                              response, distFunc,
                              gapInit, gapExtend,
                              factorDiag, factorGap,
                              localAlignment, initPenalty)

        if(length(obj2@scantime) > valscantime2) {
            object@rt$corrected[[s]] <- c(rtimecor[[s]],
                                          obj2@scantime[(max(rtup2)+1):length(obj2@scantime)])
        } else {
            object@rt$corrected[[s]] <- rtimecor[[s]]
        }

        rtdevsmo[[s]] <- round(rtcor[[s]]-object@rt$corrected[[s]],2)

        rm(obj2)
        gc()

        ## updateProgressInfo
        object@progressInfo$retcor.obiwarp <-  si / length(idx)
        xcms:::progressInfoUpdate(object)

    }

    cat("\n")
    rtdevsmo[[center]] <- round(rtcor[[center]] - object@rt$corrected[[center]], 2)

    if (plottype == "deviation") {

        ## define the colors and line types and returns a list of
        ## mypal, col and ty. Uses the original code if no colors are
        ## submitted. Supports manually selected colors (e.g. in hex)
        vals <- defineColAndTy(col, ty, classlabel)
        col <- vals$col
        mypal <- vals$mypal
        ty <- vals$ty

        rtrange <- range(do.call("c", rtcor))
        devrange <- range(do.call("c", rtdevsmo))

        layout(matrix(c(1, 2),ncol=2,  byrow=F),widths=c(1,0.3))
        par(mar=c(4,4,2,0))

        plot(0, 0, type="n", xlim = rtrange, ylim = devrange,
             main = "Retention Time Deviation vs. Retention Time",
             xlab = "Retention Time", ylab = "Retention Time Deviation")

        for (i in 1:N) {
            points(rtcor[[i]], rtdevsmo[[i]], type="l", col = mypal[col[i]], lty = ty[i])
        }

        plot.new() ;  par(mar= c(2, 0, 2, 0))
        plot.window(c(0,1), c(0,1))
        legend(0,1.04, basename(samples), col = mypal[col], lty = ty)
    }

    for (i in 1:N) {
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

############################################################
## plotrt
setMethod("plotrt", "xcmsSet", function(object, col = NULL, ty = NULL, leg = TRUE, densplit = FALSE) {

    samples <- sampnames(object)
    classlabel <- as.vector(unclass(sampclass(object)))
    n <- length(samples)
    rtuncor <- object@rt$raw
    rtcor <- object@rt$corrected

    ## define the colors and line types and returns a list of
    ## mypal, col and ty. Uses the original code if no colors are
    ## submitted. Supports manually selected colors (e.g. in hex)
    vals <- defineColAndTy(col, ty, classlabel)
    col <- vals$col
    mypal <- vals$mypal
    ty <- vals$ty

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
        close.screen(all.screens = TRUE)
    }
})

############################################################
## fillPeaks.chrom
## New version using BiocParallel instead of nSlaves and manual setup.
setMethod("fillPeaks.chrom", "xcmsSet", function(object, nSlaves = NULL,
                                                  expand.mz = 1,expand.rt = 1,
                                                  BPPARAM = bpparam()) {
  ## development mockup:
  if (FALSE) {
    library(xcms)
    library(faahKO)
    object <- group(faahko)
    gf <- fillPeaks(object)
    pkgEnv = getNamespace("xcms")
    attach(pkgEnv)
  }

    if (!is.null(nSlaves))
        warning("Use of argument 'nSlaves' is deprecated!",
                " Please use 'BPPARAM' instead.")

    peakmat <- peaks(object)
    groupmat <- groups(object)
    if (length(groupmat) == 0)
        stop("No group information found")
    files <- filepaths(object)
    samp <- sampnames(object)
    classlabel <- as.vector(unclass(sampclass(object)))
    prof <- profinfo(object)
    rtcor <- object@rt$corrected

    ## Remove groups that overlap with more "well-behaved" groups
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
    lastpeakOrig <- lastpeak

    ##    peakmat <- rbind(peakmat, matrix(nrow = sum(is.na(gvals)), ncol = ncol(peakmat)))

    cnames <- colnames(object@peaks)

    ## Making gvals environment so that when it is repeated for each file it only uses the memory one time
    gvals_env <- new.env(parent=baseenv())
    assign("gvals", gvals, envir = gvals_env)

    ft <- cbind(file=files,id=1:length(files))
    argList <- apply(ft,1,function(x) {
        ## Add only those samples which actually have NA in them
        if (!any(is.na(gvals[,as.numeric(x["id"])]))) {
        ## nothing to do.
            list()
        } else {
            list(file=x["file"],id=as.numeric(x["id"]),
                 params=list(method="chrom",
                             gvals=gvals_env,
                             prof=prof,
                             dataCorrection=object@dataCorrection,
                             polarity=object@polarity,
                             rtcor=object@rt$corrected[[as.numeric(x["id"])]],
                             peakrange=peakrange,
                             expand.mz=expand.mz,
                             expand.rt=expand.rt))
        }
    })

    nonemptyIdx <- (sapply(argList, length) > 0)

    if (!any(nonemptyIdx)) {
        ## Nothing to do
        return(invisible(object))
    }

    argList <- argList[nonemptyIdx]

    ## Use BiocParallel for parallel computing.
    newpeakslist <- bplapply(argList, fillPeaksChromPar, BPPARAM = BPPARAM)

    o <- order(sapply(newpeakslist, function(x) x$myID))
    newpeaks <- do.call(rbind, lapply(newpeakslist[o], function(x) x$newpeaks))

    ## Make sure colnames are compatible
    newpeaks <- newpeaks[, match(cnames, colnames(newpeaks)), drop = FALSE]
    colnames(newpeaks) <- cnames

    peakmat <- rbind(peakmat, newpeaks)

    for (i in seq(along = files)) {
        naidx <- which(is.na(gvals[,i]))

        for (j in seq(along = naidx))
            groupindex[[naidx[j]]] <- c(groupindex[[naidx[j]]], lastpeak+j)

        lastpeak <- lastpeak + length(naidx)
    }

    peaks(object) <- peakmat
    object@filled <- seq((lastpeakOrig+1),nrow(peakmat))
    groups(object) <- groupmat
    groupidx(object) <- groupindex

    invisible(object)
})

############################################################
## fillPeaks.MSW
setMethod("fillPeaks.MSW", "xcmsSet", function(object, mrange=c(0,0), sample=NULL) {

    peakmat <- peaks(object)
    groupmat <- groups(object)
    if (length(groupmat) == 0)
        stop("No group information found")
    files <- filepaths(object)
    samp <- sampnames(object)
    classlabel <- as.vector(unclass(sampclass(object)))
    prof <- profinfo(object)
    rtcor <- object@rt$corrected
    ## Remove groups that overlap with more "well-behaved" groups
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

############################################################
## fillPeaks
setMethod("fillPeaks", "xcmsSet", function(object, method=getOption("BioC")$xcms$fillPeaks.method,...) {
    method <- match.arg(method, getOption("BioC")$xcms$fillPeaks.methods)
    if (is.na(method))
        stop("unknown method : ", method)
    method <- paste("fillPeaks", method, sep=".")
    invisible(do.call(method, alist(object, ...)))
})

############################################################
## getEIC
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
        stop("rtrange must be a matrix or single number")
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
        ## getXcmsRaw takes care of rt correction, susetting to scanrage and other
        ## stuff.
        lcraw <- getXcmsRaw(object, sampleidx=i, rt=rt)
        currenteic <- getEIC(lcraw, mzrange, rtrange, step = prof$step)
        eic[[i]] <- currenteic@eic[[1]]
        rm(lcraw)
        gc()
    }
    cat("\n")

    invisible(new("xcmsEIC", eic = eic, mzrange = mzrange, rtrange = rtrange,
                  rt = rt, groupnames = gnames))
})

############################################################
## peakTable
setMethod("peakTable", "xcmsSet", function(object, filebase = character(), ...) {

    if (length(sampnames(object)) == 1) {
        return(object@peaks)
    }

    if (nrow(object@groups) < 1) {
        stop ('First argument must be an xcmsSet with group information or contain only one sample.')
    }

    groupmat <- groups(object)


    if (! "value" %in% names(list(...))) {
        ts <- data.frame(cbind(groupmat,groupval(object, value="into",  ...)), row.names = NULL)
    } else {
        ts <- data.frame(cbind(groupmat,groupval(object, ...)), row.names = NULL)
    }

    cnames <- colnames(ts)

    if (cnames[1] == 'mzmed') {
        cnames[1] <- 'mz'
    } else {
        stop ('mzmed column missing')
    }
    if (cnames[4] == 'rtmed') {
        cnames[4] <- 'rt'
    } else {
        stop ('mzmed column missing')
    }

    colnames(ts) <- cnames

    if (length(filebase))
        write.table(ts, paste(filebase, ".tsv", sep = ""), quote = FALSE, sep = "\t", col.names = NA)

    ts
})

############################################################
## diffreport
setMethod("diffreport", "xcmsSet", function(object, class1 = levels(sampclass(object))[1],
                                            class2 = levels(sampclass(object))[2],
                                            filebase = character(), eicmax = 0, eicwidth = 200,
                                            sortpval = TRUE, classeic = c(class1,class2),
                                            value = c("into","maxo","intb"), metlin = FALSE,
                                            h = 480, w = 640, mzdec=2, ...) {

    if ( nrow(object@groups)<1 || length(object@groupidx) <1) {
        stop("No group information. Use group().")
    }

    ## require(multtest) || stop("Couldn't load multtest")

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

    ## Check against missing Values
    if (any(is.na(values[,c(c1,c2)]))) {
        stop("NA values in xcmsSet. Use fillPeaks()")
    }

    mean1 <- rowMeans(values[,c1,drop=FALSE], na.rm = TRUE)
    mean2 <- rowMeans(values[,c2,drop=FALSE], na.rm = TRUE)

    ## Calculate fold change.
    ## For foldchange <1 set fold to 1/fold
    ## See tstat to check which was higher
    fold <- mean2 / mean1
    fold[!is.na(fold) & fold < 1] <- 1/fold[!is.na(fold) & fold < 1]

    testval <- values[,c(c1,c2)]
    testclab <- c(rep(0,length(c1)),rep(1,length(c2)))

    if (min(length(c1), length(c2)) >= 2) {
        tstat <- mt.teststat(testval, testclab, ...)
        pvalue <- pval(testval, testclab, tstat)
    } else {
        message("Too few samples per class, skipping t-test.")
        tstat <- pvalue <- rep(NA,nrow(testval))
    }
    stat <- data.frame(fold = fold, tstat = tstat, pvalue = pvalue)
    if (length(levels(sampclass(object))) >2) {
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
            plot(eics, object, rtrange = eicwidth, mzdec=mzdec)

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

############################################################
## progressCallback
setMethod("progressCallback", "xcmsSet", function(object) object@progressCallback)
setReplaceMethod("progressCallback", "xcmsSet", function(object, value) {
    object@progressCallback <- value
    object
})

############################################################
## progressInfoUpdate
setMethod("progressInfoUpdate", "xcmsSet", function(object)
          do.call(object@progressCallback, list(progress=object@progressInfo))
          )

############################################################
## $
## eSet/ExpressionSet like methods to access columns in @phenoData
setMethod("$", "xcmsSet", function(x, name) {
    eval(substitute(phenoData(x)$NAME_ARG, list(NAME_ARG=name)))
})
setReplaceMethod("$", "xcmsSet", function(x, name, value) {
  phenoData(x)[[name]] = value
  x
})

############################################################
## getXcmsRaw
## read the raw data for a xset.
## argument which allows to specify which file from the xcmsSet should be read, if length > 1
## a list of xcmsRaw is returned
setMethod("getXcmsRaw", "xcmsSet", function(object, sampleidx = 1,
                                            profmethod = profMethod(object),
                                            profstep = profStep(object),
                                            profparam = profinfo(object),
                                            mslevel = NULL,
                                            scanrange = NULL,
                                            rt = c("corrected", "raw"),
                                            BPPARAM = bpparam()) {
              if (is.numeric(sampleidx))
                  sampleidx <- sampnames(object)[sampleidx]
              sampidx <- match(sampleidx, sampnames(object)) ## numeric
              if (length(sampidx) == 0)
                  stop("submitted value for sampleidx outside of the available files!")
              fn <- filepaths(object)[sampidx]
              rt <- match.arg(rt)
              if (rt == "corrected" & !any(names(object@rt) == "corrected")) {
                  message("No RT correction has been performed, thus returning raw",
                          " retention times.")
                  rt <- "raw"
              }
              if (missing(mslevel)) {
                  msl <- mslevel(object)
              } else {
                  msl <- mslevel
              }
              if (missing(scanrange)) {
                  srange <- scanrange(object)
              } else {
                  srange <- scanrange
              }
              ## include MSn?
              includeMsn <- FALSE
              if (!is.null(msl)) {
                  if(msl > 1)
                      includeMsn <- TRUE
              }
              ret <- bplapply(as.list(fn), FUN = function(z) {
                                  raw <- xcmsRaw(z, profmethod=profmethod, profstep=profstep,
                                                 mslevel=msl, scanrange=srange,
                                                 includeMSn=includeMsn)
                                  return(raw)
                            }, BPPARAM=BPPARAM)
              ## do corrections etc.
              for(i in 1:length(ret)){
                  if(length(object@dataCorrection) > 1){
                      if(object@dataCorrection[[sampidx[i]]] == 1){
                          ret[[i]] <- stitch(ret[[i]], AutoLockMass(ret[[i]]))
                          message(paste0("Applying lock Waters mass correction to ", fn[i]))
                      }
                  }
                  if(rt == "corrected"){
                      ## check if there is any need to apply correction...
                      ## This includes fix for the bug reported by Aleksandr (issue 44)
                      if(all(object@rt$corrected[[sampidx[i]]] == object@rt$raw[[sampidx[i]]])){
                          message("No need to perform retention time correction,",
                                  " raw and corrected rt are identical for ", fn[i])
                          ret[[i]]@scantime <- object@rt$raw[[sampidx[i]]]
                      }else{
                          message(paste0("Applying retention time correction to ", fn[i]))
                          ret[[i]]@scantime <- object@rt$corrected[[sampidx[i]]]
                      }
                  }
              }

              ## what's missing?
              ## + consider calibration(s)?
              ## + what with polarity?

              if(length(ret)==1)
                  return(ret[[1]])
              return(ret)
          })

############################################################
## levelplot
setMethod("levelplot", "xcmsSet",
          function(x, log=TRUE, sampleidx=1,
                   col.regions=colorRampPalette(brewer.pal(9, "YlOrRd"))(256),
                   highlight.peaks=FALSE, highlight.col="#377EB880",
                   rt="raw", ...){
              if (is.numeric(sampleidx))
                  sampleidx <- sampnames(x)[sampleidx]
              sampidx <- match(sampleidx, sampnames(x))
              if(length(sampidx) > 1)
                  stop("A levelplot can only be created for one file at a time!")
              xraw <- getXcmsRaw(x, sampleidx=sampidx, rt=rt)
              if(highlight.peaks){
                  pks <- peaks(x)
                  pks <- pks[pks[, "sample"]==sampidx, ]
                  ## call the levelplot AND specify the panel...
                  ## some code taken from plotSurf...
                  sel <- profRange(xraw, ...)
                  zvals <- xraw@env$profile[sel$massidx, sel$scanidx]
                  if(log){
                      zvals <- log(zvals+max(c(-min(zvals), 1)))
                  }
                  ## y axis is time
                  yvals <- xraw@scantime[sel$scanidx]
                  yrange <- range(yvals)
                  ## x is m/z
                  xvals <- profMz(xraw)[sel$massidx]
                  ## now i have to match the x and y to the z.
                  ## as.numeric of z returns values by column(!), i.e. the first nrow(z) correspond to
                  ## the x of 1.
                  xvals <- rep(xvals, ncol(zvals))
                  yvals <- rep(yvals, each=nrow(zvals))
                  zvals <- as.numeric(zvals)
                  ## get the file name
                  fileName <- sampnames(x)[sampidx]
                  plt <- levelplot(zvals~xvals*yvals,
                                   panel=function(...){
                                       panel.levelplot(...)
                                       panel.rect(xleft=pks[, "mzmin"], xright=pks[, "mzmax"],
                                                ybottom=pks[, "rtmin"], ytop=pks[, "rtmax"],
                                                border=highlight.col)
                                   },
                                   xlab="m/z", ylab="Time",
                                   colorkey=list(height=1, width=0.7),
                                   main=list(fileName, side=1, line=0.5), col.regions=col.regions)
                  ## plt <- plt + layer(panel.rect(xleft=pks[, "mzmin"], xright=pks[, "mzmax"],
                  ##                               ybottom=pks[, "rtmin"], ytop=pks[, "rtmax"],
                  ##                               border="#377EB820")
                  ##                   )
              }else{
                  plt <- levelplot(xraw, col.regions = col.regions, log=log)
              }
              plt
          })

############################################################
## mslevel
## getter methods for the slots mslevel and scanrange of the xcmsSet object.
setMethod("mslevel", "xcmsSet", function(object){
              ## for xcmsSet objects that don't have (yet) the slot...
              if(!.hasSlot(object, "mslevel")){
                  message("No slot mslevel available, consider updating the",
                          " object with the 'updateObject' method.")
                  return(NULL)
              }else{
                  mlevel <- object@mslevel
                  if(length(mlevel) == 0){
                      ## for compatibility with the xcmsSet and xcmsRaw functions,
                      ## which default the mslevel argument to NULL.
                      return(NULL)
                  }
                  return(mlevel)
              }
          })
setReplaceMethod("mslevel", "xcmsSet", function(object, value){
                     if(.hasSlot(object, "mslevel")){
                         object@mslevel <- value
                     }else{
                         message("Object has no slot mslevel, consider updating",
                                 " the object using 'updateObject'.")
                     }
                     object
                 })

############################################################
## scanrange
setMethod("scanrange", "xcmsSet", function(object){
              if(!.hasSlot(object, "scanrange")){
                  message("No slot scanrange available, consider updating the",
                          " object with the 'updateObject' method.")
                  return(NULL)
              }else{
                  srange <- object@scanrange
                  if(length(srange) == 0){
                      ## for compatibility with the xcmsSet and xcmsRaw functions,
                      ## which default the scanrange argument to NULL.
                      return(NULL)
                  }
                  return(srange)
              }
          })
setReplaceMethod("scanrange", "xcmsSet", function(object, value){
                     if(.hasSlot(object, "scanrange")){
                         object@scanrange <- value
                     }else{
                         message("Object has no slot scanrange,  consider updating",
                          " the object with the 'updateObject' method.")
                     }
                     object
                 })

############################################################
## profMethod
## getter methods for the prof method and step.
setMethod("profMethod", "xcmsSet", function(object) {
              return(profinfo(object)$method)
})
setMethod("profStep", "xcmsSet", function(object) {
              return(profinfo(object)$step)
})

## sub setting an xcmsSet object...
setMethod("[", "xcmsSet", function(x, i, j, ..., drop = FALSE) {
    if (missing(drop))
        drop <- FALSE
    if (missing(i) && missing(j)) {
        if (length(list(...))!=0)
            stop("specify samples to subset; use '",
                 substitute(x), "$", names(list(...))[[1]],
                 "' to access phenoData variables")
        return(x)
    }
    ## don't allow i, but allow j to be: numeric or logical. If
    ## it's a character vector <- has to fit to sampnames(x)
    if(!missing(i))
        stop("Subsetting to rows is not supported!")
    if(missing(j))
        j <- 1:length(sampnames(x))
    if(class(j)=="character"){
        ## check if these match to the sampnames.
        matches <- match(j, sampnames(x))
        if(length(matches)!=length(j))
            stop("All provided sample names have to match the",
                 " sample names in the xcmsSet!")
        j <- matches
    }
    if(class(j)=="logical"){
        if(length(j) != length(sampnames(x)))
            stop("If j is a logical its length has to match",
                 " the number of samples in the xcmsSet!")
        j <- which(j)
    }
    if(class(j)=="numeric")
        j <- as.integer(j)
    if(class(j)!="integer")
        stop("j has to be a numeric vector specifying the",
             " index of the samples for which the data has to be extracted")
    ## check if j is within the range of 1:length(sampnames)
    if(any(!j %in% (1:length(sampnames(x)))))
        stop("j has to be a numeric with values between 1 and ",
             length(sampnames(x)), "!")
    ## OK, j is now an integer vector...
    ## "copy" the xcmsSet, that way we keep parameters mslevel, scanrange,
    ## profinfo, polarity.
    xsub <- x
    ## first of all, subset the phenoData
    phenoData(xsub) <- droplevels(phenoData(xsub)[j,, drop=FALSE])
    ## then the file paths
    filepaths(xsub) <- filepaths(x)[j]
    ## now starting to subset data:
    ## 1) @rt$raw, @rt$corrected
    xsub@rt$raw <- x@rt$raw[j]
    xsub@rt$corrected <- x@rt$corrected[j]
    ## 2) @peaks
    keep.peaks <- x@peaks
    rownames(keep.peaks) <- as.character(1:nrow(keep.peaks))
    ## subsetting the peaks. Since we want to also allow reverse ordering we have to
    ## do it a little more complicated...
    keep.peaks <- split(data.frame(keep.peaks), f=keep.peaks[, "sample"])
    keep.peaks <- keep.peaks[as.character(j)]
    names(keep.peaks) <- NULL
    keep.peaks <- as.matrix(do.call(rbind, keep.peaks))
    ##keep.peaks <- keep.peaks[as.character(keep.peaks[, "sample"]) %in% as.character(j), ]
    ## have to replace the sample index.
    newsample <- numeric(nrow(keep.peaks))
    for(idx in 1:length(j)){
        newsample[keep.peaks[, "sample"]==j[idx]] <- idx
    }
    keep.peaks[, "sample"] <- newsample
    xsub@peaks <- keep.peaks
    rownames(xsub@peaks) <- NULL
    ## 3) groupidx if present. subset this to indices of peaks which we will keep.
    ##    will use the rownames of the keep.peaks matrix for that (represents the
    ##    original index).
    if(length(x@groupidx) > 0){
        keep.groupidx <- lapply(x@groupidx, function(z){
            newidx <- match(as.character(z), rownames(keep.peaks))
            return(newidx[!is.na(newidx)])
            ## that way I just return the indices as they are!
            ##return(z[as.character(z) %in% rownames(keep.peaks)])
        })
        xsub@groupidx <- keep.groupidx
        ## 4) groups have to be re-calculated.
        keep.groups <- x@groups
        keep.groups <- keep.groups[, -(which(colnames(keep.groups)=="npeaks"):ncol(keep.groups))]
        sampclasses <- sampclass(xsub)
        peakCounts <- matrix(ncol=length(levels(sampclasses)), nrow=nrow(keep.groups), 0)
        colnames(peakCounts) <- levels(sampclasses)
        ## loop throught the peakcounts
        for(idx in 1:nrow(keep.groups)){
            groupidx <- keep.groupidx[[idx]]
            if(length(groupidx) == 0)
                next
            tab <- keep.peaks[groupidx, "sample"]
            tab <- table(sampclasses[tab])
            peakCounts[idx, names(tab)] <- as.numeric(tab)
        }
        keep.groups <- cbind(keep.groups, npeaks=rowSums(peakCounts), peakCounts)
        xsub@groups <- keep.groups
    }
    ## 5) filled
    if(length(x@filled) > 0){
        xsub@filled <- (1:nrow(keep.peaks))[rownames(keep.peaks) %in% as.character(x@filled)]
    }
    ## 6) dataCorrection
    if(length(x@dataCorrection)>0)
        xset@dataCorrection <- x@dataCorrection[j]
    return(xsub)
})

############################################################
## present
setMethod("present", "xcmsSet", function(object, class, minfrac) {
    if ( nrow(object@groups)<1 || length(object@groupidx) <1) {
        stop("No group information. Use group().")
    }

    classlabel <- sampclass(object)
    classlabel <- levels(classlabel)[as.vector(unclass(classlabel))]

    sampidx <- which(classlabel %in% class)

    if (length(sampidx) == 0) {
        stop("Class ", class, "not found")
    }

    classnum <- length(sampidx)
    minpresent <- classnum * minfrac

    filled <- rep(FALSE, nrow(peaks(object)))

    ## exists(object@filled) always returns FALSE ??
    sloti <- try(slot(object, "filled"), silent = TRUE)
    if (class(sloti) != "try-error") {
        filled[object@filled] <- TRUE
    }

    apply (groupval(object), 1, function(x) {
        length(which(  (!(is.na(x[sampidx]) | is.nan(x[sampidx])))
                     & !filled[x[sampidx]])) >= minpresent
    })
})

############################################################
## absent
setMethod("absent", "xcmsSet", function(object, class, minfrac) {
    if ( nrow(object@groups)<1 || length(object@groupidx) <1) {
        stop("No group information. Use group().")
    }

    classlabel <- sampclass(object)
    classlabel <- levels(classlabel)[as.vector(unclass(classlabel))]

    sampidx <- which(classlabel %in% class)

    if (length(sampidx) == 0) {
        stop("Class ", class, "not found")
    }

    classnum <- length(sampidx)
    minabsent <- classnum * minfrac

    filled <- rep(FALSE, nrow(peaks(object)))

    ## exists(object@filled) always returns FALSE ??
    sloti <- try(slot(object, "filled"), silent = TRUE)
    if (class(sloti) != "try-error") {
        filled[object@filled] <- TRUE
    }

    apply (groupval(object), 1, function(x) {
        length(which(is.na(x[sampidx]) | is.nan(x[sampidx]) | filled[x[sampidx]])) >= minabsent
    })
})

############################################################
## specDist
setMethod("specDist", signature(object="xcmsSet"),
          function(object, peakIDs1, peakIDs2,
                   method=getOption("BioC")$xcms$specDist.method,
                   ...) {
              if (missing(peakIDs1)) {
                  stop("missing argument peakIDs1")
              }
              if (missing(peakIDs2)) {
                  stop("missing argument peakIDs2")
              }

              peaks <- object@peaks
              peakTable1 <- peaks[peakIDs1,c("mz","into")]
              peakTable2 <- peaks[peakIDs2,c("mz","into")]

              method <- match.arg(method, getOption("BioC")$xcms$specDist.methods)
              if (is.na(method))
                  stop("unknown method : ", method)
              method <- paste("specDist", method, sep=".")
              distance <- do.call(method, alist<-list(peakTable1, peakTable2, ...))
              distance
          })
